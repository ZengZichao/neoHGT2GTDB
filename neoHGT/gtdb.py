"""
////////////////////////////////////////////////////////////
│
│	Filename:	gtdb.py
│	Comment:	GTDB data source provider: local-layout probing, input
│			self-consistency checks, GTDB -> NCBI-compatible taxdump
│			conversion with stable synthetic taxIds, and the sequence
│			ingestion channels used by "neoHGT database --source gtdb".
│
│	=======================================================
│	Copyright:	Distributed under the terms of the Modified BSD License.
│
////////////////////////////////////////////////////////////
"""
#!/usr/bin/env python3

import re
import os
import io
import json
import gzip
import bz2
import lzma
import hashlib
import statistics
import sys
from bisect import bisect_right
from os import listdir, makedirs
from os.path import (
    join, isfile, isdir, basename, dirname, realpath, exists, splitext,
    getsize)
from datetime import datetime
from collections import OrderedDict

import pandas as pd

from neoHGT.util import (
    GTDB_RANKS, GTDB_RANK_PREFIXES, NCBI_RANKS, timestamp,
    strip_rank_prefix, rank_of_label, rank_normalize, read_taxdump,
    taxid_at_rank, find_lca, find_config_file)


# ---------------------------------------------------------------------------
# Constants (implementation contract, see plan §3.3 and §5.3)
# ---------------------------------------------------------------------------

SERIAL_ANCHOR = '_SerialNumber_'          # local protein-Id anchor (§3.3d)
STEM_ANCHOR_RE = re.compile(r'_(d|p|c|o|f|g)_')
LOCAL_LAYOUTS = ('local-renamed', 'official-plain', 'official-gz')
# Layout label written by builds made before the label stopped encoding a
# release number; still accepted everywhere a renamed layout is expected.
LEGACY_RENAMED_LAYOUT = 'gtdb-r220-renamed'
# Start of the reserved GTDB taxId range. Plan ADR-2 proposed 1e9, but spike
# S1 measured that DIAMOND sizes db.dmnd by the highest taxId value, not by the
# number of nodes: the same 13 MB test database became 5.7 GB at 1e9 and 620 MB
# at 1e8. 1e8 keeps a ~25x margin over the largest NCBI taxId measured in the
# reference dumps this project sees (3,905,377 in a 2026 r220-era dump; the
# "~3.4e6" figure quoted by the first comment revision is stale). Re-measure on
# every release: if NCBI crosses 1e7 the whole synthetic range must move.
DEFAULT_ID_BASE = 100_000_000
MIN_SAFE_ID_BASE = 10_000_000              # re-check against real dumps per release


def tool_version():
    """Version string recorded in meta.json / provenance lines.

    Reads the package version so the fork identity cannot drift from
    neoHGT.__version__ (two hardcoded literals used to say '2.0b3+gtdb1'
    while the package still said '2.0b3').
    """
    try:
        from neoHGT import __version__ as v
    except Exception:
        return 'unknown'
    return v


def check_id_base(id_base):
    """Refuse synthetic taxId ranges that could collide with real NCBI taxIds.

    Every entry point that mints GTDB taxIds must call this: the validation
    used to live only on the `database` stage, so `gtdb taxdump
    --gtdb-taxid-base 1000` could silently emit ids inside the real NCBI
    range and scramble any mixed-source comparison.

    Parameters
    ----------
    id_base : int
        first taxId of the synthetic range

    Raises
    ------
    ValueError
        when the base falls below MIN_SAFE_ID_BASE
    """
    id_base = int(id_base)
    if id_base < MIN_SAFE_ID_BASE:
        raise ValueError(
            f'--gtdb-taxid-base {id_base} is too low: it would collide with '
            f'real NCBI taxIds. Use at least {MIN_SAFE_ID_BASE}.')
    return id_base
ID_SLOT = 1_000_000                        # ids reserved per rank
MAX_TAXID = 2 ** 31 - 1                    # int32 ceiling assumed by 3rd parties
ROOT_TID = '1'
CELLULAR_ROOT_TID = '131567'
RANK_ORDINAL = {r: i for i, r in enumerate(GTDB_RANKS)}
PREFIX_BY_CHAR = {'d': 'domain', 'p': 'phylum', 'c': 'class', 'o': 'order',
                  'f': 'family', 'g': 'genus', 's': 'species'}
LINEAGE_SEP = ';'
DEPTH = len(GTDB_RANKS)                    # GTDB lineages are always complete

# Process-local caches: these inputs are immutable during a build.
_METADATA_CACHE = {}
_TAXONOMY_CACHE = {}


def provenance(info=None, release=None, extra=None):
    """Provenance header line for derived files (plan §5 principle 4).

    Parameters
    ----------
    info : dict, optional
        layout information from detect_layout (supplies release and layout)
    release : str, optional
        release label to record when info is not given
    extra : dict, optional
        additional ``key=value`` fields, written in insertion order

    Returns
    -------
    str
        a single ``#``-comment line ending in a newline

    Notes
    ----
    Every text artifact written under ``gtdb/`` starts with this line, so a file
    can always be traced back to the code version, the GTDB release and the
    moment it was produced. The two maps handed to external tools verbatim
    (``taxon.map.gz``, ``genome.map.gz``) are the deliberate exception: DIAMOND,
    BLAST and ``read_prot2taxid`` consume them byte-for-byte.
    """
    from neoHGT import __version__
    rel = release if release is not None else (info or {}).get('release')
    parts = [f'neoHGT version={__version__}']
    if rel is not None:
        parts.append(f'gtdb_release={rel}')
    if info and info.get('layout'):
        parts.append(f'layout={info["layout"]}')
    parts.append(f'generated={timestamp()}')
    for key, value in (extra or {}).items():
        parts.append(f'{key}={value}')
    return '# ' + '; '.join(parts) + '\n'


# Columns neoHGT needs from GTDB metadata. Resolved BY NAME, never by index.
REQUIRED_COLS = ('accession', 'gtdb_taxonomy', 'gtdb_representative',
                 'gtdb_genome_representative', 'checkm2_completeness',
                 'checkm2_contamination', 'ncbi_assembly_level',
                 'ncbi_refseq_category', 'ncbi_taxid', 'ncbi_taxonomy',
                 'protein_count')
OPTIONAL_COLS = ('gtdb_type_species_of_genus', 'gtdb_type_designation_ncbi_taxa',
                 'ncbi_genbank_assembly_accession', 'ncbi_species_taxid',
                 'ncbi_organism_name', 'ncbi_type_material_designation',
                 'mimag_high_quality', 'mimag_medium_quality',
                 'mimag_low_quality', 'checkm2_model', 'contig_count',
                 'n50_contigs', 'genome_size')

# Files that record how this local copy was produced (plan §3.3f / ADR-9).
PROVENANCE_PATTERNS = ('00_rename_', '01_sequence_prefix_',
                       'Merge_GTDB207_', '.tar.sh')

WORKSPACE_ENVS = ('NEOHGT_ROOT', 'NEOHGT_CODE', 'NEOHGT_CACHE',
                  'NEOHGT_TMP')
DEFAULT_WORKSPACE = None


def is_renamed_layout(label):
    """Whether a layout label denotes a lineage-rewritten per-genome FASTA set.

    The label used to hardcode the release it was first seen on
    ("gtdb-r220-renamed"), which made every other release claim to be r220.
    Layouts are now named by category ("local-renamed"); both spellings are
    accepted so existing meta.json artifacts keep working.
    """
    return str(label).endswith('-renamed')


class SelfCheckError(RuntimeError):
    """Raised when the local dataset fails a consistency check."""


class GuardError(RuntimeError):
    """Raised when a write target escapes the allowed workspace."""


# ---------------------------------------------------------------------------
# Workspace guards (safety rules 10, 11 and 12)
# ---------------------------------------------------------------------------

def workspace_root():
    """Get the workspace root that outputs should stay inside, if one is set.

    Returns
    -------
    str or None
        real path of the workspace root, or None when nothing is configured

    Notes
    ----
    The value is taken from $NEOHGT_ROOT when set, else from $NEOHGT_CODE's
    parent, else from the config file's "workspace.code_dir" ancestry. There is
    deliberately **no baked-in default**: the previous constant pointed at one
    developer's machine, so every other machine either failed on its first
    write or, worse, silently wrote into a same-named directory. When nothing
    is configured the containment check is skipped (guard_write still enforces
    the read-only-input rule) and a NOTE is printed once.
    """
    env = os.environ.get('NEOHGT_ROOT')
    if env:
        return realpath(expand(env))
    code = os.environ.get('NEOHGT_CODE')
    if code:
        return dirname(realpath(expand(code)))
    cfg = load_config_file()
    ws = ((cfg or {}).get('workspace') or {})
    if ws.get('code_dir'):
        return dirname(realpath(expand(ws['code_dir'])))
    return None


def load_config_file():
    """The shipped config.yml, or None when it cannot be read.

    Kept separate so workspace_root() has one decision path per source, and so
    a test can substitute it without touching the environment.
    """
    import yaml
    path = find_config_file()
    if not path:
        return None
    try:
        with open(path, encoding='utf-8') as fh:
            return yaml.safe_load(fh)
    except OSError:
        return None


def expand(path):
    """Expand '~' and environment variables in a path (boundary handling)."""
    return os.path.expandvars(os.path.expanduser(str(path)))


def is_within(path, parent):
    """Check if a real path is inside another real path.

    Parameters
    ----------
    path : str
        candidate path (any form)
    parent : str
        container path

    Returns
    -------
    bool
        whether path resolves to a location inside parent
    """
    path, parent = realpath(path), realpath(parent)
    return path == parent or path.startswith(parent + os.sep)


def guard_write(path, root=None, must_not_contain=None):
    """Validate that an output path is allowed, and return it.

    Parameters
    ----------
    path : str
        intended output path
    root : str, optional
        workspace root (default: workspace_root()); when neither is set the
        containment check is skipped and only the read-only input rule applies
    must_not_contain : str, optional
        a read-only input directory that must never be written to

    Returns
    -------
    str
        the real, absolute path

    Raises
    ------
    GuardError
        If no path was given, if the path escapes the workspace, is a symlink
        pointing outside it, or lands inside the read-only input dataset.
    """
    raw = expand(path)
    if path is None or not str(raw).strip() or str(raw).strip() == 'None':
        # a missing -o used to reach this point as the literal string "None"
        # and created a ./None/ tree in the current directory, exit code 0
        raise GuardError(
            'No output path given: pass -o/--output explicitly (got '
            f'{path!r}).')
    if root is None:
        root = workspace_root()
        if root is None:
            print('  NOTE: no workspace root is configured (set NEOHGT_ROOT '
                  'or source tools/env.sh); writes are not contained to a '
                  'workspace, but the read-only input dataset still is.')
    root = realpath(root)
    real = realpath(raw)
    # A symlinked output must not smuggle writes outside of the workspace.
    if not is_within(real, root):
        raise GuardError(
            f'Refusing to write "{raw}" -> "{real}": outside workspace '
            f'"{root}". To write elsewhere, set NEOHGT_ROOT (or pass '
            '--workspace-out/--workspace-tmp on the database stage); writing '
            'into the read-only input dataset stays forbidden either way.')
    if must_not_contain:
        inp = realpath(expand(must_not_contain))
        if is_within(real, inp):
            raise GuardError(
                f'Refusing to write "{real}": the GTDB input dataset at '
                f'"{inp}" is read-only.')
    return real


def assert_input_read_only(gtdb_dir):
    """Warn (never fail) when the input dataset is technically writable.

    Parameters
    ----------
    gtdb_dir : str
        GTDB dataset root

    Returns
    -------
    bool
        whether the directory is writable by the current user
    """
    writable = os.access(realpath(expand(gtdb_dir)), os.W_OK)
    if writable:
        print('WARNING: input dataset is writable by this user; neoHGT will '
              f'never write to it: {gtdb_dir}')
    return writable


# ---------------------------------------------------------------------------
# Layout probing (plan §3.2/§3.3, P0-S0)
# ---------------------------------------------------------------------------

def _first_existing(cands):
    """Return the first candidate path that exists (file or directory)."""
    for c in cands:
        if c and exists(c):
            return c
    return None


def _list_faa(dir_, exts=('.faa', '.faa.gz')):
    """List genome FASTA files in a directory (read-only)."""
    if not isdir(dir_):
        return []
    return sorted(f for f in listdir(dir_) if f.endswith(exts))


def find_input_files(gtdb_dir):
    """Locate the taxonomy, metadata, sequence and list files of a dataset.

    Parameters
    ----------
    gtdb_dir : str
        GTDB dataset root (read-only)

    Returns
    -------
    dict
        discovered paths, keyed by role
    """
    gtdb_dir = realpath(expand(gtdb_dir))
    res = {'gtdb_dir': gtdb_dir}

    res['taxonomy'] = [p for p in (
        _first_existing([join(gtdb_dir, 'tree_bac120', f) for f in
                         _safe_list(join(gtdb_dir, 'tree_bac120'))
                         if f.endswith(('.tsv', '.tsv.gz'))
                         and 'taxonomy' in f]),
        _first_existing([join(gtdb_dir, 'tree_ar53', f) for f in
                         _safe_list(join(gtdb_dir, 'tree_ar53'))
                         if f.endswith(('.tsv', '.tsv.gz'))
                         and 'taxonomy' in f])) if p] or (
                             _root_table_files(gtdb_dir, 'taxonomy'))
    res['metadata'] = [p for p in (
        _first_existing([join(gtdb_dir, 'metadata_bac120', f) for f in
                         _safe_list(join(gtdb_dir, 'metadata_bac120'))
                         if f.endswith(('.tsv', '.tsv.gz'))
                         and 'metadata' in f and 'Sample' not in f]),
        _first_existing([join(gtdb_dir, 'metadata_ar53', f) for f in
                         _safe_list(join(gtdb_dir, 'metadata_ar53'))
                         if f.endswith(('.tsv', '.tsv.gz'))
                         and 'metadata' in f and 'Sample' not in f])) if p] or (
                             _root_table_files(gtdb_dir, 'metadata'))

    # sequence roots, in order of preference: an '*_aa' directory (the local
    # r220 copy), a bare 'aa'/'genomic_files_reps' directory, or the root
    aa_cands = [join(gtdb_dir, d) for d in _safe_list(gtdb_dir)
                if d.endswith('_aa') or d in ('aa', 'nt')]
    # one level deeper, which is how the official tarballs are unpacked
    for d in _safe_list(gtdb_dir):
        p = join(gtdb_dir, d)
        if isdir(p):
            for e in _safe_list(p):
                if e.endswith('_aa') or e in ('aa', 'genomic_files_reps'):
                    aa_cands.append(join(p, e))
    aa = _first_existing(aa_cands) or gtdb_dir
    reps = _first_existing([join(aa, 'protein_faa_reps'),
                            join(gtdb_dir, 'protein_faa_reps'),
                            join(aa, 'genomic_files_reps'), aa])
    sample = _first_existing([join(aa, 'protein_faa_reps_sample'),
                              join(gtdb_dir, 'protein_faa_reps_sample')])
    res['reps_dir'] = reps
    res['sample_dir'] = sample
    res['aa_dir'] = aa if isdir(aa) else gtdb_dir

    domains = {}
    if reps:
        for d in ('bacteria', 'archaea'):
            if isdir(join(reps, d)):
                files = _list_faa(join(reps, d))
                domains[d] = {'dir': join(reps, d), 'n_files': len(files),
                              'ext': '.faa.gz' if (
                                  files and files[0].endswith('.gz'))
                                  else '.faa',
                              'first': files[0] if files else None}
    res['domains'] = domains

    # per-rank sample lists (plan §3.3e), keyed '<domain>_<rank>' so that the
    # bacteria and archaea lists of the same rank do not overwrite each other
    idlists, samples = {}, {}
    for md in res['metadata']:
        dom = 'bacteria' if 'bac120' in md else 'archaea'
        for f in _safe_list(dirname(md)):
            if '_Sample_' in f and f.endswith('.IDlist'):
                rank = f.split('_Sample_')[1][:-len('.IDlist')].lower()
                idlists[f'{dom}_{rank}'] = join(dirname(md), f)
    if sample:
        for d in _safe_list(sample):
            m = re.match(r'^(bacteria|archaea)_(phylum|class|order|family|'
                         r'genus|species)$', d)
            if m:
                samples[f'{m.group(1)}_{m.group(2)}'] = {
                    'dir': join(sample, d),
                    'n_files': len(_list_faa(join(sample, d)))}
    res['idlists'] = idlists
    res['sample_sets'] = samples

    # provenance clues (plan §3.3f)
    prov = []
    for cand in (join(reps or '', 'script'), join(res['aa_dir'],
                                                  'script_decompress')):
        if cand and isdir(cand):
            for f in sorted(listdir(cand)):
                if any(p in f for p in PROVENANCE_PATTERNS):
                    prov.append(join(cand, f))
    res['provenance_scripts'] = prov

    # marker-gene / msa directories are not used by neoHGT, but recording their
    # absence documents the dataset scope
    res['absent_official'] = [name for name in (
        'MD5SUM.txt', 'VERSION.txt', 'FILE_DESCRIPTIONS.txt', 'sp_clusters.tsv',
        'gtdbtk210_taxonomy_file.stable.tsv', 'gtdb.dic')
        if not _first_existing([join(gtdb_dir, name),
                                join(gtdb_dir, 'auxillary_files', name)])]
    return res


def _safe_list(dir_):
    """listdir that tolerates missing directories."""
    try:
        return sorted(listdir(dir_))
    except OSError:
        return []


def _root_table_files(gtdb_dir, kind):
    """Find release tables that sit directly in the dataset root.

    An unpacked official GTDB release is flat: 'bac120_taxonomy_r220.tsv.gz'
    next to 'bac120_metadata_r220.tsv.gz' and the sequence directories, with no
    'tree_bac120/' wrapper. Discovery looked only inside those wrappers, so the
    official layout -- which is what the fetch routes deliver -- was invisible
    unless the user handed over every file through --gtdb-taxonomy.

    Parameters
    ----------
    gtdb_dir : str
        dataset root
    kind : str
        'taxonomy' or 'metadata'

    Returns
    -------
    list of str
        one path per family (bac120, ar53, ...), preferring the plain file
        over its compressed twin
    """
    pat = re.compile(rf'^(?:bac|ar)\d+_{kind}_r\d+(?:\.tsv|\.tsv\.gz)$')
    by_fam = {}
    for f in _safe_list(gtdb_dir):
        if not pat.match(f):
            continue
        fam = f.split('_')[0]
        prev = by_fam.get(fam)
        if prev is None or (prev.endswith('.gz') and not f.endswith('.gz')):
            by_fam[fam] = f
    return [join(gtdb_dir, by_fam[k]) for k in sorted(by_fam)]


def open_reader(path):
    """Open a plain or compressed text file for reading.

    Parameters
    ----------
    path : str
        file path ('.gz', '.bz2', '.xz' or plain)

    Returns
    -------
    file object
        text-mode reader
    """
    ext = '.gz' if path.endswith('.gz') else (
        '.bz2' if path.endswith('.bz2') else (
            '.xz' if path.endswith(('.xz', '.lz')) else ''))
    if not ext:
        return open(path, 'rt', encoding='utf-8', errors='replace')
    mod = {'.gz': gzip, '.bz2': bz2, '.xz': lzma, '.lz': lzma}[ext]
    return mod.open(path, 'rt', encoding='utf-8', errors='replace')


def apply_table_overrides(info, metadata=None, taxonomy=None):
    """Honour explicit --gtdb-metadata/--gtdb-taxonomy paths (plan §7.1).

    Parameters
    ----------
    info : dict
        layout information from detect_layout()
    metadata, taxonomy : iterable of str, optional
        files to use instead of the discovered ones

    Returns
    -------
    dict
        info with the requested tables substituted

    Raises
    ------
    ValueError
        When a file is missing, or when the two sides describe different
        releases or different taxonomic families. The row-level agreement is
        still proven by selfcheck(), which runs right after this.
    """
    if not metadata and not taxonomy:
        return info

    def norm(values, what):
        out = []
        for path in values:
            p = realpath(expand(str(path).strip()))
            if not isfile(p):
                raise ValueError(
                    f'--gtdb-{what} points at a file that does not exist: '
                    f'{path}')
            out.append(p)
        return sorted(out)

    def tokens(files):
        """(release token, set of family tokens) of a GTDB table file set."""
        rels, fams = set(), set()
        def rel(name):
            # '_r220' in the file name, not the 'r53' inside 'ar53_'
            m = re.search(r'_r(\d+)', name)
            return m.group(1) if m else None

        for p in files:
            name = basename(p)
            r = rel(name)
            if r:
                rels.add(r)
            f = re.match(r'((?:bac|ar)\d+)', name)
            if f:
                fams.add(f.group(1))
        return rels, fams

    new = dict(info)
    if metadata:
        new['metadata'] = norm(metadata, 'metadata')
    if taxonomy:
        new['taxonomy'] = norm(taxonomy, 'taxonomy')
    rel_m, fam_m = tokens(new['metadata'])
    rel_t, fam_t = tokens(new['taxonomy'])
    if rel_m and rel_t and rel_m != rel_t:
        raise ValueError(
            f'--gtdb-metadata (release {sorted(rel_m)}) and --gtdb-taxonomy '
            f'(release {sorted(rel_t)}) describe different GTDB releases.')
    if fam_m and fam_t and fam_m != fam_t:
        raise ValueError(
            f'--gtdb-metadata covers {sorted(fam_m)} but --gtdb-taxonomy '
            f'covers {sorted(fam_t)}; the two sides must describe the same '
            'families.')
    return new


def empty_genome_table():
    """The genome-table shape for a lineage-only run (no metadata supplied).

    A taxdump is defined by lineage strings, so the conversion does not need a
    metadata table at all; the derived genome maps simply come out empty.
    """
    return pd.DataFrame(columns=('accession', 'genome', 'domain',
                                'gtdb_taxonomy', 'gtdb_species', 'is_rep',
                                'proteins'))


def metadata_dependent_checks(info):
    """Names of the §3.6 checks a lineage-only dataset cannot answer.

    Built the way `selfcheck()` names them, so a metadata-less run reports
    precisely what it skipped. A hand-copied list here is not an option: the
    first version invented three check names that never existed and printed a
    wrong count.
    """
    names = [f'faa_vs_metadata_reps[{dom}]'
             for dom in ('bacteria', 'archaea')]
    names += [f'taxonomy_vs_metadata[{fam}]' for fam in ('bac120', 'ar53')]
    names += ['metadata_accessions_in_taxonomy',
              'rep_genomes_have_local_fasta']
    return names


def selfcheck_lineages(info, strict=True):
    """The consistency checks that a metadata-less dataset can still pass.

    Only the lineage files themselves are verifiable without a metadata table,
    so this runs the two checks that do not need one and says explicitly which
    §3.6 items were skipped. The rank-depth check is what rejects a lineage
    file whose strings were truncated or rejoined by a third-party script.
    """
    checks = []
    n_rows = 0
    depth_ok = True
    err = ''
    try:
        tax = read_taxonomy(info['taxonomy'])
        n_rows = len(tax)
    except ValueError as exc:
        depth_ok, err = False, str(exc)
    checks.append({'name': 'taxonomy_files_readable',
                   'expected': 'every lineage row has 7 prefixed ranks',
                   'observed': f'{n_rows:,} genome rows' if depth_ok else err,
                   'status': 'OK' if depth_ok else 'FAIL', 'detail': err})
    skipped = metadata_dependent_checks(info)
    report = {'status': 'OK' if depth_ok else 'FAIL', 'checks': checks,
              'n_checks': len(checks), 'n_failed': 0 if depth_ok else 1,
              'n_skipped': len(skipped), 'skipped': skipped}
    print('Self-consistency checks (lineage-only dataset):')
    for c in checks:
        print(f'  [{"ok" if c["status"] == "OK" else "FAIL"}] '
              f'{c["name"]}: {c["observed"]}')
    print(f'  [skipped] {len(skipped)} of the plan §3.6 items need a metadata '
          'table: ' + ', '.join(skipped))
    if not depth_ok and strict:
        raise SelfCheckError(err)
    return report


def detect_layout(gtdb_dir, release=None, layout='auto', metadata=None,
                  taxonomy=None, require_metadata=True,
                  require_sequences=True):
    """Probe a GTDB dataset and classify its layout.

    Parameters
    ----------
    gtdb_dir : str
        dataset root (read-only)
    release : str, optional
        expected release number; inferred from file names when omitted
    layout : str, optional
        force a layout instead of detecting it ('auto' by default)
    metadata, taxonomy : list of str, optional
        explicit table files, i.e. `--gtdb-metadata` / `--gtdb-taxonomy`. They
        are applied *before* the presence checks: previously a dataset whose
        files did not sit in the discovered directories raised here, so these
        options could never take effect in exactly the case they exist for.
    require_metadata, require_sequences : bool, optional
        a database build needs both; converting lineages into a taxdump needs
        neither. Defaults keep the historical behaviour.

    Returns
    -------
    dict
        layout information: label, release, discovered files, id convention,
        sequence counts and coverage statistics
    """
    info = apply_table_overrides(find_input_files(gtdb_dir), metadata,
                                 taxonomy)
    if not info['taxonomy']:
        raise ValueError(
            f'No GTDB taxonomy files found under {gtdb_dir}. Discovery looks '
            "for 'tree_bac120/*_taxonomy*.tsv' and 'tree_ar53/*_taxonomy*.tsv'"
            "; point at the files with --gtdb-taxonomy if they live elsewhere, "
            'as they do in an unpacked official release.')
    if not info['metadata'] and require_metadata:
        raise ValueError(
            f'No GTDB metadata files found under {gtdb_dir}. Discovery looks '
            "for 'metadata_bac120/*_metadata*.tsv' and "
            "'metadata_ar53/*_metadata*.tsv'; pass them with --gtdb-metadata. "
            'Only the taxdump conversion can run without a metadata table.')

    # release: from file names such as bac120_metadata_r220.tsv
    if release is None:
        for p in info['metadata'] + info['taxonomy']:
            m = re.search(r'r(\d{3})', basename(p))
            if m:
                release = m.group(1)
                break
    info['release'] = release or 'unknown'

    # sequence presence and id convention
    exts = {d: v['ext'] for d, v in info['domains'].items()}
    plain = all(e == '.faa' for e in exts.values())
    gzipped = all(e == '.faa.gz' for e in exts.values())

    anchor_hits = 0
    total_ids = 0
    max_id_len = 0
    for d, v in info['domains'].items():
        if not v['first']:
            continue
        for stem in (v['first'],):
            with open_reader(join(v['dir'], stem)) as f:
                n = 0
                for line in f:
                    if not line.startswith('>'):
                        continue
                    sid = line[1:].split('\t')[0].split(' ')[0]
                    total_ids += 1
                    max_id_len = max(max_id_len, len(line.rstrip('\r\n')))
                    if SERIAL_ANCHOR in sid:
                        anchor_hits += 1
                    n += 1
                    if n >= 25:
                        break
    coverage = anchor_hits / total_ids if total_ids else 0.0
    info['serial_anchor_coverage'] = round(coverage, 4)
    info['max_first_line_len'] = max_id_len

    if layout == 'auto':
        if info['domains'] and plain and coverage >= 0.9:
            # renamed local copy: plain .faa, lineage in file name, ids
            # rewritten. The label carries the detected release: the previous
            # three-way labelled every release "gtdb-r220-renamed", which made
            # the tag (and the provenance lines derived from it) meaningless.
            label = f'gtdb-r{info["release"] or "unknown"}-renamed'
        elif info['domains'] and gzipped:
            label = 'official-gz'
        elif info['domains'] and plain:
            label = 'official-plain'
        elif not info['domains']:
            if not require_sequences:
                # converting lineages into a taxdump is defined by the lineage
                # files alone; the per-genome FASTA tree is irrelevant here
                label = ('official-gz' if any(p.endswith('.gz')
                                              for p in info['taxonomy'])
                         else 'official-plain')
                info['lineage_only'] = True
            else:
                raise ValueError(
                    f'Cannot detect a per-genome protein FASTA directory '
                    f'under {gtdb_dir}; checked '
                    f'{join(info["aa_dir"], "protein_faa_reps")}. Candidate '
                    'layouts: ' + ', '.join(LOCAL_LAYOUTS) + '.')
        else:
            raise ValueError(
                f'Undecided layout under {gtdb_dir}. Candidate layouts: '
                + ', '.join(LOCAL_LAYOUTS) + '.')
    else:
        if layout not in LOCAL_LAYOUTS + (LEGACY_RENAMED_LAYOUT,):
            raise ValueError(
                f'Unknown --local-layout "{layout}". Choices: '
                + ', '.join(LOCAL_LAYOUTS) + '.')
        label = layout

    info['layout'] = label
    info['id_convention'] = (
        'stem_SerialNumber_contig.ver_locus' if is_renamed_layout(label)
        else 'contig.ver_locus')
    info['n_genome_files'] = sum(v['n_files'] for v in info['domains'].values())
    return info


# ---------------------------------------------------------------------------
# Lineage and identifier parsing
# ---------------------------------------------------------------------------

def split_lineage(lineage):
    """Split a GTDB lineage string into its seven prefixed segments.

    Parameters
    ----------
    lineage : str
        lineage such as 'd__Bacteria;p__...;s__Escherichia coli'

    Returns
    -------
    list of tuple
        (rank, label) pairs, high to low

    Raises
    ------
    ValueError
        If the lineage is incomplete or uses an unknown rank prefix.
    """
    segs = [s for s in str(lineage).split(LINEAGE_SEP)]
    if len(segs) != DEPTH:
        raise ValueError(
            f'GTDB lineage must have {DEPTH} ranks, got {len(segs)}: {lineage}')
    out = []
    for rank, seg in zip(GTDB_RANKS, segs):
        char = seg[:1].lower()
        if PREFIX_BY_CHAR.get(char) != rank or not seg.startswith(char):
            raise ValueError(
                f'Segment "{seg}" does not carry the expected rank prefix '
                f'"{GTDB_RANK_PREFIXES[rank]}" for rank "{rank}".')
        out.append((rank, seg))
    return out


def lineage_paths(lineage):
    """Get all ancestral lineage paths of a lineage string.

    Parameters
    ----------
    lineage : str
        full GTDB lineage

    Returns
    -------
    list of str
        paths from domain to species, each a ';'-joined prefix
    """
    segs = [s for s in str(lineage).split(LINEAGE_SEP)]
    return [LINEAGE_SEP.join(segs[:i + 1]) for i in range(len(segs))]


def parent_path(path):
    """Get the parent lineage path of a node key ('' for the root level)."""
    i = path.rfind(LINEAGE_SEP)
    return path[:i] if i > 0 else ''


def path_rank(path):
    """Get the rank of a lineage path (its deepest segment)."""
    seg = path.split(LINEAGE_SEP)[-1]
    rank = rank_of_label(seg)
    if rank is None:
        raise ValueError(f'Cannot infer rank of lineage segment "{seg}".')
    return rank


def path_label(path):
    """Get the prefixed label of a lineage path."""
    return path.split(LINEAGE_SEP)[-1]


def strip_faa_ext(name):
    """Strip a FASTA file extension (one helper for all callers).

    ``.faa.gz`` must be tested before ``.faa``: stripping the inner extension
    first (the previous per-callsite behaviour) left a stale ``.faa`` on
    ``*.faa.gz`` names and mis-parsed the trailing lineage anchor.

    Parameters
    ----------
    name : str
        file name or stem, with or without extension

    Returns
    -------
    str
        name without ``.faa.gz`` / ``.faa`` / ``.gz`` suffix
    """
    for ext in ('.faa.gz', '.faa', '.gz'):
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


def parse_stem(stem):
    """Parse a renamed local FASTA file name into its lineage segments.

    Parameters
    ----------
    stem : str
        file name without extension, e.g.
        'GB_GCA_000008885.1_d_Bacteria_p_Pseudomonadota_c_Gammaproteobacteria'
        '_o_Enterobacterales_f_Enterobacteriaceae_A_g_Wigglesworthia'

    Returns
    -------
    dict
        'accession' plus one key per rank ('domain'..'genus'), each holding the
        raw segment including any cluster letter suffix

    Notes
    ----
    Local file names stop at genus level (plan §3.3c), so no species is
    invented; only rank anchors are used, because names may themselves contain
    hyphens and digits (e.g. 'o_ML615J-28', 'g_336-3').
    """
    stem = strip_faa_ext(stem)
    matches = list(STEM_ANCHOR_RE.finditer(stem))
    res = OrderedDict()
    expected = ['d', 'p', 'c', 'o', 'f', 'g']   # domain..genus, no species
    chain, acc = [], None
    for m in matches:
        want = expected[len(chain)] if len(chain) < len(expected) else None
        if m.group(1) == want:
            if acc is None:
                acc = stem[:m.start()]
            elif chain:
                res[res_key(chain[-1])] = stem[res_end:m.start()]
            chain.append(m.group(1))
            res_end = m.end()
    if len(chain) == len(expected):
        res[res_key(chain[-1])] = stem[res_end:]
        res['accession'] = acc
        return res
    raise ValueError(
        f'File name stem does not follow the local rank-anchor convention: '
        f'{stem}')


def res_key(char):
    """Map a rank anchor character to its rank name."""
    return PREFIX_BY_CHAR[char]


def parse_seqid(seqid, layout='local-renamed'):
    """Extract the genome identifier from a protein sequence Id.

    Parameters
    ----------
    seqid : str
        protein sequence Id
    layout : str, optional
        layout label, which decides the parsing anchor

    Returns
    -------
    dict or None
        'genome' (GTDB or nucleotide accession) plus 'stem' when available
    """
    seqid = str(seqid).strip().lstrip('>').split('\t')[0].split(' ')[0]
    if SERIAL_ANCHOR in seqid:
        stem = seqid.split(SERIAL_ANCHOR, 1)[0]
        try:
            parsed = parse_stem(stem)
        except ValueError:
            return {'stem': stem, 'genome': None}
        return {'stem': stem, 'genome': parsed['accession']}
    # official GTDB convention: <nucleotide acc.ver>_<locus number>
    m = re.match(r'^(?P<contig>.+)\.(?P<ver>\d+)_(?P<locus>\d+)$', seqid)
    if m:
        return {'contig': m.group('contig') + '.' + m.group('ver'),
                'contig_nover': m.group('contig'), 'genome': None,
                'locus': m.group('locus')}
    m = re.match(r'^(?P<contig>.+)_(?P<locus>\d+)$', seqid)
    if m:
        return {'contig': m.group('contig'), 'genome': None,
                'locus': m.group('locus')}
    return None


# ---------------------------------------------------------------------------
# Metadata and taxonomy loading
# ---------------------------------------------------------------------------

def read_taxonomy(files):
    """Read GTDB lineage files.

    Parameters
    ----------
    files : list of str
        paths of '*_taxonomy*.tsv' files (2 columns, no header)

    Returns
    -------
    dict
        accession -> lineage string
    """
    res = {}
    dup = 0
    for path in files:
        with open_reader(path) as f:
            for ln, line in enumerate(f, 1):
                line = line.rstrip('\r\n')
                if not line:
                    continue
                x = line.split('\t')
                if len(x) < 2:
                    raise ValueError(
                        f'{path}:{ln}: expected 2 tab-delimited fields, got '
                        f'{len(x)}.')
                acc, lin = x[0], x[1]
                if lin.count(LINEAGE_SEP) != DEPTH - 1:
                    raise ValueError(
                        f'{path}:{ln}: lineage of {acc} does not have {DEPTH} '
                        'ranks (provenance check failed).')
                if acc in res:
                    dup += 1
                    if res[acc] != lin:
                        raise ValueError(
                            f'{path}:{ln}: accession {acc} has two different '
                            'lineages; the dataset is inconsistent.')
                res[acc] = lin
    if dup:
        print(f'WARNING: {dup} duplicated accession rows in taxonomy files.')
    return res


def load_metadata(files, release=None):
    """Load GTDB metadata tables into a normalized frame.

    Parameters
    ----------
    files : list of str
        paths of '*_metadata*.tsv' files
    release : str, optional
        release label recorded on every row

    Returns
    -------
    pandas.DataFrame
        one row per genome, using neoHGT internal column names

    Raises
    ------
    ValueError
        If a required column is missing (columns are matched by name only).

    Notes
    ----
    Results are cached per (file set, release) because a build reads the same
    780 MB table from several entry points; the frame is never mutated in
    place by callers, who copy it before selecting.
    """
    key = (tuple(sorted(realpath(expand(p)) for p in files)), release)
    if key in _METADATA_CACHE:
        return _METADATA_CACHE[key].copy()
    frames = []
    for path in files:
        header = next(open_reader(path)).rstrip('\r\n').split('\t')
        miss = [c for c in REQUIRED_COLS if c not in header]
        if miss:
            raise ValueError(
                f'{basename(path)}: missing required metadata column(s): '
                + ', '.join(miss) + '. Columns are resolved by name, so this '
                'release is not supported by the current column map.')
        usecols = [c for c in tuple(REQUIRED_COLS) + tuple(OPTIONAL_COLS)
                   if c in header]
        df = pd.read_csv(path, sep='\t', usecols=usecols, dtype=str,
                         keep_default_na=False, na_filter=False)
        df['_source_file'] = realpath(path)
        frames.append(df)
        print(f'  loaded {len(df):,} rows from {basename(path)} '
              f'({len(usecols)} of {len(header)} columns)')
    df = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]

    # normalize into internal columns (plan §5.4 table)
    dom = df['gtdb_taxonomy'].str.split(LINEAGE_SEP).str[0]
    df['domain'] = dom.map({'d__Bacteria': 'bacteria',
                            'd__Archaea': 'archaea'}).fillna('other')
    df['accession'] = df['accession'].astype(str)
    df['accnov'] = df['accession'].str.replace(r'\.\d+$', '', regex=True)
    # NB: str.lstrip would strip the character set {R,S,G,B,_} and destroy
    # 'GCA_'/'GCF_'; an anchored regex is required here (plan §5.4).
    df['genome'] = df['accnov'].str.replace(r'^(RS|GB)_', '', regex=True)
    df['genome'] = df['domain'].map({'bacteria': 'BAC', 'archaea': 'ARC'}) \
        .fillna('OTH') + '_' + df['genome']
    df['gtdb_species'] = df['gtdb_taxonomy'].str.split(
        LINEAGE_SEP).str[DEPTH - 1]
    df['organism_name'] = df['gtdb_species'].map(strip_rank_prefix)
    df['species'] = df['gtdb_species']
    df['is_rep'] = df['gtdb_representative'].astype(str).str.strip() == 't'
    df['release'] = release or 'unknown'
    df['proteins'] = pd.to_numeric(df.get('protein_count', 0), errors='coerce')\
        .fillna(0).astype('int64')
    df['completeness'] = pd.to_numeric(
        df.get('checkm2_completeness', 'nan'), errors='coerce')
    df['contamination'] = pd.to_numeric(
        df.get('checkm2_contamination', 'nan'), errors='coerce')
    df['quality_score'] = df['completeness'] - 5 * df['contamination']
    df['assembly_level'] = df.get('ncbi_assembly_level', '')
    df['ncbi_taxid'] = df.get('ncbi_taxid', '')
    df['category'] = df['gtdb_representative']
    df['tax_source'] = 'gtdb'
    if 'gtdb_type_species_of_genus' in df:
        df['type_material'] = (df['gtdb_type_species_of_genus'] == 't')
    else:
        df['type_material'] = False
    _METADATA_CACHE[key] = df
    return df.copy()


# ---------------------------------------------------------------------------
# Stable synthetic taxIds and taxdump conversion (plan §5.3)
# ---------------------------------------------------------------------------

def collect_nodes(lineages, keep=None):
    """Expand lineage strings into a unique node set.

    Parameters
    lineages : iterable of str
        full GTDB lineage strings
    keep : set of str, optional
        restrict the tree to nodes whose lineage path is a prefix of any
        accession lineage in this set (taxdump trimming)

    Returns
    -------
    dict
        node path -> {'rank', 'label', 'name', 'cluster_letter', 'parent_key'}
    """
    nodes = {}
    for lin in lineages:
        for path in lineage_paths(lin):
            if path in nodes:
                continue
            label = path_label(path)
            rank = path_rank(path)
            nodes[path] = {
                'rank': rank,
                'label': label,
                'name': strip_rank_prefix(label),
                'cluster_letter': _cluster_letter(label),
                'parent_key': parent_path(path)}
    if keep is not None:
        nodes = {k: v for k, v in nodes.items() if k in keep}
    return nodes


def _cluster_letter(label):
    """Get the species-cluster letter of a label, if any."""
    bare = strip_rank_prefix(label, strip_cluster=False)
    stripped = strip_rank_prefix(label)
    if bare == stripped:
        return ''
    m = re.search(r'_([A-Z])(\s.*|)$', bare)
    return m.group(1) if m else ''


def reserved_ids_path(manifest_path):
    """Sidecar next to taxid.map.gz listing taxIds no longer in use but never reusable.

    A taxId that some earlier release assigned stays reserved forever: if a
    later release handed it to a different lineage, historical `scores.tsv`
    and `hgts/*.txt` files would silently mean another organism (plan RK6).
    """
    if not manifest_path:
        return None
    return join(dirname(manifest_path), 'reserved_ids.txt')


def load_reserved(manifest_path):
    """Read the reserved taxId sidecar (empty when absent).

    Returns
    -------
    dict
        taxId -> node path it used to belong to (the path is audit
        information and used to be dropped, which erased the "who owned this
        id" trail at the second rebuild)
    """
    path = reserved_ids_path(manifest_path)
    out = {}
    if not path or not exists(path):
        if manifest_path and exists(manifest_path):
            print(f'  NOTE: no reserved-ids sidecar next to '
                  f'"{manifest_path}": historical taxIds from earlier builds '
                  f'cannot be carried forward this round.')
        return out
    for line in open(path):
        if line.startswith('#') or line.startswith('taxid\t'):
            continue
        x = line.rstrip('\r\n').split('\t')
        if x and x[0].isdigit():
            out[x[0]] = x[1] if len(x) > 1 else ''
    return out


def write_reserved(path, ids, info=None):
    """Write the reserved taxId sidecar (sorted, one id per line)."""
    path = guard_write(path)
    with open(path, 'w') as f:
        f.write(provenance(info, extra={'count': len(ids)}))
        f.write('# taxIds used by an earlier release; never reassigned '
                '(plan §5.3.3 rule 2 / RK6)\n')
        f.write('taxid\tnode_path\n')
        for tid, key in sorted(ids, key=lambda kv: int(kv[0])):
            f.write(f'{tid}\t{key}\n')
    return path


def load_manifest(path):
    """Read a previous taxid manifest for ID stability.

    Parameters
    ----------
    path : str
        path of 'gtdb/taxid.map.gz' (missing file yields an empty manifest)

    Returns
    -------
    dict
        node path -> {'taxid', 'rank', 'label', 'name', 'ncbi_taxid'}
    """
    man = {}
    if not path or not exists(path):
        return man
    with open_reader(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            x = line.rstrip('\r\n').split('\t')
            if len(x) < 2 or x[0] == 'lineage':
                continue
            man[x[0]] = {'taxid': x[1], 'rank': x[2] if len(x) > 2 else '',
                         'label': x[3] if len(x) > 3 else '',
                         'name': x[4] if len(x) > 4 else '',
                         'ncbi_taxid': x[7] if len(x) > 7 else ''}
    return man


def assign_ids(nodes, manifest=None, id_base=DEFAULT_ID_BASE, reserve=()):
    """Assign stable numeric taxIds to GTDB nodes.

    Parameters
    ----------
    nodes : dict
        node path -> node record (from collect_nodes)
    manifest : dict, optional
        previous taxid.map contents; matching nodes keep their old Id
    id_base : int, optional
        start of the reserved Id range
    reserve : iterable of str, optional
        taxIds that were used by an older release and must never be handed out
        again, even though their node is gone (plan RK6)

    Returns
    -------
    tuple of dict
        (node path -> taxId, node path -> taxId) for superseded nodes

    Notes
    ----
    Fresh allocation uses rank_ordinal * 1e6 + lexicographic index, which is
    deterministic and readable. When a manifest is supplied, existing Ids are
    reused verbatim and only new nodes are appended from a free cursor, so Ids
    never change between releases (plan §5.3.3, rule 1).
    """
    manifest = manifest or {}
    taxid_map, used = {}, set(reserve)

    for key, rec in manifest.items():
        if not rec['taxid'].isdigit():
            continue
        if key in nodes:
            taxid_map[key] = rec['taxid']
        # ids of nodes that vanished from this release stay reserved: otherwise a
        # renamed taxon would inherit a historical taxId and old scores.tsv
        # would silently mean a different organism (plan RK6)
        used.add(rec['taxid'])

    cursor = id_base
    for rank in GTDB_RANKS:
        keys = [k for k in nodes if nodes[k]['rank'] == rank]
        keys.sort(key=lambda k: (taxid_map.get(parent_path(k), '~'),
                                 nodes[k]['name'], k))
        slot = id_base + RANK_ORDINAL[rank] * ID_SLOT
        for i, key in enumerate(keys):
            if key in taxid_map:
                continue
            tid = slot + i
            if tid - slot >= ID_SLOT or str(tid) in used or tid > MAX_TAXID:
                # rank slot exhausted or occupied: append from the cursor
                while str(cursor) in used or cursor > MAX_TAXID:
                    cursor += 1
                tid = cursor
                cursor += 1
            if tid > MAX_TAXID:
                raise ValueError(
                    f'Cannot assign taxId for "{key}": int32 range exceeded.')
            taxid_map[key] = str(tid)
            used.add(str(tid))

    # bijection checks (UT-1.2)
    if len(set(taxid_map.values())) != len(taxid_map):
        raise AssertionError('taxId assignment is not a bijection.')
    for tid in taxid_map.values():
        if not (id_base <= int(tid) <= MAX_TAXID):
            raise AssertionError(
                f'taxId {tid} is outside the reserved range [{id_base}, '
                f'{MAX_TAXID}].')
    superseded = {k: v['taxid'] for k, v in manifest.items()
                  if k not in nodes}
    return taxid_map, superseded


def build_taxdump_from_gtdb(lineages, manifest=None, id_base=DEFAULT_ID_BASE,
                            root_parent=None, manifest_path=None):
    """Convert GTDB lineages into an NCBI-compatible taxdump.

    Parameters
    ----------
    lineages : iterable of str
        full GTDB lineage strings
    manifest : dict, optional
        previous taxid manifest for ID stability
    id_base : int, optional
        start of the reserved GTDB taxId range
    root_parent : str, optional
        taxId that GTDB domains attach to (default: the 'root' node '1';
        hybrid mode uses the NCBI cellular root '131567')

    Returns
    -------
    tuple
        (taxdump, taxid_map, stats) where taxdump has the same structure as
        util.read_taxdump() output
    """
    root_parent = root_parent or ROOT_TID
    nodes = collect_nodes(lineages)
    taxid_map, superseded = assign_ids(nodes, manifest, id_base,
                                       reserve=load_reserved(manifest_path))

    taxdump = {ROOT_TID: {'parent': ROOT_TID, 'rank': 'no rank',
                          'name': 'root'}}
    for key, rec in nodes.items():
        tid = taxid_map[key]
        parent_key = rec['parent_key']
        parent = taxid_map[parent_key] if parent_key else root_parent
        if parent_key and parent_key not in taxid_map:
            raise AssertionError(
                f'Node "{key}" has no resolvable parent "{parent_key}".')
        taxdump[tid] = {'parent': parent, 'rank': rec['rank'],
                        'name': rec['name']}

    # integrity assertions (plan 附录 A.2 step 7)
    for tid, t in taxdump.items():
        if t['parent'] not in taxdump:
            raise AssertionError(
                f'Node {tid} references missing parent {t["parent"]}.')
    if find_lca(list(taxdump), taxdump) != ROOT_TID:
        raise AssertionError(
            'GTDB conversion failed: the tree is not single-rooted, so '
            'find_lca() would be unusable.')
    stats = {
        'n_nodes': len(taxdump),
        'n_gtdb_nodes': len(nodes),
        'by_rank': {r: sum(1 for v in nodes.values() if v['rank'] == r)
                    for r in GTDB_RANKS},
        'n_superseded': len(superseded),
        'id_base': id_base,
        'id_min': min(int(v) for v in taxid_map.values()) if taxid_map else 0,
        'id_max': max(int(v) for v in taxid_map.values()) if taxid_map else 0,
        'tree': nodes,
        'superseded': superseded,
    }
    return taxdump, taxid_map, stats


# Ranks as stored in nodes.dmp. DIAMOND validates the rank column against the
# NCBI vocabulary and aborts with "Invalid taxonomic rank: domain" (spike S1),
# so GTDB's 'domain' is written as NCBI's 'superkingdom' while the official
# name is preserved in gtdb/labels.tsv and gtdb/taxid.map.gz.
DMP_RANK_ALIASES = {'domain': 'superkingdom'}

# The rank vocabulary DIAMOND 2.1.10 actually accepts. Hybrid mode grafts NCBI
# nodes in verbatim, and current NCBI dumps spell taxId 131567 'cellular root'
# and use 'clade' above the superkingdoms -- both make makedb abort with
# "Error: Invalid taxonomic rank: cellular root". Anything outside this set is
# written as 'no rank', which is what those same nodes were called before NCBI
# replaced the rank. Only the file is affected: the merged tree, the lineage
# tables and gtdb/labels.tsv keep the true rank, and neoHGT never branches on
# these ranks (it uses taxIds and labels).
DIAMOND_RANKS = frozenset((
    'no rank', 'superkingdom', 'kingdom', 'subkingdom', 'branch', 'fleet',
    'subbranch', 'superphylum', 'phylum', 'subphylum', 'infraphylum',
    'microphylum', 'parvphylum', 'superclass', 'class', 'subclass',
    'infraclass', 'magnus', 'ordo', 'superorder', 'order', 'suborder',
    'infraorder', 'parvorder', 'cfamily', 'superfamily', 'family',
    'subfamily', 'tribe', 'subtribe', 'species_group', 'species_subgroup',
    'genus', 'subgenus', 'section', 'subsection', 'series', 'subseries',
    'species', 'subspecies', 'varietas', 'subvarietas', 'forma', 'subforma',
    'strain', 'serogroup', 'serotype', 'pathogroup', 'biotype', 'clone',
    'biovar', 'chemovar', 'morphovar', 'phagovar', 'serovar', 'vargroup'))
DMP_RANK_FALLBACK = 'no rank'


def dmp_rank(rank):
    """Get the rank name to store in nodes.dmp for a GTDB rank."""
    rank = DMP_RANK_ALIASES.get(rank, rank)
    return rank if rank in DIAMOND_RANKS else DMP_RANK_FALLBACK


def illegal_dmp_ranks(taxdump):
    """Return the ranks of a tree that DIAMOND would reject, with their taxIds.

    Parameters
    ----------
    taxdump : dict
        taxonomy database

    Returns
    -------
    dict
        rank -> list of taxIds, empty if the tree can be written as is

    Notes
    ----
    The check runs after DMP_RANK_ALIASES and before the 'no rank' fallback,
    i.e. it lists exactly the nodes write_taxdump() has to flatten.
    """
    bad = {}
    for tid, rec in taxdump.items():
        rank = rec.get('rank')
        if DMP_RANK_ALIASES.get(rank, rank) not in DIAMOND_RANKS:
            bad.setdefault(rank, []).append(tid)
    return bad


def write_taxdump(taxdump, outdir):
    """Write nodes.dmp and names.dmp in NCBI format.

    Parameters
    ----------
    taxdump : dict
        taxonomy database
    outdir : str
        target directory (created if missing)

    Returns
    -------
    dict
        paths written

    Notes
    ----
    Two details matter for external tools:
    - nodes.dmp carries the full NCBI column layout (13 fields). A 3-field
      line is fine for util.read_taxdump(), but DIAMOND splits on the
      three-character delimiter "\\t|\\t", so a short line makes the rank the
      last field and reads back as "no rank\\t|", which aborts makedb.
    - the rank vocabulary is validated by DIAMOND, hence dmp_rank() below.
    """
    makedirs(outdir, exist_ok=True)
    nodes = join(outdir, 'nodes.dmp')
    names = join(outdir, 'names.dmp')
    order = sorted(taxdump, key=lambda x: (int(x) if x.isdigit() else 0, x))
    with open(nodes, 'w') as f:
        for tid in order:
            t = taxdump[tid]
            f.write(f'{tid}\t|\t{t["parent"]}\t|\t{dmp_rank(t["rank"])}\t|\t'
                    f'\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n')
    with open(names, 'w') as f:
        for tid in order:
            name = taxdump[tid].get('name', '')
            f.write(f'{tid}\t|\t{name}\t|\t\t|\tscientific name\t|\n')
    return {'nodes.dmp': nodes, 'names.dmp': names}


# ---------------------------------------------------------------------------
# Input self-consistency checks (plan §3.6, P0-S0)
# ---------------------------------------------------------------------------

def _count_lines(path):
    n = 0
    with open_reader(path) as f:
        for _ in f:
            n += 1
    return n


def selfcheck(info, strict=True, deep=True):
    """Run the input consistency checks required before building.

    The number of checks is a property of the dataset (one per metadata
    table, taxonomy file and idlist), so no fixed count is promised here; the
    report carries n_checks and callers should print that, not a literal.


    Parameters
    ----------
    info : dict
        layout information from detect_layout()
    strict : bool, optional
        raise SelfCheckError when any check fails
    deep : bool, optional
        also verify per-accession alignment, which re-reads the metadata table

    Returns
    -------
    dict
        report with one entry per check and an overall status
    """
    checks = []
    df = load_metadata(info['metadata'], release=info['release'])
    tax = read_taxonomy(info['taxonomy'])

    def add(name, expected, observed, ok, detail=''):
        checks.append({'name': name, 'expected': expected,
                       'observed': observed,
                       'status': 'OK' if ok else 'FAIL', 'detail': detail})
        print(f'  [{"ok" if ok else "FAIL"}] {name}: expected {expected}, '
              f'observed {observed}{(" - " + detail) if detail else ""}')

    print('Self-consistency checks (plan §3.6):')

    # 1 & 2: per-genome FASTA file counts vs representative counts
    for dom in ('bacteria', 'archaea'):
        got = info['domains'].get(dom, {}).get('n_files', 0)
        want = int((df[(df['domain'] == dom) & df['is_rep']]).shape[0])
        add(f'faa_vs_metadata_reps[{dom}]', want, got, want == got)

    # 3: taxonomy rows vs metadata rows
    by_dom = df['domain'].value_counts().to_dict()
    tax_dom = {'bac120': 'bacteria', 'ar53': 'archaea'}
    for path in info['taxonomy']:
        key = next((k for k in tax_dom if k in basename(path)), None)
        if key is None:
            continue
        want, got = by_dom.get(tax_dom[key], -1), _count_lines(path)
        add(f'taxonomy_vs_metadata[{key}]', want, got, want == got)

    # 4: lineage depth is always complete
    bad = [a for a, lin in tax.items() if lin.count(LINEAGE_SEP) != DEPTH - 1]
    add('lineage_depth_is_7', 0, len(bad), not bad,
        f'{len(tax):,} lineages scanned')

    # 5: rank prefixes follow the closed vocabulary
    wrong = 0
    for lin in tax.values():
        try:
            split_lineage(lin)
        except ValueError:
            wrong += 1
    add('rank_prefix_vocabulary', 0, wrong, wrong == 0)

    # 6: per-rank sample lists match the number of unique taxa
    for name, path in sorted(info['idlists'].items()):
        dom, _, rank = name.partition('_')
        if rank not in RANK_ORDINAL:
            continue
        col = GTDB_RANK_PREFIXES[rank]
        uniq = set()
        for lin in tax.values():
            if not lin.startswith('d__Bacteria' if dom == 'bacteria'
                                  else 'd__Archaea'):
                continue
            segs = lin.split(LINEAGE_SEP)
            uniq.add(LINEAGE_SEP.join(segs[:RANK_ORDINAL[rank] + 1]))
        got = _count_lines(path)
        add(f'idlist_vs_taxa[{name}]', len(uniq), got, got == len(uniq))

    # 7: every domain is one of the two GTDB domains (no eukaryotes/viruses).
    # "Exactly two" used to hard-reject single-domain datasets even though
    # --domain bacteria|archaea is a supported selection.
    doms = {lin.split(LINEAGE_SEP)[0] for lin in tax.values()}
    add('domains_are_bacteria_archaea', 'subset of {d__Bacteria, d__Archaea}',
        len(doms),
        bool(doms) and doms <= {'d__Bacteria', 'd__Archaea'},
        ', '.join(sorted(doms)))

    # 8: accession alignment between the three sources actually used
    if deep:
        md_acc = set(df['accession'])
        missing_seq = df[df['is_rep'] & ~df['accession'].isin(
            _index_accessions(info))]
        add('metadata_accessions_in_taxonomy', len(df), len(md_acc & set(tax)),
            md_acc == set(tax),
            f'{len(md_acc - set(tax)):,} not in taxonomy, '
            f'{len(set(tax) - md_acc):,} in taxonomy but not in metadata')
        add('rep_genomes_have_local_fasta', int(df['is_rep'].sum()),
            int(df['is_rep'].sum()) - len(missing_seq), len(missing_seq) == 0,
            f'{len(missing_seq):,} representative genomes lack a local FASTA')
        if len(missing_seq):
            info['missing_fasta'] = missing_seq['accession'].tolist()[:1000]

    failed = [c for c in checks if c['status'] == 'FAIL']
    res = {'checks': checks, 'n_checks': len(checks),
           'n_failed': len(failed), 'status': 'OK' if not failed else 'FAIL',
           'started': str(timestamp())}
    if failed and strict:
        raise SelfCheckError(
            'Local dataset self-consistency check failed: '
            + '; '.join(f'{c["name"]} (expected {c["expected"]}, observed '
                        f'{c["observed"]})' for c in failed)
            + '. Refusing to build on an inconsistent input (ADR-9).')
    return res


def _index_accessions(info):
    """Collect the GTDB accessions present as per-genome FASTA files."""
    out = set()
    for dom, v in info['domains'].items():
        for fname in _list_faa(v['dir']):
            stem = strip_faa_ext(fname)
            out.add(stem.split('_d_')[0] if '_d_' in stem
                    else stem.replace('_protein', ''))
    return out


def reconcile(gtdb_dir=None, info=None, out=None, release=None,
              layout='auto'):
    """Write the full accession-level difference lists of the three sources.

    Parameters
    ----------
    gtdb_dir : str, optional
        read-only dataset root (ignored when info is given)
    info : dict, optional
        already-detected layout information
    out : str, optional
        destination directory for the TSV lists
    release, layout : str, optional
        passed to detect_layout() when info is not given

    Returns
    -------
    dict
        per-list counts plus the written paths

    Notes
    ----
    Plan DoD-2 requires that "counts must agree one by one and the mismatch
    list must be written to disk". selfcheck only proves the counts, so this
    function emits the actual accession sets that differ; on a consistent copy
    every file holds only its header.
    """
    if info is None:
        if not gtdb_dir:
            raise ValueError('reconcile() needs gtdb_dir or info')
        info = detect_layout(gtdb_dir, release=release, layout=layout)
    df = load_metadata(list(info['metadata']), release=info['release'])
    tax = read_taxonomy(info['taxonomy'])
    local = _index_accessions(info)
    reps = {str(a) for a in df.loc[df['is_rep'], 'accession']}
    allmd = {str(a) for a in df['accession']}
    nonreps = allmd - reps

    lists = {
        'missing_faa.tsv': sorted(reps - local),        # rep w/o local FASTA
        'extra_faa.tsv': sorted(local - allmd),         # FASTA unknown to GTDB
        'nonrep_faa.tsv': sorted(local & nonreps),      # member files present
        'metadata_not_in_taxonomy.tsv': sorted(allmd - set(tax)),
        'taxonomy_not_in_metadata.tsv': sorted(
            {str(k) for k in tax} - allmd),
    }
    summary = {'gtdb_dir': realpath(expand(info['gtdb_dir'])),
               'release': info['release'], 'layout': info['layout'],
               'n_metadata_genomes': len(allmd),
               'n_metadata_reps': len(reps),
               'n_local_faa': len(local),
               'n_taxonomy_rows': len(tax),
               'lists': {k: len(v) for k, v in lists.items()},
               'files': {}}
    if out:
        out = guard_write(out, must_not_contain=info['gtdb_dir'])
        makedirs(out, exist_ok=True)
        for name, items in lists.items():
            path = join(out, name)
            with open(path, 'w') as f:
                f.write(provenance(info, extra={
                    'source': info['gtdb_dir'], 'list': name,
                    'count': len(items)}))
                f.write(f'# {name}: {len(items)} accession(s)\n')
                f.write('accession\n')
                for acc in items:
                    f.write(f'{acc}\n')
            summary['files'][name] = path
        summary['consistent'] = not any(lists.values())
    return summary


def sha256_file(path, block=1 << 20):
    """Hash a file with sha256 (streamed, so 780 MB inputs are safe)."""
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_text(text, block=1 << 20):
    """Hash a string with sha256 (used for the resume selection fingerprint)."""
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


GTDB_HOSTS = {'ecogenomic': 'data.gtdb.ecogenomic.org',
              'aau': 'data.gtdb.aau.ecogenomic.org'}


def gtdb_host(mirror=None):
    """Resolve the GTDB data host for a remote channel (plan RK13).

    Parameters
    ----------
    mirror : str, optional
        'ecogenomic', 'aau' or 'auto' (None/'' behave like 'auto')

    Returns
    -------
    str
        hostname to build URLs from

    Notes
    ----
    'auto' means "the host we can actually verify": the secondary AAU mirror
    presented an expired TLS certificate when last measured (2026-11-01), so
    auto and ecogenomic resolve to the primary host and choosing 'aau'
    explicitly prints that warning instead of failing inside a downloader.
    """
    key = (mirror or 'auto').lower()
    if key in ('auto', 'ecogenomic'):
        return GTDB_HOSTS['ecogenomic']
    if key == 'aau':
        print('  WARNING: --gtdb-mirror aau uses data.gtdb.aau.ecogenomic.org, '
              'whose TLS certificate was measured expired on 2026-11-01; '
              'downloads may fail on certificate validation (plan RK13).')
        return GTDB_HOSTS['aau']
    raise ValueError(
        f'Unknown --gtdb-mirror "{mirror}"; expected one of '
        + ', '.join(sorted(set(list(GTDB_HOSTS) + ['auto']))) + '.')


def fingerprint(files, out, with_lines=True, info=None):
    """Build an input fingerprint table (there is no local MD5SUM.txt).

    Parameters
    ----------
    files : iterable of str
        input files used by the build
    out : str
        destination TSV path (must be inside the workspace)
    with_lines : bool, optional
        also count lines, which requires a second read pass
    info : dict, optional
        layout information, recorded in the provenance header

    Returns
    -------
    dict
        path -> digest
    """
    out = guard_write(out)
    makedirs(dirname(out), exist_ok=True)
    res = {}
    with open(out, 'w') as f:
        f.write(provenance(info))
        f.write('# input fingerprints (local GTDB copies ship no MD5SUM.txt; '
                'see plan §3.4)\n')
        f.write('path\tsha256\tbytes\tlines\tmtime\n')
        for path in files:
            if not isfile(path):
                continue
            dig = sha256_file(path)
            res[realpath(path)] = dig
            n = _count_lines(path) if with_lines else -1
            st = os.stat(path)
            f.write(f'{realpath(path)}\t{dig}\t{st.st_size}\t{n}\t'
                    f'{datetime.fromtimestamp(st.st_mtime):%Y-%m-%dT%H:%M:%S}'
                    f'\n')
            print(f'  sha256 {dig[:16]}… {basename(path)}')
    return res


# ---------------------------------------------------------------------------
# Derived files and runtime annotation routes (plan §5.5)
# ---------------------------------------------------------------------------

def stem_index(info, domain=None):
    """Index local FASTA files by GTDB accession.

    Parameters
    ----------
    info : dict
        layout information
    domain : str, optional
        restrict to one domain ('bacteria' or 'archaea')

    Returns
    -------
    dict
        accession -> {'stem', 'path', 'domain'}
    """
    idx = {}
    for dom, v in info['domains'].items():
        if domain and dom != domain:
            continue
        for fname in _list_faa(v['dir']):
            stem = strip_faa_ext(fname)
            acc = stem.split('_d_')[0] if '_d_' in stem else stem
            idx[acc] = {'stem': stem, 'path': join(v['dir'], fname),
                        'domain': dom}
    return idx


def taxid_of_seqid(seqid, layout='local-renamed', stem2taxid=None,
                   prefix2genome=None, genome2taxid=None):
    """Resolve a hit sequence Id to a GTDB taxId (route C).

    Parameters
    ----------
    seqid : str
        subject sequence Id reported by the search engine
    layout : str, optional
        layout label, decides which anchor is used
    stem2taxid : dict, optional
        file-name stem -> taxId (local layout)
    prefix2genome : dict, optional
        nucleotide accession -> GTDB accession (official layout)
    genome2taxid : dict, optional
        GTDB accession -> species taxId

    Returns
    -------
    str or None
        taxId, or None when the Id cannot be resolved
    """
    parsed = parse_seqid(seqid, layout)
    if not parsed:
        return None
    if parsed.get('stem') and stem2taxid:
        # exact stem lookup only. The previous suffix-scan fallback (stems
        # longer than 512 characters matched by their last 64) was doubly
        # wrong: the real max_seqid_len measured on r220 builds is ~160, so
        # the branch was unreachable, and suffix matching cannot be unique --
        # same-genus stems share their tail, so the first dict hit could hand
        # another genome's taxId to a protein. An unresolvable id must surface
        # as "unresolved" (it is counted and listed), not as a wrong taxId.
        return stem2taxid.get(parsed['stem'])
    acc = parsed.get('genome')
    if acc and genome2taxid:
        return genome2taxid.get(acc)
    if parsed.get('contig') and prefix2genome and genome2taxid:
        acc = prefix2genome.get(parsed['contig']) or prefix2genome.get(
            parsed['contig_nover'] or '')
        return genome2taxid.get(acc) if acc else None
    return None


def write_side_files(db_dir, info, df, taxid_map, taxdump, stats,
                     stem2taxid=None, prefix2genome=None, manifest_path=None,
                     previous_manifest=None):
    """Write the gtdb/ derived files of a database directory.

    Parameters
    ----------
    db_dir : str
        database output directory (inside the workspace)
    info : dict
        layout information
    df : pandas.DataFrame
        genome table (all rows, not only the selected subset)
    taxid_map : dict
        node path -> taxId
    taxdump : dict
        taxonomy database
    stats : dict
        conversion statistics from build_taxdump_from_gtdb()
    stem2taxid : dict, optional
        file-name stem -> taxId for the genomes in this build
    prefix2genome : dict, optional
        nucleotide accession -> GTDB accession (collected while extracting)
    manifest_path : str, optional
        where to write taxid.map.gz (default: <db>/gtdb/taxid.map.gz)
    previous_manifest : str, optional
        the taxid.map.gz this build inherited its ids from; its reserved-ids
        sidecar is carried forward. Defaults to the manifest about to be
        written, which preserves history for same-directory rebuilds but must
        be passed explicitly when the ids come from an older build elsewhere
        (the database stage passes it; the historical default silently wiped
        the reserved set on the main build path).

    Returns
    -------
    dict
        role -> written path
    """
    db_dir = guard_write(db_dir, must_not_contain=info['gtdb_dir'])
    gtdb_dir = join(db_dir, 'gtdb')
    makedirs(gtdb_dir, exist_ok=True)
    tree = stats['tree']
    written = {}

    # One pass over the genomes aggregates every per-node NCBI hint, so the
    # cost stays linear in the number of genomes instead of nodes x rows.
    agg = {key: [0, set(), set()] for key in taxid_map}
    has_ncbi = 'ncbi_taxonomy' in df.columns and 'ncbi_taxid' in df.columns
    lins = df['gtdb_taxonomy'].tolist()
    ncbi_taxids = df['ncbi_taxid'].tolist() if has_ncbi else [''] * len(lins)
    ncbi_lins = df['ncbi_taxonomy'].tolist() if has_ncbi else [''] * len(lins)
    for lin, ntid, nlin in zip(lins, ncbi_taxids, ncbi_lins):
        segs = lin.split(LINEAGE_SEP)
        for depth in range(1, DEPTH + 1):
            rec = agg.get(LINEAGE_SEP.join(segs[:depth]))
            if rec is None:
                continue
            rec[0] += 1
            if ntid and ntid != 'na':
                rec[1].add(ntid)
            if nlin and nlin != 'na':
                nsegs = nlin.split(LINEAGE_SEP)
                # index 5 = the genus segment of the 7-rank NCBI lineage
                # (measured on r220: the column is always 7 prefixed ranks)
                if len(nsegs) > 5:
                    genus = nsegs[5]
                    # a bare 'g__' is GTDB's "unassigned", not a genus name;
                    # counting it as a distinct genus inflated nconflict on
                    # ~30% of the real r220 rows
                    if genus and genus != 'g__' and len(genus) > 3:
                        rec[2].add(genus)

    path = manifest_path or join(gtdb_dir, 'taxid.map.gz')
    with gzip.open(path, 'wt', encoding='utf-8') as f:
        f.write(provenance(info, extra={'id_base': stats['id_base']}))
        f.write('# taxId manifest: append-only, one row per lineage node\n')
        f.write('lineage\ttaxid\trank\tlabel\tname\tcluster_letter\t'
                'parent_taxid\tncbi_taxid_hint\n')
        for key in sorted(taxid_map, key=lambda k: int(taxid_map[k])):
            rec = tree[key]
            parent = (ROOT_TID if not rec['parent_key']
                      else taxid_map.get(rec['parent_key'], ''))
            ncbi = ','.join(sorted(agg[key][1])[:5]) or 'na'
            f.write('\t'.join((
                key, taxid_map[key], rec['rank'], rec['label'], rec['name'],
                rec['cluster_letter'], parent, ncbi)) + '\n')
    written['taxid.map.gz'] = path

    path = join(gtdb_dir, 'labels.tsv')
    with open(path, 'w') as f:
        f.write(provenance(info))
        f.write('taxid\tlabel\trank\tname\tcluster_letter\tn_genomes\n')
        for key, tid in sorted(taxid_map.items(), key=lambda kv: int(kv[1])):
            rec = tree[key]
            f.write(f'{tid}\t{rec["label"]}\t{rec["rank"]}\t{rec["name"]}\t'
                    f'{rec["cluster_letter"]}\t{agg[key][0]}\n')
    written['labels.tsv'] = path

    path = join(gtdb_dir, 'lineages.tsv')
    with open(path, 'w') as f:
        f.write(provenance(info))
        f.write('accession\tdomain\tgtdb_taxonomy\tspecies_taxid\n')
        for row in df.itertuples():
            f.write(f'{row.accession}\t{row.domain}\t{row.gtdb_taxonomy}\t'
                    f'{taxid_map.get(row.gtdb_taxonomy, "")}\n')
    written['lineages.tsv'] = path

    path = join(gtdb_dir, 'genome2taxid.tsv.gz')
    with gzip.open(path, 'wt', encoding='utf-8') as f:
        f.write(provenance(info, extra={'n_genomes': len(df)}))
        f.write('accession\tgenome\ttaxid\tdomain\tspecies_label\tproteins\t'
                'completeness\tcontamination\tquality_score\trepresentative\t'
                'type_species\tmimag_high\tassembly_level\tncbi_taxid\n')
        for row in df.itertuples():
            f.write('\t'.join(map(str, (
                row.accession, row.genome,
                taxid_map.get(row.gtdb_taxonomy, ''), row.domain,
                row.gtdb_species, row.proteins,
                getattr(row, 'completeness', ''),
                getattr(row, 'contamination', ''),
                getattr(row, 'quality_score', ''),
                row.gtdb_representative,
                getattr(row, 'gtdb_type_species_of_genus', ''),
                getattr(row, 'mimag_high_quality', ''),
                getattr(row, 'assembly_level', ''),
                getattr(row, 'ncbi_taxid', '')))) + '\n')
    written['genome2taxid.tsv.gz'] = path

    if stem2taxid:
        path = join(gtdb_dir, 'stem2taxid.tsv.gz')
        with gzip.open(path, 'wt', encoding='utf-8') as f:
            f.write(provenance(info, extra={'route': 'C'}))
            f.write('# route C: local file-name stem -> species taxId\n')
            f.write('stem\ttaxid\n')
            for stem in sorted(stem2taxid):
                f.write(f'{stem}\t{stem2taxid[stem]}\n')
        written['stem2taxid.tsv.gz'] = path

    if prefix2genome:
        acc2tid = {row.accession: taxid_map.get(row.gtdb_taxonomy, '')
                   for row in df.itertuples()}
        path = join(gtdb_dir, 'prefix2genome.tsv.gz')
        with gzip.open(path, 'wt', encoding='utf-8') as f:
            f.write(provenance(info, extra={'route': 'C-fallback'}))
            f.write('# route C fallback: nucleotide accession -> genome\n')
            f.write('nucleotide_accession\tgtdb_accession\ttaxid\n')
            for pacc in sorted(prefix2genome):
                acc = prefix2genome[pacc]
                f.write(f'{pacc}\t{acc}\t{acc2tid.get(acc, "")}\n')
        written['prefix2genome.tsv.gz'] = path

    path = join(gtdb_dir, 'ncbi_map.tsv.gz')
    with gzip.open(path, 'wt', encoding='utf-8') as f:
        f.write(provenance(info))
        f.write('gtdb_taxid\tgtdb_label\tgtdb_rank\tn_genomes\tnconflict\t'
                'ncbi_taxids\n')
        for key, tid in sorted(taxid_map.items(), key=lambda kv: int(kv[1])):
            rec = tree[key]
            n, taxids, genera = agg[key]
            f.write('\t'.join(map(str, (
                tid, rec['label'], rec['rank'], n,
                len(genera) if len(genera) > 1 else 0,
                ','.join(sorted(taxids)[:20]) or 'na'))) + '\n')
    written['ncbi_map.tsv.gz'] = path

    if stats['superseded']:
        path = join(gtdb_dir, 'delnodes.txt')
        with open(path, 'w') as f:
            f.write(provenance(info))
            f.write('# nodes present in the previous manifest but absent from '
                    'this release; retained so historical scores.tsv stay '
                    'interpretable (plan §5.3.3 rule 2)\n')
            for key, tid in sorted(stats['superseded'].items()):
                f.write(f'{tid}\t{key}\tsuperseded_by=removed\n')
        written['delnodes.txt'] = path

    # taxIds that this release no longer uses stay reserved for the next one
    # (plan §5.3.3 rule 2 / RK6): carrying them forward is what stops a renamed
    # taxon from inheriting a historical id and invalidating old scores.tsv
    reserved = dict(load_reserved(previous_manifest or manifest_path))
    reserved.update({tid: key for key, tid in
                     (stats.get('superseded') or {}).items()})
    if reserved:
        written['reserved_ids.txt'] = write_reserved(
            join(gtdb_dir, 'reserved_ids.txt'), sorted(
                reserved.items(), key=lambda kv: int(kv[0])), info)

    return written


def taxids_in_taxonmap(path):
    """Collect the taxIds referenced by a streamed taxonmap.

    Parameters
    ----------
    path : str
        DIAMOND taxonmap written during extraction

    Returns
    -------
    set of str
        taxIds in use (bounded by the number of taxa, not of proteins)
    """
    tids = set()
    with open(path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('accession'):
                continue
            x = line.rstrip('\r\n').split('\t')
            if len(x) >= 3:
                tids.add(x[2])
    return tids


def write_compact_taxdump(taxdump, tids, outdir):
    """Write only the taxa a taxonmap actually references.

    Parameters
    ----------
    taxdump : dict
        full taxonomy database
    tids : set of str
        taxIds in use
    outdir : str
        destination directory

    Returns
    -------
    tuple of (dict, dict)
        (written paths, refined taxdump)

    Notes
    ----
    DIAMOND loads every row of names.dmp at makedb time and spends about a
    minute per 146k nodes (spike S1), while a subset build only needs the
    handful of taxa it references. The full tree still goes to the database's
    taxdump/ directory, so taxId meanings never change.
    """
    from neoHGT.util import refine_taxdump
    used = {t for t in tids if t in taxdump}
    if not used:
        used = {ROOT_TID}
    sub = refine_taxdump(used, taxdump)
    sub.setdefault(ROOT_TID, {'parent': ROOT_TID, 'rank': 'no rank',
                              'name': 'root'})
    return write_taxdump(sub, outdir), sub


def load_labels(db_dir):
    """Read gtdb/labels.tsv from a database directory.

    Parameters
    ----------
    db_dir : str
        database directory

    Returns
    -------
    dict
        taxId -> official prefixed label (empty when absent)
    """
    path = join(db_dir, 'gtdb', 'labels.tsv')
    res = {}
    if not isfile(path):
        return res
    with open(path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('taxid\t'):
                continue
            x = line.rstrip('\r\n').split('\t')
            if len(x) >= 2:
                res[x[0]] = x[1]
    return res


_LABEL_INDEX = (None, {})


def label_index(labels):
    """taxId -> label reversed to label -> taxId, built once per mapping.

    The reverse lookup used to be spelled inline at three call sites as
    ``next((k for k, v in labels.items() if v == ref), None)``, which is linear
    in the number of taxa -- once per donor row. On a 146k-node database that is
    quadratic overall, so the index is cached; the first element of the cache
    tuple is the mapping it was built from, so a different (or rebuilt)
    database invalidates it. Ties keep the first taxId, as the inline scan did.
    """
    global _LABEL_INDEX
    src, idx = _LABEL_INDEX
    if src is not labels:
        idx = {}
        for tid, label in labels.items():
            idx.setdefault(label, tid)
        _LABEL_INDEX = (labels, idx)
    return idx


def resolve_label(ref, labels, taxid_map=None):
    """Resolve a taxon reference to its synthetic taxId (plan 附录 B).

    Parameters
    ----------
    ref : str
        a numeric taxId, or an official GTDB label such as 'g__Escherichia'
    labels : dict
        taxId -> label, as returned by load_labels()
    taxid_map : dict, optional
        node path -> taxId, consulted when the label is absent from `labels`
        (a subset build keeps only the nodes it selected)

    Returns
    -------
    str or None
        the taxId, or None when the reference is unknown here
    """
    ref = str(ref).strip()
    if not ref:
        return None
    if ref.isdigit():
        return ref
    tid = label_index(labels).get(ref)
    if tid is None and taxid_map is not None:
        tid = taxid_map.get(ref)
    return None if tid is None else str(tid)


def find_db_root(path):
    """Locate a GTDB-flavored neoHGT database root from any path inside it.

    Parameters
    ----------
    path : str
        a taxdump directory, db.faa, meta.json or the database root itself

    Returns
    -------
    str or None
        the database root, i.e. the directory whose gtdb/ subfolder holds a
        labels.tsv or taxid.map.gz and whose meta.json (or taxdump/) says
        this really is a neoHGT GTDB build

    Notes
    ----
    A real GTDB table is required as the marker. Keying on meta.json alone
    would misdetect an ordinary NCBI database, or an unrelated output
    directory, and silently switch it into GTDB mode. The marker alone is not
    enough either: an NCBI database nested under a directory that happens to
    hold a gtdb/ side folder used to flip the NCBI run into GTDB semantics.
    The root must therefore also carry meta.json whose source is gtdb/hybrid,
    or its own taxdump/ -- every build this project writes has both.
    """
    if not path:
        return None
    cur = realpath(expand(path))
    if isfile(cur):
        cur = dirname(cur)
    for _ in range(4):
        marker = join(cur, 'gtdb')
        if isdir(marker) and (isfile(join(marker, 'labels.tsv'))
                              or isfile(join(marker, 'taxid.map.gz'))):
            meta = read_meta(cur) or {}
            if str(meta.get('source', '')) in ('gtdb', 'hybrid') or                     isfile(join(cur, 'taxdump', 'nodes.dmp')):
                return cur
            # marker without corroboration: keep walking, this may be an
            # unrelated directory that merely contains a gtdb/ folder
        parent = dirname(cur)
        if parent == cur:
            break
        cur = parent
    return None


def load_genome2taxid(db_dir):
    """Read gtdb/genome2taxid.tsv.gz into an accession -> taxId dict."""
    path = join(db_dir, 'gtdb', 'genome2taxid.tsv.gz')
    res = {}
    if not isfile(path):
        return res
    with open_reader(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            header = line.rstrip('\r\n').split('\t')
            break
        else:
            return res
        if 'accession' not in header or 'taxid' not in header:
            raise ValueError(
                f'{path}: expected an "accession...taxid" header, got '
                f'{header[:3]}')
        ia, it = header.index('accession'), header.index('taxid')
        for line in f:
            if line.startswith('#'):
                continue
            x = line.rstrip('\r\n').split('\t')
            if len(x) > max(ia, it) and x[it]:
                res[x[ia]] = x[it]
    return res


def load_prefix2genome(db_dir):
    """Read gtdb/prefix2genome.tsv.gz (nucleotide accession -> genome)."""
    path = join(db_dir, 'gtdb', 'prefix2genome.tsv.gz')
    res = {}
    if not isfile(path):
        return res
    with open_reader(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            x = line.rstrip('\r\n').split('\t')
            if len(x) >= 2:
                res[x[0]] = x[1]
    return res


def read_meta(db_dir):
    """Read meta.json of a database directory ({} when absent).

    Only this directory's meta.json is read: the previous fallback to the
    parent directory could adopt a *neighbouring* database's metadata and
    silently pick the wrong seqid parser. A file that exists but does not
    parse raises -- a truncated meta.json used to become "no metadata" one
    level up and a bare traceback here, depending on which hit first.
    """
    path = join(db_dir, 'meta.json')
    if not isfile(path):
        return {}
    try:
        with open(path) as f:
            return json.load(f)
    except ValueError as exc:
        raise ValueError(
            f'{path} is not valid JSON ({exc}); the database metadata is '
            'damaged.') from None


def load_stem2taxid(db_dir):
    """Read gtdb/stem2taxid.tsv.gz (route C table) from a database directory."""
    path = join(db_dir, 'gtdb', 'stem2taxid.tsv.gz')
    res = {}
    if not isfile(path):
        return res
    with open_reader(path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('stem\t'):
                continue
            x = line.rstrip('\r\n').split('\t')
            if len(x) >= 2:
                res[x[0]] = x[1]
    return res


# ---------------------------------------------------------------------------
# Phylogeny-assisted grouping inputs (plan P7, §2.4 A4): the local release also
# ships per-rank collapsed trees (tree_{bac120,ar53}/*_{phylum,...,genus}.tree).
# Patristic distances from those trees give an objective definition of "how many
# taxa are close", which is what grouping.profiles.*.closesize tries to proxy.
# ---------------------------------------------------------------------------

class NewickTree(object):
    """Small Newick reader: labels, branch lengths, patristic distances.

    Notes
    -----
    Supports what the GTDB collapsed trees use: nested parentheses, internal
    labels that may be support values (")100.0:0.05"), and branch lengths in
    scientific notation. No third-party dependency is introduced.
    """

    def __init__(self, text):
        self.children = {}      # node -> [child node]
        self.parent = {}        # node -> parent node
        self.branch = {}        # node -> length of the branch to its parent
        self.name = {}          # node -> label
        self._n = 0
        text = text.strip()
        if text.endswith(';'):
            text = text[:-1]
        self.root, pos = self._parse(text, 0)
        rest = text[pos:].strip()
        if rest:
            raise ValueError(f'unparsed trailing Newick content: "{rest[:40]}"')

    def _new(self):
        self._n += 1
        return self._n

    def _parse(self, s, i):
        """Parse one (sub)tree starting at s[i]; return (node, next index)."""
        node = self._new()
        if s[i] == '(':
            i += 1
            while True:
                child, i = self._parse(s, i)
                self.children.setdefault(node, []).append(child)
                self.parent[child] = node
                if i < len(s) and s[i] == ',':
                    i += 1
                    continue
                if i < len(s) and s[i] == ')':
                    i += 1
                    break
                raise ValueError(f'unbalanced parentheses at offset {i}')
        j = i
        while j < len(s) and s[j] not in '(),:':
            j += 1
        label = s[i:j].strip()
        if label:
            self.name[node] = label
        if j < len(s) and s[j] == ':':
            k = j + 1
            while k < len(s) and s[k] not in '(),:':
                k += 1
            try:
                self.branch[node] = float(s[j + 1:k])
            except ValueError:
                raise ValueError(f'invalid branch length "{s[j + 1:k]}"')
            i = k
        else:
            i = j
        return node, i

    def leaves(self):
        """Node ids without children."""
        return [n for n in range(1, self._n + 1) if n not in self.children]

    def by_name(self):
        """Label -> node id (labels are unique within one rank tree)."""
        return {v: k for k, v in self.name.items()}

    def _path_to_root(self, node):
        """[(node, distance from root)] from the root down to this node."""
        chain, total, cur = [], 0.0, node
        steps = []
        while cur in self.parent:
            steps.append(cur)
            total += self.branch.get(cur, 0.0)
            cur = self.parent[cur]
        steps.append(cur)
        # 'total' is the distance from the root to the original node
        return list(reversed(steps)), total

    def distance(self, a, b):
        """Patristic distance between two nodes."""
        pa, da = self._path_to_root(a)
        pb, db = self._path_to_root(b)
        common = self.root
        for x, y in zip(pa, pb):
            if x != y:
                break
            common = x
        dc = self._path_to_root(common)[1]
        return da + db - 2 * dc

    def distances_from(self, node):
        """Patristic distance from one node to every labelled leaf."""
        out = {}
        for leaf in self.leaves():
            label = self.name.get(leaf)
            if not label or leaf == node:
                continue
            out[label] = self.distance(node, leaf)
        return out


class TaxidIndex(object):
    """Dict-like random access to a sorted protein-to-taxId TSV without loading it.

    Parameters
    ----------
    path : str
        plain (not gzipped) file of 'protein<TAB>taxId' lines, byte-sorted by
        the first column - the order in which both database routes write it
    spacing : int
        number of lines between two index anchors

    Notes
    -----
    This is plan P7's "feasible route B": a full r223-scale map has ~600M rows,
    which the upstream 'read_prot2taxid' dict cannot hold. An anchor table of
    one key per 'spacing' lines costs 0.2% of the data and is cached next to the
    map as '<path>.idx', so a lookup is one binary search plus a short scan.
    """

    def __init__(self, path, spacing=2048):
        if not isfile(path):
            raise ValueError(f'{path}: file not found.')
        if path.endswith('.gz'):
            raise ValueError(
                f'{path} is gzipped: the index needs byte offsets. Decompress '
                'it once (gunzip -k) and point --taxmap at the plain file.')
        self.path = path
        self.spacing = int(spacing)
        self.index_path = path + '.idx'
        self.keys, self.offsets, self.n_lines = [], [], 0
        self._fh = None
        if isfile(self.index_path) and self._load_index():
            return
        try:
            self._build_index()
        except OSError as exc:
            raise OSError(
                f'cannot build the --taxmap index next to the data file '
                f'({exc}); the map\'s directory may be read-only. Copy the '
                'map somewhere writable, or pass --taxmap through a writable '
                'path.') from None

    def _identity(self):
        """Cheap identity probe of the data file (size + head checksum)."""
        size = getsize(self.path)
        h = hashlib.sha256()
        with open(self.path, 'rb') as f:
            h.update(f.read(65536))
        return f'{size}:{h.hexdigest()[:16]}'

    def _write_index(self):
        # the header records the data file's identity: without it, pointing
        # --taxmap at a different (or truncated) map while the old .idx
        # survived silently served wrong offsets and therefore wrong taxIds
        with open(self.index_path, 'w') as f:
            f.write(f'#taxid-index\t{self.spacing}\t{self.n_lines}\t'
                    f'{self._identity()}\n')
            for k, o in zip(self.keys, self.offsets):
                f.write(f'{k}\t{o}\n')

    def _load_index(self):
        try:
            with open(self.index_path) as f:
                head = f.readline().rstrip('\n').split('\t')
                if len(head) < 3 or int(head[1]) != self.spacing:
                    return False
                self.n_lines = int(head[2])
                if len(head) >= 4 and head[3] != self._identity():
                    print(f'  NOTE: {self.index_path} was built for a '
                          f'different {basename(self.path)}; rebuilding.')
                    return False
                if len(head) < 4:
                    # legacy index without an identity header: rebuild once,
                    # the size check alone has already caught the worst case
                    return False
                self.keys, self.offsets = [], []
                for line in f:
                    k, _, o = line.rstrip('\n').rpartition('\t')
                    self.keys.append(k)
                    self.offsets.append(int(o))
            return bool(self.keys)
        except (ValueError, OSError):
            return False

    def _build_index(self):
        keys, offsets = [], []
        n = 0
        with open(self.path, 'rb') as f:
            prev = None
            while True:
                off = f.tell()
                line = f.readline()
                if not line:
                    break
                if not line.strip():
                    continue
                key = line.split(b'\t', 1)[0]
                if prev is not None and key < prev:
                    raise ValueError(
                        f'{self.path} is not sorted at line {n + 1}; sort it '
                        'with LC_ALL=C sort.')
                prev = key
                if n % self.spacing == 0:
                    keys.append(key.decode('utf-8', 'replace'))
                    offsets.append(off)
                n += 1
        self.keys, self.offsets, self.n_lines = keys, offsets, n
        if not keys:
            raise ValueError(f'{self.path}: no usable protein-to-taxId rows.')
        # the index lives next to the data file, which may sit inside the
        # read-only GTDB dataset; go through the guard so a read-only input
        # tree never gets a stray .idx written into it
        self.index_path = guard_write(self.index_path)
        self._write_index()

    def _file(self):
        if self._fh is None:
            self._fh = open(self.path, 'rb')
        return self._fh

    def close(self):
        if self._fh is not None:
            self._fh.close()
            self._fh = None

    def get(self, key, default=None):
        i = bisect_right(self.keys, key) - 1
        if i < 0:
            i = 0
        f = self._file()
        f.seek(self.offsets[i])
        for _ in range(self.spacing + 1):
            line = f.readline()
            if not line:
                break
            k, _, v = line.rstrip(b'\r\n').partition(b'\t')
            k = k.decode('utf-8', 'replace')
            if k == key:
                return v.decode('utf-8', 'replace')
            if k > key:
                break
        return default

    def __getitem__(self, key):
        got = self.get(key)
        if got is None:
            raise KeyError(key)
        return got

    def __contains__(self, key):
        return self.get(key) is not None

    def __len__(self):
        return self.n_lines

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


def read_newick(path):
    """Read a Newick tree file (plain or gzipped)."""
    with open_reader(path) as f:
        text = f.read()
    if not text.strip():
        raise ValueError(f'{path}: empty Newick file')
    return NewickTree(text)


def rank_trees(gtdb_dir):
    """Locate the per-rank collapsed trees shipped with a GTDB release.

    Returns
    -------
    dict
        rank -> [paths]; 'all' holds the whole-tree files
    """
    res = {}
    ranks = ('phylum', 'class', 'order', 'family', 'genus')
    if not isdir(gtdb_dir):
        raise ValueError(
            f'rank_trees: not a directory: {gtdb_dir}. Pass --gtdb-dir '
            'pointing at a GTDB release root that holds tree_*/ folders.')
    for sub in sorted(os.listdir(gtdb_dir)):
        d = join(gtdb_dir, sub)
        if not sub.startswith('tree_') or not isdir(d):
            continue
        for fname in sorted(os.listdir(d)):
            if not fname.endswith('.tree'):
                continue
            stem = fname[:-len('.tree')]
            hit = [r for r in ranks if stem.endswith('_' + r)]
            if hit:
                res.setdefault(hit[0], []).append(join(d, fname))
            else:
                res.setdefault('all', []).append(join(d, fname))
    return res


def tree_close_taxa(tree_path, focus, max_distance=None, top_k=None):
    """Rank taxa by phylogenetic distance to a focus taxon in one tree.

    Parameters
    ----------
    tree_path : str
        Newick file, typically a GTDB per-rank collapsed tree
    focus : str
        label of a node in that tree, e.g. 'g__Halorubrum'
    max_distance : float, optional
        keep only taxa at or below this patristic distance
    top_k : int, optional
        keep only the k nearest taxa

    Returns
    -------
    tuple
        (rows, summary): rows is [(label, distance)] sorted by distance,
        summary reports the distance distribution and how many taxa fall below
        common cut-offs
    """
    tree = read_newick(tree_path)
    by_name = tree.by_name()
    if focus not in by_name:
        sample = ', '.join(sorted(by_name)[:3])
        raise ValueError(
            f'"{focus}" is not a labelled node of {basename(tree_path)}. '
            f'Example labels: {sample}.')
    dists = tree.distances_from(by_name[focus])
    rows = sorted(dists.items(), key=lambda x: (x[1], x[0]))
    all_rows = rows
    if max_distance is not None:
        rows = [r for r in rows if r[1] <= float(max_distance)]
    if top_k is not None:
        rows = rows[:int(top_k)]
    d = [v for _, v in all_rows]
    summary = {'tree': tree_path, 'focus': focus, 'n_taxa': len(all_rows),
               'min': min(d) if d else None, 'max': max(d) if d else None,
               'median': statistics.median(d) if d else None,
               'within': {str(c): sum(1 for v in d if v <= c)
                          for c in (0.01, 0.05, 0.1, 0.2, 0.3, 0.5)},
               'returned': len(rows)}
    return rows, summary


# ---------------------------------------------------------------------------
# Lineage subset expansion beyond representatives (plan P7). The official
# sp_clusters.tsv is absent from local copies, so the member -> representative
# relation is taken from the metadata column 'gtdb_genome_representative'.
# ---------------------------------------------------------------------------

def member_index(df):
    """Representative accession -> member accessions, from metadata.

    Parameters
    ----------
    df : pandas.DataFrame
        genome table as loaded by load_metadata()

    Returns
    -------
    dict
        representative -> list of member accessions (the representative itself
        is included, which is how GTDB metadata encodes it)
    """
    col = 'gtdb_genome_representative'
    if col not in df.columns:
        raise ValueError(
            f'cannot expand to members: metadata has no "{col}" column.')
    idx = {}
    for rep, acc in zip(df[col].tolist(), df['accession'].tolist()):
        idx.setdefault(str(rep), []).append(str(acc))
    return idx


def expand_to_members(df, sel, seq_index=None):
    """Expand a genome selection to the member genomes it stands for.

    Parameters
    ----------
    df : pandas.DataFrame
        all genomes (the metadata table)
    sel : pandas.DataFrame
        the current selection
    seq_index : dict, optional
        accession -> local sequence path; when given, members are marked with
        whether their sequence is actually available

    Returns
    -------
    tuple
        (members, report) where 'members' is a DataFrame holding every member
        genome (columns of df plus 'representative' and, when seq_index was
        given, 'has_local_fasta') and 'report' summarises the expansion
    """
    idx = member_index(df)
    selected = {str(a) for a in sel['accession']}
    pairs = [(rep, m) for rep in sorted(selected) for m in idx.get(rep, [])
             if m not in selected]
    report = {'representatives': len(selected), 'members': len(pairs)}
    if not pairs:
        report.update({'members_in_metadata': 0, 'members_with_sequence': 0})
        return df.iloc[0:0].assign(representative=[]), report
    want = {m for _, m in pairs}
    members = df[df['accession'].astype(str).isin(want)].copy()
    rep_of = {m: rep for rep, m in pairs}
    members['representative'] = [rep_of.get(str(a), '')
                                for a in members['accession']]
    if seq_index is not None:
        members['has_local_fasta'] = [bool(seq_index.get(str(a)))
                                      for a in members['accession']]
        report['members_with_sequence'] = int(members['has_local_fasta'].sum())
    report['members_in_metadata'] = int(members.shape[0])
    return members, report


# ---------------------------------------------------------------------------
# GTDB-Tk results as an input format (plan P7, N2: no runtime gtdbtk dependency)
# ---------------------------------------------------------------------------

GTDBTK_GENOME_COLS = ('user_genome', 'Genome', 'genome')
GTDBTK_LINEAGE_COLS = ('classification', 'classify_info', 'gtdbtk.lineage',
                       'lineage')
GTDBTK_RANK_COLS = ('d__Domain', 'p__Phylum', 'c__Class', 'o__Order',
                    'f__Family', 'g__Genus', 's__Species')
RANK_DEPTH = {r: i + 1 for i, r in enumerate(GTDB_RANKS)}


def parse_gtdbtk_lineage(text):
    """Split a GTDB-Tk lineage into its leading well-formed rank segments.

    Parameters
    ----------
    text : str
        lineage such as 'd__Bacteria;p__X;c__Y;NA;NA;NA;NA'

    Returns
    -------
    list of str
        segments down to the last assigned rank (may be shorter than 7)

    Raises
    ------
    ValueError
        If a segment carries the wrong rank prefix (i.e. is not a GTDB lineage)
    """
    segs = [str(s).strip() for s in str(text).split(LINEAGE_SEP)]
    out = []
    for rank, seg in zip(GTDB_RANKS, segs):
        bare = seg[len(GTDB_RANK_PREFIXES[rank]):].strip() if seg.startswith(
            GTDB_RANK_PREFIXES[rank]) else seg
        if not seg or bare.upper() in ('NA', 'N/A', 'NAN', 'UNASSIGNED', ''):
            break
        char = seg[:1].lower()
        if PREFIX_BY_CHAR.get(char) != rank or not seg.startswith(char):
            raise ValueError(
                f'Segment "{seg}" does not carry the expected rank prefix '
                f'"{GTDB_RANK_PREFIXES[rank]}" for rank "{rank}".')
        out.append(seg)
    return out


def load_gtdbtk_summary(path):
    """Read a GTDB-Tk classify or de novo summary table.

    Returns
    -------
    tuple
        (rows, meta) where rows is a list of (genome, lineage text) and meta
        records which columns were used

    Notes
    ----
    Supported shapes: v2 'classify.max_summary' ('classification'), v1
    'classify.summary' ('classify_info'), de novo summaries ('gtdbtk.lineage'
    or 'lineage'), and summaries that carry the seven rank columns instead of
    one lineage column. v2 banners ('#####') are skipped but the '#user_genome'
    header line is kept and only its leading '#' is dropped.
    """
    with open_reader(path) as f:
        text = '\n'.join(l for l in f.read().splitlines()
                         if not l.startswith('#####'))
    df = pd.read_csv(io.StringIO(text), sep='\t', dtype=str,
                     keep_default_na=False)
    df.columns = [str(c).lstrip('#').strip() for c in df.columns]
    gcol = next((c for c in GTDBTK_GENOME_COLS if c in df.columns),
                df.columns[0])
    lcols = [c for c in GTDBTK_LINEAGE_COLS if c in df.columns]
    rows = []
    if not lcols and all(c in df.columns for c in GTDBTK_RANK_COLS):
        used = '<rank columns>'
        for rec in df[[gcol, *GTDBTK_RANK_COLS]].itertuples(index=False):
            rows.append((rec[0], LINEAGE_SEP.join(rec[1:])))
    elif lcols:
        # v2 leaves 'classification' empty whenever a rank is unassigned, so
        # each row falls back to the next lineage column that has a value
        used = ','.join(lcols)
        sub = df[[gcol, *lcols]]
        for rec in sub.itertuples(index=False):
            lineage = next((v for v in rec[1:] if str(v).strip()), '')
            rows.append((rec[0], lineage))
    else:
        raise ValueError(
            f'{path}: no GTDB-Tk lineage column found. Expected one of '
            f'{", ".join(GTDBTK_LINEAGE_COLS)} or the seven rank columns '
            f'{", ".join(GTDBTK_RANK_COLS)}.')
    return rows, {'genome_col': gcol, 'lineage_col': used, 'n_rows': len(df)}


def gtdbtk_to_input_tax(summary, db_dir, out=None, min_rank=None):
    """Map GTDB-Tk lineages onto the synthetic taxIds of a GTDB database.

    Parameters
    ----------
    summary : str
        path of a GTDB-Tk summary table
    db_dir : str
        database built with --source gtdb (its gtdb/taxid.map.gz is the map)
    out : str, optional
        destination 'genome<TAB>taxId' file for neoHGT analyze --input-tax
    min_rank : str, optional
        drop genomes whose deepest assigned rank is below this rank

    Returns
    -------
    dict
        report with counts, the rank at which each genome matched, and every
        unmapped genome with its reason
    """
    man = load_manifest(join(db_dir, 'gtdb', 'taxid.map.gz'))
    if not man:
        raise ValueError(
            f'{db_dir} has no gtdb/taxid.map.gz: build the database with '
            'neoHGT --source gtdb first, it is the taxId source here.')
    rows, meta = load_gtdbtk_summary(summary)
    want_depth = RANK_DEPTH.get(min_rank) if min_rank else None
    if min_rank and want_depth is None:
        raise ValueError(f'Unknown --min-rank "{min_rank}".')

    mapped, unmapped = {}, []
    by_rank = {}
    for genome, lineage in rows:
        try:
            segs = parse_gtdbtk_lineage(lineage)
        except ValueError as e:
            unmapped.append({'genome': genome, 'lineage': lineage,
                             'reason': str(e)})
            continue
        if not segs:
            unmapped.append({'genome': genome, 'lineage': lineage,
                             'reason': 'no rank assigned'})
            continue
        hit = None
        for depth in range(len(segs), 0, -1):
            rec = man.get(LINEAGE_SEP.join(segs[:depth]))
            if rec:
                hit = (rec['taxid'], GTDB_RANKS[depth - 1], segs[depth - 1])
                break
        if hit is None:
            unmapped.append({'genome': genome, 'lineage': lineage,
                             'reason': 'lineage not in this database'})
            continue
        if want_depth is not None and RANK_DEPTH[hit[1]] < want_depth:
            unmapped.append({'genome': genome, 'lineage': lineage,
                             'reason': f'deepest rank {hit[1]} is below '
                                       f'--min-rank {min_rank}'})
            continue
        mapped[genome] = hit
        by_rank[hit[1]] = by_rank.get(hit[1], 0) + 1

    report = {'summary': summary, 'db': db_dir, 'columns': meta,
              'min_rank': min_rank, 'n_rows': len(rows),
              'n_mapped': len(mapped), 'n_unmapped': len(unmapped),
              'matched_rank': by_rank, 'unmapped': unmapped}
    if out:
        path = guard_write(out, must_not_contain=db_dir)
        makedirs(dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            for genome in sorted(mapped):
                f.write(f'{genome}\t{mapped[genome][0]}\n')
        report['input_tax'] = path
        report['deepest_label'] = {g: mapped[g][2] for g in sorted(mapped)}
    return report


# ---------------------------------------------------------------------------
# Source providers (plan §5.4)
# ---------------------------------------------------------------------------

class GenomeFetch(tuple):
    """One genome's sequence location."""

    __slots__ = ()
    _fields = ('accession', 'genome', 'taxid', 'path', 'url', 'kind')

    def __new__(cls, accession, genome=None, taxid=None, path=None, url=None,
                kind='local'):
        # genome/taxid default to None so copy/pickle's no-arg __new__ probe
        # works (the positional-only form broke copy.copy/deepcopy)
        return super().__new__(cls, (accession, genome, taxid, path, url, kind))

    def __getattr__(self, name):
        try:
            return self[self._fields.index(name)]
        except ValueError:
            # AttributeError is what hasattr()/copy/pickle probe for; letting
            # ValueError escape broke every protocol that asks for a missing
            # attribute politely
            raise AttributeError(
                f'{type(self).__name__} has no field {name!r}') from None


def read_idlist(path):
    """Read a GTDB accession list (one per line, '#' comments allowed)."""
    out = []
    with open_reader(path) as f:
        for line in f:
            line = line.split('#')[0].strip()
            if line and not line.startswith('accession'):
                out.append(line)
    return out


class SourceProvider:
    """Genome list, taxonomy and sequence retrieval for one data source."""

    name = 'ncbi'
    ranks = NCBI_RANKS
    rank_prefixes = {r: '' for r in NCBI_RANKS}
    default_taxonomy_opts = None
    # filled in by select(); see the plan §7 note on reporting actual counts
    quality_classes = {}

    def __init__(self, opts=None):
        self.opts = opts or {}
        self.df = None
        self.taxdump = None
        self.taxid_map = None
        self.stats = None
        self.info = None

    def load_taxonomy(self):
        raise NotImplementedError

    def load_genomes(self):
        raise NotImplementedError

    def select(self, df, opts=None):
        raise NotImplementedError

    def sample(self, df, opts=None):
        raise NotImplementedError

    def protein_sources(self, df):
        raise NotImplementedError


def selection_opts(args):
    """Map a parsed CLI namespace onto the opts GtdbProvider.select() reads.

    This used to be an inline dict in the build path, which is how `--complete`
    and `--typemater` ended up parsed by argparse but never handed to the
    selector: the flags looked implemented from the outside.

    Parameters
    ----------
    args : argparse.Namespace
        parsed `database` arguments

    Returns
    -------
    dict
        opts for GtdbProvider.select() / sample()
    """
    return {
        'gtdb_taxa': getattr(args, 'gtdb_taxa', None),
        'reps': getattr(args, 'reps', False),
        'genome_list': getattr(args, 'genome_list', None),
        'taxids': getattr(args, 'taxids', None),
        'genoids': getattr(args, 'genoids', None),
        'exclude': getattr(args, 'exclude', False),
        'min_completeness': getattr(args, 'min_completeness', None),
        'max_contamination': getattr(args, 'max_contamination', None),
        'mimag': getattr(args, 'mimag', None),
        'type_material': getattr(args, 'type_material', None),
        'typemater': bool(getattr(args, 'typemater', False)),
        'complete': bool(getattr(args, 'complete', False)),
        'domain': None if getattr(args, 'domain', 'both') == 'both'
        else getattr(args, 'domain', None),
        'rank_sample': getattr(args, 'rank_sample', None),
        'rank': getattr(args, 'rank', 'species'),
        'sample': getattr(args, 'sample', None),
        'above': bool(getattr(args, 'above', False)),
    }


class GtdbProvider(SourceProvider):
    """GTDB data source backed by a local read-only dataset.

    Parameters
    ----------
    gtdb_dir : str
        dataset root (read-only)
    release : str, optional
        GTDB release label
    layout : str, optional
        forced layout ('auto' by default)
    id_base : int, optional
        start of the synthetic taxId range
    manifest_path : str, optional
        previous taxid.map.gz to reuse for ID stability
    selfcheck : bool, optional
        run the §3.6 consistency checks before anything else
    """

    name = 'gtdb'
    ranks = GTDB_RANKS
    rank_prefixes = GTDB_RANK_PREFIXES
    default_taxonomy_opts = {'capital': False, 'latin': False, 'block': []}

    def __init__(self, gtdb_dir, release=None, layout='auto',
                 id_base=DEFAULT_ID_BASE, manifest_path=None, selfcheck=True,
                 strict_selfcheck=True, opts=None, metadata_files=None,
                 taxonomy_files=None, require_metadata=True,
                 require_sequences=True):
        super().__init__(opts)
        self.gtdb_dir = realpath(expand(gtdb_dir))
        self.id_base = int(id_base)
        self.manifest_path = manifest_path
        self.want_selfcheck = selfcheck
        self.strict_selfcheck = strict_selfcheck
        self.require_metadata = require_metadata
        self.require_sequences = require_sequences
        self.report = None
        self.override_layout = layout
        self.override_release = release
        self.override_metadata = metadata_files
        self.override_taxonomy = taxonomy_files

    def prepare(self):
        """Probe the dataset, check consistency and load the tables."""
        assert_input_read_only(self.gtdb_dir)
        self.info = detect_layout(self.gtdb_dir, release=self.override_release,
                                  layout=self.override_layout,
                                  metadata=self.override_metadata,
                                  taxonomy=self.override_taxonomy,
                                  require_metadata=self.require_metadata,
                                  require_sequences=self.require_sequences)
        self.info = apply_table_overrides(self.info, self.override_metadata,
                                           self.override_taxonomy)
        self.lineage_only = not self.info['metadata']
        if self.lineage_only:
            print('NOTICE: no metadata table under', self.gtdb_dir, '-- the '
                  'taxdump is built from the lineage files alone, so the '
                  'genome maps stay empty and a database build (--source gtdb) '
                  'is not possible from this dataset.')
        if self.want_selfcheck:
            if self.lineage_only:
                self.report = selfcheck_lineages(
                    self.info, strict=self.strict_selfcheck)
            else:
                self.report = selfcheck(self.info,
                                        strict=self.strict_selfcheck)
        else:
            print('WARNING: input self-consistency checks were skipped by '
                  'request; the build is not guaranteed to be reproducible.')
        if self.lineage_only:
            self.df_all = empty_genome_table()
            self.lineages = sorted(set(read_taxonomy(
                self.info['taxonomy']).values()))
        else:
            self.df_all = load_metadata(self.info['metadata'],
                                        release=self.info['release'])
            self.lineages = list(self.df_all['gtdb_taxonomy'].unique())
        return self.info

    def load_taxonomy(self, trim=None):
        """Build the taxdump for the (optionally trimmed) lineage set.

        Parameters
        ----------
        trim : set of str, optional
        keep only these node paths plus their ancestors

        Returns
        -------
        tuple
            (taxdump, taxid_map, stats)
        """
        lins = self.lineages
        if trim:
            keep = set()
            for lin in lins:
                if lin in trim:
                    keep.update(lineage_paths(lin))
            nodes = collect_nodes(lins, keep=keep)
            taxid_map, superseded = assign_ids(
                nodes, load_manifest(self.manifest_path), self.id_base,
                reserve=load_reserved(self.manifest_path))
            taxdump = {ROOT_TID: {'parent': ROOT_TID, 'rank': 'no rank',
                                  'name': 'root'}}
            for key, rec in nodes.items():
                taxdump[taxid_map[key]] = {
                    'parent': ROOT_TID if not rec['parent_key'] else
                    taxid_map.get(rec['parent_key'], ROOT_TID),
                    'rank': rec['rank'], 'name': rec['name']}
            self.taxdump, self.taxid_map = taxdump, taxid_map
            self.stats = {'n_nodes': len(taxdump), 'n_gtdb_nodes': len(nodes),
                          'by_rank': {r: sum(1 for v in nodes.values()
                                             if v['rank'] == r)
                                       for r in GTDB_RANKS},
                          'n_superseded': len(superseded),
                          'id_base': self.id_base, 'tree': nodes,
                          'superseded': superseded}
            return taxdump, taxid_map, self.stats

        manifest = load_manifest(self.manifest_path)
        self.taxdump, self.taxid_map, self.stats = build_taxdump_from_gtdb(
            lins, manifest=manifest, id_base=self.id_base,
            manifest_path=self.manifest_path)
        return self.taxdump, self.taxid_map, self.stats

    def load_genomes(self):
        """Get the normalized genome table."""
        if self.df is None:
            self.df = self.df_all.copy()
        return self.df

    def node_key_of(self, ref):
        """Resolve a taxon reference to a lineage node path.

        Parameters
        ----------
        ref : str
            'f__Enterobacteriaceae', a taxId, or a bare name

        Returns
        -------
        str
            node path
        """
        ref = str(ref).strip()
        if ref.isdigit():
            for key, tid in self.taxid_map.items():
                if tid == ref:
                    return key
            raise ValueError(f'GTDB taxId {ref} is not in this database.')
        rank = rank_of_label(ref)
        bare = strip_rank_prefix(ref)
        if rank:
            want = GTDB_RANK_PREFIXES[rank] + bare
            for key in self.taxid_map:
                if path_label(key) == want or key.endswith(LINEAGE_SEP + want):
                    return key
            # tolerate local underscore form ('f_Enterobacteriaceae')
            want2 = GTDB_RANK_PREFIXES[rank] + bare.replace(' ', '_')
            for key in self.taxid_map:
                if path_label(key) == want2:
                    return key
        cand = [key for key in self.taxid_map
                if path_label(key) in (bare, f'{bare}')
                or strip_rank_prefix(path_label(key)) == bare]
        if len(cand) == 1:
            return cand[0]
        if len(cand) > 1:
            raise ValueError(
                f'Taxon "{ref}" matches {len(cand)} GTDB nodes; use a prefixed '
                f'label (e.g. {GTDB_RANK_PREFIXES.get(rank_of_label(cand[0]),"")} '
                f'{bare}), a taxId, or a full lineage string.')
        raise ValueError(f'Taxon "{ref}" is not found in the GTDB taxonomy.')

    def select(self, df=None, opts=None):
        """Apply lineage, quality and list filters to the genome table.

        Parameters
        ----------
        df : pandas.DataFrame, optional
            genome table (defaults to the provider's table)
        opts : dict, optional
            selection options (plan §7.1)

        Returns
        -------
        pandas.DataFrame
            selected genomes
        """
        df = df if df is not None else self.load_genomes()
        o = dict(self.opts)
        o.update(opts or {})
        sel = df
        keys, labels = [], []

        # plan §7 asks the build report to state how many genomes actually carry
        # each quality class instead of leaving the reader to guess; count them
        # on the unfiltered table so the numbers describe the dataset, not the
        # selection that happened to be requested.
        classes = {}
        if 'gtdb_type_species_of_genus' in df:
            classes['gtdb_type_species'] = int(
                (df['gtdb_type_species_of_genus'] == 't').sum())
        if 'ncbi_type_material_designation' in df:
            classes['ncbi_type_material'] = int(
                (~_designation_is_empty(
                    df['ncbi_type_material_designation'])).sum())
        if 'mimag_high_quality' in df:
            classes['mimag_high_quality'] = int(
                (df['mimag_high_quality'] == 't').sum())
        self.quality_classes = classes

        for ref in _as_list(o.get('gtdb_taxa')):
            keys.append(self.node_key_of(ref))
        if keys:
            # Every requested node must participate in the mask. Keying a dict
            # by rank ordinal (the previous implementation) silently collapsed
            # same-rank requests to one survivor -- and which one survived
            # depended on set() iteration order, i.e. on PYTHONHASHSEED.
            # Semantics: same-rank requests are a union; a request at a deeper
            # rank is NOT restricted to a shallower one (each prefix matches on
            # its own, so cross-rank requests are also a union).
            wanted = set(keys)
            if len(wanted) != len(keys):
                print(f'  NOTE: {len(keys) - len(wanted)} duplicate '
                      '--gtdb-taxa entries ignored.')
            mask = sel['gtdb_taxonomy'].apply(
                lambda lin: any(
                    LINEAGE_SEP.join(lin.split(LINEAGE_SEP)[:d + 1]) in wanted
                    for d in range(DEPTH)))
            sel = sel[mask]
            labels.append(f'lineage({len(wanted)} taxa)')

        if o.get('reps'):
            sel = sel[sel['is_rep']]
            labels.append('representatives')

        if o.get('genome_list'):
            accs = set()
            for path in _as_list(o['genome_list']):
                for item in read_idlist(path):
                    accs.add(item)
            # `_expand_acc` documents that a bare 'GCA_000000003' is accepted,
            # but matching it against the versioned `accession` column would
            # then find nothing. Compare unversioned, prefixed forms on both
            # sides so every spelling in the docstring actually selects.
            def acc_key(a):
                return re.sub(r'\.\d+$', '', _expand_acc(a))
            want = {acc_key(a) for a in accs}
            col = (sel['accession'].map(acc_key).isin(want)
                   | sel['accnov'].map(acc_key).isin(want))
            sel = sel[col]
            labels.append(f'genome-list({len(accs)})')

        if o.get('taxids'):
            tids = set()
            for ref in _as_list(o['taxids']):
                try:
                    tids.add(self.taxid_map[self.node_key_of(ref)])
                except (ValueError, KeyError):
                    if str(ref).isdigit():
                        tids.add(str(ref))
                    else:
                        raise
            col = sel['gtdb_taxonomy'].apply(
                lambda lin: any(self.taxid_map.get(
                    LINEAGE_SEP.join(lin.split(LINEAGE_SEP)[:d + 1])) in tids
                    for d in range(DEPTH)))
            sel = sel[col]
            labels.append(f'taxids({len(tids)})')

        if o.get('genoids'):
            accs = {_expand_acc(a) for a in _as_list(o['genoids'])}
            sel = sel[sel['accession'].isin(accs) | sel['accnov'].isin(accs)
                      | sel['genome'].isin(accs)]
            labels.append('genome-ids')

        if o.get('exclude'):
            # Never silently treat an exclusion request as an inclusion (the
            # CLI rejects -e/--exclude before reaching this point; this raise
            # covers direct library use of select()).
            raise NotImplementedError(
                'exclude is not implemented for --source gtdb: select the '
                'genomes to keep with --gtdb-taxa/--genome-list instead.')

        if o.get('min_completeness') is not None:
            sel = sel[sel['completeness'] >= float(o['min_completeness'])]
            labels.append(f'completeness>={o["min_completeness"]}')
        if o.get('max_contamination') is not None:
            sel = sel[sel['contamination'] <= float(o['max_contamination'])]
            labels.append(f'contamination<={o["max_contamination"]}')
        if o.get('mimag'):
            tiers = {'high': ('mimag_high_quality',),
                     'medium': ('mimag_medium_quality',),
                     'low': ('mimag_low_quality',),
                     # 'any' is a plan §7 choice: require that a genome carries
                     # *some* MIMAG tier instead of a specific one.
                     'any': ('mimag_high_quality', 'mimag_medium_quality',
                             'mimag_low_quality')}
            if o['mimag'] not in tiers:
                raise ValueError(f'--mimag {o["mimag"]} is not one of '
                                 'high / medium / low / any.')
            cols = tiers[o['mimag']]
            miss = [c for c in cols if c not in sel]
            if miss:
                raise ValueError(
                    '--mimag ' + o['mimag'] + ' needs the ' +
                    ', '.join(miss) + ' column(s), which this metadata table '
                    'does not carry; keeping every genome would be worse than '
                    'refusing. Not every GTDB release publishes MIMAG quality '
                    'classes.')
            keep = sel[cols[0]] == 't'
            for c in cols[1:]:
                keep = keep | (sel[c] == 't')
            sel = sel[keep]
            labels.append('MIMAG-any-tier' if o['mimag'] == 'any'
                          else f'MIMAG-{o["mimag"]}')
        which = o.get('type_material') or ('ncbi' if o.get('typemater') else None)
        if which:
            before = sel.shape[0]
            if which == 'ncbi':
                col = 'ncbi_type_material_designation'
                if col not in sel:
                    raise ValueError(
                        f'--type-material ncbi needs the {col} column; the '
                        'loaded metadata table does not carry it.')
                keep = ~_designation_is_empty(sel[col])
                sel = sel[keep]
                labels.append(f'NCBI type material ({before} -> {sel.shape[0]} '
                              'genomes)')
            elif which == 'gtdb':
                col = 'gtdb_type_species_of_genus'
                if col not in sel:
                    raise ValueError(
                        f'--type-material gtdb needs the {col} column; the '
                        'loaded metadata table does not carry it, and keeping '
                        'every genome would silently ignore the request.')
                sel = sel[sel[col] == 't']
                labels.append(f'GTDB type species ({before} -> '
                              f'{sel.shape[0]} genomes)')
            else:
                raise ValueError(f'--type-material {which} is not one of '
                                 'gtdb / ncbi.')
        if o.get('complete'):
            # Upstream means "assembly level is Complete Genome or Chromosome".
            # Measured on r220 (缓存测试/reports/R29_quality_class_counts.json):
            # that keeps 8,593 of 113,104 representatives, and the 104,511 it
            # drops have a median GTDB completeness of 94-96 % -- assembly level
            # describes how a sequence was submitted, not whether the genome is
            # complete. Plan §7 R5 replaces the flag with the quality
            # thresholds, so it must not pass silently.
            raise ValueError(
                '--complete is not a GTDB filter: NCBI assembly level would '
                'drop 92.4% of the r220 representatives whose median '
                'completeness is 94-96% (see plan §7 R5). Use '
                '--min-completeness/--max-contamination instead.')

        if not o.get('all_domains', True):
            # `--all-domains no` has no GTDB implementation either; refusing
            # beats pretending the flag did something.
            raise NotImplementedError(
                'all_domains=no is not implemented for --source gtdb; use '
                '--domain bacteria|archaea instead.')
        if o.get('domain'):
            sel = sel[sel['domain'].isin(_as_list(o['domain']))]
            labels.append(f'domain={",".join(_as_list(o["domain"]))}')

        self.selection_labels = labels
        return sel

    def sample(self, df, opts=None):
        """Sample genomes per taxonomic group (plan §6.4 task 2.13).

        Parameters
        ----------
        df : pandas.DataFrame
            selected genomes
        opts : dict, optional
            'rank'/'sample'/'above' and 'rank_sample'

        Returns
        -------
        pandas.DataFrame
            sampled genomes
        """
        o = dict(self.opts)
        o.update(opts or {})

        if o.get('rank_sample'):
            rank = rank_normalize(o['rank_sample'], 'gtdb')
            if rank not in RANK_ORDINAL:
                raise ValueError(
                    f'--rank-sample rank "{o["rank_sample"]}" is not a GTDB '
                    'rank: ' + ', '.join(GTDB_RANKS))
            depth = RANK_ORDINAL[rank] + 1
            key = df['gtdb_taxonomy'].str.split(LINEAGE_SEP).str[:depth] \
                .str.join(LINEAGE_SEP)
            # prefer the highest-quality genome of each group (plan R5)
            df = df.assign(_key=key, _q=df['quality_score'].fillna(-1))
            if 'is_rep' in df:
                df = df.assign(_rep=df['is_rep'].astype(int))
            df = df.sort_values(['_key', '_rep', '_q', 'accession'],
                                ascending=[True, False, False, True])
            df = df.groupby('_key', as_index=False, sort=False).head(1)
            df = df.drop(columns=['_key', '_q', '_rep'])
            self.sample_note = (
                f'one genome per {rank} ({df.shape[0]:,} groups)')
            return df

        n = o.get('sample')
        if n is None:
            return df
        try:
            n = int(n)
        except (TypeError, ValueError):
            raise ValueError(
                f'--sample takes an integer genome count, got {o.get("sample")!r}.'
                ) from None
        if n <= 0:
            return df
        rank = rank_normalize(o.get('rank') or 'species', 'gtdb')
        if rank not in RANK_ORDINAL:
            raise ValueError(
                f'--rank "{o.get("rank")}" is not a GTDB rank: '
                + ', '.join(GTDB_RANKS))
        depth = RANK_ORDINAL[rank] + 1
        lin = df['gtdb_taxonomy'].str.split(LINEAGE_SEP)
        key = lin.str[:depth].str.join(LINEAGE_SEP)
        # a lineage that stops above the sampled rank is not identified to it,
        # exactly like an NCBI genome with a NaN in that rank column: it is left
        # out of the base sample and only --above can pick it up
        key = key.where(lin.str.len() >= depth)
        n = int(n)
        df = df.assign(_key=key)
        base = df.dropna(subset=['_key'])
        if base.shape[0] < df.shape[0]:
            n_short = df.shape[0] - base.shape[0]
            print(f'  {n_short:,} genome(s) are not identified to {rank} '
                  'and were left out of the sampling'
                  + ('' if o.get('above') else
                     ' (add --above to take them at a coarser rank)'))
        out = []
        for k, grp in base.groupby('_key', sort=True):
            grp = grp.sort_values('quality_score', ascending=False,
                                  na_position='last')
            out.append(grp.head(n))
        res = pd.concat(out) if out else df.head(0)
        res = res.drop(columns=['_key'])
        self.sample_note = f'{n} genome(s) per {rank}'
        if not o.get('above'):
            return res
        # --above (plan §6.4 task 2.5) tops every coarser group up to n genomes,
        # in the same order as the NCBI route: from the rank just above the
        # sampled one down to phylum. The domain level is not sampled, because
        # upstream's rank list stops at phylum too, and topping up to "n per
        # domain" would silently drop whole phyla.
        chosen = set(res['genome'])
        more = []
        parts = []
        for depth in range(RANK_ORDINAL[rank] - 1, 0, -1):
            sub = df.assign(_k=lin.str[:depth + 1].str.join(LINEAGE_SEP))
            done = sub[sub['genome'].isin(chosen)].groupby('_k').size()
            before = len(chosen)
            for k, grp in sub.groupby('_k', sort=True):
                left = n - int(done.get(k, 0))
                if left <= 0:
                    continue
                cand = grp[~grp['genome'].isin(chosen)].sort_values(
                    'quality_score', ascending=False, na_position='last')
                for g in cand.head(left)['genome']:
                    chosen.add(g)
                    more.append(g)
            parts.append(f'+{len(chosen) - before} at {GTDB_RANKS[depth]}')
        self.sample_note = f'{n} genome(s) per {rank} and above (' \
            + ', '.join(parts) + ')'
        if not more:
            return res
        extra = df[df['genome'].isin(more)].drop(columns=['_key'])
        return pd.concat([res, extra]).drop_duplicates(subset=['genome'])

    def protein_sources(self, df):
        """Yield the sequence location of every selected genome.

        Parameters
        ----------
        df : pandas.DataFrame
            selected genomes

        Returns
        -------
        list of GenomeFetch
            fetch plans; kind is 'local' for on-disk files
        """
        if df is None or df.empty:
            return []
        idx = stem_index(self.info, domain=None)
        if self.taxid_map is None:
            self.load_taxonomy()
        tid_of = {row.accession: self.taxid_map.get(row.gtdb_taxonomy)
                  for row in df.itertuples()}
        accs = list(df['accession'])
        for acc in accs:
            hit = idx.get(acc) or idx.get(acc.split('.')[0])
            if hit:
                yield GenomeFetch(acc, None, tid_of.get(acc), path=hit['path'],
                                  kind='local')
            else:
                yield GenomeFetch(acc, None, tid_of.get(acc), path=None,
                                  kind='missing')

EMPTY_DESIGNATIONS = ('', 'na', 'none', '-', 'nan')


def _designation_is_empty(series):
    """True where an NCBI designation column carries no information.

    GTDB copies leave 'na' in `ncbi_type_material_designation` for most genomes
    (measured on the local r220 metadata: 385,255 of the first 400,000 rows),
    so 'has a designation' cannot be tested with a plain truthiness check.
    """
    return series.astype(str).str.strip().str.lower().isin(EMPTY_DESIGNATIONS)


def _as_list(val):
    """Normalize a comma-separated or list argument into a list."""
    if val is None:
        return []
    if isinstance(val, str):
        return [x.strip() for x in val.split(',') if x.strip()]
    if isinstance(val, (list, tuple, set)):
        out = []
        for x in val:
            out.extend(_as_list(x))
        return out
    return [str(val)]


def _expand_acc(acc):
    """Accept 'RS_GCF_000123456.1', 'GCF_000123456.1' or 'GCF_000123456'."""
    acc = str(acc).strip()
    if acc.startswith(('RS_', 'GB_')):
        return acc
    if acc.startswith(('GCF_', 'GCA_')):
        return 'RS_' + acc if acc.startswith('GCF_') else 'GB_' + acc
    return acc


# ---------------------------------------------------------------------------
# Streaming sequence ingestion (plan §6.4 tasks 2.7/2.8, fixes RK2)
# ---------------------------------------------------------------------------

DIAMOND_TAXONMAP_HEADER = 'accession\taccession.version\ttaxid\tgi\n'
BLAST_TAXIDMAP_HEADER = 'seqid\ttaxid\n'


def extract_genomes_stream(fetches, out_faa, out_taxonmap, layout,
                           blast_taxidmap=None, genome_map=None,
                           on_progress=None, limit=None, resolve=None):
    """Concatenate per-genome FASTA files into db.faa, streaming.

    Parameters
    ----------
    fetches : iterable of GenomeFetch
        plan produced by GtdbProvider.protein_sources()
    out_faa : str
        destination db.faa path (inside the workspace)
    out_taxonmap : str
        destination DIAMOND taxonmap path
    layout : str
        layout label, which selects the sequence-Id anchor
    blast_taxidmap : str, optional
        also write a makeblastdb '-taxid_map' table here
    genome_map : str, optional
        also write a protein-to-genome table here
    on_progress : callable, optional
        called as on_progress(done, total) after each genome
    limit : int, optional
        stop after this many genomes (used by the smoke tests)
    resolve : callable, optional
        given the (seqid, taxId) candidates of one shared protein, return the
        pairs to keep (`--hybrid-lca-policy`); without it the first taxId met
        during extraction silently wins

    Returns
    -------
    dict
        counters plus the stem2taxid, prefix2genome and per-genome tables

    Notes
    ----
    Nothing is kept in memory except the per-genome counters and the route C
    tables, plus the dedup set `seen_ids`, which is bounded by the number of
    *proteins*: at ~200 B/id the full 333M-protein representative set projects
    to tens of GB of RAM (the project's own stream benchmark extrapolates
    98.9 GB), so a full-release build on a 128 GB machine is not established.
    The sqlite-dedup pattern of the NCBI stream path (database.py
    extract_genomes_stream_ncbi) is the drop-in fix if that build is ever
    attempted; until then treat RK2 as "db.faa never sits in memory", not
    "memory is independent of protein count".
    """
    out_faa = guard_write(out_faa)
    out_taxonmap = guard_write(out_taxonmap)
    makedirs(dirname(out_faa), exist_ok=True)
    makedirs(dirname(out_taxonmap), exist_ok=True)

    stem2taxid, prefix2genome = {}, {}
    per_genome = {}
    by_domain = {}
    seen_ids = set()
    dups = {}
    stats = {'genomes_planned': 0, 'genomes_done': 0, 'genomes_missing': 0,
             'genomes_unreadable': 0,
             'proteins': 0, 'proteins_dup': 0, 'proteins_no_taxid': 0,
             'proteins_shared': 0, 'proteins_shared_reassigned': 0,
             'proteins_shared_dropped': 0,
             'residues': 0, 'bytes': 0, 'max_id_len': 0, 'lines': 0}
    unresolved = []

    fetches = list(fetches)
    if limit is not None:
        fetches = fetches[:limit]
    total = len(fetches)
    stats['genomes_planned'] = total

    # write db.faa under a .partial name and rename on success: an interrupt
    # then leaves a stale-but-complete db.faa plus a partial, never a truncated
    # file that a later resume would mistake for the real one
    partial_faa = out_faa + '.partial'
    with open(partial_faa, 'w') as faa, open(out_taxonmap, 'w') as tmap, \
            (open(blast_taxidmap, 'w') if blast_taxidmap
             else io.StringIO()) as bmap, \
            (open(genome_map, 'w') if genome_map else io.StringIO()) as gmap:
        tmap.write(DIAMOND_TAXONMAP_HEADER)
        if blast_taxidmap:
            bmap.write(BLAST_TAXIDMAP_HEADER)
        for fetch in fetches:
            if not fetch.path or not isfile(fetch.path):
                stats['genomes_missing'] += 1
                unresolved.append(fetch.accession)
                continue
            stem = strip_faa_ext(basename(fetch.path))
            if fetch.taxid:
                stem2taxid[stem] = fetch.taxid
            n_prot = n_res = 0
            seqid = None
            seq = []
            tid = fetch.taxid

            def emit(sid, s, t):
                """Write one record and its lookup rows."""
                nonlocal n_prot, n_res
                if not s or not sid:
                    return
                if sid in seen_ids:
                    stats['proteins_dup'] += 1
                    dups.setdefault(sid, []).append(t)
                    return
                if not t:
                    stats['proteins_no_taxid'] += 1
                    return
                seen_ids.add(sid)
                stats['proteins'] += 1
                stats['residues'] += len(s)
                stats['max_id_len'] = max(stats['max_id_len'], len(sid))
                n_prot += 1
                n_res += len(s)
                rec = f'>{sid}\n{s}\n'
                faa.write(rec)
                stats['bytes'] += len(rec)
                stats['lines'] += 2
                # both key columns hold the identical full Id (plan 附录 A.3)
                tmap.write(f'{sid}\t{sid}\t{t}\tna\n')
                if blast_taxidmap:
                    bmap.write(f'{sid}\t{t}\n')
                if genome_map:
                    gmap.write(f'{sid}\t{fetch.accession}\n')

            try:
                with open_reader(fetch.path) as f:
                    for line in f:
                        if line.startswith('>'):
                            if seqid is not None:
                                emit(seqid, ''.join(seq), tid)
                            raw = line[1:].rstrip('\r\n')
                            seqid = raw.split('\t')[0].split(' ')[0]
                            seq = []
                            tid = taxid_of_seqid(seqid, layout, stem2taxid) \
                                or fetch.taxid
                            parsed = parse_seqid(seqid, layout)
                            if parsed and parsed.get('contig'):
                                prefix2genome.setdefault(
                                    parsed['contig'], fetch.accession)
                        elif seqid is not None:
                            # newline first, then the stop codon: rstrip('*')
                            # on a line that still ends in '\n' strips
                            # nothing (the original order left '*' on ~99.6%
                            # of proteins)
                            seq.append(line.rstrip('\r\n').rstrip('*'))
                    if seqid is not None:
                        emit(seqid, ''.join(seq), tid)
            except (OSError, EOFError, UnicodeError) as exc:
                # a corrupt archive member used to abort hours into the build
                # with half of db.faa written and no name for the culprit; the
                # genome is skipped, counted, and the completeness gate
                # refuses to publish the build unless --force
                print(f'  WARNING: unreadable sequence file for '
                      f'{fetch.accession} ({basename(fetch.path)}): {exc}')
                stats['genomes_unreadable'] += 1
                unresolved.append(fetch.accession)
                seqid, seq = None, []
                continue
            per_genome[fetch.accession] = (n_prot, n_res)
            dom = (parse_stem(stem).get('domain', '')
                   if is_renamed_layout(layout) else '')
            if dom:
                d = by_domain.setdefault(dom.lower(), {'genomes': 0,
                                                       'proteins': 0})
                d['genomes'] += 1
                d['proteins'] += n_prot
            stats['genomes_done'] += 1
            if on_progress:
                on_progress(stats['genomes_done'], total)
    os.replace(partial_faa, out_faa)
    if dups and resolve:
        # the first taxId met during streaming was written already, so the
        # policy is applied in one substitution pass over the maps: memory
        # stays bounded by the number of *shared* proteins, not of proteins
        maps = [(out_taxonmap, 2)]
        if blast_taxidmap:
            maps.append((blast_taxidmap, 1))
        for path, taxid_col in maps:
            tmp = guard_write(path + '.resolved')
            with open(path) as src, open(tmp, 'w') as dst:
                for line in src:
                    if line.startswith('#') or not line.strip():
                        dst.write(line)
                        continue
                    cells = line.rstrip('\n').split('\t')
                    sid = cells[0]
                    if sid not in dups:
                        dst.write(line)
                        continue
                    first = cells[taxid_col]
                    chosen = resolve([(sid, first)]
                                     + [(sid, t) for t in dups[sid]])
                    stats['proteins_shared'] += 1
                    if not chosen:
                        stats['proteins_shared_dropped'] += 1
                        keep = 'na'
                    else:
                        keep = str(chosen[0][1])
                        if keep != str(first):
                            stats['proteins_shared_reassigned'] += 1
                    cells[taxid_col] = keep
                    dst.write('\t'.join(cells) + '\n')
            os.replace(tmp, path)
    stats['shared_proteins'] = len(dups)
    stats['unresolved'] = unresolved
    stats['stem2taxid'] = stem2taxid
    stats['prefix2genome'] = prefix2genome
    stats['per_genome'] = per_genome
    stats['by_domain'] = by_domain
    return stats


def reconcile_protein_counts(df, stats):
    """Compare the proteins the metadata promises with those actually extracted.

    Plan §10.2 UT-2.x asks for the release-wide protein total to be a checked
    quantity rather than a number quoted from the paper. A build that silently
    loses a genome (skipped file, unreadable record) then shows up in the
    manifest instead of surfacing weeks later as unexplained HGT counts.

    Parameters
    ----------
    df : pandas.DataFrame
        the genome table of this build (selected subset, or all rows)
    stats : dict
        extraction counters, 'proteins' being the records written to db.faa

    Returns
    -------
    tuple
        (planned, extracted, verdict) with verdict 'ok' when they agree
    """
    planned = int(df['proteins'].sum()) if 'proteins' in df.columns else 0
    extracted = int(stats.get('proteins') or 0)
    if extracted == planned:
        verdict = 'ok'
    elif not extracted:
        verdict = 'not extracted yet'
    else:
        verdict = f'{planned - extracted:,} short'
    return planned, extracted, verdict


def proteins_reconciliation(df, stats):
    """The manifest block comparing promised with extracted protein counts.

    Returned as a dict on purpose: the first version unpacked
    `reconcile_protein_counts()` positionally and bound `proteins_planned` to
    the *extracted* count, which made a 4,896-genome build look both complete
    and 1.94 billion proteins short in the same file.
    """
    planned, extracted, verdict = reconcile_protein_counts(df, stats)
    return {'proteins_planned': planned,
            'proteins_extracted_vs_planned': verdict}


def write_meta_json(path, provider, args, stats, taxonmap_sha=None,
                   diamond_version=None, started=None, finished=None,
                   config_snapshot=None, selection=None):
    """Write the reproducibility manifest (plan §5.6).

    Parameters
    ----------
    path : str
        destination meta.json
    provider : GtdbProvider
        prepared provider holding layout and provenance information
    args : dict
        effective build options
    stats : dict
        extraction counters
    taxonmap_sha : str, optional
        sha256 of the written taxonmap
    diamond_version : str, optional
        'diamond --version' output
    started, finished : str, optional
        ISO timestamps
    config_snapshot : dict, optional
        loaded config.yml

    Returns
    -------
    dict
        the manifest that was written
    """
    path = guard_write(path, must_not_contain=provider.gtdb_dir)
    info = provider.info
    df = provider.df_all
    reps = df[df['is_rep']]
    # the subset actually requested (--gtdb-taxa, --rank-sample, ...) is what
    # the extraction counters describe, so reconcile against that. Note that
    # provider.df is the *loaded* table (load_genomes() caches every genome of
    # the release there), not the selection: reconciling against it reported a
    # 4,896-genome build as "1,944,720,987 proteins short", which is the whole
    # release minus the build. The caller must pass the selection explicitly.
    sel = selection
    if sel is None:
        sel = getattr(provider, 'df_selected', None)
    if sel is None:
        sel = df
    recon = proteins_reconciliation(sel, stats)
    by_dom = {}
    for dom, grp in reps.groupby('domain'):
        by_dom[dom] = {'genomes': int(grp.shape[0]),
                       'proteins': int(grp['proteins'].sum())}
    taxid_map = provider.taxid_map or {}
    manifest = {
        'tool': f'neoHGT {tool_version()}',
        'source': provider.name,
        'gtdb_release': info['release'],
        'seq_source': args.get('seq_source', 'local'),
        'local_layout': info['layout'],
        'id_convention': info['id_convention'],
        # provenance: whether this run reused an existing db.faa (IT-6)
        'extraction_resumed': bool(stats.get('resumed', False)),
        'input_provenance': {
            'gtdb_dir': info['gtdb_dir'],
            'source_writable_by_user': os.access(info['gtdb_dir'], os.W_OK),
            'read_only_policy': ('the dataset is treated as read-only by '
                                 'guard_write plus the post-build mtime audit '
                                 '(IT-7); filesystem permissions are not the '
                                 'enforcement point'),
            'taxonomy_files': [realpath(p) for p in info['taxonomy']],
            'metadata_files': [realpath(p) for p in info['metadata']],
            'sequence_dirs': {d: v['dir'] for d, v in info['domains'].items()},
            'third_party_scripts': [realpath(p)
                                    for p in info['provenance_scripts']],
            'absent_official_files': info['absent_official'],
            'serial_anchor_coverage': info['serial_anchor_coverage'],
            'input_sha256': stats.get('input_sha256', {}),
        },
        'diamond_version': diamond_version,
        'taxid_map_sha256': taxonmap_sha,
        'counts': {
            'genomes': stats.get('genomes_done', 0),
            'genomes_planned': stats.get('genomes_planned', 0),
            'genomes_missing': stats.get('genomes_missing', 0),
            'proteins': stats.get('proteins', 0),
            'proteins_planned': recon['proteins_planned'],
            'proteins_extracted_vs_planned':
                recon['proteins_extracted_vs_planned'],
            'proteins_duplicated': stats.get('proteins_dup', 0),
            'proteins_without_taxid': stats.get('proteins_no_taxid', 0),
            'taxa': len(taxdump_of(provider)),
            'db_faa_bytes': stats.get('bytes', 0),
            'max_seqid_len': stats.get('max_id_len', 0),
            'by_domain_selected': stats.get('by_domain', {}),
            'representatives_in_dataset': {
                'total': int(reps.shape[0]),
                'proteins': int(reps['proteins'].sum()),
                'by_domain': by_dom},
        },
        'taxonomy': {
            'id_base': provider.id_base,
            'rank_counts': (provider.stats or {}).get('by_rank', {}),
            'n_superseded': (provider.stats or {}).get('n_superseded', 0),
        },
        'selfcheck': (provider.report or {}).get('status', 'skipped'),
        'selfcheck_items': [
            {k: c[k] for k in ('name', 'expected', 'observed', 'status')}
            for c in (provider.report or {}).get('checks', [])],
        'selection': {'labels': getattr(provider, 'selection_labels', []),
                      'sampling': getattr(provider, 'sample_note', None),
                      'options': {k: str(v) for k, v in args.items()
                                  if v is not None}},
        'build_started': started or str(timestamp()),
        'build_finished': finished or str(timestamp()),
        'config_snapshot': config_snapshot or {},
    }
    if provider.name == 'hybrid':
        st = provider.stats or {}
        manifest['taxonomy']['hybrid'] = {
            'ncbi_taxdump': getattr(provider, 'ncbi_dir', None),
            'keep_ranks': list(provider.keep_ranks),
            'lca_policy': st.get('lca_policy'),
            'n_ncbi_nodes': st.get('n_ncbi_nodes', 0),
            'n_nodes_merged': st.get('n_nodes_merged', 0),
            'ranks_flattened_for_diamond':
                st.get('ranks_flattened_for_diamond', {}),
        }
    with open(path, 'w') as f:
        json.dump(manifest, f, indent=2, sort_keys=True, default=str)
    return manifest


def taxdump_of(provider):
    """Get the provider's taxdump ({} when not built yet)."""
    return provider.taxdump or {}


def public_stats(stats):
    """Strip internal structures from conversion stats before serializing.

    Parameters
    ----------
    stats : dict
        conversion statistics from build_taxdump_from_gtdb()

    Returns
    -------
    dict
        JSON-safe subset ('tree' holds 146k node records and 'superseded' may
        hold many Ids; both would bloat meta.json by tens of MB)
    """
    res = {k: v for k, v in (stats or {}).items()
           if k not in ('tree', 'superseded')}
    res['superseded_nodes'] = sorted((stats or {}).get('superseded', {}).items())
    return res


# ---------------------------------------------------------------------------
# NCBI provider (plan 附录 B: a thin wrapper over the existing logic)
# ---------------------------------------------------------------------------


class NcbiProvider(SourceProvider):
    """The plan's Appendix B adapter around the untouched NCBI pipeline.

    `--source ncbi` still runs upstream code; this class exists so a caller that
    wants one interface across the three sources can get it -- same method set
    as `GtdbProvider`, behaviour deliberately unchanged: no ranking, no
    sampling, and no sequence list (upstream downloads per genome during the
    build, so there is nothing to hand out beforehand).

    `wrap` is the object the pipeline already uses (a `Database` or `Search`
    instance); when given, `load_taxonomy`/`load_genomes` defer to it instead of
    reading paths, which is what keeps upstream behaviour byte-for-byte.
    """

    name = 'ncbi'
    ranks = NCBI_RANKS
    default_taxonomy_opts = {'capital': True, 'latin': False, 'block': []}

    def __init__(self, taxdump_dir=None, df=None, opts=None, wrap=None):
        super().__init__(opts)
        self.taxdump_dir = taxdump_dir
        self.df = df
        self.wrap = wrap

    def prepare(self):
        self.info = {'release': 'ncbi', 'layout': 'ncbi',
                     'gtdb_dir': self.taxdump_dir or ''}
        return self.info

    def load_taxonomy(self, trim=None):
        if self.taxdump is not None:
            return self.taxdump, {}, {'n_nodes': len(self.taxdump)}
        if self.wrap is not None:
            self.taxdump = getattr(self.wrap, 'taxdump', None)
            return self.taxdump, {}, {'n_nodes': len(self.taxdump or {})}
        if not self.taxdump_dir:
            raise ValueError(
                'NcbiProvider needs a taxdump directory, a genome table, or '
                'the upstream object to wrap.')
        self.taxdump = read_taxdump(self.taxdump_dir)
        return self.taxdump, {}, {'n_nodes': len(self.taxdump)}

    def load_genomes(self):
        if self.wrap is not None:
            return getattr(self.wrap, 'df', self.df)
        return self.df

    def select(self, df=None, opts=None):
        return df if df is not None else self.load_genomes()

    def sample(self, df, opts=None):
        return df

    def protein_sources(self, df):
        return []


# ---------------------------------------------------------------------------
# Hybrid provider (plan §5.7)
# ---------------------------------------------------------------------------

HYBRID_LCA_POLICIES = ('root', 'domain', 'refuse')


class HybridProvider(SourceProvider):
    """GTDB sequences under an NCBI-compatible root.

    Semantics, stated plainly (ADR-7): a hybrid database contains **GTDB
    proteins only**. The NCBI taxdump contributes a small classification
    skeleton -- the cellular root plus the `ncbi_superkingdoms` anchors and
    their ancestors -- so the tree has a defined root in NCBI coordinates. No
    NCBI protein is ingested: eukaryote/virus genomes are **not** searchable
    here, and searching a hybrid database returns the same hits as the
    underlying GTDB database.

    Parameters
    ----------
    gtdb_provider : GtdbProvider
        prepared GTDB provider
    ncbi_taxdump_dir : str
        directory holding nodes.dmp/names.dmp for the NCBI part
    keep_ranks : iterable of str, optional
        NCBI ranks to graft; the cellular root is always kept
    lca_policy : str, optional
        'root', 'domain' or 'refuse' (plan §5.7 item 3)
    """

    name = 'hybrid'
    ranks = GTDB_RANKS
    rank_prefixes = GTDB_RANK_PREFIXES

    def __init__(self, gtdb_provider, ncbi_taxdump_dir, keep_ranks=('superkingdom',),
                 lca_policy='root', opts=None, ncbi_superkingdoms=('Eukaryota',)):
        super().__init__(opts)
        self.gtdb = gtdb_provider
        self.ncbi_dir = ncbi_taxdump_dir
        self.keep_ranks = tuple(keep_ranks)
        if lca_policy not in HYBRID_LCA_POLICIES:
            raise ValueError(
                f'Unknown --hybrid-lca-policy "{lca_policy}": '
                + ', '.join(HYBRID_LCA_POLICIES))
        self.lca_policy = lca_policy
        self.ncbi_superkingdoms = tuple(ncbi_superkingdoms)
        self.ncbi_taxdump = None
        self.x_source_proteins = []

    def prepare(self):
        self.info = self.gtdb.info
        # the genome tables the build reads come from the wrapped provider; only
        # the taxonomy step is overridden by the merge below
        self.df_all = self.gtdb.df_all
        self.lineages = self.gtdb.lineages
        return self.info

    def __getattr__(self, name):
        # the hybrid wrapper replaces only the taxonomy step; everything else the
        # build reads off a provider (selection labels, id base, notes, ...) is
        # the wrapped GTDB provider's
        inner = self.__dict__.get('gtdb')
        if inner is None:
            raise AttributeError(name)
        return getattr(inner, name)

    def load_taxonomy(self, trim=None):
        """Merge the GTDB subtree under the NCBI cellular root."""
        inner, taxid_map, stats = self.gtdb.load_taxonomy(trim=trim)
        # copy before merging: the graft re-parents the GTDB domains, and
        # writing through to the wrapped provider would leave a "GTDB-only"
        # tree that silently contains NCBI nodes
        taxdump = {tid: dict(rec) for tid, rec in inner.items()}
        ncbi = read_taxdump(self.ncbi_dir)
        self.ncbi_taxdump = ncbi
        cr = ncbi.get(CELLULAR_ROOT_TID)
        if cr is None:
            raise ValueError(
                f'NCBI taxdump at {self.ncbi_dir} has no cellular root '
                f'({CELLULAR_ROOT_TID}); hybrid mode needs a complete NCBI '
                'taxdump (it is not part of the local GTDB copy, plan §3.4).')
        # keep the NCBI path from cellular root down to the chosen kingdoms
        # (current NCBI dumps spell the top rank 'domain', older ones
        # 'superkingdom', so compare through the alias table)
        keep_ranks = {rank_normalize(r, 'ncbi') for r in self.keep_ranks}
        keep = {CELLULAR_ROOT_TID}
        for tid, rec in ncbi.items():
            if rank_normalize(rec.get('rank') or '', 'ncbi') in keep_ranks \
                    and _in_ncbi_lineage(tid, self.ncbi_superkingdoms, ncbi):
                for anc in _ncbi_lineage_tids(tid, ncbi):
                    keep.add(anc)
        for tid in sorted(keep, key=int):
            rec = ncbi[tid]
            parent = rec['parent']
            if tid == CELLULAR_ROOT_TID:
                # the cellular root keeps its own NCBI parent (taxId 1); making
                # it self-referential would create a second root and break
                # find_lca
                target = parent if parent in ncbi else ROOT_TID
            elif parent in keep:
                target = parent
            else:
                target = CELLULAR_ROOT_TID
            if tid in taxdump and tid != ROOT_TID:
                # taxId 1 is the shared root by construction; anything else
                # colliding means two different taxa would wear one id
                raise ValueError(
                    f'NCBI taxId {tid} ({rec.get("name")}) collides with a '
                    'GTDB node of the same id: the synthetic range '
                    f'(id_base={getattr(self.gtdb, "id_base", "?")}) is not '
                    'safe against this NCBI taxdump. Raise --gtdb-taxid-base '
                    'instead of merging two taxa into one id.')
            taxdump[tid] = {'parent': target, 'rank': rec['rank'],
                            'name': rec['name']}
        # re-attach GTDB domains under the cellular root
        for key, tid in taxid_map.items():
            if gtdb_is_domain(key):
                taxdump[tid]['parent'] = CELLULAR_ROOT_TID
        if find_lca(list(taxdump), taxdump) != ROOT_TID:
            raise AssertionError(
                'Hybrid tree is not single-rooted; check that the NCBI '
                'taxdump contains taxId 1 and 131567.')
        stats['n_nodes_merged'] = len(taxdump)
        stats['n_ncbi_nodes'] = len(keep)
        stats['lca_policy'] = self.lca_policy
        # report the NCBI ranks that nodes.dmp cannot carry verbatim, so a
        # hybrid database records exactly how much of the graft was flattened
        collapsed = illegal_dmp_ranks(taxdump)
        stats['ranks_flattened_for_diamond'] = {
            r: len(t) for r, t in collapsed.items()}
        self.taxdump, self.taxid_map, self.stats = taxdump, taxid_map, stats
        return taxdump, taxid_map, stats

    def load_genomes(self):
        return self.gtdb.load_genomes()

    def select(self, df=None, opts=None):
        # the build sets `provider.opts` on the wrapper, so forward it: without
        # this a hybrid run silently keeps every genome (--genome-list ignored)
        return self.gtdb.select(df, opts if opts is not None else self.opts)

    def sample(self, df, opts=None):
        return self.gtdb.sample(df, opts if opts is not None else self.opts)

    def protein_sources(self, df):
        return self.gtdb.protein_sources(df)

    def domain_ancestor(self, tid):
        """The domain-level node above this taxId in the merged tree.

        A GTDB domain hangs directly off the cellular root, and an NCBI
        superkingdom is spelled 'domain' in current dumps and 'superkingdom'
        in older ones, so both tests are needed.
        """
        cur, seen = str(tid), set()
        while cur in self.taxdump and cur not in seen:
            seen.add(cur)
            rec = self.taxdump[cur]
            par = str(rec.get('parent') or '')
            if par == cur or not par:
                return cur
            if par == CELLULAR_ROOT_TID or par == ROOT_TID:
                return cur
            prec = self.taxdump.get(par)
            if prec and rank_normalize(prec.get('rank') or '', 'ncbi') in (
                    'domain', 'superkingdom'):
                return par
            cur = par
        return None

    def resolve_shared_protein(self, pairs):
        """Apply `--hybrid-lca-policy` to proteins shared between sources.

        Parameters
        ----------
        pairs : iterable of tuple
            (seqid, taxId) candidates for one and the same protein, in the
            order the extraction met them

        Returns
        -------
        list of tuple
            the pairs to keep; 'root' keeps every candidate, 'domain' labels
            the protein with the donor's own domain node and marks it
            ``x_source``, 'refuse' drops non-GTDB candidates
        """
        pairs = list(pairs)
        if self.lca_policy == 'root':
            return pairs
        if self.lca_policy == 'domain':
            by_seq = {}
            for sid, tid in pairs:
                by_seq.setdefault(sid, []).append(tid)
            out = []
            for sid, tids in by_seq.items():
                doms = {self.domain_ancestor(t) for t in tids}
                doms.discard(None)
                if len(doms) > 1 or len(set(tids)) > 1:
                    # plan §5.7 item 3: cross-source proteins carry an
                    # explicit mark instead of pretending to be one taxon
                    self.x_source_proteins.append(sid)
                pick = tids[0]
                if len(doms) > 1:
                    gtdb_first = [t for t in tids
                                  if int(t) >= int(self.gtdb.id_base)]
                    pick = (gtdb_first or tids)[0]
                d = self.domain_ancestor(pick)
                out.append((sid, d or ROOT_TID))
            return out
        # 'refuse': drop any protein whose genome is not GTDB-owned
        return [(sid, tid) for sid, tid in pairs
                if int(tid) >= int(self.gtdb.id_base)]


def gtdb_is_domain(key):
    """Check whether a node path is a GTDB domain node."""
    return key.count(LINEAGE_SEP) == 0


def _ncbi_lineage_tids(tid, ncbi):
    out, cid = [], tid
    seen = set()
    while cid in ncbi and cid not in seen:
        seen.add(cid)
        out.append(cid)
        pid = ncbi[cid]['parent']
        if pid == cid:
            break
        cid = pid
    return out


def _in_ncbi_lineage(tid, names, ncbi):
    return any(ncbi.get(x, {}).get('name') in names
               for x in _ncbi_lineage_tids(tid, ncbi))


# ---------------------------------------------------------------------------
# Reporting helpers
# ---------------------------------------------------------------------------

def load_ncbi_map(db_dir):
    """Read gtdb/ncbi_map.tsv.gz from a database directory."""
    path = join(db_dir, 'gtdb', 'ncbi_map.tsv.gz')
    res = {}
    if not isfile(path):
        return res
    with open_reader(path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('gtdb_taxid'):
                continue
            x = line.rstrip('\r\n').split('\t')
            if len(x) >= 6:
                res[x[0]] = {'label': x[1], 'rank': x[2], 'n': x[3],
                             'conflict': x[4], 'ncbi_taxids': x[5]}
    return res


def load_taxon_synonyms(path):
    """Read a user-supplied table of deprecated -> current GTDB names.

    Parameters
    ----------
    path : str or None
        TSV with `alias<TAB>current` (a third column of notes is ignored).
        `current` may be a prefixed label (`p__Pseudomonadota`) or a taxId.
        Lines starting with '#' and blank lines are skipped. None/'' returns {}.

    Returns
    -------
    dict
        alias -> current reference

    Notes
    ----
    GTDB renames higher taxa between releases (e.g. `p__Proteobacteria` ->
    `p__Pseudomonadota`) but the official `synonyms.*.tsv` is not part of this
    local copy (plan §3.4), so the alias table has to come from the user. This is
    the mitigation for risk RK11: filters written against the old name keep
    working, and analyze prints a note when one is rewritten.
    """
    if not path:
        return {}
    if not isfile(path):
        raise ValueError(f'--taxon-synonyms is not a file: {path}.')
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            x = line.rstrip('\r\n').split('\t')
            if len(x) < 2 or not x[1].strip():
                raise ValueError(
                    f'--taxon-synonyms line needs "alias<TAB>current": '
                    f'{line.rstrip()}')
            out[x[0].strip()] = x[1].strip()
    return out


def map_report(db_dir, in_list, out, to='ncbi'):
    """Map donor taxIds in an HGT list to another classification system.

    Parameters
    ----------
    db_dir : str
        database directory holding gtdb/ncbi_map.tsv.gz
    in_list : str
        HGT list (one taxId or label per line, as written by analyze)
    out : str
        destination TSV
    to : str, optional
        target system ('ncbi' is the only one implemented in v1)

    Returns
    -------
    dict
        summary counters
    """
    if to != 'ncbi':
        raise ValueError(f'Unsupported mapping target "{to}"; use "ncbi".')
    out = guard_write(out)
    makedirs(dirname(out), exist_ok=True)
    ncbi = load_ncbi_map(db_dir)
    labels = load_labels(db_dir)
    stats = {'rows': 0, 'mapped': 0, 'unmapped': 0}
    with open(out, 'w') as f:
        f.write('donor_taxid\tdonor_label\tdonor_rank\tncbi_taxids\t'
                'nconflict\tstatus\n')
        with open(in_list) as g:
            for line in g:
                ref = line.split('\t')[0].strip()
                if not ref or ref.startswith('#'):
                    continue
                tid = resolve_label(ref, labels)
                stats['rows'] += 1
                rec = ncbi.get(tid) if tid else None
                if rec and rec['ncbi_taxids'] != 'na':
                    stats['mapped'] += 1
                    status = 'ok'
                else:
                    stats['unmapped'] += 1
                    status = 'unmapped'
                f.write('\t'.join((
                    tid or ref, (labels.get(tid) or rec['label'])
                    if (tid and (labels.get(tid) or rec)) else 'na',
                    rec['rank'] if rec else 'na',
                    rec['ncbi_taxids'] if rec else 'na',
                    rec['conflict'] if rec else '0', status)) + '\n')
    return stats


# ---------------------------------------------------------------------------
# Command line (plan §6.3: python -m neoHGT.gtdb <action>)
# ---------------------------------------------------------------------------

description = ('GTDB data source tools: probe a local dataset, build a '
               'taxdump, or map donors back to NCBI names')

SUBCOMMANDS = ('probe', 'taxdump', 'map-report', 'gtdbtk', 'reconcile',
               'tree-dist')

arguments = [
    'basic',
    ['--action',   'action to perform (legacy flat form; the full interface '
                   'is "neoHGT gtdb <subcommand> ..." with one of '
                   + ', '.join(SUBCOMMANDS) + ')',
     {'choices': ['probe', 'taxdump', 'map-report'], 'required': True}],
    ['--dir',      'GTDB dataset root (read-only)'],
    ['-o|--output', 'output file or directory'],
    ['--db',       'database directory (for map-report)'],
    ['--in',       'input HGT list (for map-report)'],
    ['--to',       'mapping target for map-report', {'default': 'ncbi'}],
    ['--gtdb-release', 'GTDB release label (inferred from file names if '
                      'omitted)'],
    ['--local-layout', 'force a layout instead of detecting it',
     {'choices': list(LOCAL_LAYOUTS) + ['auto'], 'default': 'auto'}],
    ['--gtdb-taxid-base', 'start of the synthetic taxId range',
     {'type': int, 'default': DEFAULT_ID_BASE}],
    ['--gtdb-id-manifest', 'previous taxid.map.gz to keep IDs stable'],
    ['--gtdb-metadata', 'explicit GTDB metadata TSV path(s) (probe/taxdump)'],
    ['--gtdb-taxonomy', 'explicit GTDB taxonomy TSV path(s) '
                        '(probe/taxdump)'],
    ['--no-selfcheck', 'skip the input self-consistency checks (not '
                       'recommended)', {'action': 'store_true'}],
    ['--deep', 'run the per-accession alignment checks too',
     {'action': 'store_true'}],
    ['--threads', 'number of threads', {'type': int}],
    ['--tmpdir',  'temporary directory (defaults to workspace.tmp_dir)'],
]


def _split_paths(value):
    """Split a comma-delimited path option into a list (None stays None)."""
    if not value:
        return None
    return [x.strip() for x in str(value).split(',') if x.strip()]


class Gtdb(object):
    """CLI module for GTDB utilities (usable as 'neoHGT gtdb ...')."""

    def __init__(self):
        self.arguments = arguments
        self.description = description

    def __call__(self, args):
        print(f'GTDB tool started at {timestamp()}.')
        self.cfg = load_configs_safe()
        action = args.action
        if action == 'probe':
            return cmd_probe(args.dir, args.output, args.gtdb_release,
                             args.local_layout, not args.no_selfcheck,
                             args.deep,
                             metadata=_split_paths(getattr(args,
                                                           'gtdb_metadata',
                                                           None)),
                             taxonomy=_split_paths(getattr(args,
                                                           'gtdb_taxonomy',
                                                           None)))
        if action == 'taxdump':
            return cmd_taxdump(args.dir, args.output, args.gtdb_release,
                               args.local_layout, args.gtdb_taxid_base,
                               args.gtdb_id_manifest, not args.no_selfcheck,
                               metadata=_split_paths(getattr(args,
                                                             'gtdb_metadata',
                                                             None)),
                               taxonomy=_split_paths(getattr(args,
                                                             'gtdb_taxonomy',
                                                             None)))
        if action == 'map-report':
            return map_report(args.db, getattr(args, 'in', None), args.output,
                              args.to)
        raise ValueError(f'Unknown action: {action}')


def load_configs_safe():
    """Load config.yml, tolerating an absent file."""
    try:
        from neoHGT.util import load_configs
        return load_configs()
    except Exception:
        return None


def cmd_probe(gtdb_dir, out=None, release=None, layout='auto', do_check=True,
              deep=True, metadata=None, taxonomy=None):
    """Probe a dataset and write layout_probe.json (IT-0).

    probe is the pre-flight command, so it must run on the datasets that most
    need probing: one without a metadata table, or one whose tables live under
    unexpected names. Missing metadata is reported, not fatal, and
    --gtdb-metadata/--gtdb-taxonomy point at explicit tables.
    """
    info = detect_layout(gtdb_dir, release=release, layout=layout,
                         metadata=metadata, taxonomy=taxonomy,
                         require_metadata=False)
    if do_check and info.get('metadata'):
        report = selfcheck(info, strict=False, deep=deep)
    elif do_check:
        report = {'status': 'SKIPPED', 'n_checks': 0, 'n_failed': 0,
                  'detail': 'no metadata table found; pass --gtdb-metadata '
                            'to check one explicitly'}
    else:
        report = None
    payload = {'layout': info['layout'], 'release': info['release'],
               'gtdb_dir': info['gtdb_dir'],
               'id_convention': info['id_convention'],
               'serial_anchor_coverage': info['serial_anchor_coverage'],
               'max_first_line_len': info['max_first_line_len'],
               'n_genome_files': info['n_genome_files'],
               'domains': {d: {k: v for k, v in rec.items() if k != 'first'}
                           for d, rec in info['domains'].items()},
               'first_files': {d: rec['first']
                               for d, rec in info['domains'].items()},
               'taxonomy_files': info['taxonomy'],
               'metadata_files': info['metadata'],
               'sample_sets': info['sample_sets'],
               'idlists': info['idlists'],
               'provenance_scripts': info['provenance_scripts'],
               'absent_official_files': info['absent_official'],
               'selfcheck': report,
               'writable_input': assert_input_read_only(gtdb_dir),
               'workspace_root': workspace_root(),
               'generated': str(timestamp())}
    if out:
        out = guard_write(out)
        makedirs(dirname(out), exist_ok=True)
        with open(out, 'w') as f:
            json.dump(payload, f, indent=2, sort_keys=True, default=str)
        print(f'Wrote {out}')
    print(json.dumps({k: payload[k] for k in
                      ('layout', 'release', 'id_convention',
                       'serial_anchor_coverage', 'n_genome_files',
                       'max_first_line_len')}, indent=2))
    if report:
        print(f"selfcheck: {report['status']} "
              f"({report['n_failed']}/{report['n_checks']} failed)")
    return payload


def cmd_taxdump(gtdb_dir, out, release=None, layout='auto',
                id_base=DEFAULT_ID_BASE, manifest=None, do_check=True,
                metadata=None, taxonomy=None):
    """Build a full taxdump from a local GTDB dataset (IT-1).

    A taxdump is a function of the lineage strings, so this step also runs on a
    dataset that ships no metadata table (an unpacked official release, or an
    older release kept only for its taxonomy files).
    """
    check_id_base(id_base)
    prov = GtdbProvider(gtdb_dir, release=release, layout=layout,
                        id_base=id_base, manifest_path=manifest,
                        selfcheck=do_check, require_metadata=False,
                        require_sequences=False,
                        metadata_files=metadata, taxonomy_files=taxonomy)
    prov.prepare()
    started = timestamp()
    taxdump, taxid_map, stats = prov.load_taxonomy()
    out = guard_write(out, must_not_contain=gtdb_dir)
    makedirs(out, exist_ok=True)
    written = write_taxdump(taxdump, join(out, 'taxdump'))
    side = write_side_files(out, prov.info, prov.df_all, taxid_map, taxdump,
                            stats, previous_manifest=manifest)
    finished = timestamp()
    secs = (finished - started).total_seconds()
    print(f'Built {stats["n_nodes"]:,} nodes in {secs:.1f}s -> {out}')
    print('rank counts: ' + ', '.join(
        f'{r}={stats["by_rank"][r]:,}' for r in GTDB_RANKS))
    meta = {'tool': f'neoHGT {tool_version()}', 'action': 'taxdump',
            'gtdb_release': prov.info['release'],
            'local_layout': prov.info['layout'],
            'lineage_only': bool(prov.lineage_only),
            'skipped_selfcheck_items': (
                (prov.report or {}).get('skipped', [])),
            'id_base': id_base, 'counts': public_stats(stats),
            'seconds': secs, 'files': {**written, **side},
            'selfcheck': (prov.report or {}).get('status'),
            'build_finished': str(finished)}
    with open(join(out, 'meta.json'), 'w') as f:
        json.dump(meta, f, indent=2, sort_keys=True, default=str)
    return meta


def main(argv=None):
    """Entry point for 'python -m neoHGT.gtdb'."""
    import argparse
    ap = argparse.ArgumentParser(
        prog='python -m neoHGT.gtdb', description=description)
    sub = ap.add_subparsers(dest='action', required=True)

    p = sub.add_parser('probe', help='detect layout and check consistency')
    p.add_argument('--dir', required=True, help='GTDB dataset root (read-only)')
    p.add_argument('--out', help='destination layout_probe.json')
    p.add_argument('--gtdb-release', help='expected release label')
    p.add_argument('--local-layout', default='auto',
                   choices=list(LOCAL_LAYOUTS) + ['auto'])
    p.add_argument('--gtdb-metadata',
                   help='explicit metadata TSV path(s), comma-delimited')
    p.add_argument('--gtdb-taxonomy',
                   help='explicit taxonomy TSV path(s), comma-delimited')
    p.add_argument('--no-selfcheck', action='store_true')
    p.add_argument('--no-deep', action='store_true',
                   help='skip the per-accession alignment checks')

    t = sub.add_parser('taxdump', help='build a taxdump from GTDB lineages')
    t.add_argument('--dir', required=True)
    t.add_argument('--out', required=True)
    t.add_argument('--gtdb-release')
    t.add_argument('--local-layout', default='auto',
                   choices=list(LOCAL_LAYOUTS) + ['auto'])
    t.add_argument('--gtdb-taxid-base', type=int, default=DEFAULT_ID_BASE)
    t.add_argument('--gtdb-id-manifest')
    t.add_argument('--gtdb-metadata',
                   help='explicit metadata TSV path(s), comma-delimited')
    t.add_argument('--gtdb-taxonomy',
                   help='explicit taxonomy TSV path(s), comma-delimited')
    t.add_argument('--no-selfcheck', action='store_true')

    m = sub.add_parser('map-report', help='map donor taxIds to NCBI names')
    m.add_argument('--db', required=True)
    m.add_argument('--in', dest='in', required=True)
    m.add_argument('--out', required=True)
    m.add_argument('--to', default='ncbi', choices=['ncbi'])

    g = sub.add_parser('gtdbtk', help='convert a GTDB-Tk summary into '
                       '"--input-tax" taxIds of a GTDB database')
    g.add_argument('--in', dest='in', required=True,
                   help='GTDB-Tk classify / de novo summary table')
    g.add_argument('--db', required=True,
                   help='database built with neoHGT --source gtdb')
    g.add_argument('--out', required=True,
                   help='destination genome<TAB>taxId file')
    g.add_argument('--report', help='destination JSON report')
    g.add_argument('--min-rank', choices=list(GTDB_RANKS),
                   help='drop genomes assigned below this rank')

    r = sub.add_parser('reconcile', help='write the accession-level '
                       'difference lists between metadata, taxonomy and the '
                       'local FASTA copy (plan DoD-2)')
    r.add_argument('--dir', required=True, help='GTDB dataset root (read-only)')
    r.add_argument('--out', required=True, help='destination directory')
    r.add_argument('--gtdb-release')
    r.add_argument('--local-layout', default='auto',
                   choices=list(LOCAL_LAYOUTS) + ['auto'])
    r.add_argument('--report', help='destination JSON summary of the counts')

    t2 = sub.add_parser('tree-dist', help='patristic distances in a GTDB '
                        'per-rank tree around one taxon (phylogeny-assisted '
                        'grouping evidence, plan P7)')
    t2.add_argument('--gtdb-dir', help='GTDB release root (uses tree_*/)')
    t2.add_argument('--tree', help='a specific Newick tree file')
    t2.add_argument('--rank', choices=['phylum', 'class', 'order', 'family',
                                      'genus'], default='genus',
                    help='rank tree to use when --gtdb-dir is given')
    t2.add_argument('--focus', required=True,
                    help='label of the focal taxon, e.g. g__Halorubrum')
    t2.add_argument('--max-distance', type=float)
    t2.add_argument('--top-k', type=int)
    t2.add_argument('--out', help='destination TSV of label<TAB>distance')
    t2.add_argument('--report', help='destination JSON summary')
    t2.add_argument('--db', help='GTDB database directory: also resolve each '
                    'label to its synthetic taxId (gtdb/labels.tsv)')
    t2.add_argument('--close-tax-out', help='write the resolved taxIds as a '
                    'comma-separated list, ready for analyze --close-tax')

    a = ap.parse_args(argv)
    if a.action == 'probe':
        res = cmd_probe(a.dir, a.out, a.gtdb_release, a.local_layout,
                        not a.no_selfcheck, not a.no_deep,
                        metadata=_split_paths(a.gtdb_metadata),
                        taxonomy=_split_paths(a.gtdb_taxonomy))
        failed = (res.get('selfcheck') or {}).get('n_failed', 0)
        return 1 if failed else 0
    if a.action == 'taxdump':
        cmd_taxdump(a.dir, a.out, a.gtdb_release, a.local_layout,
                    a.gtdb_taxid_base, a.gtdb_id_manifest,
                    not a.no_selfcheck,
                    metadata=_split_paths(a.gtdb_metadata),
                    taxonomy=_split_paths(a.gtdb_taxonomy))
        return 0
    if a.action == 'map-report':
        stats = map_report(a.db, getattr(a, 'in'), a.out, a.to)
        print(stats)
        return 0
    if a.action == 'gtdbtk':
        rep = gtdbtk_to_input_tax(getattr(a, 'in'), a.db, a.out, a.min_rank)
        print(f'{rep["n_mapped"]}/{rep["n_rows"]} genomes mapped to taxIds of '
              f'{a.db}; wrote {a.out}.')
        if rep['matched_rank']:
            print('  deepest assigned rank: ' + ', '.join(
                f'{k}={v}' for k, v in sorted(rep['matched_rank'].items(),
                                              key=lambda x: RANK_DEPTH[x[0]])))
        if rep['n_unmapped']:
            print(f'  WARNING: {rep["n_unmapped"]} genomes unmapped '
                  f'(first: {rep["unmapped"][0]["genome"]}: '
                  f'{rep["unmapped"][0]["reason"]}).')
        if a.report:
            path = guard_write(a.report, must_not_contain=a.db)
            makedirs(dirname(path), exist_ok=True)
            with open(path, 'w') as f:
                json.dump(rep, f, indent=1, default=str)
            print(f'  report written to {a.report}.')
        return 0 if rep['n_mapped'] else 1
    if a.action == 'reconcile':
        res = reconcile(gtdb_dir=a.dir, out=a.out, release=a.gtdb_release,
                        layout=a.local_layout)
        print(f'metadata {res["n_metadata_genomes"]:,} genomes '
              f'({res["n_metadata_reps"]:,} representatives), '
              f'{res["n_local_faa"]:,} local FASTA, '
              f'{res["n_taxonomy_rows"]:,} taxonomy rows.')
        for name, n in sorted(res['lists'].items()):
            print(f'  {name:34} {n:>8,}')
        print('  -> lists written to ' + a.out +
              (' (all empty: the copy is consistent)'
               if res['consistent'] else ''))
        if a.report:
            path = guard_write(a.report, must_not_contain=a.dir)
            makedirs(dirname(path), exist_ok=True)
            with open(path, 'w') as f:
                json.dump(res, f, indent=1, default=str)
            print(f'  summary written to {a.report}.')
        return 0 if res['consistent'] else 1
    if a.action == 'tree-dist':
        if not a.tree and not a.gtdb_dir:
            raise ValueError('tree-dist needs --tree or --gtdb-dir.')
        trees = [a.tree] if a.tree else rank_trees(a.gtdb_dir).get(a.rank, [])
        if not trees:
            raise ValueError(f'no {a.rank}-rank tree found under {a.gtdb_dir}')
        chosen, rows, summary = None, None, None
        errors = []
        for cand in trees:
            try:
                rows, summary = tree_close_taxa(
                    cand, a.focus, a.max_distance, a.top_k)
                chosen = cand
                break
            except ValueError as e:
                errors.append(f'{basename(cand)}: {e}')
        if chosen is None:
            raise ValueError('; '.join(errors))
        if a.db:
            tid_of = {v: k for k, v in load_labels(a.db).items()}
            if not tid_of:
                raise ValueError(
                    f'{a.db} has no gtdb/labels.tsv; build the database first.')
            rows = [(l, d, tid_of.get(l, '')) for l, d in rows]
            summary['db'] = a.db
            summary['taxids_unmapped'] = sum(1 for r in rows if not r[2])
        print(f'{summary["n_taxa"]:,} taxa in {basename(chosen)}; focus '
              f'{a.focus}: median {summary["median"]:.4f}, '
              + ', '.join(f'<= {k}: {v}' for k, v in
                          sorted(summary['within'].items(),
                                 key=lambda x: float(x[0])))
              + f'; returned {summary["returned"]}.')
        if a.db:
            print(f'  {sum(1 for r in rows if not r[2])} of {len(rows)} '
                  f'returned labels did not resolve to a taxId in {a.db}.')
        if a.out:
            path = guard_write(a.out)
            makedirs(dirname(path), exist_ok=True)
            with open(path, 'w') as f:
                f.write('label\tdistance' + ('\ttaxid' if a.db else '') + '\n')
                for r in rows:
                    f.write('\t'.join(str(x) if not isinstance(x, float)
                                      else f'{x:.6f}' for x in r) + '\n')
            summary['table'] = path
        if a.close_tax_out:
            ids = [r[2] for r in rows if len(r) > 2 and r[2]]
            path = guard_write(a.close_tax_out)
            makedirs(dirname(path), exist_ok=True)
            with open(path, 'w') as f:
                f.write(','.join(ids) + '\n')
            summary['close_tax'] = path
            summary['close_tax_n'] = len(ids)
            print(f'  phylogenetic close group: {len(ids)} taxIds written to '
                  f'{a.close_tax_out}.')
        if a.report:
            path = guard_write(a.report)
            makedirs(dirname(path), exist_ok=True)
            summary['nearest'] = [list(r[:3]) for r in rows[:20]]
            with open(path, 'w') as f:
                json.dump(summary, f, indent=1)
            print(f'  report written to {a.report}.')
        return 0
    return 2


def cli_main(argv=None):
    """Console entry: report input-contract failures instead of dumping a traceback.

    `neoHGT gtdb probe|taxdump|reconcile` is the command a user runs *before*
    trusting a dataset, so the most likely outcome is a refused dataset. The
    reason is a one-line sentence (`main()` raises ValueError / SelfCheckError
    with it), and IT-0 asks for "non-zero exit plus the differences listed" --
    which a stack trace satisfies only in the exit code.
    """
    try:
        return main(argv)
    except (ValueError, SelfCheckError, GuardError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(cli_main())
