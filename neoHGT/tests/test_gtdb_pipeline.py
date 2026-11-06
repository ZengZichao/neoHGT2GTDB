"""Offline integration tests for the GTDB pipeline (UT-0.1/0.3, UT-3.x,
UT-3.6 no-network, UT-5.x workspace, DoD-1/3/4).

Everything here runs from the committed KB-scale fixtures under
neoHGT/tests/data/GTDB; the read-only GTDB dataset is never needed, and no
network access is attempted (asserted explicitly below).
"""

import gzip
import hashlib
import inspect
import io
import json
import ast
import contextlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
import urllib.request
import pandas as pd
from os.path import (basename, exists, getsize, isdir, isfile,
                     islink, join)

from neoHGT import gtdb
from neoHGT import util
from neoHGT.search import Search, guard_precomp_knobs, \
    parse_def_table_checked, taxid_wf_gtdb
from neoHGT.analyze import write_hgt_list_gtdb
from neoHGT.util import find_lca

HERE = os.path.dirname(os.path.abspath(__file__))
FX = join(HERE, 'data', 'GTDB')
FIXTURES = json.load(open(join(FX, 'fixtures.json'))) if isfile(
    join(FX, 'fixtures.json')) else None


def tmpdir():
    # the workspace guard (and NEOHGT_ROOT) is the boundary these tests run
    # under, so the scratch dir is created inside it -- under NEOHGT_TMP (the
    # designated scratch area) when there is one, so leaks never land in the
    # workspace root; TMPDIR is the fallback for a session with no workspace
    root = os.environ.get('NEOHGT_TMP') or os.environ.get('NEOHGT_ROOT')
    if root and isdir(root):
        os.makedirs(root, exist_ok=True)
        d = tempfile.mkdtemp(prefix='neohgt_test_', dir=root)
    else:
        d = tempfile.mkdtemp(prefix='neohgt_test_',
                             dir=os.environ.get('TMPDIR') or None)
    os.chmod(d, 0o700)
    return d


def renamed_tree(d):
    """Lay the committed fixtures out as the verified local ('renamed') layout.

    metadata_bac120/ + tree_bac120/ + 00_*_aa/protein_faa_reps/<domain>/, which
    is what both the layout checks and the offline build test need.
    """
    for name, dst in (('metadata_bac120_sample.tsv',
                       'metadata_bac120/bac120_metadata_r220.tsv'),
                      ('metadata_ar53_sample.tsv',
                       'metadata_ar53/ar53_metadata_r220.tsv')):
        shutil.copy(join(FX, name), join(d, *dst.split('/')))
    rows = []
    for name in ('metadata_bac120_sample.tsv', 'metadata_ar53_sample.tsv'):
        with open(join(FX, name)) as f:
            head = f.readline().rstrip('\n').split('\t')
            i, j = head.index('accession'), head.index('gtdb_taxonomy')
            for line in f:
                x = line.rstrip('\n').split('\t')
                rows.append((x[i], x[j]))
    for sub, fam in (('tree_bac120', 'bac120'), ('tree_ar53', 'ar53')):
        dom = 'bacteria' if fam == 'bac120' else 'archaea'
        with open(join(d, sub, f'{fam}_taxonomy_r220.tsv'), 'w') as g:
            for acc, lin in rows:
                if lin.startswith('d__' + dom.capitalize()):
                    g.write(f'{acc}\t{lin}\n')
    for rec in FIXTURES['genomes']:
        dst = join(d, '00_GTDB_r220_aa/protein_faa_reps', rec['domain'],
                   os.path.basename(rec['kb_excerpt']))
        shutil.copy(rec['kb_excerpt'], dst)
    return rows


def cli_args(command, argv):
    """Parse command-line arguments through the real scripts/neoHGT parser.

    Using the shipped parser (instead of a Namespace built by hand) keeps the
    CLI surface covered by the integration tests: an option that is declared
    but not wired, or wired but not declared, fails here.
    """
    import importlib.util
    from importlib.machinery import SourceFileLoader
    from importlib import import_module
    script = join(HERE, os.pardir, os.pardir, 'scripts', 'neoHGT')
    # the CLI is an extension-less script, so the loader has to be named
    loader = SourceFileLoader('neoHGT_cli', script)
    spec = importlib.util.spec_from_file_location('neoHGT_cli', script,
                                                  loader=loader)
    cli = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cli)
    module = getattr(import_module(f'neoHGT.{command}'),
                     command.capitalize())()
    old = sys.argv
    sys.argv = ['neoHGT', command] + list(argv)
    try:
        return cli.parse_args([(command, module)])
    finally:
        sys.argv = old


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


class TestFixtures(unittest.TestCase):
    """The fixtures must be verifiable copies of a read-only source."""

    @unittest.skipIf(FIXTURES is None, 'fixtures not built')
    def test_hashes_recorded_and_copies_match_manifest(self):
        genomes = FIXTURES['genomes']
        self.assertLessEqual(len(genomes), 10,
                             'safety cap: at most 10 genomes may be copied')
        self.assertGreaterEqual(len(genomes), 4)
        # recorded paths are workspace-agnostic: the excerpt is repo-relative,
        # the source is release-relative, the cache copy is cache-relative --
        # resolve each against its own anchor instead of the current directory
        repo = join(HERE, os.pardir, os.pardir)
        cache = join(repo, os.pardir, 'neoHGT2GTDB-缓存测试')
        release = os.environ.get('GTDB_DIR') or FIXTURES['gtdb_dir']
        for rec in genomes:
            kb = rec['kb_excerpt']
            if not os.path.isabs(kb):
                kb = join(repo, kb)
            self.assertTrue(isfile(kb), kb)
            # the digests are compared against the *recorded* values (the
            # previous form hashed the excerpt and compared it with itself,
            # which no substitution could ever fail)
            copy = rec['fixture_copy']
            if not os.path.isabs(copy):
                copy = join(cache, copy)
            if isfile(copy):
                self.assertEqual(sha256(copy), rec['fixture_sha256'],
                                 'fixture copy must match its recorded digest')
                self.assertEqual(sha256(copy), rec['source_sha256'],
                                 'fixture copy must equal the read-only source')
            source = rec['source_path']
            if not os.path.isabs(source):
                source = join(release, source)
            # the excerpt must be the head of the same file, byte for byte
            if isfile(source):
                self.assertEqual(sha256(source),
                                 rec['source_sha256'],
                                 'source must match its recorded digest')
                self.assertEqual(getsize(source),
                                 rec['source_bytes'])
                with open(source, 'rb') as f:
                    head = f.read(getsize(kb))
                with open(kb, 'rb') as f:
                    self.assertEqual(f.read(), head,
                                     'excerpt must be a byte-prefix of the '
                                     'source (first kb_records records)')
            self.assertTrue(not os.path.isabs(rec['source_path'])
                            and not rec['source_path'].startswith('..'),
                            'source must be recorded release-relative and '
                            'never escape the declared read-only dataset')

    def test_layout_label_recorded(self):
        self.assertTrue(gtdb.is_renamed_layout(
            open(join(FX, 'layout.txt')).read().strip()))

    def test_metadata_slices_are_column_named(self):
        df = gtdb.load_metadata([join(FX, 'metadata_bac120_sample.tsv'),
                                 join(FX, 'metadata_ar53_sample.tsv')],
                                release='220')
        for col in ('accession', 'genome', 'ncbi_taxid', 'gtdb_taxonomy',
                    'proteins', 'domain', 'quality_score'):
            self.assertIn(col, df.columns)
        # the RS_/GB_ prefix must be removed with an anchored regex, not
        # str.lstrip (which would also eat GCA_/GCF_)
        for acc, gen in zip(df['accession'], df['genome']):
            self.assertRegex(gen, r'^(BAC|ARC)_(GCF|GCA)_\d+$')


class FakeInput(unittest.TestCase):
    """Synthesize a minimal GTDB dataset directory for layout checks."""

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        for sub in ('tree_bac120', 'tree_ar53', 'metadata_bac120',
                    'metadata_ar53',
                    '00_GTDB_r220_aa/protein_faa_reps/bacteria',
                    '00_GTDB_r220_aa/protein_faa_reps/archaea'):
            os.makedirs(join(self.d, sub), exist_ok=True)


class TestLayoutDetection(FakeInput):

    def test_renamed_layout_recognized(self):
        shutil.copy(join(FX, 'metadata_bac120_sample.tsv'),
                    join(self.d, 'metadata_bac120',
                         'bac120_metadata_r220.tsv'))
        shutil.copy(join(FX, 'metadata_ar53_sample.tsv'),
                    join(self.d, 'metadata_ar53', 'ar53_metadata_r220.tsv'))
        # taxonomy slices without header
        rows = []
        for name in ('metadata_bac120_sample.tsv', 'metadata_ar53_sample.tsv'):
            with open(join(FX, name)) as f:
                head = f.readline().rstrip('\n').split('\t')
                i = head.index('accession')
                j = head.index('gtdb_taxonomy')
                for line in f:
                    x = line.rstrip('\n').split('\t')
                    rows.append((x[i], x[j]))
        for path, dom in ((join(self.d, 'tree_bac120',
                                'bac120_taxonomy_r220.tsv'), 'bacteria'),
                          (join(self.d, 'tree_ar53',
                                'ar53_taxonomy_r220.tsv'), 'archaea')):
            with open(path, 'w') as g:
                for acc, lin in rows:
                    if lin.startswith('d__' + dom.capitalize()):
                        g.write(f'{acc}\t{lin}\n')
        for rec in FIXTURES['genomes']:
            dst = join(self.d, '00_GTDB_r220_aa/protein_faa_reps',
                       rec['domain'], os.path.basename(rec['kb_excerpt']))
            shutil.copy(rec['kb_excerpt'], dst)
        info = gtdb.detect_layout(self.d)
        self.assertTrue(gtdb.is_renamed_layout(info['layout']),
                        info['layout'])
        self.assertEqual(info['release'], '220')
        self.assertEqual(info['serial_anchor_coverage'], 1.0)

    def test_missing_taxonomy_is_reported(self):
        with self.assertRaises(ValueError) as ctx:
            gtdb.detect_layout(self.d)
        self.assertIn('tree_bac120', str(ctx.exception))

    def _flat_release_tables(self, gz=False, truncate=False):
        """Write the fixture lineages the way an unpacked release ships them.

        Official GTDB releases are flat -- 'bac120_taxonomy_r220.tsv.gz' sits in
        the dataset root, next to 'metadata' and the sequence directories, with
        no 'tree_bac120/' wrapper. The local r214 copy on disk is laid out that
        way, which is what exposed that discovery only looked inside wrappers.
        """
        rows = []
        for line in open(join(FX, 'taxonomy_sample.tsv')):
            acc, lin = line.rstrip('\n').split('\t')
            if truncate:
                # a third-party rejoin: ranks glued with '_' and the species
                # rank dropped (this is the shape of the local r226 copy)
                lin = '_'.join(lin.split(';')[:-1])
            rows.append((acc, lin))
        written = []
        for fam, dom in (('bac120', 'd__Bacteria'), ('ar53', 'd__Archaea')):
            name = f'{fam}_taxonomy_r999' + ('.tsv.gz' if gz else '.tsv')
            path = join(self.d, name)
            body = ''.join(f'{a}\t{l}\n' for a, l in rows
                           if l.startswith(dom))
            if gz:
                with gzip.open(path, 'wt', encoding='utf-8') as f:
                    f.write(body)
            else:
                with open(path, 'w') as f:
                    f.write(body)
            written.append(path)
        return written

    def test_official_flat_layout_is_discovered(self):
        for gz in (False, True):
            d = tmpdir()
            self.addCleanup(shutil.rmtree, d, ignore_errors=True)
            saved, self.d = self.d, d
            try:
                self._flat_release_tables(gz=gz)
                info = gtdb.detect_layout(d, require_metadata=False,
                                          require_sequences=False)
            finally:
                self.d = saved
            self.assertEqual(info['layout'],
                             'official-gz' if gz else 'official-plain',
                             f'gz={gz}')
            self.assertEqual(info['release'], '999')
            self.assertTrue(info['lineage_only'])

    def test_wrapper_directories_take_precedence_over_root_files(self):
        # the verified local copy lives in tree_bac120/; a stray root-level
        # table of another release must never silently replace it
        self._flat_release_tables()
        shutil.copy(join(FX, 'taxonomy_sample.tsv'),
                    join(self.d, 'tree_bac120', 'bac120_taxonomy_r220.tsv'))
        info = gtdb.find_input_files(self.d)
        self.assertTrue(info['taxonomy'])
        for p in info['taxonomy']:
            self.assertNotIn('_r999', p)

    def test_lineage_only_dataset_builds_a_taxdump(self):
        paths = self._flat_release_tables()
        prov = gtdb.GtdbProvider(self.d, require_metadata=False,
                                 require_sequences=False)
        prov.prepare()
        self.assertTrue(prov.lineage_only)
        self.assertEqual(prov.report['status'], 'OK')
        self.assertEqual(prov.report['skipped'],
                         gtdb.metadata_dependent_checks(prov.info))
        self.assertEqual(prov.report['n_skipped'], 6)
        self.assertEqual(len(prov.df_all), 0)
        lins = sorted(set(gtdb.read_taxonomy(paths).values()))
        self.assertEqual(prov.lineages, lins)
        _, tmap, stats = prov.load_taxonomy()
        nodes = set()
        for x in lins:
            nodes.update(gtdb.lineage_paths(x))
        self.assertEqual(stats['n_gtdb_nodes'], len(nodes))
        self.assertEqual(stats['by_rank']['species'], len(lins))
        self.assertEqual(len(tmap), len(nodes))

    def test_database_refuses_a_dataset_without_metadata(self):
        self._flat_release_tables()
        with self.assertRaises(ValueError) as ctx:
            gtdb.detect_layout(self.d)
        self.assertIn('--gtdb-metadata', str(ctx.exception))
        self.assertIn('taxdump', str(ctx.exception))

    def test_table_overrides_are_applied_before_discovery_fails(self):
        # files that no discovery rule can name must still be usable through
        # --gtdb-taxonomy / --gtdb-metadata
        elsewhere = join(self.d, 'elsewhere')
        os.makedirs(elsewhere, exist_ok=True)
        tax = join(elsewhere, 'lineages.tsv')
        meta = join(elsewhere, 'table.tsv')
        shutil.copy(join(FX, 'taxonomy_sample.tsv'), tax)
        shutil.copy(join(FX, 'metadata_bac120_sample.tsv'), meta)
        prov = gtdb.GtdbProvider(self.d, taxonomy_files=[tax],
                                 metadata_files=[meta], selfcheck=False,
                                 require_sequences=False)
        prov.prepare()
        self.assertFalse(prov.lineage_only)
        self.assertEqual([p for p in prov.info['taxonomy']], [tax])
        self.assertGreater(len(prov.df_all), 0)

    def test_rejoined_lineages_are_refused_before_any_write(self):
        # the local r226 copy glues ranks with '_' and stops at genus; turning
        # that into a taxdump would silently invent a shallower tree
        self._flat_release_tables(truncate=True)
        prov = gtdb.GtdbProvider(self.d, require_metadata=False,
                                 require_sequences=False)
        with self.assertRaises(gtdb.SelfCheckError) as ctx:
            prov.prepare()
        self.assertIn('does not have 7 ranks', str(ctx.exception))

    def _full_renamed_tree(self):
        """Materialise the verified local layout: metadata + taxonomy + faa."""
        renamed_tree(self.d)

    def test_reported_skips_are_real_check_names(self):
        """The lineage-only skip list must name checks that actually exist.

        A hand-copied list drifts: the first version invented three names that
        `selfcheck()` never produced and printed the wrong count. Build a
        dataset that *does* have a metadata table, run the full §3.6 battery,
        and require every reported skip to be one of its item names.
        """
        self._full_renamed_tree()
        info = gtdb.find_input_files(self.d)
        info['release'] = '220'
        rep = gtdb.selfcheck(info, strict=False, deep=True)
        ran = {i['name'] for i in rep['checks']}
        skipped = gtdb.metadata_dependent_checks(info)
        self.assertGreaterEqual(len(ran), len(skipped),
                                'the battery must actually run these checks')
        for name in skipped:
            self.assertIn(name, ran, f'{name} is not a real check name')
        self.assertEqual(
            set(gtdb.selfcheck_lineages(info)['skipped']), set(skipped))

    def test_refused_dataset_reports_one_line_instead_of_a_traceback(self):
        """IT-0 asks for "non-zero exit with the differences listed".

        A stack trace has a non-zero exit code but is not a report, and this is
        the command people run *before* trusting a dataset. main() keeps raising
        for library callers; only the console entry turns it into a message.
        """
        bad = join(self.d, 'elsewhere')
        os.makedirs(bad, exist_ok=True)
        with open(join(bad, 'bac120_taxonomy_r999.tsv'), 'w') as f:
            f.write('GB_GCA_000000000.1\td__Bacteria_p__X_c__Y_o__Z_f__F'
                    '_g__G\n')
        out = join(self.d, 'should_not_exist')
        with contextlib.redirect_stderr(io.StringIO()) as err:
            rc = gtdb.cli_main(['taxdump', '--dir', bad, '--out', out])
        self.assertEqual(rc, 2, rc)
        self.assertIn('ERROR:', err.getvalue())
        self.assertIn('does not have 7 ranks', err.getvalue())
        self.assertNotIn('Traceback', err.getvalue())
        self.assertFalse(isdir(out), 'nothing may be written for a refused '
                                     'dataset')

    def test_undecidable_layout_lists_candidates(self):
        # taxonomy + metadata present, sequences absent
        shutil.copy(join(FX, 'metadata_bac120_sample.tsv'),
                    join(self.d, 'metadata_bac120',
                         'bac120_metadata_r220.tsv'))
        shutil.copy(join(FX, 'metadata_ar53_sample.tsv'),
                    join(self.d, 'metadata_ar53', 'ar53_metadata_r220.tsv'))
        for sub, name in (('tree_bac120', 'bac120_taxonomy_r220.tsv'),
                          ('tree_ar53', 'ar53_taxonomy_r220.tsv')):
            shutil.copy(join(FX, 'taxonomy_sample.tsv'), join(self.d, sub,
                                                                name))
        # remove the empty per-genome directories the fixture scaffold made,
        # so that no sequence source exists at all (test-owned scratch dirs)
        for d in ('00_GTDB_r220_aa/protein_faa_reps/archaea',
                  '00_GTDB_r220_aa/protein_faa_reps/bacteria',
                  '00_GTDB_r220_aa/protein_faa_reps', '00_GTDB_r220_aa'):
            os.rmdir(join(self.d, d))
        with self.assertRaises(ValueError) as ctx:
            gtdb.detect_layout(self.d)
        msg = str(ctx.exception)
        self.assertIn('protein_faa_reps', msg)
        for layout in gtdb.LOCAL_LAYOUTS:
            self.assertIn(layout, msg)

    def test_forced_layout_validated(self):
        with self.assertRaises(ValueError):
            gtdb.detect_layout(self.d, layout='nonsense')


class TestStreamExtraction(unittest.TestCase):
    """DoD-3 / UT-3.1: 100% of emitted sequence Ids resolve to a taxId."""

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        df = gtdb.load_metadata([join(FX, 'metadata_bac120_sample.tsv'),
                                 join(FX, 'metadata_ar53_sample.tsv')],
                                release='220')
        self.lineages = list(df['gtdb_taxonomy'].unique())
        self.td, self.tmap, self.stats = gtdb.build_taxdump_from_gtdb(
            self.lineages)

    def _fetches(self):
        out = []
        for rec in FIXTURES['genomes']:
            out.append(gtdb.GenomeFetch(rec['accession'], None, None,
                                        path=rec['kb_excerpt']))
        return out

    def test_every_protein_gets_a_taxid(self):
        fetches = self._fetches()
        # attach the species taxId of each genome's lineage
        acc2lin = {}
        for name in ('metadata_bac120_sample.tsv',
                    'metadata_ar53_sample.tsv'):
            with open(join(FX, name)) as f:
                head = f.readline().rstrip('\n').split('\t')
                ia, il = head.index('accession'), head.index('gtdb_taxonomy')
                for line in f:
                    x = line.rstrip('\n').split('\t')
                    acc2lin[x[ia]] = x[il]
        fetches = [gtdb.GenomeFetch(f.accession, None,
                                    self.tmap.get(acc2lin.get(f.accession)),
                                    path=f.path) for f in fetches]
        faa = join(self.d, 'db.faa')
        tmap = join(self.d, 'taxon.map')
        stats = gtdb.extract_genomes_stream(
            fetches, faa, tmap, 'gtdb-r220-renamed',
            blast_taxidmap=join(self.d, 'seq2tid.txt'),
            genome_map=join(self.d, 'gmap.txt'))
        self.assertGreater(stats['proteins'], 0)
        self.assertEqual(stats['proteins_no_taxid'], 0)
        self.assertEqual(stats['genomes_missing'], 0)
        n_seq = sum(1 for line in open(faa) if line.startswith('>'))
        n_map = sum(1 for _ in open(tmap)) - 1
        self.assertEqual(n_seq, n_map,
                         'UT-3.1: taxonmap rows must equal db.faa records')
        self.assertEqual(stats['proteins'], n_seq)
        # artifact-level invariant (review A 4.1): GTDB's prodigal output ends
        # every CDS in '*', and the old rstrip order left that star on ~99.6%
        # of the extracted proteins. The published db.faa must be clean.
        with open(faa) as f:
            dirty = [line.split()[0][1:] for line in f
                     if not line.startswith('>') and line.rstrip().endswith('*')]
        self.assertEqual(dirty, [],
                         f'db.faa must not contain stop-codon tails: '
                         f'{dirty[:5]}')
        with open(faa) as f:
            partial = [line for line in f if line.startswith('>')]
        self.assertFalse(any(isfile(faa + '.partial') for _ in [0]),
                         'the .partial scratch file must be renamed away')

    def test_taxonmap_key_columns_are_identical(self):
        fetches = [gtdb.GenomeFetch(r['accession'], None, '1006000001',
                                    path=r['kb_excerpt'])
                   for r in FIXTURES['genomes']]
        tmap = join(self.d, 'taxon.map')
        gtdb.extract_genomes_stream(fetches, join(self.d, 'db.faa'), tmap,
                                    'gtdb-r220-renamed')
        with open(tmap) as f:
            self.assertEqual(f.readline().rstrip('\n'),
                             'accession\taccession.version\ttaxid\tgi')
            for line in f:
                x = line.rstrip('\n').split('\t')
                self.assertEqual(x[0], x[1],
                                 'both key columns must carry the full Id')
                self.assertIn(gtdb.SERIAL_ANCHOR, x[0])

    def test_resume_counts_and_unresolved(self):
        missing = gtdb.GenomeFetch('GB_GCA_999999999.1', None, '1',
                                   path=join(self.d, 'absent.faa'))
        stats = gtdb.extract_genomes_stream([missing], join(self.d, 'a.faa'),
                                            join(self.d, 'a.map'),
                                            'gtdb-r220-renamed')
        self.assertEqual(stats['genomes_missing'], 1)
        self.assertEqual(stats['unresolved'], ['GB_GCA_999999999.1'])



class _NetworkGuard(object):
    """Refuse every network exit the pipeline could take.

    Patching ``urllib.request.urlopen`` and ``neoHGT.util.run_command`` guarded
    nothing: search.py/database.py do ``from ... import urlopen/run_command``,
    which binds the name into their own namespace at import time. The patch has
    to land on the using modules, and a socket-level backstop catches anything
    the two named exits miss.
    """

    def __init__(self, calls):
        self.calls = calls
        self._saved = {}

    def _boom(self, *a, **k):
        self.calls.append(a)
        raise AssertionError(f'network access attempted: {a}')

    def __enter__(self):
        import neoHGT.search as _search_mod
        import neoHGT.database as _db_mod
        import socket
        targets = {
            (_search_mod, 'urlopen'): getattr(_search_mod, 'urlopen', None),
            (_search_mod, 'run_command'): getattr(_search_mod, 'run_command',
                                                  None),
            (_db_mod, 'run_command'): getattr(_db_mod, 'run_command', None),
            (urllib.request, 'urlopen'): urllib.request.urlopen,
            (util, 'run_command'): getattr(util, 'run_command', None),
        }
        for (mod, name), original in targets.items():
            if original is not None:
                self._saved[(mod, name)] = original
                setattr(mod, name, self._boom)
        self._saved[(socket, 'socket')] = socket.socket
        socket.socket = self._boom
        return self

    def __exit__(self, *exc):
        import socket
        for (mod, name), original in self._saved.items():
            setattr(mod, name, original)
        return False


class TestNoNetwork(unittest.TestCase):
    """DoD-4 / UT-3.6: a local GTDB build must not touch the network."""

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.calls = []
        self._guard = _NetworkGuard(self.calls)
        self._guard.__enter__()

    def tearDown(self):
        self._guard.__exit__()

    def test_local_ingestion_and_conversion_are_offline(self):
        df = gtdb.load_metadata([join(FX, 'metadata_bac120_sample.tsv'),
                                 join(FX, 'metadata_ar53_sample.tsv')],
                                release='220')
        td, tmap, stats = gtdb.build_taxdump_from_gtdb(
            list(df['gtdb_taxonomy'].unique()))
        gtdb.write_taxdump(td, join(self.d, 'taxdump'))
        fetches = [gtdb.GenomeFetch(r['accession'], None,
                                    tmap.get(df.iloc[i]['gtdb_taxonomy'], '1'),
                                    path=r['kb_excerpt'])
                   for i, r in enumerate(FIXTURES['genomes'])]
        st = gtdb.extract_genomes_stream(
            fetches, join(self.d, 'db.faa'), join(self.d, 'taxon.map'),
            'gtdb-r220-renamed')
        self.assertGreater(st['proteins'], 0)
        self.assertEqual(self.calls, [], 'local path must be fully offline')


class TestHitTableValidation(unittest.TestCase):
    """UT-3.2 / RK5: a missing staxids column must fail loudly."""

    def test_six_columns_rejected(self):
        s = Search.__new__(Search)
        s.source = 'auto'
        s.gtdb_root = None
        s.evalue = s.identity = s.coverage = s.maxhits = 0
        lines = ['q1\ts1\t95.0\t1e-50\t200.0\t100']
        with self.assertRaises(ValueError) as ctx:
            parse_def_table_checked(s, lines)
        self.assertIn('staxids', str(ctx.exception))

    def test_seven_columns_ok(self):
        s = Search.__new__(Search)
        s.source = 'auto'
        s.gtdb_root = None
        s.evalue = s.identity = s.coverage = s.maxhits = 0
        got = parse_def_table_checked(
            s, ['q1\ts1\t95.0\t1e-50\t200.0\t100\t1006000001'])
        self.assertEqual(got['q1'][0]['taxid'], '1006000001')


class TestWorkspaceGuards(unittest.TestCase):
    """UT-5.x: outputs stay inside the workspace and off the input dataset."""

    def setUp(self):
        self.root = os.environ.get(
            'NEOHGT_ROOT',
            os.path.dirname(os.path.dirname(os.path.dirname(
                os.path.abspath(__file__)))))

    def test_symlinked_output_cannot_escape(self):
        # Only the symlink itself is created, and it is created inside the
        # workspace; nothing is ever written to the outside target.
        base = join(self.root, 'neoHGT2GTDB-临时文件')
        os.makedirs(base, exist_ok=True)   # CI checkouts ship no scratch dir
        link = join(base, 'escape_link')
        if exists(link) or islink(link):
            os.remove(link)
        os.symlink('/etc', link)
        try:
            with self.assertRaises(gtdb.GuardError):
                gtdb.guard_write(join(link, 'x.txt'), root=self.root)
        finally:
            os.remove(link)

    def test_input_dir_is_never_a_write_target(self):
        fake_input = join(self.root, 'neoHGT2GTDB-缓存测试', 'fixtures')
        with self.assertRaises(gtdb.GuardError):
            gtdb.guard_write(join(fake_input, 'x'), root=self.root,
                             must_not_contain=fake_input)


class TestDonorReporting(unittest.TestCase):
    """DoD-6: donors must be expressible as official GTDB labels, and the
    offline NCBI mapping side table must be produced.

    The KDE-based prediction itself needs a large reference set, so this test
    drives the reporting layer directly with a prepared score table.
    """

    def setUp(self):
        import pandas as pd
        from neoHGT.analyze import Analyze
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        df = gtdb.load_metadata([join(FX, 'metadata_bac120_sample.tsv'),
                                 join(FX, 'metadata_ar53_sample.tsv')],
                                release='220')
        self.td, self.tmap, self.stats = gtdb.build_taxdump_from_gtdb(
            list(df['gtdb_taxonomy'].unique()))
        info = {'gtdb_dir': os.environ.get('GTDB_DIR', '/nonexistent'),
                'release': '220',
                'layout': 'gtdb-r220-renamed',
                'provenance_scripts': [], 'absent_official': [],
                'serial_anchor_coverage': 1.0, 'max_first_line_len': 270,
                'domains': {}, 'taxonomy': [], 'metadata': [],
                'idlists': {}, 'sample_sets': {}}
        gtdb.write_side_files(self.d, info, df, self.tmap, self.td,
                              self.stats)
        # a species taxId inside Enterobacteriaceae, used as the donor
        self.donor = next(t for k, t in self.tmap.items()
                          if k.endswith(';f__Enterobacteriaceae;g__Escherichia'))
        self.sp = next(t for k, t in self.tmap.items()
                       if 'g__Escherichia;s__' in k)
        an = Analyze.__new__(Analyze)
        an.taxdump = self.td
        an.donor_name = False
        an.donor_rank = 'genus'
        an.donor_label = 'prefixed'
        an.donor_ncbi_map = 'yes'
        an.gtdb_root = self.d
        an.gtdb_labels = gtdb.load_labels(self.d)
        an.output = self.d
        os.makedirs(join(self.d, 'hgts'), exist_ok=True)
        an.df = pd.DataFrame([
            {'sample': 's1', 'protein': 'p1', 'silh': 0.9, 'match': self.sp,
             'hgt': True},
            {'sample': 's1', 'protein': 'p2', 'silh': 0.8, 'match': self.sp,
             'hgt': True}])
        self.an = an

    def test_prefixed_labels_and_ncbi_map(self):
        self.an.__class__.write_hgt_list(self.an, 's1')
        path = join(self.d, 'hgts', 's1.txt')
        self.assertTrue(isfile(path))
        rows = [l.rstrip('\n').split('\t') for l in open(path) if l.strip()]
        self.assertEqual(len(rows), 2)
        for r in rows:
            self.assertEqual(r[2], 'g__Escherichia',
                             'donor must be the official prefixed label')
        m = join(self.d, 'hgts', 's1.donor_ncbi.tsv')
        self.assertTrue(isfile(m))
        lines = open(m).read().splitlines()
        self.assertEqual(
            lines[0], 'protein\tsilh\tdonor_taxid\tdonor_label\tdonor_rank\t'
                      'ncbi_taxids\tnconflict')
        self.assertIn('g__Escherichia', lines[1])

    def test_taxid_mode_leaves_ids(self):
        self.an.donor_label = 'taxid'
        write_hgt_list_gtdb(self.an, 's1')
        rows = [l.rstrip('\n').split('\t')
                for l in open(join(self.d, 'hgts', 's1.txt')) if l.strip()]
        # donor_rank='genus' raises the species donor to its genus taxId
        self.assertEqual(rows[0][2], str(self.donor))


class TestRankAliasesInDump(unittest.TestCase):
    """UT-4.x: 'domain' must work even though the dump stores superkingdom."""

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        df = gtdb.load_metadata([join(FX, 'metadata_bac120_sample.tsv'),
                                 join(FX, 'metadata_ar53_sample.tsv')])
        td, tmap, _ = gtdb.build_taxdump_from_gtdb(
            list(df['gtdb_taxonomy'].unique()))
        gtdb.write_taxdump(td, join(self.d, 'taxdump'))
        self.td = util.read_taxdump(join(self.d, 'taxdump'))
        self.tmap = tmap

    def test_domain_and_superkingdom_both_resolve(self):
        dom = self.tmap['d__Bacteria']
        self.assertEqual(
            util.taxid_at_rank(dom, util.rank_for_dump('domain', self.td),
                               self.td), dom)
        self.assertEqual(
            util.taxid_at_rank(dom, util.rank_for_dump('superkingdom',
                                                        self.td), self.td),
            dom)

    def test_genus_rank_lookup_on_real_gtdb_taxa(self):
        tid = next(v for k, v in self.tmap.items()
                   if k.endswith(';g__Escherichia'))
        # 'domain' must be resolved through the alias actually stored in the
        # dump ('superkingdom'), which is what rank_for_dump does
        self.assertEqual(
            util.taxid_at_rank(
                tid, util.rank_for_dump('domain', self.td), self.td),
            self.tmap['d__Bacteria'])
        self.assertEqual(
            util.taxid_at_rank(tid, util.rank_for_dump('phylum', self.td),
                               self.td),
            self.tmap['d__Bacteria;p__Pseudomonadota'])
        # genus and species are readable under their GTDB rank names
        self.assertEqual(util.taxid_at_rank(tid, 'genus', self.td), tid)


class StubGtdbProvider:
    """Minimal stand-in exposing what HybridProvider needs from GTDB."""

    def __init__(self, lineages, id_base=gtdb.DEFAULT_ID_BASE):
        self.id_base = id_base
        self.lineages = lineages
        self.taxdump = self.taxid_map = self.stats = None
        self.info = {'gtdb_dir': '/nonexistent', 'release': '220',
                     'layout': 'gtdb-r220-renamed'}

    def load_taxonomy(self, trim=None):
        self.taxdump, self.taxid_map, self.stats = \
            gtdb.build_taxdump_from_gtdb(self.lineages)
        return self.taxdump, self.taxid_map, self.stats


class TestHybrid(unittest.TestCase):
    """P5 / IT-5: GTDB prokaryotes grafted under the NCBI cellular root.

    The NCBI subtree is synthesized here because downloading NCBI taxdump is
    not permitted in this environment (project rule 13); the merge logic under
    test is the same that runs against a real NCBI taxdump.
    """

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        df = gtdb.load_metadata([join(FX, 'metadata_bac120_sample.tsv'),
                                 join(FX, 'metadata_ar53_sample.tsv')])
        self.gtdb_prov = StubGtdbProvider(
            list(df['gtdb_taxonomy'].unique()))
        ncbi = {
            '1': {'parent': '1', 'rank': 'no rank', 'name': 'root'},
            '131567': {'parent': '1', 'rank': 'no rank',
                       'name': 'cellular organisms'},
            '2759': {'parent': '131567', 'rank': 'superkingdom',
                     'name': 'Eukaryota'},
            '4932': {'parent': '2759', 'rank': 'genus',
                     'name': 'Saccharomyces'},
            '4930': {'parent': '4932', 'rank': 'species',
                     'name': 'Saccharomyces cerevisiae'},
            '10239': {'parent': '1', 'rank': 'no rank', 'name': 'Viruses'}}
        gtdb.write_taxdump(ncbi, join(self.d, 'ncbi_taxdump'))

    def _hybrid(self, policy='root'):
        return gtdb.HybridProvider(self.gtdb_prov,
                                   join(self.d, 'ncbi_taxdump'),
                                   lca_policy=policy)

    def test_single_root_and_domain_reattachment(self):
        h = self._hybrid()
        td, tmap, stats = h.load_taxonomy()
        self.assertEqual(find_lca(list(td), td), gtdb.ROOT_TID)
        bacteria = tmap['d__Bacteria']
        self.assertEqual(td[bacteria]['parent'], gtdb.CELLULAR_ROOT_TID)
        self.assertEqual(td[gtdb.CELLULAR_ROOT_TID]['parent'], '1')
        self.assertIn('2759', td)   # Eukaryota is grafted
        self.assertNotIn('4930', td)  # deeper NCBI ranks are not pulled in
        self.assertGreaterEqual(stats['n_ncbi_nodes'], 3)

    def test_modern_ncbi_ranks_survive_diamond(self):
        # NCBI's current dumps spell taxId 131567 'cellular root' and use
        # 'clade' above the superkingdoms; DIAMOND 2.1.10 aborts makedb on
        # both, so the hybrid graft must flatten them in nodes.dmp only.
        ncbi = {
            '1': {'parent': '1', 'rank': 'no rank', 'name': 'root'},
            '131567': {'parent': '1', 'rank': 'cellular root',
                       'name': 'cellular organisms'},
            '33154': {'parent': '131567', 'rank': 'clade',
                      'name': 'Opisthokonta'},
            '2759': {'parent': '33154', 'rank': 'superkingdom',
                     'name': 'Eukaryota'}}
        modern = join(self.d, 'ncbi_modern')
        os.makedirs(modern, exist_ok=True)
        with open(join(modern, 'nodes.dmp'), 'w') as f:
            for tid, rec in ncbi.items():
                f.write(f'{tid}\t|\t{rec["parent"]}\t|\t{rec["rank"]}\t|\t'
                        f'\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n')
        with open(join(modern, 'names.dmp'), 'w') as f:
            for tid, rec in ncbi.items():
                f.write(f'{tid}\t|\t{rec["name"]}\t|\t\t|\t'
                        'scientific name\t|\n')
        h = gtdb.HybridProvider(self.gtdb_prov, modern)
        td, _, stats = h.load_taxonomy()
        self.assertEqual(td['131567']['rank'], 'cellular root')
        self.assertEqual(stats['ranks_flattened_for_diamond'],
                         {'cellular root': 1, 'clade': 1})
        # what a compile step would hand to diamond makedb
        _, sub = gtdb.write_compact_taxdump(td, set(td), join(self.d, 'cmp'))
        with open(join(self.d, 'cmp', 'nodes.dmp')) as f:
            written = {line.split('\t|\t')[2] for line in f}
        self.assertTrue(written <= set(gtdb.DIAMOND_RANKS), written)
        self.assertIn('no rank', written)
        self.assertEqual(sub['131567']['rank'], 'cellular root')

    def test_merge_does_not_mutate_the_gtdb_tree(self):
        # HybridProvider re-parents the GTDB domains. Writing that through the
        # tree object the wrapped provider had just cached would turn a
        # "GTDB-only" tree into a hybrid one, depending only on who is asked
        # after the merge -- so the wrapped provider must stay untouched.
        h = self._hybrid()
        merged, _, _ = h.load_taxonomy()
        inner = self.gtdb_prov.taxdump
        bacteria = self.gtdb_prov.taxid_map['d__Bacteria']
        self.assertNotIn(gtdb.CELLULAR_ROOT_TID, inner)
        self.assertLess(len(inner), len(merged))
        self.assertEqual(inner[bacteria]['parent'], gtdb.ROOT_TID)
        self.assertEqual(merged[bacteria]['parent'], gtdb.CELLULAR_ROOT_TID)

    def test_gtdb_ranks_are_never_flattened(self):
        # the fallback must not touch the GTDB side: every GTDB rank has to be
        # DIAMOND-legal after the domain -> superkingdom alias
        self.assertEqual(gtdb.dmp_rank('domain'), 'superkingdom')
        for rank in util.GTDB_RANKS:
            self.assertIn(gtdb.dmp_rank(rank), gtdb.DIAMOND_RANKS)
        self.assertEqual({gtdb.dmp_rank(r) for r in util.GTDB_RANKS}
                         & {gtdb.DMP_RANK_FALLBACK}, set())
        td, _, _ = gtdb.build_taxdump_from_gtdb(self.gtdb_prov.lineages)
        self.assertEqual(gtdb.illegal_dmp_ranks(td), {})

    def test_cross_source_lca_sits_above_domains(self):
        h = self._hybrid()
        td, tmap, _ = h.load_taxonomy()
        # a GTDB domain and a grafted NCBI superkingdom meet at the cellular
        # root, which is exactly the RK7 semantic risk documented in the plan
        self.assertEqual(find_lca([tmap['d__Archaea'], '2759'], td),
                         gtdb.CELLULAR_ROOT_TID)
        # ...and a bacterium with a eukaryote meets one level higher, at root
        self.assertEqual(find_lca([tmap['d__Bacteria'], '2759'], td),
                         gtdb.CELLULAR_ROOT_TID)

    def test_missing_cellular_root_is_rejected(self):
        bad = join(self.d, 'no_cr')
        os.makedirs(bad, exist_ok=True)
        with open(join(bad, 'nodes.dmp'), 'w') as f:
            f.write("1\t|\t1\t|\tno rank\t|\t\t|\t0\t|\t0\t|\t0\t|\t0\t|"
                    "\t0\t|\t0\t|\t0\t|\t\t|\n")
        with open(join(bad, 'names.dmp'), 'w') as f:
            f.write('1\t|\troot\t|\t\t|\tscientific name\t|\n')
        with self.assertRaises(ValueError) as ctx:
            gtdb.HybridProvider(self.gtdb_prov, bad).load_taxonomy()
        self.assertIn('131567', str(ctx.exception))

    def test_lca_policies(self):
        h = self._hybrid('refuse')
        h.load_taxonomy()
        pairs = [('a', str(gtdb.DEFAULT_ID_BASE)), ('b', '4930')]
        self.assertEqual(h.resolve_shared_protein(pairs), [pairs[0]])
        keep_all = self._hybrid('root')
        keep_all.load_taxonomy()
        self.assertEqual(len(keep_all.resolve_shared_protein(pairs)), 2)

    def test_unknown_policy_rejected(self):
        with self.assertRaises(ValueError):
            self._hybrid('nonsense')

    def test_the_three_lca_policies_produce_three_answers(self):
        """Every declared value of --hybrid-lca-policy must assign differently.

        'root' keeps every candidate, 'domain' re-labels a shared protein with
        the donor's own domain node and marks it x_source (plan §5.7 item 3),
        'refuse' drops what is not GTDB-owned. Two of these used to be
        indistinguishable for any taxId the tree knows, and none of them ran
        during a build (see TestSharedProteinResolution).
        """
        out = {}
        for policy in ('root', 'domain', 'refuse'):
            args = cli_args('database', ['-o', 'db', '--source', 'hybrid',
                                         '--gtdb-dir', 'x',
                                         '--hybrid-ncbi-taxdump', 'y',
                                         '--hybrid-lca-policy', policy])
            self.assertEqual(vars(args)['hybrid_lca_policy'], policy)
            h = self._hybrid(policy)
            h.load_taxonomy()
            bac = str(h.taxid_map['d__Bacteria'])
            # p1 is shared between a GTDB genome and an NCBI one, p2 is GTDB
            # only, p3 has only an NCBI candidate that the graft did not pull in
            pairs = [('p1', bac), ('p1', '4930'),
                     ('p2', bac), ('p3', '4930')]
            out[policy] = sorted(h.resolve_shared_protein(pairs))
            if policy == 'domain':
                self.assertIn('p1', h.x_source_proteins,
                              'a cross-source protein must be marked')
                self.assertNotIn('p2', h.x_source_proteins,
                                 'a single-candidate protein is not shared')
        self.assertEqual(out['root'], sorted(pairs))
        self.assertEqual([s for s, _t in out['refuse']], ['p1', 'p2'])
        self.assertEqual([s for s, _t in out['domain']], ['p1', 'p2', 'p3'])
        bac = out['domain'][0][1]
        self.assertEqual({s: t for s, t in out['domain']}['p3'], gtdb.ROOT_TID,
                         'an unplaceable candidate lands at the root, marked '
                         'in the report rather than guessed')
        self.assertEqual({s: t for s, t in out['refuse']}['p1'], bac)
        self.assertNotEqual(out['root'], out['domain'])
        self.assertNotEqual(out['domain'], out['refuse'])


class TestSharedProteinResolution(unittest.TestCase):
    """The build has to consult --hybrid-lca-policy, not just print it.

    extract_genomes_stream let the first taxId met win and counted the rest as
    duplicates, so the policy -- validated at CLI level, echoed into the build
    report -- changed nothing about the database. The resolver now receives the
    colliding candidates and the taxonmap row is rewritten in one streaming
    pass, bounded by the number of *shared* proteins.
    """

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.fetches = []
        for i, (acc, tid) in enumerate([('GB_GCA_000000001.1', '100000001'),
                                        ('GB_GCA_000000002.1', '205')], 1):
            p = join(self.d, f'{acc}.faa')
            with open(p, 'w') as f:
                f.write('>NP_SHARED_1\nMKK\n')
                if i == 2:
                    f.write('>NP_ONLY_2\nMMM\n')
            self.fetches.append(
                gtdb.GenomeFetch(acc, acc, tid, path=p))

    def _run(self, resolve):
        out = join(self.d, f'out{resolve is None}')
        stats = gtdb.extract_genomes_stream(
            self.fetches, join(out, 'db.faa'), join(out, 'taxon.map'),
            'official-plain', resolve=resolve)
        rows = [l.rstrip('\n').split('\t') for l in open(join(out, 'taxon.map'))
                if not l.startswith('#')]
        return stats, {r[0]: r[2] for r in rows}

    def test_without_a_resolver_the_first_taxid_wins(self):
        stats, m = self._run(None)
        self.assertEqual(m['NP_SHARED_1'], '100000001')
        self.assertEqual(stats['proteins_dup'], 1)
        self.assertEqual(stats['proteins_shared'], 0,
                         'no policy ran, so nothing was "resolved"')

    def test_the_policy_rewrites_the_colliding_row(self):
        stats, m = self._run(lambda pairs: [('NP_SHARED_1', pairs[-1][1])])
        self.assertEqual(m['NP_SHARED_1'], '205',
                         'the taxonmap must carry the resolved taxId')
        self.assertEqual(stats['proteins_shared'], 1)
        self.assertEqual(stats['proteins_shared_reassigned'], 1)

    def test_refused_proteins_are_left_unassigned_not_guessed(self):
        stats, m = self._run(lambda pairs: [])
        self.assertEqual(m['NP_SHARED_1'], 'na')
        self.assertEqual(stats['proteins_shared_dropped'], 1)

    def test_only_the_shared_rows_change(self):
        stats, m = self._run(lambda pairs: [pairs[0]])
        self.assertEqual(m['NP_ONLY_2'], '205')
        self.assertEqual(stats['proteins'], 2)
        self.assertEqual(stats['proteins_shared_reassigned'], 0)


class TestResumeState(unittest.TestCase):
    """IT-6: reusing an existing db.faa must reproduce a fresh build."""

    def setUp(self):
        from neoHGT import database
        from types import SimpleNamespace
        self.mod = database
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.out = join(self.d, 'db')
        self.tmp = join(self.d, 'tmp')
        os.makedirs(join(self.out, 'gtdb'))
        os.makedirs(self.tmp)
        self.self = SimpleNamespace(output=self.out, tmpdir=self.tmp)
        with open(join(self.out, 'db.faa'), 'w') as f:
            f.write('>p1\nMAST\n>p2\nMAST\n')
        with gzip.open(join(self.out, 'taxon.map.gz'), 'wt') as f:
            f.write('GB_X_d_Bacteria_p_Cyanobacteriota_1\t1006000002\n')
            f.write('GB_X_d_Bacteria_p_Cyanobacteriota_2\t1006000002\n')
        with gzip.open(join(self.out, 'genome.map.gz'), 'wt') as f:
            f.write('GB_X_d_Bacteria_p_Cyanobacteriota_1\tGB_X\tCP000001.1\n')
        self.ext = {
            'genomes_planned': 1, 'genomes_done': 1, 'genomes_missing': 0,
            'proteins': 2, 'proteins_dup': 0, 'proteins_no_taxid': 0,
            'bytes': 18, 'max_id_len': 35, 'unresolved': ['nope'],
            'stem2taxid': {'GB_X': '1006000002'},
            'prefix2genome': {'CP000001.1': 'GB_X'},
            'per_genome': {'GB_X': (2, 240)}}
        self.fetches = [gtdb.GenomeFetch('GB_X', 'BAC_GCA_000000001',
                                        '1006000002', path=None)]

    def fingerprint(self):
        return self.mod.selection_fingerprint(self.self)

    def taxonmap(self, fetches=None):
        p = join(self.tmp, 'prot.accession2taxid')
        ext = self.mod.resume_extraction_state(self.self, join(
            self.out, 'db.faa'), p, fetches or self.fetches)
        if ext is None:
            return None, []
        return ext, open(p).read().splitlines()

    def test_state_roundtrip_restores_counters_and_maps(self):
        p = self.mod.write_extraction_state(self.out, self.ext, self.fetches,
                                            fingerprint=self.fingerprint())
        self.assertTrue(isfile(p))
        ext, lines = self.taxonmap()
        self.assertTrue(ext['resumed'])
        self.assertEqual(ext['proteins'], 2)
        self.assertEqual(ext['stem2taxid'], {'GB_X': '1006000002'})
        self.assertEqual(ext['per_genome'], {'GB_X': (2, 240)})
        # both key columns of the DIAMOND map must hold the identical full Id
        self.assertEqual(lines[0], gtdb.DIAMOND_TAXONMAP_HEADER.rstrip('\n'))
        self.assertEqual(lines[1],
                         'GB_X_d_Bacteria_p_Cyanobacteriota_1\t'
                         'GB_X_d_Bacteria_p_Cyanobacteriota_1\t'
                         '1006000002\tna')
        self.assertEqual(
            open(join(self.tmp, 'gtdb_seq2taxid.txt')).read().splitlines()[1],
            'GB_X_d_Bacteria_p_Cyanobacteriota_1\t1006000002')
        self.assertTrue(isfile(join(self.tmp, 'gtdb_genome.map')))

    def test_unresolved_list_is_capped(self):
        self.ext['unresolved'] = [str(i) for i in range(5000)]
        p = self.mod.write_extraction_state(self.out, self.ext, self.fetches,
                                            fingerprint=self.fingerprint())
        self.assertEqual(len(json.load(open(p))['unresolved']), 2000)

    def test_refuses_changed_selection(self):
        self.mod.write_extraction_state(self.out, self.ext, self.fetches,
                                        fingerprint=self.fingerprint())
        other = [gtdb.GenomeFetch('GB_Y', 'BAC_GCA_000000002',
                                 '1006000003', path=None)]
        ext, lines = self.taxonmap(fetches=other)
        self.assertIsNone(ext)
        self.assertEqual(lines, [])

    def test_refuses_without_state(self):
        self.mod.write_extraction_state(self.out, self.ext, self.fetches,
                                        fingerprint=self.fingerprint())
        os.remove(join(self.out, 'db.faa'))
        ext, _ = self.taxonmap()
        self.assertIsNone(ext)

    def test_refuses_truncated_db_faa(self):
        self.mod.write_extraction_state(self.out, self.ext, self.fetches,
                                        fingerprint=self.fingerprint())
        with open(join(self.out, 'db.faa'), 'w') as f:
            f.write('>p1\nMAS')          # 9 of the recorded 18 bytes
        ext, _ = self.taxonmap()
        self.assertIsNone(ext, 'a truncated db.faa must not be resumed')

    def test_refuses_without_taxonmap(self):
        self.mod.write_extraction_state(self.out, self.ext, self.fetches,
                                        fingerprint=self.fingerprint())
        os.remove(join(self.out, 'taxon.map.gz'))
        with self.assertRaises(ValueError) as ctx:
            self.taxonmap()
        self.assertIn('taxon.map.gz', str(ctx.exception))


class TestArchaeaProfile(unittest.TestCase):
    """UT-4.x: the archaeal close-group default must be picked by domain."""

    TAXDUMP = {
        '1': {'parent': '1', 'rank': 'no rank', 'name': 'root'},
        '2': {'parent': '1', 'rank': 'superkingdom', 'name': 'Archaea'},
        '21': {'parent': '1', 'rank': 'superkingdom', 'name': 'Bacteria'},
        '20': {'parent': '2', 'rank': 'phylum', 'name': 'Halobacteriota'},
        '200': {'parent': '20', 'rank': 'class', 'name': 'Halobacteria'},
        '2000': {'parent': '200', 'rank': 'order', 'name': 'Halobacteriales'},
        '20000': {'parent': '2000', 'rank': 'family', 'name': 'Halorubraceae'},
        '200000': {'parent': '20000', 'rank': 'genus', 'name': 'Halorubrum'}}

    def _probe(self, input_tax, lca):
        from types import SimpleNamespace
        from neoHGT.analyze import domain_of_input
        return domain_of_input(SimpleNamespace(
            input_tax=input_tax, data={}, taxdump=self.TAXDUMP))

    def _hits(self, tids):
        import pandas as pd
        from types import SimpleNamespace
        from neoHGT.analyze import domain_of_input
        prots = [{'hits': pd.DataFrame({'taxid': [t]})} for t in tids]
        return domain_of_input(SimpleNamespace(
            input_tax=None, data={'q': prots}, taxdump=self.TAXDUMP))

    def test_domain_of_input(self):
        self.assertEqual(self._probe({'p': '200000'}, '200000'), 'archaea')
        self.assertEqual(self._probe({'p': '21'}, '21'), 'bacteria')
        self.assertIsNone(self._probe({}, None))
        self.assertIsNone(self._probe(None, None))
        self.assertIsNone(self._probe({'p': '999999'}, None))

    def test_domain_from_best_hits(self):
        self.assertEqual(self._hits(['200000'] * 3 + ['21']), 'archaea')
        self.assertEqual(self._hits(['21', '21', '21', '200000']), 'bacteria')
        self.assertIsNone(self._hits([]))

    def test_profiles_declared_in_config(self):
        import yaml
        cfg = yaml.safe_load(open(join(HERE, '..', 'config.yml')))
        prof = cfg['grouping']['profiles']
        self.assertEqual(prof['gtdb']['closesize'], 200)
        # r220 archaeal representatives are 18x rarer than bacterial (RK16)
        self.assertLess(prof['gtdb_archaea']['closesize'],
                        prof['gtdb']['closesize'])


class TestConfigDefaults(unittest.TestCase):
    """config.yml database.gtdb supplies defaults; the CLI always wins."""

    def _params(self, **over):
        import argparse
        import yaml
        from neoHGT.database import Database, set_source_parameters
        d = Database.__new__(Database)
        d.cfg = yaml.safe_load(open(join(HERE, '..', 'config.yml')))
        d.output = join(os.environ['NEOHGT_CACHE'], 'db', '_never_created_')
        d.input = None
        d.gtdb_taxa = over.get('gtdb_taxa')
        a = argparse.Namespace(source='gtdb', **over)
        # the shipped config carries no machine paths any more, so the dataset
        # root comes from $GTDB_DIR: provide it for the duration of the call
        # unless the test names one explicitly
        saved = os.environ.get('GTDB_DIR')
        if not getattr(a, 'gtdb_dir', None):
            # a real directory: set_source_parameters validates isdir()
            os.environ['GTDB_DIR'] = tmpdir()
            self.addCleanup(shutil.rmtree, os.environ['GTDB_DIR'],
                            ignore_errors=True)
            if saved:
                self.addCleanup(os.environ.__setitem__, 'GTDB_DIR', saved)
            else:
                self.addCleanup(os.environ.pop, 'GTDB_DIR', None)
        try:
            set_source_parameters(d, a)
        finally:
            if saved is not None and 'GTDB_DIR' not in os.environ:
                os.environ['GTDB_DIR'] = saved
        return d

    def test_config_supplies_omitted_options(self):
        d = self._params()
        self.assertEqual(d.seq_source, 'local')
        self.assertEqual(d.local_layout, 'auto')
        self.assertEqual(d.id_base, gtdb.DEFAULT_ID_BASE)
        self.assertTrue(d.fingerprints)      # verify_input_fingerprints: yes
        self.assertFalse(d.no_selfcheck)     # selfcheck: required
        self.assertIsNone(d.gtdb_release)    # release: auto -> infer

    def test_command_line_overrides_config(self):
        d = self._params(seq_source='ena', local_layout='official-gz',
                         gtdb_taxid_base=200000000, input_fingerprints='no',
                         gtdb_no_selfcheck=True, gtdb_release='232')
        self.assertEqual(d.seq_source, 'ena')
        self.assertEqual(d.local_layout, 'official-gz')
        self.assertEqual(d.id_base, 200000000)
        self.assertFalse(d.fingerprints)
        self.assertTrue(d.no_selfcheck)
        self.assertEqual(d.gtdb_release, '232')

    def test_low_id_base_rejected(self):
        with self.assertRaises(ValueError):
            self._params(gtdb_taxid_base=9999999)

    def test_nothing_written(self):
        self._params()
        self.assertFalse(exists(join(os.environ['NEOHGT_CACHE'], 'db',
                                     '_never_created_')))


    def test_the_shipped_config_is_the_one_that_loads_from_any_cwd(self):
        """util.find_config_file() tries ./config.yml first, so a stray file in
        the working directory would silently replace the configuration that all
        the GTDB profile defaults live in."""
        cwd0 = os.getcwd()
        other = tmpdir()
        self.addCleanup(shutil.rmtree, other, ignore_errors=True)
        try:
            os.chdir(other)
            path = util.find_config_file()
            self.assertTrue(path.endswith(join('neoHGT', 'config.yml')), path)
            self.assertEqual(os.path.realpath(os.path.dirname(path)),
                             os.path.realpath(join(HERE, os.pardir)))
            cfg = util.load_configs()
            self.assertIn('workspace', cfg)
            self.assertIn('gtdb', ((cfg.get('taxonomy') or {})
                                   .get('profiles') or {}))
            self.assertIn('gtdb_archaea', ((cfg.get('grouping') or {})
                                           .get('profiles') or {}))
        finally:
            os.chdir(cwd0)


class TestGtdbtkInput(unittest.TestCase):
    """P7: a GTDB-Tk summary must become --input-tax taxIds of the GTDB tree."""

    LIN = ('d__Archaea;p__Halobacteriota;c__Halobacteria;o__Halobacteriales;'
           'f__Haloferacaceae;g__Halorubrum;s__Halorubrum vacuolatum')
    PARTIAL = ('d__Archaea;p__Halobacteriota;c__Halobacteria;o__Halobacteriales;'
               'f__Haloferacaceae;g__Halorubrum;NA')

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        db = join(self.d, 'db')
        os.makedirs(join(db, 'gtdb'))
        self.db = db
        segs = self.LIN.split(';')
        with gzip.open(join(db, 'gtdb', 'taxid.map.gz'), 'wt') as f:
            f.write('# manifest\nlineage\ttaxid\trank\tlabel\tname\t'
                    'cluster_letter\tparent_taxid\tncbi_taxid_hint\n')
            for i in range(1, len(segs) + 1):
                f.write(';'.join(segs[:i]) + f'\t{i * 100}\t'
                        f'{gtdb.GTDB_RANKS[i - 1]}\t{segs[i - 1]}\tx\t\t\t\n')

    def _summary(self, text, name='tk.summary'):
        p = join(self.d, name)
        with open(p, 'w') as f:
            f.write(text)
        return p

    def _run(self, path, **kw):
        out = join(self.d, 'input_tax.tsv')
        return gtdb.gtdbtk_to_input_tax(path, self.db, out, **kw), out

    def test_classify_info_and_partial_lineage(self):
        p = self._summary(
            '#user_genome\tclassification\tclassify_info\n'
            f'BIN1\t{self.LIN}\t{self.LIN}\n'
            f'BIN2\t\t{self.PARTIAL}\n')
        rep, out = self._run(p)
        self.assertEqual(rep['n_mapped'], 2)
        self.assertEqual(rep['matched_rank']['species'], 1)
        self.assertEqual(rep['matched_rank']['genus'], 1)
        self.assertEqual(open(out).read().splitlines(), ['BIN1\t700',
                                                         'BIN2\t600'])

    def test_v2_banner_lines_are_skipped(self):
        p = self._summary(
            '##### GTDB-Tk v2.3.2 (gtdbtkv1)\n#####\n'
            '#user_genome\tclassification\n'
            f'BIN4\t{self.LIN}\n', 'max.summary')
        rep, _ = self._run(p)
        self.assertEqual(rep['n_mapped'], 1)
        self.assertEqual(rep['columns']['genome_col'], 'user_genome')

    def test_rank_columns_shape(self):
        p = self._summary(
            'user_genome\td__Domain\tp__Phylum\tc__Class\to__Order\t'
            'f__Family\tg__Genus\ts__Species\n'
            'BIN3\td__Archaea\tp__Halobacteriota\tc__Halobacteria\t'
            'o__Halobacteriales\tf__Haloferacaceae\tg__Halorubrum\t'
            's__Halorubrum vacuolatum\n', 'denovo.summary')
        rep, _ = self._run(p)
        self.assertEqual(rep['n_mapped'], 1)
        self.assertEqual(rep['columns']['lineage_col'], '<rank columns>')

    def test_min_rank_gate_and_unknown(self):
        p = self._summary(
            'user_genome\tclassify_info\n'
            f'B_A\t{self.PARTIAL}\n'
            'B_D\td__Archaea;NA;NA;NA;NA;NA;NA\n'
            'B_Z\td__Bacteria;p__X;c__Y;o__Z;f__Q;g__W;s__W sp9\n')
        rep, _ = self._run(p, min_rank='family')
        self.assertEqual(rep['n_mapped'], 1)
        reasons = {u['genome']: u['reason'] for u in rep['unmapped']}
        self.assertIn('below --min-rank family', reasons['B_D'])
        self.assertEqual(reasons['B_Z'], 'lineage not in this database')

    def test_missing_manifest_rejected(self):
        empty = join(self.d, 'empty_db')
        os.makedirs(join(empty, 'gtdb'))
        with self.assertRaises(ValueError) as ctx:
            gtdb.gtdbtk_to_input_tax(self._summary('user_genome\nc\n'),
                                     empty, join(self.d, 'o.tsv'))
        self.assertIn('taxid.map.gz', str(ctx.exception))


class TestMemberExpansion(unittest.TestCase):
    """P7-2: representatives can be expanded through the metadata column."""

    def setUp(self):
        import pandas as pd
        self.df = pd.DataFrame({
            'accession': ['RS_A', 'RS_B', 'GB_a1', 'GB_a2', 'GB_b1'],
            'genome': ['GCA_1', 'GCA_2', 'GCA_11', 'GCA_12', 'GCA_21'],
            'gtdb_taxonomy': ['d__Bacteria;p__X;c__Y;o__Z;f__F;g__G;s__G g1'] * 5,
            'gtdb_genome_representative':
                ['RS_A', 'RS_B', 'RS_A', 'RS_A', 'RS_B'],
            'checkm2_completeness': [90.0, 91.0, 88.0, 89.0, 92.0],
            'checkm2_contamination': [1.0, 2.0, 3.0, 4.0, 5.0]})
        self.sel = self.df[self.df['accession'].isin(['RS_A'])]

    def test_member_index(self):
        idx = gtdb.member_index(self.df)
        self.assertEqual(sorted(idx['RS_A']), ['GB_a1', 'GB_a2', 'RS_A'])
        self.assertEqual(idx['RS_B'], ['RS_B', 'GB_b1'])

    def test_expand_lists_members_and_locality(self):
        index = {'GB_a1': '/x/GB_a1.faa'}          # only one member on disk
        members, rep = gtdb.expand_to_members(self.df, self.sel,
                                              seq_index=index)
        self.assertEqual(sorted(members['accession']), ['GB_a1', 'GB_a2'])
        self.assertEqual(rep['members_in_metadata'], 2)
        self.assertEqual(rep['members_with_sequence'], 1)
        self.assertEqual(list(members['representative']), ['RS_A', 'RS_A'])
        self.assertEqual(sorted(members.loc[members['has_local_fasta'],
                                            'accession']), ['GB_a1'])

    def test_expand_without_sequence_index(self):
        members, rep = gtdb.expand_to_members(self.df, self.sel)
        self.assertNotIn('has_local_fasta', members.columns)
        self.assertEqual(rep['members_in_metadata'], 2)
        self.assertNotIn('members_with_sequence', rep)

    def test_no_members_returns_empty(self):
        sel = self.df[self.df['accession'] == 'GB_a1']   # a member, not a rep
        members, rep = gtdb.expand_to_members(self.df, sel, seq_index={})
        self.assertEqual(members.shape[0], 0)
        self.assertEqual(rep['members'], 0)

    def test_missing_column_rejected(self):
        with self.assertRaises(ValueError) as ctx:
            gtdb.member_index(self.df.drop(
                columns=['gtdb_genome_representative']))
        self.assertIn('gtdb_genome_representative', str(ctx.exception))


class TestPhyloDistances(unittest.TestCase):
    """P7: patristic distances from the per-rank trees group the close set."""

    NEWICK = '((A:0.1,B:0.2)AB:0.3,(C:0.4,D:0.5)CD:0.6)root;'

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.tree = join(self.d, 'x_genus.tree')
        with open(self.tree, 'w') as f:
            f.write(self.NEWICK)

    def test_structure(self):
        t = gtdb.NewickTree(self.NEWICK)
        self.assertEqual(sorted(t.name[n] for n in t.leaves()),
                         ['A', 'B', 'C', 'D'])
        self.assertEqual(t.branch[t.by_name()['A']], 0.1)
        self.assertEqual(t.by_name()['AB'] in t.children, True)

    def test_distances(self):
        t = gtdb.NewickTree(self.NEWICK)
        n = t.by_name()
        self.assertAlmostEqual(t.distance(n['A'], n['B']), 0.3, places=9)
        self.assertAlmostEqual(t.distance(n['A'], n['C']), 1.4, places=9)
        self.assertAlmostEqual(t.distance(n['AB'], n['CD']), 0.9, places=9)
        got = t.distances_from(n['A'])
        self.assertNotIn('A', got)
        self.assertAlmostEqual(got['D'], 1.5, places=9)

    def test_close_taxa_with_cutoffs(self):
        rows, s = gtdb.tree_close_taxa(self.tree, 'A', top_k=2)
        self.assertEqual([r[0] for r in rows], ['B', 'C'])
        self.assertEqual(s['n_taxa'], 3)
        self.assertAlmostEqual(s['within']['0.5'], 1, places=9)
        rows, _ = gtdb.tree_close_taxa(self.tree, 'A', max_distance=0.35)
        self.assertEqual([r[0] for r in rows], ['B'])

    def test_unknown_focus_lists_examples(self):
        with self.assertRaises(ValueError) as ctx:
            gtdb.tree_close_taxa(self.tree, 'g__Nope')
        self.assertIn('not a labelled node', str(ctx.exception))

    def test_rank_trees_discovery(self):
        root = join(self.d, 'rel')
        for sub, fname in (('tree_bac120', 'bac120_r220_genus.tree'),
                           ('tree_ar53', 'ar53_r220_genus.tree'),
                           ('tree_ar53', 'ar53_r220.tree'),
                           ('tree_ar53', 'ar53_taxonomy_r220.tsv')):
            os.makedirs(join(root, sub), exist_ok=True)
            open(join(root, sub, fname), 'w').write(self.NEWICK)
        found = gtdb.rank_trees(root)
        self.assertEqual(sorted(found['genus']),
                         [join(root, 'tree_ar53/ar53_r220_genus.tree'),
                          join(root, 'tree_bac120/bac120_r220_genus.tree')])
        self.assertEqual(found['all'],
                         [join(root, 'tree_ar53/ar53_r220.tree')])
        self.assertNotIn('tsv', str(found))

    def test_malformed_newick(self):
        with self.assertRaises(ValueError):
            gtdb.NewickTree('((A:0.1,B:0.2)')
        bad = join(self.d, 'bad.tree')
        open(bad, 'w').write('')
        with self.assertRaises(ValueError):
            gtdb.read_newick(bad)


    def test_cli_with_db_resolution(self):
        # a synthetic release + a synthetic database whose labels map back
        rel = join(self.d, 'rel')
        os.makedirs(join(rel, 'tree_demo'))
        open(join(rel, 'tree_demo', 'demo_r1_genus.tree'), 'w').write(
            '((A:0.1,B:0.2)AB:0.3,(C:0.4,D:0.5)CD:0.6)root;')
        db = join(self.d, 'db')
        os.makedirs(join(db, 'gtdb'))
        with open(join(db, 'gtdb', 'labels.tsv'), 'w') as f:
            f.write('taxid\tlabel\trank\tname\tcluster_letter\tn_genomes\n')
            for tid, lab in (('101', 'A'), ('102', 'B'), ('103', 'D')):
                f.write(f'{tid}\t{lab}\tgenus\t{lab[1:]}\t\t1\n')
        out = join(self.d, 'close.tsv')
        lst = join(self.d, 'close.txt')
        rep = join(self.d, 'close.json')
        rc = gtdb.main(['tree-dist', '--gtdb-dir', rel, '--rank', 'genus',
                        '--focus', 'A', '--max-distance', '1.5', '--db', db,
                        '--out', out, '--close-tax-out', lst,
                        '--report', rep])
        self.assertEqual(rc, 0)
        lines = open(out).read().splitlines()
        self.assertEqual(lines[0], 'label\tdistance\ttaxid')
        self.assertEqual(lines[1].split('\t')[0], 'B')
        self.assertEqual(open(lst).read().strip(), '102,103')
        summary = json.load(open(rep))
        self.assertEqual(summary['close_tax_n'], 2)
        self.assertEqual(summary['taxids_unmapped'], 1)          # C has no row
        near = summary['nearest'][0]
        self.assertEqual([near[0], near[2]], ['B', '102'])
        self.assertAlmostEqual(near[1], 0.3, places=9)

    def test_cli_requires_tree_or_dir(self):
        with self.assertRaises(ValueError):
            gtdb.main(['tree-dist', '--focus', 'A'])


class TestCloseGroupReachability(unittest.TestCase):
    """A phylogeny-derived close group must survive taxonomy refinement."""

    def setUp(self):
        import pandas as pd
        from types import SimpleNamespace
        from neoHGT.analyze import _keep_visitable_close_taxa
        from neoHGT.util import add_children
        self.fn = _keep_visitable_close_taxa
        self.taxdump = {
            '1': {'parent': '1', 'rank': 'no rank', 'name': 'root'},
            '2': {'parent': '1', 'rank': 'superkingdom', 'name': 'Archaea'},
            '20': {'parent': '2', 'rank': 'phylum', 'name': 'Halobacteriota'},
            '200': {'parent': '20', 'rank': 'class', 'name': 'Halobacteria'},
            '2000': {'parent': '200', 'rank': 'order', 'name': 'Halcanales'},
            '20000': {'parent': '2000', 'rank': 'family', 'name': 'Halcaceae'},
            'g1': {'parent': '20000', 'rank': 'genus', 'name': 'Hitgenus'},
            's1': {'parent': 'g1', 'rank': 'species', 'name': 'Hitgenus sp1'},
            'g2': {'parent': '20000', 'rank': 'genus', 'name': 'Emptygenus'}}
        add_children(self.taxdump)
        hits = pd.DataFrame({'taxid': ['s1', 's1']})
        self.me = SimpleNamespace(
            data={'q': [{'hits': hits}, {'hits': hits}]}, taxdump=self.taxdump,
            gtdb_labels={'g1': 'g__Hitgenus', 'g2': 'g__Emptygenus'})

    def test_keeps_taxa_reached_by_hits(self):
        # g1 is an ancestor of the present s1, so it stays; the class node too
        out = self.fn(self.me, ['g1', '200'])
        self.assertEqual(out, ['g1', '200'])

    def test_drops_unreachable_taxa_with_warning(self):
        out = self.fn(self.me, ['g1', 'g2'])
        self.assertEqual(out, ['g1'])

    def test_all_unreachable_raises(self):
        from types import SimpleNamespace
        import pandas as pd
        me = SimpleNamespace(
            data={'q': [{'hits': pd.DataFrame({'taxid': ['s1']})}]},
            taxdump=self.taxdump, gtdb_labels={})
        with self.assertRaises(ValueError) as ctx:
            self.fn(me, ['g2'])
        self.assertIn('--close-tax', str(ctx.exception))

    def test_no_hits_passes_through(self):
        from types import SimpleNamespace
        me = SimpleNamespace(data={}, taxdump=self.taxdump, gtdb_labels={})
        self.assertEqual(self.fn(me, ['g2']), ['g2'])


class TestTaxmapIndex(unittest.TestCase):
    """P7 route B: --taxmap must work at release scale without loading it."""

    N = 5000

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.path = join(self.d, 'map.tsv')
        with open(self.path, 'w') as f:
            for i in range(self.N):
                f.write(f'SP_{i:08d}.1\t{2000 + i}\n')

    def test_lookups(self):
        ix = gtdb.TaxidIndex(self.path, spacing=512)
        self.addCleanup(ix.close)
        self.assertEqual(len(ix), self.N)
        self.assertEqual(ix.get('SP_00000000.1'), '2000')       # first anchor
        self.assertEqual(ix.get('SP_00000002.1'), '2002')       # after anchor
        self.assertEqual(ix.get(f'SP_{self.N - 1:08d}.1'),
                         str(2000 + self.N - 1))                # last line
        self.assertIsNone(ix.get('SP_99999999.1'))              # past the end
        self.assertIsNone(ix.get('no-such-id'))
        self.assertIn('SP_00000123.1', ix)
        self.assertNotIn('SP_00000123.2', ix)
        with self.assertRaises(KeyError):
            ix['no-such-id']

    def test_index_is_small_and_cached(self):
        ix = gtdb.TaxidIndex(self.path, spacing=512)
        anchors = len(ix.keys)
        ix.close()
        self.assertLessEqual(anchors, self.N / 512 + 1)
        self.assertTrue(isfile(self.path + '.idx'))
        t0 = time.time()
        again = gtdb.TaxidIndex(self.path, spacing=512)
        again.close()
        self.assertEqual(len(again.keys), anchors)
        self.assertLess(time.time() - t0, 0.5)     # reuses the cached index

    def test_rejects_gzipped_and_unsorted(self):
        gz = self.path + '.gz'
        with gzip.open(gz, 'wt') as f:
            f.write('A\t1\n')
        with self.assertRaises(ValueError) as ctx:
            gtdb.TaxidIndex(gz)
        self.assertIn('Decompress', str(ctx.exception))
        bad = join(self.d, 'bad.tsv')
        with open(bad, 'w') as f:
            f.write('B\t1\nA\t2\n')
        with self.assertRaises(ValueError) as ctx:
            gtdb.TaxidIndex(bad)
        self.assertIn('not sorted', str(ctx.exception))

    def test_search_wrapper_caches_and_validates(self):
        from types import SimpleNamespace
        from neoHGT.search import open_taxmap_index
        me = SimpleNamespace(taxmap=self.path)
        idx = open_taxmap_index(me)
        self.assertIsNotNone(idx)
        self.assertIs(open_taxmap_index(me), idx)     # cached on the instance
        idx.close()
        self.assertIsNone(open_taxmap_index(SimpleNamespace(taxmap=None)))
        bad_gz = join(self.d, 'wrapper.tsv.gz')
        with gzip.open(bad_gz, 'wt') as f:
            f.write('A\t1\n')
        me2 = SimpleNamespace(taxmap=bad_gz)
        with self.assertRaises(ValueError) as ctx:
            open_taxmap_index(me2)
        self.assertIn('--taxmap could not be indexed', str(ctx.exception))

    def test_route_b_after_route_c(self):
        # A hit table without staxids whose subject ids cannot be parsed by
        # route C must be rescued by --taxmap read through the disk index.
        db = join(self.d, 'db')
        os.makedirs(join(db, 'gtdb'))
        with gzip.open(join(db, 'gtdb', 'stem2taxid.tsv.gz'), 'wt') as f:
            f.write('stem\ttaxid\n')          # nothing resolvable by route C
        map_path = join(self.d, 'map.tsv')
        with open(map_path, 'w') as f:
            f.write('SP_00000000.1\t2000\nSP_00000001.1\t2001\n')
        s = Search.__new__(Search)
        s.source, s.gtdb_root = 'gtdb', db
        s.gtdb_layout = 'gtdb-r220-renamed'
        s.gtdb_stems, s.gtdb_p2g, s.gtdb_g2t = {}, {}, {}
        s.prot2tid, s.method = {}, 'precomp'
        s.taxmap, s.keep_unresolved = map_path, 'no'
        s.tmpdir = s.output = self.d
        s.n_unresolved = 0
        prots = {'q1': [{'id': 'SP_00000001.1', 'taxid': '', 'identity': 90.0,
                         'evalue': 1e-40, 'score': 300.0, 'coverage': 95.0}]}
        taxid_wf_gtdb(s, prots)
        self.assertEqual(prots['q1'][0]['taxid'], '2001')
        self.assertEqual(s.n_unresolved, 0)
        self.assertFalse(isfile(join(self.d, 'unresolved_taxids.txt')))


class TestCapacityPreflight(unittest.TestCase):
    """UT-5.x: the build must stop before writing when the volume is too small.

    The estimate uses the measured 466 B/protein constant and the 3x peak
    factor, so the thresholds below are derived from those constants rather
    than hardcoded bytes.
    """

    def setUp(self):
        import pandas as pd
        from types import SimpleNamespace
        from unittest import mock
        from neoHGT import database
        self.mock = mock
        self.db = database
        self.df = pd.DataFrame({'proteins': [1000]})
        self.need = (1000 * self.db.GTDB_BYTES_PER_PROTEIN
                     * self.db.GTDB_PEAK_FACTOR)
        self.Self = lambda force=False: SimpleNamespace(
            output='/nowhere', force=force)

    def _free(self, nbytes, force=False):
        usage = type('U', (), {'free': nbytes})
        with self.mock.patch.object(self.db._shutil, 'disk_usage',
                                    lambda path: usage):
            return self.db.check_capacity(self.Self(force), self.df)

    def test_roomy_volume_passes(self):
        res = self._free(self.need * 10)
        self.assertEqual(res['estimated_bytes'],
                         1000 * self.db.GTDB_BYTES_PER_PROTEIN)
        self.assertEqual(res['plan_proteins'], 1000)

    def test_tight_volume_stops_with_an_actionable_message(self):
        with self.assertRaises(ValueError) as ctx:
            self._free(self.need // 2)
        msg = str(ctx.exception)
        for hint in ('--rank-sample', '--keep-faa', '--force'):
            self.assertIn(hint, msg)
        self.assertIn('Insufficient free space', msg)

    def test_force_downgrades_to_a_warning(self):
        res = self._free(1024, force=True)      # must not raise
        self.assertLess(res['free_bytes'], res['estimated_bytes'])

    def test_boundary_exactly_at_the_peak_passes(self):
        self.assertEqual(self._free(self.need)['free_bytes'], self.need)


class TestCloseGroupSizing(unittest.TestCase):
    """UT-4.x: `infer_close_group` must reach `close_size` on a balanced tree.

    A balanced synthetic tree makes the expectation exact: the walk climbs from
    the self-group LCA until the clade (minus the self group) is large enough,
    and stops there.
    """

    def setUp(self):
        from neoHGT.analyze import Analyze
        from neoHGT.util import add_children
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        seq, self.td = [0], {'1': {'parent': '1', 'rank': 'no rank',
                                   'name': 'root'}}

        def node(parent, rank, name):
            seq[0] += 1
            x = str(1000 + seq[0])
            self.td[x] = {'parent': parent, 'rank': rank, 'name': name}
            return x

        self.dom = node('1', 'domain', 'd__X')
        self.genuses = [node(self.dom, 'genus', f'g__{i}') for i in range(4)]
        self.species = {i: [node(g, 'species', f's__{i}_{j}')
                            for j in range(10)]
                        for i, g in enumerate(self.genuses)}
        add_children(self.td)
        self.an = Analyze.__new__(Analyze)
        self.an.taxdump = self.td
        self.an.groups = {'self': set(), 'close': set(), 'distal': set()}

    def _run(self, self_tax, close_size):
        self.an.self_tax = set(self_tax)
        self.an.groups['self'] = set(self_tax)
        self.an.close_size = close_size
        self.an.infer_close_group()
        return self.an.groups['close'], self.an.close_tax

    def test_climbs_out_of_the_genus_until_close_size_is_met(self):
        focus = self.species[0][:3]
        # one genus holds 10 species + itself = 8 candidates after excluding the
        # self group, which is below 25, so the walk must rise to the domain
        close, close_tax = self._run(focus, 25)
        self.assertEqual(close_tax, [self.dom])
        self.assertGreaterEqual(len(close), 25)
        self.assertFalse(set(focus) & close)

    def test_small_close_size_stops_at_the_own_genus(self):
        focus = self.species[0][:4]
        close, close_tax = self._run(focus, 4)
        self.assertEqual(close_tax, [self.genuses[0]])
        self.assertEqual(close, {self.genuses[0], *self.species[0][4:]})

    def test_two_genera_escalate_to_the_domain_at_once(self):
        focus = [self.species[0][0], self.species[1][0]]
        close, close_tax = self._run(focus, None)
        self.assertEqual(close_tax, [self.dom])
        self.assertEqual(len(close), 4 + 40 - 2 + 1)   # genera + species + dom

    def test_distal_group_is_what_is_left(self):
        focus = self.species[0][:3]
        close, _ = self._run(focus, 8)
        self.assertEqual(len(close), 8)


class TestRiskMitigations(unittest.TestCase):
    """RK6 / RK13 / RK8: the mitigations promised in the risk table, tested.

    Each case drives the mitigation itself, not a comment claiming it exists.
    """

    ENT = ('d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;'
           'o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;'
           's__Escherichia coli')
    SALM = ('d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;'
            'o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;'
            's__Salmonella enterica')
    ARCH = ('d__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Nanoarchaeales;'
            'f__Nanoarchaeaceae;g__Nanoarchaeum;s__Nanoarchaeum equitans')
    RENAMED = ('d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;'
               'o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;'
               's__Escherichia fenusiae')      # new species, genus kept
    SALM_GENUS = SALM.rsplit(';', 1)[0]

    # --- RK6: cross-release Id stability -------------------------------
    def test_ids_are_reused_and_dropped_taxa_are_recorded(self):
        base = gtdb.DEFAULT_ID_BASE
        td_a, map_a, st_a = gtdb.build_taxdump_from_gtdb(
            [self.ENT, self.SALM, self.ARCH], id_base=base)
        manifest = {k: {'taxid': v} for k, v in map_a.items()}

        # release B: one species gone, another added, rest untouched
        td_b, map_b, st_b = gtdb.build_taxdump_from_gtdb(
            [self.ENT, self.RENAMED, self.ARCH], manifest=manifest,
            id_base=base)
        for key in (self.ENT,
                    'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;'
                    'o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia',
                    'd__Bacteria', self.ARCH):
            self.assertEqual(map_a[key], map_b[key],
                             f'{key} must keep its taxId across releases')
        self.assertNotIn(self.SALM, map_b)
        self.assertIn(self.RENAMED, map_b)
        # the new node may not steal the id of a node that just disappeared
        self.assertNotEqual(map_b[self.RENAMED], map_a[self.SALM])
        superseded = st_b.get('superseded') or {}
        self.assertIn(self.SALM, superseded)
        self.assertEqual(superseded[self.SALM], map_a[self.SALM])

    def test_delnodes_file_lists_the_superseded_nodes(self):
        import pandas as pd
        base = gtdb.DEFAULT_ID_BASE
        out = tmpdir()
        self.addCleanup(shutil.rmtree, out, ignore_errors=True)
        os.makedirs(join(out, 'gtdb'), exist_ok=True)
        td_a, map_a, _ = gtdb.build_taxdump_from_gtdb(
            [self.ENT, self.SALM], id_base=base)
        manifest = {k: {'taxid': v} for k, v in map_a.items()}
        td_b, map_b, st_b = gtdb.build_taxdump_from_gtdb(
            [self.ENT], manifest=manifest, id_base=base)
        df = pd.DataFrame([{
            'accession': 'RS_GCF_000026225.1', 'genome': 'RS_GCF_000026225.1',
            'gtdb_taxonomy': self.ENT, 'gtdb_species': self.ENT,
            'gtdb_representative': 't', 'is_rep': True,
            'domain': 'bacteria', 'proteins': 4}])
        stats = dict(st_b)
        stats.setdefault('id_base', base)
        gtdb.write_side_files(out, {'release': '220', 'layout': 'auto',
                                    'gtdb_dir': join(out, 'src'),
                                    'metadata': [],
                                    'taxonomy': [], 'domains': {},
                                    'idlists': {}}, df, map_b, td_b, stats)
        path = join(out, 'gtdb', 'delnodes.txt')
        self.assertTrue(isfile(path))
        rows = {tuple(l.rstrip('\n').split('\t'))
                for l in open(path) if not l.startswith('#')
                and not l.startswith('accession')}
        # dropping the Salmonella species also drops its genus node, so both
        # historical ids must be reserved and listed
        self.assertEqual(rows, {
            (map_a[self.SALM], self.SALM, 'superseded_by=removed'),
            (map_a[self.SALM_GENUS], self.SALM_GENUS,
             'superseded_by=removed')})
        for tid, _path, _tag in rows:
            self.assertNotIn(tid, set(map_b.values()))

    def test_reserved_ids_carry_across_chained_rebuilds(self):
        """RK6, second hop: a sidecar id must not be handed out again."""
        base = gtdb.DEFAULT_ID_BASE
        n1 = gtdb.collect_nodes([self.ENT, self.SALM])
        m1 = {}
        ids1, sup1 = gtdb.assign_ids(n1, m1, base)
        n2 = gtdb.collect_nodes([self.ENT])                    # Salmonella gone
        man2 = {k: {'taxid': v} for k, v in ids1.items()}
        ids2, sup2 = gtdb.assign_ids(n2, man2, base)
        freed = sup2[self.SALM]
        n3 = gtdb.collect_nodes([self.ENT, self.RENAMED])
        man3 = {k: {'taxid': v} for k, v in ids2.items()}
        leaked = gtdb.assign_ids(n3, man3, base)[0][self.RENAMED]
        held = gtdb.assign_ids(n3, man3, base, reserve={freed})[0][self.RENAMED]
        self.assertEqual(leaked, freed)       # demonstrates the hazard
        self.assertNotEqual(held, freed)      # reserving closes it

    def test_sidecar_round_trip(self):
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'gtdb'), exist_ok=True)
        manifest = join(d, 'gtdb', 'taxid.map.gz')
        self.assertEqual(gtdb.load_reserved(manifest), {})
        gtdb.write_reserved(join(d, 'gtdb', 'reserved_ids.txt'),
                            [('106000002', self.SALM), ('106000007', '')],
                            {'release': '220', 'layout': 'auto'})
        self.assertEqual(gtdb.load_reserved(manifest),
                         {'106000002': self.SALM, '106000007': ''})

    # --- RK13: the expired AAU certificate must not sink a build -------
    def test_mirror_resolution_prefers_the_verifiable_host(self):
        self.assertEqual(gtdb.gtdb_host('auto'),
                         'data.gtdb.ecogenomic.org')
        self.assertEqual(gtdb.gtdb_host(None), 'data.gtdb.ecogenomic.org')
        self.assertEqual(gtdb.gtdb_host('ecogenomic'),
                         'data.gtdb.ecogenomic.org')

    def test_aau_mirror_warns_instead_of_failing_later(self):
        import contextlib
        import io
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            host = gtdb.gtdb_host('aau')
        self.assertEqual(host, 'data.gtdb.aau.ecogenomic.org')
        self.assertIn('expired', buf.getvalue())
        with self.assertRaises(ValueError):
            gtdb.gtdb_host('mirror-of-convenience')

    # --- RK8: a missing DIAMOND must be an actionable error ------------
    def test_demo_tool_refuses_cleanly_without_diamond(self):
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        from unittest import mock
        path = join(HERE, os.pardir, os.pardir, 'tools', 'run_hgt_demo.py')
        spec = spec_from_file_location('rhd', path,
                                       loader=SourceFileLoader('rhd', path))
        rhd = module_from_spec(spec)
        spec.loader.exec_module(rhd)
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        argv = ['--db', join(d, 'db'), '--accession', 'RS_GCF_000000000.1',
                '--out-dir', join(d, 'out'), '--gtdb-dir', d,
                '--diamond', 'definitely-not-installed-here']
        with mock.patch.object(rhd.shutil, 'which', lambda x: None):
            with self.assertRaises(SystemExit) as ctx:
                rhd.main(argv)
        msg = str(ctx.exception)
        for hint in ('conda install', '--diamond', '--skip-search'):
            self.assertIn(hint, msg)
        # the check runs before the dataset is opened, so a bogus --gtdb-dir
        # above never gets far enough to be read


class TestRemoteChannelOptions(unittest.TestCase):
    """§7.1: `--gtdb-mirror`, `--downloader` and `--verify-md5` must be honoured.

    The remote channels may not download in this environment, but the options
    used to be parsed and then dropped on the floor: the plan they shape has to
    carry them, and the CLI has to win over config.yml.
    """

    def setUp(self):
        from types import SimpleNamespace
        from neoHGT import database
        self.db = database
        self.NS = SimpleNamespace
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        os.makedirs(join(self.d, 'gtdb'), exist_ok=True)
        self.fetch = [gtdb.GenomeFetch('RS_GCF_000000001.1', 'G1',
                                       '100000001', path=None, url=None,
                                       kind='local')]

    def _stub(self, **kw):
        vals = dict(output=self.d, seq_source='gtdb-tar',
                    gtdb_release='r220', gtdb_mirror='ecogenomic',
                    downloader='curl', verify_md5=True, allow_network=False)
        vals.update(kw)
        return self.NS(**vals)

    def _plan(self, **kw):
        me = self._stub(**kw)
        with self.assertRaises(ValueError) as ctx:
            self.db.plan_remote_fetches(me, list(self.fetch))
        self.assertIn('--allow-network', str(ctx.exception))
        rows = [l.rstrip('\n').split('\t')
                for l in open(join(self.d, 'gtdb', 'fetch_plan.tsv'))
                if not l.startswith('#')]
        return rows[0], dict(zip(rows[0][1:], rows[1][1:]))

    def test_plan_carries_the_downloader_and_md5_policy(self):
        head, row = self._plan()
        self.assertEqual(head, ['accession', 'url', 'kind', 'downloader',
                                'md5_check', 'command'])
        self.assertEqual(row['downloader'], 'curl')
        self.assertEqual(row['md5_check'], 'yes')
        self.assertIn('curl -L', row['command'])
        self.assertIn('tar -xzO', row['command'],
                      'the gtdb-tar channel streams the release bundle')
        self.assertIn('data.gtdb.ecogenomic.org', row['url'])

    def test_wget_and_no_md5_are_reflected_in_the_plan(self):
        _head, row = self._plan(downloader='wget', verify_md5=False)
        self.assertEqual(row['downloader'], 'wget')
        self.assertEqual(row['md5_check'], 'no')
        self.assertTrue(row['command'].startswith('wget -q -O - '),
                        row['command'])

    def test_aau_mirror_is_used_and_warned_about(self):
        import contextlib
        import io
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            _head, row = self._plan(gtdb_mirror='aau')
        self.assertIn('data.gtdb.aau.ecogenomic.org', row['url'])
        self.assertIn('expired', buf.getvalue())

    def test_command_line_beats_config_for_the_remote_knobs(self):
        class Loose(self.NS):
            """Test double: unset option names read as None, like argparse would
            leave them when the option is not in that subcommand's group."""

            def __getattr__(self, name):
                return None

        args = Loose(
            source='gtdb', gtdb_dir=FX, seq_source=None, local_layout=None,
            gtdb_release=None, gtdb_taxid_base=None, gtdb_mirror='aau',
            downloader='wget', verify_md5='no', allow_network='no',
            keep_faa=None, max_genomes=None, limit=None, force=False,
            input_fingerprints=None, gtdb_no_selfcheck=None,
            workspace_out=None, workspace_tmp=None, gtdb_taxa=None)
        me = Loose(cfg={'database': {'gtdb': {'mirror': 'ecogenomic',
                                               'taxid_base': 200000000}}},
                   output=self.d, source='gtdb')
        self.db.set_source_parameters(me, args)
        self.assertEqual(me.gtdb_mirror, 'aau')     # CLI wins over config.yml
        self.assertEqual(me.downloader, 'wget')
        self.assertFalse(me.verify_md5)
        self.assertEqual(int(me.id_base), 200000000)  # config fills what is unset

    def test_ncbi_and_ena_channels_keep_their_urls(self):
        _head, row = self._plan(seq_source='ncbi')
        self.assertIn('ftp.ncbi.nlm.nih.gov', row['url'])
        self.assertIn('curl -L -C -', row['command'],
                      'per-file channels keep resumable downloads; only '
                      'gtdb-tar streams a bundle')
        _head, row = self._plan(seq_source='ena')
        self.assertIn('www.ebi.ac.uk', row['url'])


class TestHybridCliWiring(unittest.TestCase):
    """§5.7 / §7.1: `--source hybrid` must actually reach HybridProvider.

    Before this wiring the two hybrid options were parsed and dropped, and a
    hybrid build silently produced a plain GTDB database. The merge logic
    itself is covered by TestHybrid; this drives the selection code.
    """

    def setUp(self):
        from types import SimpleNamespace
        from neoHGT import database
        self.db = database
        self.NS = SimpleNamespace

        class Loose(SimpleNamespace):
            """Test double whose unset attributes read as None."""

            def __getattr__(self, name):
                return None
        self.Loose = Loose

    def _self(self, source='hybrid', **kw):
        vals = dict(source=source, gtdb_dir='/nowhere', gtdb_release='220',
                    local_layout='auto', id_base=100000000,
                    no_selfcheck=True, gtdb_metadata=None, gtdb_taxonomy=None,
                    output='/tmp-out', hybrid_ncbi_taxdump='/ncbi/dump',
                    hybrid_lca_policy='root')
        vals.update(kw)
        return self.NS(**vals)

    def _make(self, me):
        from unittest import mock

        class RealGtdb:
            def __init__(self, *a, **k):
                pass

            def prepare(self):
                return {}

        captured = {}

        class SpyHybrid:
            def __init__(self, inner, ncbi_dir, lca_policy='root'):
                captured.update(inner=inner, ncbi_dir=ncbi_dir,
                                lca_policy=lca_policy)

            def prepare(self):
                captured['prepared'] = True
                return {}
        with mock.patch.object(self.db._gtdb, 'GtdbProvider', RealGtdb), \
                mock.patch.object(self.db._gtdb, 'HybridProvider', SpyHybrid):
            prov = self.db.make_provider(me, self.NS(gtdb_id_manifest=None))
        return prov, captured

    def test_hybrid_source_wraps_the_gtdb_provider(self):
        prov, seen = self._make(self._self())
        self.assertEqual(seen['ncbi_dir'], '/ncbi/dump')
        self.assertEqual(seen['lca_policy'], 'root')
        self.assertTrue(seen['prepared'])
        self.assertIsNotNone(prov)

    def test_gtdb_source_does_not_wrap(self):
        _prov, seen = self._make(self._self(source='gtdb'))
        self.assertEqual(seen, {})

    def test_selection_options_reach_the_wrapped_provider(self):
        # the build assigns provider.opts on the wrapper; a hybrid run that
        # forgets to forward them keeps every genome (silent whole-dataset build)
        seen = {}

        class Inner:
            def select(self, df, opts=None):
                seen['select'] = opts
                return df

            def sample(self, df, opts=None):
                seen['sample'] = opts
                return df
        h = gtdb.HybridProvider.__new__(gtdb.HybridProvider)
        h.gtdb = Inner()
        opts = {'genome_list': ['GB_GCA_000008085.1'], 'reps': True}
        h.opts = opts
        h.select('DF')
        h.sample('DF')
        self.assertEqual(seen, {'select': opts, 'sample': opts})
        # an explicit argument still wins over the stored options
        h.select('DF', {'reps': False})
        self.assertEqual(seen['select'], {'reps': False})

    def test_hybrid_requires_an_ncbi_taxdump(self):
        args = self.NS(cfg={}, source='hybrid', gtdb_dir=FX,
                       seq_source='local', local_layout=None,
                       gtdb_release=None, gtdb_taxid_base=None,
                       gtdb_mirror=None, downloader=None, verify_md5=None,
                       allow_network=None, keep_faa=None, max_genomes=None,
                       limit=None, force=False, input_fingerprints=None,
                       gtdb_no_selfcheck=None, workspace_out=None,
                       workspace_tmp=None, gtdb_taxa=None,
                       hybrid_ncbi_taxdump=None, hybrid_lca_policy=None)
        me = self.Loose(cfg={}, output='/x')
        with self.assertRaises(ValueError) as ctx:
            self.db.set_source_parameters(me, args)
        self.assertIn('--hybrid-ncbi-taxdump', str(ctx.exception))
        self.assertIn('--source gtdb', str(ctx.exception))

    def test_unknown_lca_policy_is_rejected(self):
        args = self.NS(cfg={}, source='hybrid', gtdb_dir=FX,
                       seq_source='local', local_layout=None,
                       gtdb_release=None, gtdb_taxid_base=None,
                       gtdb_mirror=None, downloader=None, verify_md5=None,
                       allow_network=None, keep_faa=None, max_genomes=None,
                       limit=None, force=False, input_fingerprints=None,
                       gtdb_no_selfcheck=None, workspace_out=None,
                       workspace_tmp=None, gtdb_taxa=None,
                       hybrid_ncbi_taxdump=FX, hybrid_lca_policy='greedy')
        me = self.Loose(cfg={}, output='/x')
        with self.assertRaises(ValueError) as ctx:
            self.db.set_source_parameters(me, args)
        self.assertIn('root', str(ctx.exception))


class TestExplicitTablePaths(unittest.TestCase):
    """§7.1: --gtdb-metadata / --gtdb-taxonomy must be honoured and checked."""

    def test_noop_without_overrides(self):
        info = {'metadata': ['a'], 'taxonomy': ['b']}
        self.assertIs(gtdb.apply_table_overrides(info), info)

    def test_missing_file_is_reported(self):
        with self.assertRaises(ValueError) as ctx:
            gtdb.apply_table_overrides({'metadata': [], 'taxonomy': []},
                                       metadata=['/nope/bac120_metadata_r220.tsv'])
        self.assertIn('--gtdb-metadata', str(ctx.exception))

    def test_cross_release_combination_is_refused(self):
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        md = join(d, 'bac120_metadata_r220.tsv')
        tx = join(d, 'bac120_taxonomy_r232.tsv')
        for p in (md, tx):
            open(p, 'w').write('accession\n')
        with self.assertRaises(ValueError) as ctx:
            gtdb.apply_table_overrides({'metadata': [], 'taxonomy': []},
                                       metadata=[md], taxonomy=[tx])
        self.assertIn('different GTDB releases', str(ctx.exception))

    def test_family_mismatch_is_refused(self):
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        md = join(d, 'bac120_metadata_r220.tsv')
        tx = join(d, 'ar53_taxonomy_r220.tsv')
        for p in (md, tx):
            open(p, 'w').write('accession\n')
        with self.assertRaises(ValueError) as ctx:
            gtdb.apply_table_overrides({'metadata': [], 'taxonomy': []},
                                       metadata=[md], taxonomy=[tx])
        self.assertIn('same families', str(ctx.exception))

    def test_matching_pair_is_substituted(self):
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        md = join(d, 'bac120_metadata_r220.tsv')
        tx = join(d, 'bac120_taxonomy_r220.tsv')
        for p in (md, tx):
            open(p, 'w').write('accession\n')
        got = gtdb.apply_table_overrides(
            {'metadata': ['discovered'], 'taxonomy': ['discovered'],
             'release': '220'}, metadata=[md], taxonomy=[tx])
        self.assertEqual([basename(p) for p in got['metadata']],
                         [basename(md)])
        self.assertEqual([basename(p) for p in got['taxonomy']],
                         [basename(tx)])


class TestPlanVerificationIds(unittest.TestCase):
    """UT-2.1 / UT-0.3 / UT-3.3 / UT-3.4: the §10.3 ids that had no test.

    A traceability sweep of the plan's verification tables found these four
    named requirements without a test behind them, so each is pinned here.
    """

    STEM = ('RS_GCF_000026225.1_d_Bacteria_p_Pseudomonadota_c_'
            'Gammaproteobacteria_o_Enterobacterales_f_Enterobacteriaceae'
            '_g_Escherichia')

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.files = [join(FX, 'metadata_bac120_sample.tsv'),
                      join(FX, 'metadata_ar53_sample.tsv')]

    # --- UT-2.1 / UT-2.x -------------------------------------------------
    def test_metadata_is_parsed_by_column_name_not_position(self):
        df = gtdb.load_metadata(self.files, release='220')
        # reorder the columns of the first file: the result must be identical,
        # which is only true if the mapping is by name (plan §10 UT-2.1)
        shuffled = join(self.d, 'shuffled.tsv')
        with gtdb.open_reader(self.files[0]) as f:
            lines = f.read().splitlines()
        head = lines[0].split('\t')
        order = sorted(range(len(head)), key=lambda i: (i % 7, i))
        with open(shuffled, 'w') as f:
            f.write('\t'.join(head[i] for i in order) + '\n')
            for line in lines[1:]:
                cells = line.split('\t')
                f.write('\t'.join(cells[i] for i in order) + '\n')
        got = gtdb.load_metadata([shuffled] + self.files[1:], release='220')
        self.assertEqual(sorted(df.columns), sorted(got.columns))
        self.assertEqual(df['accession'].tolist(), got['accession'].tolist())
        self.assertEqual(df['gtdb_taxonomy'].tolist(),
                         got['gtdb_taxonomy'].tolist())

    def test_rows_are_preserved_and_keys_unique(self):
        df = gtdb.load_metadata(self.files, release='220')
        want = sum(sum(1 for _ in gtdb.open_reader(p)) - 1 for p in self.files)
        self.assertEqual(len(df), want)          # no silent filtering at load
        self.assertTrue(df['accession'].is_unique)
        self.assertTrue(df['genome'].is_unique)
        # the plan's "protein_count sum == planned proteins" invariant: the
        # build's capacity plan is exactly this column sum
        self.assertGreater(int(df.loc[df['is_rep'], 'proteins'].astype(int)
                                 .sum()), 0)

    def test_missing_column_names_the_column(self):
        src = self.files[0]
        bad = join(self.d, 'broken.tsv')
        with gtdb.open_reader(src) as f:
            lines = f.read().splitlines()
        head = lines[0].split('\t')
        drop = head.index('gtdb_taxonomy')
        with open(bad, 'w') as f:
            for line in lines:
                cells = line.split('\t')
                f.write('\t'.join(c for i, c in enumerate(cells)
                                  if i != drop) + '\n')
        with self.assertRaises(ValueError) as ctx:
            gtdb.load_metadata([bad], release='220')
        msg = str(ctx.exception)
        self.assertIn('gtdb_taxonomy', msg)
        self.assertIn('by name', msg)

    # --- UT-0.3 ----------------------------------------------------------
    def test_selfcheck_fails_loudly_on_one_inconsistent_file_count(self):
        df = gtdb.load_metadata(self.files, release='220')
        n_bac = int(((df['domain'] == 'bacteria') & df['is_rep']).sum())
        info = {'gtdb_dir': self.d, 'release': '220',
                'layout': 'gtdb-r220-renamed', 'metadata': self.files,
                'taxonomy': [join(FX, 'taxonomy_sample.tsv')],
                'idlists': {}, 'provenance_scripts': [],
                'absent_official': [], 'serial_anchor_coverage': 1.0,
                'max_first_line_len': 270,
                'domains': {'bacteria': {'dir': join(FX, 'protein'),
                                         'n_files': n_bac + 1},
                            'archaea': {'dir': join(FX, 'protein'),
                                        'n_files': int(((df['domain'] ==
                                                         'archaea') &
                                                        df['is_rep']).sum())}}}
        res = gtdb.selfcheck(info, strict=False, deep=False)
        failed = [c for c in res['checks'] if c['status'] == 'FAIL']
        self.assertEqual([c['name'] for c in failed],
                         ['faa_vs_metadata_reps[bacteria]'])
        self.assertEqual(failed[0]['expected'], n_bac)
        self.assertEqual(failed[0]['observed'], n_bac + 1)
        with self.assertRaises(gtdb.SelfCheckError) as ctx:
            gtdb.selfcheck(info, strict=True, deep=False)
        self.assertIn('faa_vs_metadata_reps[bacteria]', str(ctx.exception))
        self.assertIn('Refusing to build', str(ctx.exception))

    # --- UT-3.3 ----------------------------------------------------------
    def _search_stub(self, keep=False):
        s = Search.__new__(Search)
        s.source, s.method = 'gtdb', 'precomp'
        s.gtdb_root = join(self.d, 'db')
        s.gtdb_layout = 'gtdb-r220-renamed'
        # None (not {}) is what makes _load_route_c() read the tables
        s.gtdb_stems, s.gtdb_p2g, s.gtdb_g2t = None, None, None
        s.gtdb_stem2taxid = None
        s.prot2tid, s.taxmap = {}, None
        s.keep_unresolved, s.n_unresolved = keep, 0
        s.tmpdir = s.output = self.d
        os.makedirs(join(s.gtdb_root, 'gtdb'), exist_ok=True)
        with gzip.open(join(s.gtdb_root, 'gtdb', 'stem2taxid.tsv.gz'),
                       'wt') as f:
            f.write(f'{self.STEM}\t1006000001\n')
        return s

    @staticmethod
    def _hit(subject, taxid=''):
        return {'id': subject, 'taxid': taxid, 'identity': 90.0,
                'evalue': 1e-40, 'score': 300.0, 'coverage': 95.0}

    def test_route_a_then_c_then_unresolved_keeps_hits(self):
        s = self._search_stub(keep=True)
        a = self._hit('GB_GCA_000008085.1_d_Archaea_x', '1000000000')
        c = self._hit(f'{self.STEM}{gtdb.SERIAL_ANCHOR}NC_000001.1_7')
        u = self._hit('totally_unparseable_subject')
        prots = {'q1': [a], 'q2': [c], 'q3': [u]}
        taxid_wf_gtdb(s, prots)
        self.assertEqual(a['taxid'], '1000000000')        # route A untouched
        self.assertEqual(c['taxid'], '1006000001')       # route C resolved
        self.assertEqual(u['taxid'], '')                 # RK5: not guessed
        self.assertEqual(s.n_unresolved, 1)
        unresolved = open(join(self.d, 'unresolved_taxids.txt')).read()
        self.assertIn('totally_unparseable_subject', unresolved)
        hits = open(join(self.d, 'unresolved_hits.tsv')).read()
        self.assertIn('totally_unparseable_subject', hits)
        self.assertEqual(len(prots['q3']), 1,
                         '--keep-unresolved yes must retain the hit record')

    def test_unresolved_hit_dropped_when_not_kept(self):
        s = self._search_stub(keep=False)
        prots = {'q3': [self._hit('totally_unparseable_subject')]}
        taxid_wf_gtdb(s, prots)
        self.assertEqual(prots['q3'], [])
        self.assertEqual(s.n_unresolved, 1)

    # --- UT-3.4 ----------------------------------------------------------
    def test_gtdb_mode_never_touches_update_dmp_files(self):
        s = self._search_stub()

        def boom(*a, **k):
            raise AssertionError(
                'GTDB mode must not extend the taxdump at search time')
        s.update_dmp_files = boom
        prots = {'q1': [self._hit('GB_GCA_000008085.1_d_Archaea_x',
                                 '1000000000')]}
        taxid_wf_gtdb(s, prots)                       # no exception => not called
        self.assertEqual(prots['q1'][0]['taxid'], '1000000000')


    def test_protein_total_survives_column_normalization(self):
        """UT-2.x: the protein total must be a checked quantity, not a number
        quoted from the paper. Reading the raw TSV independently has to give the
        same total as the normalized frame, with no row silently dropped."""
        files = [join(FX, 'metadata_bac120_sample.tsv'),
                 join(FX, 'metadata_ar53_sample.tsv')]
        declared, rows = 0, 0
        for path in files:
            with open(path) as f:
                head = f.readline().rstrip('\n').split('\t')
                ip = head.index('protein_count')
                for line in f:
                    x = line.rstrip('\n').split('\t')
                    if not x or not x[0]:
                        continue
                    rows += 1
                    declared += int(x[ip])
        df = gtdb.load_metadata(files, release='220')
        self.assertEqual(int(df['proteins'].sum()), declared)
        self.assertEqual(df.shape[0], rows,
                         'no genome row may vanish during column mapping')
        self.assertFalse(df['proteins'].isna().any(),
                         'a protein count that will not parse must not become '
                         'a silent NaN in the total')
        self.assertEqual(len(set(df['accession'])), df.shape[0])

    def test_proteins_reconciliation_labels_both_sides_correctly(self):
        """`proteins_planned` is the promised count, never the extracted one.

        The first implementation unpacked the helper positionally and bound the
        planned field to the extracted count, so a real build reported itself as
        both complete and "1,944,720,987 proteins short".
        """
        selection = pd.DataFrame({'proteins': [10, 20]})
        release = pd.DataFrame({'proteins': [10, 20, 1957679562]})
        ok = gtdb.proteins_reconciliation(selection, {'proteins': 30})
        self.assertEqual(ok['proteins_planned'], 30)
        self.assertEqual(ok['proteins_extracted_vs_planned'], 'ok')
        short = gtdb.proteins_reconciliation(selection, {'proteins': 25})
        self.assertEqual(short['proteins_planned'], 30,
                         'planned must not track the extracted count')
        self.assertEqual(short['proteins_extracted_vs_planned'], '5 short')
        wrong = gtdb.proteins_reconciliation(release, {'proteins': 30})
        self.assertEqual(wrong['proteins_planned'], 1957679592)
        self.assertIn('short', wrong['proteins_extracted_vs_planned'],
                      'the release table is the wrong basis and must scream')

    def test_protein_counts_are_reconciled_against_extraction(self):
        """The manifest must be able to say that a build extracted fewer
        proteins than the metadata promised."""
        df = gtdb.load_metadata([join(FX, 'metadata_bac120_sample.tsv'),
                                 join(FX, 'metadata_ar53_sample.tsv')],
                                release='220')
        total = int(df['proteins'].sum())
        self.assertEqual(gtdb.reconcile_protein_counts(df, {'proteins': total}),
                         (total, total, 'ok'))
        planned, got, verdict = gtdb.reconcile_protein_counts(
            df, {'proteins': total - 1000})
        self.assertEqual((planned, got), (total, total - 1000))
        self.assertIn('short', verdict)
        self.assertEqual(gtdb.reconcile_protein_counts(df, {})[2],
                         'not extracted yet')

    def test_meta_publishes_the_reconciliation(self):
        """The manifest must publish both sides of the reconciliation and the
        caller must hand it the *selection*, not the loaded release table."""
        body = inspect.getsource(gtdb.write_meta_json)
        self.assertIn('recon = proteins_reconciliation(sel, stats)', body)
        self.assertIn("'proteins_planned': recon['proteins_planned']", body)
        self.assertIn("'proteins_extracted_vs_planned':", body)
        self.assertNotIn(
            'reconcile_protein_counts', body,
            'unpacking the helper positionally is what bound the extracted '
            'count to the planned field')

        dbsrc = open(join(HERE, '..', 'database.py'), encoding='utf-8').read()
        self.assertIn('def write_gtdb_meta(self, provider, args, ext, '
                      'started, df=None)', dbsrc)
        self.assertIn('selection=df)', dbsrc,
                      'the build must reconcile against the genomes it was '
                      'asked to extract')
        self.assertIn('write_gtdb_meta(self, provider, args, ext, started, df)',
                      dbsrc)


class TestReconcile(unittest.TestCase):
    """DoD-2: the accession-level mismatch lists must be written out."""

    LIN = 'd__Bacteria;p__X;c__Y;o__Z;f__F;g__G;s__G g1'

    def setUp(self):
        import pandas as pd
        from unittest import mock
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.info = {
            'gtdb_dir': join(self.d, 'src'), 'release': '999',
            'layout': 'official-plain',
            'metadata': [], 'taxonomy': [], 'idlists': {},
            'domains': {'bacteria': {'dir': join(self.d, 'src'),
                                    'n_files': 0}},
            'provenance_scripts': [], 'absent_official': [],
            'serial_anchor_coverage': 1.0, 'max_first_line_len': 40}
        self.pd = pd
        self.mock = mock

    def _patch(self, md_rows, tax_rows, files):
        import pandas as pd

        def load_metadata(files, release=None):
            return pd.DataFrame(md_rows, columns=[
                'accession', 'gtdb_taxonomy', 'gtdb_representative',
                'is_rep', 'domain'])

        def read_taxonomy(files):
            return dict(tax_rows)

        def list_faa(dir_):
            return files
        for name, fn in (('load_metadata', load_metadata),
                         ('read_taxonomy', read_taxonomy),
                         ('_list_faa', list_faa)):
            patcher = self.mock.patch.object(gtdb, name, fn)
            patcher.start()
            self.addCleanup(patcher.stop)

    def _faa(self, acc):
        return f'{acc}_d_Bacteria_p_X_c_Y_o_Z_f_F_g_G.faa'

    def test_lists_capture_every_kind_of_difference(self):
        self._patch(
            [('RS_A.1', self.LIN, 't', True, 'bacteria'),
             ('RS_B.1', self.LIN, 't', True, 'bacteria'),
             ('GB_C.1', self.LIN, '.', False, 'bacteria')],
            {'RS_A.1': self.LIN, 'RS_B.1': self.LIN, 'GB_C.1': self.LIN,
             'RS_TAXONLY.1': self.LIN},
            [self._faa('RS_A.1'), self._faa('RS_B.1'), self._faa('GB_C.1'),
             self._faa('RS_GHOST.1')])
        out = join(self.d, 'rec')
        res = gtdb.reconcile(info=self.info, out=out)
        self.assertFalse(res['consistent'])
        self.assertEqual(res['lists'], {
            'missing_faa.tsv': 0, 'extra_faa.tsv': 1, 'nonrep_faa.tsv': 1,
            'metadata_not_in_taxonomy.tsv': 0,
            'taxonomy_not_in_metadata.tsv': 1})
        body = [l.strip() for l in open(join(out, 'extra_faa.tsv'))
                if not l.startswith('#')]
        self.assertEqual(body, ['accession', 'RS_GHOST.1'])
        missing = [l.strip() for l in open(join(out, 'missing_faa.tsv'))
                   if not l.startswith('#')]
        self.assertEqual(missing, ['accession'])     # empty but present

    def test_consistent_copy_is_flagged(self):
        self._patch([('RS_A.1', self.LIN, 't', True, 'bacteria')],
                    {'RS_A.1': self.LIN}, [self._faa('RS_A.1')])
        res = gtdb.reconcile(info=self.info, out=join(self.d, 'rec2'))
        self.assertTrue(res['consistent'])
        self.assertTrue(all(v == 0 for v in res['lists'].values()))

    def test_missing_local_file_is_listed(self):
        self._patch([('RS_A.1', self.LIN, 't', True, 'bacteria'),
                     ('RS_GONE.1', self.LIN, 't', True, 'bacteria')],
                    {'RS_A.1': self.LIN, 'RS_GONE.1': self.LIN},
                    [self._faa('RS_A.1')])
        res = gtdb.reconcile(info=self.info, out=join(self.d, 'rec3'))
        self.assertEqual(res['lists']['missing_faa.tsv'], 1)
        self.assertFalse(res['consistent'])

    def test_needs_dir_or_info(self):
        with self.assertRaises(ValueError):
            gtdb.reconcile(out=join(self.d, 'nope'))

    def test_cli_writes_json_report_and_exits(self):
        import json
        from unittest import mock

        self._patch([('RS_A.1', self.LIN, 't', True, 'bacteria')],
                    {'RS_A.1': self.LIN}, [self._faa('RS_A.1')])
        out, rep = join(self.d, 'rec_cli'), join(self.d, 'report.json')
        for name, fn in (('detect_layout', lambda *a, **k: dict(self.info)),
                         ('load_metadata', lambda *a, **k: self.pd.DataFrame(
                             [('RS_A.1', self.LIN, 't', True, 'bacteria')],
                             columns=['accession', 'gtdb_taxonomy',
                                      'gtdb_representative', 'is_rep',
                                      'domain'])),
                         ('read_taxonomy', lambda *a, **k:
                             {'RS_A.1': self.LIN}),
                         ('_list_faa', lambda d_: [self._faa('RS_A.1')])):
            patcher = self.mock.patch.object(gtdb, name, fn)
            patcher.start()
            self.addCleanup(patcher.stop)
        code = gtdb.main(['reconcile', '--dir', join(self.d, 'src'),
                         '--out', out, '--report', rep])
        self.assertEqual(code, 0)
        data = json.load(open(rep))
        self.assertTrue(data['consistent'])
        self.assertEqual(data['lists']['missing_faa.tsv'], 0)
        self.assertTrue(isfile(data['files']['extra_faa.tsv']))


class TestCliWrapper(unittest.TestCase):
    """'neoHGT gtdb <subcommand>' must reach the real GTDB interface."""

    @staticmethod
    def _load_wrapper():
        from importlib.machinery import SourceFileLoader
        from importlib.util import spec_from_file_location, module_from_spec
        path = join(HERE, os.pardir, os.pardir, 'scripts', 'neoHGT')
        spec = spec_from_file_location('neoHGT_cli_wrapper', path,
                                       loader=SourceFileLoader(
                                           'neoHGT_cli_wrapper', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_all_subcommands_have_a_parser(self):
        for name in gtdb.SUBCOMMANDS:
            with self.assertRaises(SystemExit) as ctx:   # --help exits 0
                gtdb.main([name, '--help'])
            self.assertEqual(ctx.exception.code, 0)

    def test_wrapper_delegates_and_propagates_exit_code(self):
        from unittest import mock
        cli = self._load_wrapper()
        seen = {}

        def fake_main(argv):
            seen['argv'] = argv
            return 7
        with mock.patch.object(gtdb, 'main', fake_main), \
                mock.patch.object(sys, 'argv',
                                  ['neoHGT', 'gtdb', 'reconcile',
                                   '--dir', 'D', '--out', 'O']):
            with self.assertRaises(SystemExit) as ctx:
                cli.main()
        self.assertEqual(ctx.exception.code, 7)
        self.assertEqual(seen['argv'],
                         ['reconcile', '--dir', 'D', '--out', 'O'])

    def test_wrapper_keeps_the_legacy_flat_form(self):
        from types import SimpleNamespace
        from unittest import mock
        cli = self._load_wrapper()
        calls = []

        def fake_main(argv):
            calls.append(argv)
            return 0

        parsed = SimpleNamespace(command='gtdb', action='probe', dir='D',
                                 output='O', gtdb_release=None,
                                 local_layout='auto', gtdb_taxid_base=1,
                                 gtdb_id_manifest=None, no_selfcheck=True,
                                 deep=False, threads=1, tmpdir=None, to='ncbi')
        with mock.patch.object(gtdb, 'main', fake_main), \
                mock.patch.object(gtdb, 'cmd_probe',
                                  lambda *a, **k: 'probed'), \
                mock.patch.object(cli, 'parse_args', lambda modules: parsed), \
                mock.patch.object(sys, 'argv',
                                  ['neoHGT', 'gtdb', '--action', 'probe',
                                   '--dir', 'D', '-o', 'O']):
            self.assertEqual(cli.main(), None)
        self.assertEqual(calls, [])      # no delegation for the flat form


class TestGtdbActionChoices(unittest.TestCase):
    """Every legal `--action` value must reach its own implementation.

    The value-coverage rule read only three of the four CLI surfaces, so
    `neoHGT gtdb`'s own table was never asked: `--action taxdump` and
    `--action map-report` were values that no test handed to the option which
    declares them, even though `--action probe` was covered. Each value is
    parsed by the shipped parser and must arrive at its function with its own
    arguments in the right order.
    """

    def _dispatch(self, argv, target):
        from unittest import mock
        args = cli_args('gtdb', argv)
        seen = []

        def recorder(*a, **k):
            seen.append((a, k))
            return 'dispatched'

        with mock.patch.object(gtdb, target, recorder):
            out = gtdb.Gtdb()(args)
        return args, seen, out

    def test_probe_action_arrives_at_the_prober(self):
        args, seen, out = self._dispatch(
            ['--action', 'probe', '--dir', 'D', '-o', 'O', '--deep',
             '--gtdb-release', '220'], 'cmd_probe')
        self.assertEqual(out, 'dispatched')
        self.assertEqual(seen[0][0][:3], ('D', 'O', '220'), seen)

    def test_taxdump_action_arrives_at_the_converter(self):
        args, seen, out = self._dispatch(
            ['--action', 'taxdump', '--dir', 'D', '-o', 'O',
             '--gtdb-taxid-base', '123456', '--no-selfcheck'], 'cmd_taxdump')
        self.assertEqual(out, 'dispatched')
        self.assertEqual(seen[0][0][:2], ('D', 'O'), seen)
        self.assertEqual(seen[0][0][4], 123456,
                         'the id base has to reach the converter, not the '
                         'manifest argument')
        self.assertEqual(seen[0][0][6], False,
                         '--no-selfcheck must arrive as "do not check"')

    def test_map_report_action_arrives_at_the_mapper(self):
        args, seen, out = self._dispatch(
            ['--action', 'map-report', '--db', 'DB', '--in', 'I', '-o', 'O',
             '--to', 'ncbi'], 'map_report')
        self.assertEqual(out, 'dispatched')
        self.assertEqual(seen[0][0], ('DB', 'I', 'O', 'ncbi'), seen)

    def test_an_unknown_action_is_refused_by_the_parser(self):
        with self.assertRaises(SystemExit):
            cli_args('gtdb', ['--action', 'report-maps', '--dir', 'D'])

    def test_the_three_actions_are_three_different_destinations(self):
        # a single dispatch would show up here: all three values landing in the
        # same function is what an inert option looks like from the outside
        from unittest import mock
        called = []
        for name, argv in (('cmd_probe', ['--action', 'probe']),
                           ('cmd_taxdump', ['--action', 'taxdump']),
                           ('map_report', ['--action', 'map-report'])):
            args = cli_args('gtdb', argv + ['--dir', 'D', '--db', 'D',
                                            '--in', 'I', '-o', 'O'])
            with mock.patch.object(gtdb, 'cmd_probe',
                                   lambda *a, **k: called.append('probe')), \
                    mock.patch.object(gtdb, 'cmd_taxdump',
                                      lambda *a, **k: called.append('taxdump')), \
                    mock.patch.object(gtdb, 'map_report',
                                      lambda *a, **k: called.append('map')):
                gtdb.Gtdb()(args)
        self.assertEqual(called, ['probe', 'taxdump', 'map'], called)


    def test_the_value_rule_reads_all_four_cli_surfaces(self):
        # the rule that asks "did some test hand this value to its option?" used
        # to loop over three stages, so `neoHGT gtdb`'s own table was never asked
        v = TestCliSurfaceAndDocExamples._tool()
        code = join(HERE, os.pardir, os.pardir)
        entries = v.cli_choice_options(code)
        self.assertEqual({e['stage'] for e in entries} - set(v.CLI_STAGES), set(),
                         'a CLI surface is in the table but not in CLI_STAGES')
        gtdb_entry = [e for e in entries if e['stage'] == 'gtdb']
        self.assertEqual([e['flags'] for e in gtdb_entry], [['--action']],
                         'the gtdb table has no option with choices?')
        gated, backlog = v.choice_values_without_tests(
            code, v.declared_choice_values(
                join(code, 'docs', '交付物清单.md'), code))
        self.assertEqual(backlog, [], 'a stage was left unjudged')
        self.assertEqual(gated, [], gated)

    def test_the_value_rule_bites_on_a_surface_it_has_just_gained(self):
        # a fresh CLI surface with two values and no test handing either must
        # produce two problems, not silence
        v = TestCliSurfaceAndDocExamples._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'neoHGT', 'tests'))
        with open(join(d, 'neoHGT', 'gtdb.py'), 'w', encoding='utf-8') as fh:
            fh.write("arguments = [['--mystery', 'x',\n"
                     "  {'choices': ['alpha', 'beta']}]]\n")
        for stage in ('database', 'search', 'analyze'):
            with open(join(d, 'neoHGT', f'{stage}.py'), 'w') as fh:
                fh.write('arguments = []\n')
        with open(join(d, 'neoHGT', 'tests', 'test_anything.py'), 'w') as fh:
            fh.write('import unittest\n\n\nclass T(unittest.TestCase):\n'
                     '    def test_one(self):\n        assert True\n')
        gated, backlog = v.choice_values_without_tests(d, ())
        self.assertEqual(gated, ['--mystery=alpha', '--mystery=beta'],
                         (gated, backlog))
        self.assertEqual(backlog, [])


class TestFixtureProvenance(unittest.TestCase):
    """Every committed fixture sequence has to trace back to the read-only release.

    `fixtures.json` recorded hashes once, when the generator wrote it, and CI has
    claimed "KB-scale excerpts of GTDB r220 genomes" ever since -- but nothing
    re-read the data, so an edited, invented or unlisted fixture would have stayed
    invisible. The controls build their own tiny trees, so they run in CI where
    the dataset is not mounted.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_fixture_provenance.py')
        spec = spec_from_file_location('vfp', path,
                                       loader=SourceFileLoader('vfp', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _tree(self, v, tamper=None):
        """A manifest, a two-record source .faa and a one-record excerpt."""
        d = tmpdir()
        self.root = d
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'release', 'bacteria'))
        os.makedirs(join(d, 'fx', 'protein'))
        src = join(d, 'release', 'bacteria', 'GB_GCA_1.faa')
        with open(src, 'w') as fh:
            fh.write('>p1_1\nMKV\n>p1_2\nLLL\n')
        kb = join(d, 'fx', 'protein', 'GB_GCA_1.faa')
        with open(kb, 'w') as fh:
            fh.write('>p1_1\nMKV\n')
        with open(join(d, 'fx', 'db.faa'), 'w') as fh:
            fh.write('>p1_1\nMKV\n')

        def write(man):
            with open(join(d, 'fx', 'fixtures.json'), 'w',
                      encoding='utf-8') as fh:
                fh.write(json.dumps(man))

        man = {'release': '220', 'gtdb_dir': join(d, 'release'),
               'genomes': [{'accession': 'GB_GCA_1', 'label': 'x',
                            'domain': 'bacteria', 'source_path': src,
                            'source_sha256': v.sha256(src),
                            'source_bytes': os.path.getsize(src),
                            'kb_excerpt': kb, 'fixture_copy': kb,
                            'fixture_sha256': v.sha256(kb)}],
               'derived': {'db_faa_records': 1}}
        write(man)
        if tamper:
            tamper(d, man, src, kb, write)
        return v.check(join(d, 'fx', 'fixtures.json'))

    def test_a_clean_tree_reports_nothing(self):
        v = self._tool()
        _report, problems = self._tree(v)
        self.assertEqual(problems, [])

    def test_each_kind_of_drift_is_caught(self):
        v = self._tool()

        def edit_source(d, man, src, kb, write):
            with open(src, 'a') as fh:
                fh.write('>p1_3\nXXX\n')

        def invent(d, man, src, kb, write):
            with open(kb, 'w') as fh:
                fh.write('>p9_9\nZZZ\n')

        def alter(d, man, src, kb, write):
            with open(kb, 'w') as fh:
                fh.write('>p1_1\nMMM\n')

        def unlist(d, man, src, kb, write):
            with open(join(d, 'fx', 'protein', 'GB_ORPHAN.faa'), 'w') as fh:
                fh.write('>q1\nAAA\n')

        def stale_hash(d, man, src, kb, write):
            # the source stays as it is but the recorded hash is a lie
            man['genomes'][0]['source_sha256'] = 'dead' + 'beef'
            write(man)

        def gone(d, man, src, kb, write):
            os.remove(src)

        def count_drift(d, man, src, kb, write):
            with open(join(d, 'fx', 'db.faa'), 'a') as fh:
                fh.write('>p1_2\nLLL\n')

        for name, fn, needle in (('source changed under us', edit_source,
                                  'sha256'),
                                 ('invented sequence', invent, '不存在'),
                                 ('edited sequence', alter, '被改过'),
                                 ('unlisted fixture file', unlist,
                                  '没有对应条目'),
                                 ('stale recorded hash', stale_hash, '不符'),
                                 ('missing source', gone, '已不存在'),
                                 ('published count drifted', count_drift,
                                  'db.faa')):
            _report, problems = self._tree(v, fn)
            self.assertTrue(problems, f'{name} passed unnoticed')
            self.assertTrue(any(needle in p for p in problems),
                            f'{name}: {problems} (wanted {needle!r})')

    def test_the_ten_genome_limit_is_enforced(self):
        v = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'fx'))
        man = {'gtdb_dir': d, 'genomes': [{'accession': str(i)}
                                          for i in range(v.GENOME_LIMIT + 1)]}
        path = join(d, 'fx', 'fixtures.json')
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(json.dumps(man))
        _report, problems = v.check(path)
        self.assertTrue(any('10' in p for p in problems), problems)

    def test_a_release_mounted_elsewhere_still_verifies(self):
        # --gtdb-dir re-bases the recorded paths instead of distrusting them:
        # the same data under another root has to stay traceable
        v = self._tool()
        _report, problems = self._tree(v)
        self.assertEqual(problems, [])
        alt = join(self.root, 'release_alt')
        shutil.copytree(join(self.root, 'release', 'bacteria'),
                        join(alt, 'bacteria'))
        report2, problems2 = v.check(join(self.root, 'fx', 'fixtures.json'),
                                     alt)
        self.assertEqual(problems2, [], (report2, problems2))
        self.assertEqual(report2['sources_found'], 1)
        # ... and a root without the data says so, rather than passing silently
        _report3, problems3 = v.check(join(self.root, 'fx', 'fixtures.json'),
                                      join(self.root, 'release_empty'))
        self.assertTrue(any('已不存在' in p for p in problems3), problems3)

    def test_the_shipped_manifest_traces_back_to_the_dataset(self):
        v = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        path = join(code, 'neoHGT', 'tests', 'data', 'GTDB', 'fixtures.json')
        # pass the mounted release so the strict check runs where it can;
        # on a checkout without the dataset this stays None and the test
        # takes its no-dataset skipTest path below
        report, problems = v.check(path, gtdb_dir=os.environ.get('GTDB_DIR'))
        self.assertEqual(report['genomes'], 6, report)
        self.assertLessEqual(report['genomes'], v.GENOME_LIMIT)
        if report['sources_found'] == 0:
            # CI without the dataset: the only complaints allowed are the ones
            # that say so, never a provenance mismatch
            self.assertTrue(all('已不存在' in p for p in problems), problems)
            self.skipTest('the read-only GTDB release is not mounted here')
        self.assertEqual(report['sources_found'], report['genomes'], report)
        self.assertEqual(problems, [])
        self.assertEqual(report['db_faa_records'], 30,
                         'the committed fixture DB has 6 genomes x 5 proteins')


class TestProvenanceHeaders(unittest.TestCase):
    """Plan §5 principle 4: every derived file carries release, time, version.

    The two maps handed to external tools verbatim (taxon.map.gz, genome.map.gz)
    are the documented exception, so the test asserts their absence as well.
    """

    STAMPED = ('taxid.map.gz', 'labels.tsv', 'lineages.tsv',
               'genome2taxid.tsv.gz', 'ncbi_map.tsv.gz', 'stem2taxid.tsv.gz',
               'prefix2genome.tsv.gz')

    def setUp(self):
        import pandas as pd
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        df = gtdb.load_metadata([join(FX, 'metadata_bac120_sample.tsv')],
                                release='220')
        self.td, self.tmap, self.stats = gtdb.build_taxdump_from_gtdb(
            list(df['gtdb_taxonomy'].unique()))
        self.info = {'gtdb_dir': join(FX, 'protein'), 'release': '220',
                     'layout': 'gtdb-r220-renamed', 'provenance_scripts': [],
                     'absent_official': [], 'serial_anchor_coverage': 1.0,
                     'max_first_line_len': 270, 'domains': {}, 'taxonomy': [],
                     'metadata': [], 'idlists': {}, 'sample_sets': {}}
        stems = {f'{acc}_d_Bacteria_p_X_c_X_o_X_f_X_g_X': tid
                 for acc, tid in
                 [(df.iloc[0].accession, self.tmap[df.iloc[0].gtdb_taxonomy])]}
        gtdb.write_side_files(self.d, self.info, df, self.tmap, self.td,
                              self.stats, stem2taxid=stems,
                              prefix2genome={'NC_000001.1': df.iloc[0].accession})
        self.df = df

    def _first_line(self, name):
        path = join(self.d, 'gtdb', name)
        with gtdb.open_reader(path) as f:
            return f.readline()

    def test_every_derived_table_is_stamped(self):
        from neoHGT import __version__
        for name in self.STAMPED:
            head = self._first_line(name)
            self.assertTrue(head.startswith('# '), f'{name}: no header')
            self.assertIn(f'version={__version__}', head, name)
            self.assertIn('gtdb_release=220', head, name)
            self.assertIn('generated=', head, name)

    def test_loaders_still_round_trip_after_stamping(self):
        self.assertTrue(gtdb.load_labels(self.d))
        self.assertTrue(gtdb.load_manifest(join(self.d, 'gtdb',
                                                'taxid.map.gz')))
        self.assertTrue(gtdb.load_genome2taxid(self.d))
        self.assertTrue(gtdb.load_stem2taxid(self.d))
        self.assertTrue(gtdb.load_prefix2genome(self.d))
        self.assertTrue(gtdb.load_ncbi_map(self.d))

    def test_fingerprint_table_is_stamped(self):
        fp = join(self.d, 'gtdb', 'input_fingerprints.tsv')
        gtdb.fingerprint([join(FX, 'metadata_bac120_sample.tsv')], fp,
                         info=self.info)
        with open(fp) as f:
            head = f.readline()
        self.assertIn('gtdb_release=220', head)
        self.assertIn('version=', head)

    def test_external_tool_maps_stay_plain(self):
        # DIAMOND, BLAST and read_prot2taxid consume the protein maps
        # byte-for-byte, so they are root-level artifacts and never stamped.
        for name in ('taxon.map.gz', 'genome.map.gz'):
            self.assertFalse(isfile(join(self.d, 'gtdb', name)),
                             f'{name} must stay outside the stamped gtdb/ set')
            self.assertFalse(isfile(join(self.d, name)))


class TestWorkspaceChecker(unittest.TestCase):
    """IT-7: the compliance checker must catch each kind of violation."""

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools', 'verify_workspace.py')
        spec = spec_from_file_location('vw', path,
                                       loader=SourceFileLoader('vw', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _tree(self, stray=False, escape=False, missing=False):
        root = tmpdir()
        self.addCleanup(shutil.rmtree, root, ignore_errors=True)
        for name in self._tool().PROJECT_DIRS:
            os.makedirs(join(root, name))
        cache = join(root, 'neoHGT2GTDB-缓存测试')
        tmp = join(root, 'neoHGT2GTDB-临时文件')
        for sub in self._tool().CACHE_LAYOUT:
            if missing and sub == 'benchmarks':
                continue
            os.makedirs(join(cache, sub), exist_ok=True)
        os.makedirs(join(cache, 'db', 'd1', 'gtdb'), exist_ok=True)
        with open(join(cache, 'db', 'd1', 'meta.json'), 'w') as f:
            json.dump({'tool': 'neoHGT', 'workspace': {'tmp_dir': tmp}}, f)
        if stray:
            os.makedirs(join(root, 'someplace_else'))
        if escape:
            os.symlink('/etc', join(tmp, 'outside'))
        return root, cache, tmp

    def test_clean_tree_passes(self):
        vw = self._tool()
        root, cache, tmp = self._tree()
        seen, at_root, links = vw.check_project_files(root)
        self.assertEqual((at_root, links), ([], []))
        self.assertGreaterEqual(seen, 1)
        self.assertEqual(vw.check_cache_layout(cache), [])
        refs = vw.check_tmp_disposable(root, cache, tmp)
        self.assertEqual([r['path'] for r in refs], [tmp])   # inventory only

    def test_stray_directory_at_root_is_reported(self):
        vw = self._tool()
        root, cache, tmp = self._tree(stray=True)
        _, at_root, _ = vw.check_project_files(root)
        self.assertEqual(len(at_root), 1)
        self.assertIn('someplace_else', at_root[0]['path'])
        self.assertFalse(at_root[0]['predates_project'])

    def test_pre_existing_file_at_root_is_not_a_violation(self):
        vw = self._tool()
        root, cache, tmp = self._tree(stray=True)
        cutoff = time.time() + 3600          # everything predates this cutoff
        _, at_root, _ = vw.check_project_files(root, cutoff)
        self.assertTrue(at_root[0]['predates_project'])

    def test_symlink_escaping_the_workspace_is_reported(self):
        vw = self._tool()
        root, cache, tmp = self._tree(escape=True)
        _, _, links = vw.check_project_files(root)
        self.assertTrue(any('outside' in p for p in links), links)

    def test_missing_cache_subdirectory_is_reported(self):
        vw = self._tool()
        root, cache, tmp = self._tree(missing=True)
        self.assertEqual(vw.check_cache_layout(cache), ['benchmarks'])


class TestBlastRoute(unittest.TestCase):
    """ADR-8: makeblastdb cannot index 151-270 character GTDB ids."""

    def test_blast_refused_unless_forced(self):
        import argparse
        from types import SimpleNamespace
        from neoHGT.database import compile_gtdb_database
        me = SimpleNamespace(force=False, gtdb_max_id_len=151)
        for compile_ in ('blast', 'both'):
            with self.assertRaises(ValueError) as ctx:
                compile_gtdb_database(me, argparse.Namespace(
                    compile=compile_))
            self.assertIn('50 characters', str(ctx.exception))
        self.assertIsNone(compile_gtdb_database(
            me, argparse.Namespace(compile='none')))


class TestGtdbStyleFiltersKeepPrefixedNames(unittest.TestCase):
    """DoD-5: with the style filters on, a GTDB rank prefix must not cost hits.

    GTDB labels arrive with their `s__`/`g__` prefix in several places (user
    tables, `--input-tax`, dumps converted from the taxonomy files), and the
    upstream predicates judge the raw string, so `s__Escherichia coli` reads as
    "not capitalized" and every hit under it would vanish. The GTDB branch
    judges the prefix-stripped name instead.
    """

    def _predicates(self, keep=True, capital=True):
        from types import SimpleNamespace
        from unittest import mock
        from neoHGT import search as search_mod
        seen = {}

        def fake(self_, prots):
            seen['is_capital'] = search_mod.is_capital
            seen['is_latin'] = search_mod.is_latin
            return prots

        me = SimpleNamespace(gtdb_root='/db', tax_unirank=None,
                             tax_capital=capital, tax_latin=capital,
                             keep_placeholder_names=keep)
        with mock.patch.object(search_mod, '_orig_taxfilt_wf', fake):
            search_mod.taxfilt_wf_gtdb(me, {})
        return seen

    def test_prefixed_labels_survive_the_style_filters(self):
        seen = self._predicates()
        self.assertTrue(seen['is_capital']('s__Escherichia coli'))
        self.assertTrue(seen['is_latin']('s__Escherichia coli'))
        self.assertTrue(seen['is_capital']('g__Halobacterium'))
        self.assertTrue(seen['is_latin']('MK2026 sp010841923'))

    def test_plain_non_names_are_still_rejected(self):
        seen = self._predicates()
        self.assertFalse(seen['is_capital']('lowercase thing'))
        self.assertFalse(seen['is_latin']('Escherichia'))       # not a binomial
        self.assertFalse(seen['is_latin']('s__lowercase thing'))

    def test_gtdb_profile_defaults_the_filters_off(self):
        import yaml
        cfg = yaml.safe_load(open(join(HERE, os.pardir, 'config.yml')))
        prof = cfg['taxonomy']['profiles']['gtdb']
        self.assertFalse(prof['capital'])
        self.assertFalse(prof['latin'])
        self.assertTrue(prof['keep_placeholder'])


class TestPlaceholderExemption(unittest.TestCase):
    """UT-2.6/UT-3.5: --tax-capital/--tax-latin must not eat GTDB labels.

    The upstream filter loop looks the two predicates up in the search module,
    so the GTDB wrapper rebinds them for the duration of the call. These tests
    observe what the wrapped filter sees, rather than re-testing the filter.
    """

    def _predicates_seen(self, keep, capital=True):
        from types import SimpleNamespace
        from unittest import mock
        from neoHGT import search as search_mod
        seen = {}

        def fake(self_, prots):
            seen['is_capital'] = search_mod.is_capital
            seen['is_latin'] = search_mod.is_latin
            return prots

        me = SimpleNamespace(gtdb_root='/db', tax_unirank=None,
                             tax_capital=capital, tax_latin=False,
                             keep_placeholder_names=keep)
        with mock.patch.object(search_mod, '_orig_taxfilt_wf', fake):
            search_mod.taxfilt_wf_gtdb(me, {'p': []})
        return seen

    def test_exemption_is_active_when_requested(self):
        seen = self._predicates_seen(keep=True)
        # r220 has 1,658 taxa whose label starts with a digit; the plain
        # is_capital() would throw all of their hits away
        self.assertTrue(seen['is_capital']('2-12-FULL-45-22'))
        self.assertTrue(seen['is_capital']('Halobacterium'))
        self.assertFalse(seen['is_capital']('lowercase thing'))

    def test_strict_mode_still_tolerates_prefixes_but_not_placeholders(self):
        # strict means "do not exempt coded placeholder labels", not "reject
        # official GTDB names that carry a rank prefix" (plan DoD-5)
        seen = self._predicates_seen(keep=False)
        self.assertFalse(seen['is_capital']('2-12-FULL-45-22'),
                         '--keep-placeholder-names no must stop exempting')
        self.assertTrue(seen['is_capital']('s__Halobacterium'),
                        'a prefixed official name is still capitalized')
        self.assertFalse(seen['is_capital']('lowercase thing'))

    def test_predicates_are_restored_afterwards(self):
        from neoHGT import util
        from neoHGT import search as search_mod
        self._predicates_seen(keep=True)
        self.assertIs(search_mod.is_capital, util.is_capital)
        self.assertIs(search_mod.is_latin, util.is_latin)


class TestQualitySelectionOptions(unittest.TestCase):
    """`--type-material gtdb|ncbi`, `--typemater`, and why `--complete` refuses.

    Plan §7 gives the two type-material definitions different columns, and says
    the GTDB quality thresholds *replace* `--complete` rather than inheriting
    its assembly-level meaning. Before round 29 the `ncbi` branch read the GTDB
    column, `--typemater` did nothing in GTDB mode, and `--complete` was parsed
    and then silently ignored.
    """

    LIN = ('d__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G1;s__G1 alpha')

    def setUp(self):
        rows = [
            {'genome': 'BAC_1', 'accession': 'GB_GCA_1.1', 'accnov': 'GCA_1',
             'gtdb_taxonomy': self.LIN, 'domain': 'bacteria', 'is_rep': True,
             'completeness': 99.0, 'contamination': 0.3,
             'gtdb_type_species_of_genus': 't',
             'ncbi_type_material_designation': 'assembly from type material',
             'mimag_high_quality': 't'},
            {'genome': 'BAC_2', 'accession': 'GB_GCA_2.1', 'accnov': 'GCA_2',
             'gtdb_taxonomy': self.LIN, 'domain': 'bacteria', 'is_rep': True,
             'completeness': 60.0, 'contamination': 9.0,
             'gtdb_type_species_of_genus': 'f',
             'ncbi_type_material_designation': 'na',
             'mimag_high_quality': 'f'},
            {'genome': 'BAC_3', 'accession': 'GB_GCA_3.1', 'accnov': 'GCA_3',
             'gtdb_taxonomy': self.LIN, 'domain': 'bacteria', 'is_rep': False,
             'completeness': 95.0, 'contamination': 1.0,
             'gtdb_type_species_of_genus': 'f',
             'ncbi_type_material_designation': 'assembly from synonym type '
                                                'material',
             'mimag_high_quality': 'f'},
            {'genome': 'BAC_4', 'accession': 'GB_GCA_4.1', 'accnov': 'GCA_4',
             'gtdb_taxonomy': self.LIN, 'domain': 'bacteria', 'is_rep': False,
             'completeness': 95.0, 'contamination': 1.0,
             'gtdb_type_species_of_genus': 't',
             'ncbi_type_material_designation': '',
             'mimag_high_quality': 't'},
        ]
        self.df = pd.DataFrame(rows)
        p = gtdb.GtdbProvider.__new__(gtdb.GtdbProvider)
        p.taxid_map, _ = gtdb.assign_ids(gtdb.collect_nodes([self.LIN]))
        p.opts, p.selection_labels, p.sample_note = {}, [], None
        self.p = p

    def _genomes(self, opts):
        return list(self.p.select(self.df, opts)['genome'])

    def test_gtdb_definition_uses_the_gtdb_column(self):
        self.assertEqual(self._genomes({'type_material': 'gtdb'}),
                         ['BAC_1', 'BAC_4'])
        self.assertIn('GTDB type species', self.p.selection_labels[0])

    def test_ncbi_definition_uses_the_ncbi_column(self):
        got = self._genomes({'type_material': 'ncbi'})
        self.assertEqual(got, ['BAC_1', 'BAC_3'],
                         '"na" and empty must not count as a designation')
        self.assertIn('NCBI type material', self.p.selection_labels[0])

    def test_typemater_is_the_upstream_alias_for_the_ncbi_definition(self):
        self.assertEqual(self._genomes({'typemater': True}),
                         self._genomes({'type_material': 'ncbi'}))

    def test_complete_is_refused_with_the_replacements_named(self):
        with self.assertRaises(ValueError) as ctx:
            self._genomes({'complete': True})
        msg = str(ctx.exception)
        self.assertIn('--min-completeness', msg)
        self.assertIn('--max-contamination', msg)
        self.assertIn('assembly level', msg,
                      'the refusal must say why the upstream meaning does '
                      'not transfer, not just that it is unsupported')

    def test_quality_class_counts_are_reported(self):
        self.p.select(self.df, {})
        self.assertEqual(self.p.quality_classes, {
            'gtdb_type_species': 2, 'ncbi_type_material': 2,
            'mimag_high_quality': 2})

    def test_missing_column_says_so(self):
        slim = self.df.drop(columns=['ncbi_type_material_designation'])
        with self.assertRaises(ValueError) as ctx:
            self.p.select(slim, {'type_material': 'ncbi'})
        self.assertIn('ncbi_type_material_designation', str(ctx.exception))

    def test_gtdb_definition_missing_column_refuses(self):
        slim = self.df.drop(columns=['gtdb_type_species_of_genus'])
        with self.assertRaises(ValueError) as ctx:
            self.p.select(slim, {'type_material': 'gtdb'})
        self.assertIn('gtdb_type_species_of_genus', str(ctx.exception))

    def test_unknown_definition_refuses(self):
        with self.assertRaises(ValueError) as ctx:
            self._genomes({'type_material': 'refseq'})
        self.assertIn('not one of', str(ctx.exception))

    def test_mimag_missing_column_refuses(self):
        # Not every release publishes MIMAG classes (GTDB-Tk input tables do
        # not either). Before round 29 the filter was `if col in sel`, so
        # --mimag high on such a table quietly returned every genome.
        slim = self.df.drop(columns=['mimag_high_quality'])
        self.assertEqual(self._genomes({'mimag': 'high'}), ['BAC_1', 'BAC_4'],
                         'with the column present it must filter')
        with self.assertRaises(ValueError) as ctx:
            self.p.select(slim, {'mimag': 'high'})
        self.assertIn('mimag_high_quality', str(ctx.exception))

    def test_mimag_unknown_class_refuses(self):
        with self.assertRaises(ValueError) as ctx:
            self._genomes({'mimag': 'ultra'})
        self.assertIn('high / medium / low', str(ctx.exception))

    def test_any_tier_keeps_genomes_carrying_one_tier(self):
        wide = self.df.assign(mimag_high_quality=['t', 'f', 'f', 'f'],
                              mimag_medium_quality=['f', 't', 'f', 'f'],
                              mimag_low_quality=['f', 'f', 't', 'f'])
        self.assertEqual(list(self.p.select(wide, {'mimag': 'any'})['genome']),
                         ['BAC_1', 'BAC_2', 'BAC_3'],
                         '"any" must accept a genome carrying *some* tier and '
                         'still drop the unclassified one')
        self.assertEqual(self.p.selection_labels, ['MIMAG-any-tier'])

    def test_any_tier_without_the_other_columns_refuses(self):
        # the fixture table only publishes the high column
        with self.assertRaises(ValueError) as ctx:
            self._genomes({'mimag': 'any'})
        self.assertIn('mimag_medium_quality', str(ctx.exception))


class TestResolveLabel(unittest.TestCase):
    """The plan's Appendix B helper `resolve_label(label, taxid_map)`.

    It replaced three inline `next(k for k, v in labels.items() if v == ref)`
    scans -- one per donor row against a 146k-node label table. The cached
    reverse index must never go stale when a different database is passed in.
    """

    def test_numeric_ids_pass_through_and_labels_resolve(self):
        labels = {'100001': 'g__Escherichia', '100002': 'g__Salmonella'}
        self.assertEqual(gtdb.resolve_label('100001', labels), '100001')
        self.assertEqual(gtdb.resolve_label(' g__Salmonella ', labels),
                         '100002')
        self.assertIsNone(gtdb.resolve_label('g__Nobody', labels))
        self.assertIsNone(gtdb.resolve_label('', labels))
        self.assertIsNone(gtdb.resolve_label(None, labels))

    def test_taxid_map_is_the_fallback(self):
        labels = {'100001': 'g__Escherichia'}
        tmap = {'d__Bacteria;g__Wigglesworthia': '100500'}
        self.assertEqual(
            gtdb.resolve_label('d__Bacteria;g__Wigglesworthia', labels,
                               taxid_map=tmap), '100500')
        self.assertIsNone(gtdb.resolve_label('d__Bacteria;g__Missing', labels,
                                            taxid_map=tmap))

    def test_duplicate_labels_keep_the_first_taxid(self):
        # same tie-breaking the inline scan had
        labels = {'100001': 'g__Same', '100002': 'g__Same'}
        self.assertEqual(gtdb.resolve_label('g__Same', labels), '100001')

    def test_the_index_follows_a_different_mapping(self):
        a = {'1': 'g__A'}
        b = {'2': 'g__A'}
        self.assertEqual(gtdb.resolve_label('g__A', a), '1')
        self.assertEqual(gtdb.resolve_label('g__A', b), '2',
                         'a cached index from another database must not leak')
        self.assertEqual(gtdb.resolve_label('g__A', a), '1')

    def test_the_call_sites_use_it(self):
        # map_report and the donor mapping in analyze go through this helper
        import inspect
        from neoHGT import analyze
        self.assertIn('resolve_label',
                      inspect.getsource(gtdb.map_report))
        self.assertIn('resolve_label',
                      inspect.getsource(analyze.write_hgt_list_gtdb))
        self.assertNotIn('labels.items() if v ==',
                         inspect.getsource(gtdb.map_report))


class TestNcbiProviderAdapter(unittest.TestCase):
    """The plan's Appendix B promises an `NcbiProvider` beside the GTDB one.

    It is a wrapper, not a reimplementation: whatever the upstream pipeline
    already holds must come back unchanged, and nothing may be read twice. A
    round-28 cleanup deleted it and the register gate rejected that, which is
    the right outcome -- a promised symbol has to stay, so this pins what
    "thin wrapper" actually means.
    """

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        lins = sorted(set(gtdb.read_taxonomy(
            [join(FX, 'taxonomy_sample.tsv')]).values()))
        td, tmap, stats = gtdb.build_taxdump_from_gtdb(lins)
        gtdb.write_taxdump(td, self.d)
        self.td = td
        self.p = gtdb.NcbiProvider

    def test_it_offers_the_same_interface_as_the_gtdb_provider(self):
        for meth in ('prepare', 'load_taxonomy', 'load_genomes', 'select',
                     'sample', 'protein_sources'):
            self.assertTrue(callable(getattr(self.p, meth, None)), meth)
        prov = self.p(taxdump_dir=self.d)
        self.assertEqual(prov.name, 'ncbi')
        self.assertEqual(prov.ranks, gtdb.NCBI_RANKS)
        self.assertTrue(prov.default_taxonomy_opts['capital'],
                        'NCBI labels are capitalised; GTDB ones are not')

    def test_taxonomy_is_read_once_and_returned_as_is(self):
        prov = self.p(taxdump_dir=self.d)
        tax, tmap, stats = prov.load_taxonomy()
        self.assertEqual(set(tax), set(self.td))
        self.assertEqual(tmap, {})
        self.assertEqual(stats['n_nodes'], len(self.td))
        again, _, s2 = prov.load_taxonomy()
        self.assertIs(again, tax, 'a second call must not re-read the files')
        self.assertEqual(s2['n_nodes'], len(self.td))

    def test_wrapping_the_upstream_object_keeps_behaviour_unchanged(self):
        from types import SimpleNamespace
        df = object()
        inner = SimpleNamespace(taxdump={'1': {'parent': '1'}}, df=df)
        prov = self.p(wrap=inner)
        tax, _, _ = prov.load_taxonomy()
        self.assertIs(tax, inner.taxdump)
        self.assertIs(prov.load_genomes(), df)
        # the tables pass straight through: no ranking, no sampling
        self.assertIs(prov.select(df), df)
        self.assertIs(prov.sample(df), df)
        self.assertEqual(list(prov.protein_sources(df)), [])

    def test_a_source_must_be_given(self):
        with self.assertRaises(ValueError) as ctx:
            self.p().load_taxonomy()
        self.assertIn('taxdump', str(ctx.exception))


class TestSearchGtdbModeDetection(unittest.TestCase):
    """`gtdb_db_dir` / `is_gtdb_mode` decide every GTDB-only behaviour in search.

    `is_gtdb_mode` cannot be deleted: search.py is an upstream file held to the
    pure-append rule, so an unused helper has to be made *tested* instead of
    removed. Round 28's execution tracer found it unexecuted, which is also how
    a wrong answer here would silently disable the whole GTDB branch.
    """

    @staticmethod
    def _ns(**kw):
        from types import SimpleNamespace
        base = dict(source='auto', gtdb_stem2taxid=None, taxdump=None,
                    taxmap=None, input=None, output=None, db=None)
        base.update(kw)
        return SimpleNamespace(**base)

    def test_a_gtdb_database_is_detected_from_any_hint(self):
        from neoHGT.search import gtdb_db_dir, is_gtdb_mode
        inside = join(FX, 'gtdb', 'labels.tsv')
        for attr in ('taxdump', 'input', 'output', 'db'):
            me = self._ns(**{attr: inside})
            self.assertEqual(gtdb_db_dir(me), FX, f'{attr} hint')
            # detection happens once, in set_params_gtdb; everything after that
            # asks the predicate, which reads the resolved value
            self.assertTrue(is_gtdb_mode(self._ns(gtdb_root=gtdb_db_dir(me))),
                            attr)

    def test_a_plain_directory_is_not_a_gtdb_database(self):
        from neoHGT.search import gtdb_db_dir, is_gtdb_mode
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        me = self._ns(output=d, taxdump=d)
        self.assertIsNone(gtdb_db_dir(me))
        self.assertFalse(is_gtdb_mode(self._ns(gtdb_root=gtdb_db_dir(me))))

    def test_mode_checks_go_through_the_one_predicate(self):
        # four places used to re-ask `getattr(self, 'gtdb_root', None)` by hand;
        # the reachability audit found `is_gtdb_mode` had no caller at all, which
        # is how a named concept and its duplicated checks drift apart.
        # Matched on the syntax tree, not on source text: the fifth site wrapped
        # its call across two lines and a literal-string assertion had sailed
        # straight past it -- the one failure mode this rule must not have.
        with open(join(HERE, os.pardir, 'search.py'),
                  encoding='utf-8') as fh:
            src = fh.read()
        self.assertIn('if not is_gtdb_mode(self):', src)
        self.assertIn('if is_gtdb_mode(self) and self.tax_unirank:', src)
        tree = ast.parse(src)
        pred = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                    and n.name == 'is_gtdb_mode')
        asks = [n.lineno for n in ast.walk(tree)
                if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                and n.func.id == 'getattr' and len(n.args) > 1
                and isinstance(n.args[0], ast.Name) and n.args[0].id == 'self'
                and isinstance(n.args[1], ast.Constant)
                and n.args[1].value == 'gtdb_root']
        dup = [ln for ln in asks if not pred.lineno <= ln <= pred.end_lineno]
        self.assertEqual(dup, [], f'is_gtdb_mode() bypassed at {dup}')

    def test_the_wrapped_call_is_still_a_bypass(self):
        # the shape that fooled the original text-matching rule: a second
        # function calls getattr with the option string on one line and the
        # default on the next. The AST rule must flag it; the text rule must
        # miss it -- that difference is the whole reason for the rewrite.
        wrapped = ("def is_gtdb_mode(self):\n"
                   "    return bool(getattr(self, 'gtdb_root', None))\n"
                   "\n"
                   "def other(self):\n"
                   "    if getattr(self, 'gtdb_root',\n"
                   "                   None) and self.x:\n"
                   "        return 1\n")
        tree = ast.parse(wrapped)
        pred = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                    and n.name == 'is_gtdb_mode')
        asks = [n.lineno for n in ast.walk(tree)
                if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                and n.func.id == 'getattr' and len(n.args) > 1
                and isinstance(n.args[0], ast.Name) and n.args[0].id == 'self'
                and isinstance(n.args[1], ast.Constant)
                and n.args[1].value == 'gtdb_root']
        self.assertEqual(asks, [2, 5], asks)
        self.assertEqual([ln for ln in asks
                          if not pred.lineno <= ln <= pred.end_lineno], [5],
                         'the AST scanner cannot see a wrapped call')
        self.assertNotIn("getattr(self, 'gtdb_root', None)",
                         wrapped.replace(
                             "return bool(getattr(self, 'gtdb_root', None))",
                             ''),
                         'the old text rule would have found nothing either, '
                         'so it cannot be the guard')

    def test_the_ncbi_source_short_circuits(self):
        from neoHGT.search import gtdb_db_dir, is_gtdb_mode
        me = self._ns(source='ncbi', db=FX)
        self.assertIsNone(gtdb_db_dir(me),
                          'an --source ncbi run must never take the GTDB branch')


class TestGtdbGroupingProfile(unittest.TestCase):
    """The per-domain grouping profile has to reach the object.

    Round 40 found it never did, twice over: `apply_gtdb_grouping` looked up
    `grouping.profiles['closesize']` (the map, not the selected profile), and
    even after that, `set_parameters` had already merged the generic
    `grouping.closesize: 10` into the object -- so the documented bacteria 200 /
    archaea 50 defaults could not apply. The measured consequence: runs without
    an explicit --close-size silently used 10.
    """

    def _cfg(self):
        import yaml
        return yaml.safe_load(open(join(HERE, os.pardir, 'config.yml')))

    def _analyse(self, domain, close_size_cli=None, input_cov_cli=None):
        from types import SimpleNamespace
        from neoHGT import analyze as A
        me = SimpleNamespace(
            cfg=self._cfg(), close_size=close_size_cli,
            input_cov=input_cov_cli, self_rank=None, donor_rank=None,
            tax_unirank=None,
            taxdump={'1': {'parent': '1', 'rank': 'no rank',
                           'name': 'root'}},
            cli_grouping={'close_size': close_size_cli,
                          'input_cov': input_cov_cli})
        orig, A.domain_of_input = A.domain_of_input, lambda self: domain
        try:
            A.apply_gtdb_grouping(me)
        finally:
            A.domain_of_input = orig
        return me

    def test_the_archaeal_profile_applies(self):
        self.assertEqual(self._analyse('archaea').close_size, 50)

    def test_the_bacterial_profile_applies(self):
        self.assertEqual(self._analyse('bacteria').close_size, 200)

    def test_an_explicit_option_beats_the_profile(self):
        # the falsifying half: without a CLI-vs-config distinction the profile
        # would simply overwrite whatever the user asked for
        me = self._analyse('archaea', close_size_cli=7)
        self.assertEqual(me.close_size, 7)
        self.assertNotEqual(me.close_size, 50)

    def test_the_generic_grouping_value_is_not_mistaken_for_the_profile(self):
        # config.yml says grouping.closesize = 10; a run that gets 10 has not
        # had its profile applied
        self.assertEqual(self._cfg()['grouping']['closesize'], 10)
        self.assertNotEqual(self._analyse('bacteria').close_size, 10)

    def test_set_parameters_records_what_the_command_line_actually_said(self):
        from types import SimpleNamespace
        from unittest import mock
        from neoHGT import analyze as A
        me = SimpleNamespace()
        args = SimpleNamespace(close_size=42, input_cov=None)
        called = []
        with mock.patch.object(A, '_orig_set_parameters',
                               lambda self, a: called.append((self, a))):
            A.set_parameters_gtdb(me, args)
        self.assertEqual(len(called), 1, 'the upstream call has to still happen')
        self.assertEqual(me.cli_grouping, {'close_size': 42, 'input_cov': None},
                         'what the user passed must survive the config merge')

    # a two-domain taxonomy tree in the shape read_taxdump() returns
    DUMP = {
        '1': {'parent': '1', 'rank': 'no rank', 'name': 'root'},
        '2': {'parent': '1', 'rank': 'superkingdom', 'name': 'Archaea'},
        '21': {'parent': '2', 'rank': 'phylum', 'name': 'Halobacteriota'},
        '22': {'parent': '21', 'rank': 'family', 'name': 'Halobacteriaceae'},
        '100': {'parent': '1', 'rank': 'superkingdom', 'name': 'Bacteria'},
        '101': {'parent': '100', 'rank': 'phylum', 'name': 'Actinomycetota'},
        '102': {'parent': '101', 'rank': 'family', 'name': 'Strepto_f'},
    }

    def _base(self, **kw):
        from types import SimpleNamespace
        me = SimpleNamespace(cfg=self._cfg(), close_size=None, input_cov=None,
                             self_rank=None, donor_rank=None,
                             tax_unirank=None, taxdump=self.DUMP,
                             cli_grouping={'close_size': None,
                                           'input_cov': None},
                             data=None, input_tax=None, df=None,
                             from_scores=False, output='out')
        for key, val in kw.items():
            setattr(me, key, val)
        return me

    def _run(self, me):
        """apply_gtdb_grouping() and everything it printed."""
        import contextlib
        import io
        from neoHGT import analyze as A
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            A.apply_gtdb_grouping(me)
        return me, buf.getvalue()

    def _both_routes(self, donors):
        """The same taxa through the two analyze routes: `data` (full run) and
        a reused `df` (--from-scores)."""
        prots = [{'id': f'p{i}', 'hits': pd.DataFrame({'taxid': [d]})}
                 for i, d in enumerate(donors)]
        full = self._run(self._base(data={'query': prots}))
        reuse = self._run(self._base(
            df=pd.DataFrame(dict(match=donors)), from_scores=True))
        return full, reuse

    def test_both_routes_pick_the_same_profile_for_the_same_taxa(self):
        # this is the whole point of the round: --from-scores used to answer
        # "bacterial, close-size 200" for the archaeal sample a full run had
        # just reported as archaea/50
        full, reuse = self._both_routes(['22'] * 5 + ['21'] * 3)
        self.assertEqual(full[0].close_size, 50)
        self.assertEqual(reuse[0].close_size, 50, reuse[1])
        self.assertIn('using "gtdb_archaea"', reuse[1])
        self.assertEqual(reuse[0].domain_hint, full[0].domain_hint)

    def test_the_reuse_route_reads_a_numeric_donor_column(self):
        # scores.tsv holds donors as bare numbers, so what comes back is int64
        # -- or float64 as soon as one protein has no donor at all
        for donors, dtype in (([22, 22, 21], None),
                              ([22, 22, None], 'float64')):
            reuse = self._run(self._base(df=pd.DataFrame({'match': donors}),
                                         from_scores=True))
            self.assertEqual(reuse[0].close_size, 50, donors)
            if dtype:
                self.assertEqual(reuse[0].df['match'].dtype, dtype)

    def test_a_bacterial_score_table_still_gets_the_bacterial_profile(self):
        full, reuse = self._both_routes(['102'] * 4 + ['101'] * 2)
        self.assertEqual(full[0].close_size, 200)
        self.assertEqual(reuse[0].close_size, 200)
        self.assertNotIn('gtdb_archaea', reuse[1])

    def test_the_reuse_route_says_the_grouping_knobs_are_inert(self):
        # the values above the note are true but they change nothing here:
        # --from-scores skips define_groups(), the only reader of close_size
        full, reuse = self._both_routes(['22'] * 3)
        self.assertIn('--from-scores reuses the score table', reuse[1])
        self.assertNotIn('--from-scores reuses the score table', full[1])
        self.assertEqual(full[0].close_size, reuse[0].close_size)

    def test_an_unresolvable_domain_is_reported_as_a_fallback(self):
        for donors in ([None, None], ['999999', '999998']):
            _, log = self._run(self._base(
                df=pd.DataFrame({'match': donors}), from_scores=True))
            self.assertIn('input domain = not determined', log, donors)
            self.assertIn('falling back to the default profile', log, donors)

    def test_a_half_and_half_input_claims_neither_profile(self):
        # a tie used to be resolved by whatever the first protein happened to
        # be, which silently picked the 4x smaller archaeal close group
        for donors in (['22', '102'], ['22', '102', '21', '101']):
            _, log = self._run(self._base(
                df=pd.DataFrame({'match': donors}), from_scores=True))
            self.assertIn('not determined', log, donors)
            self.assertNotIn('gtdb_archaea', log, donors)

    def test_the_donor_column_is_what_makes_the_reuse_route_work(self):
        # falsifying half: strip the score table and the route is blind again,
        # so the agreement above is earned by the vote and not by luck
        from neoHGT import analyze as A
        me = self._base(df=pd.DataFrame({'match': ['22', '21', '22']}),
                        from_scores=True)
        self.assertEqual(A.domain_of_input(me), 'archaea')
        me.df = pd.DataFrame({'match': []})
        self.assertIsNone(A.domain_of_input(me))


class TestSourceHonouredOnBothRoutes(unittest.TestCase):
    """`--source` cannot mean two things depending on which route ran.

    Round 41's inert-options sweep (`缓存测试/reports/R41_inert_options.json`)
    found `source` read in exactly one place -- `read_input()`, the step
    `--from-scores` skips. The reuse route therefore still held the taxdump
    *path* (the dict swap lives in that same skipped step), which is precisely
    the string `find_db_root()` walks up to the GTDB marker: measured on the
    archaeal sample, `--source ncbi` gave NCBI mode on the full run and GTDB
    mode on the reused one (`临时文件/r49_src_ncbi_full.log` vs
    `r49_src_ncbi_reuse.log`).
    """

    def _db(self, marker=True):
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        if marker:
            os.makedirs(join(d, 'gtdb'), exist_ok=True)
            with open(join(d, 'gtdb', 'labels.tsv'), 'w') as f:
                f.write('100001\tg__Example\n')
            # a real build root always carries meta.json naming its source;
            # find_db_root requires that corroboration since the nested-NCBI
            # misdetection (review C 4.11)
            with open(join(d, 'meta.json'), 'w') as f:
                f.write('{"source": "gtdb"}')
        else:
            with open(join(d, 'db.faa'), 'w') as f:
                f.write('>x\n')
        return d

    def _me(self, source, taxdump):
        from types import SimpleNamespace
        from neoHGT import analyze as A
        return A, SimpleNamespace(source=source, taxdump=taxdump,
                                  input=None, output=None)

    def test_ncbi_source_never_looks_for_the_marker(self):
        A, me = self._me('ncbi', self._db())
        self.assertIsNone(A.detect_gtdb_root(me))
        self.assertIsNone(getattr(me, 'gtdb_root', None),
                          'the option has to decide the mode, not the layout')

    def test_auto_source_finds_the_marker_the_same_way_on_both_routes(self):
        db = self._db()
        A, reuse = self._me('auto', db)        # reuse route: still a path
        A, full = self._me('auto', {'1': {'parent': '1', 'rank': 'no rank',
                                          'name': 'root'}})
        # by the time the full route reaches prediction, the path is gone
        self.assertEqual(A.detect_gtdb_root(reuse), db)
        self.assertEqual(A.detect_gtdb_root(full, [db]), db)

    def test_explicit_gtdb_source_refuses_to_silently_downgrade(self):
        A, me = self._me('gtdb', self._db(marker=False))
        with self.assertRaises(ValueError) as ctx:
            A.detect_gtdb_root(me)
        msg = str(ctx.exception)
        self.assertIn('gtdb/labels.tsv', msg)
        # the message must name the alternative, not just the failure
        self.assertIn('--source', msg)

    def test_hybrid_source_refuses_to_silently_downgrade(self):
        # the search stage raises for a missing marker; analyse used to print
        # a WARNING and run the NCBI path, so one flag meant two behaviours
        A, me = self._me('hybrid', self._db(marker=False))
        with self.assertRaises(ValueError) as ctx:
            A.detect_gtdb_root(me)
        self.assertIn('--source', str(ctx.exception))

    def test_defaulted_option_counts_as_auto(self):
        A, me = self._me(None, self._db())
        self.assertEqual(A.source_of(me), 'auto')
        self.assertEqual(A.detect_gtdb_root(me), me.taxdump)

    def _predict(self, source, db, marker=True):
        """Run the prediction wrapper and report what it did to the object."""
        import contextlib
        import io
        from types import SimpleNamespace
        from unittest import mock
        from neoHGT import analyze as A
        me = SimpleNamespace(source=source, taxdump=db, input=None,
                             output=None, bandwidth=None,
                             df=pd.DataFrame({'match': [22, 21],
                                              'close': [1.0, 2.0],
                                              'distal': [3.0, 4.0],
                                              'self': [0.0, 0.0]}),
                             gtdb_labels=None)
        calls = []
        dump = {'22': {'parent': '2', 'rank': 'family', 'name': 'f1'}}
        buf = io.StringIO()
        with mock.patch.object(A, '_orig_predict_hgt',
                               lambda self: calls.append(self) or 0), \
                mock.patch.object(A, 'read_taxdump', lambda p: dump), \
                mock.patch.object(A, 'apply_gtdb_grouping',
                                  lambda self: calls.append('grouping')):
            with contextlib.redirect_stdout(buf):
                A.predict_hgt_gtdb(me)
        return me, calls, buf.getvalue(), (marker and db)

    def test_explicit_source_repairs_the_skipped_read_input(self):
        # --from-scores skips the step that turns the taxdump path into the
        # taxonomy dict; a run that named its mode explicitly gets that repair
        me, calls, log, _ = self._predict('ncbi', self._db())
        self.assertIsInstance(me.taxdump, dict, 'the path became the taxonomy')
        self.assertEqual(list(me.taxdump), ['22'])
        self.assertEqual(list(me.df['match']), ['22', '21'],
                         'donors come back as the string keys the dict uses')
        self.assertIn('Reuse route: reading the taxonomy', log)
        self.assertNotIn('grouping', calls,
                         'NCBI mode must not take the GTDB grouping branch')

    def test_no_source_means_upstream_behaviour_unchanged(self):
        # the RT-1 boundary: without an explicit --source, a non-GTDB database on
        # the reuse route is left exactly as upstream leaves it (crash included)
        me, calls, log, _ = self._predict('auto', self._db(marker=False),
                                          marker=False)
        self.assertIsInstance(me.taxdump, str)
        self.assertNotIn('Reuse route: reading the taxonomy', log)
        self.assertEqual(calls, [me], 'only the upstream prediction ran')


class TestReuseRouteInertOptions(unittest.TestCase):
    """The options `--from-scores` cannot honour have to be named, and checked.

    Round 41 fixed the two grouping knobs and then measured the whole class:
    `neoHGT/analyze.py::REUSE_INERT` is what the run prints, and
    `tools/verify_plan_traceability.py::reuse_inert_problems` re-derives it
    from the call graph, so the list cannot rot in either direction.
    """

    @staticmethod
    def _tool(path=None):
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        if path is None:
            path = join(HERE, os.pardir, os.pardir, 'tools',
                        'verify_plan_traceability.py')
        tail = 'verify_plan_traceability.py'
        name = 'vpt7' if path.endswith(tail) else 'vpt7x'
        spec = spec_from_file_location(name, path,
                                       loader=SourceFileLoader(name, path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _me(self, **kw):
        """A reuse-route-shaped object: no hits, only a score table."""
        import yaml
        from types import SimpleNamespace
        cfg = yaml.safe_load(open(join(HERE, os.pardir, 'config.yml')))
        me = SimpleNamespace(cfg=cfg, close_size=None, input_cov=None,
                             self_rank=None, donor_rank=None,
                             tax_unirank=None,
                             taxdump=TestGtdbGroupingProfile.DUMP,
                             cli_grouping={'close_size': None,
                                           'input_cov': None},
                             cli_inert={}, data=None, input_tax=None,
                             df=pd.DataFrame({'match': ['22', '21']}),
                             from_scores=True, output='out')
        for key, val in kw.items():
            setattr(me, key, val)
        return me

    def test_the_shipped_declaration_is_what_the_graph_derives(self):
        v = self._tool()
        report, problems = v.reuse_inert_problems(
            join(HERE, os.pardir, os.pardir))
        self.assertEqual(problems, [], report)
        self.assertGreaterEqual(len(report['declared']), 10, report)
        self.assertEqual(set(report['declared']),
                         set(report['derived']) | set(report['exceptions']),
                         'every declared-dead option is either derived or '
                         'excused, and nothing is excused without needing it')
        # `--source` was in this list until round 42 made both routes read it
        self.assertNotIn('source', report['declared'])
        # and the option that report-time still needs is not swept in either
        self.assertNotIn('donor_rank', report['declared'])

    def test_an_unDeclared_dead_option_is_named(self):
        v = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        for sub in ('neoHGT', 'tools'):
            shutil.copytree(join(HERE, os.pardir, os.pardir, sub),
                            join(d, sub))
        path = join(d, 'neoHGT', 'analyze.py')
        text = open(path, encoding='utf-8').read()
        # drop one honest entry: the graph still says the option is dead, so
        # the user would be left silent about a knob that cannot work
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text.replace("'taxon_synonyms', ", '', 1))
        report, problems = v.reuse_inert_problems(d)
        self.assertTrue(any('taxon_synonyms' in p for p in problems), problems)

    def test_a_declared_option_that_is_actually_live_is_named(self):
        v = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        for sub in ('neoHGT', 'tools'):
            shutil.copytree(join(HERE, os.pardir, os.pardir, sub),
                            join(d, sub))
        path = join(d, 'neoHGT', 'analyze.py')
        text = open(path, encoding='utf-8').read()
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text.replace(
                "'weighted')", "'weighted', 'donor_rank')", 1))
        report, problems = v.reuse_inert_problems(d)
        self.assertTrue(any('donor_rank' in p for p in problems), problems)

    def test_a_stale_exemption_is_named(self):
        v = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        for sub in ('neoHGT', 'tools'):
            shutil.copytree(join(HERE, os.pardir, os.pardir, sub),
                            join(d, sub))
        tool = join(d, 'tools', 'verify_plan_traceability.py')
        text = open(tool, encoding='utf-8').read()
        # an exemption for an option that is inert anyway excuses nothing --
        # and would quietly cover up the next change to that option
        with open(tool, 'w', encoding='utf-8') as fh:
            fh.write(text.replace("    'self_rank': 'only rewritten",
                                  "    'weighted': 'excused for no reason',\n"
                                  "    'self_rank': 'only rewritten", 1))
        mod = self._tool(tool)
        report, problems = mod.reuse_inert_problems(d)
        self.assertTrue(any('weighted' in p for p in problems), problems)

    def test_a_blind_graph_is_not_a_green_gate(self):
        # the first traversal walked only module-level functions, never saw a
        # method, and so "proved" that nothing is inert. An empty walk and a
        # fixed product have to be distinguishable.
        v = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'neoHGT'), exist_ok=True)
        target = join(d, 'neoHGT', 'analyze.py')
        with open(target, 'w', encoding='utf-8') as fh:
            fh.write('arguments = []\n')
        report, problems = v.reuse_inert_problems(d)
        self.assertTrue(problems, report)
        self.assertTrue(any('一个选项都没读到' in p for p in problems), problems)
        self.assertTrue(any('没看见方法' in p for p in problems), problems)

    def test_the_note_names_the_options_the_user_passed(self):
        import contextlib
        import io
        from neoHGT import analyze as A
        me = self._me(cli_inert={'weighted': True, 'self_rank': 'genus',
                                 'evalue': 1e-5, 'donor_rank': 'genus',
                                 'taxon_synonyms': None})
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            A.apply_gtdb_grouping(me)
        log = buf.getvalue()
        self.assertIn('NOTE: --from-scores reuses the score table', log)
        for flag in ('--weighted', '--self-rank', '--evalue'):
            self.assertIn(flag, log, log)
        # options this route can still honour, and ones that were never
        # passed, must not be swept into the warning
        self.assertNotIn('--donor-rank', log, log)
        self.assertNotIn('--taxon-synonyms', log, log)

    def test_the_full_route_does_not_print_the_reuse_note(self):
        import contextlib
        import io
        from neoHGT import analyze as A
        hits = pd.DataFrame({'taxid': ['22']})
        me = self._me(from_scores=False, df=None,
                      data={'query': [{'id': 'p1', 'hits': hits}]},
                      cli_inert={'weighted': True})
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            A.apply_gtdb_grouping(me)
        self.assertNotIn('reuses the score table', buf.getvalue())

    def test_input_tax_is_read_but_only_as_a_dict(self):
        # the exemption's reason, tested rather than asserted in a comment: on
        # this route the value arrives as the path string the CLI carried
        from neoHGT import analyze as A
        me = self._me()
        me.input_tax = join('somewhere', 'input_tax.tsv')
        self.assertEqual(A.domain_of_input(me), 'archaea',
                         'the donor vote is what decided it, as it would with '
                         'no --input-tax at all')
        me.input_tax = {'query': '102'}        # bacteria, if the dict were read
        self.assertEqual(A.domain_of_input(me), 'bacteria')

class TestDocExampleValues(unittest.TestCase):
    """Rule 18: a documented value must exist on the stage it is documented for.

    The flag-level rule already asks "does this sub-command have `--source`?",
    which is blind to the fact that database accepts ncbi|gtdb|hybrid while
    search and analyze also accept `auto`: `neoHGT database --source auto` is a
    plausible example that argparse rejects.
    """

    NEED = ('neoHGT', 'tools', 'scripts', 'docs')

    def _tool(self, path=None):
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        if path is None:
            path = join(HERE, os.pardir, os.pardir, 'tools',
                        'verify_plan_traceability.py')
        name = 'vpt8' if path.endswith('verify_plan_traceability.py') else 'vpt8x'
        spec = spec_from_file_location(name, path,
                                       loader=SourceFileLoader(name, path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _tree(self):
        code = join(HERE, os.pardir, os.pardir)
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        for sub in self.NEED:
            shutil.copytree(join(code, sub), join(d, sub))
        shutil.copy(join(code, 'README.md'), join(d, 'README.md'))
        return d

    def test_shipped_examples_use_values_their_stage_accepts(self):
        v = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        report, problems = v.doc_value_problems(code, v.cli_subcommands(code))
        self.assertEqual(problems, [], report)
        self.assertEqual(report['unbalanced_fences'], [], report)
        self.assertGreaterEqual(report['judged'], 30, report)

    def test_a_value_the_stage_does_not_accept_is_named(self):
        v = self._tool()
        d = self._tree()
        path = join(d, 'README.md')
        text = open(path, encoding='utf-8').read()
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text + '\n```bash\nneoHGT database --source auto -o db\n'
                         '```\n')
        report, problems = v.doc_value_problems(d, v.cli_subcommands(d))
        self.assertTrue(any('--source auto' in p for p in problems), problems)

    def test_an_alternation_is_judged_member_by_member(self):
        v = self._tool()
        d = self._tree()
        path = join(d, 'README.md')
        text = open(path, encoding='utf-8').read()
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text + '\n```bash\nneoHGT gtdb --action probe|nonsense'
                             '\n```\n')
        report, problems = v.doc_value_problems(d, v.cli_subcommands(d))
        joined = ' '.join(problems)
        self.assertIn('nonsense', joined, problems)
        self.assertNotIn('`neoHGT gtdb --action probe`', joined, problems)

    def test_placeholders_are_counted_not_judged(self):
        v = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        report, problems = v.doc_value_problems(code, v.cli_subcommands(code))
        self.assertGreaterEqual(len(report['placeholders']), 5, report)
        self.assertTrue(any('$' in p[3] for p in report['placeholders']),
                        report['placeholders'][:5])
        self.assertEqual(problems, [])

    def test_prose_mentioning_a_value_is_not_an_example(self):
        # the fence scoping is what keeps this rule from crying wolf: the probe
        # version that scanned every line "found" six violations in clean docs
        v = self._tool()
        d = self._tree()
        path = join(d, 'docs', 'GTDB-数据库.md')
        text = open(path, encoding='utf-8').read()
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text + '\n也可以写 `neoHGT database --source wrongvalue`，'
                            '但那是散文不是示例。\n')
        report, problems = v.doc_value_problems(d, v.cli_subcommands(d))
        self.assertEqual(problems, [], problems)

    def test_an_unbalanced_fence_is_reported_as_blindness(self):
        v = self._tool()
        d = self._tree()
        path = join(d, 'docs', '交付物清单.md')
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write('```bash\nneoHGT analyze --source gtdb -i in -o out\n')
        report, problems = v.doc_value_problems(d, v.cli_subcommands(d))
        self.assertIn('交付物清单.md', report['unbalanced_fences'], report)
        self.assertTrue(any('围栏不成对' in p for p in problems), problems)

    def test_the_rule_refuses_an_empty_walk(self):
        v = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'docs'), exist_ok=True)
        with open(join(d, 'README.md'), 'w', encoding='utf-8') as fh:
            fh.write('# nothing here\n')
        report, problems = v.doc_value_problems(d, ('database', 'search',
                                                    'analyze', 'gtdb'))
        self.assertTrue(any('空跑' in p for p in problems), problems)

    def test_the_plan_document_is_judged_by_the_same_reader(self):
        # the plan is the source the code is measured against; a wrong value in
        # it is inherited by every copy, and round 43 only ever measured it by
        # hand
        v = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        with open(join(d, 'PLAN.md'), 'w', encoding='utf-8') as fh:
            fh.write('# 方案\n\n```bash\n'
                     'neoHGT database --source automatic -o db\n'
                     'neoHGT database --seq-source local -o db\n```\n')
        report, problems = v.doc_value_problems(
            join(HERE, os.pardir, os.pardir), v.cli_subcommands(
                join(HERE, os.pardir, os.pardir)), join(d, 'PLAN.md'))
        self.assertIn('automatic', ' '.join(problems), problems)
        self.assertNotIn('local', ' '.join(problems), problems)
        self.assertGreaterEqual(report['plan_judged'], 2, report)
        self.assertEqual(report['plan_state'], 'judged')

    def test_a_missing_plan_is_reported_as_not_checked(self):
        # CI has no plan file; "nothing was judged" must never read as "the plan
        # is fine"
        v = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        report, problems = v.doc_value_problems(code, v.cli_subcommands(code))
        self.assertEqual(report['plan_state'], 'not checked (no --plan)')
        self.assertEqual(report['plan_judged'], 0, report)
        self.assertEqual(problems, [], report)


class TestExternalWriteGate(unittest.TestCase):
    """IT-7 has to read the runner's own external-write report back.

    `tools/run.sh` writes `EXTERNAL` lines into commands.log since round 39 and
    nothing consumed them; and a command run outside the runner leaves no trace
    there at all (round 42's `/tmp` scratch file is absent from the log for that
    reason). So the gate also samples the watched prefixes -- attributing by
    *namespace*, because round 44 measured that our own toolchain creates 0 files
    under $HOME/.cache while another application's `fontconfig` cache appears
    there on any busy day.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools', 'verify_workspace.py')
        spec = spec_from_file_location('vw1', path,
                                       loader=SourceFileLoader('vw1', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _dirs(self):
        root = tmpdir()
        self.addCleanup(shutil.rmtree, root, ignore_errors=True)
        cache = join(root, 'cache')
        os.makedirs(join(cache, 'logs'), exist_ok=True)
        return root, cache

    def _log(self, cache, entries):
        with open(join(cache, 'logs', 'commands.log'), 'w',
                  encoding='utf-8') as f:
            f.write(''.join(entries))

    def test_an_external_write_named_in_the_log_is_found(self):
        v = self._tool()
        root, cache = self._dirs()
        self._log(cache, [
            '2026-11-03T10:00:00Z\t10:00:01Z\tintent\texit=0\tcwd=x\tcmd\n',
            '\t\tEXTERNAL\t%s/neoHGT2GTDB-临时文件/ok.json\n' % root,
            '\t\tEXTERNAL\t/home/someone/.cache/matplotlib/fontList.json\n',
        ])
        report, problems = v.check_external_writes(root, cache, 0,
                                                   quiet=(), noisy=())
        self.assertEqual(report['external_in_workspace'], 1, report)
        self.assertEqual(len(report['external_outside']), 1, report)
        self.assertTrue(any('matplotlib/fontList.json' in p for p in problems),
                        problems)

    def test_a_missing_or_empty_log_is_not_a_pass(self):
        v = self._tool()
        root, cache = self._dirs()
        _, problems = v.check_external_writes(root, cache, 0, quiet=(),
                                              noisy=())
        self.assertTrue(any('不存在' in p for p in problems), problems)
        self._log(cache, ['not a command line\n'])
        _, problems = v.check_external_writes(root, cache, 0, quiet=(),
                                              noisy=())
        self.assertTrue(any('一条命令都没有' in p for p in problems), problems)

    def test_our_own_cache_namespace_fails_the_gate(self):
        v = self._tool()
        root, cache = self._dirs()
        self._log(cache, ['2026-11-03T10:00:00Z\tcmd\n'])
        watch = join(root, 'watch')
        os.makedirs(join(watch, 'matplotlib'), exist_ok=True)
        with open(join(watch, 'matplotlib', 'fontList.cache'), 'w') as f:
            f.write('x')
        report, problems = v.check_external_writes(
            root, cache, time.time() - 60, quiet=(watch,), noisy=(),
            own_namespaces=('matplotlib',))
        self.assertEqual(report['quiet'][watch]['n_ours'], 1, report)
        self.assertTrue(any('matplotlib' in p for p in problems), problems)

    def _declared_watch(self, sub='pip/selfcheck', evidence=True, stated='now',
                        scope=True):
        """A watched prefix with one fresh file, plus the declaration's evidence.

        The synthetic root has to contain the code the exemption is measured over
        (round 53's freshness rule reads `newest_in_scope`): an evidence file over
        an empty scope is refused, so these fixtures lay down a `tools/` file too.
        """
        v = self._tool()
        root, cache = self._dirs()
        self._log(cache, ['2026-11-03T10:00:00Z\tcmd\n'])
        if scope:
            for rel in ('neoHGT2GTDB-项目代码/tools', 'neoHGT2GTDB-项目代码/neoHGT',
                        'neoHGT2GTDB-项目代码/neoHGT/tests'):
                os.makedirs(join(root, rel), exist_ok=True)
            with open(join(root, 'neoHGT2GTDB-项目代码/tools/v.py'), 'w') as f:
                f.write('#\n')
        scope_mtime, scope_path = v.newest_in_scope(root)
        watch = join(root, 'watch')
        os.makedirs(join(watch, sub), exist_ok=True)
        with open(join(watch, sub, 'afile'), 'w') as f:
            f.write('x')
        decl = v.DECLARED_FOREIGN_WRITERS[0]
        if evidence:
            ev = join(root, decl['evidence'])
            os.makedirs(join(root, os.path.dirname(decl['evidence'])), exist_ok=True)
            stated_mtime = (scope_mtime + 5 if stated == 'now'
                            else scope_mtime - 3600 if stated == 'before' else stated)
            with open(ev, 'w', encoding='utf-8') as f:
                json.dump({'our_stack_wrote': 0,
                           'commands_log_pip_invocations': 0,
                           'scope_measured': scope_path,
                           'newest_measured_code_mtime': stated_mtime}, f)
        return v, root, cache, watch

    def test_a_declared_foreign_writer_with_evidence_is_not_blamed(self):
        # round 52: pip's version self-check is written by whatever else on this
        # machine uses /opt/anaconda3; the exemption exists only while the probe
        # file says our own commands wrote nothing there.
        v, root, cache, watch = self._declared_watch()
        report, problems = v.check_external_writes(
            root, cache, time.time() - 60, quiet=(watch,), noisy=(),
            own_namespaces=('pip',))
        self.assertEqual(problems, [], (report, problems))
        entry = report['quiet'][watch]
        self.assertEqual(entry['n_declared_foreign'], 1, entry)
        self.assertEqual(entry['n_ours'], 0, entry)
        self.assertIn('IT7_pip_attribution', entry['recent_declared'][0], entry)

    def test_evidence_older_than_the_code_it_measures_is_refused(self):
        # round 53's freshness rule: a measurement cannot vouch for code written
        # after it -- which is not hypothetical, because round 52's own evidence
        # was three minutes older than verify_workspace.py when the rule was added
        v, root, cache, watch = self._declared_watch(stated='before')
        _, problems = v.check_external_writes(
            root, cache, time.time() - 60, quiet=(watch,), noisy=(),
            own_namespaces=('pip',))
        self.assertTrue(any('stale:' in p for p in problems), problems)
        self.assertTrue(any('r119_pip_attribution.py' in p for p in problems),
                        'a red gate has to say what to re-run: ' + str(problems))

    def test_touching_the_evidence_file_does_not_refresh_it(self):
        # the stated measurement wins over the file's own mtime: mtime is free to
        # anyone who runs `touch`, the number inside is part of what was measured
        v, root, cache, watch = self._declared_watch(stated='before')
        ev = join(root, v.DECLARED_FOREIGN_WRITERS[0]['evidence'])
        os.utime(ev, (time.time() + 10, time.time() + 10))
        _, problems = v.check_external_writes(
            root, cache, time.time() - 60, quiet=(watch,), noisy=(),
            own_namespaces=('pip',))
        self.assertTrue(any('stale:' in p for p in problems),
                        'a touch must not be able to renew an exemption')

    def test_an_exemption_over_an_empty_scope_is_refused(self):
        v, root, cache, watch = self._declared_watch(scope=False)
        _, problems = v.check_external_writes(
            root, cache, time.time() - 60, quiet=(watch,), noisy=(),
            own_namespaces=('pip',))
        self.assertTrue(any('no scope' in p for p in problems),
                        'a measurement over nothing is not a measurement: '
                        + str(problems))

    def test_the_pip_detector_reads_the_command_field_not_the_prose(self):
        # round 52's probe matched the words "pip redirection" inside its own
        # intent line and reported a pip invocation that never happened
        v = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        log = join(d, 'commands.log')
        with open(log, 'w', encoding='utf-8') as f:
            f.write('2026-11-04T07:00:00Z\t07:00:01Z\tround52 C pass: env.sh pip '
                    'redirection, 4 new tests\texit=0\tcwd=x\tbash gates.sh\n')
            f.write('\t\tOUT\tsome/file\n')
            f.write('2026-11-04T07:01:00Z\t07:01:02Z\tintent\texit=0\tcwd=x\t'
                    '/opt/anaconda3/bin/python -m pip install requests\n')
            f.write('2026-11-04T07:02:00Z\t07:02:01Z\tintent2\texit=0\tcwd=x\t'
                    'pip3 list\n')
        hits = v.pip_command_hits(log)
        self.assertEqual(len(hits), 2, hits)
        self.assertTrue(all('pip' in h for h in hits), hits)
        self.assertEqual(v.pip_command_hits(join(d, 'absent.log')), [])

    def test_the_freshness_scope_is_the_code_the_gates_execute(self):
        v = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        code = join(d, 'neoHGT2GTDB-项目代码')
        os.makedirs(join(code, 'neoHGT', 'tests'), exist_ok=True)
        os.makedirs(join(code, 'tools'), exist_ok=True)
        os.makedirs(join(code, 'docs'), exist_ok=True)
        for rel, when in (('neoHGT/a.py', 100), ('neoHGT/tests/b.py', 300),
                          ('tools/c.py', 200), ('tools/run.sh', 400),
                          ('docs/d.md', 900)):          # docs are not executed code
            p = join(code, rel)
            with open(p, 'w') as f:
                f.write('x\n')
            os.utime(p, (time.time() - (1000 - when),) * 2)
        mtime, path = v.newest_in_scope(d)
        self.assertEqual(path, 'neoHGT2GTDB-项目代码/tools/run.sh', path)
        _, outside = v.newest_in_scope(join(d, 'empty'))
        self.assertIsNone(outside, 'a missing scope has to be reported, not agreed with')

    def test_a_stated_measurement_survives_an_old_file_mtime(self):
        # the evidence's own number is what counts; a copy that resets mtime must
        # not invalidate a measurement that is still about the code it saw
        v, root, cache, watch = self._declared_watch(stated='now')
        ev = join(root, v.DECLARED_FOREIGN_WRITERS[0]['evidence'])
        old = time.time() - 86400
        os.utime(ev, (old, old))
        _, problems = v.check_external_writes(
            root, cache, time.time() - 60, quiet=(watch,), noisy=(),
            own_namespaces=('pip',))
        self.assertEqual(problems, [], (root, problems))

    def test_the_declaration_without_its_evidence_is_a_problem(self):
        v, root, cache, watch = self._declared_watch(evidence=False)
        report, problems = v.check_external_writes(
            root, cache, time.time() - 60, quiet=(watch,), noisy=(),
            own_namespaces=('pip',))
        self.assertTrue(any('豁免必须有证据' in p for p in problems),
                        (report, problems))

    def test_a_declaration_whose_evidence_says_otherwise_is_a_problem(self):
        v, root, cache, watch = self._declared_watch()
        ev = v.DECLARED_FOREIGN_WRITERS[0]['evidence']
        with open(join(root, ev), 'w', encoding='utf-8') as f:
            json.dump({'our_stack_wrote': 3, 'commands_log_pip_invocations': 0}, f)
        _, problems = v.check_external_writes(
            root, cache, time.time() - 60, quiet=(watch,), noisy=(),
            own_namespaces=('pip',))
        self.assertTrue(any('豁免必须有证据' in p for p in problems), problems)

    def test_only_the_declared_sub_path_is_excused_not_the_namespace(self):
        # excusing 'pip' wholesale would turn the rule off for the one tool that
        # can install packages, which is the opposite of what an exemption is for
        v, root, cache, watch = self._declared_watch(sub='pip/http')
        _, problems = v.check_external_writes(
            root, cache, time.time() - 60, quiet=(watch,), noisy=(),
            own_namespaces=('pip',))
        self.assertTrue(any('缓存' in p and 'pip/http' in p for p in problems),
                        problems)

    def test_another_app_s_namespace_is_reported_not_judged(self):
        v = self._tool()
        root, cache = self._dirs()
        self._log(cache, ['2026-11-03T10:00:00Z\tcmd\n'])
        watch = join(root, 'watch')
        os.makedirs(join(watch, 'fontconfig'), exist_ok=True)
        with open(join(watch, 'fontconfig', 'abc-le64.cache-9'), 'w') as f:
            f.write('x')
        report, problems = v.check_external_writes(
            root, cache, time.time() - 60, quiet=(watch,), noisy=(),
            own_namespaces=('matplotlib',))
        self.assertEqual(report['quiet'][watch]['n_other'], 1, report)
        self.assertEqual(problems, [], report)

    def test_not_judged_prefixes_are_named_without_being_walked(self):
        # round 44 walked /tmp (200k+ entries) every run just to print "not
        # attributable"; the caveat is the point, the count was noise
        v = self._tool()
        root, cache = self._dirs()
        self._log(cache, ['2026-11-03T10:00:00Z\tcmd\n'])
        report, problems = v.check_external_writes(
            root, cache, time.time(), quiet=(), noisy=('/definitely-absent',))
        self.assertEqual(problems, [], report)
        entry = report['noisy']['/definitely-absent']
        self.assertFalse(entry['judged'], entry)
        self.assertNotIn('entries', entry, entry)
        self.assertIn('not attributable', entry['note'], entry)

    def test_the_real_workspace_passes_on_its_own_window(self):
        # this machine is not idle: the shipped tree must still come out clean
        # when judged over a window that is only the current round
        v = self._tool()
        root = os.environ.get('NEOHGT_ROOT')
        cache = os.environ.get('NEOHGT_CACHE')
        if not root or not cache:
            self.skipTest('NEOHGT_ROOT/NEOHGT_CACHE not set')
        report, problems = v.check_external_writes(root, cache, time.time())
        if not report.get('log_exists'):
            # the command log is a workspace-history artifact: a session that
            # never ran the pipeline has no window to judge (the fresh-cache
            # branch used to fail with "这条核对什么都没看")
            self.skipTest('no commands.log in this cache: nothing to judge')
        self.assertEqual(problems, [], report)
        self.assertGreater(report['commands'], 100, report)


class TestBandwidthValueTypes(unittest.TestCase):
    """A value is valid or refused -- not valid on one route and refused on another.

    Measured in round 45 (`缓存测试/reports/R45_bandwidth_routes.txt`, from
    `临时文件/r54_bandwidth_routes.sh`): a copied ./config.yml carrying
    `bandwidth: 1` -- YAML's integer, and the maximum the help text allows --
    stopped the run with `ValueError: Invalid bandwidth: 1.`, while the same
    command with `--bandwidth 1` predicted 222 HGTs. The reason is two partial
    normalisations: `perform_kde()` asks for `isinstance(bw, float)`, and the
    GTDB wrapper only understood the *string* form the parser hands over.
    """

    def _run(self, bandwidth):
        from types import SimpleNamespace
        from unittest import mock
        from neoHGT import analyze as A
        me = SimpleNamespace(source='gtdb', gtdb_root='/tmp-ish',
                             taxdump={'1': {'parent': '1', 'rank': 'no rank',
                                            'name': 'root'}},
                             gtdb_labels={'1': 'g__x'}, bandwidth=bandwidth,
                             df=pd.DataFrame({'close': [1.0, 2.0],
                                              'distal': [3.0, 4.0],
                                              'self': [0.0, 0.0]}))
        with mock.patch.object(A, 'apply_gtdb_grouping', lambda s: None), \
                mock.patch.object(A, 'check_score_groups', lambda s: None), \
                mock.patch.object(A, '_orig_predict_hgt', lambda s: 0):
            A.predict_hgt_gtdb(me)
        return me.bandwidth

    def test_the_documented_maximum_works_from_both_routes(self):
        for value, want in ((1, 1.0), ('1', 1.0), (1.0, 1.0), ('0.5', 0.5)):
            got = self._run(value)
            self.assertEqual(got, want, value)
            # the type is the point: perform_kde() branches on isinstance(bw,
            # float), and an int that merely compares equal to 1.0 still misses
            # that branch -- so an assertion on the value alone cannot see the bug.
            self.assertIsInstance(got, float, value)

    def test_keywords_are_left_to_the_branches_that_read_them(self):
        for kw in ('auto', 'grid', 'silverman'):
            self.assertEqual(self._run(kw), kw, kw)

    def test_an_out_of_range_number_is_refused_with_its_domain_and_origin(self):
        for value, origin in ((2, 'predict.bandwidth'), ('0.05', '--bandwidth 0.05')):
            with self.assertRaises(ValueError) as ctx:
                self._run(value)
            msg = str(ctx.exception)
            self.assertIn(origin, msg, msg)
            self.assertIn('between 0.1 and 1.0', msg, msg)

    def test_a_boolean_is_not_a_bandwidth(self):
        # python counts True as an int; it must not be floated into 1.0
        self.assertIs(self._run(True), True)


class TestFetchCountCoercion(unittest.TestCase):
    """`--fetch-queries 50` used to reach a slice as a string.

    Measured in round 46 (`缓存测试/reports/R46_fetch_route_types.json`, before
    `缓存测试/reports/R46_fetch_after_fix.json`): in the shipped `arguments` table
    `--fetch-delay` and `--fetch-timeout` declare `{'type': int}` while their two
    neighbours `--fetch-queries` and `--fetch-retries` declare nothing, so argparse
    hands those two over as strings -- and the consumers count with them
    (`cids[i:i + cq]`, `range(self.fetch_retries)`). The CLI route died with
    `TypeError: 'str' object cannot be interpreted as an integer` before it ever
    reached the network; the same value in config.yml (an int) ran fine.
    """

    def _obj(self, **kw):
        from neoHGT import search as S
        obj = S.Search.__new__(S.Search)   # no __init__: no filesystem, no db
        obj.fetch_queries = 50
        obj.fetch_retries = 3
        obj.fetch_delay = 0
        obj.fetch_timeout = 1
        obj.fetch_server = 'http://127.0.0.1:1/x'
        for k, v in kw.items():
            setattr(obj, k, v)
        S.coerce_fetch_counts(obj)
        return obj

    def test_both_routes_end_up_with_the_same_object(self):
        for attr, cli, cfg in (('fetch_queries', '50', 50),
                               ('fetch_retries', '3', 3)):
            a = self._obj(**{attr: cli})
            b = self._obj(**{attr: cfg})
            self.assertEqual(getattr(a, attr), getattr(b, attr))
            self.assertIsInstance(getattr(a, attr), int, attr)

    def test_a_non_numeric_count_is_refused_naming_the_option(self):
        from neoHGT import search as S
        for attr, flag in (('fetch_queries', '--fetch-queries'),
                           ('fetch_retries', '--fetch-retries')):
            obj = S.Search.__new__(S.Search)
            setattr(obj, attr, 'many')
            with self.assertRaises(ValueError) as ctx:
                S.coerce_fetch_counts(obj)
            self.assertIn(flag, str(ctx.exception), str(ctx.exception))

    def test_a_count_below_one_is_refused_rather_than_silently_nooping(self):
        # `range(0)` and `range(-1)` are both empty loops: a fetch that was asked
        # for and did nothing is worse than a refusal.
        from neoHGT import search as S
        for attr in ('fetch_queries', 'fetch_retries'):
            obj = S.Search.__new__(S.Search)
            setattr(obj, attr, '0')
            with self.assertRaises(ValueError) as ctx:
                S.coerce_fetch_counts(obj)
            self.assertIn('at least 1', str(ctx.exception))

    def test_a_boolean_is_not_a_count(self):
        from neoHGT import search as S
        obj = S.Search.__new__(S.Search)
        obj.fetch_retries = True            # python counts True as int 1
        with self.assertRaises(ValueError):
            S.coerce_fetch_counts(obj)

    def test_an_unset_count_is_left_unset(self):
        # `None` means neither the flag nor the config supplied it; the upstream
        # default has already been applied by then, so inventing 1 here would
        # overwrite a value the program chose on purpose.
        obj = self._obj(fetch_queries=None)
        self.assertIsNone(obj.fetch_queries)

    def test_the_consumers_no_longer_die_on_the_cli_route(self):
        # The fix has to be seen by the code that had the bug, not just by the
        # helper: with urlopen replaced by a recorder, the CLI route now reaches
        # the request instead of raising TypeError on the way there.
        from unittest import mock
        from neoHGT import search as S
        for attr, value, method, callargs in (
                ('fetch_queries', '50', 'remote_fetches',
                 (['P1', 'P2'], 'query={}')),
                ('fetch_retries', '3', 'remote_fetch', ('x=1',))):
            obj = self._obj(**{attr: value})
            called = []

            def fake_urlopen(*a, **k):
                called.append(a)
                raise OSError('stub')

            with mock.patch.object(S, 'urlopen', fake_urlopen):
                with self.assertRaises(OSError) as ctx:
                    getattr(obj, method)(*callargs)
            self.assertEqual(str(ctx.exception), 'stub', method)
            self.assertTrue(called, f'{method} never reached the request: {obj}')


    def test_the_wrapper_runs_the_coercion_on_every_route(self):
        # The helper being correct is not the fix: `set_params_gtdb` has to call
        # it, and it has to do so for NCBI runs too (the wrapper replaces
        # `args_wf` unconditionally). Everything the wrapper does *around* the
        # coercion is stubbed out, so this exercises the call site and nothing
        # else -- a mutant that deletes the line fails here and nowhere else.
        from types import SimpleNamespace
        from unittest import mock
        from neoHGT import search as S
        obj = S.Search.__new__(S.Search)
        obj.fetch_queries = '50'
        obj.fetch_retries = '3'
        args = SimpleNamespace(source='auto', gtdb_stem2taxid=None,
                               keep_unresolved=None, require_staxids=None,
                               fetch_enable='no')
        with mock.patch.object(S, '_orig_set_params', lambda self, a: None), \
                mock.patch.object(S, 'gtdb_db_dir', lambda self: None), \
                mock.patch.object(S, 'guard_precomp_multisample',
                                  lambda self: None):
            S.set_params_gtdb(obj, args)
        self.assertEqual(obj.fetch_queries, 50)
        self.assertEqual(obj.fetch_retries, 3)


class TestRouteEquivalenceRule(unittest.TestCase):
    """Register rule 19 has to be able to *see* the thing it claims to check.

    Every gate in this project has eventually been caught reading nothing: a walk
    that missed a call form, a grep whose pattern went stale, a ledger nobody
    re-derives. So these tests plant the violations the rule exists to catch and
    require it to name them -- and require it to refuse when it can't look.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt19', path,
                                       loader=SourceFileLoader('vpt19', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _tree(self, keys, config_lines, table, reads):
        d = tmpdir()
        self.addCleanup(__import__('shutil').rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'neoHGT'))
        os.makedirs(join(d, 'docs'))
        with open(join(d, 'neoHGT', 'analyze.py'), 'w',
                  encoding='utf-8') as fh:
            fh.write(f'arguments = ["g", {table}]\n\n\n'
                     f'def set_parameters(self):\n'
                     f'    {reads}\n')
        for other in ('search', 'database'):
            with open(join(d, 'neoHGT', f'{other}.py'), 'w',
                      encoding='utf-8') as fh:
                fh.write('arguments = []\n')
        with open(join(d, 'neoHGT', 'config.yml'), 'w',
                  encoding='utf-8') as fh:
            fh.write(config_lines)
        ledger = {'keys': keys, 'unresolved_call_sites': []}
        path = join(d, 'docs', '路由等价性清单.json')
        with open(path, 'w', encoding='utf-8') as fh:
            json.dump(ledger, fh, ensure_ascii=False)
        return d, path

    def test_an_unexplained_divergence_is_named(self):
        v = self._tool()
        keys = {'predict.someint': {'attr': 'some_int', 'converter': None,
                                    'cli_type': 'int', 'yaml_type': 'str',
                                    'in_config': True, 'divergent': True}}
        d, ledger = self._tree(
            keys, 'predict:\n  someint: "12"\n',
            "['--some-int', 'some count', {'type': int}]",
            "get_config(self, 'some_int', 'predict.someint')")
        report, problems = v.route_equivalence_problems(d, ledger)
        self.assertTrue(any('predict.someint' in p for p in problems), problems)
        self.assertEqual(report['divergent_undeclared'], ['predict.someint'])

    def test_a_stale_exemption_is_also_a_problem(self):
        # An exemption outlives the bug it defended, and then it is quietly
        # covering nothing -- worse than no list, because it reads like care.
        v = self._tool()
        key = 'predict.weighted'
        self.assertIn(key, v.ROUTE_NEUTRALIZERS)
        keys = {key: {'attr': 'weighted', 'converter': None, 'cli_type': 'bool',
                      'yaml_type': 'bool', 'in_config': True,
                      'divergent': False}}
        d, ledger = self._tree(
            keys, 'predict:\n  weighted: true\n',
            "['--weighted', 'use weighted scores', {'action': 'store_true'}]",
            "get_config(self, 'weighted', 'predict.weighted')")
        report, problems = v.route_equivalence_problems(d, ledger)
        self.assertTrue(any('豁免' in p and key in p for p in problems), problems)

    def test_a_missing_ledger_refuses_instead_of_passing(self):
        v = self._tool()
        d = tmpdir()
        report, problems = v.route_equivalence_problems(
            d, join(d, 'docs', '路由等价性清单.json'))
        self.assertEqual(report['state'], 'ledger missing')
        self.assertTrue(problems, 'a missing ledger must not read as a pass')

    def test_a_blind_walk_refuses(self):
        v = self._tool()
        keys = {'a.b': {'attr': 'b', 'converter': None, 'cli_type': 'int',
                        'yaml_type': 'str', 'in_config': True,
                        'divergent': True}}
        d, ledger = self._tree(keys, 'a:\n  b: "1"\n',
                               "['--b', 'x', {'type': int}]",
                               "get_config(self, 'b', 'a.b')")
        report, problems = v.route_equivalence_problems(d, ledger)
        self.assertTrue(any('绿灯不算数' in p for p in problems),
                        f'{problems}: a walk that found {report["reads"]} reads '
                        'should not present itself as a clean audit')

    def test_the_shipped_tree_passes_with_every_divergence_explained(self):
        v = self._tool()
        code = os.environ.get('NEOHGT_CODE')
        if not code:
            self.skipTest('NEOHGT_CODE is not set (bare checkout)')
        ledger = join(code, 'docs', '路由等价性清单.json')
        self.assertTrue(isfile(ledger), f'{ledger} must be committed')
        report, problems = v.route_equivalence_problems(code, ledger)
        self.assertEqual(problems, [], report)
        self.assertEqual(len(report['divergent_explained']),
                         report['divergent'],
                         'every divergent key must be explained -- and no more')
        self.assertEqual(report['divergent_undeclared'], [])
        self.assertEqual(report['ledger_unlisted_keys'], [])


    def test_a_declared_slot_expands_and_an_undeclared_one_is_reported(self):
        # Both halves matter: expanding `f'{self.method}.{key}'` over the choices
        # the option declares is reading the code, but if the choices stop being
        # declared the call site has to surface as unresolved -- silently dropping
        # it is how rule 19 came to audit 58 keys while the census saw 67.
        v = self._tool()
        d = tmpdir()
        self.addCleanup(__import__('shutil').rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'neoHGT'))
        src = ('arguments = ["g",\n'
               "    ['-m|--method', 'search method', "
               "{'choices': ['auto', 'diamond']}],\n"
               "    ['--queries', 'q', {'type': int}],\n"
               "]\n\n\n"
               'def args_wf(self, args):\n'
               "    for key in 'queries',:\n"
               "        get_config(self, key, f'{self.method}.{key}')\n")
        for mod in ('search', 'database', 'analyze'):
            with open(join(d, 'neoHGT', f'{mod}.py'), 'w',
                      encoding='utf-8') as fh:
                fh.write(src if mod == 'search' else 'arguments = []\n')
        rows, unres = v.route_reads(d)
        entries = {e for (_m, _a, e, _c, _b) in rows}
        self.assertIn('auto.queries', entries)
        self.assertIn('diamond.queries', entries)
        self.assertEqual(unres, [],
                         'a slot whose choices are declared must not be refused')
        # now the same file with the choices removed from the option
        with open(join(d, 'neoHGT', 'search.py'), 'w', encoding='utf-8') as fh:
            fh.write(src.replace("{'choices': ['auto', 'diamond']}", "{}"))
        rows, unresolved = v.route_reads(d)
        self.assertEqual({e for (_m, _a, e, _c, _b) in rows}, set(),
                         'without declared choices there is nothing to expand to')
        self.assertTrue(unresolved, 'a refused call site must be listed')
        # the line number is whatever the snippet happens to be; the claim under
        # test is that the site is *reported*, so assert that and not a coordinate
        self.assertIn('search.py:', unresolved[0], unresolved)
        self.assertIn('not resolvable', unresolved[0], unresolved)

    def test_the_shipped_ledger_refuses_nothing_and_marks_its_expansions(self):
        v = self._tool()
        code = os.environ.get('NEOHGT_CODE')
        if not code:
            self.skipTest('NEOHGT_CODE is not set (bare checkout)')
        ledger = join(code, 'docs', '路由等价性清单.json')
        with open(ledger, encoding='utf-8') as fh:
            data = json.load(fh)
        self.assertEqual(data['unresolved_call_sites'], [],
                         'every config read is in scope today; a reappearing '
                         'refusal means the walk narrowed again')
        expanded = [k for k, val in data['keys'].items()
                    if val.get('expansion_basis')]
        self.assertTrue(expanded,
                        'the method-family keys must say why they are in scope')
        self.assertTrue(all(val['in_config'] is False or val['yaml_type'] != 'str'
                            for val in data['keys'].values()
                            if val.get('divergent')),
                        'divergent means the two routes disagree, nothing more')


class TestRouteClaimVerification(unittest.TestCase):
    """Rule 19b: an exemption must name a normaliser that exists, covers the key,
    and runs after the read -- not merely a sentence that sounds like one. The
    places that sentence cites have to exist too.

    Rule 19 asks whether a divergence is *explained*; the explanation was prose, so
    a stale or wrong one still read as care. These build minimal trees rather than
    mutating the shipped ones: the point is which feature of the code the check
    reacts to, and a synthetic tree answers that in both directions (a claim that
    should pass does, so a failure means the claim and not the checker).
    """

    HEAD = ('arguments = ["g",\n'
            "    ['-m|--method', 'search method', "
            "{'choices': ['auto', 'diamond']}],\n"
            "    ['--capital', 'capitalised', {}],\n"
            "]\n\n\n")

    def _code(self, body):
        d = tmpdir()
        self.addCleanup(__import__('shutil').rmtree, d, ignore_errors=True)
        for mod in ('search', 'database', 'analyze'):
            os.makedirs(join(d, 'neoHGT'), exist_ok=True)
            with open(join(d, 'neoHGT', f'{mod}.py'), 'w',
                      encoding='utf-8') as fh:
                fh.write(self.HEAD + body if mod == 'search'
                         else 'arguments = []\n')
        os.makedirs(join(d, 'neoHGT', 'tests'), exist_ok=True)
        with open(join(d, 'neoHGT', 'tests', 'test_gtdb_pipeline.py'), 'w',
                  encoding='utf-8') as fh:
            fh.write('')
        return d

    def _claims(self, d, key):
        v = self._tool()
        return v.route_claim_problems(d, {key})

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt19b', path,
                                       loader=SourceFileLoader('vpt19b', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_a_sound_claim_passes_and_says_which_lines_it_read(self):
        v = self._tool()
        key = next(k for k, c in v.ROUTE_CLAIMS.items() if c[0] == 'search')
        attr = v.ROUTE_CLAIMS[key][2]
        body = ('def args_wf(self, args):\n'
                "    for key in 'capital',:\n"
                f"        get_config(self, f'{attr[:4]}{{key}}', "
                "    f'taxonomy.{key}')\n"
                "    for key in 'tax_capital', 'tax_latin':\n"
                "        setattr(self, key, arg2bool(getattr(self, key, None)))\n")
        d = self._code(body)
        report, problems = self._claims(d, key)
        self.assertEqual(problems, [], report)
        got = report['verdicts'][key]
        self.assertEqual(got['function'], 'args_wf')
        self.assertGreater(got['line'], got['read_line'])

    def test_a_claim_naming_nothing_real_is_a_problem(self):
        v = self._tool()
        key = 'taxonomy.capital'
        d = self._code('def args_wf(self, args):\n'
                       "    for key in 'capital',:\n"
                       "        get_config(self, f'tax_{key}', f'taxonomy.{key}')\n")
        report, problems = self._claims(d, key)
        self.assertTrue(any('找不到' in p and key in p for p in problems),
                        f'{problems}: a reason naming nothing must not pass')
        self.assertNotIn(key, report['verdicts'])

    def test_normalising_before_the_read_is_a_problem(self):
        # Same call, same key list, one line earlier: the value is converted, then
        # overwritten by the config read. Line order is the only thing that tells
        # these apart, so the check has to look at it.
        v = self._tool()
        key = 'taxonomy.capital'
        d = self._code('def args_wf(self, args):\n'
                       "    for key in 'tax_capital',:\n"
                       "        setattr(self, key, arg2bool(getattr(self, key)))\n"
                       "    for key in 'capital',:\n"
                       "        get_config(self, f'tax_{key}', f'taxonomy.{key}')\n")
        report, problems = self._claims(d, key)
        self.assertTrue(any('之前' in p for p in problems), problems)

    def test_an_unorderable_claim_needs_a_test_that_exists(self):
        # Loop-form reads can hide which iteration wrote which attribute. When the
        # line order cannot be decided, the claim must point at a test that pins the
        # behaviour by execution -- and that test must exist.
        v = self._tool()
        key = 'fetch.queries'
        d = self._code('def args_wf(self, args):\n'
                       "    for key in 'queries',:\n"
                       "        get_config(self, f'fetch_{key}', f'fetch.{key}')\n"
                       'def set_params_gtdb(self, args):\n'
                       '    coerce_fetch_counts(self)\n'
                       'def coerce_fetch_counts(self):\n'
                       "    for attr in FETCH_COUNTS:\n"
                       "        setattr(self, attr, int(getattr(self, attr)))\n")
        with open(join(d, 'neoHGT', 'search.py'), 'a', encoding='utf-8') as fh:
            fh.write('\nFETCH_COUNTS = {"fetch_queries": "--fetch-queries"}\n')
        report, problems = self._claims(d, key)
        self.assertTrue(any('不在测试里' in p for p in problems),
                        f'{problems}: no proof test exists in this synthetic tree')

    def test_the_shipped_claims_are_all_checked(self):
        v = self._tool()
        code = os.environ.get('NEOHGT_CODE')
        if not code:
            self.skipTest('NEOHGT_CODE is not set (bare checkout)')
        ledger = join(code, 'docs', '路由等价性清单.json')
        with open(ledger, encoding='utf-8') as fh:
            data = json.load(fh)
        explained = {k for k, val in data['keys'].items()
                     if val.get('divergent') and not val.get('converter')}
        report, problems = v.route_claim_problems(code, explained)
        self.assertEqual(problems, [], report)
        self.assertEqual(report['claims_checked'], len(explained),
                         'every written reason must be checked, none skipped')
        self.assertTrue(explained, 'no exemptions left to verify means the ledger '
                                  'changed shape and this test stopped testing')

    def test_a_reason_naming_a_function_that_is_not_there_is_a_problem(self):
        # The reason is prose, and prose rots quietly. Three shipped reasons said
        # "search.py set_parameters" while that module's argument handler is
        # `args_wf`; the normaliser lookup cannot see that, because it follows the
        # call rather than the sentence. So the sentence is read too: a location a
        # reason names must be defined, and it must be the *only* thing that fails
        # here -- the normaliser is found, ordered and looped exactly as before.
        v = self._tool()
        key = 'taxonomy.capital'          # its reason cites `search.py args_wf`
        d = self._code('def args_wf_renamed(self, args):\n'
                       "    for key in 'capital',:\n"
                       "        get_config(self, f'tax_{key}', f'taxonomy.{key}')\n"
                       "    for key in 'tax_capital',:\n"
                       "        setattr(self, key, arg2bool(getattr(self, key)))\n")
        report, problems = self._claims(d, key)
        self.assertIn(key, report['verdicts'],
                      'the normaliser itself was found, so this failure is the '
                      'citation and not a missing converter')
        self.assertTrue(any('没有定义' in p and key in p for p in problems), problems)
        self.assertTrue(any('search.py args_wf' in p for p in problems), problems)

    def test_the_shipped_reasons_cite_only_things_that_exist(self):
        v = self._tool()
        code = os.environ.get('NEOHGT_CODE')
        if not code:
            self.skipTest('NEOHGT_CODE is not set (bare checkout)')
        with open(join(code, 'docs', '路由等价性清单.json'), encoding='utf-8') as fh:
            data = json.load(fh)
        explained = {k for k, val in data['keys'].items()
                     if val.get('divergent') and not val.get('converter')}
        report, problems = v.route_claim_problems(code, explained)
        self.assertEqual(problems, [], report)
        # `>= len(explained)`, not a fixed number: a reason with no location in it
        # at all would leave the citation check reading nothing while still
        # reporting no problems.
        self.assertGreaterEqual(report['citations_seen'], len(explained),
                                'some reason stopped naming a place, so there is '
                                'nothing left for this check to verify')
        self.assertEqual(report['citations_unjudged'], [],
                         'a citation this checkout cannot look up is not a pass')
        self.assertEqual(report['citations_missing_modules'], [], report)


class TestDefaultPromises(unittest.TestCase):
    """Rule 20: a documented default must equal the value the code delivers.

    Round 50 found the first real instance by executing rather than reading: a
    default GTDB run wrote the synthetic taxId `105000000` into hgts/*.txt while
    the deliverables table promised prefixed labels, because `--donor-label` had
    no table default, no config key, and a consumer that only acted on the literal
    'prefixed'. These build one synthetic tree carrying every shape the rule
    models, then mutate one thing at a time: the claims themselves are shipped
    data, so a rule that only ever sees the real tree proves it can agree, not that
    it can disagree.
    """

    PY = {
        'database.py': (
            "def set_source_parameters(self, args):\n"
            "    self.seq_source = (getattr(args, 'seq_source', None)\n"
            "                       or gcfg.get('seq_source') or 'local')\n"
            "    self.id_base = (getattr(args, 'gtdb_taxid_base', None)\n"
            "                    or gcfg.get('taxid_base')\n"
            "                    or _gtdb.DEFAULT_ID_BASE)\n"
            "    self.allow_network = arg2bool(\n"
            "        getattr(args, 'allow_network', None)) is True\n"),
        'search.py': (
            "def set_params_gtdb(self, args):\n"
            "    prof = {}\n"
            "    if getattr(args, 'tax_capital', None) is None:\n"
            "        self.keep_placeholder_names = bool(prof.get(\n"
            "            'keep_placeholder', True))\n"
            "    if getattr(args, 'require_staxids', None) is None:\n"
            "        self.require_staxids = is_gtdb_mode(self)\n"
            "    else:\n"
            "        self.require_staxids = arg2bool(args.require_staxids)\n"),
        'analyze.py': (
            "def write_hgt_list_gtdb(self, sample):\n"
            "    label = getattr(self, 'donor_label', None)\n"
            "    if label is None:\n"
            "        label = self.donor_label = 'prefixed'\n"
            "\n"
            "def apply_gtdb_grouping(self):\n"
            "    for attr, name in (('close_size', 'closesize'),):\n"
            "        setattr(self, attr, 200)\n"),
        'gtdb.py': ("DEFAULT_ID_BASE = 100_000_000\n"
                    "arguments = ['g',\n"
                    "             ['--gtdb-taxid-base', 'start of the synthetic "
                    "taxId range',\n"
                    "              {'type': int, 'default': DEFAULT_ID_BASE}]]\n"),
        'util.py': "def arg2bool(val):\n    return bool(val)\n",
    }
    DOCS = {
        'docs/决策落点对照.md': "默认路线：`--seq-source local` 为默认值。\n",
        'docs/方案修订记录核对.md': "`--seq-source local` 为默认，远端补数。\n",
        'docs/GTDB-数据库.md': (
            "`--gtdb-taxid-base` 默认 `100000000`。\n"
            "`--keep-placeholder-names`（默认 `yes`）豁免编码标签。\n"
            "`--require-staxids` 默认立刻报错（RK5）。\n"
            "`--close-size` 默认按域取档（细菌 200、古菌 50）。\n"
            "`--seq-source` 默认 `local`。\n"
            "`--donor-label` 默认 `prefixed`（带前缀标签）。\n"
            "`--allow-network` 默认 no（联网开关）。\n"),
        'docs/计划任务对照表.md': "`--donor-label` 三档输出，默认前缀标签 prefixed。\n",
        'docs/验收标准对照.md': "`--allow-network` 默认关。\n",
        'docs/风险应对登记.md': ("`--keep-placeholder-names`（默认 `yes`）；"
                                 "`--allow-network` 默认关。\n"),
    }
    CONFIG = ("grouping:\n"
              "  profiles:\n"
              "    gtdb:\n"
              "      closesize: 200\n"
              "    gtdb_archaea:\n"
              "      closesize: 50\n"
              # the real config also carries taxid_base; without it the chain's
              # `_gtdb.DEFAULT_ID_BASE` operand is an alias this fixture does not
              # resolve, and the promise would read as unjudgeable rather than kept
              "database:\n"
              "  gtdb:\n"
              "    taxid_base: 100000000\n")

    def _tool(self):
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vptr20', path,
                                       loader=SourceFileLoader('vptr20', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _tree(self, edits=None):
        """The fixture tree, with `edits` = [(path, old, new)] applied per file."""
        d = tmpdir()
        self.addCleanup(__import__('shutil').rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'neoHGT'))
        os.makedirs(join(d, 'docs'))
        for name, src in self.PY.items():
            with open(join(d, 'neoHGT', name), 'w', encoding='utf-8') as fh:
                fh.write('arguments = []\n\n\n' + src)
        for rel, src in self.DOCS.items():
            with open(join(d, rel), 'w', encoding='utf-8') as fh:
                fh.write(src)
        with open(join(d, 'neoHGT', 'config.yml'), 'w', encoding='utf-8') as fh:
            fh.write(self.CONFIG)
        for rel, old, new in (edits or []):
            full = join(d, rel)
            text = open(full, encoding='utf-8').read()
            assert text.count(old) == 1, (rel, old, text.count(old))
            with open(full, 'w', encoding='utf-8') as fh:
                fh.write(text.replace(old, new, 1))
        return d

    def _judge(self, d):
        v = self._tool()
        return v.default_promise_problems(d)

    def test_every_shipped_promise_agrees_and_none_is_unjudgeable(self):
        v = self._tool()
        report, problems = v.default_promise_problems(
            join(HERE, os.pardir, os.pardir))
        self.assertEqual(problems, [], report)
        self.assertGreaterEqual(report['promises'], 7, report)
        self.assertEqual(len(report['agree']), report['promises'],
                         'a promise that could not be judged must not be counted '
                         'as agreement')
        self.assertEqual(report['unjudgeable'], [], report)
        self.assertEqual(report['config_state'], 'read')

    def test_a_delivery_that_moved_is_reported_as_disagreement(self):
        # the docs still promise 'local'; the code now delivers 'remote'
        report, problems = self._judge(self._tree([
            ('neoHGT/database.py', "or 'local')", "or 'remote')")]))
        self.assertTrue(any('--seq-source' in p and 'remote' in p
                            for p in problems), problems)
        self.assertTrue(any('seq_source' in a for a in report['disagree']), report)
        self.assertFalse(any(a.startswith('seq_source=') for a in report['agree']),
                         report['agree'])
        self.assertTrue(any('taxid_base' in a for a in report['agree']),
                        'the other promises must still pass: one changed delivery '
                        'must not silence the whole rule')

    def test_a_promise_deleted_from_the_docs_is_reported(self):
        report, problems = self._judge(self._tree([
            ('docs/GTDB-数据库.md',
             '`--donor-label` 默认 `prefixed`（带前缀标签）。', '`--donor-label` 见上。')]))
        self.assertTrue(any('--donor-label' in p and '默认' in p
                            for p in problems), problems)

    def test_a_table_default_written_as_a_constant_is_resolved(self):
        """`default=DEFAULT_ID_BASE` is a delivery, not a gap.

        Rule 20 shipped without this shape -- `gtdb --gtdb-taxid-base` writes the
        constant rather than a literal, so the claim was unreadable on the day the
        rule landed. Moving the constant must move the verdict.
        """
        report, problems = self._judge(self._tree())
        self.assertEqual(problems, [], report)
        self.assertEqual(report['unjudgeable'], [], report)
        self.assertTrue(any(a.startswith('gtdb_taxid_base=100000000')
                            for a in report['agree']), report['agree'])
        bad, _p = self._judge(self._tree([
            ('neoHGT/gtdb.py', 'DEFAULT_ID_BASE = 100_000_000',
             'DEFAULT_ID_BASE = 900_000_000')]))
        self.assertIn('gtdb_taxid_base', bad['disagree'], bad)

    def test_an_establishing_function_that_is_gone_is_not_a_pass(self):
        _report, problems = self._judge(self._tree([
            ('neoHGT/analyze.py', 'def write_hgt_list_gtdb(self, sample):',
             'def write_hgt_list_renamed(self, sample):')]))
        self.assertTrue(any('donor_label' in p and '不在了' in p
                            for p in problems), problems)

    def test_a_profile_default_without_config_is_unjudgeable_not_agreed(self):
        # remove the per-domain values: the rule must say "I could not look"
        report, _problems = self._judge(self._tree([
            ('neoHGT/config.yml', self.CONFIG, 'grouping:\n  closesize: 10\n')]))
        self.assertTrue(any('close_size' in u for u in report['unjudgeable']),
                        report)

    def test_a_tree_where_nothing_can_be_judged_says_so(self):
        d = self._tree([('neoHGT/database.py',
                         'def set_source_parameters(self, args):',
                         'def other_name(self, args):'),
                        ('neoHGT/search.py', 'def set_params_gtdb(self, args):',
                         'def other_search(self, args):'),
                        ('neoHGT/analyze.py', 'def write_hgt_list_gtdb(self, sample):',
                         'def other_a(self, sample):'),
                        ('neoHGT/analyze.py', 'def apply_gtdb_grouping(self):',
                         'def other_b(self):'),
                        ('neoHGT/gtdb.py', "{'type': int, 'default': DEFAULT_ID_BASE}",
                         "{'type': int}")])
        _report, problems = self._judge(d)
        self.assertTrue(any('一条都没判定成' in p for p in problems), problems)


class TestNumericNoneGuards(unittest.TestCase):
    """Rule 21: a value that can reach code as None must never land in a number slot.

    The three real bugs of this project's history were all this shape -- a string
    into `isinstance(bw, float)` (round 45), a string into `cids[i:i+cq]` and
    `range()` (round 46), None into a comparison that silently skips a rewrite
    (round 50) -- and each was found by running code, not by reading it. This rule
    watches the shape. These fixtures are synthetic because the shipped tree is
    clean, and a test suite that only ever asserts "no problems" cannot tell a
    fixed rule from a deaf one.
    """

    HEAD = ('arguments = ["g",\n'
            "    ['-n|--narrow', 'how many', {'type': int}],\n"
            "]\n\n\n")

    def _tool(self):
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt21', path,
                                       loader=SourceFileLoader('vpt21', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _tree(self, body, config=None):
        d = tmpdir()
        self.addCleanup(__import__('shutil').rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'neoHGT'))
        with open(join(d, 'neoHGT', 'search.py'), 'w', encoding='utf-8') as fh:
            fh.write(self.HEAD + body)
        for mod in ('database', 'analyze', 'gtdb', 'util'):
            with open(join(d, 'neoHGT', f'{mod}.py'), 'w', encoding='utf-8') as fh:
                fh.write('arguments = []\n')
        if config is not None:
            with open(join(d, 'neoHGT', 'config.yml'), 'w', encoding='utf-8') as fh:
                fh.write(config)
        return d

    def _run(self, d):
        return self._tool().numeric_none_problems(d)

    def test_an_unguarded_range_argument_is_reported(self):
        report, problems = self._run(self._tree(
            'def run(self, args):\n'
            '    self.narrow = args.narrow\n'
            '    for i in range(self.narrow):\n'
            '        pass\n'))
        self.assertTrue(any('range' in p or 'bound' in p for p in problems),
                        problems)
        self.assertEqual(report['unguarded'][0]['name'], 'narrow', report)
        self.assertGreaterEqual(report['can_be_none'], 1, report)

    def test_a_truthiness_guard_makes_the_same_code_pass(self):
        _report, problems = self._run(self._tree(
            'def run(self, args):\n'
            '    self.narrow = args.narrow\n'
            '    if self.narrow:\n'
            '        for i in range(self.narrow):\n'
            '            pass\n'))
        self.assertEqual(problems, [],
                         'the read is the same, the guard is the difference')

    def test_a_message_is_not_arithmetic(self):
        # `'... %s' % self.x` is formatting. Reading it as a numeric use is how the
        # first run of this rule reported six false alarms on --taxdump.
        _report, problems = self._run(self._tree(
            'def run(self, args):\n'
            '    self.narrow = args.narrow\n'
            "    print('narrow is %s' % self.narrow)\n"))
        self.assertEqual(problems, [], problems)

    def test_a_config_value_counts_as_delivery_and_stops_the_report(self):
        # The delivery has to be reachable the way the real code reaches it: an
        # option that argparse leaves as None is filled from config.yml by
        # get_config(), which is why 69 of the shipped options are not None-capable.
        _report, problems = self._run(self._tree(
            'def run(self, args):\n'
            '    self.narrow = args.narrow\n'
            "    get_config(self, 'narrow', 'search.narrow')\n"
            '    for i in range(self.narrow):\n'
            '        pass\n',
            config='search:\n  narrow: 5\n'))
        self.assertEqual(problems, [], problems)

    def test_an_attribute_reassigned_from_disk_is_not_the_options_value(self):
        # `self.taxdump` starts as args.taxdump and later holds a loaded dict; the
        # name is dropped from the option's population and the drop is counted.
        report, problems = self._run(self._tree(
            'def run(self, args):\n'
            '    self.narrow = args.narrow\n'
            '    self.narrow = load_narrow(self.narrow)\n'
            '    for i in range(self.narrow):\n'
            '        pass\n'))
        self.assertEqual(problems, [], problems)
        self.assertIn('narrow', report['names_dropped_as_reassigned'], report)

    def test_a_tree_where_every_option_has_a_default_is_refused_as_blind(self):
        # Vacuity: if every one of 24 options magically has a table default, the
        # rule is not looking at the shapes it claims to look at. The guard is
        # silent on small trees on purpose (see VACUITY_MIN_OPTIONS), so this
        # fixture has to be bigger than one option -- the first version of this
        # test asserted the guard on a single-option tree and proved nothing.
        d = self._tree('def run(self, args):\n    pass\n')
        table = ['arguments = ["g",']
        for i in range(24):
            table.append("    ['-x%d|--x%d', 'n%d', {'type': int, 'default': %d}],"
                         % (i, i, i, i))
        with open(join(d, 'neoHGT', 'search.py'), 'w', encoding='utf-8') as fh:
            fh.write('\n'.join(table) + '\n]\n')
        report, problems = self._run(d)
        self.assertGreaterEqual(report['options'], 20, report)
        self.assertEqual(report['can_be_none'], 0, report)
        self.assertTrue(any('判定形态没认全' in p for p in problems), problems)

    def test_the_shipped_tree_is_clean_but_not_because_it_was_unexamined(self):
        v = self._tool()
        code = os.environ.get('NEOHGT_CODE')
        if not code:
            self.skipTest('NEOHGT_CODE is not set (bare checkout)')
        report, problems = v.numeric_none_problems(code)
        self.assertEqual(problems, [], problems)
        self.assertGreaterEqual(report['options'], 150, report)
        self.assertGreaterEqual(report['can_be_none'], 50, report)
        # the detector itself must be demonstrably alive on this tree
        self.assertGreaterEqual(report['numeric_reads_in_all_cli_attrs'], 1, report)
        self.assertEqual(report['config_state'], 'read', report)


class TestConfigKnobs(unittest.TestCase):
    """A config key nothing reads is a knob that only looks like a setting.

    Round 40's probe found `database.gtdb.keep_placeholder_names` in
    config.yml: the search stage reads
    `taxonomy.profiles.gtdb.keep_placeholder` instead, so editing it did nothing.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt6', path,
                                       loader=SourceFileLoader('vpt6', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_every_shipped_gtdb_knob_is_read(self):
        v = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        report, problems = v.config_key_problems(
            code, join(code, 'docs', '交付物清单.md'))
        self.assertGreaterEqual(len(report['declared_own_keys']), 15,
                                'almost no keys were classified (wrong path?)')
        self.assertEqual(report['never_read'], [], report)
        self.assertEqual(report['read_but_absent'], {}, report)
        self.assertEqual(problems, [])

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_a_dynamic_section_is_reported_only_until_rule_23_owns_it(self):
        """The handoff between the config-key rule and rule 23, judged both ways.

        `precomp.queries/maxchars/extrargs` are reads the code really performs and
        config.yml really lacks. Rule 22 would call that a defect unless the
        disposition table says why it is not -- so with the ledger they belong to
        rule 23 and must not be counted twice, and *without* the ledger they must
        come back. A missing ledger must never look like a clean run.
        """
        v = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        doc = join(code, 'docs', '交付物清单.md')
        self.assertEqual(v.config_key_problems(code, doc)[0]['read_but_absent'],
                         {})
        without = v.config_key_problems(code, None)[0]
        self.assertEqual(sorted(without['read_but_absent']),
                         ['precomp.extrargs', 'precomp.maxchars',
                          'precomp.queries'], without['read_but_absent'])
        self.assertTrue(without['problems'], 'an unread ledger is not a pass')
        # and rule 23 must be the one that owns them, with the shipped dispositions
        report, problems = v.dynamic_section_problems(code, doc)
        self.assertEqual(sorted(report['candidates']), ['auto', 'precomp'],
                         report)
        self.assertEqual(problems, [], problems)

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_an_unread_knob_is_named(self):
        v = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        for sub in ('neoHGT', 'tools'):
            shutil.copytree(join(HERE, os.pardir, os.pardir, sub),
                            join(d, sub))
        path = join(d, 'neoHGT', 'config.yml')
        text = open(path, encoding='utf-8').read()
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text.replace('    taxid_base:',
                                  '    a_knob_nobody_reads: 7\n    taxid_base:',
                                  1))
        report, problems = v.config_key_problems(d, None)
        self.assertIn('database.gtdb.a_knob_nobody_reads',
                      report['never_read'], problems)
        self.assertTrue(any('a_knob_nobody_reads' in p for p in problems),
                        f'the rule reported something, but not this knob: '
                        f'{problems}')


class TestPrecompMultiSampleGuard(unittest.TestCase):
    """F6: a GTDB search must not read the wrong pre-compacted table.

    Upstream hands `list(pcmap.values())[0]` to every batch, so a directory with
    several per-sample tables silently returns zero hits for all samples after
    the first. The GTDB branch refuses instead.
    """

    def _guard(self, method, pcmap):
        from types import SimpleNamespace
        from neoHGT.search import guard_precomp_multisample
        return guard_precomp_multisample(
            SimpleNamespace(method=method, pcmap=pcmap))

    def test_single_table_is_allowed(self):
        # "it did not raise" is only evidence when the case that must raise is
        # shown in the same test: a guard that never raised at all would pass a
        # test that only calls it and asserts nothing (round 37 found exactly
        # that shape here)
        self.assertIsNone(self._guard('precomp', {'sampleA': '/t/A.tsv'}))
        self.assertIsNone(self._guard('precomp', {}))
        with self.assertRaises(ValueError):
            self._guard('precomp', {'a': '/t/A', 'b': '/t/B'})

    def test_other_methods_are_untouched(self):
        for method in ('diamond', 'blastp', 'selfaln'):
            self.assertIsNone(self._guard(method, {'a': '/t/A', 'b': '/t/B'}),
                              f'{method} must not be refused')
        with self.assertRaises(ValueError):
            self._guard('precomp', {'a': '/t/A', 'b': '/t/B'})

    def test_several_tables_are_refused(self):
        with self.assertRaises(ValueError) as ctx:
            self._guard('precomp', {'a': '/t/A', 'b': '/t/B'})
        msg = str(ctx.exception)
        self.assertIn('appendix F6', msg)
        self.assertIn('zero hits', msg)
        # the refusal must say what to do instead
        self.assertIn('-s', msg)


class TestMetaWithoutDiamond(unittest.TestCase):
    """UT-2.12: `--compile none` must still publish a meta.json.

    That path never locates DIAMOND, so `self.diamond` is None; the meta writer
    used to pass it straight to shutil.which(), which raised TypeError *after*
    db.faa, the taxdump and the tables were already on disk -- leaving a
    finished-looking database with no manifest.
    """

    def test_missing_diamond_does_not_abort_meta_writing(self):
        import argparse
        from types import SimpleNamespace
        from unittest import mock
        from neoHGT import database as db
        out = tmpdir()
        self.addCleanup(shutil.rmtree, out, ignore_errors=True)
        me = SimpleNamespace(fingerprints=False, output=out,
                             gtdb_taxonmap=join(out, 'absent.gz'),
                             diamond=None, cfg={})
        seen = {}

        def fake_write_meta(path, provider, args, stats, **kw):
            seen['kw'] = kw
            with open(path, 'w') as f:
                f.write('{}')

        provider = SimpleNamespace(info={}, df_all=None)
        with mock.patch.object(db.gtdb, 'write_meta_json', fake_write_meta):
            db.write_gtdb_meta(me, provider, argparse.Namespace(), {}, 'ts')
        self.assertTrue(isfile(join(out, 'meta.json')))
        self.assertIsNone(seen['kw']['diamond_version'])

    def test_protein_reconciliation_uses_the_selection_not_the_release(self):
        """`proteins_planned` has to be the sum over the genomes actually built.

        provider.df is the *loaded* release table (load_genomes() caches every
        genome there), so reconciling against it reported the 4,896-genome
        bacterial build as "1,944,720,987 proteins short" -- the whole release
        minus the build. The build path must hand the selection to the writer.
        """
        import argparse
        from types import SimpleNamespace
        from unittest import mock
        from neoHGT import database as db
        out = tmpdir()
        self.addCleanup(shutil.rmtree, out, ignore_errors=True)
        me = SimpleNamespace(fingerprints=False, output=out,
                             gtdb_taxonmap=join(out, 'absent.gz'),
                             diamond=None, cfg={})
        release = pd.DataFrame({'proteins': [10, 20, 30, 40]})
        selection = pd.DataFrame({'proteins': [10, 20]})
        provider = SimpleNamespace(info={}, df_all=release, df=release,
                                   gtdb_dir=out)
        seen = {}

        def fake_write_meta(path, prov, args, stats, **kw):
            seen['selection'] = kw.get('selection')
            with open(path, 'w') as f:
                f.write('{}')

        with mock.patch.object(db.gtdb, 'write_meta_json', fake_write_meta):
            db.write_gtdb_meta(me, provider, argparse.Namespace(),
                               {'proteins': 30}, 'ts', selection)
        self.assertIs(seen['selection'], selection,
                      'the manifest must reconcile against the selection')
        self.assertEqual(
            gtdb.reconcile_protein_counts(seen['selection'],
                                          {'proteins': 30})[2], 'ok')
        self.assertEqual(gtdb.reconcile_protein_counts(
            release, {'proteins': 30})[2], '70 short',
            'the release table is the wrong basis and must not be used')

    def test_diamond_version_still_recorded_when_present(self):
        import argparse
        from types import SimpleNamespace
        from unittest import mock
        from neoHGT import database as db
        out = tmpdir()
        self.addCleanup(shutil.rmtree, out, ignore_errors=True)
        me = SimpleNamespace(fingerprints=False, output=out,
                             gtdb_taxonmap=join(out, 'absent.gz'),
                             diamond='sh', cfg={})
        seen = {}
        with mock.patch.object(db.gtdb, 'write_meta_json',
                               lambda p, *a, **kw: seen.update(kw) or
                               open(p, 'w').write('{}')):
            with mock.patch.object(db, 'run_command',
                                   lambda cmd: (0, ['diamond v2.1.10'])):
                db.write_gtdb_meta(me, SimpleNamespace(info={}, df_all=None),
                                   argparse.Namespace(), {}, 'ts')
        self.assertIn('2.1.10', seen['diamond_version'])


class TestScratchDirectoryStaysInWorkspace(unittest.TestCase):
    """UT-5.x: an omitted --tmpdir must not scatter scratch into /tmp.

    Upstream defaults it to `mkdtemp()` wherever the ambient TMPDIR points, so
    the workspace rule (§5.8) held only because tools/env.sh exports TMPDIR
    inside the project. The GTDB/hybrid branch now relocates the directory, and
    only that branch: the NCBI route keeps the upstream default untouched.
    """

    def setUp(self):
        import argparse
        from types import SimpleNamespace
        self.argparse, self.ns = argparse, SimpleNamespace
        self.root = tmpdir()
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.code = join(self.root, 'neoHGT2GTDB-项目代码')
        self.cache = join(self.root, 'neoHGT2GTDB-缓存测试')
        self.tmp = join(self.root, 'neoHGT2GTDB-临时文件')
        self.dataset = join(self.root, 'GTDB_r220')
        for d in (self.code, self.cache, self.tmp, self.dataset,
                  join(self.cache, 'db')):
            os.makedirs(d, exist_ok=True)
        self.out = join(self.cache, 'db', 'x')
        self.argv = argparse.Namespace(tmpdir=None)

    def _configure(self, argv_tmpdir=None):
        from tempfile import mkdtemp
        from unittest import mock
        from neoHGT.database import set_source_parameters
        me = self.ns(cfg={'workspace': {'cache_dir': self.cache,
                                         'tmp_dir': self.tmp}},
                     source='gtdb', seq_source='local',
                     gtdb_dir=self.dataset, gtdb_taxa=None,
                     output=self.out, tmpdir=None)
        # emulate what upstream set_parameters() does when --tmpdir is omitted
        upstream = None if argv_tmpdir else mkdtemp()
        me.tmpdir = argv_tmpdir or upstream
        if upstream:
            setattr(me, 'mkdtemp', True)
        argv = self.argparse.Namespace(tmpdir=argv_tmpdir)
        with mock.patch.dict(os.environ, {'NEOHGT_ROOT': self.root}):
            set_source_parameters(me, argv)
        return me, upstream

    def test_relocated_under_workspace_tmp_dir(self):
        me, upstream = self._configure()
        self.assertTrue(os.path.realpath(me.tmpdir).startswith(
            os.path.realpath(self.tmp)), me.tmpdir)
        self.assertNotEqual(me.tmpdir, upstream)
        # the empty directory upstream created is gone, and only because rmdir
        # cannot remove anything that still has contents
        self.assertFalse(exists(upstream))

    def test_explicit_tmpdir_is_left_alone(self):
        chosen = join(self.tmp, 'chosen')
        os.makedirs(chosen)
        me, _ = self._configure(argv_tmpdir=chosen)
        self.assertEqual(me.tmpdir, chosen)


class TestTaskLedgerGate(unittest.TestCase):
    """UT: the §6 numbered-task ledger must be mechanically complete."""

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt', path,
                                       loader=SourceFileLoader('vpt', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_every_plan_option_is_handed_to_the_parser(self):
        """A planned option that no test parses is an unverified option.

        'Exists in the CLI' and 'is read by the code' cannot see a flag that
        parses into the wrong dest or whose value never reaches the selection;
        only handing it to the real parser does. The rule is checked both ways:
        that it bites on a made-up flag, and that the shipped suite really
        covers every §7 option.
        """
        tool = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        # built at runtime: a flag written out as a literal here would appear in
        # the very file the rule scans, and the control would prove nothing
        ghost = '--' + 'an-option' + '-nobody-tests'
        self.assertEqual(tool.flags_without_tests([ghost], code), [ghost],
                         'the coverage rule must actually bite')
        plan = join(code, os.pardir, 'neoHGT2.GTDB-ADAPTATION-PLAN-CHN.md')
        if not isfile(plan):
            self.skipTest('the plan document is not in this checkout')
        flags = tool.section7_flags(open(plan, encoding='utf-8').read())
        self.assertGreater(len(flags), 20, flags)
        self.assertEqual(tool.flags_without_tests(flags, code), [],
                         '方案 §7 里某个选项没有任何用例把它交给解析器')

    def setUp(self):
        self.root = tmpdir()
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.code = join(self.root, 'neoHGT2GTDB-项目代码')
        os.makedirs(join(self.code, 'neoHGT'))
        with open(join(self.code, 'neoHGT', 'gtdb.py'), 'w') as f:
            f.write('# a real file the ledger may cite\nline two\nline three\n'
                    'def here_ok(x):\n    return x\n')
        os.makedirs(join(self.root, 'neoHGT2GTDB-缓存测试', 'reports'),
                    exist_ok=True)
        with open(join(self.root, 'neoHGT2GTDB-缓存测试', 'reports', 'e.json'),
                  'w') as f:
            f.write('{}')
        self.plan = '\n'.join((
            '## 6. 分阶段实施计划', '', '| 编号 | 任务 |', '| --- | --- |',
            '| S0 | spike |', '| 2.1 | wiring |', '| 2.13 | lists |',
            '', '## 7. CLI 与配置变更清单', ''))

    def _check(self, rows, statuses=None, ids=None):
        v = self._tool()
        docs = join(self.code, 'docs')
        os.makedirs(docs, exist_ok=True)
        path = join(docs, 'ledger.md')
        with open(path, 'w', encoding='utf-8') as f:
            f.write('| 编号 | 状态 | 实现位置 | 证据 |\n'
                    '| --- | --- | --- | --- |\n' + rows)
        kwargs = {} if statuses is None else {'statuses': statuses}
        return v.ledger_problems(v.ledger_rows(path),
                                 v.plan_task_ids(self.plan) if ids is None
                                 else ids, self.code, self.root, **kwargs)

    def test_plan_ids_are_extracted_from_section_6_only(self):
        v = self._tool()
        self.assertEqual(v.plan_task_ids(self.plan), ['S0', '2.1', '2.13'])

    def test_task_missing_from_ledger_is_reported(self):
        unlisted, extra, problems = self._check(
            '| S0 | 已实现 | `neoHGT/gtdb.py` | `缓存测试/reports/e.json` |\n')
        self.assertEqual(unlisted, ['2.1', '2.13'])
        self.assertEqual(extra, [])
        self.assertEqual(problems, [])

    def test_invented_ledger_id_is_reported(self):
        _, extra, _ = self._check(
            '| S0 | 已实现 | `neoHGT/gtdb.py` | x |\n'
            '| 9.9 | 已实现 | `neoHGT/gtdb.py` | x |\n')
        self.assertEqual(extra, ['9.9'])

    def test_unknown_status_is_reported(self):
        _, _, problems = self._check(
            '| S0 | 大概做完了 | `neoHGT/gtdb.py` | 证据足够长的一条说明文字内容 |\n'
            '| 2.1 | 已实现 | `neoHGT/gtdb.py` | x |\n'
            '| 2.13 | 已实现 | `neoHGT/gtdb.py` | x |\n')
        self.assertTrue(any('状态' in p for p in problems), problems)

    def test_partial_needs_an_explicit_reason(self):
        _, _, problems = self._check(
            '| S0 | 部分实现 | `neoHGT/gtdb.py` | 嗯 |\n'
            '| 2.1 | 已实现 | `neoHGT/gtdb.py` | x |\n'
            '| 2.13 | 已实现 | `neoHGT/gtdb.py` | x |\n')
        self.assertTrue(any('缺的是哪一块' in p for p in problems), problems)

    def test_ghost_path_is_reported(self):
        _, _, problems = self._check(
            '| S0 | 已实现 | `neoHGT/nope.py` | `缓存测试/reports/e.json` |\n'
            '| 2.1 | 已实现 | `neoHGT/gtdb.py` | `缓存测试/reports/ghost.json` |\n'
            '| 2.13 | 已实现 | `neoHGT/gtdb.py` | x |\n')
        self.assertEqual(len(problems), 2, problems)

    def test_workspace_shorthand_resolves(self):
        _, _, problems = self._check(
            '| S0 | 已实现 | `neoHGT/gtdb.py` | `缓存测试/reports/e.json` |\n'
            '| 2.1 | 已实现 | `neoHGT/gtdb.py` | x |\n'
            '| 2.13 | 已实现 | `neoHGT/gtdb.py` | x |\n')
        self.assertEqual(problems, [])

    def test_line_citations_are_checked_too(self):
        # a citation may name a line, and the gate notices when the file no
        # longer has it -- this is how the defect register points at code
        _, _, problems = self._check(
            '| S0 | 已实现 | `neoHGT/gtdb.py:3` | x |\n'
            '| 2.1 | 已实现 | `gtdb.py:2` | x |\n'      # bare name resolves
            '| 2.13 | 已实现 | `neoHGT/gtdb.py:99` | x |\n')
        self.assertEqual(
            [p for p in problems if '行号' in p],
            ['2.13: `neoHGT/gtdb.py:99` 的行号超出文件长度'], problems)

    def test_defect_register_has_its_own_status_vocabulary(self):
        v = self._tool()
        _, _, ok = self._check(
            '| S0 | 已修复 | `neoHGT/gtdb.py` | 说明写得足够长所以不该被挑出来 |\n',
            statuses=v.DEFECT_STATUSES, ids=['S0'])
        self.assertEqual(ok, [])
        _, _, bad = self._check(
            '| S0 | 已实现 | `neoHGT/gtdb.py` | x |\n',
            statuses=v.DEFECT_STATUSES, ids=['S0'])
        self.assertTrue(any('不在' in p for p in bad), bad)
        # 登记未改 is one of the statuses that must explain itself
        _, _, short = self._check(
            '| S0 | 登记未改 | `neoHGT/gtdb.py` | 懒 |\n',
            statuses=v.DEFECT_STATUSES, ids=['S0'])
        self.assertTrue(any('没有写清' in p for p in short), short)

    def test_appendix_flags_must_be_explained_somewhere(self):
        v = self._tool()
        plan = ('### 附录 C:示例命令\n'
                '```sh\nneoHGT analyze -i x --to ncbi --frobnicate yes\n'
                'rsync --list-only rsync://h/releases/\n```\n### 附录 E\n')
        flags = v.appendix_flags(plan)
        self.assertEqual(flags, ['--frobnicate', '--list-only', '--to'])
        # a third-party option is not our CLI's business
        self.assertIn('--list-only', v.THIRD_PARTY_FLAGS)
        docs = join(self.code, 'docs')
        os.makedirs(docs, exist_ok=True)
        with open(join(docs, 'd.md'), 'w', encoding='utf-8') as f:
            f.write('## 附录 C/D 命令形态偏差\n\n'
                    '| 方案写法 | 实际形态 |\n| --- | --- |\n'
                    '| `--to ncbi` | `--donor-ncbi-map yes` |\n')
        self.assertEqual(v.deviation_flags(join('docs', 'd.md'), self.code),
                         {'--to'})

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_deliverable_register_is_structurally_complete(self):
        # the shipped register must cover the plan's 12 deliverables with the
        # right vocabulary; path existence is checked by the gate itself, which
        # needs the cache artifacts and therefore cannot run in CI
        v = self._tool()
        rows = v.ledger_rows(join(HERE, os.pardir, os.pardir, 'docs',
                                  '交付物清单.md'))
        d = {k: r for k, r in rows.items() if k.startswith('D')}
        self.assertEqual(sorted(d, key=lambda x: int(x[1:])),
                         [f'D{i}' for i in range(1, 13)])
        self.assertTrue(all(r['status'] in v.DELIVER_STATUSES
                            for r in d.values()),
                        {k: r['status'] for k, r in d.items()})
        # 部分交付/未交付 must say what is missing
        self.assertTrue(all(len(r['evidence']) >= v.MIN_REASON
                            for r in d.values()
                            if r['status'] != v.DELIVER_STATUSES[0]))


    def test_appendix_b_symbols_resolve_or_are_declared(self):
        v = self._tool()
        plan = ('### 附录 B:签名\n```\nneoHGT/gtdb.py\n'
                '  def here_ok(x)\n  MISSING_SYMBOL = 1\n'
                '  DECLARED_ELSEWHERE = 1\n'
                'neoHGT/util.py\n  def in_other_module(x)\n```\n'
                '### 附录 C\n')
        with open(join(self.code, 'neoHGT', 'util.py'), 'w') as f:
            f.write('def in_other_module(x):\n    return x\n')
        docs = join(self.code, 'docs')
        os.makedirs(docs, exist_ok=True)
        with open(join(docs, 'd.md'), 'w', encoding='utf-8') as f:
            f.write('## 附录 B 符号偏差\n\n| 方案符号 | 实际实现 |\n'
                    '| --- | --- |\n'
                    '| `DECLARED_ELSEWHERE = (1)` | 换名了 |\n')
        syms = v.appendix_b_symbols(plan)
        self.assertEqual([n for _, n in syms],
                         ['here_ok', 'MISSING_SYMBOL', 'DECLARED_ELSEWHERE',
                          'in_other_module'])
        problems = v.appendix_b_problems(
            syms, self.code, join('docs', 'd.md'))
        self.assertEqual(problems,
                         ['附录 B 承诺的 MISSING_SYMBOL（neoHGT/gtdb.py）'
                          '既不存在也未登记偏差'])

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_acceptance_register_covers_every_dod(self):
        v = self._tool()
        plan = ('### 1.4 成功标准\n| 编号 | 判据 |\n| --- | --- |\n'
                '| DoD-1 | a |\n| DoD-9 | b |\n\n## 2. x\n'
                '### 10.7 验收 Checklist\n- [ ] DoD-0 本地布局探测\n'
                '- [ ] DoD-1 建库\n\n## 11. 预算\n')
        self.assertEqual(v.plan_dod_ids(plan),
                         ['DoD-0', 'DoD-1', 'DoD-9'])
        rows = v.ledger_rows(join(HERE, os.pardir, os.pardir, 'docs',
                                  '验收标准对照.md'))
        dod = {k: r for k, r in rows.items() if k.startswith('DoD-')}
        self.assertEqual(sorted(dod, key=lambda x: int(x.split('-')[1])),
                         [f'DoD-{i}' for i in range(10)])
        self.assertTrue(all(r['status'] in v.DOD_STATUSES
                            for r in dod.values()),
                        {k: r['status'] for k, r in dod.items()})
        # 部分达标/被阻塞 must name the missing half
        self.assertTrue(all(len(r['evidence']) >= v.MIN_REASON
                            for r in dod.values()
                            if r['status'] != v.DOD_STATUSES[0]))

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_decision_register_covers_every_adr(self):
        v = self._tool()
        plan = ('## 8. 关键设计决策(ADR)\n| 编号 | 决策 |\n| --- | --- |\n'
                '| ADR-1 | a |\n| ADR-4b | b |\n| ADR-10 | c |\n\n'
                '## 9. 风险登记册\n')
        self.assertEqual(v.plan_adr_ids(plan), ['ADR-1', 'ADR-4b', 'ADR-10'])
        rows = v.ledger_rows(join(HERE, os.pardir, os.pardir, 'docs',
                                  '决策落点对照.md'))
        adr = {k: r for k, r in rows.items() if k.startswith('ADR-')}
        self.assertEqual(sorted(adr, key=lambda x: (int(x[4:].rstrip('b')), x)),
                         ['ADR-1', 'ADR-2', 'ADR-3', 'ADR-4', 'ADR-4b',
                          'ADR-5', 'ADR-6', 'ADR-7', 'ADR-8', 'ADR-9',
                          'ADR-10'])
        self.assertTrue(all(r['status'] in v.DECISION_STATUSES
                            for r in adr.values()),
                        {k: r['status'] for k, r in adr.items()})
        self.assertTrue(all(len(r['evidence']) >= v.MIN_REASON
                            for r in adr.values()
                            if r['status'] != v.DECISION_STATUSES[0]))

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_risk_register_covers_every_plan_risk(self):
        v = self._tool()
        plan = ('## 9. 风险登记册\n| ID | 风险 | 缓解 |\n| --- | --- | --- |\n'
                '| RK1 | a | m |\n| RK12 | b | m |\n\n## 10. 测试与验证方案\n')
        self.assertEqual(v.plan_risk_ids(plan), ['RK1', 'RK12'])
        rows = v.ledger_rows(join(HERE, os.pardir, os.pardir, 'docs',
                                  '风险应对登记.md'))
        rk = {k: r for k, r in rows.items() if k.startswith('RK')}
        self.assertEqual(sorted(rk, key=lambda x: int(x[2:])),
                         [f'RK{i}' for i in range(1, 20)])
        self.assertTrue(all(r['status'] in v.RISK_STATUSES
                            for r in rk.values()),
                        {k: r['status'] for k, r in rk.items()})
        # 部分实现/登记接受 must say which half is missing
        self.assertTrue(all(len(r['evidence']) >= v.MIN_REASON
                            for r in rk.values()
                            if r['status'] != v.RISK_STATUSES[0]))

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_verification_register_covers_the_plan_matrix(self):
        v = self._tool()
        sect = '\n'.join((
            '### 10.2 单元 / 属性测试', '| 编号 | 对象 |', '| --- | --- |',
            '| UT-0.1 | a |', '| UT-3.1b | b |', '| UT-2.x | c |',
            '### 10.3 集成测试', '| 编号 | 场景 |', '| --- | --- |',
            '| IT-0 | x |', '| IT-2b | y |',
            '### 10.4 科学有效性验证', '1. **一致性**: a', '2. 参考: b',
            '### 10.5 回归与兼容', '- RT-1: a', '- RT-3: b', '### 10.6 基准'))
        self.assertEqual(v.plan_verification_ids(sect),
                         ['UT-0.1', 'UT-3.1b', 'UT-2.x', 'IT-0', 'IT-2b',
                          '10.4-1', '10.4-2', 'RT-1', 'RT-3'])
        plan_file = join(HERE, os.pardir, os.pardir, os.pardir,
                         'neoHGT2.GTDB-ADAPTATION-PLAN-CHN.md')
        if not isfile(plan_file):
            self.skipTest('the plan document is not next to the code tree')
        ids = v.plan_verification_ids(open(plan_file, encoding='utf-8').read())
        self.assertEqual(len(ids), 36, ids)
        rows = v.ledger_rows(join(HERE, os.pardir, os.pardir, 'docs',
                                  '测试用例对照.md'))
        self.assertEqual(sorted(rows), sorted(ids))
        self.assertTrue(all(r['status'] in v.VERIF_STATUSES
                            for r in rows.values()),
                        {k: r['status'] for k, r in rows.items()})
        # 部分覆盖/未覆盖 must name the missing half
        self.assertTrue(all(len(r['evidence']) >= v.MIN_REASON
                            for r in rows.values()
                            if r['status'] != v.VERIF_STATUSES[0]))
        self.assertIn('test_gtdb_pipeline.py', rows['IT-4']['where'],
                      'the offline end-to-end row must cite the test itself')

    def test_citation_drift_into_a_function_body_is_caught(self):
        v = self._tool()
        tdir = join(self.code, 'neoHGT', 'tests')
        os.makedirs(tdir, exist_ok=True)
        with open(join(tdir, 'test_x.py'), 'w') as f:
            f.write('def test_a():\n    body_line\n')
        rows = {'S1': {'id': 'S1', 'status': v.LEDGER_STATUSES[0],
                       'where': '`neoHGT/tests/test_x.py:1`',
                       'evidence': 'points at the definition itself'}}
        self.assertEqual(
            v.ledger_problems(rows, ['S1'], self.code, self.root,
                              v.LEDGER_STATUSES)[2], [])
        rows['S1']['where'] = '`neoHGT/tests/test_x.py:2`'
        problems = v.ledger_problems(rows, ['S1'], self.code, self.root,
                                     v.LEDGER_STATUSES)[2]
        self.assertTrue(any('漂移' in p for p in problems), problems)

    def test_frozen_plan_digest_equals_the_plan_it_came_from(self):
        """CI reads the digest instead of the plan, so the digest must not rot."""
        v = self._tool()
        digest = join(HERE, os.pardir, os.pardir, 'docs', '方案编号清单.json')
        plan_file = join(HERE, os.pardir, os.pardir, os.pardir,
                         'neoHGT2.GTDB-ADAPTATION-PLAN-CHN.md')
        if not (isfile(digest) and isfile(plan_file)):
            self.skipTest('needs both the digest and the plan document')
        saved = json.load(open(digest, encoding='utf-8'))
        plan_text = open(plan_file, encoding='utf-8').read()
        self.assertEqual(saved['inventory'], v.plan_inventory(plan_text),
                         'the committed id digest no longer matches the plan; '
                         'regenerate it with --emit-ids')
        self.assertEqual(saved['counts']['verification'], 36)
        self.assertEqual(saved['counts']['tasks'], 28)

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_the_assumption_and_revision_registers_are_complete(self):
        """§2.3 (S1–S9) and §0.5 (REV-01–19) are numbered and truth-apt."""
        v = self._tool()
        docs = join(HERE, os.pardir, os.pardir, 'docs')
        plan_file = join(HERE, os.pardir, os.pardir, os.pardir,
                         'neoHGT2.GTDB-ADAPTATION-PLAN-CHN.md')
        if not isfile(plan_file):
            self.skipTest('the plan document is not next to the code tree')
        text = open(plan_file, encoding='utf-8').read()
        ids, revs = v.plan_assumption_ids(text), v.plan_revision_ids(text)
        self.assertEqual(ids, [f'S{i}' for i in range(1, 10)])
        self.assertEqual(revs, [f'REV-{i:02d}' for i in range(1, 20)])
        a = v.ledger_rows(join(docs, '搜索层假设核对.md'))
        r = v.ledger_rows(join(docs, '方案修订记录核对.md'))
        self.assertEqual(sorted(a), sorted(ids))
        self.assertEqual(sorted(r), sorted(revs))
        self.assertTrue(all(x['status'] in v.ASSUMPTION_STATUSES
                            for x in a.values()))
        self.assertTrue(all(x['status'] in v.REVISION_STATUSES
                            for x in r.values()))
        self.assertIn('makeblastdb', a['S6']['evidence'],
                      'S6 must record that the plan judged BLAST feasible and '
                      'measurement said otherwise')
        self.assertEqual(r['REV-12']['status'], '部分采纳')
        self.assertGreaterEqual(len(r['REV-12']['evidence']), v.MIN_REASON)

    def test_wildcard_citations_are_rejected(self):
        """`dir/*.json` cannot be opened, so it must not count as evidence."""
        v = self._tool()
        rows = {'S1': {'id': 'S1', 'status': v.LEDGER_STATUSES[0],
                       'where': '`neoHGT/gtdb.py`',
                       'evidence': '见 `缓存测试/reports/e*.json`'}}
        problems = v.ledger_problems(rows, ['S1'], self.code, self.root,
                                     v.LEDGER_STATUSES)[2]
        self.assertTrue(any('通配符' in p for p in problems), problems)
        rows['S1']['evidence'] = '见 `缓存测试/reports/e.json`'
        self.assertEqual(
            v.ledger_problems(rows, ['S1'], self.code, self.root,
                              v.LEDGER_STATUSES)[2], [])

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_milestone_register_covers_the_plan_table(self):
        v = self._tool()
        plan_file = join(HERE, os.pardir, os.pardir, os.pardir,
                         'neoHGT2.GTDB-ADAPTATION-PLAN-CHN.md')
        if not isfile(plan_file):
            self.skipTest('the plan document is not next to the code tree')
        ids = v.plan_milestone_ids(open(plan_file, encoding='utf-8').read())
        self.assertEqual(len(ids), 8, ids)
        rows = v.ledger_rows(join(HERE, os.pardir, os.pardir, 'docs',
                                  '里程碑核对.md'))
        self.assertEqual(sorted(rows), sorted(ids))
        self.assertTrue(all(r['status'] in v.MILESTONE_STATUSES
                            for r in rows.values()))
        # the two genuinely open milestones must name what is missing
        for mid in ('W7', 'W8'):
            self.assertEqual(rows[mid]['status'], '部分达成', mid)
            self.assertGreaterEqual(len(rows[mid]['evidence']), v.MIN_REASON)

    def test_report_paths_outside_the_workspace_are_refused(self):
        """Rule 11: this tool must not write outside the workspace either."""
        v = self._tool()
        with self.assertRaises(SystemExit) as ctx:
            v.guarded_out('/tmp/definitely_outside_workspace.json', self.code)
        self.assertIn('outside the workspace', str(ctx.exception))
        inside = join(self.root, 'reports', 'x.json')
        self.assertEqual(v.guarded_out(inside, self.code),
                         os.path.realpath(inside))


class TestProbeCommandOffline(FakeInput):
    """`neoHGT gtdb probe` has to run offline too.

    IT-0 was only ever executed against the read-only dataset, whose evidence
    lives in the workspace tree, so no test executed `cmd_probe` at all -- the
    subcommand a user runs *first*.
    """

    def test_probe_writes_a_machine_readable_report(self):
        renamed_tree(self.d)
        out = join(self.d, 'layout_probe.json')
        rc = gtdb.main(['probe', '--dir', self.d, '--out', out])
        self.assertIn(rc, (0, 1), f'probe returned {rc}')
        rep = json.load(open(out))
        self.assertTrue(gtdb.is_renamed_layout(rep['layout']), rep['layout'])
        self.assertEqual(str(rep['release']), '220')
        self.assertEqual(rep['serial_anchor_coverage'], 1.0)
        checks = rep['selfcheck']['checks']
        self.assertGreaterEqual(len(checks), 8)
        names = {c['name'] for c in checks}
        for want in ('lineage_depth_is_7', 'rank_prefix_vocabulary',
                     'faa_vs_metadata_reps[bacteria]',
                     'faa_vs_metadata_reps[archaea]'):
            self.assertIn(want, names)

    def test_probe_names_the_unusable_dataset_instead_of_crashing(self):
        ghost = join(self.d, 'not-a-gtdb')
        os.makedirs(ghost)
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rc = gtdb.cli_main(['probe', '--dir', ghost,
                                '--out', join(self.d, 'x.json')])
        self.assertEqual(rc, 2, err.getvalue())
        self.assertIn('ERROR:', err.getvalue())
        self.assertFalse(isfile(join(self.d, 'x.json')),
                         'a refused dataset must write nothing')


class TestBuildArtifactChoices(unittest.TestCase):
    """Choice options whose promise is about the *artifacts* of a build.

    `--keep-faa`, `--input-fingerprints` and `--seq-source` are not filters on
    the genome table but switches on what the build leaves on disk. Each value
    is run through a real (fixture-scale, `--compile none`) build so the claim
    is about files, not about a code path being reachable.
    """

    @classmethod
    def setUpClass(cls):
        from neoHGT.database import Database
        cls.d = tmpdir()
        src = join(cls.d, 'src')
        for sub in ('tree_bac120', 'tree_ar53', 'metadata_bac120',
                    'metadata_ar53',
                    '00_GTDB_r220_aa/protein_faa_reps/bacteria',
                    '00_GTDB_r220_aa/protein_faa_reps/archaea'):
            os.makedirs(join(src, sub), exist_ok=True)
        renamed_tree(src)
        cls.src, cls.db = src, Database

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, ignore_errors=True)

    def _build(self, tag, *extra):
        out = join(self.d, tag)
        scratch = join(self.d, tag + '_tmp')
        os.makedirs(scratch, exist_ok=True)
        argv = ['-o', out, '--source', 'gtdb', '--gtdb-dir', self.src,
                '--reps', '--seq-source', 'local', '--tmpdir', scratch,
                '--compile', 'none'] + list(extra)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.db()(cli_args('database', argv))
        return out, buf.getvalue()

    def test_keep_faa_keep_leaves_the_fasta(self):
        out, _log = self._build('faa_keep', '--keep-faa', 'keep')
        self.assertTrue(isfile(join(out, 'db.faa')))
        self.assertFalse(exists(join(out, 'db.faa.gz')))

    def test_keep_faa_compress_replaces_it(self):
        out, log = self._build('faa_gz', '--keep-faa', 'compress')
        self.assertFalse(isfile(join(out, 'db.faa')),
                         'compress must not leave the plaintext file behind')
        self.assertTrue(isfile(join(out, 'db.faa.gz')))
        self.assertIn('--keep-faa compress', log)

    def test_keep_faa_drop_deletes_the_file_and_says_so(self):
        # 'drop' used to print a refusal while keeping the file: a tier that
        # could not do what its name says. It now honours the explicit request
        # and names the surviving artifact.
        out, log = self._build('faa_drop', '--keep-faa', 'drop')
        self.assertIn('--keep-faa drop', log)
        self.assertFalse(isfile(join(out, 'db.faa')),
                         '--keep-faa drop must remove db.faa')
        self.assertTrue(isfile(join(out, 'taxdump', 'nodes.dmp')),
                        'the compiled artifacts must survive the drop')

    def test_input_fingerprints_yes_writes_the_table(self):
        out, _log = self._build('fp_yes', '--input-fingerprints', 'yes')
        self.assertTrue(isfile(join(out, 'gtdb', 'input_fingerprints.tsv')))

    def test_stream_extraction_yes_is_what_the_gtdb_route_does(self):
        out, log = self._build('stream_yes', '--stream-extraction', 'yes')
        self.assertTrue(isfile(join(out, 'db.faa')))
        self.assertNotIn('is an NCBI-route switch', log)

    def test_local_layout_gtdb_renamed_is_accepted(self):
        args = cli_args('database', ['-o', 'db', '--source', 'gtdb',
                                     '--gtdb-dir', 'x', '--local-layout',
                                     'gtdb-r220-renamed'])
        self.assertEqual(vars(args)['local_layout'], 'gtdb-r220-renamed')
        info = gtdb.detect_layout(self.src, layout='gtdb-r220-renamed',
                                 require_sequences=False)
        self.assertEqual(info['layout'], 'gtdb-r220-renamed')

    def test_local_layout_category_name_is_accepted(self):
        # the renamed-layout label stopped encoding a release number; the
        # category name must parse and detect_layout must emit it verbatim
        args = cli_args('database', ['-o', 'db', '--source', 'gtdb',
                                     '--gtdb-dir', 'x', '--local-layout',
                                     'local-renamed'])
        self.assertEqual(vars(args)['local_layout'], 'local-renamed')
        info = gtdb.detect_layout(self.src, layout='local-renamed',
                                  require_sequences=False)
        self.assertEqual(info['layout'], 'local-renamed')
        self.assertTrue(gtdb.is_renamed_layout(info['layout']))

    def test_gtdb_mirror_auto_prefers_the_verifiable_host(self):
        args = cli_args('database', ['-o', 'db', '--source', 'gtdb',
                                     '--gtdb-dir', 'x', '--gtdb-mirror',
                                     'auto'])
        self.assertEqual(vars(args)['gtdb_mirror'], 'auto')
        self.assertIn('ecogenomic', gtdb.gtdb_host('auto'))
        self.assertNotIn('aaupol', gtdb.gtdb_host('auto'))

    def test_seq_source_gtdb_tar_keeps_its_url_and_kind(self):
        from types import SimpleNamespace
        from neoHGT import database
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'gtdb'))
        me = SimpleNamespace(output=d, seq_source='gtdb-tar',
                             gtdb_release='r220', gtdb_mirror='ecogenomic',
                             downloader='curl', verify_md5=True,
                             allow_network=False)
        fetch = [gtdb.GenomeFetch('RS_GCF_000000001.1', 'G1', '100000001',
                                  path=None, url=None, kind='gtdb-tar')]
        with self.assertRaises(ValueError) as ctx:
            database.plan_remote_fetches(me, fetch)
        self.assertIn('--allow-network', str(ctx.exception))
        rows = [l.rstrip('\n').split('\t')
                for l in open(join(d, 'gtdb', 'fetch_plan.tsv'))
                if not l.startswith('#')]
        self.assertEqual(rows[1][2], 'remote-gtdb-tar',
                         'the plan must record the channel that was asked for')
        self.assertIn('gtdb_proteins_aa_reps.tar.gz', rows[1][1],
                      'plan §5.4: the gtdb-tar channel is the release bundle, '
                      'not one HTTP request per genome')
        self.assertIn('tar -xzO', rows[1][5])
        self.assertIn('genomic_files_reps/RS_GCF_000000001.1_protein.faa.gz',
                      rows[1][5])

    def test_expand_members_reaches_the_build_decision(self):
        # the expansion itself is covered by TestExpandMembers; what used to be
        # unproven was that the CLI value arrives as a boolean at all: 'no' is
        # a *string* and any truthiness test on it would be true
        from neoHGT import util
        for value, want in (('yes', True), ('no', False)):
            args = cli_args('database', ['-o', 'db', '--source', 'gtdb',
                                         '--gtdb-dir', 'x',
                                         '--expand-members', value])
            self.assertIs(util.arg2bool(args.expand_members), want)
        with open(join(HERE, os.pardir, 'database.py'),
                  encoding='utf-8') as fh:
            src = fh.read()
        self.assertIn("arg2bool(getattr(args, 'expand_members', None))", src,
                      'the build must gate expand_selection on the parsed flag')


class TestDatabaseChoiceValues(unittest.TestCase):
    """Every declared value of a database-stage choice option must do something.

    Round 29 found `--mimag any`: argparse listed the value, the selector did
    not know it, and no test asked for it, so the mismatch stayed invisible
    until an unrelated change made it raise. Round 30's gate rule
    `tools/verify_plan_traceability.py::choice_values_without_tests` requires
    every declared value to be handed to its option; this class is the
    behavioural half -- the value has to change the outcome, not merely parse.
    """

    LINS = ('d__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G1;s__G1 a',
            'd__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G2;s__G2 b',
            'd__Bacteria;p__P1;c__C2;o__O2;f__F2;g__G3;s__G3 c',
            'd__Archaea;p__A1;c__A1;o__A1;f__A1;g__A1;s__A1 d')
    QUALITY = (99.0, 90.0, 80.0, 70.0)

    def setUp(self):
        rows = []
        for i, (lin, q) in enumerate(zip(self.LINS, self.QUALITY), 1):
            rows.append({'genome': f'ARC_1' if lin.startswith('d__Arch')
                         else f'BAC_{i}',
                         'accession': f'GB_GCA_{i:09d}.1',
                         'accnov': f'GCA_{i:09d}',
                         'gtdb_taxonomy': lin,
                         'domain': 'archaea' if lin.startswith('d__Arch')
                         else 'bacteria',
                         'is_rep': True, 'proteins': 10,
                         'completeness': q, 'contamination': 1.0,
                         'quality_score': q - 5.0,
                         'mimag_high_quality': 't' if q >= 95 else 'f',
                         'mimag_medium_quality': 't' if 80 <= q < 95 else 'f',
                         'mimag_low_quality': 't' if q < 80 else 'f',
                         'gtdb_type_species_of_genus': 'f',
                         'ncbi_type_material_designation': 'na'})
        self.df = pd.DataFrame(rows)
        p = gtdb.GtdbProvider.__new__(gtdb.GtdbProvider)
        p.taxid_map, _ = gtdb.assign_ids(gtdb.collect_nodes(list(self.LINS)))
        p.opts, p.selection_labels, p.sample_note = {}, [], None
        self.p = p

    def test_domain_choice_restricts_the_table(self):
        # the value comes from the real CLI, through the same opts mapping the
        # build installs, so 'both' is proven to be the "no restriction" case
        for value, want in (('bacteria', ['BAC_1', 'BAC_2', 'BAC_3']),
                            ('archaea', ['ARC_1']),
                            ('both', ['ARC_1', 'BAC_1', 'BAC_2', 'BAC_3'])):
            args = cli_args('database', ['-o', 'db', '--source', 'gtdb',
                                         '--gtdb-dir', 'x', '--domain', value])
            opts = gtdb.selection_opts(args)
            got = sorted(self.p.select(self.df, opts)['genome'])
            self.assertEqual(got, want, f'--domain {value} kept {got}')

    def test_magim_keeps_exactly_the_right_tier_for_each_value(self):
        tiers = {'high': ['BAC_1'], 'medium': ['BAC_2', 'BAC_3'],
                 'low': ['ARC_1'], 'any': ['ARC_1', 'BAC_1', 'BAC_2', 'BAC_3']}
        for tier, want in tiers.items():
            got = sorted(self.p.select(self.df, {'mimag': tier})['genome'])
            self.assertEqual(got, want, f'--mimag {tier} kept {got}')

    def test_rank_sample_groups_at_every_declared_rank(self):
        # 4 genomes: two share a genus boundary only at genus, one is archaeal,
        # so each rank must collapse a different number of groups
        groups = {'domain': 2, 'phylum': 2, 'class': 3, 'order': 3,
                  'family': 3, 'genus': 4}
        for rank, n in groups.items():
            got = self.p.sample(self.df, {'rank_sample': rank})
            self.assertEqual(got.shape[0], n, f'--rank-sample {rank} kept '
                             f'{got.shape[0]} genomes, expected {n}')
            self.assertIn(rank, self.p.sample_note)
        self.assertEqual(
            sorted(self.p.sample(self.df, {'rank_sample': 'domain'})['genome']),
            ['ARC_1', 'BAC_1'],
            'each group must contribute its highest-quality genome')

    def test_yes_no_values_reach_the_build_as_booleans(self):
        from types import SimpleNamespace
        from neoHGT import database
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        src = join(d, 'src')
        os.makedirs(src)
        cases = [('--capital', 'capital', 'taxonomic name style'),
                 ('--latin', 'latin', 'taxonomic name style'),
                 ('--keep-placeholder-names', 'keep_placeholder_names',
                  'encoded GTDB labels'),
                 ('--input-fingerprints', 'fingerprints',
                  'input fingerprint table'),
                 ('--verify-md5', 'verify_md5', 'remote md5 checking'),
                 ('--allow-network', 'allow_network', 'downloads')]
        for flag, attr, what in cases:
            for value, want in (('yes', True), ('no', False)):
                args = cli_args('database', ['-o', join(d, 'db'),
                                             '--source', 'gtdb',
                                             '--gtdb-dir', src, flag, value])
                # the same chain a real run goes through: Database instance,
                # set_parameters (which converts yes/no), then the GTDB one
                me = database.Database.__new__(database.Database)
                me.cfg = database.load_configs()
                me.set_parameters(args)
                database.set_source_parameters(me, args)
                got = getattr(me, attr, None)
                self.assertIs(got, want,
                              f'{flag} {value} must turn {what} '
                              f'{"on" if want else "off"}, got {got!r}')

    def test_stream_extraction_says_what_the_gtdb_route_does(self):
        from neoHGT.database import Database
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        src = join(d, 'src')
        os.makedirs(src)
        out = io.StringIO()
        args = cli_args('database', ['-o', join(d, 'db'), '--source', 'gtdb',
                                     '--gtdb-dir', src,
                                     '--stream-extraction', 'no'])
        with contextlib.redirect_stdout(out):
            try:
                Database.__call__(Database.__new__(Database), args)
            except Exception:
                pass                       # empty dataset: fails further along
        self.assertIn('always streams', out.getvalue(),
                      'a GTDB build must not quietly ignore the request')

    def test_stream_extraction_default_is_not_treated_as_a_request(self):
        from neoHGT.database import Database
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        src = join(d, 'src')
        os.makedirs(src)
        out = io.StringIO()
        args = cli_args('database', ['-o', join(d, 'db'), '--source', 'gtdb',
                                     '--gtdb-dir', src])
        self.assertIsNone(vars(args)['stream_extraction'],
                          'no default may look like an explicit choice')
        with contextlib.redirect_stdout(out):
            try:
                Database.__call__(Database.__new__(Database), args)
            except Exception:
                pass
        self.assertNotIn('always streams', out.getvalue())


class TestFixtureDatabaseBuild(unittest.TestCase):
    """Build a real GTDB database from the committed fixtures, offline.

    Everything the build path does -- `load_genomes`, `protein_sources`,
    `write_gtdb_outputs`, `write_meta_json`, the derived `gtdb/` tables and the
    `map-report` subcommand -- had only ever run against the read-only dataset,
    so the manifest's reconciliation fields and the resume bookkeeping had no
    offline test at all. The build is done once for the class: it takes ~14 s.
    """

    @classmethod
    def setUpClass(cls):
        from neoHGT.database import Database
        cls.d = tmpdir()
        # the dataset and the output must be siblings: writing a database
        # inside the read-only input root is refused on purpose
        src = join(cls.d, 'src')
        for sub in ('tree_bac120', 'tree_ar53', 'metadata_bac120',
                    'metadata_ar53',
                    '00_GTDB_r220_aa/protein_faa_reps/bacteria',
                    '00_GTDB_r220_aa/protein_faa_reps/archaea'):
            os.makedirs(join(src, sub), exist_ok=True)
        renamed_tree(src)
        cls.out = join(cls.d, 'db')
        scratch = join(cls.d, 'tmp')
        os.makedirs(scratch)
        diamond = os.environ.get('NEOHGT_DIAMOND')
        cls.have_diamond = bool(diamond) and exists(diamond)
        argv = ['-o', cls.out, '--source', 'gtdb', '--gtdb-dir', src,
                '--reps', '--seq-source', 'local', '--tmpdir', scratch]
        argv += (['--compile', 'diamond', '--diamond', diamond]
                 if cls.have_diamond else ['--compile', 'none'])
        Database()(cli_args('database', argv))
        cls.meta = json.load(open(join(cls.out, 'meta.json')))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, ignore_errors=True)

    def test_build_outputs_and_manifest(self):
        for rel in ('db.faa', 'taxdump/nodes.dmp', 'taxdump/names.dmp',
                    'genomes.tsv', 'lineages.txt',
                    'meta.json', 'gtdb/taxid.map.gz', 'gtdb/labels.tsv',
                    'gtdb/stem2taxid.tsv.gz', 'gtdb/ncbi_map.tsv.gz'):
            self.assertTrue(isfile(join(self.out, rel)), rel)
        if self.have_diamond:
            # the compressed maps are published by the compile step; without a
            # compiler on PATH (--compile none) they legitimately do not exist
            for rel in ('taxon.map.gz', 'diamond/db.dmnd', 'genome.map.gz'):
                self.assertTrue(isfile(join(self.out, rel)), rel)
        c = self.meta['counts']
        self.assertGreater(c['proteins'], 0)
        self.assertEqual(c['genomes_missing'], 0)
        self.assertFalse(self.meta['extraction_resumed'])
        self.assertEqual(self.meta['source'], 'gtdb')
        # the fixtures are KB-scale *excerpts*, so the metadata promises far
        # more proteins than db.faa holds. That gap must be reported, never
        # smoothed over: this is the same reconciliation that watches a real
        # 13M-protein build (round 23 found two bugs in it).
        n_faa = sum(1 for l in open(join(self.out, 'db.faa'))
                    if l.startswith('>'))
        self.assertEqual(c['proteins'], n_faa)
        self.assertEqual(c['proteins_planned'],
                         c['representatives_in_dataset']['proteins'])
        self.assertGreater(c['proteins_planned'], n_faa,
                           'the fixture excerpts must be shorter than promised')
        self.assertEqual(
            c['proteins_extracted_vs_planned'],
            f'{c["proteins_planned"] - n_faa:,} short',
            f'a truncated extraction must be named exactly: {c}')
        report = open(join(self.out, 'gtdb', 'build_report.md')).read()
        self.assertIn('quality classes in the dataset', report,
                      'plan §7 asks the build report to state the actual '
                      'class counts')

    def test_map_report_subcommand_runs_offline(self):
        """`gtdb map-report` was only ever executed against the real dataset."""
        labels = [l.rstrip('\n').split('\t')
                  for l in open(join(self.out, 'gtdb', 'labels.tsv'))
                  if l.strip() and not l.startswith('#')]
        header, rows = labels[0], labels[1:]
        self.assertEqual(header[:2], ['taxid', 'label'], header)
        lst = join(self.d, 'hgts.txt')
        with open(lst, 'w') as f:
            f.write('\n'.join(r[0] for r in rows[:5]) + '\n')
        out = join(self.d, 'mapped.tsv')
        rc = gtdb.main(['map-report', '--db', self.out, '--in', lst,
                        '--out', out])
        self.assertEqual(rc, 0)
        rows = [l.rstrip('\n').split('\t') for l in open(out) if l.strip()]
        head, body = rows[0], rows[1:]
        self.assertEqual(head, ['donor_taxid', 'donor_label', 'donor_rank',
                                'ncbi_taxids', 'nconflict', 'status'], head)
        self.assertEqual(len(body), 5, body)
        self.assertTrue(all(r[3] and r[3] != 'na' for r in body),
                        'every donor of the fixture tree maps to NCBI taxIds')
        self.assertTrue(all(r[5] == 'ok' for r in body), body)


class TestNoDeadHelpers(unittest.TestCase):
    """Nothing in the GTDB branch may be defined and referenced nowhere.

    Round 28 found four such helpers in `neoHGT/gtdb.py` -- one of them
    (`write_lineages_txt`) was cited in the task register as the implementation
    of a deliverable while the file was actually written elsewhere, which is how
    dead code survives a path-checking gate.
    """

    # entry points are called by the console script or by the upstream CLI,
    # and these two are reached through polymorphism (`getattr(provider, ...)`)
    ENTRY = {'main', 'cli_main'}

    def _blob(self):
        code = join(HERE, os.pardir, os.pardir)
        parts = []
        for root, dirs, files in os.walk(code):
            dirs[:] = [d for d in dirs if d not in ('__pycache__',)]
            for f in files:
                if f.endswith('.py') or f.endswith('.yml'):
                    parts.append(open(os.path.join(root, f),
                                      encoding='utf-8',
                                      errors='replace').read())
        return '\n'.join(parts)

    def test_every_public_name_in_gtdb_is_referenced(self):
        """Module-level functions, classes and methods alike.

        A method counts as referenced only if something writes `.name(` or
        `name=` somewhere in the shipped code, tools or CI -- the same shape a
        real call has. Overriding an upstream attribute name is checked too:
        the wrapper must be reachable, not just plausibly named.
        """
        path = join(HERE, os.pardir, 'gtdb.py')
        tree = ast.parse(open(path, encoding='utf-8').read())
        blob = self._blob()
        dead = []

        def check(name, kind):
            if name.startswith('_') or name in self.ENTRY:
                return
            if blob.count(name) <= 1:
                dead.append(f'{kind} {name}')

        for node in tree.body:
            if isinstance(node, ast.ClassDef):
                check(node.name, 'class')
                for m in node.body:
                    if isinstance(m, (ast.FunctionDef, ast.AsyncFunctionDef)):
                        check(m.name, f'{node.name}.')
            elif isinstance(node, ast.FunctionDef):
                check(node.name, 'function')
        self.assertEqual(dead, [],
                         'defined in gtdb.py but referenced nowhere (dead code, '
                         'or a stale register citation)')


class TestPlanOptionWiring(unittest.TestCase):
    """Every option in plan §7 has to be parseable *and* have an effect.

    The register gate proves a flag exists in the CLI and is read somewhere in
    the code. It cannot see the opposite failure: a flag that parses into the
    wrong dest, or one whose parsed value never reaches the selection. Fifteen
    §7 options were not named in a single test before this class existed, so
    the parse check is table-driven and the behaviour check runs each selection
    option against a synthetic genome table.
    """

    FLAGS = [
        ('--complete', 'complete', True),
        ('--typemater', 'typemater', True),
        ('--reps', 'reps', True),
        ('--genome-list', 'genome_list', 'ids.IDlist'),
        ('--mimag', 'mimag', 'high'),
        ('--type-material', 'type_material', 'gtdb'),
        ('--min-completeness', 'min_completeness', 90.0),
        ('--max-contamination', 'max_contamination', 5.0),
        ('--keep-placeholder-names', 'keep_placeholder_names', 'no'),
        ('--gtdb-taxa', 'gtdb_taxa', 'f__Enterobacteriaceae'),
        ('--gtdb-release', 'gtdb_release', '226'),
        ('--local-layout', 'local_layout', 'official-gz'),
        ('--gtdb-taxid-base', 'gtdb_taxid_base', 200000000),
        ('--gtdb-id-manifest', 'gtdb_id_manifest', 'm.map.gz'),
        ('--seq-source', 'seq_source', 'gtdb-tar'),
    ]

    def test_plan_options_reach_the_namespace(self):
        for flag, dest, want in self.FLAGS:
            argv = [flag] if want is True else [flag, str(want)]
            args = cli_args('database', ['-o', 'db', '--source', 'gtdb',
                                         '--gtdb-dir', 'x'] + argv)
            d = vars(args)
            self.assertIn(dest, d, f'{flag} parses to no usable dest')
            self.assertEqual(d[dest], want, f'{flag} -> {d[dest]!r}')

    LIN = ('d__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G1;s__G1 alpha')

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        rows = []
        for i, (rep, comp, cont, mag, typ) in enumerate([
                (True, 99.0, 0.5, 't', 't'),
                (True, 90.0, 5.0, 'f', 'f'),
                (False, 60.0, 30.0, 'f', 'f')], 1):
            rows.append({'genome': f'BAC_GCA_{i:09d}',
                         'accession': f'GB_GCA_{i:09d}.1',
                         'accnov': f'GCA_{i:09d}',
                         'gtdb_taxonomy': self.LIN, 'domain': 'bacteria',
                         'is_rep': rep, 'proteins': 100,
                         'completeness': comp, 'contamination': cont,
                         'quality_score': comp - 5 * cont,
                         'mimag_high_quality': mag,
                         'mimag_medium_quality': 't' if comp >= 90 else 'f',
                         'gtdb_type_species_of_genus': typ})
        self.df = pd.DataFrame(rows)
        p = gtdb.GtdbProvider.__new__(gtdb.GtdbProvider)
        p.taxid_map, _ = gtdb.assign_ids(gtdb.collect_nodes([self.LIN]))
        p.opts, p.selection_labels, p.sample_note = {}, [], None
        self.p = p

    def _sel(self, opts):
        self.p.selection_labels = []
        return self.p.select(self.df, opts)

    def test_selection_options_filter_the_genome_table(self):
        self.assertEqual(self._sel({'reps': True}).shape[0], 2)
        self.assertEqual(self.p.selection_labels, ['representatives'])
        self.assertEqual(self._sel({'min_completeness': 95.0}).shape[0], 1)
        self.assertEqual(self._sel({'max_contamination': 1.0}).shape[0], 1)
        self.assertEqual(self._sel({'mimag': 'high'}).shape[0], 1)
        self.assertEqual(self.p.selection_labels, ['MIMAG-high'])
        self.assertEqual(self._sel({'mimag': 'medium'}).shape[0], 2)
        self.assertEqual(self._sel({'type_material': 'gtdb'}).shape[0], 1)
        # plan §7: the type-material filters must say how many genomes they
        # actually kept, so the label carries both sides of the cut.
        self.assertEqual(self.p.selection_labels,
                         ['GTDB type species (3 -> 1 genomes)'])
        lst = join(self.d, 'ids.IDlist')
        with open(lst, 'w') as f:
            f.write('GB_GCA_000000002.1\nGCA_000000003\n')
        sel = self._sel({'genome_list': lst})
        self.assertEqual(sel.shape[0], 2,
                         'both a full and an unversioned accession must work')
        self.assertEqual(self.p.selection_labels, ['genome-list(2)'])

    def test_parsed_options_reach_the_build_opts(self):
        # The build path used to hand-build its provider opts dict inline, and
        # two flags that parse perfectly (--complete, --typemater) were not in
        # it: from the command line the selector never saw them at all, so the
        # round-29 alias and refusal would have been unreachable code.
        for flag, dest, want in (('--complete', 'complete', True),
                                 ('--typemater', 'typemater', True),
                                 ('--mimag', 'mimag', 'high'),
                                 ('--type-material', 'type_material', 'ncbi'),
                                 ('--min-completeness', 'min_completeness', 95.0),
                                 ('--max-contamination', 'max_contamination', 1.0),
                                 ('--reps', 'reps', True)):
            argv = [flag] if want is True else [flag, str(want)]
            args = cli_args('database', ['-o', 'db', '--source', 'gtdb',
                                         '--gtdb-dir', 'x'] + argv)
            self.assertEqual(gtdb.selection_opts(args).get(dest), want,
                             f'{flag} parses but never reaches the opts')
        with open(join(HERE, os.pardir, 'database.py'),
                  encoding='utf-8') as fh:
            self.assertIn('gtdb.selection_opts(args)', fh.read(),
                          'the build path must install the very opts the '
                          'selector is tested against')
        with self.assertRaises(ValueError):
            self._sel(gtdb.selection_opts(cli_args('database', [
                '-o', 'db', '--source', 'gtdb', '--gtdb-dir', 'x',
                '--complete'])))

    def test_placeholder_names_option_reaches_the_taxonomy_filters(self):
        args = cli_args('search', ['-i', 'x.tsv', '-t', 'taxdump',
                                   '-o', 'out', '--source', 'gtdb',
                                   '--keep-placeholder-names', 'no'])
        self.assertEqual(vars(args)['keep_placeholder_names'], 'no')
        args = cli_args('database', ['-o', 'db', '--source', 'gtdb',
                                     '--gtdb-dir', 'x'])
        self.assertNotEqual(vars(args).get('keep_placeholder_names'), 'no',
                            'the default must not silently drop placeholders')

    def test_gtdb_options_reach_the_converted_tree(self):
        # --gtdb-release / --gtdb-taxid-base / --gtdb-id-manifest / --local-layout
        # must show up in what the CLI actually writes, not just in argv
        src = join(self.d, 'flat')
        os.makedirs(src)
        body = '\n'.join(f'GB_GCA_{i:09d}.1\t{self.LIN}'
                         for i in (1, 2)) + '\n'
        for fam in ('bac120', 'ar53'):
            with open(join(src, f'{fam}_taxonomy_r999.tsv'), 'w') as f:
                f.write(body)
        out = join(self.d, 'out')
        rc = gtdb.main(['taxdump', '--dir', src, '--out', out,
                        '--gtdb-release', '999', '--local-layout',
                        'official-plain', '--gtdb-taxid-base', '200000000'])
        self.assertEqual(rc, 0)
        meta = json.load(open(join(out, 'meta.json')))
        self.assertEqual(meta['gtdb_release'], '999')
        self.assertEqual(meta['local_layout'], 'official-plain')
        self.assertEqual(meta['id_base'], 200000000)
        man = join(out, 'gtdb', 'taxid.map.gz')
        ids = sorted(int(l.split('\t')[1]) for l in
                     gzip.open(man, 'rt').read().splitlines()
                     if l and not l.startswith('#')
                     and l.split('\t')[0] != 'lineage')
        self.assertTrue(min(ids) >= 200000000,
                        f'--gtdb-taxid-base ignored: min id {min(ids)}')
        # and the manifest is honoured on a second pass: no id changes
        out2 = join(self.d, 'out2')
        gtdb.main(['taxdump', '--dir', src, '--out', out2,
                   '--gtdb-id-manifest', man, '--gtdb-taxid-base',
                   '200000000'])
        ids2 = sorted(int(l.split('\t')[1]) for l in
                      gzip.open(join(out2, 'gtdb', 'taxid.map.gz'),
                                'rt').read().splitlines()
                      if l and not l.startswith('#')
                      and l.split('\t')[0] != 'lineage')
        self.assertEqual(ids, ids2)


class TestSelectionAndSampling(unittest.TestCase):
    """UT for plan §6.4 tasks 2.5/2.13: lineage selection and -r/-s/--above.

    The genome table is synthesized so the rank geometry is exactly known; the
    read-only dataset is never touched.
    """

    LINEAGES = (
        ('d__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G1;s__G1 alpha', 3),
        ('d__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G1;s__G1 beta', 1),
        ('d__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G2', 2),   # stops at genus
        ('d__Bacteria;p__P2;c__C2;o__O2;f__F2;g__G3;s__G3 alpha', 2))

    def setUp(self):
        rows, i = [], 0
        for lin, n in self.LINEAGES:
            for j in range(n):
                i += 1
                rows.append({'genome': f'BAC_GCA_{i:09d}',
                             'accession': f'GB_GCA_{i:09d}.1',
                             'accnov': f'GCA_{i:09d}',
                             'gtdb_taxonomy': lin, 'domain': 'bacteria',
                             'is_rep': True, 'proteins': 1000,
                             'completeness': 90, 'contamination': 1,
                             # earlier rows of a group rank higher
                             'quality_score': 1000 - i})
        self.df = pd.DataFrame(rows)
        # a provider without any dataset IO: select/sample only need these
        p = gtdb.GtdbProvider.__new__(gtdb.GtdbProvider)
        p.taxid_map, _ = gtdb.assign_ids(
            gtdb.collect_nodes([l for l, _ in self.LINEAGES]))
        p.opts, p.selection_labels, p.sample_note = {}, [], None
        self.p = p

    def test_lineage_label_selection(self):
        sel = self.p.select(self.df, {'gtdb_taxa': ['f__F1']})
        self.assertEqual(sel.shape[0], 6)          # 3 + 1 + 2 genomes
        self.assertEqual(self.p.selection_labels, ['lineage(1 taxa)'])
        self.assertTrue(sel['gtdb_taxonomy'].str.contains('f__F1').all())

    def test_taxid_selection_matches_label(self):
        by_label = self.p.select(self.df, {'gtdb_taxa': ['f__F1']})
        key = self.p.node_key_of('f__F1')
        by_tid = self.p.select(self.df, {'taxids': [self.p.taxid_map[key]]})
        self.assertEqual(set(by_label['genome']), set(by_tid['genome']))

    def test_sample_excludes_lineages_shallower_than_the_rank(self):
        res = self.p.sample(self.df, {'rank': 'species', 'sample': 2})
        # 2 of 'G1 alpha', 1 of 'G1 beta', 2 of 'G3 alpha'; the two genomes that
        # GTDB only identifies to genus are left out, as in the NCBI route
        self.assertEqual(res.shape[0], 5)
        self.assertEqual(self.p.sample_note, '2 genome(s) per species')
        self.assertNotIn('_key', res.columns)
        self.assertNotIn('g__G2', ''.join(res['gtdb_taxonomy']))

    def test_above_tops_up_coarser_groups(self):
        base = set(self.p.sample(
            self.df, {'rank': 'species', 'sample': 2})['genome'])
        res = self.p.sample(self.df, {'rank': 'species', 'sample': 2,
                                      'above': True})
        self.assertEqual(res.shape[0], 7)
        self.assertTrue(base < set(res['genome']))
        self.assertIn('g__G2', ''.join(res['gtdb_taxonomy']))
        self.assertIn('and above', self.p.sample_note)
        self.assertNotIn('_key', res.columns)
        self.assertEqual(res['genome'].duplicated().sum(), 0)

    def test_rank_sample_keeps_best_per_group(self):
        res = self.p.sample(self.df, {'rank_sample': 'family'})
        self.assertEqual(res.shape[0], 2)
        # the highest-quality genome of each family
        self.assertEqual(list(res.sort_values('gtdb_taxonomy')['genome']),
                         ['BAC_GCA_000000001', 'BAC_GCA_000000007'])


class TestTaxonSynonyms(unittest.TestCase):
    """RK11: GTDB retires higher-taxon names, and the local copy has no
    synonyms.*.tsv, so the alias table must come from the user and be honoured
    by every taxon reference in analyze."""

    def setUp(self):
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)

    def _write(self, text):
        path = join(self.d, 'syn.tsv')
        with open(path, 'w') as f:
            f.write(text)
        return path

    def test_two_column_table_with_note_column(self):
        path = self._write('# comment\n'
                           'p__Proteobacteria\tp__Pseudomonad\trenamed r220\n'
                           'c__Bacilli\tc__SomeOther\t\n')
        self.assertEqual(
            gtdb.load_taxon_synonyms(path),
            {'p__Proteobacteria': 'p__Pseudomonad',
             'c__Bacilli': 'c__SomeOther'})

    def test_missing_or_malformed_is_an_explicit_error(self):
        self.assertEqual(gtdb.load_taxon_synonyms(None), {})
        with self.assertRaises(ValueError) as ctx:
            gtdb.load_taxon_synonyms(join(self.d, 'nope.tsv'))
        self.assertIn('--taxon-synonyms', str(ctx.exception))
        with self.assertRaises(ValueError):
            gtdb.load_taxon_synonyms(self._write('only_one_column\n'))

    def test_analyze_declares_the_option(self):
        from neoHGT import analyze as analyze_mod
        flags = [a[0] for a in analyze_mod.arguments if isinstance(a, (list, tuple))]
        self.assertTrue(any('--taxon-synonyms' in f for f in flags), flags)


class TestProductionReachability(unittest.TestCase):
    """Every function this adaptation added must be reachable from a CLI path.

    Coverage says "a test called this"; the product question is "can the
    command line get there". `--hybrid-lca-policy` is the case that motivated
    the check: its resolver was validated, echoed into the build report and
    called by tests, while no build path ever invoked it.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_production_reachability.py')
        spec = spec_from_file_location('vpr', path,
                                        loader=SourceFileLoader('vpr', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_no_added_function_is_out_of_reach_without_a_reason(self):
        doc = join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')
        if not isfile(doc):
            # the internal QC registers (this hub deliverables document among
            # them) live in the project workspace, not in the public release;
            # the declared-exception cross-check needs that document
            self.skipTest('internal QC registers not part of this checkout')
        tool = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        base = os.environ.get('NEOHGT_MAIN') or join(code, os.pardir,
                                                    'neoHGT-main')
        if isdir(base):
            # strict mode: fresh diff against the pristine upstream tree
            unreachable, reach, added = tool.analyse(code, base)
        else:
            # a checkout without the upstream tree (CI) cannot tell added
            # functions from inherited ones by diffing; judge the frozen
            # added-functions list instead of every definition, so inherited
            # upstream functions cannot look like additions
            added_names, _upstream = tool.frozen_lists(
                join(code, 'docs', '新增函数清单.json'))
            unreachable, reach, added = tool.analyse(code, None,
                                                     sorted(added_names))
        self.assertGreater(added, 120,
                           'the audit found almost no added functions, so it '
                           'cannot conclude anything (wrong --code path?)')
        declared = tool.declared_exceptions(
            join(code, 'docs', '交付物清单.md'))
        undeclared = {rel: [n for n in names if n not in declared
                            and n.split('.')[-1] not in declared]
                      for rel, names in unreachable.items()}
        self.assertEqual({k: v for k, v in undeclared.items() if v}, {},
                         f'{added} added functions, {reach} reachable names, '
                         f'but these are unreachable and undeclared')

    def test_the_detector_can_see_the_failure_it_claims_to_catch(self):
        # rewriting the round-30 call site into a dummy lambda must make the
        # resolver unreachable again; a control that cannot fail proves nothing
        tool = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        nc = tool.control(code)
        self.assertTrue(nc.get('applied'), nc.get('note'))
        self.assertTrue(nc.get('resolver_becomes_unreachable'),
                        'the call graph cannot see a deleted call site')


    def test_a_wrapper_reached_only_through_a_rebinding_is_reachable(self):
        # the round-40 tool fix, with its own two-way control:
        # `Class.slot = wrapper` and `from x import name as alias` are dispatch
        # edges. A walk that ignores them calls a live function unreachable --
        # which is how resolve_taxon_ref() came to be wrongly declared.
        tool = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'scripts'))
        with open(join(d, 'scripts', 'neoHGT'), 'w') as fh:
            fh.write("module_names = ('gtdb',)\n")
        os.makedirs(join(d, 'neoHGT'))

        def write(with_rebinding):
            body = ("class Gtdb(object):\n"
                    "    def __call__(self, args):\n"
                    "        self.run(args)\n"
                    "\n"
                    "    def run(self, args):\n"
                    "        return None\n"
                    "\n"
                    "def wrapper(self, args):\n"
                    "    return alias_helper(args)\n"
                    "\n"
                    "from neoHGT.util import resolve_taxon_ref as "
                    "alias_helper\n"
                    "\n"
                    "def resolve_taxon_ref(x):\n"
                    "    return x\n")
            if with_rebinding:
                body += "\nGtdb.run = wrapper\n"
            with open(join(d, 'neoHGT', 'gtdb.py'), 'w') as fh:
                fh.write(body)

        verdicts = {}
        for keep in (True, False):
            write(keep)
            unreachable, _reach, added = tool.analyse(
                d, None, {'wrapper', 'resolve_taxon_ref'})
            self.assertEqual(added, 2,
                             'the audit saw no functions, so it proves nothing')
            verdicts[keep] = unreachable.get('neoHGT/gtdb.py', [])
        self.assertEqual(verdicts[True], [],
                         'a rebound slot and an import alias must count as '
                         'entry points')
        self.assertEqual(sorted(verdicts[False]),
                         ['resolve_taxon_ref', 'wrapper'],
                         'without those edges nothing reaches them -- if this '
                         'is empty the control proves nothing')

    def test_the_frozen_list_covers_every_definition_today(self):
        # the hole: CI judges the intersection of the list with the tree, so a
        # function added after the list was frozen belonged to no set at all
        tool = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        defs = tool.defs_by_module(code)
        added, upstream = tool.frozen_lists(join(code, 'docs',
                                                 '新增函数清单.json'))
        self.assertTrue(added and upstream,
                        'the committed list is still the one-half format')
        self.assertEqual(added & upstream, set())
        self.assertEqual(tool.untracked_names(defs, added, upstream), {})
        for rel, names in defs.items():
            for n in names:
                self.assertTrue(n in added or n in upstream
                                or n.split('.')[-1] in added
                                or n.split('.')[-1] in upstream,
                                f'{rel}:{n} is in neither list')

    def test_an_untracked_function_is_named(self):
        tool = self._tool()
        defs = {'neoHGT/gtdb.py': {'known_added': 1, 'known_upstream': 2,
                                   'sneaked_in': 3}}
        self.assertEqual(tool.untracked_names(defs, {'known_added'},
                                              {'known_upstream'}),
                         {'neoHGT/gtdb.py': ['sneaked_in']},
                         'a function in neither list passed unnoticed')
        # the old one-half format cannot support this question at all: without
        # the upstream names, every inherited function looks "untracked", so the
        # gate had to skip the check -- which is how a later addition escaped it
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        legacy = join(d, 'legacy.json')
        with open(legacy, 'w') as fh:
            fh.write('{"functions": ["known_added"]}')
        a2, u2 = tool.frozen_lists(legacy)
        self.assertEqual((a2, u2), ({'known_added'}, set()))
        self.assertEqual(sorted(tool.untracked_names(defs, a2, u2)
                                 ['neoHGT/gtdb.py']),
                         ['known_upstream', 'sneaked_in'],
                         'without the upstream half even inherited code looks '
                         'new, which is why main() skipped the check instead')

    def test_the_list_is_compared_with_a_fresh_diff_locally(self):
        tool = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        for root, body in ((join(d, 'code'),
                            'def kept(self):\n    pass\n\n'
                            'def added_later(self):\n    pass\n'),
                           (join(d, 'base'),
                            'def kept(self):\n    pass\n\n'
                            'def upstream_gone(self):\n    pass\n')):
            os.makedirs(join(root, 'neoHGT'), exist_ok=True)
            with open(join(root, 'neoHGT', 'gtdb.py'), 'w') as fh:
                fh.write(body)
        defs = tool.defs_by_module(join(d, 'code'))
        self.assertIn('neoHGT/gtdb.py', defs)
        missing, stale, upstream_changed = tool.list_drift(
            defs, join(d, 'base'), {'kept', 'vanished'}, {'kept'})
        self.assertEqual(missing, ['added_later'],
                         'a function added after the freeze must show up')
        self.assertEqual(stale, ['vanished'],
                         'a name no longer defined anywhere must be dropped')
        self.assertEqual(upstream_changed, ['neoHGT/gtdb.py:upstream_gone'],
                         upstream_changed)


class TestCheckerToolUpkeep(unittest.TestCase):
    """A checker that nothing runs is evidence that can no longer be reproduced.

    Four `tools/verify_*.py` are cited in the docs as the commands that produced
    committed evidence but need things outside the repository (the read-only
    dataset, the pristine upstream tree, one-shot snapshots), so they cannot go
    into CI. They are declared instead -- and a declaration has to be paid for:
    each one must still at least *run*, or the documented command is a lie.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt5', path,
                                       loader=SourceFileLoader('vpt5', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def setUp(self):
        self.code = join(HERE, os.pardir, os.pardir)
        v = self._tool()
        self.checkers = v.checker_tools(self.code)
        self.problems, self.declared = v.checker_tool_problems(
            self.checkers, join(self.code, 'docs', '交付物清单.md'))

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_every_runner_less_tool_is_declared(self):
        unrun = sorted(k for k, who in self.checkers.items() if not who)
        self.assertEqual(self.problems, [], self.checkers)
        self.assertEqual(unrun, sorted(self.declared),
                         (unrun, self.declared))
        self.assertGreaterEqual(len(unrun), 4,
                                'the scan found almost no checker tools')

    def test_a_declared_tool_can_still_be_invoked(self):
        # this is what makes the declaration honest rather than a way of
        # forgetting about a tool: the documented command has to at least start
        for name in self.declared:
            path = join(self.code, 'tools', name)
            self.assertTrue(isfile(path), name)
            rc = subprocess.call([sys.executable, path, '--help'],
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL)
            self.assertEqual(rc, 0, f'{name} --help failed; the declared '
                                    'command can no longer be run')

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_the_rule_bites_on_an_undeclared_orphan(self):
        v = self._tool()
        problems, declared = v.checker_tool_problems(
            {'verify_ghost.py': [], 'verify_real.py': ['neoHGT/tests/x.py']},
            join(self.code, 'docs', '交付物清单.md'))
        self.assertEqual(declared, sorted(self.declared))
        self.assertEqual(len(problems), 1, problems)
        self.assertIn('verify_ghost.py', problems[0])

    def test_the_nightly_reuses_the_read_only_checker(self):
        # the nightly workflow used to carry its own `find -newermt` copy of this
        # check: two implementations of one guard, with the weaker one in CI
        path = join(self.code, '.github', 'workflows', 'nightly.yml')
        text = open(path, encoding='utf-8').read()
        self.assertIn('tools/verify_read_only.py', text)
        self.assertNotIn('find "$GTDB_DIR" -type f -newermt', text)


class TestAssertionFreeRule(unittest.TestCase):
    """A test that asserts nothing is a green light with no information in it.

    Round 37 found two such tests in this suite: they called the guard and never
    looked at the result, so a guard that had been turned into `pass` would have
    kept them green. The rule is now machine-checked, and the upstream files --
    which RT-1 keeps append-only -- are registered rather than edited.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt4', path,
                                       loader=SourceFileLoader('vpt4', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_the_shipped_suite_declares_every_vacuous_test(self):
        v = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        found = v.assertion_free_tests(code)
        declared = v.declared_vacuous_tests(join(code, 'docs',
                                                 '交付物清单.md'))
        self.assertEqual([x for x in found if x not in declared], [],
                         'a test with no observable assertion is not registered')
        self.assertEqual(sorted(declared - set(found)), [],
                         'the register lists tests that now assert something')
        self.assertGreaterEqual(len(found), 12,
                                'the scan found almost nothing (wrong path?)')
        self.assertFalse([x for x in found
                          if x.startswith(('test_gtdb_pipeline.py',
                                           'test_gtdb_core.py',
                                           'test_util_gtdb.py'))],
                         'this adaptation must not add assertion-free tests')

    def test_the_scanner_bites_on_a_synthetic_vacuous_test(self):
        v = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'neoHGT', 'tests'))
        with open(join(d, 'neoHGT', 'tests', 'test_fake.py'), 'w') as fh:
            fh.write('import unittest\n\n\nclass T(unittest.TestCase):\n'
                     '    def test_only_calls(self):\n        self.helper()\n'
                     '\n    def test_asserts(self):\n        '
                     'self.assertEqual(1, 1)\n'
                     '\n    def test_expects_a_raise(self):\n        '
                     'with self.assertRaises(ValueError):\n            '
                     'self.helper()\n')
        self.assertEqual(v.assertion_free_tests(d),
                         ['test_fake.py::test_only_calls'])


class TestCliSurfaceAndDocExamples(unittest.TestCase):
    """Our own command examples are the thing users copy, so they get a gate.

    The register gate already checked the *plan's* appendix C/D against the CLI,
    but never read README.md or docs/*.md -- and it read the `neoHGT gtdb` option
    surface out of the source text with an `add_argument('...')` regex, which
    cannot see a table entry at all: `--action`, `--deep`, `--threads`,
    `--tmpdir` and `--output` were invisible to every rule built on that set.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt3', path,
                                       loader=SourceFileLoader('vpt3', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def setUp(self):
        self.code = join(HERE, os.pardir, os.pardir)
        self.tool = self._tool()
        self.surfaces = self.tool.cli_flags(self.code)
        self.subs = self.tool.cli_subcommands(self.code)
        self.root = tmpdir()
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        for name in ('database', 'search', 'analyze', 'gtdb'):
            os.makedirs(join(self.root, 'neoHGT'), exist_ok=True)
            with open(join(self.root, 'neoHGT', f'{name}.py'), 'w') as f:
                f.write("arguments = [\n    'basic',\n    ['--action',\n"
                        "     'written on its own line',\n"
                        "     {'choices': ['probe']}],\n]\n")
        os.makedirs(join(self.root, 'scripts'), exist_ok=True)
        with open(join(self.root, 'scripts', 'neoHGT'), 'w') as f:
            f.write("module_names = ('search', 'analyze', 'database', 'gtdb')\n")
        os.makedirs(join(self.root, 'docs'), exist_ok=True)

    def _docs(self, *lines):
        with open(join(self.root, 'README.md'), 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines) + '\n')
        return self.root

    def test_the_flat_table_entries_are_part_of_the_surface(self):
        # the bug this round fixed: gtdb's options live in an `arguments` table,
        # and a regex over add_argument calls structurally cannot see them
        for flag in ('--action', '--deep', '--threads', '--tmpdir', '--output'):
            self.assertIn(flag, self.surfaces['gtdb'], flag)
        probe = self.tool.cli_flags(self.root)
        self.assertIn('--action', probe['gtdb'],
                      'a table entry with the option string on its own line '
                      'must still reach the surface')

    def test_an_example_naming_an_option_its_subcommand_lacks_is_reported(self):
        used = self.tool.doc_command_flags(
            self._docs('neoHGT analyze --mimag high -o out',
                       'conda remove -n neoHGT --all'), self.subs)
        surfaces = {'analyze': {'--output'}, 'database': {'--mimag'}}
        problems = self.tool.doc_flag_problems(used, surfaces)
        self.assertTrue(problems, used)
        self.assertIn('--mimag', problems[0])
        self.assertNotIn('--all', used, 'a third-party command line is not ours')

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_the_deviation_row_is_excused_but_the_same_flag_elsewhere_is_not(self):
        ledger = join(self.code, 'docs', '交付物清单.md')
        declared = self.tool.declared_doc_usages(ledger, self.subs,
                                                 self.surfaces)
        self.assertTrue(declared, 'the deviation table quotes no foreign flag?')
        where, stage, flag = sorted(declared)[0]
        self.assertEqual(self.tool.doc_flag_problems(
            {flag: [where]}, self.surfaces, declared), [])
        self.assertTrue(self.tool.doc_flag_problems(
            {flag: ['elsewhere.md:1 (%s)' % stage]}, self.surfaces, declared),
            'the exemption has to stay scoped to the line that declares it')

    def test_the_shipped_docs_only_use_options_their_subcommand_has(self):
        used = self.tool.doc_command_flags(self.code, self.subs)
        scanned = sum(len(v) for v in used.values())
        self.assertGreater(scanned, 60,
                           'barely any example was read, so this proves '
                           'nothing (did the scanner lose the sub-command list?)')
        self.assertEqual(self.tool.doc_flag_problems(
            used, self.surfaces,
            self.tool.declared_doc_usages(join(self.code, 'docs',
                                               '交付物清单.md'),
                                          self.subs, self.surfaces)), [],
            f'{scanned} example invocations read')


class TestSelfDiscoveryRegister(unittest.TestCase):
    """The running self-discovery list must state what it actually holds.

    Its heading carries two hand-maintained numbers -- the round span and the
    item count -- and the list grows by one every round, so the heading is the
    first thing to go stale. Round 32 came within one edit of shipping "十五个"
    above sixteen entries; a reader then cannot tell whether an item vanished or
    the count did.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt2', path,
                                       loader=SourceFileLoader('vpt2', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def setUp(self):
        self.doc = join(HERE, os.pardir, os.pardir, 'docs', 'GTDB-数据库.md')
        self.text = open(self.doc, encoding='utf-8').read()

    def _run(self, text, name):
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        path = join(d, name)
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text)
        return self._tool().self_discovery_check(path)

    def test_the_shipped_list_agrees_with_its_own_heading(self):
        # the anti-vacuous half: the real document parses, and parses clean
        res = self._tool().self_discovery_check(self.doc)
        self.assertEqual(res['problems'], [], res['heading'])
        self.assertGreaterEqual(res['items'], 15,
                                'the list was not read at all (wrong path?)')
        self.assertEqual(res['declared_count'], res['items'])
        self.assertEqual(res['rounds'][1], max(res['tagged_rounds']))

    @staticmethod
    def _cn_prev(numeral):
        """One less than a Chinese numeral up to 99, written the same way."""
        digits = '〇一二三四五六七八九'
        if numeral == '十':
            return '九'
        if numeral.endswith('十'):
            return digits[digits.index(numeral[0]) - 1] + '十'
        if '十' in numeral:
            head, _, tail = numeral.partition('十')
            return head + '十' + digits[digits.index(tail) - 1]
        return digits[digits.index(numeral) - 1]

    def _mutants(self):
        """The four drifts, spelled out of the document as it stands now.

        Nothing here may name a number or a round: the round-32 version of this
        test hard-coded "十六个" and "（第三十二轮）", and appending item 17 in
        round 33 turned that snapshot into a failing test. A control that breaks
        when the thing it guards grows is guarding the wrong thing.
        """
        v = self._tool()
        lines = self.text.splitlines()
        at = next((i for i, l in enumerate(lines)
                   if v.DISCOVERY_HEADING_RE.match(l)), None)
        self.assertIsNotNone(at, 'the list heading is no longer in the shape '
                                 'the rule reads, so this control has nothing '
                                 'to mutate')
        head = v.DISCOVERY_HEADING_RE.match(lines[at])
        line, start = lines[at], self.text.index(lines[at])
        # the rule's pattern is line-anchored without re.M because it is used
        # one line at a time; scanning the whole document needs the flag back
        # the section only: the next list in the file restarts at 1, and a scan
        # over the whole document picks up those numbers -- which made two of
        # these four "mutants" edit some other list and report nothing at all
        nxt = self.text.find('\n### ', start + len(line))
        end = nxt + 1 if nxt >= 0 else len(self.text)
        items = [int(m.group(1)) for m
                 in re.finditer(v.ITEM_START_RE.pattern,
                                self.text[start:end], re.M)]
        self.assertEqual(items, list(range(1, len(items) + 1)),
                         'this control assumes a contiguous list')
        last = self.text.index('\n%d. ' % items[-1], start)
        block = self.text[last:end]
        tag = v.ROUND_TAG_RE.search(block).group(1)
        return {
            # the heading left behind when an item is appended
            'stale_count': self.text.replace(
                line, line.replace('%s个' % head.group('count'),
                                   '%s个' % self._cn_prev(head.group('count'))),
                1),
            # an item pasted without renumbering
            'duplicate': self.text.replace('\n%d. ' % items[-1],
                                           '\n%d. ' % (items[-1] - 1), 1),
            # the newest entry does not say which round it came from
            'no_round_tag': (self.text[:last]
                             + block.replace('第%s轮' % tag, '', 1)
                             + self.text[end:]),
            # the range stops short of the newest entry
            'stale_span': self.text.replace(
                line, line.replace('至%s轮' % head.group('to'),
                                   '至%s轮' % self._cn_prev(head.group('to'))),
                1),
        }

    def test_the_counter_bites_on_every_drift_it_claims_to_catch(self):
        for name, text in self._mutants().items():
            self.assertNotEqual(text, self.text,
                                '%s mutated nothing: the document no longer '
                                'looks the way this control assumes' % name)
            res = self._run(text, '%s.md' % name)
            self.assertTrue(res['problems'],
                            '%s drifted but the rule stayed silent' % name)


class TestSearchChoiceValues(unittest.TestCase):
    """Every search-side choice value must filter, refuse, or explain itself.

    The probe in `缓存测试/reports/R31_search_value_probe.json` ran the clean
    fixture hit table through all of these values and the table came out
    unchanged every time -- which proves nothing about the switches, only that
    the fixture has nothing to filter. So the values are driven against the
    naming reality GTDB actually has: coded placeholder labels, a taxId that
    repeats within one query, and a hit table without the staxids column.
    """

    def _stub(self, **kw):
        from neoHGT.search import Search
        from neoHGT.util import arg2bool
        s = Search.__new__(Search)
        s.taxdump = {
            '0': {'name': 'root', 'parent': '0', 'rank': 'no rank'},
            '1': {'name': 'Escherichia coli', 'parent': '0', 'rank': 'species'},
            '2': {'name': '2-12-FULL-45-22', 'parent': '0', 'rank': 'species'},
        }
        s.badtaxids, s.tax_include, s.tax_exclude = set(), [], []
        s.tax_block, s.tax_unirank, s.gtdb_root = [], None, '/db'
        for name, val in dict(tax_capital=False, tax_latin=False,
                              tax_unique=False, keep_placeholder_names=True,
                              require_staxids=True).items():
            setattr(s, name, val)
        for name, val in kw.items():
            # the CLI spells these switches as the strings 'yes'/'no', and a
            # truthiness test on the string 'no' is true -- convert exactly the
            # way the pipeline does
            setattr(s, name, arg2bool(val) if isinstance(val, str) else val)
        return s

    def _hits(self):
        # three hits for one query: a proper species, a coded GTDB placeholder,
        # and a repeat of the first taxId
        return {'q1': [{'taxid': '1'}, {'taxid': '2'}, {'taxid': '1'}]}

    def test_capital_filter_keeps_placeholder_labels_only_when_exempted(self):
        for keep, want in (('yes', 3), ('no', 2)):
            s = self._stub(tax_capital='yes', keep_placeholder_names=keep)
            prots = self._hits()
            s.taxfilt_wf(prots)
            self.assertEqual(len(prots['q1']), want,
                             f'--keep-placeholder-names {keep} with '
                             '--tax-capital yes kept '
                             f'{len(prots["q1"])} hits')
        off = self._stub(tax_capital='no', keep_placeholder_names='no')
        prots = self._hits()
        off.taxfilt_wf(prots)
        self.assertEqual(len(prots['q1']), 3,
                         '--tax-capital no must drop nothing by case')

    def test_latin_filter_keeps_placeholder_labels_only_when_exempted(self):
        for keep, want in (('yes', 3), ('no', 2)):
            s = self._stub(tax_latin='yes', keep_placeholder_names=keep)
            prots = self._hits()
            s.taxfilt_wf(prots)
            self.assertEqual(len(prots['q1']), want,
                             f'--keep-placeholder-names {keep} with '
                             f'--tax-latin yes kept {len(prots["q1"])} hits')
        off = self._stub(tax_latin='no', keep_placeholder_names='no')
        prots = self._hits()
        off.taxfilt_wf(prots)
        self.assertEqual(len(prots['q1']), 3,
                         '--tax-latin no must not filter by Latinate names')

    def test_tax_unique_collapses_a_repeated_taxid(self):
        for unique, want in (('yes', 2), ('no', 3)):
            s = self._stub(tax_unique=unique)
            prots = self._hits()
            s.taxfilt_wf(prots)
            self.assertEqual(len(prots['q1']), want,
                             f'--tax-unique {unique} left '
                             f'{len(prots["q1"])} hits')

    def test_require_staxids_no_warns_instead_of_failing(self):
        rows = ['q1\ts1\t100\t1e-50\t200\t90']       # six columns, no staxids
        s = self._stub(require_staxids='yes')
        with self.assertRaises(ValueError) as ctx:
            s.parse_def_table(rows)
        self.assertIn('staxids', str(ctx.exception))
        loose = self._stub(require_staxids='no')
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            try:
                loose.parse_def_table(rows)
            except ValueError as exc:
                self.fail(f'--require-staxids no must not fail outright: {exc}')
        self.assertIn('WARNING', buf.getvalue(),
                      'waiving the requirement has to be visible in the log')

    def test_absent_flag_still_requires_staxids(self):
        # RK5's protection is the fail-fast behaviour, so "the user said
        # nothing" must mean "required": a stub that never got the attribute
        # must still refuse a table without the column
        s = self._stub()
        del s.require_staxids
        with self.assertRaises(ValueError) as ctx:
            s.parse_def_table(['q1\ts1\t100\t1e-50\t200\t90'])
        self.assertIn('staxids', str(ctx.exception))
        args = cli_args('search', ['-i', 'x.tsv', '-o', 'out',
                                   '-t', 'taxdump', '--source', 'gtdb'])
        self.assertIsNone(args.require_staxids,
                          'the CLI must not smuggle a default that reads as '
                          'an explicit "no"')

    def _fixture_search(self, tag, *extra, **kw):
        """Run the offline precomp search on a private copy of the fixtures.

        ``bogus_subject=True`` rewrites the first hit row's subject into a stem
        that no route can resolve, which is the situation --keep-unresolved is
        about.
        """
        from neoHGT.search import Search
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        db = join(d, 'db')
        shutil.copytree(FX, db)
        if kw.get('bogus_subject'):
            src = join(db, 'search', 'sample.tsv')
            rows = open(src).read().splitlines()
            head = [r for r in rows if r.startswith('#')]
            body = [r for r in rows if not r.startswith('#')]
            x = body[0].split('\t')
            x[1] = ('GB_GCA_999999999.1_d_Bacteria_p_X_c_X_o_X_f_X_g_X'
                    '_SerialNumber_UNRESOLVABLE_1')
            x[6] = 'na'
            body[0] = '\t'.join(x)
            with open(src, 'w') as fh:
                fh.write('\n'.join(head + body) + '\n')
        out = join(d, tag)
        argv = ['-i', join(db, 'protein', TestFixtureEndToEnd.GENOME + '.faa'),
                '-o', out, '-m', 'precomp',
                '-s', join(db, 'search', 'sample.tsv'),
                '-t', join(db, 'taxdump'), '--source', 'gtdb',
                '-d', db, '--threads', '1'] + list(extra)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            Search()(cli_args('search', argv))
        return out, buf.getvalue(), db

    def test_aln_method_values_resolve_or_refuse(self):
        _out, log, _db = self._fixture_search('aln_auto', '--aln-method', 'auto')
        self.assertIn('Self-alignment method: fast.', log,
                      '"auto" must decide and say what it decided')
        # upstream resolves 'native' against the search method: with a
        # pre-computed table there is no DIAMOND to reuse, so it becomes 'fast'
        _out, log, _db = self._fixture_search('aln_native', '--aln-method',
                                             'native')
        self.assertIn('Self-alignment method: fast.', log)
        _out, log, _db = self._fixture_search('aln_fast', '--aln-method', 'fast')
        self.assertIn('Self-alignment method: fast.', log)
        with self.assertRaises(ValueError) as ctx:
            self._fixture_search('aln_lookup', '--aln-method', 'lookup')
        self.assertIn('self-hit', str(ctx.exception))
        with self.assertRaises(ValueError) as ctx:
            self._fixture_search('aln_precomp', '--aln-method', 'precomp')
        self.assertIn('pre-computed', str(ctx.exception).lower())

    def test_fetch_enable_is_forced_off_in_gtdb_mode(self):
        for value in ('auto', 'yes'):
            _out, log, _db = self._fixture_search(
                f'fetch_{value}', '--fetch-enable', value)
            self.assertIn('forcing --fetch-enable no', log,
                          f'--fetch-enable {value} must report the override')
        # --fetch-enable no needs no override, and the settings summary (which
        # upstream prints before this branch) must then read 'no' unaided
        _out, log, _db = self._fixture_search('fetch_no', '--fetch-enable', 'no')
        self.assertIn('Remote fetch enabled: no.', log)
        self.assertNotIn('forcing --fetch-enable no', log,
                         'nothing to override when the user already said no')

    def test_keep_unresolved_both_values_run_and_agree_on_this_input(self):
        """Both values are handed over; the fixture cannot separate them.

        Every subject either resolves through the embedded map or comes back as
        'na', which is not the same as "unresolvable", so yes and no keep the
        same table. Recorded honestly rather than pretending a difference: the
        missing case (a subject no route can place at all) is registered in
        docs/交付物清单.md, 选项取值偏差, together with the still-unexplained
        53-vs-28 row difference between a run with and without -d.
        """
        counts = []
        for value in ('yes', 'no'):
            out, log, _db = self._fixture_search(
                f'unres_{value}', '--keep-unresolved', value)
            table = join(out, TestFixtureEndToEnd.GENOME + '.tsv')
            self.assertTrue(isfile(table), log[-300:])
            counts.append(self._rows(table))
            self.assertFalse(isfile(join(out, 'unresolved_hits.tsv')),
                             'this input has nothing to list as unresolved')
        self.assertEqual(counts[0], counts[1],
                         'both values must behave alike when nothing is '
                         'unresolved')


    @staticmethod
    def _rows(path):
        return sum(1 for l in open(path) if not l.startswith('#'))




class TestFixtureEndToEnd(unittest.TestCase):
    """IT-4 as a CI test: the real `neoHGT search -m precomp` -> `neoHGT
    analyze` chain on the committed KB-scale fixtures, with no DIAMOND, no
    dataset and no network.

    This only became possible with the derived fixtures of plan §10.1
    (`taxdump/`, `gtdb/`, `db.faa`, `search/sample.tsv`); before them the
    committed data were inputs to the converter, so the search and analyze
    paths could only be exercised against the 406 GB dataset. The fixture hit
    table plants HGT_LIKE proteins per genome whose best hit comes from the
    most distant genome, so recovering exactly those proteins with a GTDB
    donor label is a ground-truth assertion, not merely "it ran".
    """

    GENOME = ('RS_GCF_000026225.1_d_Bacteria_p_Pseudomonadota_c_'
              'Gammaproteobacteria_o_Enterobacterales_f_Enterobacteriaceae'
              '_g_Escherichia')

    def setUp(self):
        if not (isfile(join(FX, 'search', 'sample.tsv'))
                and isfile(join(FX, 'gtdb', 'stem2taxid.tsv.gz'))):
            self.skipTest('derived fixtures missing: run '
                          'tools/make_fixtures.py (plan §10.1)')
        self.d = tmpdir()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        # the committed fixtures are the read-only source of this copy
        self.db = join(self.d, 'db')
        shutil.copytree(FX, self.db)
        self.search = join(self.d, 'search')
        self.analyze = join(self.d, 'analyze')
        self.calls = []

        self._guard = _NetworkGuard(self.calls)
        self._guard.__enter__()
        self.addCleanup(self._guard.__exit__)

    def _planted(self):
        """Fixture queries whose best hit is not their own sequence."""
        best = {}
        for line in open(join(self.db, 'search', 'sample.tsv')):
            x = line.rstrip('\n').split('\t')
            if not x[0].startswith(self.GENOME):
                continue
            if x[0] not in best or float(x[4]) > float(best[x[0]][4]):
                best[x[0]] = x
        return sorted(k for k, v in best.items() if v[1] != k)

    def _run_search(self):
        args = cli_args('search', [
            '-i', join(self.db, 'protein', self.GENOME + '.faa'),
            '-o', self.search, '-m', 'precomp',
            '-s', join(self.db, 'search', 'sample.tsv'),
            '-t', join(self.db, 'taxdump'), '--source', 'gtdb',
            '--threads', '1'])
        self.assertEqual(args.source, 'gtdb')
        s = Search()
        s(args)
        self.assertEqual(s.fetch_enable, 'no',
                         'GTDB mode must switch remote fetching off')
        return join(self.search, self.GENOME + '.tsv')

    def test_the_shipped_cli_route_refuses_the_inert_knobs(self):
        """`--queries` must be refused by the parser route, not only when called.

        TestPrecompRouteRefusesBatchKnobs hands a hand-made object to
        guard_precomp_knobs(); those cases prove the guard's logic but would all
        still pass if no CLI path ever reached it -- the difference rounds 32 and
        42 kept making between "the call site exists" and "a command can arrive
        at it". This runs the shipped parser and Search.__init__ end to end.
        """
        args = cli_args('search', [
            '-i', join(self.db, 'protein', self.GENOME + '.faa'),
            '-o', self.search, '-m', 'precomp',
            '-s', join(self.db, 'search', 'sample.tsv'),
            '-t', join(self.db, 'taxdump'), '--source', 'gtdb',
            '--threads', '1', '--queries', '5'])
        self.assertEqual(args.queries, 5, 'the parser must hand it over intact')
        with self.assertRaises(ValueError) as got:
            Search()(args)
        self.assertIn('--queries', str(got.exception), str(got.exception))
        with self.assertRaises(ValueError) as got2:
            Search()(cli_args('search', [
                '-i', join(self.db, 'protein', self.GENOME + '.faa'),
                '-o', self.search, '-m', 'precomp',
                '-s', join(self.db, 'search', 'sample.tsv'),
                '-t', join(self.db, 'taxdump'), '--source', 'gtdb',
                '--threads', '1', '--extrargs=--seg']))
        self.assertIn('--extrargs', str(got2.exception), str(got2.exception))

    def test_search_resolves_every_subject_through_route_c(self):
        out = self._run_search()
        self.assertTrue(isfile(out), out)
        lines = open(out).read().splitlines()
        hits = [l.split('\t') for l in lines if not l.startswith('#')]
        want = sum(1 for l in open(join(self.db, 'search', 'sample.tsv'))
                   if l.split('\t')[0].startswith(self.GENOME))
        self.assertEqual(len(hits), want,
                         'every fixture hit of this genome must be reported')
        # the planted queries have no self hit, so they carry one row less
        nhits = [int(l.split(': ')[1]) for l in lines if l.startswith('# Hits:')]
        counts = {}
        for line in open(join(self.db, 'search', 'sample.tsv')):
            x = line.rstrip('\n').split('\t')
            if x[0].startswith(self.GENOME):
                counts[x[0]] = counts.get(x[0], 0) + 1
        self.assertEqual(sorted(nhits), sorted(counts.values()))
        self.assertTrue(all(len(x) > 5 and x[5].isdigit() for x in hits),
                        [x[:6] for x in hits][:2])
        tids = {x[5] for x in hits}
        labels = gtdb.load_labels(self.db)
        self.assertTrue(labels)
        for tid in tids:
            self.assertIn(tid, labels, f'taxId {tid} is not in the fixture tree')
        self.assertEqual(self.calls, [])

    def test_self_scores_come_from_real_self_alignment(self):
        out = self._run_search()
        scores = [float(l.split(': ')[1]) for l in
                  open(out).read().splitlines() if l.startswith('# Score:')]
        self.assertEqual(len(scores), 5)
        self.assertTrue(all(s > 0 for s in scores), scores)

    def _run_analyze(self, out, extra=(), omit=()):
        from neoHGT.analyze import Analyze
        tid = gtdb.load_stem2taxid(self.db)[self.GENOME]
        argv = [
            '-i', self.search, '-o', out, '-t', join(self.db, 'taxdump'),
            '--source', 'gtdb', '--input-tax', f'{self.GENOME}:{tid}',
            '--self-rank', 'genus', '--donor-rank', 'genus',
            '--donor-label', 'prefixed', '--close-size', '1',
            '--silhouette', '0']
        if omit:
            # drop a flag and its value: `omit` is what the run looks like when the
            # user never mentioned the option, which is the only way to measure a
            # default instead of asserting one.
            kept, skip = [], False
            for tok in argv:
                if skip and tok == omit[1]:
                    skip = False
                    continue
                if tok == omit[0]:
                    skip = True
                    continue
                kept.append(tok)
            argv = kept
        args = cli_args('analyze', argv + list(extra))
        Analyze()(args)
        return join(out, 'hgts', self.GENOME + '.txt')

    def test_the_donor_label_default_is_the_one_the_docs_promise(self):
        """`--donor-label` has no table default, so what does an ordinary run write?

        计划任务对照表 4.4 says "默认前缀标签" and 交付物清单 D6 lists the same promise,
        while the option table declares `choices` and no `default=`. The only honest
        check is a run with the flag absent: if the donor column comes out bare, the
        promise is being kept by nothing but the sentence.
        """
        self._run_search()
        hgt = self._run_analyze(self.analyze,
                               omit=('--donor-label', 'prefixed'))
        called = [l.split('\t') for l in open(hgt).read().splitlines()]
        self.assertTrue(called, hgt)
        self.assertEqual({r[2] for r in called}, {'g__Nanoarchaeum'},
                         'without --donor-label the run must still deliver the '
                         'documented default (official prefixed GTDB label)')

    def test_analyze_recovers_the_planted_hgts_with_gtdb_donors(self):
        self._run_search()
        hgt = self._run_analyze(self.analyze)

        head = open(join(self.analyze, 'scores.tsv')).readline().split('\t')
        for col in ('self', 'close', 'distal'):
            self.assertIn(col, head, head)
        rows = open(join(self.analyze, 'scores.tsv')).read().splitlines()
        self.assertEqual(len(rows) - 1, 5)

        self.assertTrue(isfile(hgt), hgt)
        called = [l.split('\t') for l in open(hgt).read().splitlines()]
        self.assertEqual([r[0] for r in called], self._planted(),
                         'the planted HGT-like proteins must be the ones found')
        self.assertEqual({r[2] for r in called}, {'g__Nanoarchaeum'},
                         'donors must be reported as official GTDB labels')
        for pdf in ('scatter.pdf', 'distal.kde.pdf', 'close.kde.pdf'):
            self.assertTrue(isfile(join(self.analyze, pdf)), pdf)
        self.assertEqual(self.calls, [])

    def _read(self, path):
        return open(path).read() if isfile(path) else None

    def test_analyze_choice_values_move_the_prediction(self):
        """Every analyze-side choice value must change the run, or be provably
        the documented default that changes nothing. Measured first with
        `neoHGT2GTDB-临时文件/r31_analyze_value_probe.py`, then pinned here:
        round 30 found a declared value (`--mimag any`) that no branch knew,
        and a value can be inert in exactly the same way.
        """
        self._run_search()
        base = self._read(self._run_analyze(join(self.d, 'an_base')))
        expect = {('--self-low', 'yes'): False,
                  ('--self-low', 'no'): True,
                  ('--outliers', 'boxplot'): False,
                  ('--outliers', 'none'): True,
                  ('--outliers', 'zscore'): True,
                  ('--donor-label', 'name'): False,
                  ('--donor-label', 'prefixed'): True,
                  ('--donor-ncbi-map', 'no'): True}
        for (flag, value), same in expect.items():
            got = self._read(self._run_analyze(
                join(self.d, f'an_{flag[2:]}_{value}'), [flag, value]))
            if same:
                self.assertEqual(got, base,
                                 f'{flag} {value} must behave like the default')
            else:
                self.assertNotEqual(got, base,
                                    f'{flag} {value} changed nothing')

    def test_weighted_no_is_refused_on_this_fixture_instead_of_silently_empty(self):
        # counting occurrences instead of weighting by bit score collapses the
        # 5-protein fixture; the GTDB guard must say so rather than predict zero
        self._run_search()
        with self.assertRaises(ValueError) as ctx:
            self._run_analyze(join(self.d, 'an_w_no'), ['--weighted', 'no'])
        self.assertIn('survived the ORFan and outlier filters',
                      str(ctx.exception))
        self._run_analyze(join(self.d, 'an_w_yes'), ['--weighted', 'yes'])

    def test_donor_ncbi_map_yes_writes_the_extra_donor_table(self):
        self._run_search()
        out = join(self.d, 'an_map_yes')
        hgt = self._run_analyze(out, ['--donor-ncbi-map', 'yes'])
        extra = hgt.replace('.txt', '.donor_ncbi.tsv')
        self.assertTrue(isfile(extra), extra)
        head = open(extra).readline().rstrip('\n').split('\t')
        self.assertEqual(head[:6], ['protein', 'silh', 'donor_taxid',
                                    'donor_label', 'donor_rank', 'ncbi_taxids'])
        out2 = join(self.d, 'an_map_no')
        hgt2 = self._run_analyze(out2, ['--donor-ncbi-map', 'no'])
        self.assertFalse(exists(hgt2.replace('.txt', '.donor_ncbi.tsv')),
                         '"no" must not write the table')

    def test_orphans_switch_decides_whether_the_removal_runs(self):
        """`--orphans` is honoured by the caller, not by the remover.

        `Analyze.remove_orphans` always drops rows whose close and distal
        scores are both 0; the option decides whether it is called at all
        (`if not self.orphans`). So the observable difference is in the run,
        and the fixture does contain such proteins -- which is why
        `--weighted no` collapses it entirely.
        """
        self._run_search()
        out = join(self.d, 'an_orphans_no')
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self._run_analyze(out, ['--orphans', 'no'])
        self.assertIn('Removed', buf.getvalue(),
                      '--orphans no must run the ORFan removal')
        self.assertIn('ORFans', buf.getvalue())

        out2 = join(self.d, 'an_orphans_yes')
        buf2 = io.StringIO()
        with contextlib.redirect_stdout(buf2):
            self._run_analyze(out2, ['--orphans', 'yes'])
        self.assertNotIn('ORFans', buf2.getvalue(),
                         '--orphans yes must skip the removal, not run it '
                         'and ignore the answer')
        # the score table itself is written before this step, so both runs
        # score the same proteins; the difference is what may be predicted
        self.assertEqual(self._read(join(out, 'scores.tsv')),
                         self._read(join(out2, 'scores.tsv')))

    def test_outliers_none_keeps_rows_the_other_methods_can_drop(self):
        from neoHGT.analyze import Analyze
        rows = [{'query': f'p{i}', 'self': 1.0, 'close': 1.0,
                 'distal': 1.0 + 0.01 * i} for i in range(12)]
        rows.append({'query': 'outlier', 'self': 1.0, 'close': 1.0,
                     'distal': 900.0})
        for method, expect_drop in (('zscore', True), ('boxplot', True),
                                    ('none', False)):
            me = Analyze.__new__(Analyze)
            me.df = pd.DataFrame(rows)
            me.outliers, me.self_low = method, False
            me.remove_outliers()
            dropped = me.df.shape[0] < len(rows)
            self.assertEqual(dropped, expect_drop,
                             f'--outliers {method} kept {me.df.shape[0]} of '
                             f'{len(rows)} rows')
            if expect_drop:
                self.assertNotIn('outlier', list(me.df['query']),
                                 'the planted extreme must be the one removed')
        def tree():
            return {os.path.relpath(join(r, x), FX): sha256(join(r, x))
                    for r, _, fs in os.walk(FX) for x in fs}

        before = tree()
        self._run_search()
        self.assertEqual(before, tree(),
                         'the test must work on a copy, not on the fixtures')

    def test_from_scores_reuse_repeats_the_same_hgt_list(self):
        """--from-scores skips read_input(), which used to leave self.taxdump as
        the directory path and kill donor-rank reporting."""
        self._run_search()
        first = self._run_analyze(self.analyze)
        again = join(self.d, 'analyze2')
        os.makedirs(again)
        shutil.copy(join(self.analyze, 'scores.tsv'), again)
        second = self._run_analyze(again, ['--from-scores'])
        a = [l.split('\t') for l in open(first).read().splitlines()]
        b = [l.split('\t') for l in open(second).read().splitlines()]
        self.assertEqual([(x[0], x[2]) for x in a], [(x[0], x[2]) for x in b],
                         'reusing the score table must not change which proteins '
                         'are HGTs, nor who the donor is')
        for x, y in zip(a, b):
            # scores.tsv is written with float_format='%g' (6 significant
            # digits), so a reused table can only reproduce the score itself
            # to the last stored digit -- compare membership and donor, not
            # bytes
            self.assertAlmostEqual(float(x[1]), float(y[1]), places=5)
        # The last-digit gap must not be mistaken for run-to-run instability:
        # the same reuse command in a second process reproduces the table
        # byte for byte (real-data check: 0 differing rows of 771).
        third = join(self.d, 'analyze3')
        os.makedirs(third)
        shutil.copy(join(self.analyze, 'scores.tsv'), third)
        again_path = self._run_analyze(third, ['--from-scores'])
        self.assertEqual(open(second, 'rb').read(),
                         open(again_path, 'rb').read(),
                         'two reuse runs on the same table must be identical')

    def test_coarse_donor_ranks_survive_the_from_scores_route(self):
        """`--from-scores` used to skip the GTDB rank translation.

        A GTDB taxdump stores the domain rank as 'superkingdom', so
        `--donor-rank domain` resolved nothing when the score table was reused
        while a full run translated the rank and answered -- the same option
        yielding two different answers depending on the route taken.
        """
        def donors_of(path):
            return [(l.split('\t')[0], l.split('\t')[2])
                    for l in open(path, encoding='utf-8').read().splitlines()]

        self._run_search()
        full = self._run_analyze(join(self.d, 'domain_full'),
                                 ['--donor-rank', 'domain'])
        reuse = join(self.d, 'domain_reuse')
        os.makedirs(reuse)
        shutil.copy(join(self.d, 'domain_full', 'scores.tsv'), reuse)
        again = self._run_analyze(reuse, ['--donor-rank', 'domain',
                                         '--from-scores'])
        a, b = donors_of(full), donors_of(again)
        self.assertEqual(a, b,
                         'donor reporting must not depend on the route taken')
        self.assertTrue(all(d != '0' for _, d in a),
                        'a domain-level donor must actually resolve')

    def test_numeric_bandwidth_is_usable_in_gtdb_mode(self):
        """The --bandwidth help promises a number in [0.1, 1.0]; the parser has
        no type, so it arrived as a string and perform_kde rejected it."""
        self._run_search()
        hgt = self._run_analyze(join(self.d, 'analyze_bw'),
                                ['--bandwidth', '0.3'])
        self.assertTrue(isfile(hgt), hgt)
        self.assertEqual(len(open(hgt).read().splitlines()), len(self._planted()))


class TestDemoToolDonorScope(unittest.TestCase):
    """tools/run_hgt_demo.py must describe the sample it just ran.

    Its donor statistics used to be read from `<cache>/analyze/halo_csN/`, a
    literal 'halo' in the middle of a template: every non-halo run therefore
    reported halo's donors in its own hgt_demo.json. The registers were not
    affected (their donor shares came from separate measurements), but the
    evidence files were, so the sample directory now comes from the run tag.
    """

    def test_reference_dir_is_built_from_the_run_tag(self):
        src = open(join(HERE, '..', '..', 'tools', 'run_hgt_demo.py'),
                   encoding='utf-8').read()
        self.assertNotIn("f'halo_cs", src)
        self.assertIn("f'{tag}_cs{ref_cs}'", src)
        # the recomputed, per-sample numbers live beside the old evidence so a
        # reader can tell which file was affected
        rep = os.environ.get('NEOHGT_CACHE')
        if not rep:
            self.skipTest('NEOHGT_CACHE not set')
        path = join(rep, 'reports', 'R25_donor_statistics_recomputed.json')
        if not isfile(path):
            self.skipTest(f'{path} is a workspace evidence file')
        d = json.load(open(path))['per_sample']
        self.assertIn('bac_bacillota', d)
        self.assertNotEqual(d['bac_bacillota']['own_prefixed_share'],
                            d['bac_bacillota']['json_prefixed_share'],
                            'the recorded share came from another sample')


class TestDegenerateGroupGuard(unittest.TestCase):
    """A group with nothing in it must be diagnosed, not crashed into.

    Upstream `predict_hgt()` only guards against a *constant* group; an empty
    score table (all proteins removed as outliers) or an all-NaN group reaches
    sklearn and dies with "Found array with 0 sample(s)". This is what a real
    bacterial sample hit when --close-tax named genera the database does not
    contain (152 of 158 had no hit; the close column became constant 0 and the
    outlier filter then removed all 4,311 proteins).
    """

    @staticmethod
    def _me(df, close_tax=None):
        from types import SimpleNamespace
        return SimpleNamespace(df=df, output='analyze/x',
                               close_tax=close_tax)

    def test_empty_score_table_is_reported(self):
        import pandas as pd
        from neoHGT.analyze import check_score_groups
        me = self._me(pd.DataFrame({'sample': [], 'protein': [], 'self': [],
                                   'close': [], 'distal': []}))
        with self.assertRaises(ValueError) as ctx:
            check_score_groups(me)
        self.assertIn('no protein survived', str(ctx.exception))

    def test_all_nan_close_group_names_the_supplied_group(self):
        import pandas as pd
        from neoHGT.analyze import check_score_groups
        me = self._me(pd.DataFrame({'self': [1.0, 2.0],
                                    'close': [float('nan')] * 2,
                                    'distal': [1.0, 1.0]}),
                      close_tax='105022102,105022277')
        with self.assertRaises(ValueError) as ctx:
            check_score_groups(me)
        msg = str(ctx.exception)
        self.assertIn('close', msg)
        self.assertIn('--close-tax', msg, 'the message must say how to fix it')

    def test_healthy_table_passes_through(self):
        import pandas as pd
        from neoHGT.analyze import check_score_groups
        me = self._me(pd.DataFrame({'self': [1.0, 2.0], 'close': [0.5, 0.4],
                                   'distal': [1.0, 1.0]}))
        self.assertIsNone(check_score_groups(me))


class TestRunnerChangeList(unittest.TestCase):
    """tools/run.sh must witness removals and renames, and run what it was given.

    The runner is the enforcement point of safety rules 14 and 16: it prints the
    list of files a command changed. Its change list came from
    `find "$root" -type f -newer "$marker"`, a question with no answer for a file
    that is gone and no answer for a renamed one (`mv` keeps the mtime), so a
    deletion and a rename both reported "(none)". It executed the command through
    `eval "$*"`, which drops every quote, so `-- python -c "print('x')"` died on a
    shell syntax error and the log recorded the mangled command rather than the
    one asked for. Both were measured before being fixed (round 52).

    These run the shipped runner in a synthetic workspace of five files, so the
    inventory diff is cheap, and they assert on what the runner *prints*: a test
    that only grepped the source for the new code would pass on a runner that
    prints nothing useful.
    """

    RUN = join(HERE, os.pardir, os.pardir, 'tools', 'run.sh')

    def setUp(self):
        self.root = tmpdir()
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.code = join(self.root, 'neoHGT2GTDB-项目代码')
        for d in (join(self.code, 'tools'),
                  join(self.root, 'neoHGT2GTDB-缓存测试', 'logs'),
                  join(self.root, 'neoHGT2GTDB-临时文件')):
            os.makedirs(d, exist_ok=True)
        shutil.copyfile(self.RUN, join(self.code, 'tools', 'run.sh'))
        self.scratch = join(self.root, 'neoHGT2GTDB-临时文件', 'scratch')
        os.makedirs(self.scratch, exist_ok=True)

    def _run(self, *argv, **kw):
        p = subprocess.run(
            ['bash', join(self.code, 'tools', 'run.sh'), 'test intent',
             'test impact', '--'] + list(argv),
            cwd=self.code, capture_output=True, text=True,
            env={**os.environ, 'NEOHGT_ROOT': self.root,
                 'NEOHGT_WATCH': kw.get('watch', '/nonexistent-watch-dir')})
        self.assertEqual(p.returncode, kw.get('rc', 0), p.stdout + p.stderr)
        return p.stdout

    def _helper(self, name, body):
        path = join(self.scratch, name)
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(body)
        return path

    def test_created_file_is_reported_with_its_category(self):
        made = join(self.scratch, 'made.txt')
        out = self._run(sys.executable,
                        self._helper('make.py', f"open(r'{made}','w').write('x\\n')"))
        rows = [l for l in out.splitlines() if 'made.txt' in l]
        self.assertTrue(rows, f'the runner listed nothing for a new file:\n{out}')
        self.assertTrue(any(l.startswith('  CREATED ') for l in rows), rows)

    def test_deleted_file_is_reported_as_gone(self):
        made = join(self.scratch, 'gone.txt')
        with open(made, 'w', encoding='utf-8') as fh:
            fh.write('x\n')
        out = self._run(sys.executable,
                        self._helper('drop.py', f"import os; os.remove(r'{made}')"))
        self.assertIn('GONE', out, 'a deletion has to be named; "(none)" is not a '
                                   'witness that nothing changed')
        self.assertIn('gone.txt', out)

    def test_rename_inside_the_tree_is_reported(self):
        src = join(self.scratch, 'a.txt')
        dst = join(self.scratch, 'b.txt')
        with open(src, 'w', encoding='utf-8') as fh:
            fh.write('payload\n')
        out = self._run('mv', src, dst)
        self.assertIn('MOVED', out, 'a rename keeps its mtime, so a -newer filter '
                                    'cannot see it at all')
        self.assertTrue(('a.txt' in out and 'b.txt' in out), out)

    def _same_signature(self, paths, text='payload\n', mtime=1700000000):
        """Files the inventory cannot tell apart: identical size and mtime."""
        for p in paths:
            with open(p, 'w', encoding='utf-8') as fh:
                fh.write(text)
            os.utime(p, (mtime, mtime))

    def test_a_bulk_move_pairs_by_name_not_by_hash_order(self):
        """Moving a tree must not print an arbitrary permutation as if it were fact.

        Round 54 moved three scratch trees in one command and the change list
        paired each vanished file with *some* created file of the same size --
        192 rows, most of them pointing at a different tree. The pairing now also
        requires the file name to match, and says so when it cannot.
        """
        src = join(self.scratch, 'treeA')
        dst = join(self.scratch, 'treeB')
        os.makedirs(src)
        script = join(self.scratch, 'bulk.sh')
        self._same_signature([join(src, n) for n in ('x.txt', 'y.txt')])
        with open(script, 'w', encoding='utf-8') as fh:
            fh.write(f'#!/bin/sh\nmkdir -p {dst}\n'
                     f'mv {src}/x.txt {dst}/x.txt\n'
                     f'mv {src}/y.txt {dst}/y.txt\n')
        out = self._run('sh', script)
        pairs = [l.split('  ', 2)[2] for l in out.splitlines()
                 if l.startswith('  MOVED ')]
        self.assertEqual(len(pairs), 2, out)
        for pair in pairs:
            old, new = pair.split(' -> ')
            self.assertEqual(os.path.basename(old), os.path.basename(new),
                             f'a same-name move must not be paired across names: '
                             f'{pair}')
        self.assertIn('moved(ambiguous)=0', out, out)

    def test_an_unpairable_move_is_declared_instead_of_guessed(self):
        src = join(self.scratch, 'treeC')
        dst = join(self.scratch, 'treeD')
        os.makedirs(src)
        script = join(self.scratch, 'bulk2.sh')
        self._same_signature([join(src, n) for n in ('c1.txt', 'c2.txt')])
        with open(script, 'w', encoding='utf-8') as fh:
            fh.write(f'#!/bin/sh\nmkdir -p {dst}\n'
                     f'mv {src}/c1.txt {dst}/d1.txt\n'
                     f'mv {src}/c2.txt {dst}/d2.txt\n')
        out = self._run('sh', script)
        self.assertNotIn('  MOVED ', out,
                         f'two files with one signature cannot be told apart: {out}')
        self.assertIn('moved(ambiguous)=2', out, out)
        # nothing is hidden by the honesty: both ends still appear
        self.assertIn('c1.txt', out)
        self.assertIn('d1.txt', out)

    def test_quoted_argument_reaches_the_program_intact(self):
        out = self._run(sys.executable, '-c', "print('quoted-arg-survived')")
        self.assertIn('quoted-arg-survived', out)
        self.assertNotIn('syntax error near unexpected token', out,
                         'the runner must not re-parse the command it was given')

    def test_a_command_that_changes_nothing_lists_nothing(self):
        # The vacuity control: without this case, "report every file every time"
        # would satisfy the four tests above. The comparison is scoped to the
        # scratch dir, because the first command a workspace ever runs also
        # creates commands.log -- that is the log's own business, not a change
        # the command made.
        before = sorted(os.listdir(self.scratch))
        out = self._run('true')
        self.assertNotIn('  CREATED ', out)
        self.assertNotIn('  MODIFIED ', out)
        self.assertNotIn('  GONE ', out)
        self.assertEqual(before, sorted(os.listdir(self.scratch)))

    def test_the_deletion_reaches_the_command_log(self):
        made = join(self.scratch, 'logged.txt')
        with open(made, 'w', encoding='utf-8') as fh:
            fh.write('x\n')
        self._run(sys.executable,
                  self._helper('drop2.py', f"import os; os.remove(r'{made}')"))
        with open(join(self.root, 'neoHGT2GTDB-缓存测试', 'logs',
                       'commands.log'), encoding='utf-8') as fh:
            log = fh.read()
        self.assertIn('GONE', log, 'the durable record has to carry what the '
                                   'terminal showed')
        self.assertIn('logged.txt', log)


class TestIntroducedConfigKeys(unittest.TestCase):
    """The keys this adaptation added to config.yml have to be read by something.

    Round 40's rule checked only the keys whose *name* contains "gtdb" (16 of the
    93 leaves), so `workspace.code_dir` and `database.source` sat outside its
    population while the summary line said "0 inert config knobs". Both were
    shipped with a sentence attached -- config.yml's own comment promises
    `database.source` selects the source, and gtdb.py's docstring promised the
    workspace root comes from `workspace.code_dir` -- and neither was read. Now
    the rule's population is the diff against the upstream config, and these
    tests hold the two product behaviors the widened rule exposed.
    """

    def _db(self, cfg_source, args_source=None):
        import argparse
        import yaml
        from neoHGT.database import Database, set_source_parameters
        cfg = yaml.safe_load(open(join(HERE, os.pardir, 'config.yml'),
                                  encoding='utf-8'))
        cfg.setdefault('database', {})['source'] = cfg_source
        d = Database.__new__(Database)
        d.cfg = cfg
        d.output = join(os.environ.get('NEOHGT_CACHE', tmpdir()), 'db',
                        '_never_created_')
        d.input = None
        d.gtdb_taxa = None
        kw = {} if args_source is None else {'source': args_source}
        # the shipped config no longer carries a machine dataset path, so the
        # env fallback must name a real directory for the isdir() validation
        if cfg_source in ('gtdb', 'hybrid') and not os.environ.get('GTDB_DIR'):
            os.environ['GTDB_DIR'] = tmpdir()
            self.addCleanup(shutil.rmtree, os.environ['GTDB_DIR'],
                            ignore_errors=True)
            self.addCleanup(os.environ.pop, 'GTDB_DIR', None)
        set_source_parameters(d, argparse.Namespace(**kw))
        return d

    def test_config_source_is_what_the_shipped_comment_promises(self):
        # config.yml says this key is the source used when `database` is started
        # without --source; a run that ignores it makes that sentence false.
        self.assertEqual(self._db('gtdb').source, 'gtdb',
                         'database.source was ignored, so the comment lies')

    def test_the_command_line_still_wins(self):
        self.assertEqual(self._db('gtdb', args_source='ncbi').source, 'ncbi')

    def test_the_shipped_config_keeps_upstream_behavior(self):
        """The value that ships must not turn an NCBI build into something else."""
        self.assertEqual(self._db('ncbi').source, 'ncbi')

    def test_search_and_analyze_do_not_read_this_key(self):
        """The second half of the config comment, held to account.

        The comment says the search and analyze stages do not read database.source
        because their own --source means "auto = detect". If someone wires it
        there without amending the comment, this fails -- a doc sentence and a
        behavior are one artifact, not two.
        """
        tool = self._register()
        reads = tool.config_reads(
            join(HERE, os.pardir, os.pardir),
            ('neoHGT/search.py', 'neoHGT/analyze.py', 'neoHGT/util.py',
             'neoHGT/gtdb.py'))
        self.assertNotIn('database.source', reads,
                         sorted(reads.get('database.source', [])))

    def test_workspace_root_falls_back_to_the_config_ancestry(self):
        from neoHGT import gtdb
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        code = join(d, 'project-code')
        os.makedirs(code, exist_ok=True)
        # workspace_root() returns the *ancestry* of code_dir, which is what
        # config.yml's comment promises: the workspace, not the code folder.
        want = os.path.dirname(os.path.realpath(code))
        saved_env = {k: os.environ.pop(k, None)
                     for k in ('NEOHGT_ROOT', 'NEOHGT_CODE')}
        saved = gtdb.load_config_file
        gtdb.load_config_file = lambda: {'workspace': {'code_dir': code}}
        try:
            self.assertEqual(gtdb.workspace_root(), want,
                             'the docstring names this source, so the code has to '
                             'use it before falling back to the baked-in default')
        finally:
            gtdb.load_config_file = saved
            for k, v in saved_env.items():
                if v is not None:
                    os.environ[k] = v

    def test_the_environment_still_wins_over_the_config(self):
        from neoHGT import gtdb
        saved = gtdb.load_config_file
        gtdb.load_config_file = lambda: {'workspace': {'code_dir': '/nope/code'}}
        try:
            env = os.environ.get('NEOHGT_ROOT')
            if env:
                self.assertEqual(gtdb.workspace_root(), os.path.realpath(env))
        finally:
            gtdb.load_config_file = saved

    def test_no_workspace_config_means_no_baked_in_default(self):
        from neoHGT import gtdb
        saved_env = {k: os.environ.pop(k, None)
                     for k in ('NEOHGT_ROOT', 'NEOHGT_CODE')}
        saved = gtdb.load_config_file
        gtdb.load_config_file = lambda: {'workspace': {}}
        try:
            # the old constant pointed at one developer's machine; without any
            # configuration the guard must degrade, not guess a root
            self.assertIsNone(gtdb.workspace_root())
        finally:
            gtdb.load_config_file = saved
            for k, v in saved_env.items():
                if v is not None:
                    os.environ[k] = v

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_an_unreachable_upstream_dir_is_not_a_vacuous_pass(self):
        """`upstream_dir=''` resolves to the shipped config, and a file diffed
        against itself has no new keys -- a rule that reports "0 unread" from an
        empty population would be greenest when it has looked at least."""
        tool = self._register()
        code = join(HERE, os.pardir, os.pardir)
        doc = join(code, 'docs', '交付物清单.md')
        report, problems = tool.config_key_problems(code, doc, '')
        self.assertIn('heuristic', report['population'], report['population'])
        self.assertGreaterEqual(len(report['declared_own_keys']), 15,
                                'an empty population is not a judgement')
        self.assertEqual(problems, [])

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_the_shipped_upstream_makes_the_population_the_diff(self):
        """The rule must judge the keys we added, not the keys whose name says so."""
        tool = self._register()
        code = join(HERE, os.pardir, os.pardir)
        up = join(code, os.pardir, 'neoHGT-main')
        if not isfile(join(up, 'neoHGT', 'config.yml')):
            self.skipTest('no upstream copy reachable from this checkout')
        report, problems = tool.config_key_problems(
            code, join(code, 'docs', '交付物清单.md'), up)
        self.assertEqual(report['population'], 'upstream diff', report)
        self.assertIn('workspace.code_dir', report['declared_own_keys'],
                      'a key we introduced escaped the population')
        self.assertIn('database.source', report['declared_own_keys'])
        self.assertEqual(problems, [], report)

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_a_key_read_only_through_an_entry_string_is_still_a_read(self):
        """The entry-string half of the population, in both directions.

        `config_reads()` follows `cfg['database']['gtdb']['dir']`; it never saw
        `get_config(self, 'gtdb', 'database.gtdb')`. Adding those paths can only
        make *more* shipped keys look read, which is a rule that risks going quiet
        for the wrong reason -- so this asserts a key the code really does not
        read is still reported when the entry strings are in the population.
        """
        tool = self._register()
        code = join(HERE, os.pardir, os.pardir)
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        for sub in ('neoHGT', 'docs', 'tools'):
            shutil.copytree(join(code, sub), join(d, sub),
                            ignore=shutil.ignore_patterns('__pycache__',
                                                           'tests'))
        path = join(d, 'neoHGT', 'config.yml')
        text = open(path, encoding='utf-8').read()
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text.replace('    taxid_base:',
                                  '    a_knob_nobody_reads: 7\n    taxid_base:',
                                  1))
        report, problems = tool.config_key_problems(
            d, join(d, 'docs', '交付物清单.md'), None)
        self.assertIn('database.gtdb.a_knob_nobody_reads',
                      report['never_read'], report['never_read'])
        self.assertTrue(any('a_knob_nobody_reads' in p for p in problems),
                        problems)
        # ...and the dispositioned dynamic sections are not this rule's finding:
        # rule 23 owns them, so a precomp read must not be counted twice.
        self.assertNotIn('precomp.queries', report['read_but_absent'],
                         sorted(report['read_but_absent']))

    def _register(self):
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt22', path,
                                       loader=SourceFileLoader('vpt22', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_the_rule_population_is_the_upstream_diff(self):
        """A knob this project introduced is judged even when its name is plain."""
        tool = self._register()
        code = join(HERE, os.pardir, os.pardir)
        up = join(code, os.pardir, 'neoHGT-main')
        if not isfile(join(up, 'neoHGT', 'config.yml')):
            self.skipTest('no upstream config to diff against here')
        own = tool.introduced_config_keys(code, up)
        self.assertIn('workspace.code_dir', own)
        self.assertIn('database.source', own)
        self.assertIn('database.gtdb.release', own)
        # an upstream key must not be in the population: judging it would blame
        # this adaptation for a decision it never made
        self.assertNotIn('local.tmpdir', own)
        never = tool.inert_introduced_keys(code, up)
        self.assertEqual(never, [], f'still unread: {never}')
        # ...and the same call must bite when a knob really is unread. A test that
        # only asserts the clean answer is satisfied by a rule muted to report
        # nothing -- the shape round 51 named, applied here to my own new test.
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        for sub in ('neoHGT', 'tools'):
            shutil.copytree(join(code, sub), join(d, sub))
        path = join(d, 'neoHGT', 'config.yml')
        with open(path, encoding='utf-8') as fh:
            text = fh.read()
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text.replace('    taxid_base:',
                                  '    a_knob_nobody_reads: 7\n    taxid_base:', 1))
        self.assertIn('database.gtdb.a_knob_nobody_reads',
                      tool.inert_introduced_keys(d, None),
                      'the rule reports nothing, so it cannot be trusted to '
                      'report the clean tree either')


class TestConfigEntryExpansion(unittest.TestCase):
    """`get_config(self, key, f'search.{key}')` has to be read as the loop means it.

    Round 54's first version collected every `for key in ...` list in a function
    into one candidate set, so `Search.args_wf` -- which reuses `key` for five
    lists -- let each f-string reach every key: 409 "entries", ~270 of them paths
    no single call can form (`grouping.name`, `diamond.threads`, ...). An
    over-generated read set is not a conservative over-approximation here: it
    judges config.yml against paths the program never asks for, and it makes the
    count meaningless, which is the failure shape round 51 named.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt23', path,
                                       loader=SourceFileLoader('vpt23', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _code(self, src):
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'neoHGT'))
        with open(join(d, 'neoHGT', 'search.py'), 'w', encoding='utf-8') as fh:
            fh.write(src)
        with open(join(d, 'neoHGT', 'config.yml'), 'w', encoding='utf-8') as fh:
            fh.write('search:\n  minsize: 10\ndiamond:\n  queries: 4\n')
        return d

    def test_two_loops_sharing_a_name_do_not_cross_multiply(self):
        v = self._tool()
        d = self._code("""arguments = ['g',
                     ['--method', 'm', {'choices': ['diamond']}]]


def args_wf(self):
    for key in 'minsize', 'maxseqs':
        get_config(self, key, f'search.{key}')
    for key in 'queries', 'maxchars':
        get_config(self, key, f'diamond.{key}')
""")
        got, unexpanded = v.get_config_entries(d, ('neoHGT/search.py',))
        self.assertEqual(sorted(got),
                         ['diamond.maxchars', 'diamond.queries',
                          'search.maxseqs', 'search.minsize'], sorted(got))
        self.assertEqual(unexpanded, [])

    def test_an_inner_loop_shadows_the_outer_one(self):
        v = self._tool()
        d = self._code("""arguments = []


def args_wf(self):
    for key in 'a', 'b':
        for key in 'c', 'd':
            get_config(self, key, f'search.{key}')
    for key in 'e', 'f':
        get_config(self, key, f'search.{key}')
""")
        got, unexpanded = v.get_config_entries(d, ('neoHGT/search.py',))
        self.assertEqual(sorted(got), ['search.c', 'search.d', 'search.e',
                                       'search.f'], sorted(got))
        self.assertEqual(unexpanded, [])

    def test_a_loop_it_cannot_read_is_reported_not_invented(self):
        """`for key in 'e':` iterates a *string*, not a tuple of one option.

        The scanner will not guess that the author meant `('e',)`, so the call
        lands in `unexpanded` and the register says so -- which is the difference
        between a blind tool and a wrong one.
        """
        v = self._tool()
        d = self._code("""arguments = []


def args_wf(self):
    for key in 'e':
        get_config(self, key, f'search.{key}')
""")
        got, unexpanded = v.get_config_entries(d, ('neoHGT/search.py',))
        self.assertEqual(got, {}, sorted(got))
        self.assertEqual(unexpanded, ['neoHGT/search.py:6'], unexpanded)

    def test_a_renamed_leaf_still_expands(self):
        v = self._tool()
        d = self._code("""arguments = []


def args_wf(self):
    for key in 'close_size', 'input_cov':
        get_config(self, key, f'grouping.{key.replace("_", "")}')
""")
        got, _ = v.get_config_entries(d, ('neoHGT/search.py',))
        self.assertEqual(sorted(got), ['grouping.closesize',
                                       'grouping.inputcov'], sorted(got))

    def test_a_call_it_cannot_expand_is_reported_not_guessed(self):
        v = self._tool()
        d = self._code("""arguments = []


def args_wf(self, fmt):
    get_config(self, 'queries', fmt)
""")
        got, unexpanded = v.get_config_entries(d, ('neoHGT/search.py',))
        self.assertEqual(got, {}, 'an invented path would be a fake read')
        self.assertEqual(unexpanded, ['neoHGT/search.py:5'], unexpanded)

    def test_the_subscript_resolver_is_call_scoped_too(self):
        """`cfg['search'][key]` went through the same function-wide merge."""
        v = self._tool()
        d = self._code("""arguments = []


def args_wf(self):
    for key in 'minsize', 'maxseqs':
        self.x = self.cfg['search'][key]
    for key in 'queries', 'maxchars':
        self.y = self.cfg['diamond'][key]
""")
        got = v.config_reads(d, ('neoHGT/search.py',))
        self.assertEqual(sorted(got),
                         ['diamond.maxchars', 'diamond.queries',
                          'search.maxseqs', 'search.minsize'], sorted(got))

    def test_the_shipped_tree_forms_no_impossible_paths(self):
        v = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        entries, unexpanded = v.get_config_entries(code,
                                                   v.CONFIG_READ_MODULES)
        self.assertEqual(unexpanded, [], 'the scanner gave up somewhere')
        self.assertIn('search.minsize', entries)
        self.assertIn('precomp.queries', entries,
                      'the dynamic-section read must stay visible')
        for invented in ('grouping.name', 'diamond.threads',
                         'analyze.donor_name', 'remote.method'):
            self.assertNotIn(invented, entries,
                             f'{invented} is a cross of two loops that never '
                             'meet, so nothing reads it')
        # the same census as 缓存测试/reports/R54_entry_expansion.json: a silent
        # change in how many reads the rule can see is a change in what it judges
        self.assertEqual(len(entries), 73, len(entries))


class TestDynamicConfigSections(unittest.TestCase):
    """Rule 23: a config section chosen by an option value needs a checked claim.

    The scanner can enumerate the candidates but cannot know which reach the read,
    so the disposition is written in the deliverables document -- and then the gate
    goes to look for the thing the row says exists: a `raise` that still names
    every option it refuses, or an `if` that really leaves no path carrying the
    value to the read.
    """

    HEAD = ("arguments = ['g',\n"
            "             ['--method', 'm', {'choices': ['remote', 'precomp',\n"
            "                                                 'blast']}]]\n"
            "\n"
            "\n"
            'def args_wf(self):\n')
    RESOLVE = ("    if self.method == 'remote':\n"
               '        if self.db:\n'
               "            self.method = 'blast'\n"
               '        else:\n'
               "            self.method = 'precomp'\n")
    READ = ("    for key in 'queries', 'maxchars':\n"
            "        get_config(self, key, f'{self.method}.{key}')\n")
    GUARD = ("\n"
             "\n"
             'def refuse(self):\n'
             "    if getattr(self, 'queries', None):\n"
             "        raise ValueError('--queries / --maxchars cannot take '\n"
             "                         'effect with -m precomp')\n")

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt24', path,
                                       loader=SourceFileLoader('vpt24', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    @staticmethod
    def _ln(src, needle):
        """1-based line of a needle that must appear exactly once."""
        hits = [i + 1 for i, l in enumerate(src.split('\n')) if needle in l]
        assert len(hits) == 1, (needle, hits)
        return hits[0]

    def _tree(self, src=None, cfg='blast:\n  queries: 4\n  maxchars: 9\n',
              rows=()):
        v = self._tool()
        src = (src or self.HEAD + self.RESOLVE + self.READ + self.GUARD)
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'neoHGT'))
        os.makedirs(join(d, 'docs'))
        with open(join(d, 'neoHGT', 'search.py'), 'w', encoding='utf-8') as fh:
            fh.write(src)
        with open(join(d, 'neoHGT', 'config.yml'), 'w', encoding='utf-8') as fh:
            fh.write(cfg)
        doc = join(d, 'docs', '交付物清单.md')
        site = self._ln(src, 'get_config(self, key')
        with open(doc, 'w', encoding='utf-8') as fh:
            fh.write('## 动态配置段的可达性（门禁规则）\n\n'
                     '| 段名 | 读取点 | 处置 | 代码证据 | 覆盖选项 | 理由 |\n'
                     '| --- | --- | --- | --- | --- | --- |\n')
            for value, kind, witness, flags in rows:
                fh.write(f'| `{value}` | `neoHGT/search.py:{site}` | {kind} '
                         f'| `{witness}` | {flags} | 测试用 |\n')
        return d, doc, v

    def _rows(self, src, kind_precomp='拒绝'):
        return [('remote', '归一',
                 f'neoHGT/search.py:{self._ln(src, "if self.method ==")}',
                 '—'),
                ('precomp', kind_precomp,
                 f'neoHGT/search.py:{self._ln(src, "raise ValueError")}',
                 '`--queries`、`--maxchars`')]

    def test_an_undeclared_runtime_section_is_a_problem(self):
        d, doc, v = self._tree(rows=())
        report, problems = v.dynamic_section_problems(d, doc)
        self.assertEqual(sorted(report['candidates']), ['precomp', 'remote'],
                         report)
        self.assertEqual(len(problems), 2, problems)
        self.assertTrue(any('`remote`' in p for p in problems), problems)
        self.assertTrue(any('`precomp`' in p for p in problems), problems)

    def test_two_checked_dispositions_clear_it(self):
        src = self.HEAD + self.RESOLVE + self.READ + self.GUARD
        d, doc, v = self._tree(src=src, rows=self._rows(src))
        report, problems = v.dynamic_section_problems(d, doc)
        self.assertEqual(problems, [], problems)
        self.assertEqual(report['sections_present'], ['precomp', 'remote'],
                         report)

    def test_a_refusal_that_stops_naming_an_option_is_a_problem(self):
        src = (self.HEAD + self.RESOLVE + self.READ
               + "\n\ndef refuse(self):\n    if getattr(self, 'queries', "
                 "None):\n        raise ValueError('--queries cannot be used"
               " with -m precomp')\n")
        d, doc, v = self._tree(src=src, rows=self._rows(src))
        _report, problems = v.dynamic_section_problems(d, doc)
        self.assertTrue(any('--maxchars' in p for p in problems), problems)

    def test_a_disposition_whose_witness_is_not_a_raise_is_a_problem(self):
        src = self.HEAD + self.RESOLVE + self.READ + self.GUARD
        rows = [('remote', '归一',
                 f'neoHGT/search.py:{self._ln(src, "if self.method ==")}',
                 '—'),
                ('precomp', '拒绝',
                 f'neoHGT/search.py:{self._ln(src, "if getattr")}',
                 '`--queries`、`--maxchars`')]
        d, doc, v = self._tree(src=src, rows=rows)
        _report, problems = v.dynamic_section_problems(d, doc)
        self.assertTrue(any('不是一行 `raise`' in p for p in problems), problems)

    def test_a_declaration_that_omits_an_option_is_a_problem(self):
        src = self.HEAD + self.RESOLVE + self.READ + self.GUARD
        rows = [('remote', '归一',
                 f'neoHGT/search.py:{self._ln(src, "if self.method ==")}',
                 '—'),
                ('precomp', '拒绝',
                 f'neoHGT/search.py:{self._ln(src, "raise ValueError")}',
                 '`--queries`')]
        d, doc, v = self._tree(src=src, rows=rows)
        _report, problems = v.dynamic_section_problems(d, doc)
        self.assertTrue(any('少了' in p and '--maxchars' in p for p in problems),
                        problems)

    def test_a_resolution_with_an_open_branch_is_a_problem(self):
        # the inner `if` loses its else, so `method` can still be 'remote'
        src = (self.HEAD
               + "    if self.method == 'remote':\n"
                 '        if self.db:\n'
                 "            self.method = 'blast'\n"
               + self.READ + self.GUARD)
        d, doc, v = self._tree(src=src, rows=self._rows(src))
        _report, problems = v.dynamic_section_problems(d, doc)
        self.assertTrue(any('仍可能是' in p for p in problems), problems)

    def test_a_value_put_back_before_the_read_voids_the_witness(self):
        src = (self.HEAD + self.RESOLVE
               + "    if self.db:\n"
                 "        self.method = 'remote'\n"
               + self.READ + self.GUARD)
        d, doc, v = self._tree(src=src, rows=self._rows(src))
        _report, problems = v.dynamic_section_problems(d, doc)
        self.assertTrue(any('赋回了' in p for p in problems), problems)

    def test_a_section_the_code_can_no_longer_form_is_a_problem(self):
        src = self.HEAD + self.RESOLVE + self.READ + self.GUARD
        rows = self._rows(src) + [
            ('zzz', '拒绝', f'neoHGT/search.py:{self._ln(src, "raise")}',
             '`--queries`')]
        d, doc, v = self._tree(src=src, rows=rows)
        _report, problems = v.dynamic_section_problems(d, doc)
        self.assertTrue(any('zzz' in p for p in problems), problems)

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_an_unreadable_ledger_is_not_a_pass(self):
        src = self.HEAD + self.RESOLVE + self.READ + self.GUARD
        d, _doc, v = self._tree(src=src, rows=())
        _report, problems = v.dynamic_section_problems(
            d, join(d, 'docs', '不存在的清单.md'))
        self.assertTrue(any('读不到' in p for p in problems), problems)

    def test_a_candidate_the_shipped_config_has_needs_no_row(self):
        """`blast` is a candidate too, but config.yml really has that section."""
        src = self.HEAD + self.RESOLVE + self.READ + self.GUARD
        d, doc, v = self._tree(src=src, rows=(), cfg='blast:\n  queries: 4\n')
        report, _problems = v.dynamic_section_problems(d, doc)
        self.assertNotIn('blast', report['candidates'], report['candidates'])

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_the_shipped_tree_is_dispositioned_and_the_witness_bites(self):
        v = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        doc = join(code, 'docs', '交付物清单.md')
        report, problems = v.dynamic_section_problems(code, doc)
        self.assertEqual(sorted(report['candidates']), ['auto', 'precomp'],
                         report)
        self.assertEqual(problems, [], problems)
        # ...and the same call must bite when the shipped refusal no longer
        # names one of the three options it says it refuses
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        shutil.copytree(join(code, 'neoHGT'), join(d, 'neoHGT'),
                        ignore=shutil.ignore_patterns('__pycache__', 'tests'))
        shutil.copytree(join(code, 'docs'), join(d, 'docs'))
        path = join(d, 'neoHGT', 'search.py')
        text = open(path, encoding='utf-8').read()
        lost = '--extrargs has no search program'
        self.assertEqual(text.count(lost), 1, text.count(lost))
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text.replace(lost, 'those flags have no search program',
                                  1))
        _report, problems = v.dynamic_section_problems(d,
                                                       join(d, 'docs',
                                                            '交付物清单.md'))
        self.assertTrue(any('--extrargs' in p for p in problems),
                        f'the refusal lost an option and nothing noticed: '
                        f'{problems}')


class TestLedgerRowDetection(unittest.TestCase):
    """A ledger is a table whose header says so, not a table four cells wide.

    Round 54 added the dynamic-section table to the same document as the D1-D12
    ledger and the register answered "交付物表里的 `auto` 不是方案编号" seven times:
    `ledger_rows()` had been treating every `|`-row with >= 4 cells as a ledger
    row, and the three prose tables that share the file only escaped because they
    were narrower.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt25', path,
                                       loader=SourceFileLoader('vpt25', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _doc(self, text):
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        path = join(d, 'ledger.md')
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text)
        return path

    def test_a_prose_table_is_not_read_as_a_ledger(self):
        v = self._tool()
        path = self._doc(
            '| ID | 状态 | 产物 | 说明 |\n| --- | --- | --- | --- |\n'
            '| D1 | 已交付 | `x` | y |\n\n'
            '## 动态配置段\n\n| 段名 | 读取点 | 处置 | 证据 | 覆盖 | 理由 |\n'
            '| --- | --- | --- | --- | --- | --- |\n'
            '| `auto` | `neoHGT/search.py:306` | 归一 | `a:1` | — | z |\n')
        self.assertEqual(sorted(v.ledger_rows(path)), ['D1'],
                         v.ledger_rows(path))

    def test_a_ledger_that_continues_after_prose_still_counts(self):
        v = self._tool()
        path = self._doc(
            '| ID | 状态 | 产物 | 说明 |\n| --- | --- | --- | --- |\n'
            '| D1 | 已交付 | `x` | y |\n\n一段说明文字\n\n'
            '| D2 | 已交付 | `z` | w |\n')
        self.assertEqual(sorted(v.ledger_rows(path)), ['D1', 'D2'],
                         v.ledger_rows(path))

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_the_shipped_deliverables_table_has_no_phantom_rows(self):
        v = self._tool()
        rows = v.ledger_rows(join(HERE, os.pardir, os.pardir, 'docs',
                                  '交付物清单.md'))
        self.assertEqual(set(rows), {f'D{i}' for i in range(1, 13)},
                         sorted(rows))


class TestPrecompRouteRefusesBatchKnobs(unittest.TestCase):
    """`-m precomp` cannot honour --queries / --maxchars / --extrargs, so it says so.

    The route reads them from config.yml as `precomp.<key>` (search.py:305-306),
    uses neither when the search table is read whole (search.py:160-161) nor when
    self-alignment is forced to the built-in method (search.py:412-414), and the
    one branch that would honour them carries a test that is always true there.
    A user setting them gets a refusal instead of a silently ignored option.
    """

    class _S(object):
        def __init__(self, **kw):
            self.method = 'precomp'
            self.queries = None
            self.maxchars = None
            self.extrargs = None
            for k, v in kw.items():
                setattr(self, k, v)

    def test_each_inert_option_is_refused_by_name(self):
        for kw, flag in (({'queries': 100}, '--queries'),
                         ({'maxchars': 5000}, '--maxchars'),
                         ({'extrargs': '--seg'}, '--extrargs')):
            with self.assertRaises(ValueError) as got:
                guard_precomp_knobs(self._S(**kw))
            self.assertIn(flag, str(got.exception), str(got.exception))

    def test_the_refusal_names_all_three_when_all_are_given(self):
        with self.assertRaises(ValueError) as got:
            guard_precomp_knobs(self._S(queries=1, maxchars=2, extrargs='--x'))
        msg = str(got.exception)
        for flag in ('--queries', '--maxchars', '--extrargs'):
            self.assertIn(flag, msg, msg)

    def test_a_route_that_can_use_them_is_left_alone(self):
        for method in ('diamond', 'blast', 'remote'):
            self.assertIsNone(
                guard_precomp_knobs(self._S(method=method, queries=100,
                                             maxchars=9, extrargs='--seg')),
                f'{method} really does subset and pass arguments')

    def test_the_values_that_mean_the_same_as_the_route_are_not_refused(self):
        # `--queries 0` is "whole sample" and `--maxchars 0` is "unlimited",
        # which is exactly what this route does; refusing them would be theatre
        self.assertIsNone(guard_precomp_knobs(self._S(queries=0, maxchars=0,
                                                      extrargs='')))
        self.assertIsNone(guard_precomp_knobs(self._S()))

    def test_it_is_wired_into_the_gtdb_argument_route(self):
        src = open(join(HERE, os.pardir, 'search.py'), encoding='utf-8').read()
        tree = ast.parse(src)
        calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                 and getattr(n.func, 'id', '') == 'guard_precomp_knobs']
        self.assertEqual(len(calls), 1, 'the guard must run exactly once')
        holder = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)
                  and n.name == 'set_params_gtdb']
        self.assertEqual(len(holder), 1)
        self.assertTrue(any(c.lineno > holder[0].lineno
                            and c.lineno < holder[0].end_lineno
                            for c in calls),
                        'the call has to be inside the GTDB parameter route')


class TestLedgerSectionHeadings(unittest.TestCase):
    """Rule 24: a heading the gate reads has to be there to be read.

    Rounds 53 and 54 each lost a document heading to one of my own edits, and no
    rule noticed: every one of those readers does `text.split(token, 1)[1] if
    token in text else ''` -- the section is gone, the rule reads an empty
    string, and it reports "nothing to check". The document is now checked against
    the list of sections the register itself depends on.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt26', path,
                                       loader=SourceFileLoader('vpt26', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _doc(self, text):
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'docs'))
        path = join(d, 'docs', '交付物清单.md')
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text)
        return d, path

    def test_a_missing_heading_is_named_with_the_rule_it_silences(self):
        v = self._tool()
        d, _p = self._doc('## 选项取值偏差\n\n内容\n')
        _report, problems = v.ledger_heading_problems(d, 'docs/交付物清单.md')
        self.assertTrue(any('动态配置段的可达性' in p for p in problems),
                        problems)
        self.assertTrue(any('appendix_b_problems' in p for p in problems),
                        'the message has to say which reader went blind')

    def test_a_prose_mention_is_not_a_section(self):
        v = self._tool()
        text = ''.join('## %s\n\n内容\n' % tok
                       for tok, _owner in v.DECLARED_SECTIONS
                       if tok != '断言缺失的既有测试')
        text += ('上面这一节提到「断言缺失的既有测试」这个词，但它不是标题。\n')
        d, _p = self._doc(text)
        _report, problems = v.ledger_heading_problems(d, 'docs/交付物清单.md')
        self.assertEqual(len(problems), 1, problems)
        self.assertTrue(any('断言缺失的既有测试' in p for p in problems), problems)

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_the_shipped_document_keeps_every_section_the_gate_reads(self):
        v = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        report, problems = v.ledger_heading_problems(
            code, join(code, 'docs', '交付物清单.md'))
        self.assertEqual(problems, [], problems)
        self.assertEqual(report['sections'], len(v.DECLARED_SECTIONS), report)

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_an_unreadable_document_refuses_to_look_clean(self):
        v = self._tool()
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        _report, problems = v.ledger_heading_problems(d, 'docs/没有这个文件.md')
        self.assertEqual(len(problems), 1, problems)
        self.assertIn('文件不在', problems[0], problems[0])

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_moving_a_round_section_does_not_lose_a_heading(self):
        """The exact accident, twice in two rounds, made mechanical.

        An edit that anchors on a heading line and forgets to write the heading
        back deletes it. Here the shipped document is copied, one heading removed,
        and rule 24 must name the reader that would have gone silent.
        """
        v = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        os.makedirs(join(d, 'docs'))
        src = join(code, 'docs', '交付物清单.md')
        text = open(src, encoding='utf-8').read()
        lost = '## 无人自动执行的核对工具'
        self.assertEqual(text.count(lost), 1, text.count(lost))
        with open(join(d, 'docs', '交付物清单.md'), 'w',
                  encoding='utf-8') as fh:
            fh.write(text.replace(lost, '', 1))
        _report, problems = v.ledger_heading_problems(
            d, join(d, 'docs', '交付物清单.md'))
        self.assertTrue(any('无人自动执行的核对工具' in p for p in problems),
                        problems)


class TestDocSectionLocator(unittest.TestCase):
    """A section found by the first mention of its title is the wrong section.

    Round 54's pass C went red on this and the defect was in the gate, not in the
    code or the docs: the D5 row names 「动态配置段的可达性」, `split(token, 1)[1]`
    started the slice there, and rule 23 reported both dispositions as missing
    while rule 22 re-reported the three `precomp` keys they cover. Two of the six
    readers had no end anchor at all, so 「命令形态偏差」 was exempting flags from
    the whole rest of the file -- wider than that section's own preamble claims.
    """

    @staticmethod
    def _tool():
        from importlib.machinery import SourceFileLoader
        from importlib.util import module_from_spec, spec_from_file_location
        path = join(HERE, os.pardir, os.pardir, 'tools',
                    'verify_plan_traceability.py')
        spec = spec_from_file_location('vpt27', path,
                                       loader=SourceFileLoader('vpt27', path))
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_prose_mention_does_not_move_the_slice(self):
        v = self._tool()
        text = ('intro line naming 动态配置段的可达性 in prose\n'
                '## 动态配置段的可达性（门禁规则 x）\n\n'
                '| 段名 | 读取点 |\n| --- | --- |\n'
                '| `precomp` | `neoHGT/search.py:306` |\n\n'
                '## 下一节\n\n| `leaked` | row outside |\n')
        sect = v.doc_section(text, '动态配置段的可达性')
        self.assertIn('`precomp`', sect)
        self.assertNotIn('leaked', sect)
        self.assertNotIn('in prose', sect)

    def test_a_slice_ends_at_its_own_heading_level(self):
        v = self._tool()
        text = ('## 甲\n\nbody one\n\n### 甲.1\n\nsubsection stays\n\n'
                '## 乙\n\nbody two\n')
        self.assertIn('subsection stays', v.doc_section(text, '甲'))
        self.assertNotIn('body two', v.doc_section(text, '甲'))
        self.assertEqual(v.doc_section(text, '乙'), '\nbody two\n')

    def test_a_token_that_never_names_a_heading_has_no_section(self):
        v = self._tool()
        text = 'prose mentions 选项取值偏差 and nothing else\n'
        self.assertEqual(v.doc_section(text, '选项取值偏差'), '',
                         'prose must not be read as the section it names')

    def test_the_stop_token_wins_on_the_plan_layout(self):
        v = self._tool()
        text = ('### 附录 B 符号\n\n| row |\n\n'
                '### 附录 B.1 变体\n\n| more |\n\n### 附录 C 命令\n')
        sect = v.doc_section(text, '附录 B', stop='### 附录 C')
        self.assertIn('| more |', sect,
                      'same-level subdivision must stay inside until its stop')
        self.assertNotIn('附录 C 命令', sect,
                         'the stop heading itself must not be inside')
        without = v.doc_section(text, '附录 B')
        self.assertNotIn('| more |', without,
                         'a locator that ignores headings reads the next appendix')

    @unittest.skipIf(not isfile(join(HERE, os.pardir, os.pardir, 'docs', '交付物清单.md')),
                         'internal QC registers are maintained in the project workspace, not in the public release')
    def test_the_shipped_ledger_hazard_is_still_there_and_still_harmless(self):
        """The D5 row really does name the section, and rule 23 still reads it.

        Asserting the hazard exists is the point: without it this test would pass
        for a document that simply stopped mentioning the title.
        """
        v = self._tool()
        code = join(HERE, os.pardir, os.pardir)
        doc = join(code, 'docs', '交付物清单.md')
        text = open(doc, encoding='utf-8').read()
        d5 = [ln for ln in text.splitlines() if ln.startswith('| D5 |')]
        self.assertEqual(len(d5), 1)
        self.assertIn(v.DYNAMIC_LEDGER_HEADING, d5[0],
                      'the regression this test exists for has been removed '
                      'from the document, so the locator is untested again')
        report, problems = v.dynamic_section_problems(code, doc)
        self.assertEqual(problems, [], problems)
        self.assertEqual(report['sections_present'], ['auto', 'precomp'], report)
        self.assertEqual(sorted(report['declarations'].values()),
                         ['归一', '拒绝'], report['declarations'])

    def test_no_reader_locates_our_ledger_by_first_mention(self):
        """Structural guard: the six readers must call `doc_section`."""
        code = join(HERE, os.pardir, os.pardir)
        src = open(join(code, 'tools', 'verify_plan_traceability.py'),
                   encoding='utf-8').read()
        literals = ["'附录 B 符号偏差'", "'命令形态偏差'", "'选项取值偏差'",
                    "'无人自动执行的核对工具'", "'断言缺失的既有测试'",
                    'DYNAMIC_LEDGER_HEADING']
        for lit in literals:
            self.assertNotIn('split(%s, 1)[1]' % lit, src,
                             '%s is still located by its first mention' % lit)
            self.assertIn('doc_section(text, %s)' % lit, src,
                          '%s no longer goes through the heading anchor' % lit)



class TestGateRoundWitness(unittest.TestCase):
    """`tools/gate_round_witness.sh` is the step that makes the witness a requirement.

    The witness itself is tested elsewhere; what had no owner was the claim that a
    round *used* it. The case that would have been easy to skip is the scoping one:
    `commands.log` is append-only, so without the start offset the step compares
    every round's writes against this round's list and is always red -- a permanently
    red check is a check someone disables.
    """

    TOOL = join(HERE, os.pardir, os.pardir, 'tools', 'gate_round_witness.sh')

    def _tree(self, round_no=91):
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        for sub in ('项目代码/tools', '缓存测试/logs', '临时文件'):
            os.makedirs(join(d, 'neoHGT2GTDB-%s' % sub))
        shutil.copyfile(self.TOOL, join(d, 'neoHGT2GTDB-项目代码', 'tools',
                                        'gate_round_witness.sh'))
        self.env = dict(os.environ, NEOHGT_ROOT=d,
                        NEOHGT_CODE=join(d, 'neoHGT2GTDB-项目代码'),
                        NEOHGT_CACHE=join(d, 'neoHGT2GTDB-缓存测试'),
                        NEOHGT_TMP=join(d, 'neoHGT2GTDB-临时文件'))
        self.round_no = round_no
        return d

    def _run(self, d):
        import subprocess
        return subprocess.run(['bash', join(d, 'neoHGT2GTDB-项目代码', 'tools',
                                            'gate_round_witness.sh'),
                               str(self.round_no)],
                              capture_output=True, text=True, env=self.env, cwd=d)

    def _log(self, d, lines):
        with open(join(d, 'neoHGT2GTDB-缓存测试', 'logs', 'commands.log'),
                  'w', encoding='utf-8') as fh:
            fh.write('\n'.join(lines) + '\n')

    def test_unbracketed_round_is_refused(self):
        d = self._tree()
        self._log(d, ['\t\tADD\tneoHGT2GTDB-缓存测试/reports/x.json'])
        r = self._run(d)
        self.assertNotEqual(r.returncode, 0, r.stdout)
        self.assertIn('never bracketed', r.stdout + r.stderr, r.stdout)

    def test_header_only_list_is_not_evidence(self):
        d = self._tree()
        open(join(d, 'neoHGT2GTDB-临时文件', '.rw_91_before'), 'w').close()
        with open(join(d, 'neoHGT2GTDB-缓存测试', 'logs',
                       'round91_files.txt'), 'w', encoding='utf-8') as fh:
            fh.write('# round 91 file witness\n')
        r = self._run(d)
        self.assertNotEqual(r.returncode, 0, r.stdout)
        self.assertIn('does not record a', r.stdout + r.stderr, r.stdout)

    def test_earlier_rounds_writes_do_not_count_against_this_one(self):
        d = self._tree()
        old = '\t\tADD\tneoHGT2GTDB-缓存测试/reports/from_a_previous_round.json\n'
        new = '\t\tADD\tneoHGT2GTDB-缓存测试/reports/this_round.json\n'
        self._log(d, [])
        path = join(d, 'neoHGT2GTDB-缓存测试', 'logs', 'commands.log')
        with open(path, 'a', encoding='utf-8') as fh:
            fh.write(old)
        offset = os.path.getsize(path)
        with open(path, 'a', encoding='utf-8') as fh:
            fh.write(new)
        open(join(d, 'neoHGT2GTDB-临时文件', '.rw_91_before'), 'w').close()
        open(join(d, 'neoHGT2GTDB-临时文件', '.rw_91_clogoffset'),
             'w').write('%d\n' % offset)
        with open(join(d, 'neoHGT2GTDB-缓存测试', 'logs',
                       'round91_files.txt'), 'w', encoding='utf-8') as fh:
            fh.write('# round 91 file witness (start 2026-01-01T00:00:00, '
                     'stop 2026-01-01T01:00:00)\n'
                     '# workspace root: x; files before=10 after=11\n'
                     'CREATED\tneoHGT2GTDB-缓存测试/reports/this_round.json\n')
        r = self._run(d)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertNotIn('from_a_previous_round', r.stdout,
                         'the step is judging the whole log, not this round')

    def test_a_write_the_list_lost_is_reported(self):
        d = self._tree()
        self._log(d, ['\t\tADD\tneoHGT2GTDB-缓存测试/reports/vanished.json'])
        open(join(d, 'neoHGT2GTDB-临时文件', '.rw_91_before'), 'w').close()
        open(join(d, 'neoHGT2GTDB-临时文件', '.rw_91_clogoffset'), 'w').write('0\n')
        with open(join(d, 'neoHGT2GTDB-缓存测试', 'logs',
                       'round91_files.txt'), 'w', encoding='utf-8') as fh:
            fh.write('# round 91 file witness (start 2026-01-01T00:00:00, '
                     'stop 2026-01-01T01:00:00)\n'
                     '# workspace root: x; files before=10 after=11\n'
                     'CREATED\tsomething_else.txt\n')
        r = self._run(d)
        self.assertNotEqual(r.returncode, 0, r.stdout)
        self.assertIn('vanished.json', r.stdout + r.stderr, r.stdout)

if __name__ == '__main__':
    unittest.main()


class TestRoundWitness(unittest.TestCase):
    """Safety rule 16 asks for the round's changed files; a per-command window
    cannot give it. `tools/run.sh` inventories before and after one command, so a
    write made between two commands is in no window -- round 54's change record had
    to be written by hand for exactly that reason. `tools/round_witness.sh` brackets
    the whole round instead, and these cases judge the three ways it could lie:
    reporting nothing because it never started, being restarted mid-round, and
    naming its own snapshots as changes.
    """

    TOOL = join(HERE, os.pardir, os.pardir, 'tools', 'round_witness.sh')

    def _root(self):
        d = tmpdir()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        for sub in ('neoHGT2GTDB-项目代码/tools', 'neoHGT2GTDB-缓存测试/logs',
                    'neoHGT2GTDB-临时文件/scratch'):
            os.makedirs(join(d, sub))
        shutil.copyfile(self.TOOL, join(d, 'neoHGT2GTDB-项目代码', 'tools',
                                        'round_witness.sh'))
        open(join(d, 'neoHGT2GTDB-临时文件', 'scratch', 'seed.txt'),
             'w').write('seed\n')
        return d

    def _run(self, d, action, number):
        import subprocess
        env = dict(os.environ,
                   NEOHGT_ROOT=d,
                   NEOHGT_CODE=join(d, 'neoHGT2GTDB-项目代码'),
                   NEOHGT_CACHE=join(d, 'neoHGT2GTDB-缓存测试'),
                   NEOHGT_TMP=join(d, 'neoHGT2GTDB-临时文件'))
        return subprocess.run(['bash', join(d, 'neoHGT2GTDB-项目代码', 'tools',
                                            'round_witness.sh'), action, number],
                              capture_output=True, text=True, env=env, cwd=d)

    def _report(self, d, number):
        p = join(d, 'neoHGT2GTDB-缓存测试', 'logs', 'round%s_files.txt' % number)
        return open(p, encoding='utf-8').read() if isfile(p) else None

    def test_writes_between_commands_are_named(self):
        d = self._root()
        self.assertEqual(self._run(d, 'start', '91').returncode, 0)
        # outside any command window, which is the case run.sh cannot see
        open(join(d, 'neoHGT2GTDB-临时文件', 'scratch', 'quiet.txt'),
             'w').write('nobody ran a command for this\n')
        os.remove(join(d, 'neoHGT2GTDB-临时文件', 'scratch', 'seed.txt'))
        stop = self._run(d, 'stop', '91')
        self.assertEqual(stop.returncode, 0, stop.stderr)
        text = self._report(d, '91')
        self.assertIn('quiet.txt', text, text)
        self.assertIn('GONE', text, 'a deletion must be named too: '
                                    'it is the case mtime scanning cannot answer')
        self.assertIn('created=1', stop.stdout, stop.stdout)

    def test_the_witness_does_not_report_its_own_snapshots(self):
        d = self._root()
        self._run(d, 'start', '92')
        self._run(d, 'stop', '92')
        text = self._report(d, '92')
        # the header may name the snapshot files it kept (that is provenance); no
        # *change row* may be one of them -- a tool that lists itself as a change
        # every run teaches nothing about changes
        rows = [ln for ln in text.splitlines()
                if ln.split('\t')[0] in ('CREATED', 'MODIFIED', 'GONE', 'MOVED')]
        self.assertEqual([ln for ln in rows if '.rw_' in ln], [], rows)

    def test_stop_without_start_refuses_and_writes_no_report(self):
        d = self._root()
        stop = self._run(d, 'stop', '93')
        self.assertNotEqual(stop.returncode, 0,
                            'an empty witness must not report "no changes"')
        self.assertIn('no start snapshot', stop.stderr + stop.stdout, stop.stderr)
        self.assertIsNone(self._report(d, '93'),
                          'it must refuse before writing, not write and then warn')

    def test_restarting_a_round_is_refused_and_the_window_survives(self):
        d = self._root()
        self._run(d, 'start', '94')
        open(join(d, 'neoHGT2GTDB-临时文件', 'scratch', 'early.txt'),
             'w').write('written before any restart\n')
        again = self._run(d, 'start', '94')
        self.assertNotEqual(again.returncode, 0,
                            'a restart would hide every write before it')
        self.assertIn('already started', again.stderr + again.stdout, again.stderr)
        self._run(d, 'stop', '94')
        self.assertIn('early.txt', self._report(d, '94'))

    def test_a_round_that_changed_says_it_in_words_not_a_blank_list(self):
        d = self._root()
        self._run(d, 'start', '95')
        self._run(d, 'stop', '95')
        self.assertIn('NOTHING CHANGED', self._report(d, '95'))

    def test_a_round_number_that_is_not_one_is_refused(self):
        d = self._root()
        for bad in ('0', 'abc', '55x', ''):
            r = self._run(d, 'start', bad)
            self.assertNotEqual(r.returncode, 0, 'accepted round number %r' % bad)
