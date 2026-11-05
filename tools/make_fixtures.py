#!/usr/bin/env python3
"""Build the GTDB test fixtures from the read-only local dataset.

Safety contract (project rules 10-12, 14):
  * the GTDB dataset is opened for reading only; nothing is written there;
  * at most MAX_GENOMES genomes are copied, and every copy is recorded in a
    manifest together with the source path and both sha256 sums;
  * outputs go only to the workspace folders (cache fixtures + KB-scale test
    data inside the code folder).

usage:
    python tools/make_fixtures.py --limit 6
"""

import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys
from os.path import basename, dirname, exists, isdir, isfile, join

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from neoHGT import gtdb
from neoHGT.util import timestamp

MAX_GENOMES = 10          # hard cap from the project safety rules
KB_RECORDS = 5            # proteins kept in the KB-scale test-data excerpts
HGT_LIKE = 2              # queries per genome given a distal best hit

# (bitscore, pident, query coverage, -log10 evalue) of the fixture hits. A
# vertically inherited protein is topped by its own genome and falls away by
# lineage distance; an HGT-like protein has no self hit at all and its best
# hit comes from the most distant genome.
GRADED = ((2000.0, 100.0, 100.0, 0), (1700.0, 95.0, 98.0, 20),
          (1400.0, 85.0, 93.0, 40), (1100.0, 75.0, 88.0, 60),
          (800.0, 65.0, 83.0, 80), (500.0, 55.0, 78.0, 100))
HGT = ((1500.0, 92.0, 95.0, 140), (1200.0, 80.0, 90.0, 100),
       (900.0, 70.0, 85.0, 70), (600.0, 60.0, 80.0, 50),
       (300.0, 50.0, 75.0, 30))

# patterns chosen to cover the plan's fixture requirements (§10.1)
WANTED = (
    ('cluster_letter_family', 'f_Enterobacteriaceae_A_g_Wigglesworthia'),
    ('archaeal_genome', 'g_Nanoarchaeum'),
    ('placeholder_genus', 'f_DUMJ01_g_DUMJ01'),
    ('hyphen_ranks', 'o_PCC-6307_f_Cyanobiaceae_g_Synechococcus_C'),
    ('reference_genus', 'f_Enterobacteriaceae_g_Escherichia'),
    ('second_genus_same_family', 'f_Enterobacteriaceae_g_Salmonella'),
)


def sha256(path, block=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(block), b''):
            h.update(b)
    return h.hexdigest()


def pick(gtdb_dir, limit):
    """Resolve the wanted patterns to real per-genome FASTA paths."""
    info = gtdb.detect_layout(gtdb_dir, layout='auto')
    files = []
    for dom, rec in info['domains'].items():
        for name in gtdb._list_faa(rec['dir']):
            files.append((dom, join(rec['dir'], name)))
    chosen, seen = [], set()
    for label, pat in WANTED:
        hits = [(dom, p) for dom, p in files if pat in basename(p)]
        if not hits:
            print(f'  WARNING: no genome matches "{pat}" ({label})')
            continue
        dom, p = sorted(hits)[0]
        if p in seen:
            continue
        seen.add(p)
        chosen.append({'label': label, 'pattern': pat, 'domain': dom,
                       'src': p})
        if len(chosen) >= min(limit, MAX_GENOMES):
            break
    if len(chosen) > MAX_GENOMES:
        raise ValueError(f'refusing to copy {len(chosen)} genomes '
                         f'(cap {MAX_GENOMES})')
    return info, chosen


def slice_table(path, header, accs, out, limit=None):
    """Copy the header plus the rows whose accession is in accs."""
    kept = []
    with gtdb.open_reader(path) as f:
        head = f.readline().rstrip('\r\n').split('\t')
        try:
            idx = head.index(header)
        except ValueError:
            idx = 0
        with open(out, 'w') as g:
            g.write('\t'.join(head) + '\n')
            for line in f:
                x = line.rstrip('\r\n').split('\t')
                if len(x) > idx and x[idx] in accs:
                    g.write(line)
                    kept.append(x[idx])
                    if limit and len(kept) >= limit:
                        break
    return kept


def excerpt(src, out, n=KB_RECORDS):
    """Keep the first n complete records (keeps the file KB-sized)."""
    with gtdb.open_reader(src) as f, open(out, 'w') as g:
        keep = 0
        for line in f:
            if line.startswith('>'):
                if keep >= n:
                    break
                keep += 1
            g.write(line)
    return keep


def write_derived(info, root, accs, taxonomy_lines):
    """Emit the converter products and a GTDB hit sample the plan asks for.

    Plan §10.1 lists, besides the copied `.faa` excerpts, the *derived* KB-scale
    fixtures: `taxdump/`, `gtdb/{taxid.map.gz,labels.tsv,lineages.tsv,
    stem2taxid.tsv.gz}`, a merged `db.faa` and `search/sample.tsv` carrying
    `staxids`. Those make `neoHGT search` testable offline (no DIAMOND, no
    dataset), which the copied inputs alone did not allow.

    Deterministic by construction: the hit table is built from the fixture
    sequences themselves, so it does not depend on the DIAMOND version
    installed here. Most queries are topped by their own genome and fall away
    by lineage distance (vertical inheritance); the last HGT_LIKE queries of
    each genome have no self hit and are topped by the most distant genome,
    which is what an HGT looks like and what lets analyze() reach the distal
    group from a fixture (plan IT-4).
    """
    prot = join(root, 'protein')
    slices = [join(root, 'metadata_bac120_sample.tsv'),
              join(root, 'metadata_ar53_sample.tsv')]
    df = gtdb.load_metadata(slices, release=info['release'])
    df = df[df['accession'].isin(accs)]
    lins = [taxonomy_lines[a] for a in sorted(accs) if a in taxonomy_lines]
    td, tmap, stats = gtdb.build_taxdump_from_gtdb(lins)

    gtdb.write_taxdump(td, join(root, 'taxdump'))
    stem2taxid = {}
    for name in sorted(os.listdir(prot)):
        stem = name[:-len('.faa')] if name.endswith('.faa') else name
        lin = taxonomy_lines.get(gtdb.parse_stem(stem).get('accession'))
        if lin and lin in tmap:
            stem2taxid[stem] = tmap[lin]
    if not stem2taxid:
        raise AssertionError('no genome stem resolved to a species taxId; '
                             'the route-C fixture would be unusable')
    gtdb.write_side_files(root, info, df, tmap, td, stats,
                          stem2taxid=stem2taxid)

    # merged db.faa: the excerpts concatenated in file order, verbatim records
    records = []
    for name in sorted(os.listdir(prot)):
        with open(join(prot, name)) as f:
            text = f.read()
        records.extend([r for r in text.split('>') if r.strip()])
    with open(join(root, 'db.faa'), 'w') as f:
        for rec in records:
            f.write('>' + rec.rstrip('\n') + '\n')

    # search/sample.tsv: 7 columns as recommended in the docs
    # (qseqid sseqid pident evalue bitscore qcovhsp staxids).
    # Hits are graded by lineage distance instead of being self-only: a table
    # in which every subject belongs to the query genome would put all hits in
    # the self group, and the close/distal split that analyze.py needs could
    # never be exercised from a fixture. Deterministic: query i takes subject
    # i % n of each donor genome, scores descend by shared-lineage depth.
    seqid_of = {}
    for rec in records:
        seqid = rec.split(None, 1)[0]
        stem = (gtdb.parse_seqid(seqid, info['layout']) or {}).get('stem')
        if stem in stem2taxid:
            seqid_of.setdefault(stem, []).append(seqid)
    lineage_of = {stem: taxonomy_lines[gtdb.parse_stem(stem)['accession']]
                  for stem in seqid_of}

    def shared_depth(a, b):
        """Number of lineage segments two stems have in common."""
        x, y = lineage_of[a].split(gtdb.LINEAGE_SEP), \
            lineage_of[b].split(gtdb.LINEAGE_SEP)
        return sum(1 for u, v in zip(x, y) if u == v)

    rows, graded, hgt_queries = [], 0, 0
    for stem in sorted(seqid_of):
        queries = seqid_of[stem]
        n_hgt = min(HGT_LIKE, max(len(queries) - 1, 0))
        # closest-first: what vertical inheritance looks like
        donors = sorted((s for s in seqid_of if s != stem),
                        key=lambda s: (-shared_depth(stem, s), s))
        # furthest-first: what a horizontally acquired gene looks like
        distal = sorted((s for s in seqid_of if s != stem),
                        key=lambda s: (shared_depth(stem, s), s))
        for i, seqid in enumerate(queries):
            if i < len(queries) - n_hgt:
                rows.append((seqid, seqid, '100.000', '0.0e+00', '2000.0',
                             '100.0', stem2taxid[stem]))
                for k, donor in enumerate(donors):
                    graded += 1
                    b, pid, cov, ev = GRADED[k + 1]
                    rows.append((seqid, seqid_of[donor][i % len(seqid_of[donor])],
                                 f'{pid:.3f}', f'{10 ** -ev:.1e}', f'{b:.1f}',
                                 f'{cov:.1f}', stem2taxid[donor]))
            else:
                hgt_queries += 1
                for k, donor in enumerate(distal):
                    graded += 1
                    b, pid, cov, ev = HGT[k]
                    rows.append((seqid,
                                 seqid_of[donor][i % len(seqid_of[donor])],
                                 f'{pid:.3f}', f'{10 ** -ev:.1e}', f'{b:.1f}',
                                 f'{cov:.1f}', stem2taxid[donor]))
    rows.sort(key=lambda r: (r[0], -float(r[4])))
    search_dir = gtdb.guard_write(join(root, 'search'))
    os.makedirs(search_dir, exist_ok=True)
    with open(join(search_dir, 'sample.tsv'), 'w') as f:
        for r in rows:
            f.write('\t'.join(r) + '\n')

    # meta.json: search/analyze read the layout and provenance from it
    meta = {'source': 'gtdb', 'release': info['release'],
            'local_layout': info['layout'], 'gtdb_dir': info['gtdb_dir'],
            'n_genomes': len(stem2taxid), 'n_proteins': len(records),
            'taxonomy': {'n_nodes': len(td), 'id_base': stats['id_base'],
                         'by_rank': stats['by_rank']},
            'fixture': True,
            'note': 'KB-scale test fixture, not a searchable reference '
                    'database (plan §10.1)'}
    with open(join(root, 'meta.json'), 'w') as f:
        json.dump(meta, f, indent=2, sort_keys=True)

    return {'taxdump_nodes': len(td), 'stem2taxid': len(stem2taxid),
            'db_faa_records': len(records), 'hit_rows': len(rows),
            'cross_genome_hits': graded, 'hgt_like_queries': hgt_queries,
            'meta_json': True}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gtdb-dir', default=os.environ.get('GTDB_DIR'))
    ap.add_argument('--cache', default=os.environ.get('NEOHGT_CACHE'))
    ap.add_argument('--code', default=os.environ.get('NEOHGT_CODE'))
    ap.add_argument('--limit', type=int, default=6)
    ap.add_argument('--dry-run', action='store_true',
                    help='only report what would be copied')
    args = ap.parse_args()

    print(f'GTDB fixture builder started at {timestamp()}')
    info, chosen = pick(args.gtdb_dir, args.limit)
    print(f'  layout={info["layout"]} release={info["release"]} '
          f'genomes to copy={len(chosen)} (cap {MAX_GENOMES})')
    for c in chosen:
        print(f'    [{c["label"]}] {basename(c["src"])} '
              f'({os.path.getsize(c["src"]):,} B)')
    if args.dry_run:
        print('  dry run: nothing written')
        return 0

    fix = gtdb.guard_write(join(args.cache, 'fixtures', 'protein'))
    small = gtdb.guard_write(join(args.code, 'neoHGT', 'tests', 'data',
                                  'GTDB', 'protein'))
    os.makedirs(fix, exist_ok=True)
    os.makedirs(small, exist_ok=True)

    accs = set()
    rows = []
    for c in chosen:
        stem = basename(c['src'])[:-len('.faa')]
        acc = gtdb.parse_stem(stem)['accession']
        accs.add(acc)
        dst = join(fix, basename(c['src']))
        shutil.copy2(c['src'], dst)          # copy only; source untouched
        n = excerpt(c['src'], join(small, basename(c['src'])))
        rows.append({
            'label': c['label'], 'accession': acc, 'domain': c['domain'],
            'source_path': gtdb.realpath(c['src']),
            'source_sha256': sha256(c['src']), 'source_bytes':
            os.path.getsize(c['src']),
            'fixture_copy': dst, 'fixture_sha256': sha256(dst),
            'kb_excerpt': join(small, basename(c['src'])),
            'kb_records': n, 'copied_at': str(timestamp())})

    man = join(args.cache, 'fixtures', 'manifest.tsv')
    with open(man, 'w') as f:
        f.write('# GTDB test fixtures copied from a read-only dataset\n')
        f.write(f'# source dataset: {info["gtdb_dir"]} (never written to)\n')
        f.write(f'# release: {info["release"]}   layout: {info["layout"]}\n')
        f.write(f'# genomes copied: {len(rows)} (cap {MAX_GENOMES})\n')
        f.write('\t'.join(rows[0].keys()) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[k]) for k in rows[0]) + '\n')

    # KB-scale table slices for the unit tests
    for path, header, out in (
            (info['metadata'][0], 'accession', 'metadata_bac120_sample.tsv'),
            (info['metadata'][1], 'accession', 'metadata_ar53_sample.tsv')):
        gtdb.guard_write(join(small, '..', out))
        slice_table(path, header, accs, join(small, '..', out))
    tax = join(small, '..', 'taxonomy_sample.tsv')
    lines = gtdb.read_taxonomy(info['taxonomy'])
    with open(tax, 'w') as g:
        for acc in sorted(accs):
            if acc in lines:
                g.write(f'{acc}\t{lines[acc]}\n')

    # the converter products + a deterministic hit table, so that the search
    # and analyze paths are testable offline (plan §10.1)
    root = join(small, '..')
    derived = write_derived(info, root, accs, lines)

    with open(join(small, '..', 'layout.txt'), 'w') as g:
        g.write(info['layout'] + '\n')
    with open(join(small, '..', 'fixtures.json'), 'w') as g:
        json.dump({'layout': info['layout'], 'release': info['release'],
                   'gtdb_dir': info['gtdb_dir'], 'genomes': rows,
                   'derived': derived}, g, indent=2)

    du = sum(os.path.getsize(join(r, x)) for r, _, fs in os.walk(root)
             for x in fs)
    print(f'  manifest: {man}')
    print(f'  full copies: {fix} ({len(rows)} genomes)')
    print(f'  derived: {derived}')
    print(f'  KB-scale test data: {small} ({du/1024:.0f} KiB total)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
