#!/usr/bin/env python3
"""D12 / plan §10.4-2: specificity control on conserved single-copy markers.

The local GTDB copy ships the ar53 phylogenomic marker families (53 files, one
per HMM, containing the marker protein of every representative genome). These
are the most strictly vertically inherited genes in an archaeon, so a sane
pipeline must almost never call them horizontal. This tool measures exactly
that on an existing GTDB database and an existing analyze run:

  1. copy the 53 read-only marker files into the workspace, renaming each
     record to '<marker>|<accession>' so the family survives concatenation;
  2. make a DIAMOND database of the query proteome (the same file the search
     used) and map every marker family onto its ortholog in that proteome;
  3. intersect those orthologs with hgts/<sample>.txt and report the false
     positive rate of markers against the background rate of all proteins.

Inputs are read-only; the copied marker file and its sha256 are recorded.

usage:
    python tools/run_marker_fp.py --db <gtdb db> --query-faa <query.faa> \
        --hgts <analyze run>/hgts/query.txt --markers <GTDB>/ar53_marker_genes_reps_r220/faa \
        --out <cache>/reports/D12_markers.json
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
from os import makedirs
from os.path import basename, isfile, join

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from neoHGT import gtdb


def run(cmd):
    """Run a command as an argument list (no shell parsing)."""
    print('# ' + ' '.join(cmd), flush=True)
    p = subprocess.run(cmd)
    if p.returncode:
        raise SystemExit(f'command failed ({p.returncode}): {cmd[0]}')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def collect_markers(src, dest):
    """Concatenate marker families into one file, prefixing record ids."""
    files = sorted(f for f in os.listdir(src) if f.endswith('.faa'))
    manifest = []
    with open(dest, 'w') as fo:
        for fname in files:
            path = join(src, fname)
            marker = basename(fname)[:-len('.faa')]
            n = 0
            with open(path) as fi:
                for line in fi:
                    if line.startswith('>'):
                        fo.write(f'>{marker}|{line[1:].strip()}\n')
                        n += 1
                    else:
                        fo.write(line)
            manifest.append({'marker': marker, 'source': path,
                             'source_sha256': sha256(path), 'records': n})
    return files, manifest


def best_per_marker(tsv):
    """(query protein, bitscore, pident, marker) with the best hit per marker."""
    best = {}
    with open(tsv) as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            if len(x) < 4:
                continue
            marker, sid, pident, bits = (x[0].split('|', 1)[0], x[1],
                                         float(x[2]), float(x[3]))
            if marker not in best or bits > best[marker][1]:
                best[marker] = (sid, bits, pident)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--db', required=True)
    ap.add_argument('--query-faa', required=True)
    ap.add_argument('--hgts', required=True)
    ap.add_argument('--markers', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-pident', type=float, default=30.0)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--diamond', default=os.environ.get('NEOHGT_DIAMOND',
                                                        'diamond'))
    args = ap.parse_args()

    work = gtdb.guard_write(join(os.environ['NEOHGT_TMP'], 'markers'))
    makedirs(work, exist_ok=True)
    rep = gtdb.guard_write(args.out)
    makedirs(os.path.dirname(rep), exist_ok=True)

    markers_faa = gtdb.guard_write(join(work, 'ar53_markers.faa'))
    files, manifest = collect_markers(args.markers, markers_faa)
    log = {'markers_dir': args.markers, 'marker_files': len(files),
           'marker_records': sum(m['records'] for m in manifest),
           'copied_to': markers_faa, 'copied_sha256': sha256(markers_faa),
           'source_writable_by_user': os.access(args.markers, os.W_OK),
           'db': args.db, 'query_faa': args.query_faa,
           'query_sha256': sha256(args.query_faa), 'hgts': args.hgts}

    # 1. DIAMOND db of the query proteome only
    qdb = gtdb.guard_write(join(work, 'query_dmnd'))
    run([args.diamond, 'makedb', '--threads', str(args.threads),
         '--in', args.query_faa, '--db', qdb, '--tmpdir', work])

    # 2. markers -> query proteome
    hits = gtdb.guard_write(join(work, 'marker_hits.tsv'))
    run([args.diamond, 'blastp', '--query', markers_faa, '--db', qdb,
         '--threads', str(args.threads),
         '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'bitscore',
         '--max-target-seqs', '1', '--sensitive',
         '--id', str(args.min_pident), '-o', hits])
    best = best_per_marker(hits)
    log['markers_found_in_query'] = len(best)
    log['query_proteins'] = sum(1 for l in open(args.query_faa)
                               if l.startswith('>'))

    # 3. intersect with the prediction
    predicted = {l.split('\t')[0] for l in open(args.hgts) if l.strip()}
    scores = join(os.path.dirname(os.path.dirname(args.hgts)), 'scores.tsv')
    if isfile(scores):
        n_scored = sum(1 for l in open(scores)) - 1
    else:
        n_scored = log['query_proteins']
    orthologs = {sid for sid, _, _ in best.values()}
    fp = sorted(orthologs & predicted)
    log['scores_table'] = scores
    log['proteins_scored'] = n_scored
    log['predicted_hgts'] = len(predicted)
    log['marker_orthologs'] = len(orthologs)
    log['marker_false_positives'] = len(fp)
    log['marker_fp_rate'] = round(len(fp) / max(1, len(orthologs)), 4)
    log['background_hgt_rate'] = round(len(predicted) / max(1, n_scored), 4)
    log['false_positive_list'] = fp
    log['per_marker'] = {m: {'query_protein': sid, 'bitscore': b,
                             'pident': p, 'predicted_hgt': sid in predicted}
                         for m, (sid, b, p) in sorted(best.items())}
    log['manifest'] = manifest

    with open(rep, 'w') as f:
        json.dump(log, f, indent=1, default=str)
    print(f'markers {log["marker_files"]}: found in query '
          f'{log["markers_found_in_query"]}, orthologs '
          f'{len(orthologs)}, false positives {len(fp)} '
          f'({log["marker_fp_rate"]:.1%}) vs background '
          f'{log["background_hgt_rate"]:.1%}')
    print(f'wrote {rep}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
