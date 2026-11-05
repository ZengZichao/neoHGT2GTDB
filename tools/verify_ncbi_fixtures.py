#!/usr/bin/env python3
"""RT-1 fixture comparison: run the NCBI analysis path in both trees.

The plan asks that the committed fixtures under ``neoHGT/tests/data/Ecoli``
produce identical results before and after the GTDB adaptation. Upstream's own
``test___call__`` asserts a hardcoded HGT count that no longer holds on current
scipy/scikit-learn (it fails identically in both trees), so comparing against
the *committed* numbers proves nothing. What RT-1 actually needs is a
byte-for-byte comparison of the artifacts produced by the pristine tree and the
adapted tree under the same interpreter and the same random seed.

Outputs are deterministic within a tree (checked by running each tree twice);
tree-to-tree identity is the gate.

usage:
    python tools/verify_ncbi_fixtures.py --baseline /path/to/neoHGT-main \
        --code /path/to/neoHGT2GTDB-项目代码 \
        --out <cache>/reports/RT1_ncbi_fixtures.json
"""

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from os.path import isfile, join

CHILD = r'''
import json, os, sys
sys.path.insert(0, os.environ['CODE'])
import numpy as np
from neoHGT.analyze import Analyze

ecoli = os.path.join(os.environ['CODE'], 'neoHGT', 'tests', 'data', 'Ecoli')
input_dir = os.path.join(ecoli, 'search')
taxdump = os.path.join(ecoli, 'taxdump')
out = sys.argv[1]
seed = int(sys.argv[2])

def args():
    return None


def run(bandwidth, from_scores):
    np.random.seed(seed)
    me = Analyze()
    args.input = None if from_scores else input_dir
    args.output = out
    args.taxdump = taxdump
    args.maxhits = args.evalue = args.identity = args.coverage = None
    args.input_tax = args.self_tax = args.close_tax = None
    args.self_rank = args.close_size = args.distal_top = None
    args.bandwidth = bandwidth
    args.from_scores = from_scores
    args.donor_name = False
    args.donor_rank = None
    me(args)
    return int(me.df['hgt'].sum())


first = run('silverman', False)
second = run('grid', True)
print(json.dumps({'hgt_silverman': first, 'hgt_grid': second}))
'''


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


# matplotlib stamps /CreationDate and /ModDate into every PDF, so two runs of
# the *same* tree already differ byte-wise. PDFs are therefore compared after
# stripping those two fields; everything else is compared byte-for-byte.
PDF_STAMP = re.compile(rb'/(CreationDate|ModDate)\s*\([^)]*\)')


def sha256_pdf(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(PDF_STAMP.sub(b'', f.read())).hexdigest()


def snapshot(root):
    """{relpath: sha256} of every file under root (sorted, stable)."""
    files = {}
    for dirpath, _, names in os.walk(root):
        for name in sorted(names):
            full = os.path.join(dirpath, name)
            rel = os.path.relpath(full, root)
            files[rel] = (sha256_pdf(full) if rel.endswith('.pdf')
                          else sha256(full))
    return files


def run_tree(tree, work, python, seed):
    """Run the analysis twice in one tree; return (snapshot, counters)."""
    counters = []
    first = None
    for rep in ('a', 'b'):
        # a fresh directory per run: no deletion, and no append-into-old-output
        out = tempfile.mkdtemp(prefix=f'ncbi_{rep}_', dir=work)
        env = dict(os.environ, CODE=tree,
                   MPLCONFIGDIR=join(work, 'mpl'),
                   XDG_CACHE_HOME=join(work, 'xdg'),
                   TMPDIR=os.environ.get('NEOHGT_TMP', work))
        p = subprocess.run([python, '-c', CHILD, out, str(seed)], env=env,
                           capture_output=True, text=True)
        if p.returncode:
            raise SystemExit(f'{tree} run {rep} failed:\n{p.stdout}\n{p.stderr}')
        line = [l for l in p.stdout.splitlines() if l.startswith('{')][-1]
        counters.append(json.loads(line))
        snap = snapshot(out)
        if rep == 'a':
            first = snap
        elif snap != first:
            raise SystemExit(f'{tree}: outputs are not deterministic '
                             f'(rep a vs rep b differ)')
    return first, counters[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--baseline', required=True)
    ap.add_argument('--code', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--python', default=sys.executable)
    args = ap.parse_args()

    work = tempfile.mkdtemp(prefix='ncbi_fixtures_',
                            dir=os.environ['NEOHGT_TMP'])
    os.makedirs(join(work, 'mpl'), exist_ok=True)
    os.makedirs(join(work, 'xdg'), exist_ok=True)

    report = {'seed': args.seed, 'scratch': work, 'trees': {}}
    snaps = {}
    for label, tree in (('baseline', args.baseline), ('code', args.code)):
        tdir = join(work, label)
        os.makedirs(tdir)
        snap, counters = run_tree(tree, tdir, args.python, args.seed)
        snaps[label] = snap
        report['trees'][label] = {
            'path': tree, 'files': len(snap), 'counters': counters,
            'file_list': sorted(snap)}

    b, c = snaps['baseline'], snaps['code']
    only_b = sorted(set(b) - set(c))
    only_c = sorted(set(c) - set(b))
    differ = sorted(k for k in set(b) & set(c) if b[k] != c[k])
    report['file_sets_identical'] = not (only_b or only_c)
    report['identical_files'] = len(set(b) & set(c)) - len(differ)
    report['missing_in_code'] = only_b
    report['extra_in_code'] = only_c
    report['differing'] = differ
    report['verdict'] = ('PASS' if report['file_sets_identical'] and not differ
                         else 'FAIL')

    # committed reference: reported, not gated. Upstream's own assertions fail
    # on a modern scipy/scikit-learn, so the committed scores.tsv is stale in
    # *both* trees; recording that separately keeps it from being mistaken for
    # a regression introduced here.
    ref = join(args.baseline, 'neoHGT', 'tests', 'data', 'Ecoli', 'analyze',
               'scores.tsv')
    if isfile(ref):
        report['committed_scores_sha256'] = sha256(ref)
        report['baseline_scores_sha256'] = next(
            (v for k, v in b.items() if k.endswith('scores.tsv')), None)
        report['code_scores_sha256'] = next(
            (v for k, v in c.items() if k.endswith('scores.tsv')), None)
        report['trees_agree_with_each_other'] = (
            report['baseline_scores_sha256'] == report['code_scores_sha256'])
        report['matches_committed_fixture'] = (
            report['code_scores_sha256'] == report['committed_scores_sha256'])
        report['upstream_matches_committed_fixture'] = (
            report['baseline_scores_sha256'] == report['committed_scores_sha256'])

    out = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(report, f, indent=1)
    text = out[:-5] + '.txt' if out.endswith('.json') else out + '.txt'
    with open(text, 'w') as f:
        f.write('# RT-1 NCBI fixture comparison (plan DoD-7)\n')
        for label in ('baseline', 'code'):
            t = report['trees'][label]
            f.write(f'{label:8} {t["path"]}\n         files={t["files"]} '
                    f'hgt(silverman)={t["counters"]["hgt_silverman"]} '
                    f'hgt(grid)={t["counters"]["hgt_grid"]}\n')
        f.write(f'file sets identical: {report["file_sets_identical"]}\n')
        f.write(f'identical files: {report["identical_files"]}\n')
        f.write(f'differing: {report["differing"] or "(none)"}\n')
        if 'committed_scores_sha256' in report:
            f.write('against the committed fixture (informational):\n')
            f.write(f'  both trees agree with each other: '
                    f'{report["trees_agree_with_each_other"]}\n')
            f.write(f'  code matches committed scores.tsv:     '
                    f'{report["matches_committed_fixture"]}\n')
            f.write(f'  upstream baseline matches committed:   '
                    f'{report["upstream_matches_committed_fixture"]}\n')
            f.write('  (a stale committed fixture fails identically in both '
                    'trees; that is a\n   pre-existing upstream condition, '
                    'not a regression from this work)\n')
        f.write(f'verdict: {report["verdict"]}\n')
    print(open(text).read())
    print(f'wrote {out}')
    return 0 if report['verdict'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
