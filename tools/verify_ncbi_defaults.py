#!/usr/bin/env python3
"""RT-2: parameter defaults of the NCBI route must be identical before/after.

The GTDB adaptation may not change how the NCBI pipeline interprets its
arguments (plan §10.5 RT-2). This tool runs a set of representative command
lines through Database.set_parameters twice - once with the current code and
once with the pristine upstream copy - and compares the resulting state.

usage:
    python tools/verify_ncbi_defaults.py --baseline <dir of pristine neoHGT> \
        --out <cache>/reports/RT2_defaults.json

Both runs are pure argument parsing: no downloads, no writes outside the
temporary output directory each run is given inside the workspace.
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from os.path import join

CHILD = r'''
import importlib.machinery, importlib.util, json, os, sys
from importlib import import_module

tree, out, cases = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
sys.path.insert(0, tree)

# build the CLI parser exactly the way scripts/neoHGT does, so this compares
# the real command line and not a re-implementation of it. The entry script
# has no .py suffix, hence the explicit source loader.
path = os.path.join(tree, 'scripts', 'neoHGT')
loader = importlib.machinery.SourceFileLoader('neohgt_cli', path)
spec = importlib.util.spec_from_loader('neohgt_cli', loader)
cli = importlib.util.module_from_spec(spec)
loader.exec_module(cli)
modules = [[n, getattr(import_module('neoHGT.' + n), n.capitalize())()]
           for n in cli.module_names]

from neoHGT.database import Database
from neoHGT.util import load_configs


def parse(case):
    sys.argv = ['neoHGT', 'database', '-o', out] + case
    args = cli.parse_args(modules)
    me = Database()
    me.cfg = load_configs()      # each tree reads its own config.yml
    me.set_parameters(args)
    keep = {}
    for attr in ('cats', 'rank', 'sample', 'manual', 'compile', 'threads',
                 'overwrite', 'capital', 'latin', 'exclude', 'genoids',
                 'taxids', 'segment', 'max_genomes', 'seq_source', 'source',
                 'gtdb_dir', 'id_base', 'keep_faa', 'fingerprints'):
        keep[attr] = getattr(me, attr, '<absent>')
    keep['n_selected'] = 0 if not hasattr(me, 'df') else int(len(me.df))
    return keep


print(json.dumps({json.dumps(c): parse(c) for c in cases}, default=str))
'''

CASES = [
    ['--cats', 'bacteria'],
    ['--cats', 'vertebrate_mammalian,viral', '--rank', 'genus',
     '--sample', '5'],
    ['--default'],
    ['--compile', 'both', '--threads', '4'],
    ['--capital', 'no', '--latin', 'yes'],
    ['--overwrite'],
    ['--cats', 'archaea', '--genbank'],
    ['--cats', 'bacteria', '--complete', '--reference'],
]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--baseline', required=True)
    ap.add_argument('--code', default=os.environ.get('NEOHGT_CODE'))
    ap.add_argument('--out', required=True)
    ap.add_argument('--python', default=sys.executable)
    args = ap.parse_args()

    for p in (args.baseline, args.code):
        if not os.path.isfile(join(p, 'neoHGT', 'database.py')):
            raise SystemExit(f'{p} does not look like a neoHGT source tree')

    root = os.environ.get('NEOHGT_TMP') or tempfile.gettempdir()
    work = join(root, 'rt2')
    os.makedirs(work, exist_ok=True)
    results = {}
    for label, tree in (('current', args.code), ('baseline', args.baseline)):
        out_dir = join(work, label)
        os.makedirs(out_dir, exist_ok=True)
        # rewrite the placeholder output path in each case
        cases = [[out_dir if x == 'OUT' else x for x in c] for c in CASES]
        p = subprocess.run([args.python, '-c', CHILD, tree, out_dir,
                            json.dumps(cases)], capture_output=True, text=True)
        if p.returncode:
            raise SystemExit(f'{label} failed:\n{p.stdout}\n{p.stderr}')
        results[label] = json.loads(
            [l for l in p.stdout.splitlines() if l.startswith('{')][-1])

    diffs = []
    new_attrs = set()
    for case, cur in results['current'].items():
        base = results['baseline'].get(case)
        if base is None:
            diffs.append({'case': case, 'current': cur, 'baseline': None})
            continue
        # attributes that only exist after the adaptation are additions, not
        # changes; the upstream ones must match exactly
        bad = {}
        for k, v in base.items():
            if v == '<absent>':          # attribute introduced by the adaptation
                if cur.get(k, '<absent>') != '<absent>':
                    new_attrs.add(k)
                continue
            if cur.get(k) != v:
                bad[k] = {'current': cur.get(k), 'baseline': v}
        for k in set(cur) - set(base):
            if cur[k] != '<absent>':
                new_attrs.add(k)
        if bad:
            diffs.append({'case': case, 'upstream_attrs_differ': bad})

    # the adapted CLI must expose everything the upstream one did, plus the
    # GTDB options; check that no upstream option disappeared
    def flags(tree):
        p = subprocess.run(
            [args.python, '-c',
             f'import sys; sys.path.insert(0, {tree!r});\n'
             'from neoHGT.database import arguments\n'
             'print("\\n".join(i[0] for i in arguments if not '
             'isinstance(i, str)))'],
            capture_output=True, text=True, env=os.environ)
        return set(x for x in p.stdout.splitlines() if x)
    f_base, f_cur = flags(args.baseline), flags(args.code)
    missing = sorted(f_base - f_cur)

    report = {'cases': len(CASES), 'identical': not diffs,
              'diffs': diffs, 'baseline_options': len(f_base),
              'current_options': len(f_cur),
              'new_attributes': sorted(new_attrs),
              'new_options': sorted(f_cur - f_base),
              'removed_options': missing,
              'verdict': 'PASS' if not diffs and not missing else 'FAIL'}
    dest = args.out
    os.makedirs(os.path.dirname(dest) or '.', exist_ok=True)
    with open(dest, 'w') as f:
        json.dump(report, f, indent=1, ensure_ascii=False)
    print(f'{report["cases"]} command lines compared; '
          f'identical: {report["identical"]}; '
          f'removed options: {missing or "(none)"}; '
          f'added options: {len(report["new_options"])} -> '
          f'{report["verdict"]}')
    for d in diffs:
        print('  DIFF', d['case'])
    return 0 if report['verdict'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
