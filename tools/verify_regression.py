#!/usr/bin/env python3
"""Generate the RT-1 regression evidence used by plan DoD-7.

Runs the upstream test suite twice — once against the untouched reference copy
(`--baseline`), once against the adapted code (`--code`) — and proves two things:

1. the *failure set* is identical (zero new failures introduced by the GTDB
   adaptation);
2. the four upstream modules are pure appends (zero deleted lines), so no NCBI
   behaviour could have been changed by editing an existing line.

Nothing outside the workspace is written; pytest caches stay inside TMPDIR.

usage:
    python tools/verify_regression.py --baseline /path/to/neoHGT-main \
        --code /path/to/neoHGT2GTDB-项目代码 \
        --out <cache>/regression/RT1_failure_set_diff.txt
"""

import argparse
import json
import os
import re
import subprocess
import sys

MODULES = ('neoHGT/database.py', 'neoHGT/search.py', 'neoHGT/analyze.py',
           'neoHGT/util.py', 'neoHGT/config.yml', 'scripts/neoHGT')
# upstream suites that are expected to fail identically on both trees
SUITES = 'neoHGT/tests'


def run_suite(root, python, tmpdir):
    """Return (counts line, sorted set of FAILED test ids) for one tree."""
    env = dict(os.environ, PYTHONPATH=root, TMPDIR=tmpdir,
               PYTHONDONTWRITEBYTECODE='1')
    p = subprocess.run(
        [python, '-m', 'pytest', SUITES, '-q', '--tb=no',
         '-p', 'no:cacheprovider', f'--rootdir={root}'],
        cwd=root, env=env, capture_output=True, text=True)
    ids = sorted(x.split(' ')[1] for x in p.stdout.splitlines()
                 if x.startswith('FAILED'))
    tail = [l for l in p.stdout.splitlines()
            if re.search(r'\d+ (failed|passed)', l)]
    if not tail:
        # No summary means the suite never ran (bad interpreter, collection
        # error, crash). Comparing two empty failure sets would look identical,
        # so surface it instead of reporting a green gate.
        raise SystemExit(
            f'{root}: pytest produced no summary line; the comparison would be '
            f'meaningless. Last output:\n{p.stdout[-800:]}\n{p.stderr[-800:]}')
    return (tail[-1], set(ids), p.returncode)


def diff_counts(base_file, cur_file):
    """(deleted, added) line counts of one file, diff-style."""
    if not os.path.isfile(base_file):
        return (None, None)
    p = subprocess.run(['diff', '-u', base_file, cur_file],
                       capture_output=True, text=True)
    removed = sum(1 for l in p.stdout.splitlines()
                  if l.startswith('-') and not l.startswith('---'))
    added = sum(1 for l in p.stdout.splitlines()
                if l.startswith('+') and not l.startswith('+++'))
    return (removed, added)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--baseline', required=True)
    ap.add_argument('--code', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--python', default=sys.executable)
    args = ap.parse_args()

    tmp = os.path.join(os.environ['NEOHGT_TMP'], 'regression_run')
    os.makedirs(tmp, exist_ok=True)

    b_line, b_set, b_rc = run_suite(args.baseline, args.python, tmp)
    c_line, c_set, c_rc = run_suite(args.code, args.python, tmp)
    new_fail = sorted(c_set - b_set)
    fixed = sorted(b_set - c_set)

    audit = {}
    for rel in MODULES:
        rem, add = diff_counts(os.path.join(args.baseline, rel),
                               os.path.join(args.code, rel))
        audit[rel] = {'deleted_lines': rem, 'added_lines': add}

    pure_append = all(v['deleted_lines'] == 0 for k, v in audit.items()
                      if k in MODULES[:5])
    verdict = 'PASS' if not new_fail and pure_append else 'FAIL'
    report = {'baseline': args.baseline, 'code': args.code,
              'baseline_summary': b_line, 'current_summary': c_line,
              'baseline_failed': len(b_set), 'current_failed': len(c_set),
              'baseline_failed_ids': sorted(b_set),
              'current_failed_ids': sorted(c_set),
              'newly_failing': new_fail, 'no_longer_failing': fixed,
              'failure_set_identical': not (new_fail or fixed),
              'deleted_line_audit': audit, 'verdict': verdict}

    out = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    lines = [
        '# RT-1: zero-regression gate (plan DoD-7)',
        f'baseline tree : {args.baseline}',
        f'adapted tree  : {args.code}',
        f'baseline suite: {b_line}  (exit {b_rc})',
        f'adapted suite : {c_line}  (exit {c_rc})',
        f'failure sets  : baseline {len(b_set)}, adapted {len(c_set)}, '
        f'newly failing {len(new_fail)}, no longer failing {len(fixed)}',
        'newly failing : ' + (', '.join(new_fail) if new_fail else '(none)'),
        'no longer failing: ' + (', '.join(fixed) if fixed else '(none)'),
        '# upstream modules must be pure appends (scripts/neoHGT is the one',
        '# known 1-line change: module_names gains "gtdb")',
    ]
    for rel, v in audit.items():
        lines.append(f'{rel:22} deleted_lines={v["deleted_lines"]} '
                     f'added_lines={v["added_lines"]}')
    lines.append(f'verdict: {verdict}')
    text = '\n'.join(lines) + '\n'
    with open(out, 'w') as f:
        f.write(text)
    with open(os.path.splitext(out)[0] + '.json', 'w') as f:
        json.dump(report, f, indent=1)
    print(text)
    print(f'wrote {out}')
    return 0 if verdict == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
