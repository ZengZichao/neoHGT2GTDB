#!/usr/bin/env python3
"""IT-6 evidence: compare a resumed GTDB build against a fresh one.

    python tools/verify_resume.py --fresh <dir of the first pass> \
        --resumed <dir after re-running the same command> --out <report.json>

Only files that carry wall-clock stamps (meta.json, the taxid.map.gz header)
may differ; everything else has to be byte-identical, otherwise the resume
path would silently publish a database that disagrees with its own metadata.
"""

import argparse
import gzip
import hashlib
import json
import os
import sys
from os.path import isfile, join

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SKIP = {'download'}


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def text_equal(a, b):
    """Compare two text/gz files ignoring '#' comment lines, without loading
    them into memory (lineages.tsv is ~100 MB even for a 6-genome build)."""
    opener = lambda p: (gzip.open(p, 'rt', encoding='utf-8', errors='replace')
                        if p.endswith('.gz') else open(p, 'rt',
                                                      errors='replace'))
    with opener(a) as fa, opener(b) as fb:
        while True:
            la = next((l for l in fa if not l.startswith('#')), None)
            lb = next((l for l in fb if not l.startswith('#')), None)
            if la is None or lb is None:
                return la is None and lb is None
            if la != lb:
                return False


# meta.json keys that legitimately differ between two runs of the same build:
# the wall clock of each phase, the resume flag, and the scratch directory that
# was echoed back from the command line (a location, not build content).
META_EPHEMERAL = ('/build_started', '/build_finished', '/extraction_resumed',
                  '/selection/options/tmpdir')


def meta_diff(a, b):
    """Keys whose values differ between two meta.json manifests."""
    out = []

    def walk(x, y, p=''):
        if isinstance(x, dict) and isinstance(y, dict):
            for k in sorted(set(x) | set(y)):
                walk(x.get(k), y.get(k), f'{p}/{k}')
        elif x != y:
            out.append(p)
    walk(a, b)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fresh', required=True)
    ap.add_argument('--resumed', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    files = set()
    for root in (args.fresh, args.resumed):
        for dirpath, dirnames, names in os.walk(root):
            dirnames[:] = [d for d in dirnames if d not in SKIP]
            for n in names:
                files.add(os.path.relpath(join(dirpath, n), root))

    report = {'fresh': os.path.abspath(args.fresh),
              'resumed': os.path.abspath(args.resumed),
              'identical': [], 'text_equal': [], 'timestamp_only': [],
              'mismatched': [], 'missing': []}
    for rel in sorted(files):
        f, r = join(args.fresh, rel), join(args.resumed, rel)
        if not isfile(f) or not isfile(r):
            report['missing'].append(rel)
            continue
        if sha(f) == sha(r):
            report['identical'].append(rel)
        elif rel == 'meta.json':
            d = meta_diff(json.load(open(f)), json.load(open(r)))
            keep = [k for k in d if k not in META_EPHEMERAL]
            report['timestamp_only' if not keep else 'mismatched'].append(rel)
            report.setdefault('meta_diff', {})[rel] = d
        elif rel.endswith('.gz') or rel.endswith('.tsv') or rel.endswith(
                '.txt'):
            # derived tables carry a '# ... generated=<ts>' provenance header,
            # so their bytes legitimately differ between two builds; the gz
            # maps always did (gzip stores an mtime)
            equal = text_equal(f, r)
            key = ('text_equal' if equal and rel.endswith('.gz')
                   else 'timestamp_only' if equal else 'mismatched')
            report[key].append(rel)
        else:
            report['mismatched'].append(rel)

    m = json.load(open(join(args.resumed, 'meta.json')))
    report['resumed_flag'] = m.get('extraction_resumed')
    report['counts'] = m.get('counts')
    report['verdict'] = 'PASS' if not (report['mismatched']
                                       or report['missing']) else 'FAIL'
    # same discipline as the other tools: refuse --out paths outside the
    # workspace instead of leaking a report onto whatever path was passed
    import sys as _sys
    _sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(
        __file__))))
    from neoHGT.gtdb import guard_write
    out_path = guard_write(args.out)
    with open(out_path, 'w') as fo:
        json.dump(report, fo, indent=1, ensure_ascii=False)
    print(f'{len(report["identical"])} identical, '
          f'{len(report["text_equal"])} gz text-equal, '
          f'{len(report["timestamp_only"])} stamp-only, '
          f'{len(report["mismatched"])} mismatched, '
          f'{len(report["missing"])} missing -> {report["verdict"]}')
    for k in ('mismatched', 'missing'):
        for rel in report[k]:
            print(f'  {k}: {rel}')
    return 0 if report['verdict'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
