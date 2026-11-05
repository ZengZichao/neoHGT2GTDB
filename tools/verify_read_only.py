#!/usr/bin/env python3
"""Verify that a read-only input dataset was not modified.

Walks the dataset and reports every entry whose mtime is newer than a cutoff.
Stat-only: no file under the dataset is opened for reading or writing, and no
symlink is followed.

usage:
    python tools/verify_read_only.py --dir "$GTDB_DIR" \
        --since "2026-11-01 17:25:00" --out <workspace>/logs/read_only_check.txt
"""

import argparse
import collections
import os
import sys
import time
from calendar import timegm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from neoHGT import gtdb


def data_suffix(name):
    """Meaningful suffix of a dataset file name ('a.faa.gz' -> 'faa').

    The plan's §3.6 reconciliation counts `.faa` / `.fna` / `.tsv` files, and in
    this dataset both plain and gzipped copies exist, so the compression layer
    has to be peeled off before counting.
    """
    base = name[:-3] if name.endswith('.gz') else name
    _, _, ext = base.rpartition('.')
    return ext if ext and ext != base else '(no extension)'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--since', required=True,
                    help='cutoff as UTC "YYYY-MM-DD HH:MM:SS"')
    ap.add_argument('--out', required=True)
    ap.add_argument('--limit', type=int, default=20)
    ap.add_argument('--no-inventory', action='store_true',
                    help='skip the per-extension file census')
    ap.add_argument('--save-hashes',
                    help='write a relpath->sha256 manifest of every regular '
                         'file (run BEFORE a build)')
    ap.add_argument('--check-hashes',
                    help='compare against a manifest written by '
                         '--save-hashes (run AFTER a build); a mismatch is a '
                         'content change that an mtime rollback would hide')
    args = ap.parse_args()

    cutoff = timegm(time.strptime(args.since, '%Y-%m-%d %H:%M:%S'))
    root = os.path.realpath(os.path.expanduser(args.dir))
    n = 0
    newer = []
    errors = 0
    by_ext = collections.Counter()
    by_dir = collections.Counter()
    hashes = {}
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        for name in dirnames + filenames:
            path = os.path.join(dirpath, name)
            n += 1
            try:
                st = os.stat(path, follow_symlinks=False)
            except OSError:
                errors += 1
                continue
            if st.st_mtime > cutoff:
                newer.append((st.st_mtime, path))
            # content evidence: an mtime rolled back below the cutoff defeats
            # the check above, so the manifest records the bytes themselves
            if (args.save_hashes or args.check_hashes) and \
                    os.path.isfile(path) and not os.path.islink(path):
                import hashlib
                h = hashlib.sha256()
                try:
                    with open(path, 'rb') as f:
                        for chunk in iter(lambda: f.read(1 << 20), b''):
                            h.update(chunk)
                    hashes[os.path.relpath(path, root)] = h.hexdigest()
                except OSError:
                    errors += 1
        if not args.no_inventory:
            rel = os.path.relpath(dirpath, root)
            top = rel.split(os.sep)[0] if rel != '.' else '(root)'
            by_dir[top] += len(filenames)
            for name in filenames:
                by_ext[data_suffix(name)] += 1

    lines = [
        '# read-only discipline verification (plan DoD-9 / IT-7)',
        f'dataset: {root}',
        f'cutoff_utc: {args.since}',
        f'checked_at_utc: {time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}',
        f'entries_seen: {n}',
        f'stat_errors: {errors}',
        f'entries_newer_than_cutoff: {len(newer)}',
    ]
    if not args.no_inventory:
        # plan §3.6 rows 6-8 are dataset-wide file counts; they are recorded
        # here (they cannot be equality-checked per build, being properties of
        # this particular third-party copy) and quoted in docs/验收标准对照.md
        lines.append(f'files_by_extension: ' + ' '.join(
            f'{k}={v}' for k, v in sorted(by_ext.items(),
                                          key=lambda kv: -kv[1])[:12]))
        lines.append('files_by_top_level_dir: ' + ' '.join(
            f'{k}={v}' for k, v in sorted(by_dir.items(),
                                          key=lambda kv: -kv[1])[:12]))
    for mt, p in newer[:args.limit]:
        lines.append(f'  NEWER {time.strftime("%F %T", time.gmtime(mt))}Z {p}')
    lines.append('verdict: ' + ('PASS (nothing written into the dataset)'
                                if not newer else
                                'FAIL (entries modified after the cutoff)'))
    hash_verdict = None
    if args.check_hashes:
        import json
        with open(args.check_hashes) as f:
            baseline = json.load(f)
        changed = sorted(k for k, v in baseline.items()
                         if hashes.get(k) != v)
        lost = sorted(k for k in baseline if k not in hashes)
        hash_verdict = 'FAIL' if (changed or lost) else 'PASS'
        lines.append(f'hash_baseline: {args.check_hashes}')
        lines.append(f'files_hashed_now: {len(hashes)}')
        lines.append(f'hash_changed: {len(changed)}'
                     + (f' (first: {changed[:5]})' if changed else ''))
        lines.append(f'hash_missing: {len(lost)}'
                     + (f' (first: {lost[:5]})' if lost else ''))
        lines.append(f'hash_verdict: {hash_verdict}')
    if args.save_hashes:
        import json
        manifest = gtdb.guard_write(os.path.expanduser(args.save_hashes))
        os.makedirs(os.path.dirname(manifest), exist_ok=True)
        with open(manifest, 'w') as f:
            json.dump(hashes, f, indent=0, sort_keys=True)
        lines.append(f'hash_manifest: {manifest} ({len(hashes)} files)')
    report = '\n'.join(lines) + '\n'
    out = gtdb.guard_write(os.path.expanduser(args.out))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        f.write(report)
    print(report)
    failed = bool(newer) or hash_verdict == 'FAIL'
    return 0 if not failed else 1


if __name__ == '__main__':
    raise SystemExit(main())
