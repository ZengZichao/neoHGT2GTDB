#!/usr/bin/env python3
"""Report the reachability of external URLs used in the user documentation.

This check is advisory: many hosts (GitHub, OneDrive, NCBI) rate-limit or
reject automated requests, so a failure here does not fail the build unless
``--strict`` is given. It exists so a release can see which links changed.
"""

import argparse
import re
import sys
import urllib.error
import urllib.request
from os.path import isfile, join

FILES = (
    'README.md', 'README.zh-CN.md', 'docs/GTDB-数据库.md',
    'docs/en/GTDB-database.md', 'docs/README.md', 'docs/cli-reference.md',
    'docs/config-reference.md', 'docs/troubleshooting.md',
)
URL_RE = re.compile(r'https?://[^\s)>\]]+')


def urls_in(path):
    out = set()
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            for url in URL_RE.findall(line):
                url = url.rstrip('.,;:`')
                url = re.sub(r'[`。，；：、）】]+$', '', url)
                out.add(url)
    return out


def check(url, timeout):
    req = urllib.request.Request(
        url, method='GET', headers={'User-Agent': 'neoHGT2GTDB-doc-check/1.0',
                                    'Range': 'bytes=0-0'})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status
    except urllib.error.HTTPError as exc:
        return exc.code
    except Exception as exc:
        return f'ERROR {type(exc).__name__}: {exc}'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--code', default='.')
    ap.add_argument('--timeout', type=float, default=8.0)
    ap.add_argument('--strict', action='store_true')
    args = ap.parse_args()
    seen = set()
    for rel in FILES:
        path = join(args.code, rel)
        if isfile(path):
            seen |= urls_in(path)
    failures = 0
    for url in sorted(seen):
        status = check(url, args.timeout)
        print(f'{status}\t{url}')
        if isinstance(status, str) or (isinstance(status, int) and status >= 400):
            failures += 1
    print(f'{len(seen)} URL(s); {failures} unreachable or unknown')
    if args.strict and failures:
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
