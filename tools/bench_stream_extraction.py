#!/usr/bin/env python3
"""Measure peak RSS and wall clock of the NCBI extraction, upstream vs streaming.

Two input modes:

  --faa-dir   : use real local FASTA files, copied into an NCBI-style download
                tree inside the workspace (capped at 10 genomes per run, which
                is what the project safety rules allow).
  --synthetic : build N genomes x M proteins from the sequences in --faa-dir,
                so the slope of memory growth can be measured without copying
                hundreds of real genomes.

Both implementations run in separate child processes; the report compares peak
RSS, elapsed time, temporary-disk use and the sha256 of the emitted db.faa and
genome.map.gz, which must be identical.

usage:
    python tools/bench_stream_extraction.py --faa-dir <dir with .faa> \
        --taxdump <db/taxdump> --out <benchmarks/stream_extraction.txt> \
        [--synthetic 200x1000]
"""

import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
from os.path import basename, join, splitext

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from neoHGT import gtdb

CHILD = r'''
import json, os, resource, sys, time
sys.path.insert(0, os.environ['CODE'])
import pandas as pd
from neoHGT.database import Database, extract_genomes_stream_ncbi
from neoHGT.util import read_taxdump

mode, root = sys.argv[1], sys.argv[2]
me = type('X', (), {})()
me.output = root
me.down = os.path.join(root, 'download')
me.tmpdir = os.path.join(root, 'tmp')
os.makedirs(me.tmpdir, exist_ok=True)
me.taxdump = read_taxdump(os.environ['TAXDUMP'])
me.df = pd.read_csv(os.path.join(root, 'df.tsv'), sep='\t', dtype=str)
me.df['proteins'] = 0
me.df['residues'] = 0
t0 = time.time()
if mode == 'stream':
    extract_genomes_stream_ncbi(me)
else:
    Database.extract_genomes(me)
print(json.dumps({
    'mode': mode,
    'seconds': round(time.time() - t0, 2),
    'max_rss_mb': round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 1),
    'db_faa_bytes': os.path.getsize(os.path.join(root, 'db.faa'))}))
'''


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def run_case(mode, root, taxdump, code, python):
    env = dict(os.environ, CODE=code, TAXDUMP=taxdump)
    p = subprocess.run([python, '-c', CHILD, mode, root], env=env,
                       capture_output=True, text=True)
    if p.returncode:
        raise SystemExit(f'{mode} failed:\n{p.stdout}\n{p.stderr}')
    return json.loads([l for l in p.stdout.splitlines()
                       if l.startswith('{')][-1])


def load_faa(path):
    """Read one FASTA file into (accession, sequence) pairs."""
    recs, acc, seq = [], None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if acc:
                    recs.append((acc, ''.join(seq)))
                acc, seq = line[1:].split(None, 1)[0], []
            else:
                seq.append(line.strip())
    if acc:
        recs.append((acc, ''.join(seq)))
    return recs


def write_input(root, genomes):
    """Write [(genome, stem, records)] as an NCBI-style download tree."""
    faa = join(root, 'download', 'faa')
    os.makedirs(faa, exist_ok=True)
    rows = []
    for g, stem, recs in genomes:
        with gzip.open(join(faa, f'{stem}_protein.faa.gz'), 'wt') as f:
            for acc, seq in recs:
                f.write(f'>{acc} benchmark protein\n')
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i + 60] + '\n')
        rows.append((g, '2', f'https://x/y/{stem}'))
    pd.DataFrame(rows, columns=['genome', 'taxid',
                               'ftp_path']).to_csv(join(root, 'df.tsv'),
                                                   sep='\t', index=False)
    return sum(len(r) for _, _, r in genomes)


def build_genomes(args):
    files = sorted(f for f in os.listdir(args.faa_dir) if f.endswith('.faa'))
    pool = []
    for fname in files:
        pool.extend(load_faa(join(args.faa_dir, fname)))
    if not pool:
        raise SystemExit(f'no sequences found under {args.faa_dir}')

    if args.synthetic:
        n, m = (int(x) for x in args.synthetic.lower().split('x'))
        genomes = [(f'G{i:05d}', f'GCF_{i:07d}.1_ASM{i}',
                    [(f'SP_{i:06d}_{j:06d}.1', pool[(i * m + j) % len(pool)][1])
                     for j in range(m)])
                   for i in range(n)]
        return genomes, (f'synthetic {n} genomes x {m} proteins from a pool '
                         f'of {len(pool)} real proteins')
    take = files[:args.max_files]
    genomes = [(f'G{i}', splitext(basename(fn))[0],
                load_faa(join(args.faa_dir, fn)))
               for i, fn in enumerate(take, start=1)]
    return genomes, (f'{len(take)} real genome files '
                     f'(copy cap {args.max_files} genomes)')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--faa-dir', required=True)
    ap.add_argument('--taxdump', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--synthetic', metavar='NXM',
                    help='build N synthetic genomes of M proteins each')
    ap.add_argument('--max-files', type=int, default=10)
    ap.add_argument('--python', default=sys.executable)
    ap.add_argument('--code', default=os.environ.get('NEOHGT_CODE'))
    args = ap.parse_args()

    scratch_root = gtdb.guard_write(join(os.environ['NEOHGT_TMP'], 'scratch'))
    os.makedirs(scratch_root, exist_ok=True)
    work = gtdb.guard_write(tempfile.mkdtemp(prefix='bench_stream_',
                                             dir=scratch_root))
    genomes, desc = build_genomes(args)

    src = join(work, 'input')
    n_proteins = write_input(src, genomes)

    report = {'input': desc, 'proteins': n_proteins,
              'taxdump': args.taxdump, 'cases': {}}
    for mode in ('upstream', 'stream'):
        case_dir = join(work, mode)
        os.makedirs(case_dir)
        subprocess.run(['cp', '-a', f'{src}/.', case_dir], check=True)
        os.makedirs(join(case_dir, 'tmp'), exist_ok=True)
        rec = run_case(mode, case_dir, args.taxdump, args.code, args.python)
        rec['db_faa_sha256'] = sha(join(case_dir, 'db.faa'))
        with gzip.open(join(case_dir, 'genome.map.gz'), 'rt') as f:
            rec['genome_map_sha256'] = gtdb.sha256_text(f.read())
        rec['tmp_bytes'] = sum(os.path.getsize(join(dp, fn))
                               for dp, _, fs in
                               os.walk(join(case_dir, 'tmp')) for fn in fs)
        report['cases'][mode] = rec

    a, b = report['cases']['upstream'], report['cases']['stream']
    report['identical_db_faa'] = a['db_faa_sha256'] == b['db_faa_sha256']
    report['identical_genome_map'] = (
        a['genome_map_sha256'] == b['genome_map_sha256'])
    per_m = n_proteins / 1e6
    report['rss_mb_per_million_proteins'] = {
        k: round(v['max_rss_mb'] / per_m, 1) for k, v in
        report['cases'].items()}
    report['projected_rss_gb_for_333m_proteins'] = {
        k: round(mb * 333.384550 / 1024.0, 1) for k, mb in
        report['rss_mb_per_million_proteins'].items()}

    out = gtdb.guard_write(args.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    lines = [f'input: {desc}; {n_proteins:,} protein records',
             f'taxdump: {args.taxdump}']
    for k, v in report['cases'].items():
        lines.append(
            f'{k:8}: rss {v["max_rss_mb"]:9} MB | extraction '
            f'{v["seconds"]:8} s | db.faa {v["db_faa_bytes"]:13,} B | '
            f'tmpdir {v["tmp_bytes"] / 2 ** 20:9,.1f} MiB')
    lines += [
        f'db.faa identical: {report["identical_db_faa"]}',
        f'genome.map.gz identical: {report["identical_genome_map"]}',
        'RSS per million proteins (MB): ' + json.dumps(
            report['rss_mb_per_million_proteins']),
        'projected RSS for 333,384,550 proteins (GB): ' + json.dumps(
            report['projected_rss_gb_for_333m_proteins']),
        f'scratch kept at: {work}']
    text = '\n'.join(lines) + '\n'
    with open(out, 'w') as f:
        f.write(text)
    with open(os.path.splitext(out)[0] + '.json', 'w') as f:
        json.dump(report, f, indent=1)
    print(text)
    print(f'wrote {out}')
    return 0 if (report['identical_db_faa']
                 and report['identical_genome_map']) else 1


if __name__ == '__main__':
    sys.exit(main())
