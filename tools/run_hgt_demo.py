#!/usr/bin/env python3
"""IT-4 / D12 evidence: a real HGT prediction plus the parameter scan that
justifies the GTDB grouping defaults (plan §10.4 items 1-3, §7.3).

Inputs are read-only; the single query genome is copied into the workspace and
recorded with its sha256 (project rule 12: at most 10 genomes per copy run).

usage:
    python tools/run_hgt_demo.py --db <built gtdb db> --accession <RS_/GB_> \
        --out-dir <workspace>/search/halo
"""

import argparse
import hashlib
import re
import json
import os
import shutil
import subprocess
import sys
import time
from os import makedirs
from os.path import basename, join

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from neoHGT import gtdb

OUTFMT = 'qseqid sseqid pident evalue bitscore qcovhsp staxids'

# Every interpolated value below comes from a CLI flag or an environment
# variable, and the commands run through a shell, so reject anything that a
# shell would treat as syntax before it reaches them.
_UNSAFE = re.compile(r'[;&|`$()<>!\'"\n\r\\]')


def safe(value):
    """Validate a path/scalar that will be interpolated into a command."""
    text = str(value)
    if _UNSAFE.search(text):
        raise SystemExit(f'refusing unsafe argument: {text!r}')
    return text


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def run(cmd):
    print('# ' + cmd, flush=True)
    p = subprocess.run(cmd, shell=True)
    if p.returncode:
        raise SystemExit(f'command failed ({p.returncode}): {cmd}')
    return p


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--db', required=True)
    ap.add_argument('--accession', required=True)
    ap.add_argument('--gtdb-dir', default=os.environ.get('GTDB_DIR'))
    ap.add_argument('--code', default=os.environ.get('NEOHGT_CODE'))
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--python', default=os.environ.get('NEOHGT_PYTHON',
                                                       sys.executable))
    ap.add_argument('--diamond', default=os.environ.get('NEOHGT_DIAMOND',
                                                        'diamond'))
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--close-sizes', default='10,50,100,200,400')
    ap.add_argument('--skip-search', action='store_true',
                    help='reuse an existing hits.tsv')
    args = ap.parse_args(argv)

    # RK8: check the external binaries before touching the dataset, so a missing
    # DIAMOND is an actionable error rather than a late OSError from subprocess
    if not args.skip_search and not shutil.which(args.diamond):
        raise SystemExit(
            f'DIAMOND executable not found: {args.diamond!r}. Install it '
            '(e.g. conda install -c bioconda diamond=2.1.10), point '
            '$NEOHGT_DIAMOND or --diamond at the binary, or reuse an existing '
            'hit table with --skip-search.')

    out = gtdb.guard_write(args.out_dir)
    makedirs(out, exist_ok=True)
    # name the run after its output directory: search tmpdir and analyze dirs
    # must not collide between different query genomes
    tag = re.sub(r'[^A-Za-z0-9._-]+', '_', basename(out)) or 'run'
    log = {'started': str(gtdb.timestamp()), 'db': args.db,
           'accession': args.accession, 'tag': tag}

    # 1. locate the query genome in the read-only dataset and copy it
    info = gtdb.detect_layout(args.gtdb_dir)
    idx = gtdb.stem_index(info)
    hit = idx.get(args.accession)
    if not hit:
        raise SystemExit(f'{args.accession} not found among local per-genome '
                         'FASTA files')
    q = join(out, 'query.faa')
    if not os.path.exists(q):
        gtdb.guard_write(q)
        with open(hit['path']) as fi, open(q, 'w') as fo:
            fo.write(fi.read())
    log['query_source'] = hit['path']
    log['query_sha256'] = sha256(q)
    log['query_source_sha256'] = sha256(hit['path'])
    n_q = sum(1 for l in open(q) if l.startswith('>'))
    log['query_proteins'] = n_q
    print(f'query: {n_q} proteins from {basename(hit["path"])}', flush=True)

    # 2. search
    hits = join(out, 'hits.tsv')
    t0 = time.time()
    if not args.skip_search or not os.path.exists(hits):
        run(safe(f'{args.diamond}') + ' blastp' + f' --query {safe(q)}' +
            f' --db {safe(join(args.db, "diamond", "db"))}' +
            f' --threads {int(args.threads)} --outfmt 6 {OUTFMT}' +
            ' --max-target-seqs 25 --sensitive --evalue 1e-5' +
            f' -o {safe(hits)}')
    with open(hits) as f:
        rows = [l for l in f if l.strip()]
    log['hits'] = len(rows)
    log['hits_with_taxid'] = sum(1 for r in rows
                                 if len(r.split('\t')) > 6
                                 and r.split('\t')[6].strip())
    log['search_seconds'] = round(time.time() - t0, 1)
    print(f'hits: {len(rows)} (staxids non-empty: {log["hits_with_taxid"]})',
          flush=True)

    # 3. neoHGT search (precomputed table)
    # Upstream 'search' appends to <sample>.tsv instead of truncating, so a
    # second run into the same -o directory would double every protein. Move
    # the previous table aside (never delete) so this run starts clean.
    sample = os.path.splitext(basename(q))[0]
    table = join(out, f'{sample}.tsv')
    tmp_search = join(os.environ['NEOHGT_TMP'], 'search', tag)
    makedirs(tmp_search, exist_ok=True)
    if os.path.exists(table):
        prev = f'{table}.prev-{int(time.time())}'
        os.rename(gtdb.guard_write(table), gtdb.guard_write(prev))
        print(f'previous search table kept as {basename(prev)}', flush=True)
    run(f'{safe(args.python)} {safe(join(args.code, "scripts", "neoHGT"))} '
        f'search -i {safe(q)} -m precomp -s {safe(hits)} '
        f'-t {safe(join(args.db, "taxdump"))} -o {safe(out)} '
        f'--tmpdir {safe(join(os.environ["NEOHGT_TMP"], "search", tag))}')
    with open(table) as f:
        n_records = sum(1 for l in f if l.startswith('# ID: '))
    log['search_records'] = n_records
    log['search_table'] = table

    # 4. analysis at several close_size values (plan §10.4-3)
    an_root = gtdb.guard_write(join(os.environ['NEOHGT_CACHE'], 'analyze'))
    tag = basename(out)                       # 'halo', 'entero', ...
    scan = []
    for cs in [int(x) for x in args.close_sizes.split(',')]:
        o = join(an_root, f'{tag}_cs{cs}')
        makedirs(o, exist_ok=True)
        t0 = time.time()
        rc = subprocess.run(
            f'{safe(args.python)} {safe(join(args.code, "scripts", "neoHGT"))}'
            f' analyze -i {safe(table)} '
            f'-t {safe(join(args.db, "taxdump"))} -o {safe(o)} '
            '--self-rank genus --donor-rank genus --donor-label prefixed '
            f'--donor-ncbi-map yes --close-size {int(cs)}',
            shell=True, capture_output=True, text=True)
        tail = (rc.stdout or '').strip().splitlines()
        hgts = join(o, 'hgts', f'{sample}.txt')
        n_hgt = sum(1 for _ in open(hgts)) if os.path.exists(hgts) else 0
        sc = join(o, 'scores.tsv')
        n_score = (sum(1 for _ in open(sc)) - 1) if os.path.exists(sc) else 0
        rec = {'close_size': cs, 'exit': rc.returncode,
               'proteins_scored': n_score, 'predicted_hgts': n_hgt,
               'scatter_pdf': os.path.exists(join(o, 'scatter.pdf')),
               'seconds': round(time.time() - t0, 1),
               'last_lines': tail[-3:]}
        scan.append(rec)
        print(f'close_size={cs}: scored={n_score} hgts={n_hgt} '
              f'scatter={rec["scatter_pdf"]} exit={rc.returncode}', flush=True)
    log['close_size_scan'] = scan

    # 5. donor summary of the reference run: the largest close_size that still
    # produced predictions (a tier too dense for this database predicts none)
    ok = [r for r in scan if r['predicted_hgts'] > 0]
    ref_cs = (ok[-1] if ok else scan[-1])['close_size']
    # the donor statistics belong to THIS run; hardcoding another sample's
    # directory here silently pasted halo's donors into every other report
    ref = join(an_root, f'{tag}_cs{ref_cs}')
    log['donor_reference_close_size'] = ref_cs
    hg = join(ref, 'hgts', 'query.txt')
    if os.path.exists(hg):
        with open(hg) as f:
            log['donor_examples'] = [l.rstrip() for l in f][:15]
        labels = gtdb.load_labels(args.db)
        mp = gtdb.load_ncbi_map(args.db)
        rows = [l.rstrip().split('\t') for l in open(hg)][:200]
        log['donor_prefixed_share'] = sum(
            1 for r in rows if r[2].startswith(tuple(
                gtdb.GTDB_RANK_PREFIXES.values()))) / max(1, len(rows))
        ncbi_tids = [gtdb.resolve_label(r[2], labels) for r in rows]
        log['donors_with_ncbi_map'] = sum(
            1 for t in ncbi_tids if t and mp.get(t, {}).get(
                'ncbi_taxids') not in (None, 'na'))
    log['finished'] = str(gtdb.timestamp())

    rep = join(out, 'hgt_demo.json')
    with open(rep, 'w') as f:
        json.dump(log, f, indent=2, default=str)
    print(f'wrote {rep}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
