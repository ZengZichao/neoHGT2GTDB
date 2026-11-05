#!/usr/bin/env python3
"""P0 spike S1/S2: does DIAMOND and BLAST annotate GTDB sequence Ids?

The plan flags this as risk RK1/RK14: local GTDB protein Ids look like
'GB_GCA_000008085.1_d_Archaea_..._g_Nanoarchaeum_SerialNumber_AE017199.1_1'
(>150 chars, containing dots). DIAMOND's --taxonmap matching and makeblastdb's
-taxid_map both have to cope with them.

Inputs are read-only; every artifact goes to a directory chosen by the caller,
which must live inside the workspace. Nothing is downloaded and no NCBI
service is contacted.

usage:
    python tools/spike_s1_s2.py --db-dir <built db> --taxonmap <file> \
        --work <dir inside workspace> [--with-blast]
"""

import argparse
import json
import os
import random
import subprocess
import sys
import time
from os import makedirs
from os.path import exists, isfile, join

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from neoHGT import gtdb
from neoHGT.util import read_taxdump, run_command

OUTFMT6 = 'qseqid sseqid pident length evalue bitscore qcovhsp staxids'


def pick_query(faa, out, n=8, seed=0):
    """Take whole records from a FASTA file, in source order."""
    random.seed(seed)
    recs = []
    sid = None
    seq = []
    with open(faa) as f:
        for line in f:
            if line.startswith('>'):
                if sid is not None and len(''.join(seq)) >= 60:
                    recs.append((sid, ''.join(seq)))
                sid = line[1:].split('\t')[0].split(' ')[0]
                seq = []
            elif sid is not None:
                seq.append(line.strip())
        if sid is not None and len(''.join(seq)) >= 60:
            recs.append((sid, ''.join(seq)))
    if not recs:
        raise ValueError(f'{faa} has no usable records')
    step = max(1, len(recs) // n)
    chosen = recs[::step][:n]
    with open(out, 'w') as f:
        for sid, seq_ in chosen:
            f.write(f'>{sid}\n{seq_}\n')
    return [s for s, _ in chosen]


def taxonmap_variant(src, dst, kind):
    """Rewrite a taxonmap into one of the three variants under test."""
    with open(dst, 'w') as fo:
        fo.write('accession\taccession.version\ttaxid\tgi\n')
        with open(src) as fi:
            next(fi)
            for line in fi:
                x = line.rstrip('\r\n').split('\t')
                if len(x) < 3:
                    continue
                sid, tid = x[1], x[2]
                if kind == 'full-both':
                    a = b = sid
                elif kind == 'rsplit-dot':
                    # what database.py does for NCBI accessions
                    a, b = sid.rsplit('.', 1)[0] if '.' in sid else sid, sid
                elif kind == 'truncate-both':
                    a = b = sid.rsplit('.', 1)[0] if '.' in sid else sid
                else:
                    raise ValueError(kind)
                fo.write(f'{a}\t{b}\t{tid}\tna\n')


def run_diamond(diamond, db, work, query, threads=4, extra=''):
    """Build and search a database, returning staxids coverage."""
    res = {}
    dmnd = join(work, 'db')
    cmd = (f'{diamond} makedb --threads {threads} --in {db} --taxonmap '
           f'{join(work, "taxonmap.tsv")} --taxonnodes '
           f'{join(work, "nodes.dmp")} --taxonnames {join(work, "names.dmp")} '
           f'--db {dmnd} --tmpdir {work} {extra}')
    t0 = time.time()
    ec, out = run_command(cmd)
    res['makedb_exit'] = ec
    res['makedb_seconds'] = round(time.time() - t0, 2)
    res['makedb_output_tail'] = _tail(out)
    if ec:
        return res
    hits = join(work, 'hits.tsv')
    cmd = (f'{diamond} blastp --threads {threads} --db {dmnd} --query {query} '
           f'--outfmt 6 {OUTFMT6} --max-target-seqs 25 --sensitive '
           f'> {hits}')
    ec, out = run_command(cmd)
    res['search_exit'] = ec
    if not exists(hits):
        return res
    rows = [l.rstrip('\r\n').split('\t') for l in open(hits) if l.strip()]
    res['hits'] = len(rows)
    res['self_hits'] = sum(1 for r in rows if len(r) > 7 and r[0] == r[1])
    res['with_staxids'] = sum(1 for r in rows if len(r) > 7
                              and r[7].strip())
    res['staxids_rate'] = (res['with_staxids'] / res['hits']
                           if res['hits'] else 0.0)
    res['distinct_staxids'] = len({r[7] for r in rows if len(r) > 7
                                  and r[7].strip()})
    res['example'] = [r[:2] + [r[7]] for r in rows[:3] if len(r) > 7]
    return res


def _tail(out, n=6):
    """Last n lines of run_command output (which returns a list)."""
    lines = out if isinstance(out, list) else (out or '').splitlines()
    return lines[-n:]


def _joined(out):
    lines = out if isinstance(out, list) else (out or '').splitlines()
    return '\n'.join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--db-dir', required=True,
                    help='a built GTDB database directory (read-only here)')
    ap.add_argument('--taxonmap', required=True,
                    help='streamed taxonmap produced by the build')
    ap.add_argument('--work', required=True,
                    help='scratch directory inside the workspace')
    ap.add_argument('--diamond', default=os.environ.get('NEOHGT_DIAMOND',
                                                        'diamond'))
    ap.add_argument('--with-blast', action='store_true')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--out', help='where to write the JSON report')
    args = ap.parse_args()

    work = gtdb.guard_write(args.work)
    makedirs(work, exist_ok=True)
    db_faa = join(args.db_dir, 'db.faa')
    report = {'db_dir': args.db_dir, 'work': work}
    ec, out = run_command(f'{args.diamond} --version')
    report['diamond_version'] = _joined(out).strip()

    # write the taxdump in the exact format the engines expect
    td = read_taxdump(join(args.db_dir, 'taxdump'))
    used = gtdb.taxids_in_taxonmap(args.taxonmap)
    gtdb.write_compact_taxdump(td, used, join(work, 'taxdump'))
    report['taxa_referenced_by_map'] = len(used)
    report['taxa_in_full_dump'] = len(td)

    query = join(work, 'query.faa')
    ids = pick_query(db_faa, query, n=8)
    report['query_ids'] = [i[:70] + ('...' if len(i) > 70 else '')
                           for i in ids]
    report['max_query_id_len'] = max(len(i) for i in ids)

    # S1: taxonmap variants, with and without DIAMOND's seqid parsing
    variants = {}
    trials = [('full-both', ''),
              ('full-both-noparse', '--no-parse-seqids'),
              ('rsplit-dot', ''),
              ('truncate-both', '')]
    for kind, extra in trials:
        w = join(work, 's1_' + kind)
        makedirs(w, exist_ok=True)
        taxonmap_variant(args.taxonmap, join(w, 'taxonmap.tsv'),
                         'full-both' if 'noparse' in kind else kind)
        for f in ('nodes.dmp', 'names.dmp'):
            with open(join(w, f), 'w') as fo, open(
                    join(work, 'taxdump', f)) as fi:
                fo.write(fi.read())
        variants[kind] = run_diamond(args.diamond, db_faa, w, query,
                                     args.threads, extra)
        report[f's1_{kind}_cmd_extra'] = extra
        print(f'S1 {kind}: makedb={variants[kind].get("makedb_exit")} '
              f'{variants[kind].get("makedb_seconds")}s '
              f'hits={variants[kind].get("hits")} '
              f'staxids_rate={variants[kind].get("staxids_rate")}')
    report['s1_variants'] = variants
    good = [k for k, v in variants.items()
            if v.get('makedb_exit') == 0 and v.get('staxids_rate', 0) >= 0.99]
    report['s1_conclusion'] = {
        'usable_variants': good,
        'recommended': 'full-both' if 'full-both' in good else (
            good[0] if good else None),
        'note': ('full-both keeps both key columns identical to the full '
                 'sequence Id, which is the rule in plan appendix A.3; '
                 'rank names are mapped to the NCBI vocabulary because '
                 'DIAMOND rejects "domain"'),
    }

    # S2: BLAST route
    if args.with_blast:
        w = join(work, 's2_blast')
        makedirs(w, exist_ok=True)
        with open(join(w, 'seq2taxid.txt'), 'w') as fo, \
                open(args.taxonmap) as fi:
            next(fi)
            for line in fi:
                x = line.rstrip('\r\n').split('\t')
                if len(x) >= 3:
                    fo.write(f'{x[1]}\t{x[2]}\n')
        cmd = (f'makeblastdb -dbtype prot -in {db_faa} -out {join(w, "db")} '
               f'-parse_seqids -taxid_map {join(w, "seq2taxid.txt")}')
        ec, out = run_command(cmd)
        s2 = {'makeblastdb_exit': ec, 'tail': _tail(out, 3)}
        if ec == 0:
            cmd = (f'blastp -query {query} -db {join(w, "db")} -outfmt '
                   f'"6 qaccver saccver pident evalue bitscore qcovhsp '
                   f'staxids" -max_target_seqs 25 -num_threads '
                   f'{args.threads} 2>&1 | head -20')
            ec, out = run_command(cmd)
            rows = [l.split('\t') for l in _joined(out).splitlines()
                    if l.strip() and not l.startswith('WARNING')]
            s2['search_exit'] = ec
            s2['rows'] = len(rows)
            s2['with_staxids'] = sum(1 for r in rows
                                     if len(r) > 6 and r[6].strip())
            s2['example'] = rows[:2]
            cmd = (f'blastdbcmd -db {join(w, "db")} -entry all -outfmt '
                   f'"%a %T" 2>&1 | head -5')
            ec, out = run_command(cmd)
            s2['blastdbcmd_sample'] = out or []
        report['s2_blast'] = s2
        print(f'S2: makeblastdb={s2.get("makeblastdb_exit")} '
              f'rows={s2.get("rows")} with_staxids={s2.get("with_staxids")}')

    path = args.out or join(work, 'spike_s1_s2.json')
    path = gtdb.guard_write(path)
    with open(path, 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f'Wrote {path}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
