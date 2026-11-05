#!/usr/bin/env python
"""Re-verify the committed GTDB fixtures against the read-only release.

Why this exists: `neoHGT/tests/data/GTDB/fixtures.json` was written once by
tools/make_fixtures.py and CI has claimed ever since that the fixtures are
"KB-scale excerpts of GTDB r220 genomes". Recording a hash is not checking one.
Until now nothing re-read the source data, so any of these would have stayed
invisible: a fixture sequence that is not in the release at all (invented), a
sequence edited after the manifest was written, a fixture file added without a
provenance row, a source file that has changed under us (the release is supposed
to be immutable), or a committed fixture DB whose record count no longer matches
the number published in the manifest.

Safety: /data2 is only ever opened for reading; the report path is guarded to
stay inside the workspace.

    usage: python tools/verify_fixture_provenance.py --code <neoHGT repo> \
               [--gtdb-dir <release root>] [--out <report.json>]
"""
import argparse
import gzip
import hashlib
import json
import os
import sys
from os.path import abspath, isfile, join

GENOME_LIMIT = 10          # project rule: test copies never exceed 10 genomes


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_records(path):
    """{sequence id: residues} from a .faa, plain or gzipped."""
    out, key, buf = {}, None, []
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if key is not None:
                    out[key] = ''.join(buf)
                key, buf = line[1:].split()[0], []
            elif key is not None:
                buf.append(line.strip())
    if key is not None:
        out[key] = ''.join(buf)
    return out


def check(manifest_path, gtdb_dir=None):
    """(report, problems) for one manifest; the dataset is only read."""
    problems = []
    with open(manifest_path, encoding='utf-8') as f:
        man = json.load(f)
    root = os.path.dirname(abspath(manifest_path))
    recorded = man.get('gtdb_dir')
    # --gtdb-dir redirects the whole tree: a release mounted somewhere else is
    # still the same data, so the recorded paths are re-based rather than distrusted
    gtdb_dir = gtdb_dir or recorded
    rows, checked = [], 0

    def source_of(g):
        # the manifest stores the source either as an absolute path (older
        # manifests) or relative to the recorded release root; --gtdb-dir
        # re-bases either form onto the mounted release
        src = g.get('source_path')
        if not src:
            return None
        if not os.path.isabs(src) and recorded:
            src = join(recorded, src)
        if gtdb_dir and recorded and abspath(src).startswith(
                abspath(recorded) + os.sep):
            return join(gtdb_dir, os.path.relpath(abspath(src),
                                              abspath(recorded)))
        return src
    if len(man.get('genomes', [])) > GENOME_LIMIT:
        problems.append(f'清单里有 {len(man["genomes"])} 个基因组，超过项目规则'
                        f'允许的 {GENOME_LIMIT} 个测试拷贝上限')
    listed = set()
    for g in man.get('genomes', []):
        row = {'accession': g.get('accession'), 'label': g.get('label')}
        recorded_src = g.get('source_path')
        if not recorded_src:
            problems.append(f'{g.get("accession")} 没有记 source_path')
            rows.append(row)
            continue
        abs_recorded_src = (recorded_src if os.path.isabs(recorded_src)
                            else (join(recorded, recorded_src)
                                  if recorded else None))
        if recorded and abs_recorded_src and not abspath(abs_recorded_src).startswith(
                abspath(recorded) + os.sep):
            problems.append(f'{g.get("accession")} 的来源 {recorded_src} 不在清单'
                            f'声明的只读数据集 {recorded} 之内')
            rows.append(row)
            continue
        src = source_of(g)
        kb = g.get('kb_excerpt')
        if kb:
            # the excerpt is the committed repo file: the recorded absolute
            # path was written by make_fixtures on the original machine and
            # need not exist on another checkout, but the same file ships at
            # protein/<basename> relative to the manifest -- prefer that so
            # every checkout can account for its own fixtures
            local = join(root, 'protein', os.path.basename(kb))
            kb = local if isfile(local) else kb
        if kb and isfile(kb):
            listed.add(abspath(kb))
        copy = g.get('fixture_copy')
        if copy:
            # cache-relative in current manifests; resolve against the
            # workspace cache so other checkouts simply find no copy
            code = abspath(join(root, os.pardir, os.pardir, os.pardir,
                                os.pardir))
            cache = join(os.path.dirname(code), 'neoHGT2GTDB-缓存测试')
            copy_abs = copy if os.path.isabs(copy) else join(cache, copy)
        else:
            copy_abs = None
        if copy_abs and isfile(copy_abs):
            listed.add(abspath(copy_abs))
        if not isfile(src):
            problems.append(f'{g.get("accession")} 的来源文件已不存在：{src}')
            rows.append(row)
            continue
        checked += 1
        if sha256(src) != g.get('source_sha256'):
            problems.append(f'{g.get("accession")} 的来源 sha256 与清单不符 —— '
                            '只读数据集被改过，或者清单是抄来的')
        if os.path.getsize(src) != g.get('source_bytes'):
            problems.append(f'{g.get("accession")} 的来源字节数与清单不符')
        if g.get('domain') and f'/{g["domain"]}/' not in src.replace(os.sep, '/'):
            problems.append(f'{g.get("accession")} 记为 {g["domain"]}，'
                            f'但来源路径里没有这一档：{src}')
        if kb and isfile(kb):
            src_recs, kb_recs = read_records(src), read_records(kb)
            row['kb_records'] = len(kb_recs)
            missing = [i for i in kb_recs if i not in src_recs]
            altered = [i for i, s in kb_recs.items()
                       if i in src_recs and src_recs[i] != s]
            row['excerpt_of_source'] = not missing and not altered
            if missing or altered:
                problems.append(
                    f'{g.get("accession")} 的夹具里有来源中不存在或被改过的序列'
                    f'（缺 {len(missing)} 条 / 改 {len(altered)} 条，例如 '
                    f'{(missing or altered)[:1]}）—— "KB 级摘录"不成立')
            if len(kb_recs) >= len(src_recs):
                problems.append(f'{g.get("accession")} 的夹具不是"摘录"：'
                                f'{len(kb_recs)} 条 vs 来源 {len(src_recs)} 条')
        if copy_abs and isfile(copy_abs):
            if sha256(copy_abs) != g.get('fixture_sha256'):
                problems.append(f'{g.get("accession")} 的工作副本与清单记录的 '
                                'fixture_sha256 不符')
        rows.append(row)

    protein = join(root, 'protein')
    if os.path.isdir(protein):
        for name in sorted(os.listdir(protein)):
            if abspath(join(protein, name)) not in listed:
                problems.append(f'夹具 protein/{name} 在清单里没有对应条目'
                                '（来源未记录 = 无法核对）')
    db = join(root, 'db.faa')
    row_db = None
    if isfile(db) and 'derived' in man:
        n = len(read_records(db))
        if man['derived'].get('db_faa_records') not in (None, n):
            problems.append(
                f'清单说 db.faa 有 {man["derived"]["db_faa_records"]} 条，'
                f'实际 {n} 条')
        row_db = n
    return {'manifest': manifest_path,
            'genomes': len(man.get('genomes', [])),
            'sources_found': checked,
            'gtdb_dir': gtdb_dir, 'release': man.get('release'),
            'db_faa_records': row_db,
            'rows': rows, 'problems': problems}, problems


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--code', default=os.environ.get('NEOHGT_CODE'),
                    help='the adaptation repository')
    ap.add_argument('--manifest', default=None,
                    help='fixtures.json (default: <code>/neoHGT/tests/data/'
                         'GTDB/fixtures.json)')
    ap.add_argument('--gtdb-dir', default=os.environ.get('GTDB_DIR'),
                    help='override the release root recorded in the manifest')
    ap.add_argument('--out', default=None, help='JSON report inside the workspace')
    args = ap.parse_args()
    path = args.manifest or join(args.code, 'neoHGT', 'tests', 'data', 'GTDB',
                                 'fixtures.json')
    if not isfile(path):
        print(f'FAIL: no manifest at {path}')
        return 2
    report, problems = check(path, args.gtdb_dir)
    if args.out:
        root = os.environ.get('NEOHGT_ROOT')
        real = os.path.realpath(os.path.expanduser(args.out))
        if root and not real.startswith(abspath(root) + os.sep):
            raise SystemExit(f'refusing to write outside the workspace: {real}')
        os.makedirs(os.path.dirname(real), exist_ok=True)
        with open(real, 'w', encoding='utf-8') as f:
            json.dump(report, f, ensure_ascii=False, indent=1)
    print(f'manifest: {path}')
    print(f'genomes: {report["genomes"]} (limit {GENOME_LIMIT}); '
          f'sources read: {report["sources_found"]}; '
          f'db.faa records: {report["db_faa_records"]}')
    print('verdict: ' + ('PASS (every fixture sequence is traceable to the '
                         'read-only release)' if not problems else 'FAIL'))
    for p in problems:
        print('  PROBLEM', p)
    return 0 if not problems else 1


if __name__ == '__main__':
    sys.exit(main())
