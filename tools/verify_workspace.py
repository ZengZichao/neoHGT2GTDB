#!/usr/bin/env python3
"""IT-7 / DoD-9: check the workspace contract mechanically.

The project rules (plan §5.8, safety rules 1-6 and 12) make four promises that
were previously only argued in prose:

1. every file created by this project lives inside the three project
   directories, and nothing escapes them through a symlink;
2. the cache directory keeps the documented layout (db/, search/, analyze/,
   fixtures/, logs/, ...);
3. the temporary directory is *disposable*: no built database or manifest
   refers to a path inside it, so clearing it cannot invalidate a result;
4. the read-only GTDB dataset is untouched (mtime audit, unless skipped).

usage:
    python tools/verify_workspace.py --out <cache>/reports/IT7_workspace.json
    # optional: --since "2026-11-01 17:25:00" to include the dataset audit
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from os.path import abspath, expanduser, isdir, isfile, join, realpath

sys.path.insert(0, os.path.dirname(os.path.dirname(abspath(__file__))))

PROJECT_DIRS = ('neoHGT2GTDB-项目代码', 'neoHGT2GTDB-临时文件',
                'neoHGT2GTDB-缓存测试')
# Pre-existing inputs that are not produced by this project.
ROOT_ALLOW = ('neoHGT-main', 'neoHGT-main.zip')
CACHE_LAYOUT = ('db', 'search', 'analyze', 'fixtures', 'logs', 'reports',
                'benchmarks', 'input')
MANIFESTS = ('meta.json', join('gtdb', 'extraction_state.json'))


def outside(root, path):
    """True when realpath(path) escapes root."""
    return not realpath(path).startswith(root + os.sep)


def check_project_files(root, cutoff=None):
    """(files_seen, unexpected_root_entries, symlink_escapes).

    A root entry that is neither a project directory nor one of the plan's own
    documents is reported. When a cutoff is given, entries that predate the
    project are reported separately as pre-existing (this project did not make
    them, so they are not a workspace violation).
    """
    seen, unexpected, links = 0, [], []
    for name in sorted(os.listdir(root)):
        if name in PROJECT_DIRS or name in ROOT_ALLOW:
            continue
        path = join(root, name)
        if isdir(path) or not name.endswith(('.md', '.pdf')):
            st = os.stat(path)
            entry = {'path': path + ('/' if isdir(path) else ''),
                     'mtime': time.strftime('%Y-%m-%d %H:%M:%S',
                                            time.localtime(st.st_mtime))}
            entry['predates_project'] = (
                st.st_mtime < cutoff if cutoff is not None else False)
            unexpected.append(entry)
    for name in PROJECT_DIRS:
        base = join(root, name)
        for dirpath, dirnames, names in os.walk(base, followlinks=False):
            for d in list(dirnames):
                p = join(dirpath, d)
                if os.path.islink(p) and outside(root, p):
                    links.append(p)
            for f in names:
                p = join(dirpath, f)
                seen += 1
                if os.path.islink(p) and outside(root, p):
                    links.append(p)
    return seen, unexpected, links


def check_cache_layout(cache):
    return [d for d in CACHE_LAYOUT if not isdir(join(cache, d))]


# Scratch locations recorded inside manifests are an inventory, not a verdict:
# "the temp dir is disposable" is proven behaviorally (verify_resume rerun with a
# brand-new empty --tmpdir, reports/IT7_tmp_disposable_probe), because a path
# string inside a manifest does not by itself say whether anything reads it back.
def check_tmp_disposable(root, cache, tmp):
    """Paths inside the temp dir referenced by a database manifest."""
    refs = []
    tmp_real = realpath(tmp)
    dbs = join(cache, 'db')
    if not isdir(dbs):
        return refs
    for db in sorted(os.listdir(dbs)):
        root_db = join(dbs, db)
        if not isdir(root_db):
            continue
        for rel in MANIFESTS:
            path = join(root_db, rel)
            if not isfile(path):
                continue
            try:
                with open(path) as f:
                    blob = f.read()
            except (OSError, UnicodeDecodeError):
                continue
            for key, token in _strings(blob):
                if not token.startswith('/'):
                    continue
                norm = realpath(token)
                if norm != tmp_real and not norm.startswith(tmp_real + os.sep):
                    continue
                refs.append({'db': db, 'manifest': rel, 'key': key,
                             'path': token,
                             'exists_now': isfile(token) or isdir(token)})
    return refs


def _strings(blob):
    """(json-path, value) pairs for every string in a JSON blob.

    Raw text that does not parse yields ('', token) instead.
    """
    out = []

    def walk(x, p):
        if isinstance(x, dict):
            for k, v in x.items():
                walk(v, f'{p}/{k}')
        elif isinstance(x, list):
            for i, v in enumerate(x):
                walk(v, f'{p}[{i}]')
        elif isinstance(x, str):
            out.append((p, x))
    try:
        walk(json.loads(blob), '')
    except ValueError:
        out.extend(('', t) for t in blob.split())
    return out


OWN_CACHE_NAMESPACES = (
    'matplotlib', 'mpl', 'pip', 'ipython', 'numba', 'Pillow', 'torch',
    'huggingface', 'fontTools', 'jupyter', 'uv', 'pandas')

# A namespace our stack uses is normally ours to redirect. Sometimes it is not:
# round 52's gates failed on a fresh `pip/selfcheck` file -- pip's "is there a new
# pip" cache, written by *whatever* process on this machine invokes that
# interpreter -- while four gate-shaped commands (pytest, the register, importing
# every module, the fixture comparison) produced zero writes there and no command
# in commands.log has ever named pip. The evidence is a file, and this list points
# at it: the exemption is only granted while the measured claim still holds, and
# only for the exact sub-path, never for a whole namespace.
DECLARED_FOREIGN_WRITERS = (
    {'path': 'pip/selfcheck',
     # a stable name on purpose: round 52 named the evidence after itself
     # (R52_...), so refreshing the measurement meant editing this declaration
     # every time -- friction that turns an exemption into a stale claim. The
     # probe writes here, and the file says which code state it measured.
     'evidence': 'neoHGT2GTDB-缓存测试/reports/IT7_pip_attribution.json',
     'must_say': {'our_stack_wrote': 0, 'commands_log_pip_invocations': 0},
     'probe': 'neoHGT2GTDB-临时文件/r119_pip_attribution.py',
     'why': 'another process sharing /opt/anaconda3 writes pip\'s version '
            'self-check; measured that our own commands write nothing there'},
)

# `pip ...`, `pip3 ...`, `python -m pip ...`, in a command field only.
PIP_CMD = re.compile(r'(?:^|\s)pip3?\s|-m\s+pip(?=\s|$)')

# What could make that measurement wrong? Our own code: a new importer of
# `neoHGT/server.py` (which pip-installs at import time), a test that shells out,
# a dropped `PIP_CACHE_DIR` in tools/env.sh. So the evidence is only usable while
# it is at least as new as every file the gates execute. Round 52 measured that
# this is not hypothetical: three minutes after the probe wrote its report, the
# session had edited `verify_workspace.py`, `tools/env.sh` and the test file, so a
# static exemption would have been defending a claim about code that no longer
# exists.
EVIDENCE_SCOPE = (
    ('neoHGT2GTDB-项目代码/neoHGT', ('.py',)),
    ('neoHGT2GTDB-项目代码/neoHGT/tests', ('.py',)),
    ('neoHGT2GTDB-项目代码/tools', ('.py', '.sh')),
)


def pip_command_hits(log_path):
    """Command lines in commands.log that invoked pip. The prose fields are not evidence.

    Round 52's probe searched whole log lines and matched the *intent* text of a
    gate run ("... env.sh pip redirection, 4 new tests ..."), which reported a pip
    invocation that never happened -- a check that reads its own labels as
    evidence will confirm whatever the author typed. A command line is the last
    tab-separated field of an entry; that field alone is searched.
    """
    hits = []
    if not isfile(log_path):
        return hits
    with open(log_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line.startswith('20'):
                continue                      # continuation rows: OUT/MOD/GONE/EXTERNAL
            cmd = line.rstrip('\n').split('\t')[-1]
            if PIP_CMD.search(cmd):
                hits.append(cmd[:120])
    return hits


def newest_in_scope(root):
    """(mtime, relative path) for the newest file that could invalidate a measurement."""
    best = (0.0, None)
    for rel, exts in EVIDENCE_SCOPE:
        base = join(root, rel)
        if not isdir(base):
            continue
        for name in sorted(os.listdir(base)):
            fp = join(base, name)
            if not isfile(fp) or os.path.splitext(name)[1] not in exts:
                continue
            try:
                mtime = os.stat(fp).st_mtime
            except OSError:
                continue
            if mtime > best[0]:
                best = (mtime, f'{rel}/{name}')
    return best


def evidence_state(root, decl):
    """'fresh' / 'stale' / 'no evidence' / 'no scope' for one declaration.

    The evidence's own statement about which code it measured wins when present:
    a report's mtime can be moved by a copy or a touch, while the field is part of
    what the probe observed. The file mtime is the fallback for evidence written
    before that field existed -- and the fallback is labelled, not silently used.
    """
    ev = join(root, decl['evidence'])
    if not isfile(ev):
        return 'no evidence', decl['evidence']
    newest, path = newest_in_scope(root)
    if path is None:
        # an exemption whose scope contains nothing was never measured over anything
        return 'no scope', 'scope 内没有任何被执行的代码文件'
    measured = None
    try:
        with open(ev, encoding='utf-8') as fh:
            measured = json.load(fh).get('newest_measured_code_mtime')
    except (OSError, ValueError, AttributeError):
        measured = None
    if measured is not None:
        stated = True
        ev_mtime = float(measured)
    else:
        stated = False
        ev_mtime = os.stat(ev).st_mtime
    if ev_mtime < newest:
        return 'stale', f'{path}' + ('' if stated else ' (evidence mtime, not self-stated)')
    return 'fresh', path


def declared_writer(rel, root):
    """(declaration, ok, reason) for a relative path under a watched prefix."""
    norm = rel.replace(os.sep, '/')
    for d in DECLARED_FOREIGN_WRITERS:
        if norm.startswith(d['path'] + '/') or norm == d['path']:
            if len(d['path'].split('/')) < 2:
                return d, False, 'whole-namespace'   # a whole-namespace exemption is not one
            import json
            try:
                with open(join(root, d['evidence']), encoding='utf-8') as fh:
                    got = json.load(fh)
            except (OSError, ValueError):
                return d, False, 'unreadable evidence'
            if not all(got.get(k) == v for k, v in d['must_say'].items()):
                return d, False, 'evidence says otherwise'
            state, detail = evidence_state(root, d)
            if state != 'fresh':
                return d, False, f'{state}: {detail}'
            return d, True, 'fresh'
    return None, False, ''
def check_external_writes(root, cache, cutoff, quiet=None, noisy=None,
                          own_namespaces=None):
    """Read the runner's external-write report back, instead of trusting it.

    `tools/run.sh` scans $NEOHGT_WATCH after every command and appends an
    `EXTERNAL` line to `缓存测试/logs/commands.log` (round 39). Until now nothing
    ever read those lines: the one enforcement point for safety rule 11 produced
    a log that no gate consumed, and a command run *outside* the runner left no
    trace at all -- round 42's `/tmp/frag.sh` slip is missing from the log for
    exactly that reason. So this check has two halves of different strength:

    1. the log, which is authoritative for anything that went through the runner;
    2. a direct sample of the watched prefixes, which is the only thing that can
       see a runner-bypassing write.

    Attribution is decided by *namespace*, not by mtime. Round 44 measured that
    our own toolchain writes nothing into $HOME/.cache when `tools/env.sh` is
    sourced (matplotlib PDFs + `import neoHGT.analyze`: 0 new files), while a
    `fontconfig` cache file from another application on this machine appeared
    there inside a 24-hour window. So a new file under a cache namespace our
    stack actually uses (matplotlib, pip, ipython, ...) is ours and fails the
    gate; a new file under someone else's namespace is reported as observed and
    not attributable -- otherwise the gate fails on this desktop session and
    gets ignored, which is worse than not having it. /tmp and /dev/shm are never
    judged: round 43 counted 215,824 and 11,894 entries, the recent ones in /tmp
    belonging to other processes and `/dev/shm/__KMP_REGISTERED_LIB_*` created by
    our own OpenMP runtime.
    """
    quiet = quiet or ('~/.cache', '~/.config', '~/.matplotlib', '/var/tmp')
    noisy = noisy or ('/tmp', '/dev/shm')
    own = own_namespaces or OWN_CACHE_NAMESPACES
    if cutoff is None:
        cutoff, window = time.time() - 86400.0, '24h default window'
    else:
        window = 'round window'
    log = join(cache, 'logs', 'commands.log')
    out = {'log': log, 'log_exists': isfile(log), 'commands': 0,
           'external_lines': 0, 'external_in_workspace': 0,
           'external_outside': [], 'quiet': {}, 'noisy': {},
           'window': window, 'own_namespaces': list(own)}
    problems = []
    if out['log_exists']:
        with open(log, encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('20'):
                    out['commands'] += 1
                if '\tEXTERNAL\t' not in line:
                    continue
                path = line.split('\tEXTERNAL\t')[-1].strip()
                out['external_lines'] += 1
                if path.startswith(root + os.sep):
                    out['external_in_workspace'] += 1
                else:
                    out['external_outside'].append(path)
    else:
        problems.append(f'{log} 不存在：这条核对什么都没看')
    if out['log_exists'] and not out['commands']:
        problems.append('commands.log 里一条命令都没有：日志被清空或格式变了')
    problems.extend(f'工作区外的写入: {p}' for p in out['external_outside'])
    for prefix in quiet:
        p = expanduser(prefix)
        if not isdir(p):
            out['quiet'][prefix] = {'exists': False, 'entries': 0,
                                    'recent': 0, 'n_ours': 0, 'n_other': 0,
                                    'recent_ours': [], 'recent_other': []}
            continue
        entries, recent_ours, recent_other = 0, [], []
        recent_declared, declared_unproven = [], []
        for base, dirs, files in os.walk(p):
            if base[len(p):].count(os.sep) >= 2:
                dirs[:] = []
            for name in files:
                entries += 1
                fp = os.path.join(base, name)
                try:
                    mtime = os.stat(fp).st_mtime
                except OSError:
                    continue
                if mtime < cutoff:
                    continue
                rel = fp[len(p):].strip(os.sep)
                d, ok, why = declared_writer(rel, root)
                if d and ok:
                    recent_declared.append(fp)
                elif d:
                    declared_unproven.append((fp, why))
                    recent_ours.append(fp)
                elif any(ns in rel for ns in own):
                    recent_ours.append(fp)
                else:
                    recent_other.append(fp)
        out['quiet'][prefix] = {
            'exists': True, 'entries': entries,
            'recent': len(recent_ours) + len(recent_other) + len(recent_declared),
            'n_ours': len(recent_ours), 'n_other': len(recent_other),
            'n_declared_foreign': len(recent_declared),
            'recent_ours': recent_ours[:5], 'recent_other': recent_other[:5],
            'recent_declared': [f'{fp} <- {declared_writer(fp[len(p):].strip(os.sep), root)[0]["evidence"]}'
                                for fp in recent_declared[:5]],
            'declared_without_evidence': [f'{fp} [{why}]' for fp, why
                                          in declared_unproven[:5]]}
        if declared_unproven:
            # the message has to say what to run, or a red gate is a riddle: a
            # stale measurement is fixed by re-measuring, not by deleting the row
            probe = DECLARED_FOREIGN_WRITERS[0]['probe']
            problems.append(
                f'{prefix} 里 {len(declared_unproven)} 个文件走的是"已声明的外部写入者"'
                f'这条豁免，但支撑它的实测不成立了：'
                + ', '.join(f'{os.path.relpath(fp, p)} [{why}]'
                            for fp, why in declared_unproven[:3])
                + f' —— 豁免必须有证据且证据不能比被测量的代码旧；'
                  f'重跑 {probe} 即可刷新（它只读 ~/.cache，不改代码）')
        if recent_ours:
            problems.append(
                f'{prefix} 里有本工具链该重定向却没重定向的缓存 '
                f'{len(recent_ours)} 个（{window}内新建）：'
                + ', '.join(os.path.relpath(x, p) for x in recent_ours[:3]))
    for prefix in noisy:
        # Not walked: round 43/44 measured these (see the register), and walking
        # 200k+ entries on every IT-7 run to print "not judged" costs seconds and
        # makes the report's numbers move with unrelated machine activity. The
        # caveat is what matters, not the count.
        out['noisy'][prefix] = {'judged': False, 'measured_elsewhere': True,
                                'note': 'shared with other processes (round 43: '
                                        '/tmp 215,828 entries, 25 recent from '
                                        'unrelated programs) or created by our own '
                                        'OpenMP runtime (/dev/shm '
                                        '__KMP_REGISTERED_LIB_*); not attributable'}
    return out, problems


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', default=os.environ.get('NEOHGT_ROOT'))
    ap.add_argument('--cache', default=os.environ.get('NEOHGT_CACHE'))
    ap.add_argument('--tmp', default=os.environ.get('NEOHGT_TMP'))
    ap.add_argument('--code', default=os.environ.get('NEOHGT_CODE'))
    ap.add_argument('--gtdb-dir', default=os.environ.get('GTDB_DIR'))
    ap.add_argument('--since', help='UTC cutoff for the dataset mtime audit')
    ap.add_argument('--out', required=True)
    ap.add_argument('--probe-report', action='append', default=[],
                    help='path of a temp-disposability probe (verify_resume '
                         'rerun with a fresh --tmpdir); repeatable, one per '
                         'database source, and every one must PASS')
    args = ap.parse_args()

    root = realpath(abspath(args.root))
    cutoff = None
    if args.since:
        from calendar import timegm
        cutoff = timegm(time.strptime(args.since, '%Y-%m-%d %H:%M:%S'))
    seen, at_root, links = check_project_files(root, cutoff)
    missing = check_cache_layout(args.cache)
    refs = check_tmp_disposable(root, args.cache, args.tmp)
    violations = [e for e in at_root if not e.get('predates_project')]
    external, external_problems = check_external_writes(root, args.cache,
                                                        cutoff)

    report = {'root': root, 'project_files': seen,
              'unexpected_at_root': at_root, 'symlink_escapes': links,
              'cache_layout_missing': missing,
              'temp_paths_in_manifests': refs,
              'external_writes': external,
              'external_write_problems': external_problems,
              'temp_disposability_probes': list(args.probe_report)}
    out = abspath(args.out)
    if args.since and args.gtdb_dir:
        # reuse the audited mtime implementation rather than duplicating it
        ds_out = os.path.splitext(out)[0] + '_dataset.txt'
        rc = subprocess.run(
            [sys.executable, join(os.path.dirname(abspath(__file__)),
                                  'verify_read_only.py'),
             '--dir', args.gtdb_dir, '--since', args.since,
             '--out', ds_out], capture_output=True, text=True).returncode
        report['dataset_audit'] = ds_out
        report['dataset_exit_code'] = rc
    else:
        report['dataset_exit_code'] = 'skipped (no --since/--gtdb-dir)'

    probes = {}
    for p in args.probe_report:
        if not isfile(p):
            probes[p] = 'missing'
            continue
        try:
            with open(p) as f:
                probes[p] = json.load(f).get('verdict', 'no verdict key')
        except (OSError, ValueError) as e:
            probes[p] = f'unreadable ({e})'
    report['probe_verdicts'] = probes
    probe_ok = bool(probes) and all(v == 'PASS' for v in probes.values())
    report['verdict'] = ('PASS' if not violations and not links and not missing
                         and not external_problems and probe_ok
                         and report['dataset_exit_code'] in
                         (0, 'skipped (no --since/--gtdb-dir)')
                         else 'FAIL')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(report, f, indent=1, ensure_ascii=False)
    text = (
        f'project files scanned: {seen:,} under {root}\n'
        f'root entries outside the three project dirs: {len(at_root)} '
        f'({len(violations)} not predating the project)\n')
    for e in at_root:
        text += f'  {e["path"]}  mtime={e["mtime"]}  ' \
                f'predates_project={e.get("predates_project")}\n'
    text += (
        f'symlinks escaping the workspace: {links or "(none)"}\n'
        f'missing cache subdirectories: {missing or "(none)"}\n'
        f'scratch paths recorded inside manifests: {len(refs)} '
        '(inventory; disposability is proven by the probes)\n')
    for p, v in probes.items():
        text += f'temp-disposability probe: {p} -> {v}\n'
    if not probes:
        text += 'temp-disposability probe: none supplied (--probe-report)\n'
    text += (
        f'commands.log: {external["commands"]} command entries, '
        f'{external["external_lines"]} EXTERNAL lines '
        f'({external["external_in_workspace"]} inside the workspace, '
        f'{len(external["external_outside"])} outside)\n'
        f'quiet prefixes sampled ({external["window"]}): '
        + ', '.join(f'{k}: {v.get("n_ours", 0)} ours / '
                    f'{v.get("n_other", 0)} other-namespace'
                    for k, v in sorted(external['quiet'].items())) + '\n'
        + 'not judged (see round 43 measurement in the register): '
        + ', '.join(f'{k}: {v["note"][:40]}…'
                    for k, v in sorted(external['noisy'].items())) + '\n')
    for problem in external_problems:
        text += f'  EXTERNAL PROBLEM {problem}\n'
    text += (
        f'dataset audit exit code: {report["dataset_exit_code"]}\n'
        f'verdict: {report["verdict"]}\n')
    with open(os.path.splitext(out)[0] + '.txt', 'w') as f:
        f.write('# IT-7 workspace compliance\n' + text)
    print(text)
    print(f'wrote {out}')
    return 0 if report['verdict'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
