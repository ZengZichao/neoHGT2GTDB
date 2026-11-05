#!/usr/bin/env python3
"""Traceability gate: every verification id in the plan must have evidence.

The plan names 30+ test identifiers (UT-*, IT-*, AT-*, RT-*). Prose claiming
that they are covered is not checkable, and two real gaps were previously found
by hand (AT-4 had no test at all; UT-0.3/UT-2.1/UT-3.3/UT-3.4 had no labelled
test). This tool makes the sweep repeatable: an id counts as traced when either

  * a test file cites the id, or
  * the id maps to an evidence artifact that actually exists on disk.

Anything else is a gap and the tool exits non-zero.

usage:
    python tools/verify_plan_traceability.py \
        --plan "$NEOHGT_ROOT/neoHGT2.GTDB-ADAPTATION-PLAN-CHN.md" \
        --out  "$NEOHGT_CACHE/reports/plan_traceability.json"
"""

import argparse
import ast
import json
import os
import re
import sys
from os.path import abspath, dirname, exists, isdir, isfile, join, splitext, basename

ID_RE = re.compile(r'\b(UT-\d+\.\d+(?:b|c)?|UT-\d+\.x|IT-\d+[ab]?|'
                   r'AT-\d+|RT-\d+)\b')

# ids whose evidence is an artifact rather than a labelled unit test
EVIDENCE = {
    'IT-0': 'logs/layout_probe.json',
    'IT-1': 'benchmarks/it1_taxdump.time',
    'IT-2': 'db/it2_smoke/meta.json',
    'IT-2b': 'db/r220_ar_family/meta.json',
    'IT-3': 'search/halo/hgt_demo.json',
    'IT-4': 'search/halo/hgt_demo.json',
    'IT-5': 'reports/IT5_hybrid_r3.json',
    'IT-6': 'reports/IT6_resume.json',
    'IT-7': 'reports/IT7_workspace.json',
    'RT-1': 'regression/RT1_failure_set_diff.txt',
    'RT-2': 'reports/RT2_defaults.json',
}
# grouped ids (UT-1.x, UT-3.x, ...) stand for their whole family: satisfied if
# any concrete member is traced
GROUP_SUFFIX = '.x'


def cite_map(test_dir):
    """{id: [test file names citing it]} from the test sources."""
    cites = {}
    for name in sorted(os.listdir(test_dir)):
        if not (name.startswith('test_') and name.endswith('.py')):
            continue
        with open(join(test_dir, name), encoding='utf-8') as f:
            text = f.read()
        for i in set(ID_RE.findall(text)):
            cites.setdefault(i, []).append(name)
    return cites


CLI_STAGES = ('database', 'search', 'analyze', 'gtdb')


def _read_text(path):
    with open(path, encoding='utf-8', errors='replace') as f:
        return f.read()


HEADING = re.compile(r'(#{1,6})\s')


def doc_section(text, start, stop=None):
    """The body of the heading named `start`, or '' when no heading names it.

    Anchored on the heading *line*, not on the token's first mention. A locator
    that splits on the first mention reads whatever sits between that mention and
    the next boundary, so a sentence that merely repeats a section title relocates
    the slice: round 54's pass C put 「动态配置段的可达性」 into a deliverables row
    338 lines above its own heading, and rule 23 then reported the dispositions as
    missing while rule 22 re-reported the three `precomp` keys they cover.

    `stop` names the heading that ends the section (used by the plan, whose
    appendices are subdivided at the same level); without it the section ends at
    the next heading of the same or a shallower level.
    """
    lines = text.splitlines(True)
    start_i = level = None
    for i, ln in enumerate(lines):
        m = HEADING.match(ln)
        if m and start in ln:
            start_i, level = i + 1, len(m.group(1))
            break
    if start_i is None:
        return ''
    for j in range(start_i, len(lines)):
        m = HEADING.match(lines[j])
        if not m:
            continue
        if stop and stop in lines[j]:
            return ''.join(lines[start_i:j])
        if not stop and len(m.group(1)) <= level:
            return ''.join(lines[start_i:j])
    return ''.join(lines[start_i:])


def option_tables(tree):
    """Every `arguments` list in a CLI module, name aliases resolved.

    Read with ast instead of imported: the gate has to run on a checkout where
    the scientific dependencies are absent, and these tables are plain data.
    """
    named, lists = {}, []
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.List):
            for t in node.targets:
                if isinstance(t, ast.Name):
                    named[t.id] = node.value
                    if t.id == 'arguments':
                        lists.append(node.value)
        elif (isinstance(node, ast.AugAssign)
              and isinstance(node.target, ast.Name)
              and node.target.id == 'arguments'):
            val = node.value
            if isinstance(val, ast.Name):
                val = named.get(val.id)
            if isinstance(val, ast.List):
                lists.append(val)
    return lists


def table_flags(lists):
    """Option strings of an `arguments` table; `-o|--output` counts twice."""
    flags = set()
    for lst in lists:
        for item in lst.elts:
            if (isinstance(item, (ast.List, ast.Tuple)) and item.elts
                    and isinstance(item.elts[0], ast.Constant)):
                flags.update(str(item.elts[0].value).split('|'))
    return {f for f in flags if f.startswith('-')}


def cli_flags(code_dir):
    """{module: set(flags)} of every CLI surface the product actually accepts.

    Three stage modules declare an `arguments` table; `neoHGT gtdb` has both a
    table (the legacy flat form) and a hand-built argparse tree (the subcommands),
    so both are read. Round 33's lesson: this function used to grep the gtdb
    source for `add_argument('--x')`, which structurally cannot see a table entry
    -- `--action`, `--dir`, `--db`, `--in`, `--to`, `--deep`, `--threads`,
    `--tmpdir` and `--output` were invisible to every rule built on this set,
    including "is this planned option handed to the parser by a test?".
    """
    out = {}
    for name in CLI_STAGES:
        path = join(code_dir, 'neoHGT', f'{name}.py')
        if not isfile(path):
            out[name] = set()
            continue
        with open(path, encoding='utf-8', errors='replace') as f:
            src = f.read()
        flags = table_flags(option_tables(ast.parse(src)))
        if name == 'gtdb':
            flags |= set(re.findall(r"add_argument\(\s*'([^']+)'", src))
        out[name] = {f for f in flags if f.startswith('-')}
    return out


def flag_dest(flag):
    return flag[2:].replace('-', '_')


# --- plan §6 numbered tasks ------------------------------------------------
# The plan's phase tables enumerate 28 numbered tasks (S0-S3, 2.1-2.13,
# 3.1-3.6, 4.1-4.5). docs/计划任务对照表.md is the ledger that says, for each of
# them, where the implementation lives and what proves it. This gate keeps the
# ledger honest: one row per plan id, a status from a fixed vocabulary, cited
# paths that exist, and an explicit reason whenever the status is not 已实现.
TASK_ID_RE = re.compile(r'^\|\s*(S\d|[1-4]\.\d+)\s*\|', re.M)


def row_cells(line):
    """Cells of a markdown table row; a pipe escaped as `\\|` stays inside."""
    return [c.strip() for c in re.split(r'(?<!\\)\|', line.strip().strip('|'))]
LEDGER_STATUSES = ('已实现', '部分实现', '未做')
# plan appendix F registers pre-existing defects that the adaptation shares a
# root cause with; the first status is the only one that needs no explanation
DEFECT_STATUSES = ('已修复', '已按方案更正', '登记未改')
DEFECT_ID_RE = re.compile(r'^\|\s*(F\d{1,2})\s*\|', re.M)
CODE_PREFIXES = ('neoHGT/', 'tools/', 'docs/', 'scripts/', 'ci/', '.github/',
                 'README')
ROOT_PREFIXES = ('缓存测试/', '临时文件/')
NOT_PATH = set('*<>{}|()')
# a bare file name in a citation is looked up in the package directory
SOURCE_EXTS = ('.py', '.yml', '.md', '.sh', '.txt', '.json', '.tsv')
MIN_REASON = 20
WORKSPACE_DIRS = {'缓存测试': 'neoHGT2GTDB-缓存测试',
                  '临时文件': 'neoHGT2GTDB-临时文件',
                  '项目代码': 'neoHGT2GTDB-项目代码'}

# §12 deliverables: same shape as the other two registers
DELIVER_STATUSES = ('已交付', '部分交付', '未交付')
DELIVER_ID_RE = re.compile(r'^\|\s*(D\d{1,2})\s*\|', re.M)
# §9 risk register: did each planned mitigation actually land in the code?
RISK_STATUSES = ('已实现缓解', '部分实现', '登记接受')
RISK_ID_RE = re.compile(r'^\|\s*(RK\d{1,2})\s*\|', re.M)
# §1.4 / §10.7 acceptance criteria: are the Definition-of-Done items really met?
DOD_STATUSES = ('达标', '部分达标', '被外部条件阻塞')
DOD_ID_RE = re.compile(r'(?<![A-Za-z-])DoD-(\d+)\b')
# §8 architecture decisions: is the shipped design still the decided one?
DECISION_STATUSES = ('已落实', '部分落实', '偏离并说明', '未落实')
ADR_ID_RE = re.compile(r'^\|\s*(ADR-\d+b?)\s*\|', re.M)
# §10.2 / §10.3 / §10.5: the plan's own verification matrix. A row here is
# only useful if it names the covering test (with a line that exists) or the
# evidence artifact, so "we tested it" cannot be claimed vaguely.
VERIF_STATUSES = ('已覆盖', '部分覆盖', '未覆盖')
UT_ID_RE = re.compile(r'^\|\s*(UT-[0-9]+\.[0-9]+b?|UT-[0-9]+\.x)\s*\|', re.M)
IT_ID_RE = re.compile(r'^\|\s*(IT-[0-9]+[a-z]?)\s*\|', re.M)
RT_ID_RE = re.compile(r'(?<![A-Za-z0-9-])(RT-[0-9]+)\b')
# §10.4 is an ordered list, not a table: its four items are the plan's
# "is the GTDB mode actually more trustworthy" evidence list.
SEC104_ID_RE = re.compile(r'^([1-9])\. ', re.M)
# §2.3 lists one row per assumption the search layer made about NCBI data
# (S1..S9), and §0.5 lists the review revisions (REV-01..REV-19). Both are
# numbered and both are truth-apt, so both belong behind a gate.
ASSUMPTION_ID_RE = re.compile(r'^\|\s*(S\d)\s*\|', re.M)
REVISION_ID_RE = re.compile(r'^\|\s*(REV-\d{2})\s*\|', re.M)
ASSUMPTION_STATUSES = ('已按GTDB改造', '确认无需改动', '与方案预期不同并已调整',
                       '未处理')
# §13 milestones: the plan's own "are we done" view (M0..M7 keyed by week row)
MILESTONE_ID_RE = re.compile(r'^\|\s*(W\d+(?:[–-]W?\d+)?\+?)\s*\|',
                             re.M)
MILESTONE_STATUSES = ('已达成', '部分达成', '未达成')
REVISION_STATUSES = ('已采纳', '部分采纳', '未采纳')


def plan_adr_ids(plan_text):
    """ADR ids of the plan's §8 decision table."""
    if '## 8. 关键设计决策' not in plan_text:
        return []
    sect = plan_text.split('## 8. 关键设计决策', 1)[1].split('## 9.', 1)[0]
    found = set(ADR_ID_RE.findall(sect))

    def order(x):
        # 'ADR-4' and 'ADR-4b' used to compare equal, which let sorted() fall
        # back to set-iteration order -- hash-randomised between processes, and
        # it made the frozen digest flap. Give the suffix a total order.
        m = re.match(r'ADR-(\d+)(.*)$', x)
        return (int(m.group(1)), m.group(2))

    return sorted(found, key=order)


def plan_verification_ids(plan_text):
    """Verification ids of plan §10.2 (UT), §10.3 (IT), §10.4 (科学有效性)
    and §10.5 (RT), in plan order.

    §10.4 is an ordered list, so its items are named `10.4-<n>` here; without a
    stable name they could not be registered or gated at all.
    """
    found = []
    for start, stop, rx, prefix in (('### 10.2', '### 10.3', UT_ID_RE, ''),
                                    ('### 10.3', '### 10.4', IT_ID_RE, ''),
                                    ('### 10.4', '### 10.5', SEC104_ID_RE,
                                     '10.4-'),
                                    ('### 10.5', '### 10.6', RT_ID_RE, '')):
        if start not in plan_text:
            continue
        sect = plan_text.split(start, 1)[1].split(stop, 1)[0]
        found += [prefix + x for x in rx.findall(sect)
                  if prefix + x not in found]
    return found


def plan_milestone_ids(plan_text):
    """W1..W9+ row keys of the plan's §13 milestone table."""
    if '## 13. 里程碑' not in plan_text:
        return []
    sect = plan_text.split('## 13. 里程碑', 1)[1].split('## 14.', 1)[0]
    found = set(MILESTONE_ID_RE.findall(sect))

    def order(x):
        head = re.match(r'W(\d+)', x)
        return (int(head.group(1)), x)

    return sorted(found, key=order)


def plan_assumption_ids(plan_text):
    """S1..S9 of the plan's §2.3 assumption table."""
    if '### 2.3' not in plan_text:
        return []
    sect = plan_text.split('### 2.3', 1)[1].split('### 2.4', 1)[0]
    found = set(ASSUMPTION_ID_RE.findall(sect))
    return sorted(found, key=lambda x: int(x[1:]))


def plan_revision_ids(plan_text):
    """REV-01.. of the plan's §0.5 review-revision table."""
    if '## 0.5' not in plan_text:
        return []
    sect = plan_text.split('## 0.5', 1)[1].split('## 1.', 1)[0]
    found = set(REVISION_ID_RE.findall(sect))
    return sorted(found)


def plan_dod_ids(plan_text):
    """DoD ids named by the plan's §1.4 table and §10.7 checklist."""
    found = set()
    for start, stop in (('### 1.4', '\n## 2.'), ('### 10.7', '\n## 11.')):
        if start not in plan_text:
            continue
        sect = plan_text.split(start, 1)[1].split(stop, 1)[0]
        found.update(int(x) for x in DOD_ID_RE.findall(sect))
    return [f'DoD-{i}' for i in sorted(found)]

# Options that belong to the tools the plan's example commands call, not to
# neoHGT itself; appendix C/D mixes both, so they need an explicit list rather
# than a guess.
THIRD_PARTY_FLAGS = {
    # diamond
    '--outfmt', '--sensitive', '--max-target-seqs', '--query', '--db',
    '--no-parse-seqids', '--taxonmap', '--threads', '--diamond', '--blastp',
    '--id', '--evalue', '--blastdbcmd', '--dbcmd',
    # blast+
    '--parse_seqids', '--taxid_map', '--max_target_seqs', '--out',
    # rsync / tar / coreutils used by the transport spikes
    '--list-only', '--wildcards', '--recursive', '--links', '--times',
    '--contimeout', '--partial', '--progress', '-L', '-t', '-q',
    # network fetches in the fetch plan
    '--location', '--retry', '--no-check-certificate',
    # python/pytest invocations in the reproduction appendix
    '--tb', '-q', '-k',
}


def plan_risk_ids(plan_text):
    """RK-ids of the plan's §9 risk register."""
    if '## 9. 风险登记册' not in plan_text:
        return []
    sect = plan_text.split('## 9. 风险登记册', 1)[1].split('## 10.', 1)[0]
    found = set(RISK_ID_RE.findall(sect))
    return sorted(found, key=lambda x: int(x[2:]))


def plan_deliverable_ids(plan_text):
    """D-ids of the plan's §12 deliverable list."""
    if '## 12. 交付物清单' not in plan_text:
        return []
    sect = plan_text.split('## 12. 交付物清单', 1)[1].split('## 13.', 1)[0]
    found = set(DELIVER_ID_RE.findall(sect))
    return sorted(found, key=lambda x: int(x[1:]))


def appendix_b_symbols(plan_text):
    """(module, name) pairs of the symbol inventory in plan appendix B.

    Appendix B is the plan's proposed API surface ("new files and function
    signatures, without touching existing ones"). The gate checks that every
    name it promises really exists in the named module, or is listed in the
    register's 附录 B 符号偏差 table (a couple of them ended up elsewhere: the
    rank tables live in util.py, and one constant became a regex).
    """
    if '### 附录 B' not in plan_text:
        return []
    sect = plan_text.split('### 附录 B', 1)[1].split('### 附录 C', 1)[0]
    blocks = re.findall(r'```\n(.*?)\n```', sect, re.S)
    out, module = [], None
    for line in (blocks[0] if blocks else '').splitlines():
        stripped = line.strip()
        if not stripped.startswith(('def ', 'class ')) \
                and re.match(r'^(neoHGT/\w+\.py|工作区)', stripped):
            module = stripped.split('(')[0].split()[0] if '.py' in stripped \
                else None
            continue
        m = re.match(r'(?:def|class)\s+(\w+)', stripped)
        if not m:
            m = re.match(r'([A-Z][A-Z0-9_]*)\s*[:=]', stripped)
        if m and module:
            out.append((module, m.group(1)))
    return out


def appendix_b_problems(symbols, code_dir, ledger):
    """Missing appendix-B symbols, minus the ones the register declares."""
    pkg = join(code_dir, 'neoHGT')
    sources = {}
    for name in sorted(os.listdir(pkg)) if isdir(pkg) else []:
        if not (name.endswith('.py') and isfile(join(pkg, name))):
            continue
        with open(join(pkg, name), encoding='utf-8') as f:
            sources[join(code_dir, 'neoHGT', name)] = f.read()
    declared = set()
    path = ledger if os.path.isabs(ledger) else join(code_dir, ledger)
    if isfile(path):
        with open(path, encoding='utf-8') as f:
            text = f.read()
        sect = doc_section(text, '附录 B 符号偏差')
        for line in sect.splitlines():
            first = row_cells(line)[0] if line.strip().startswith('|') else ''
            m = re.match(r'`(\w+)', first)
            if m:
                declared.add(m.group(1))
    problems = []
    for module, name in symbols:
        path = join(code_dir, module)
        if path not in sources:
            try:
                with open(path, encoding='utf-8') as f:
                    sources[path] = f.read()
            except OSError:
                sources[path] = ''
        present = re.search(rf'^(?:def|class)\s+{name}\b|^{name}\s*[=:]',
                            sources[path], re.M)
        if present:
            continue
        # a symbol the plan put in gtdb.py may simply live in util.py instead
        elsewhere = [p for p in sources
                     if p != path and re.search(
                         rf'^(?:def|class)\s+{name}\b|^{name}\s*[=:]',
                         sources[p], re.M)]
        if elsewhere:
            continue
        if name in declared:
            continue
        problems.append(f'附录 B 承诺的 {name}（{module}）既不存在也未登记偏差')
    return problems


def appendix_flags(plan_text):
    """`--flags` used by the example commands in appendix C and D."""
    if '### 附录 C' not in plan_text:
        return []
    sect = plan_text.split('### 附录 C', 1)[1]
    sect = sect.split('### 附录 E', 1)[0] if '### 附录 E' in sect else sect
    return sorted(set(re.findall(r'(?<![\w-])(--[a-z][a-z0-9-]+)', sect)))


def deviation_flags(ledger, code_dir):
    """`--flags` the deliverables document registers as deliberately renamed."""
    path = join(code_dir, ledger) if not os.path.isabs(ledger) else ledger
    if not isfile(path):
        return set()
    with open(path, encoding='utf-8') as f:
        text = f.read()
    sect = doc_section(text, '命令形态偏差')
    # only the plan-side column names a plan flag; the cell may show it with its
    # argument (`--to ncbi`), so keep the flag part only
    out = set()
    for line in sect.splitlines():
        line = line.strip()
        if not line.startswith('|'):
            continue
        cells = row_cells(line)
        if not cells or set(cells[0]) <= set('-: ') or cells[0] == '方案写法':
            continue
        m = re.search(r'`(--[a-z][a-z0-9-]+)', cells[0])
        if m:
            out.add(m.group(1))
    return out


def plan_task_ids(plan_text):
    """Numbered task ids of the plan's §6 phase tables."""
    if '## 6. 分阶段实施计划' not in plan_text:
        return []
    sect = plan_text.split('## 6. 分阶段实施计划', 1)[1].split('## 7.', 1)[0]
    found = set(TASK_ID_RE.findall(sect))

    def order(x):
        head, _, tail = x.partition('.')
        # the P0 spikes lead, then the phase tables in phase order
        return (0 if head.startswith('S') else 1, head,
                int(tail) if tail.isdigit() else 0)

    return sorted(found, key=order)


def plan_defect_ids(plan_text):
    """F-ids of the plan's appendix F (pre-existing defects register)."""
    if '### 附录 F' not in plan_text:
        return []
    sect = plan_text.split('### 附录 F', 1)[1]
    sect = re.split(r'\n## |\n### 附录 [A-E]', sect)[0]
    found = set(DEFECT_ID_RE.findall(sect))
    return sorted(found, key=lambda x: int(x[1:]))


LEDGER_HEADERS = ('编号', 'ID')


def ledger_rows(path):
    """{id: row} from the markdown *ledger* tables in a document.

    A ledger table is one whose header names its first column 编号 or ID; the rows
    under it are read by position (status / where / evidence). Round 54 hit the
    difference: these documents also carry prose tables -- command-shape
    deviations, the fixture coverage census, and now the dynamic-config-section
    declarations -- and the old rule was only a column count (`len(cells) < 4`),
    which happened to keep three-column tables out. The first four-column prose
    table became seven phantom "deliverables that are not plan ids", so the check
    that is meant to catch a drifting ledger was deciding what a ledger is by how
    wide its rows happened to be. Prose between two fragments of the same ledger
    does not end it, because the 用例表 continues after a heading.
    """
    with open(path, encoding='utf-8') as f:
        lines = [l.strip() for l in f]
    cells = [row_cells(l) if l.startswith('|') else None for l in lines]

    def is_sep(c):
        return bool(c) and set(c[0]) <= set('-: ')

    rows, in_ledger = {}, False
    for i, c in enumerate(cells):
        if c is None or is_sep(c):
            continue
        nxt = cells[i + 1] if i + 1 < len(cells) else None
        if is_sep(nxt):
            in_ledger = c[0] in LEDGER_HEADERS
            continue
        if in_ledger and len(c) >= 4:
            rows[c[0]] = {'id': c[0], 'status': c[1], 'where': c[2],
                          'evidence': c[3]}
    return rows


def ledger_path(token, code_dir, root_dir):
    """Resolve a backticked ledger citation to (path, line) or (None, None).

    Accepted spellings: `<code-relative>` (neoHGT/, tools/, docs/, scripts/,
    ci/, .github/, README*), `<root-relative>` (缓存测试/, 临时文件/ -- also spelled
    neoHGT2GTDB-*), a bare file name that exists in the package directory, and
    any of those with a `:NNN` or `:NNN-NNN` line citation.
    """
    tok, sep, tail = token.partition(':')
    line = None
    if sep and tail[:1].isdigit():
        line = int(tail.split('-')[0])
    elif sep:
        return None, None
    if any(c in tok for c in NOT_PATH):
        return None, None
    head, slash, rest = tok.partition('/')
    if slash and head in WORKSPACE_DIRS:
        real = join(root_dir, WORKSPACE_DIRS[head])
        return (join(real, rest) if isdir(real) else join(root_dir, tok), line)
    if tok.startswith(CODE_PREFIXES):
        return join(code_dir, tok), line
    if tok.startswith(ROOT_PREFIXES):
        return join(root_dir, tok), line
    if '/' not in tok and splitext(tok)[1] in SOURCE_EXTS:
        cand = join(code_dir, 'neoHGT', tok)
        return (cand if exists(cand) else None), line
    return None, None


def _line_count(path, cache={}):
    if path not in cache:
        try:
            with open(path, 'rb') as f:
                cache[path] = sum(1 for _ in f)
        except OSError:
            cache[path] = -1
    return cache[path]


def _read_lines(path, cache={}):
    if path not in cache:
        try:
            with open(path, encoding='utf-8') as f:
                cache[path] = f.read().splitlines()
        except OSError:
            cache[path] = []
    return cache[path]


def _line_at(path, n):
    """Text of line n of path ('' when out of range)."""
    lines = _read_lines(path)
    return lines[n - 1] if 0 < n <= len(lines) else ''


def ledger_problems(rows, ids, code_dir, root_dir, statuses=LEDGER_STATUSES,
                    tolerate_missing_workspace=False, ws_skips=None):
    """(unlisted plan ids, ids with no plan counterpart, per-row problems).

    tolerate_missing_workspace: a repository checkout (CI) has the code and
    the registers but not the workspace evidence tree, so citations of
    `缓存测试/` / `临时文件/` artifacts cannot resolve there. When the mapped
    workspace directory is absent, such citations are skipped and counted in
    ws_skips instead of failing the gate; when the workspace IS present the
    check stays strict. Default False keeps the historical strict behaviour
    for direct callers (tests)."""
    ws_absent = (tolerate_missing_workspace
                 and not isdir(join(root_dir, WORKSPACE_DIRS['缓存测试'])))
    unlisted = [i for i in ids if i not in rows]
    extra = [i for i in sorted(rows) if i not in ids]
    problems = []
    for tid, r in sorted(rows.items()):
        if r['status'] not in statuses:
            problems.append(f'{tid}: 状态 "{r["status"]}" 不在 '
                            + '/'.join(statuses))
        if r['status'] != statuses[0] and len(r['evidence']) < MIN_REASON:
            problems.append(f'{tid}: {r["status"]} 却没有写清缺的是哪一块')
        for token in re.findall(r'`([^`]+)`', f'{r["where"]} {r["evidence"]}'):
            if '::' in token:
                head, _, sym = token.partition('::')
                spath, _sline = ledger_path(head, code_dir, root_dir)
                if spath is None:
                    continue
                if not exists(spath):
                    if ws_absent and head.partition('/')[0] in WORKSPACE_DIRS:
                        if ws_skips is not None:
                            ws_skips.append(token)
                        continue
                    problems.append(f'{tid}: 引用了不存在的路径 `{token}`')
                elif not any(ln.lstrip().startswith(('def ' + sym,
                                                     'class ' + sym))
                             for ln in _read_lines(spath)):
                    problems.append(f'{tid}: `{token}` 里找不到 {sym}')
                continue
            path, line = ledger_path(token, code_dir, root_dir)
            if path is None:
                # a glob is not a citation: nobody can open it, so a wildcard
                # silently proves nothing while looking like evidence
                if any(ch in token for ch in '*?[]') and '/' in token:
                    problems.append(f'{tid}: `{token}` 用了通配符，无法核对'
                                    '（请写具体路径或 `路径::符号`）')
                continue                      # not a path-like citation
            if not exists(path):
                if ws_absent and token.partition('/')[0] in WORKSPACE_DIRS:
                    if ws_skips is not None:
                        ws_skips.append(token)
                    continue
                problems.append(f'{tid}: 引用了不存在的路径 `{token}`')
            elif line and _line_count(path) < line:
                problems.append(f'{tid}: `{token}` 的行号超出文件长度')
            elif line and '/tests/' in path.replace(os.sep, '/'):
                # test citations must land on the definition, not inside a body:
                # inserting a test shifts every following line number, and a
                # silently shifted citation is worse than no citation at all
                txt = _line_at(path, line)
                if 'def ' not in txt and 'class ' not in txt:
                    problems.append(
                        f'{tid}: `{token}` 已漂移到函数体内部（该行应为 def/class '
                        '定义行）')
    return unlisted, extra, problems


def check_ledger(label, ids, path, code_dir, root_dir, statuses,
                 tolerate_missing_workspace=False, ws_skips=None):
    """Run one ledger against its plan ids; returns a compact dict."""
    if not isfile(path):
        return {'label': label, 'path': path, 'rows': 0, 'ids': len(ids),
                'unlisted': list(ids), 'extra': [],
                'problems': [f'{label}: 缺少登记文件 {path}'],
                'counts': {}, 'ws_skipped': 0}
    rows = ledger_rows(path)
    unlisted, extra, problems = ledger_problems(
        rows, ids, code_dir, root_dir, statuses,
        tolerate_missing_workspace=tolerate_missing_workspace,
        ws_skips=ws_skips)
    problems += [f'任务 {i} 没有出现在{label}里' for i in unlisted]
    problems += [f'{label}里的 {i} 不是方案编号' for i in extra]
    counts = {s: sum(1 for r in rows.values() if r['status'] == s)
              for s in statuses}
    return {'label': label, 'path': path, 'rows': len(rows), 'ids': len(ids),
            'unlisted': unlisted, 'extra': extra, 'problems': problems,
            'counts': counts,
            'ws_skipped': len(ws_skips) if ws_skips is not None else 0}


def dead_gtdb_options(code_dir, flags):
    """GTDB-era options that are parsed but never read anywhere.

    Only flags of the GTDB surface are checked: upstream NCBI options have
    consumers in places this heuristic would not see.
    """
    pkg = join(code_dir, 'neoHGT')
    sources = {}
    for name in ('database.py', 'search.py', 'analyze.py', 'util.py',
                 'gtdb.py'):
        with open(join(pkg, name), encoding='utf-8') as f:
            sources[name] = f.read()
    dead = []
    for module, module_flags in flags.items():
        for flag in sorted(module_flags):
            dest = flag_dest(flag)
            if not (dest.startswith('gtdb') or dest in GTDB_SURFACE):
                continue
            # the declaration spells the option with hyphens, so every
            # underscore occurrence is a real read or assignment of the value
            uses = sum(src.count(dest) for src in sources.values())
            if uses == 0:
                dead.append(f'{module}:{flag}')
    return dead


# GTDB-era options that do not carry the gtdb_ prefix
GTDB_SURFACE = {'seq_source', 'keep_faa', 'keep_placeholder_names', 'mimag',
                'type_material', 'rank_sample', 'genome_list', 'max_genomes',
                'stream_extraction', 'expand_members', 'taxdump_trim',
                'allow_network', 'downloader', 'verify_md5', 'input_fingerprints',
                'hybrid_ncbi_taxdump', 'hybrid_lca_policy', 'workspace_out',
                'workspace_tmp', 'keep_unresolved', 'require_staxids',
                'gtdb_stem2taxid', 'donor_label', 'donor_ncbi_map'}


def section7_flags(plan_text):
    """The `--flag` tokens listed in the plan's §7 option tables."""
    sect, keep = None, []
    for line in plan_text.splitlines():
        if line.startswith('## '):
            sect = line
        elif sect and sect.startswith('## 7.') and line.startswith('|'):
            keep.append(line)
    return sorted(set(re.findall(r'`(--[a-z0-9-]+)`', '\n'.join(keep))))


def flags_without_tests(plan_flags, code_dir):
    """§7 flags that no test ever passes.

    'Exists in the CLI' and 'is read by the code' are both checkable from the
    source, and neither catches a flag that parses into the wrong dest or whose
    parsed value never reaches the selection. The only way to catch that is for
    something to hand the flag to the real parser, so every planned option has
    to appear in a test file -- otherwise it is reported here.
    """
    tests_dir = join(code_dir, 'neoHGT', 'tests')
    blob = ''
    if not isdir(tests_dir):
        return []
    for n in sorted(os.listdir(tests_dir)):
        if n.startswith('test_') and n.endswith('.py'):
            with open(join(tests_dir, n), encoding='utf-8',
                      errors='replace') as f:
                blob += f.read()
    return [x for x in plan_flags if x not in blob]


def cli_choice_options(code_dir):
    """Options with a `choices` list, as (stage, flags, dest, choices, default).

    Read with ast rather than imported: the gate also runs on a checkout where
    the scientific dependencies are absent, and these tables are plain data.
    All four CLI surfaces count -- `neoHGT gtdb` carries its own `--action`
    choices, and round 33 found that the option surface rule had been reading
    that stage out of the source text with a regex, which missed it entirely.
    """
    out = []
    for stage in CLI_STAGES:
        path = join(code_dir, 'neoHGT', f'{stage}.py')
        if not exists(path):
            continue
        with open(path, encoding='utf-8', errors='replace') as f:
            tree = ast.parse(f.read())
        for lst in option_tables(tree):
            for item in lst.elts:
                if (not isinstance(item, (ast.List, ast.Tuple))
                        or len(item.elts) < 3
                        or not isinstance(item.elts[0], ast.Constant)
                        or not isinstance(item.elts[2], ast.Dict)):
                    continue
                spec = {}
                for k, v in zip(item.elts[2].keys, item.elts[2].values):
                    try:
                        spec[k.value] = ast.literal_eval(v)
                    except Exception:
                        continue
                choices = spec.get('choices')
                if (not isinstance(choices, (list, tuple)) or not choices
                        or not all(isinstance(c, (str, int, float, bool))
                                   for c in choices)):
                    continue
                flags = re.findall(r'--?[a-z][\w-]*', str(item.elts[0].value))
                if not flags:
                    continue
                out.append({
                    'stage': stage, 'flags': flags,
                    'dest': str(spec.get('dest') or max(flags, key=len)
                                .lstrip('-').replace('-', '_')),
                    'choices': [str(c) for c in choices],
                    'default': str(spec.get('default'))})
    return out


def test_function_bodies(code_dir):
    """Source of every `test_*` function: the scope in which a value is handed.

    Function level rather than file level, because the suite loops over the
    values of an option (`for tier in ('high', 'medium', 'low', 'any')`) -- a
    file-wide match would let an unrelated test claim credit for a value it
    never passes, and a purely textual adjacency match would fail the loops.
    """
    tests_dir = join(code_dir, 'neoHGT', 'tests')
    out = []
    if not isdir(tests_dir):
        return out
    for n in sorted(os.listdir(tests_dir)):
        if not (n.startswith('test_') and n.endswith('.py')):
            continue
        with open(join(tests_dir, n), encoding='utf-8',
                  errors='replace') as f:
            text = f.read()
        try:
            tree = ast.parse(text)
        except SyntaxError:
            out.append(text)
            continue
        for node in ast.walk(tree):
            if (isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and node.name.startswith('test')):
                seg = ast.get_source_segment(text, node)
                if seg:
                    out.append(seg)
    return out


def value_handed_to_option(text, flags, dest, value):
    """Does this test body name *this option* and *this value*?

    Deliberately a scope test, not data flow: it catches the failure mode that
    matters -- a declared value no test thinks to pass -- while leaving the
    behavioural proof to the test itself.
    """
    named = (any(f'"{f}"' in text or f"'{f}'" in text for f in flags)
             or re.search(rf'\b{re.escape(dest)}\b', text) is not None)
    lit = re.search(r'[\'"]' + re.escape(str(value)) + r'[\'"]', text)
    return named and lit is not None


def declared_choice_values(ledger, code_dir):
    """`--flag=value` pairs the deliverables document registers as covered by
    an artifact (a benchmark or a report) instead of a unit test."""
    path = join(code_dir, ledger) if not os.path.isabs(ledger) else ledger
    if not isfile(path):
        return set()
    with open(path, encoding='utf-8') as f:
        text = f.read()
    sect = doc_section(text, '选项取值偏差')
    return {f'{m.group(1)}={m.group(2)}'
            for m in re.finditer(r'(-{1,2}[a-z][a-z0-9-]*)=([A-Za-z0-9_.-]+)',
                                 sect)}


def choice_values_without_tests(code_dir, declared=(),
                                gated_stages=CLI_STAGES):
    """Choice values no test hands to the option that declares them.

    Returns (gated, backlog). All four CLI surfaces are judged: `gated` fails
    the gate, so a value that cannot be exercised offline has to be declared in
    the deliverables document's 选项取值偏差 section instead of quietly staying
    untested. `backlog` holds any stage left out of `gated_stages`, which is
    empty as long as all of CLI_STAGES are judged.
    """
    bodies = test_function_bodies(code_dir)
    if not bodies:
        return [], []
    gated, backlog = [], []
    for ent in cli_choice_options(code_dir):
        for val in ent['choices']:
            key = f'{ent["flags"][0]}={val}'
            if key in declared:
                continue
            if any(value_handed_to_option(b, ent['flags'], ent['dest'], val)
                   for b in bodies):
                continue
            (gated if ent['stage'] in gated_stages else backlog).append(key)
    return sorted(set(gated)), sorted(set(backlog))


# the config keys this adaptation owns: everything under a `gtdb` named branch
# (config.yml is append-only for upstream keys, so `taxonomy.profiles.ncbi.*`
# stays upstream's business and is not judged here)
def _own_key(path):
    parts = path.split('.')
    return 'gtdb' in parts or 'gtdb_archaea' in parts


def _upstream_config(upstream_dir, ship):
    """The upstream config.yml to diff against, or None.

    A relative or empty `upstream_dir` is refused on purpose: with cwd inside the
    adapted tree, `join('', 'neoHGT', 'config.yml')` is the SHIPPED file, and
    diffing a file against itself yields an empty population -- a rule that
    passes by having looked at nothing. Same-file identity is checked too.
    """
    if not upstream_dir or not os.path.isabs(upstream_dir):
        return None
    up = join(upstream_dir, 'neoHGT', 'config.yml')
    if not isfile(up) or os.path.realpath(up) == os.path.realpath(ship):
        return None
    return up


def introduced_config_keys(code_dir, upstream_dir=None):
    """The config.yml leaves this adaptation added, judged by diff not by name.

    Round 40's rule asked whether a key's NAME contains "gtdb". That defined its
    population by accident: `workspace.code_dir` and `database.source` are keys we
    introduced, ship values for, and attach promises to, while their names look
    like upstream's -- so they were never candidates for "inert", and the summary
    line said "0 inert config knobs" while never having looked at them.

    The honest population is the difference between the shipped config and the
    upstream one. When no upstream copy is reachable (a CI job that only checks
    this repository out) the caller says so in the report instead of quietly
    judging a subset.
    """
    import yaml
    ship = join(code_dir, 'neoHGT', 'config.yml')
    if not isfile(ship):
        return []
    mine = set(config_leaf_paths(yaml.safe_load(open(ship, encoding='utf-8'))))
    up = _upstream_config(upstream_dir, ship)
    if up is None:
        return sorted(k for k in mine if _own_key(k))
    theirs = set(config_leaf_paths(yaml.safe_load(open(up, encoding='utf-8'))))
    return sorted(mine - theirs)


def inert_introduced_keys(code_dir, upstream_dir=None):
    """Introduced config keys that no shipped code reads.

    A profile leaf is picked by a computed name (`profiles[key]['closesize']`),
    which no static scan resolves past that level, so "same subtree, same leaf
    name" counts as a read -- but only when both halves match.
    """
    own = set(introduced_config_keys(code_dir, upstream_dir))
    reads = config_reads(code_dir, CONFIG_READ_MODULES)
    read_paths = set(reads)
    read = {k for k in own
            if k in read_paths
            or any(r.startswith('.'.join(k.split('.')[:-2]) + '.')
                   and r.split('.')[-1] == k.split('.')[-1]
                   for r in read_paths)}
    return sorted(own - read)


# The files a config read may legitimately live in. `workspace.*` is read by
# tools/env.sh, so tools count too -- a key only readable from a shell wrapper is
# still read, and excluding the wrapper would invent a defect.
CONFIG_READ_MODULES = ('neoHGT/database.py', 'neoHGT/search.py',
                       'neoHGT/analyze.py', 'neoHGT/util.py', 'neoHGT/gtdb.py',
                       'tools/make_fixtures.py', 'tools/run_hgt_demo.py',
                       'tools/verify_workspace.py')


# Every top-level section of config.yml. The rule that judges "a key nobody
# reads" used to collect reads only from three sections, which silently defined
# its population: `workspace.code_dir` was introduced by this adaptation and read
# by nothing, yet it could never be reported because `workspace` was not in the
# list of places a read is recognised.
CFG_SECTIONS = ('database', 'taxonomy', 'grouping', 'workspace', 'local',
                'program', 'search', 'analyze', 'predict', 'donor', 'diamond',
                'blast', 'remote', 'fetch', 'download', 'server', 'precomp')


def _any_config_branch(path):
    return path.split('.')[0] in CFG_SECTIONS


def _lost_a_hop(path, declared):
    """True when the scan knows the leaf and the grandparent but not the level
    between them, because that level is a computed name.

    `analyze.py` picks a profile with `prof = prof.get(key)` where key is
    'gtdb' or 'gtdb_archaea'; the resolver cannot carry that hop, so it reports
    `grouping.profiles.closesize` -- a path that is not a claim the code makes,
    just the same read with one segment missing. Calling it an invented key would
    be a false defect, and a rule that cries wolf is read as noise.
    """
    parts = path.split('.')
    if len(parts) < 2:
        return False
    parent, leaf = '.'.join(parts[:-1]), parts[-1]
    return any(d.startswith(parent + '.') and d.split('.')[-1] == leaf
               for d in declared)


def _cfg_resolve(expr, aliases, consts):
    """Config path read by an expression, or None; '' root is ().

    Aliases come from assignments inside the same function, and `consts` maps a
    loop variable to the string literals it iterates over -- analyze.py picks
    its profile keys with `for attr, name in (('close_size', 'closesize'), ...)`,
    which a scanner that only understands literal subscripts would call unread.
    """
    if expr is None:
        return None
    if isinstance(expr, ast.Name):
        if expr.id in aliases:
            return aliases[expr.id]
        return None
    if isinstance(expr, ast.Attribute):
        text = ast.unparse(expr)
        if re.search(r'\bcfg\b|\bconfig\b', text):
            return ()
        return _cfg_resolve(expr.value, aliases, consts)
    if isinstance(expr, ast.Call) and isinstance(expr.func, ast.Attribute):
        if expr.func.attr != 'get' or not expr.args:
            return None
        base = _cfg_resolve(expr.func.value, aliases, consts)
        if base is None:
            return None
        first = expr.args[0]
        if isinstance(first, ast.Constant) and isinstance(first.value, str):
            return base + (first.value,)
        if isinstance(first, ast.Name) and first.id in consts:
            return [base + (c,) for c in consts[first.id]]
        return None
    # A function whose own name says "config" returns the config: gtdb.py reads
    # its workspace defaults through `load_config_file()` rather than `self.cfg`,
    # and a scanner that only knows attribute roots would call that key unread.
    if isinstance(expr, ast.Call):
        fname = (getattr(expr.func, 'id', None)
                 or getattr(expr.func, 'attr', None))
        if fname and 'config' in fname:
            return ()
    if isinstance(expr, ast.Subscript):
        base = _cfg_resolve(expr.value, aliases, consts)
        if base is None:
            return None
        sl = expr.slice
        if isinstance(sl, ast.Constant) and isinstance(sl.value, str):
            return base + (sl.value,)
        if isinstance(sl, ast.Name) and sl.id in consts:
            return [base + (c,) for c in consts[sl.id]]
        return None
    for child in ast.iter_child_nodes(expr):
        got = _cfg_resolve(child, aliases, consts)
        if got is not None:
            return got
    return None


def _normalise(got):
    if got is None:
        return []
    return [got] if isinstance(got, tuple) else list(got)


def _loops_in(fn):
    """[(for node, {variable: [literal strings]}, {id(node) in its body})].

    One entry per loop, never merged: `search.py`'s `__init__` reuses the name
    `key` for five different lists, and a function-wide {name: literals} map lets
    every list reach every f-string -- which is how round 54's first census
    invented `grouping.name` and `diamond.threads` out of two loops that never
    meet. Keeping the loops apart and matching each call to the loop that
    *encloses* it is the scoping Python actually applies.
    """
    out = []
    for node in ast.walk(fn):
        if not isinstance(node, ast.For):
            continue
        targets = ([node.target] if isinstance(node.target, ast.Name)
                   else getattr(node.target, 'elts', []))
        names = [t.id for t in targets if isinstance(t, ast.Name)]
        elems = (node.iter.elts if isinstance(node.iter, (ast.Tuple, ast.List))
                 else [])
        lits = [e.value for e in elems
                if isinstance(e, ast.Constant) and isinstance(e.value, str)]
        d = {}
        if lits and len(names) == 1:
            d[names[0]] = lits
        else:
            # zip by position: 'attr' only ever takes the first element of each
            # pair, so its candidates are those
            for i, name in enumerate(names):
                vals = [e.elts[i].value for e in elems
                        if isinstance(e, (ast.Tuple, ast.List))
                        and len(getattr(e, 'elts', [])) > i
                        and isinstance(e.elts[i], ast.Constant)
                        and isinstance(e.elts[i].value, str)]
                if vals:
                    d[name] = vals
        if d:
            inside = {id(n) for stmt in node.body for n in ast.walk(stmt)}
            out.append((node, d, inside))
    return out


def _keys_at(call, loops):
    """Loop variables visible to this call, innermost loop winning."""
    consts = {}
    for node, d, inside in sorted(loops, key=lambda x: -x[0].lineno):
        if id(call) in inside:
            for k, v in d.items():
                consts.setdefault(k, list(v))
    return consts


def _assigned_values(tree, attr, before):
    """String constants assigned to `self.<attr>` above a line (e.g. method='remote')."""
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign) or node.lineno >= before:
            continue
        for t in node.targets:
            if (isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name)
                    and t.value.id == 'self' and t.attr == attr
                    and isinstance(node.value, ast.Constant)
                    and isinstance(node.value.value, str)):
                out.append(node.value.value)
    return out


def _entry_paths(arg, consts, value_of):
    """Concrete dotted paths an entry argument can form, or None if unresolvable.

    `value_of(attr)` answers what a `self.<attr>` used as a config *section* can
    hold; the caller decides whether that is knowable and labels the answer. A
    path is only returned when every part of the entry resolves -- an entry the
    scanner cannot expand is reported, not guessed at, because a guessed read is
    how round 51's census produced six different numbers from one tree.
    """
    if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
        return [arg.value]
    if not isinstance(arg, ast.JoinedStr):
        return None
    groups = []
    for v in arg.values:
        if isinstance(v, ast.Constant) and isinstance(v.value, str):
            groups.append([v.value])
            continue
        if not isinstance(v, ast.FormattedValue):
            return None
        e = v.value
        if isinstance(e, ast.Name) and e.id in consts:
            groups.append(consts[e.id])
            continue
        if (isinstance(e, ast.Call) and isinstance(e.func, ast.Attribute)
                and e.func.attr == 'replace'
                and isinstance(e.func.value, ast.Name)
                and e.func.value.id in consts
                and len(e.args) == 2
                and all(isinstance(x, ast.Constant) for x in e.args)):
            old, new = e.args[0].value, e.args[1].value
            groups.append([x.replace(old, new) for x in consts[e.func.value.id]])
            continue
        if isinstance(e, ast.Attribute) and isinstance(e.value, ast.Name) \
                and e.value.id == 'self':
            vals = value_of(e.attr)
            if not vals:
                return None
            groups.append(vals)
            continue
        return None
    out = ['']
    for g in groups:
        out = [base + piece for base in out for piece in g]
    return out


def get_config_entries(code_dir, modules, stages=None):
    """{dotted path: {file:line}} for every entry string a get_config() call forms.

    The other half of how config.yml is read: `config_reads()` follows dict
    subscripts (`cfg['database']['gtdb']['dir']`), while the upstream style is
    `get_config(self, key, f'search.{key}')`. Round 54 found the difference
    mattered: `for key in 'queries','maxchars','extrargs': get_config(...,
    f'{self.method}.{key}')` (search.py:306) has a *section* that comes from an
    option value, so nothing was recorded for it and `precomp.queries/maxchars/
    extrargs` -- entries the code can genuinely form -- stayed invisible to the
    "reads a key config.yml does not have" check, the same check that caught
    `database.hybrid.ncbi_taxdump` last round.
    """
    out, unexpanded = {}, set()
    options = {}
    for o in cli_choice_options(code_dir):
        options.setdefault(o['dest'], set()).update(o['choices'])
    for rel in modules:
        path = join(code_dir, rel)
        if not isfile(path):
            continue
        with open(path, encoding='utf-8', errors='replace') as fh:
            tree = ast.parse(fh.read())
        for fn in [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)]:
            loops = _loops_in(fn)
            for node in ast.walk(fn):
                if not (isinstance(node, ast.Call)
                        and getattr(node.func, 'id', '') == 'get_config'
                        and len(node.args) >= 3):
                    continue

                def value_of(attr, _line=node.lineno):
                    # an option's declared choices, plus every value the same
                    # module assigns to it before this line
                    return sorted(set(options.get(attr, set()))
                                  | set(_assigned_values(tree, attr, _line)))
                paths = _entry_paths(node.args[2], _keys_at(node, loops),
                                     value_of)
                if paths is None:
                    unexpanded.add(f'{rel}:{node.lineno}')
                    continue
                for p in paths:
                    if '.' in p:
                        out.setdefault(p, set()).add(f'{rel}:{node.lineno}')
    return out, sorted(unexpanded)


DYNAMIC_LEDGER_HEADING = '动态配置段的可达性'
DYNAMIC_KINDS = ('拒绝', '归一')


def _ledger_path(code_dir, ledger):
    """Absolute path of the deliverables document, or None when there is none.

    The rules that read prose out of the documents are also run by tests that hand
    in `ledger=None` to judge only the code; a None must mean "no declarations",
    not "join(code_dir, None)" -- and it must not quietly look like a clean pass,
    which is why `dynamic_section_problems` reports an unreadable ledger as a
    problem rather than as an empty table.
    """
    if not ledger:
        return None
    return ledger if os.path.isabs(ledger) else join(code_dir, ledger)


def _self_assigns(stmts, attr):
    """[(line, value)] for every `self.<attr> = <str>` inside these statements."""
    out = []
    for st in stmts:
        for node in ast.walk(st):
            if not isinstance(node, ast.Assign):
                continue
            if any(isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name)
                   and t.value.id == 'self' and t.attr == attr
                   for t in node.targets) \
                    and isinstance(node.value, ast.Constant) \
                    and isinstance(node.value.value, str):
                out.append((node.lineno, node.value.value))
    return out


def _definitely_rebinds(stmts, attr, value):
    """True when `self.<attr>` is certain to hold some *other* constant.

    Only two shapes are accepted, deliberately: an unconditional assignment of a
    different constant, and an `if/else` whose every branch certainly rebinds.
    A loop, a `try`, or a bare `if` without `else` can skip the assignment, so no
    witness is granted from them. Anything in the subtree that can still put the
    value back -- including a branch that assigns it later -- refuses the claim,
    because the question is what reaches the read, not what could once have been.
    """
    if any(v == value for _, v in _self_assigns(stmts, attr)):
        return False
    for st in stmts:
        if isinstance(st, ast.Assign) \
                and any(isinstance(t, ast.Attribute)
                        and isinstance(t.value, ast.Name)
                        and t.value.id == 'self' and t.attr == attr
                        for t in st.targets) \
                and isinstance(st.value, ast.Constant) \
                and st.value.value != value:
            return True
        if isinstance(st, ast.If) and st.orelse \
                and _definitely_rebinds(st.body, attr, value) \
                and _definitely_rebinds(st.orelse, attr, value):
            return True
    return False


def _witness_problems(code_dir, kind, value, attr, witness, flags, site):
    """Why a dynamic-section witness does not witness, or [] when it does."""
    rel, sep, tail = witness.partition(':')
    path = join(code_dir, rel)
    if not (sep and tail.isdigit() and isfile(path)):
        return [f'动态配置段 `{value}` 的证据 `{witness}` 不是可打开的 文件:行']
    with open(path, encoding='utf-8', errors='replace') as fh:
        src = fh.read()
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return [f'动态配置段 `{value}` 的证据 `{witness}` 所在文件解析不了']
    n = int(tail)
    if kind == '拒绝':
        hit = [x for x in ast.walk(tree) if isinstance(x, ast.Raise)
               and x.lineno == n]
        if not hit:
            return [f'`{witness}` 不是一行 `raise`：登记说这段的代码会拒绝，'
                    '但拒绝它的语句已经不在了（或挪了行）']
        text = ast.get_source_segment(src, hit[0]) or ''
        miss = [f for f in flags if f not in text]
        bad = []
        if value not in text:
            bad.append(f'`{witness}` 的拒绝理由里没有出现段名 `{value}`')
        if miss:
            bad.append(f'`{witness}` 的拒绝理由没有点名 {", ".join(miss)} '
                       '-- 用户看不出是哪个选项被拒')
        return bad
    tests = [x for x in ast.walk(tree) if isinstance(x, ast.If)
             and x.lineno == n]
    if not tests:
        return [f'`{witness}` 不是一行 `if`：登记说 `{value}` 在读取点之前已被归一，'
                '但归一它的分支已经不在了（或挪了行）']
    node = tests[0]
    st = node.test
    if not (isinstance(st, ast.Compare) and isinstance(st.left, ast.Attribute)
            and isinstance(st.left.value, ast.Name)
            and st.left.value.id == 'self' and st.left.attr == attr
            and len(st.ops) == 1 and isinstance(st.ops[0], ast.Eq)
            and isinstance(st.comparators[0], ast.Constant)
            and st.comparators[0].value == value):
        return [f'`{witness}` 判的不是 `{attr} == {value!r}`，归一声明与代码不符']
    site_line = int(site.rsplit(':', 1)[1]) if (
        ':' in site and site.rsplit(':', 1)[1].isdigit()) else None
    if site_line is not None and rel == site.split(':')[0]:
        if node.lineno >= site_line:
            return [f'`{witness}` 的归一分支在读取点 {site} 之后，管不到它']
        back = [x.lineno for x in ast.walk(tree) if isinstance(x, ast.Assign)
                and node.lineno < x.lineno < site_line
                and any(isinstance(t, ast.Attribute)
                        and isinstance(t.value, ast.Name)
                        and t.value.id == 'self' and t.attr == attr
                        for t in x.targets)
                and isinstance(x.value, ast.Constant)
                and x.value.value == value]
        if back:
            return [f'{rel}:{back[0]} 又把 `{attr}` 赋回了 `{value}`，'
                    f'`{witness}` 的归一到不了读取点 {site}']
        # only the `then` branch needs to rebind: when the test is false the
        # attribute already holds something else at that line, and the `back`
        # scan above says nothing puts the value back before the read.
        if not _definitely_rebinds(node.body, attr, value):
            return [f'`{witness}` 成立时 `{attr}` 仍可能是 `{value}`：分支体没有'
                    '在所有路径上把它改成别的取值（需要无条件赋值，或有 else 的 if）']
    return []


def dynamic_section_sites(code_dir, modules):
    """[{site, attr, values, paths}] for get_config() calls that pick a section at runtime.

    `f'{self.method}.{key}'` (search.py:306) asks config.yml for a subtree whose
    name is an *option value*, so which keys the code reads depends on what the
    user passed. The scanner can enumerate the candidates -- they are the option's
    choices plus what the module assigns -- but it cannot know which of them
    survive to the read, and a rule that pretended to would be grading its own
    guess. So the candidates are handed to rule 23, which requires a written
    disposition for every one that names no section in config.yml.
    """
    out = []
    options = {}
    for o in cli_choice_options(code_dir):
        options.setdefault(o['dest'], set()).update(o['choices'])
    for rel in modules:
        path = join(code_dir, rel)
        if not isfile(path):
            continue
        with open(path, encoding='utf-8', errors='replace') as fh:
            tree = ast.parse(fh.read())
        for fn in [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)]:
            loops = _loops_in(fn)
            for node in ast.walk(fn):
                arg = node.args[2] if (
                    isinstance(node, ast.Call)
                    and getattr(node.func, 'id', '') == 'get_config'
                    and len(node.args) >= 3) else None
                if not (isinstance(arg, ast.JoinedStr) and arg.values
                        and isinstance(arg.values[0], ast.FormattedValue)):
                    continue
                e = arg.values[0].value
                if not (isinstance(e, ast.Attribute)
                        and isinstance(e.value, ast.Name)
                        and e.value.id == 'self'):
                    continue
                attr = e.attr
                vals = sorted(set(options.get(attr, set()))
                              | set(_assigned_values(tree, attr, node.lineno)))

                def value_of(a, _line=node.lineno):
                    return sorted(set(options.get(a, set()))
                                  | set(_assigned_values(tree, a, _line)))
                paths = _entry_paths(arg, _keys_at(node, loops), value_of)
                if not vals or paths is None:
                    continue
                out.append({'site': f'{rel}:{node.lineno}', 'attr': attr,
                            'values': vals, 'paths': sorted(paths)})
    return out


def dynamic_section_declarations(ledger):
    """{section value: row} from the deliverables document's dynamic-section table."""
    if not ledger or not isfile(ledger):
        return {}
    text = _read_text(ledger)
    if DYNAMIC_LEDGER_HEADING not in text:
        return {}
    sect = doc_section(text, DYNAMIC_LEDGER_HEADING)
    rows = {}
    for line in sect.split('\n'):
        line = line.strip()
        if not line.startswith('|'):
            continue
        cells = row_cells(line)
        if len(cells) < 5 or cells[0] in ('段名',) or set(cells[0]) <= set('-: '):
            continue
        value = cells[0].strip('`')
        rows[value] = {'value': value, 'site': cells[1].strip('`'),
                       'kind': cells[2], 'witness': cells[3].strip('`'),
                       'flags': re.findall(r'--?[a-z][\w-]*', cells[4]),
                       'why': cells[5] if len(cells) > 5 else ''}
    return rows


def dynamic_section_problems(code_dir, ledger):
    """Rule 23: every config section the code chooses at runtime must have a home.

    A candidate section value that config.yml does not have means one of two
    things, and they look identical to a user: the read can never reach the file
    because the value is resolved away first (`auto`), or the read lands on a
    subtree nobody can write because the options it would configure do nothing on
    that route (`precomp`). Both are claims about behaviour, so both have to be
    written down in the deliverables document -- and the gate checks the written
    claim against the code that is supposed to make it true: a `raise` that still
    names every option it refuses, or an `if` that still rebinds the attribute on
    every path before the read.
    """
    import yaml
    try:
        sections = set(yaml.safe_load(
            open(join(code_dir, 'neoHGT', 'config.yml'))) or {})
    except (OSError, ValueError):
        sections = set()
    sites = dynamic_section_sites(code_dir, CONFIG_READ_MODULES)
    candidates = {}
    for s in sites:
        for p in s['paths']:
            head, _, leaf = p.partition('.')
            if head in sections or '.' not in p:
                continue
            candidates.setdefault(head, {}).setdefault(
                s['site'], {'leaves': set(), 'attr': s['attr'],
                            'values': s['values']})['leaves'].add(leaf)
    decls = dynamic_section_declarations(_ledger_path(code_dir, ledger))
    problems = []
    if candidates and not isfile(_ledger_path(code_dir, ledger) or ''):
        problems.append(
            f'{len(candidates)} 个按取值选择的配置段无人处置，而登记表 '
            f'{ledger!r} 读不到 —— 缺文件不等于没问题')
    for value in sorted(set(candidates) - set(decls)):
        where = ', '.join(sorted(candidates[value]))
        leaves = sorted({l for m in candidates[value].values()
                         for l in m['leaves']})
        problems.append(
            f'代码会按取值选择配置段 `{value}`（{where} 读 {", ".join(leaves)}），'
            f'但 config.yml 没有 `{value}` 段，登记表里也没有这一行 —— '
            '要么给它一个段，要么写清楚这些读取为什么永远落不了地')
    for value in sorted(set(decls) - set(candidates)):
        problems.append(
            f'登记表申明了动态配置段 `{value}`，但代码里已没有任何 get_config() '
            '调用能形成它 —— 段改名或读取点被删之后，这一行不再核对任何东西')
    for value, d in sorted(decls.items()):
        if value not in candidates:
            continue
        if d['kind'] not in DYNAMIC_KINDS:
            problems.append(f'动态配置段 `{value}` 的处置 `{d["kind"]}` 不是 '
                            f"{' 或 '.join(DYNAMIC_KINDS)} 之一")
            continue
        if d['site'] not in candidates[value]:
            problems.append(
                f'动态配置段 `{value}` 登记的读取点 `{d["site"]}` 不是代码实际形成'
                f'它的那一处（可选：{", ".join(sorted(candidates[value]))}）')
            continue
        info = candidates[value].get(d['site'], {})
        if d['kind'] == '归一' and d['flags']:
            problems.append(f'动态配置段 `{value}` 归一后即不可达，登记里却点名了 '
                            f'{", ".join(d["flags"])} —— 拒绝的写法用错了处置')
        if d['kind'] == '拒绝':
            want = {f'--{l.replace("_", "-")}' for l in info.get('leaves', ())}
            missing = sorted(want - set(d['flags']))
            if missing:
                problems.append(
                    f'动态配置段 `{value}` 在 {d["site"]} 形成 '
                    f"{', '.join(sorted(want))} 这些读取，登记的覆盖选项少了 "
                    f'{", ".join(missing)} —— 少登记一个就等于那个旋钮没人管')
        problems += _witness_problems(code_dir, d['kind'], value,
                                      info.get('attr', ''), d['witness'],
                                      d['flags'], d['site'])
    return {'candidates': {v: sorted(candidates[v]) for v in candidates},
            'candidate_leaves': {
                v: sorted({l for m in candidates[v].values()
                           for l in m['leaves']}) for v in candidates},
            'sections_present': sorted(set(candidates) & set(decls)),
            'declarations': {v: decls[v]['kind'] for v in sorted(decls)}}, problems



# The prose sections this register slices out of our own documents. A heading
# that disappears does not fail the rule that reads it -- the rule reads an empty
# string and reports "nothing to check". Round 53 lost a `###` heading to a bad
# edit and round 54 lost the `## 第五十二轮新增交付物` heading the same way, twice
# in a row, so "the document still has the section the gate reads" is now itself
# a checked claim instead of something to remember.
DECLARED_SECTIONS = (
    ('附录 B 符号偏差', 'appendix_b_problems'),
    ('命令形态偏差', 'deviation_flags'),
    ('选项取值偏差', 'declared_choice_values'),
    ('无人自动执行的核对工具', 'checker_tool_problems'),
    ('断言缺失的既有测试', 'declared_vacuous_tests'),
    ('可达性偏差', 'verify_production_reachability.py'),
    (DYNAMIC_LEDGER_HEADING, 'dynamic_section_problems'),
)


def ledger_heading_problems(code_dir, ledger):
    """Rule 24: every section the register reads must still be a heading."""
    path = _ledger_path(code_dir, ledger)
    if not path or not isfile(path):
        return {'sections': 0}, [f'登记门禁要读 {path!r}，但文件不在 —— '
                                 '读不到不等于没问题']
    text = _read_text(path)
    titled = {l.lstrip('# ').strip() for l in text.split('\n')
              if l.startswith('## ') or l.startswith('### ')}
    problems, seen = [], []
    for token, owner in DECLARED_SECTIONS:
        hit = [t for t in titled if token in t]
        if not hit:
            problems.append(
                f'`{owner}` 要从文档里读「{token}」一节，但该标题已不存在 —— '
                '这一节的核对会读到空串并交出"没问题"，这是最坏的一种绿')
        else:
            seen.append(token)
    return {'sections': len(seen), 'read': sorted(seen)}, problems


def config_reads(code_dir, modules):
    """{(subtree, key): {file}} for every config path the code touches."""
    out = {}
    for rel in modules:
        path = join(code_dir, rel)
        if not isfile(path):
            continue
        with open(path, encoding='utf-8', errors='replace') as f:
            tree = ast.parse(f.read())
        for fn in [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)]:
            loops = _loops_in(fn)
            aliases = {}
            for node in ast.walk(fn):
                if isinstance(node, ast.Assign) and len(node.targets) == 1 \
                        and isinstance(node.targets[0], ast.Name):
                    got = _normalise(_cfg_resolve(node.value, aliases,
                                                 _keys_at(node, loops)))
                    if len(got) == 1:
                        aliases[node.targets[0].id] = got[0]
            for node in ast.walk(fn):
                if not isinstance(node, (ast.Subscript, ast.Call)):
                    continue
                for got in _normalise(_cfg_resolve(node, aliases,
                                                  _keys_at(node, loops))):
                    if len(got) >= 2 and _any_config_branch('.'.join(got)):
                        out.setdefault('.'.join(got), set()).add(rel)
    return out


def config_leaf_paths(cfg):
    """{dotted.path: value} for the leaves of a parsed config dict."""
    out = {}

    def walk(node, prefix):
        if isinstance(node, dict):
            for k, v in node.items():
                walk(v, prefix + (str(k),))
        else:
            out['.'.join(prefix)] = node

    walk(cfg, ())
    return out


def config_key_problems(code_dir, ledger, upstream_dir=None):
    """Every key this adaptation introduced in config.yml must be read; every
    key the code reads under those subtrees must exist.

    The rule found `database.gtdb.keep_placeholder_names` in round 40: a knob
    shipped in config.yml, named after a command-line option, read by nobody --
    the search stage reads `taxonomy.profiles.gtdb.keep_placeholder` instead.
    A user could set it to anything and nothing would change.

    Round 52 widened the population from "name contains gtdb" to the diff against
    the upstream config, which is what turned up `workspace.code_dir` (a
    docstring promised the workspace root came from it) and `database.source`
    (config.yml's own comment promised it selects the source).

    Round 54 added the other half of how config.yml is read: `config_reads()` only
    follows dict subscripts, while the upstream style is
    `get_config(self, key, f'search.{key}')` -- measured as two disjoint sets (26
    subscript paths, 73 entry paths, none in common), so the entry-string reads
    had never been judged at all. Paths whose *section* is an option value are not
    judged here; rule 23 owns them, because "can the code reach this section" is a
    behaviour question a static scan cannot answer.
    """
    import yaml
    declared = config_leaf_paths(
        yaml.safe_load(open(join(code_dir, 'neoHGT', 'config.yml'))))
    own = set(introduced_config_keys(code_dir, upstream_dir))
    population = ('upstream diff' if _upstream_config(
        upstream_dir, join(code_dir, 'neoHGT', 'config.yml'))
        else 'name heuristic (no upstream config here)')
    sub = config_reads(code_dir, CONFIG_READ_MODULES)
    subscripts = len(sub)
    entries, unexpanded = get_config_entries(code_dir, CONFIG_READ_MODULES)
    reads = dict(sub)
    dyn = dynamic_section_declarations(_ledger_path(code_dir, ledger))
    for p, where in entries.items():
        reads.setdefault(p, set()).update(where)
    # a section rule 23 has already dispositioned is not this rule's finding
    owned = {p for p in entries if p.split('.')[0] in dyn}
    # a profile is selected by a computed name (`profiles[key]['closesize']`),
    # which no static scan can resolve past that level; treat "same subtree,
    # same leaf name" as a read, but require both halves to match
    read_paths = set(reads)
    read_leaves = {k for k in own
                   if k in read_paths
                   or any(r.startswith('.'.join(k.split('.')[:-2]) + '.')
                          and r.split('.')[-1] == k.split('.')[-1]
                          for r in read_paths)}
    prefixes = {p for k in declared for p in
                ['.'.join(k.split('.')[:i]) for i in range(1, len(k.split('.')))]}
    invented = {k: sorted(v) for k, v in reads.items()
                if k not in declared and k not in prefixes and k not in owned
                and _any_config_branch(k) and not _lost_a_hop(k, declared)}
    never = sorted(own - read_leaves)
    problems = [f'config.yml 声明了 {k}，但没有任何代码读它 —— 一个改了没反应的旋钮'
                for k in never]
    problems += [f'{", ".join(where)} 读了 config.yml 里没有的键 {k}'
                 for k, where in sorted(invented.items())]
    problems += [f'{site} 的 get_config() 调用展开不出具体配置键 —— '
                 '扫描器看不见这次读取，既不能说它合法也不能说它越界，改掉写法或在'
                 '登记表里说明' for site in unexpanded]
    return {'declared_own_keys': sorted(own), 'read_keys': sorted(read_leaves),
            'never_read': never, 'read_but_absent': invented,
            'entry_paths': len(entries), 'entry_unexpanded': unexpanded,
            'subscript_paths': len(reads) - len(entries),
            'population': population, 'problems': problems}, problems


def checker_tools(code_dir):
    """{tool: who runs it}. A checker nobody runs is documentation that rots.

    README and the registers cite `python tools/verify_*.py ...` as the command
    that produced a committed piece of evidence. If nothing executes that file
    any more, the evidence can no longer be reproduced, and the tool can break
    without any signal -- round 28 found the same shape inside the package.
    """
    tools_dir = join(code_dir, 'tools')
    if not isdir(tools_dir):
        return {}
    runners = []
    for rel in ('.github/workflows/ci.yml', '.github/workflows/nightly.yml'):
        path = join(code_dir, rel)
        runners.append((rel, _read_text(path)))
    for folder in ('neoHGT/tests', 'tools'):
        base = join(code_dir, folder)
        if not isdir(base):
            continue
        for name in sorted(os.listdir(base)):
            if name.endswith('.py') and not name.startswith('__'):
                runners.append((join(folder, name),
                                _read_text(join(base, name))))
    out = {}
    for name in sorted(os.listdir(tools_dir)):
        # shell tools carry the same risk this rule exists for: run.sh, env.sh,
        # round_witness.sh and gate_round_witness.sh enforce the workspace rules,
        # and nothing notices when they stop being executed or start failing
        # quietly -- a checker that only covers Python checkers has a blind spot
        # shaped exactly like the tools that guard this workspace
        if not ((name.startswith('verify_') and name.endswith('.py'))
                or name.endswith('.sh')):
            continue
        who = [f'{src}' for src, text in runners
               if name in text and not src.endswith(name)]
        out[name] = who
    return out


def checker_tool_problems(tools, ledger):
    """Every verify_*.py must be run by something, or be declared."""
    declared = set()
    if isfile(ledger):
        with open(ledger, encoding='utf-8') as f:
            text = f.read()
        sect = doc_section(text, '无人自动执行的核对工具')
        if sect:
            declared = set(re.findall(
                r'`((?:verify_\w+\.py)|(?:[A-Za-z_][\w-]*\.sh))`', sect))
    problems = []
    for name, who in sorted(tools.items()):
        if not who and name not in declared:
            problems.append(f'{name} 没有任何 CI 步骤、用例或工具引用它 —— '
                            '文档里那条命令再也无人执行（要跑它就跑，否则在登记表里申明）')
    return problems, sorted(declared)


ASSERTION_MARK_RE = re.compile(
    r'self\.assert|\bassert\s|self\.fail|pytest\.raises|assertRaises|'
    r'assertRegex|assertNotRegex|skipTest|skipIf|skipUnless')


def assertion_free_tests(code_dir):
    """`file::test_name` for every test function with no observable assertion.

    A test whose body only calls something passes just as happily when the
    called code is a no-op, which is the same "declared but inert" shape the
    option and reachability rules chase -- one level lower down. Round 37 found
    two of this adaptation's own tests written that way.
    """
    tests = join(code_dir, 'neoHGT', 'tests')
    out = []
    if not isdir(tests):
        return out
    for name in sorted(os.listdir(tests)):
        if not (name.startswith('test_') and name.endswith('.py')):
            continue
        path = join(tests, name)
        with open(path, encoding='utf-8', errors='replace') as f:
            src = f.read()
        try:
            tree = ast.parse(src)
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name.startswith('test_'):
                seg = ast.get_source_segment(src, node) or ''
                if not ASSERTION_MARK_RE.search(seg):
                    out.append(f'{name}::{node.name}')
    return sorted(out)


def declared_vacuous_tests(ledger):
    """Upstream tests registered as "asserts nothing" in the deliverables doc."""
    if not isfile(ledger):
        return set()
    with open(ledger, encoding='utf-8') as f:
        text = f.read()
    if '断言缺失的既有测试' not in text:
        return set()
    sect = doc_section(text, '断言缺失的既有测试')
    return {f'{m.group(1)}::{m.group(2)}' for m in
            re.finditer(r'`([\w.]+\.py)::(test_\w+)`', sect)}


# --- reuse route (--from-scores) option reachability ------------------------
# `--from-scores` jumps from argument validation to prediction. Any option whose
# only consumer lives in a skipped step is parsed, validated, and then cannot
# take effect -- round 41 measured 13 of them while fixing the two grouping
# knobs it collided with. The register in `neoHGT/analyze.py::REUSE_INERT` is
# the user-facing half (it is what the run prints), and this rule is the half
# that keeps the register from rotting: the list is re-derived from the call
# graph on every run.
REUSE_SKIPPED_STEPS = ('read_input', 'assign_taxonomy', 'define_groups',
                       'calc_scores', 'make_score_table')
REUSE_VALIDATION = ('set_parameters', 'set_parameters_gtdb')
# Options the register declares dead on this route even though the graph finds a
# live reader for them. Every entry has to be justified, and the justification
# is itself tested (`TestReuseRouteInertOptions`), because "the graph cannot see
# this one" is the easiest place to hide a mistake:
#   input_tax -- `domain_of_input()` does read it, but only as a *dict*, and the
#                step that turns `--input-tax` into that dict is the skipped
#                `read_input()`. On this route the value arrives as a path
#                string, is not a dict, and is passed over for the donor vote.
#   self_rank -- `apply_gtdb_grouping()` reads it to rewrite 'domain' into
#                'superkingdom', and nothing else on this route consumes a self
#                rank: the group itself is built by the skipped `define_groups()`.
#                (`--donor-rank`/`--tax-unirank` are *not* here: their ranks are
#                used at report time, which this route does run.)
REUSE_ROUTE_EXCEPTIONS = {
    'input_tax': 'read as a dict by domain_of_input(), but the string -> dict '
                 'parse happens in the skipped read_input()',
    'self_rank': 'only rewritten by apply_gtdb_grouping(); the group that would '
                 'use it is built by the skipped define_groups()',
}


def reuse_route_graph(code_dir):
    """(live functions, blocked functions, {function: names it reads}).

    Deciding this needs the same two rebind idioms the reachability auditor
    needs: `Analyze.read_input = read_input_gtdb` (the CLI arrives under the old
    attribute name) and `_orig_read_input = Analyze.read_input` (the wrapper
    still calls the upstream body).

    Traversal refuses to *enter* a skipped step rather than subtracting its
    subtree afterwards. The subtraction version deleted any helper both routes
    reach -- `detect_gtdb_root()` is exactly that, called from the reuse route's
    prediction as well as from `read_input()` -- and so reported a finished fix
    as still missing: a graph rule that cannot see a fix is not a check.
    """
    path = join(code_dir, 'neoHGT', 'analyze.py')
    if not isfile(path):
        return set(), set(), {}
    with open(path, encoding='utf-8', errors='replace') as f:
        tree = ast.parse(f.read())
    fns = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            fns.setdefault(node.name, node)      # methods live in the class body
    rebound, orig = {}, {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        tgt, val = node.targets[0], node.value
        if isinstance(tgt, ast.Attribute) and isinstance(val, ast.Name):
            rebound[tgt.attr] = val.id
        if isinstance(tgt, ast.Name) and isinstance(val, ast.Attribute):
            orig[tgt.id] = val.attr
    blocked = set(REUSE_SKIPPED_STEPS)
    blocked |= {rebound.get(n, n) for n in REUSE_SKIPPED_STEPS}

    def edges(fn):
        out = set()
        for sub in ast.walk(fn):
            if not isinstance(sub, ast.Call):
                continue
            fun = sub.func
            if (isinstance(fun, ast.Attribute)
                    and isinstance(fun.value, ast.Name)
                    and fun.value.id == 'self'):
                out.add(rebound.get(fun.attr, fun.attr))
            if (sub.args and isinstance(sub.args[0], ast.Name)
                    and sub.args[0].id == 'self'):
                if isinstance(fun, ast.Name):
                    out.add(orig.get(fun.id, fun.id))
                elif isinstance(fun, ast.Attribute):
                    out.add(fun.attr)
        return {x for x in out if x in fns}

    eg = {name: edges(fn) for name, fn in fns.items()}
    seen, stack = set(), [r for r in ('__call__',) if r not in blocked]
    while stack:
        cur = stack.pop()
        for nxt in eg.get(cur, ()):
            if nxt in blocked or nxt in seen:
                continue
            seen.add(nxt)
            stack.append(nxt)
    seen.add('__call__')

    reads = {}
    for name, fn in fns.items():
        par = {}
        for sub in ast.walk(fn):
            for kid in ast.iter_child_nodes(sub):
                par[id(kid)] = sub
        # `for k in ('taxdump', 'input', 'output'): getattr(self, k)` -- the
        # indirect form. Without this the graph reads `--input` as dead, because
        # the only string naming it sits in a tuple, not in the getattr call.
        looped = {}
        for sub in ast.walk(fn):
            if isinstance(sub, ast.For) and isinstance(sub.target, ast.Name):
                src, var = sub.iter, sub.target.id
            elif (isinstance(sub, ast.comprehension)
                  and isinstance(sub.target, ast.Name)):
                src, var = sub.iter, sub.target.id
            else:
                continue
            if isinstance(src, (ast.Tuple, ast.List)):
                names = {e.value for e in src.elts
                         if isinstance(e, ast.Constant)
                         and isinstance(e.value, str)}
                if names:
                    looped.setdefault(var, set()).update(names)
        got = set()
        for sub in ast.walk(fn):
            p = par.get(id(sub))
            if (isinstance(sub, ast.Attribute)
                    and isinstance(sub.value, ast.Name)
                    and sub.value.id == 'self'):
                # a write (`self.gtdb_root = root`) is not a consumption
                if isinstance(p, ast.Assign) and any(t is sub for t in p.targets):
                    continue
                if isinstance(p, ast.AugAssign):
                    continue
                got.add(sub.attr)
            if (isinstance(sub, ast.Constant) and isinstance(sub.value, str)
                    and isinstance(p, ast.Call)
                    and any(isinstance(a, ast.Name) and a.id == 'self'
                            for a in p.args)):
                # getattr(self, 'dest') / get_config(self, 'dest', ...): the
                # string form of the same read
                got.add(sub.value)
            if (isinstance(sub, ast.Name) and sub.id in looped
                    and isinstance(p, ast.Call)
                    and p.args and isinstance(p.args[0], ast.Name)
                    and p.args[0].id == 'self'):
                got |= looped[sub.id]
        reads[name] = got
    return seen, blocked, reads


def declared_reuse_inert(path):
    """The REUSE_INERT tuple in analyze.py, read with ast (not imported)."""
    if not isfile(path):
        return []
    with open(path, encoding='utf-8', errors='replace') as f:
        tree = ast.parse(f.read())
    for node in tree.body:
        if not (isinstance(node, ast.Assign) and isinstance(node.value, ast.Tuple)
                and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id == 'REUSE_INERT'):
            continue
        return [e.value for e in node.value.elts
                if isinstance(e, ast.Constant) and isinstance(e.value, str)]
    return []


def reuse_inert_problems(code_dir):
    """Rule: the declared dead-option list == the call graph's, and it is shown.

    Three directions, because all three have happened while writing this rule:
    an exemption that exempts nothing (my first version excused `input`,
    `output`, `taxdump` and `from_scores` as "route plumbing" -- all four have
    live readers and are honoured), a declared-dead option that the graph says
    is read (`input_tax`, which needs a *reason*, not a shrug), and a dead
    option nobody declared (the user gets silence). A stale exemption is the
    worst of the three: it looks like a considered decision and hides the next
    real change.
    """
    live, blocked, reads = reuse_route_graph(code_dir)
    path = join(code_dir, 'neoHGT', 'analyze.py')
    declared = set(declared_reuse_inert(path))
    dests = cli_flags(code_dir).get('analyze', set())
    problems = []
    if not dests:
        problems.append('analyze 表面一个选项都没读到：判定无法进行（判据不能空跑）')
    if len(live) < 10:
        # the first version of this traversal walked only module-level functions,
        # never saw a single method, and so "proved" that nothing is inert --
        # a blind graph is indistinguishable from a fixed product unless the
        # rule refuses to be satisfied by an empty walk
        problems.append(f'调用图只走到 {len(live)} 个函数：多半是没看见方法，'
                        '不是真的没有跳过任何东西')
    derived_all, consumers = set(), {}
    for flag in sorted(dests):
        dest = flag_dest(flag)
        who = {n for n, got in reads.items() if dest in got}
        users = who - set(REUSE_VALIDATION) - blocked
        live_users = users & live
        consumers[dest] = sorted(live_users or users)
        if who and not live_users:
            derived_all.add(dest)
    derived = derived_all - set(REUSE_ROUTE_EXCEPTIONS)
    for dest in sorted(REUSE_ROUTE_EXCEPTIONS):
        if dest not in declared:
            problems.append(f'豁免 {dest} 没有对应的 REUSE_INERT 申明，'
                            '等于凭空放行一个选项')
        elif dest in derived_all:
            problems.append(f'豁免 {dest} 已经用不上：调用图本来就没找到'
                            '它的可达读者，留着它会掩盖下一次的真实变化')
    for dest in sorted(derived_all - set(REUSE_ROUTE_EXCEPTIONS) - declared):
        problems.append(f'调用图判定 {dest} 在复用路径无人消费，'
                        '但它既不在 REUSE_INERT 里也没有豁免理由')
    for dest in sorted(declared - derived_all - set(REUSE_ROUTE_EXCEPTIONS)):
        problems.append(f'REUSE_INERT 申明 {dest} 在复用路径不生效，'
                        f'但调用图显示它还有可达的读者 '
                        f'{consumers.get(dest) or ""}')
    text = _read_text(path)
    if 'NOTE: --from-scores reuses the score table' not in text:
        problems.append('复用路径的提示不见了：REUSE_INERT 必须被某条打印读出来')
    if 'REUSE_INERT' not in text:
        problems.append('neoHGT/analyze.py 里找不到 REUSE_INERT 的申明')
    return {'declared': sorted(declared),
            'derived': sorted(derived),
            'exceptions': sorted(REUSE_ROUTE_EXCEPTIONS),
            'live_functions': len(live),
            'blocked_functions': sorted(blocked),
            }, problems


CN_DIGITS = {'〇': 0, '零': 0, '一': 1, '二': 2, '三': 3, '四': 4, '五': 5,

             '六': 6, '七': 7, '八': 8, '九': 9}
DISCOVERY_HEADING_RE = re.compile(
    r'^###\s*第(?P<from>[〇零一二三四五六七八九十]+)至'
    r'(?P<to>[〇零一二三四五六七八九十]+)轮的自我发现：'
    r'(?P<count>[〇零一二三四五六七八九十百\d]+)个')
ROUND_TAG_RE = re.compile(r'第([〇零一二三四五六七八九十百]+)轮')
LEAD_RE = re.compile(r'\*\*(.+?)\*\*', re.S)
ITEM_START_RE = re.compile(r'^\s*(\d+)\.\s')


def cn_int(text):
    """Chinese numeral up to 999 ('三十二' -> 32); None when it is not one.

    The list's round numbers and item count are written as Chinese numerals
    because that is how the surrounding prose reads; parsing them is what lets
    the heading be checked against the list instead of trusted.
    """
    if text.isdigit():
        return int(text)
    if '十' in text:
        head, _, tail = text.partition('十')
        if head and head not in CN_DIGITS:
            return None
        if tail and tail not in CN_DIGITS:
            return None
        return CN_DIGITS.get(head, 1) * 10 + CN_DIGITS.get(tail, 0)
    if '百' in text:
        head, _, rest = text.partition('百')
        if head not in CN_DIGITS:
            return None
        return CN_DIGITS[head] * 100 + (cn_int(rest) or 0) if rest \
            else CN_DIGITS[head] * 100
    return CN_DIGITS.get(text)


def self_discovery_check(path):
    """The self-discovery list has to say how many items it holds.

    Its heading carries two hand-maintained numbers -- the round span and the
    item count -- while the list itself grows one item per round. Round 32 came
    within one edit of shipping "十五个" above a list of sixteen, which is the
    shape of mistake that makes a running log stop being citable: the reader
    cannot tell whether the missing item was deleted or never written down.
    """
    out = {'path': path, 'heading': None, 'rounds': None, 'declared_count': None,
           'items': 0, 'numbering': [], 'tagged_rounds': [], 'problems': []}
    if not isfile(path):
        out['problems'] = [f'找不到自我发现清单文件 {path}']
        return out
    with open(path, encoding='utf-8') as f:
        lines = f.read().splitlines()
    at = next((i for i, line in enumerate(lines)
               if DISCOVERY_HEADING_RE.match(line)), None)
    if at is None:
        out['problems'] = ['自我发现清单的标题不再写明轮次区间与条数'
                           '（应为「### 第X至Y轮的自我发现：N个…」）']
        return out
    m = DISCOVERY_HEADING_RE.match(lines[at])
    out['heading'] = lines[at].lstrip('#').strip()
    lo, hi, count = (cn_int(m.group(k)) for k in ('from', 'to', 'count'))
    out['rounds'] = [lo, hi]
    out['declared_count'] = count
    if None in (lo, hi, count):
        out['problems'] = [f'自我发现清单标题的编号无法解析: {out["heading"]}']
        return out
    if lo > hi:
        out['problems'] = [f'自我发现清单标题的轮次区间是反的: 第{lo}至{hi}轮']
    body, nums = [], []
    for line in lines[at + 1:]:
        if line.startswith('### '):
            break
        start = ITEM_START_RE.match(line)
        if start:
            nums.append(int(start.group(1)))
            body.append([line])
        elif body and (line.startswith('    ') or line.startswith('\t')):
            body[-1].append(line)
        elif line.strip() == '':
            continue
        else:
            if body:
                body[-1].append(line)     # a wrapped item continues flush left
    # the round is read out of each item's bold lead, not out of the whole item:
    # prose in an entry regularly names *other* rounds (the round a defect was
    # first found in), and those references must not drive the range check
    leads = []
    for item in body:
        text = '\n'.join(item)
        bold = LEAD_RE.search(text)
        found = [cn_int(t) for t in ROUND_TAG_RE.findall(bold.group(1))] \
            if bold else []
        leads.append([r for r in found if r is not None])
    out['numbering'] = nums
    out['items'] = len(nums)
    out['tagged_rounds'] = [r for lead in leads for r in lead]
    if nums != list(range(1, len(nums) + 1)):
        out['problems'].append(
            f'自我发现清单的编号不连续或重复: {nums}')
    if count != len(nums):
        out['problems'].append(
            f'自我发现清单标题写着 {count} 个，实际有 {len(nums)} 条')
    stray = [r for r in out['tagged_rounds'] if r < lo or r > hi]
    if stray:
        out['problems'].append(
            f'清单里小标题标着第{"、".join(str(r) for r in sorted(set(stray)))}轮的'
            f'条目落在标题声明的区间（第{lo}至{hi}轮）之外')
    if not leads or not leads[-1]:
        out['problems'].append(
            f'自我发现清单最后一条（第 {nums[-1] if nums else "?"} 条）的加粗小标题'
            '没有写明它出自哪一轮，标题的轮次区间因此无从核对')
    elif leads[-1][0] != hi:
        out['problems'].append(
            f'标题的轮次区间止于第{hi}轮，但清单最新一条出自第{leads[-1][0]}轮')
    return out


def invocation_flags(line, subs):
    """[(sub-command, flags)] for every `neoHGT <subcommand>` on one line.

    Only text following a known sub-command, and only up to the next one:
    `conda remove -n neoHGT --all` names the product without being our command
    line, and a rule that scanned whole lines would have to memorise every
    third-party tool on earth to stay quiet. A register row may hold two
    invocations (the plan's spelling beside ours), so each tail must be
    attributed to the invocation that starts it.
    """
    hits = [(m, m.group(1)) for m in re.finditer(r'\bneoHGT\s+([a-z][a-z-]*)',
                                                 line)
            if m.group(1) in subs]
    out = []
    for i, (m, stage) in enumerate(hits):
        end = hits[i + 1][0].start() if i + 1 < len(hits) else len(line)
        tail = line[m.end():end]
        out.append((stage, set(re.findall(r'(?<![\w-])(--[a-z][a-z0-9-]*)',
                                          tail))))
    return out


FENCE_RE = re.compile(r'^```[A-Za-z0-9_+#.\-]*$')


def invocation_values(line, subs):
    r"""[(sub-command, flag, raw value)] for one line of a *fenced* example.

    Values are only judged inside ``` blocks: prose says things like
    "`--source gtdb`）仍然**要求** metadata 表", and a rule that scanned it would
    report `gtdb\`）仍然**要求**` as a rejected value -- round 33's shape again, a
    text rule reading structure it cannot see. The first version of this probe did
    exactly that and "found" six violations in clean documents.

    `--flag value` and `--flag=value` are both read; a token that is another flag
    means the option takes no value here, and a token with trailing punctuation
    (`--close-tax\`,`) is prose mentioning a flag, not a command line.
    """
    hits = [(m, m.group(1)) for m in re.finditer(r'\bneoHGT\s+([a-z][a-z-]*)',
                                                 line)
            if m.group(1) in subs]
    out = []
    for i, (m, stage) in enumerate(hits):
        end = hits[i + 1][0].start() if i + 1 < len(hits) else len(line)
        tokens = line[m.end():end].split()
        for j, tok in enumerate(tokens):
            if not tok.startswith('--'):
                continue
            flag = re.match(r'--[a-z][a-z0-9-]*', tok).group(0) if \
                re.match(r'--[a-z][a-z0-9-]*', tok) else None
            if flag is None:
                continue
            if tok != flag:
                if not tok.startswith(flag + '='):
                    continue      # trailing punctuation: prose, not a command
                out.append((stage, flag, tok[len(flag) + 1:]))
                continue
            nxt = tokens[j + 1] if j + 1 < len(tokens) else ''
            if nxt.startswith('-') or nxt == '':
                continue
            out.append((stage, flag, nxt))
    return out


def doc_command_values(code_dir, subs, extra_file=None):
    """[(file:line, stage, flag, value)] for every example inside a code fence.

    Returns (pairs, judged, placeholders, unbalanced). `extra_file` lets the same
    reader judge the plan document itself: its command tables and examples are
    what the code is measured against, and a wrong value there is a plan bug that
    every downstream copy inherits.
    """
    files = [join(code_dir, 'README.md')]
    docs = join(code_dir, 'docs')
    if isdir(docs):
        files += [join(docs, f) for f in sorted(os.listdir(docs))
                  if f.endswith('.md')]
    if extra_file and isfile(extra_file):
        files.append(extra_file)
    pairs, judged, placeholders, unbalanced = [], [], [], []
    for path in files:
        if not isfile(path):
            continue
        with open(path, encoding='utf-8') as f:
            text = f.read().replace('\\\n', ' ')
        fence = False
        for n, line in enumerate(text.splitlines(), 1):
            # a fence boundary is a line that holds *only* the marker (with an
            # optional language): upstream writes `​```server.py`​`` inline at the
            # start of a sentence, and toggling on any line that begins with
            # backticks flips the whole file's parity -- the first version of
            # this rule therefore "found" an unbalanced README and would have
            # blindened every check after it
            if FENCE_RE.match(line.strip()):
                fence = not fence
                continue
            if not fence:
                continue
            for stage, flag, raw in invocation_values(line, subs):
                where = f'{basename(path)}:{n}'
                val = raw.strip('.,;:()"\'`')
                # `--action probe|taxdump|map-report` documents the whole set;
                # every alternative has to be a real value
                parts = [p for p in re.split(r'\||/| 或 ', val) if p]
                pairs.append((where, stage, flag, val, parts))
                if re.match(r'^[\$<>~*\[]', val) or re.match(r'^[-_.]+$', val):
                    placeholders.append((where, stage, flag, val))
                else:
                    judged.append((where, stage, flag, val, parts))
        if fence:
            # an odd number of fences means the parser believes half the file is
            # outside a code block, so those examples are never judged -- the
            # blindness has to be reported, not quietly skipped
            unbalanced.append(basename(path))
    return pairs, judged, placeholders, unbalanced


def doc_value_problems(code_dir, subs=None, plan=None):
    """Rule 18: a documented example's value must exist on *that* sub-command.

    `--source` is the case this rule was written for: the database stage accepts
    ncbi|gtdb|hybrid (default ncbi) while search and analyze accept
    auto|ncbi|gtdb|hybrid -- so `neoHGT database --source auto …` is a
    plausible-looking line that argparse rejects, and the flag-level rule is
    blind to it because the flag itself does exist everywhere.

    The plan document is judged too when it is available, and its absence is
    reported as `not checked` rather than as a pass: CI runs without the plan
    file, so a silent skip would look identical to a clean plan.
    """
    subs = subs or cli_subcommands(code_dir)
    choices = {}
    for entry in cli_choice_options(code_dir):
        for flag in entry['flags']:
            choices.setdefault(entry['stage'], {})[flag] = set(entry['choices'])
    pairs, judged, placeholders, unbalanced = doc_command_values(
        code_dir, subs, plan)
    problems = []
    plan_label = basename(plan) if plan else ''
    for where, stage, flag, val, parts in judged:
        opts = choices.get(stage, {})
        if flag not in opts:
            continue          # free-form value, or the flag rule's business
        for part in parts:
            if part not in opts[flag]:
                problems.append(
                    f'{where}: `neoHGT {stage} {flag} {part}` 里的取值不存在'
                    f'（该子命令只接受 {"|".join(sorted(opts[flag]))}）')
    for name in unbalanced:
        problems.append(f'{name} 的代码围栏不成对：围栏判定会漏掉一半示例，'
                        '取值核对对这些行是瞎的')
    if not pairs:
        problems.append('文档里一个带取值的示例命令都没读到：判据不能空跑')
    elif len(judged) < 10:
        problems.append(f'只判得了 {len(judged)} 个字面取值（共 {len(pairs)} 个）：'
                        '多半是围栏或分词读丢了，不是文档真的没写值')
    plan_judged = [j for j in judged if j[0].startswith(plan_label + ':')] \
        if plan_label else []
    return {'pairs': len(pairs), 'judged': len(judged),
            'placeholders': [list(p) for p in placeholders],
            'unbalanced_fences': unbalanced,
            'plan': plan or '', 'plan_state': ('judged' if plan_label else
                                               'not checked (no --plan)'),
            'plan_judged': len(plan_judged),
            'values': [[w, s, f, v] for w, s, f, v, _ in judged],
            'choice_stages': sorted(choices),
            }, sorted(set(problems))


def cli_subcommands(code_dir):
    """The sub-command names of the entry script (its own dispatch table)."""
    path = join(code_dir, 'scripts', 'neoHGT')
    if not isfile(path):
        return ()
    with open(path, encoding='utf-8') as f:
        src = f.read()
    m = re.search(r'module_names\s*=\s*\(([^)]*)\)', src)
    return tuple(re.findall(r"'([a-z-]+)'", m.group(1))) if m else ()


def doc_command_flags(code_dir, subs):
    """{flag: [where]} for every flag our shipped docs use on a neoHGT command.

    `surfaces` maps a sub-command to the flags that sub-command accepts; a flag
    that exists somewhere else in the CLI is still a broken example, because
    argparse builds one sub-parser per stage.
    """
    files = [join(code_dir, 'README.md')]
    docs = join(code_dir, 'docs')
    if isdir(docs):
        files += [join(docs, f) for f in sorted(os.listdir(docs))
                  if f.endswith('.md')]
    used = {}
    for path in files:
        if not isfile(path):
            continue
        with open(path, encoding='utf-8') as f:
            text = f.read().replace('\\\n', ' ')
        for n, line in enumerate(text.splitlines(), 1):
            for stage, flags in invocation_flags(line, subs):
                for flag in flags:
                    used.setdefault(flag, []).append(
                        f'{basename(path)}:{n} ({stage})')
    return used


def declared_doc_usages(ledger, subs, surfaces):
    """Rows the deviation table quotes on purpose, as (where, stage, flag).

    附录 C/D 命令形态偏差 exists to show what we deliberately did not copy, so a
    row there names options our CLI has never accepted -- that text *is* the
    evidence, and rewriting it to satisfy a checker would destroy the point.
    The exemption is tied to the exact line: excusing `analyze --out` because the
    deviation table quotes it must not also forgive the same mistake in a real
    example two paragraphs away.
    """
    if not isfile(ledger):
        return set()
    with open(ledger, encoding='utf-8') as f:
        lines = f.read().splitlines()
    marker = '## 附录 C/D 命令形态偏差'
    if marker not in lines:
        return set()
    start = lines.index(marker)
    end = next((i for i in range(start + 1, len(lines))
                if lines[i].startswith('## ')), len(lines))
    declared = set()
    for n in range(start, end):
        for stage, flags in invocation_flags(lines[n], subs):
            for flag in flags:
                if flag not in surfaces.get(stage, set()):
                    declared.add((f'{basename(ledger)}:{n + 1} ({stage})',
                                  stage, flag))
    return declared


def doc_flag_problems(used, surfaces, declared=()):
    """Examples naming an option the sub-command in that very example lacks.

    Judged per (sub-command, flag) pair: a flag used correctly under `gtdb` and
    wrongly under `analyze` must still be reported, which an aggregation by flag
    alone would forgive.
    """
    pairs = {}
    for flag, where in used.items():
        for w in where:
            stage = w.rsplit('(', 1)[1].rstrip(')')
            if flag in surfaces.get(stage, set()) or (w, stage, flag) in declared:
                continue
            pairs.setdefault((stage, flag), []).append(w)
    return [
        f'{where[0].rsplit(":", 1)[0]} 里的 `neoHGT {stage} … {flag}` '
        f'这个选项在 {stage} 这一档并不存在（{len(where)} 处，例如 {where[0]}）'
        for (stage, flag), where in sorted(pairs.items())]


def plan_inventory(plan_text):
    """Every plan-derived id list, in one dict.

    The registers can only be checked if the id inventory is known. That used to
    mean the plan document had to be present, which it is not in a repository
    checkout -- so the inventory can also be frozen to a committed digest.
    """
    return {
        'test_ids': sorted(set(ID_RE.findall(plan_text))),
        'tasks': plan_task_ids(plan_text),
        'defects': plan_defect_ids(plan_text),
        'deliverables': plan_deliverable_ids(plan_text),
        'risks': plan_risk_ids(plan_text),
        'dods': plan_dod_ids(plan_text),
        'decisions': plan_adr_ids(plan_text),
        'verification': plan_verification_ids(plan_text),
        'assumptions': plan_assumption_ids(plan_text),
        'revisions': plan_revision_ids(plan_text),
        'milestones': plan_milestone_ids(plan_text),
        'section7_flags': section7_flags(plan_text),
        'appendix_flags': appendix_flags(plan_text),
        'appendix_b_symbols': [list(x) for x in appendix_b_symbols(plan_text)],
    }


def write_inventory(path, plan_text, plan_name):
    """Freeze plan_inventory() to a committed digest."""
    import time
    data = {'source': basename(str(plan_name)),
            'generated_utc': time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime()),
            'counts': {}, 'inventory': plan_inventory(plan_text)}
    for k, v in data['inventory'].items():
        data['counts'][k] = len(v)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2, sort_keys=True)
    return path


def load_inventory(path):
    """Read a digest written by write_inventory()."""
    with open(path, encoding='utf-8') as f:
        data = json.load(f)
    inv = data['inventory']
    inv['appendix_b_symbols'] = [tuple(x) for x in inv['appendix_b_symbols']]
    return inv, data


def guarded_out(path, code_dir):
    """Refuse report/digest paths outside the project workspace.

    This tool is itself run by contributors, and project rule 11 says nothing may
    be written outside the workspace root. Command-line redirection is caught by
    tools/run.sh, but an `--out /tmp/x` argument is not a redirection, so the
    writers guard their own targets.
    """
    root = os.environ.get('NEOHGT_ROOT') or dirname(abspath(code_dir))
    real = os.path.realpath(os.path.expanduser(path))
    if not real.startswith(abspath(root) + os.sep):
        raise SystemExit(
            f'refusing to write outside the workspace: {path} -> {real}\n'
            f'workspace root: {abspath(root)}')
    return real


# --- rule 19: config route vs CLI route -------------------------------------
# `util.get_config(obj, attr, entry, func=None)` stores the YAML loader's object
# verbatim unless the call site passes a converter, while argparse always converts
# the *string* it received with the option's `type=`. So a value can arrive as two
# different objects depending on how it was supplied -- round 45 met this as
# `bandwidth: 1` (config) refusing while `--bandwidth 1` ran, and round 46 met its
# mirror: `--fetch-retries 3` reached `range(...)` as a string and raised
# TypeError, while the config file's int worked.
#
# Most key pairs are harmless because something downstream normalises them, and
# that is a legitimate answer -- but it has to be *said*. Each exempt key names the
# coercion that saves it; a key with neither a converter nor a stated neutralizer
# is a problem.
#
# What counts as divergent is a modelling choice, stated here because it is not
# obvious: `int` and `float` are treated as the *same* class (both are numbers, and
# `/ 100`, `<=`, `if val and val < 1` give the same answer either way), while
# `bool`/`str` versus a numeric option counts as divergent. Round 45's bandwidth is
# the honest counter-example to the first half -- there an int versus float did crash,
# because `perform_kde()` asked `isinstance(bw, float)` outright. That is a property
# of the consumer, not of the option, and the ledger cannot see it; it is covered
# instead by `_bandwidth_in_range()`, which normalises every incoming form, and by
# the keyword pass-through it leaves alone. A future key whose consumer really does
# branch on int-versus-float has to be handled the same way: fixed in code, not
# excused here.
# An exemption is a *claim about code*, so it is written as one: which module, which
# callable normalises it, which attribute that callable must actually touch, and the
# prose a human reads. Rule 19b then tries to falsify the claim -- the callable has to
# exist, has to name this attribute, and has to run after the `get_config` read.
#
# Ordering cannot be located by line for the two `fetch.*` keys, because their read is
# loop-form (`f'fetch.{key}'`): the AST knows a key list, not which iteration wrote
# which attribute. Rather than let "unknown" pass, those claims carry `proven_by`: the
# name of a test that pins the behaviour by execution, which this rule also checks
# exists in the suite. Round 48's probe measured all of this first
# (缓存测试/reports/R49_neutralizer_soundness.json), including the two ways a check of
# this shape goes wrong: reading only the call's arguments calls a real normaliser
# missing, while reading every string constant in the callee lets a claim pass on its
# own error-message wording.
ROUTE_CLAIMS = {
    'predict.weighted': ('analyze', 'arg2bool', 'weighted',
                         'arg2bool() over (weighted, orphans, self_low, donor_name), analyze.py set_parameters'),
    'predict.orphans': ('analyze', 'arg2bool', 'orphans',
                        'arg2bool() over (weighted, orphans, self_low, donor_name), analyze.py set_parameters'),
    'predict.selflow': ('analyze', 'arg2bool', 'self_low',
                        'arg2bool() over (weighted, orphans, self_low, donor_name), analyze.py set_parameters'),
    'taxonomy.capital': ('search', 'arg2bool', 'tax_capital',
                         'arg2bool() over (tax_unique, tax_capital, tax_latin), search.py args_wf'),
    'taxonomy.latin': ('search', 'arg2bool', 'tax_latin',
                       'arg2bool() over (tax_unique, tax_capital, tax_latin), search.py args_wf'),
    'taxonomy.unique': ('search', 'arg2bool', 'tax_unique',
                        'arg2bool() over (tax_unique, tax_capital, tax_latin), search.py args_wf'),
    'fetch.queries': ('search', 'coerce_fetch_counts', 'fetch_queries',
                      'coerce_fetch_counts() in the GTDB region of search.py (round 46: this was a real TypeError, now converted)'),
    'fetch.retries': ('search', 'coerce_fetch_counts', 'fetch_retries',
                      'coerce_fetch_counts() in the GTDB region of search.py (round 46: this was a real TypeError, now converted)'),
}
# Which claims rest on an executed test instead of a line order, and which test.
ROUTE_ORDER_PROOF = {
    'fetch.queries': 'test_the_wrapper_runs_the_coercion_on_every_route',
    'fetch.retries': 'test_the_wrapper_runs_the_coercion_on_every_route',
}
ROUTE_NEUTRALIZERS = {k: v[3] for k, v in ROUTE_CLAIMS.items()}
NO_DEST = '<no CLI option>'

# Every module a written reason may cite by name. The list is wider than
# ROUTE_CLAIMS' own modules because the reasons point at helpers (`arg2bool()`
# lives in util.py, not in the stage that reads the key).
CLAIM_MODULES = ('search', 'database', 'analyze', 'util', 'gtdb')
CITE_MODULE_RE = re.compile(r'\b([a-z_]+)\.py\s+([A-Za-z_]\w*)')
CITE_CALL_RE = re.compile(r'\b([A-Za-z_]\w*)\(\)')


def _citation_problems(note, defs_by_module, missing_modules):
    """Do the locations a reason names exist?

    A reason is prose, and prose rots quietly: three of the eight shipped reasons
    said "search.py set_parameters" while search.py's argument handler is
    `args_wf` -- nothing the normaliser lookup would notice, because the lookup
    follows the call, not the sentence. So the sentence is read too: `mod.py name`
    must be defined in that module, and a bare `name()` must be defined somewhere
    in CLAIM_MODULES. A module that is not in the tree at all is reported as
    unjudged rather than as a clean pass -- "could not look" is not "looked and
    was fine".
    """
    bad, unjudged = [], []
    all_defs = {n for defs in defs_by_module.values() for n in defs}
    seen = 0
    for match in CITE_MODULE_RE.finditer(note):
        seen += 1
        mod, name = match.group(1), match.group(2)
        if mod in defs_by_module:
            if name not in defs_by_module[mod]:
                bad.append(f'{mod}.py {name}')
        else:
            unjudged.append(f'{mod}.py {name}')
    for match in CITE_CALL_RE.finditer(note):
        seen += 1
        name = match.group(1)
        if name in all_defs:
            continue
        if missing_modules:
            unjudged.append(f'{name}()')
        else:
            bad.append(f'{name}()')
    return bad, unjudged, seen


def _slot_parts(node, target):
    """Template parts: strings, None for a `{target}` slot, ('dyn', attr) otherwise."""
    if isinstance(node, ast.Constant):
        return [str(node.value)]
    if not isinstance(node, ast.JoinedStr):
        return None
    parts = []
    for v in node.values:
        if isinstance(v, ast.Constant):
            parts.append(str(v.value))
            continue
        expr = v.value if isinstance(v, ast.FormattedValue) else None
        if isinstance(expr, ast.Name) and expr.id == target:
            parts.append(None)
        elif (isinstance(expr, ast.Call)
              and isinstance(expr.func, ast.Attribute)
              and expr.func.attr == 'replace'
              and isinstance(expr.func.value, ast.Name)
              and expr.func.value.id == target
              and len(expr.args) == 2
              and all(isinstance(a, ast.Constant) for a in expr.args)):
            parts.append(('replace', str(expr.args[0].value),
                          str(expr.args[1].value)))
        elif (isinstance(expr, ast.Attribute)
              and isinstance(expr.value, ast.Name)
              and expr.value.id == 'self'):
            parts.append(('dyn', expr.attr))
        else:
            parts.append(('dyn', None))
    return parts


def _render(parts, key):
    """A template with one key substituted, or None if a runtime slot is inside."""
    out = ''
    for p in parts:
        if p is None:
            out += str(key)
        elif isinstance(p, tuple) and p[0] == 'dyn':
            return None
        elif isinstance(p, tuple):
            out += str(key).replace(p[1], p[2])
        else:
            out += p
    return out


def _render_entries(parts, key, choices):
    """Every entry a template denotes for `key`, expanding runtime slots.

    `f'{self.method}.{key}'` (search.py:305) names a config *section* only at
    runtime, so an earlier draft refused the call site and rule 19 silently
    audited 58 keys while the census enumerated 67 -- the 9 missing ones being
    `diamond./blast./remote.` + queries/maxchars/extrargs. The refusal was not
    necessary: the set of sections that slot can take is written down in the
    option that sets it (`-m|--method`'s `{'choices': [...]}`), which the AST
    already reads. Expanding over the declared choices is reading the code, not
    guessing it; and if the option ever stops declaring them, the call site goes
    back to `unresolved` instead of being invented.
    """
    slots = [i for i, p in enumerate(parts)
             if isinstance(p, tuple) and p[0] == 'dyn']
    if not slots:
        one = _render(parts, key)
        return [one] if one else []
    attr = parts[slots[0]][1]
    values = (choices or {}).get(attr) or []
    if not values:
        return []
    out, basis = [], []
    for value in values:
        replaced = list(parts)
        replaced[slots[0]] = str(value)
        entry = _render(replaced, key)
        if entry:
            out.append(entry)
            basis.append(f'self.{attr} in {sorted(values)}')
    return out


def route_cli_choices(code_dir):
    """{module: {dest: [declared choices]}} from the same `arguments` tables."""
    out = {}
    for module in ('database', 'search', 'analyze'):
        path = join(code_dir, 'neoHGT', f'{module}.py')
        if not isfile(path):
            continue
        with open(path, encoding='utf-8', errors='replace') as fh:
            tree = ast.parse(fh.read())
        table = {}
        for node in ast.walk(tree):
            if not (isinstance(node, (ast.List, ast.Tuple)) and node.elts):
                continue
            for item in node.elts:
                if not (isinstance(item, (ast.List, ast.Tuple)) and item.elts):
                    continue
                first = item.elts[0]
                if not (isinstance(first, ast.Constant)
                        and isinstance(first.value, str)
                        and first.value.startswith('-')):
                    continue
                longs = [f.strip() for f in first.value.split('|')
                         if f.strip().startswith('--')]
                if not longs:
                    continue
                dest = max(longs, key=len)[2:].replace('-', '_')
                for sub in item.elts[1:]:
                    if not isinstance(sub, ast.Dict):
                        continue
                    for k, v in zip(sub.keys, sub.values):
                        if (isinstance(k, ast.Constant)
                                and k.value == 'choices'
                                and isinstance(v, (ast.List, ast.Tuple))):
                            vals = [e.value for e in v.elts
                                    if isinstance(e, ast.Constant)]
                            if vals:
                                table[dest] = vals
        out[module] = table
    return out


def route_reads(code_dir):
    """(rows, unresolved) for every `get_config(...)` call in the stage modules.

    Rows are (module, attr, entry, converter_name). Attr and entry are paired
    **per loop key**: an earlier draft paired every attr with every entry, which
    is precisely the mistake round 46's census made and then documented -- it
    reported `predict.bandwidth` under the attribute `weighted`, and a rule built
    on that pairing would accuse the wrong option.
    """
    rows, unresolved = [], []
    all_choices = route_cli_choices(code_dir)
    for module in ('database', 'search', 'analyze'):
        path = join(code_dir, 'neoHGT', f'{module}.py')
        mod_choices = all_choices.get(module, {})
        if not isfile(path):
            continue
        with open(path, encoding='utf-8', errors='replace') as fh:
            tree = ast.parse(fh.read())
        # innermost constant-list loop around each node: (target, keys, assigns)
        loop_of = {}
        for node in ast.walk(tree):
            if not isinstance(node, ast.For):
                continue
            it = node.iter
            if isinstance(it, ast.Call) and it.args:
                it = it.args[0]
            if isinstance(it, (ast.List, ast.Tuple)):
                keys = [e.value for e in it.elts if isinstance(e, ast.Constant)]
                target = (node.target.id if isinstance(node.target, ast.Name)
                          else None)
                assigns = {}
                for stmt in node.body:
                    if (isinstance(stmt, ast.Assign)
                            and len(stmt.targets) == 1
                            and isinstance(stmt.targets[0], ast.Name)):
                        name = stmt.targets[0].id
                        parts = _slot_parts(stmt.value, target)
                        # assigned twice in one body is not resolvable per key
                        assigns[name] = (parts if name not in assigns
                                         or assigns[name] is None else None) \
                            if parts is not None else None
                for sub in ast.walk(node):
                    loop_of[id(sub)] = (target, keys, assigns)
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name)
                    and node.func.id == 'get_config'):
                continue
            args = node.args
            attr = args[1] if len(args) > 1 else None
            entry = args[2] if len(args) > 2 else None
            func = args[3] if len(args) > 3 else None
            conv = (func.id if isinstance(func, ast.Name)
                    else 'lambda' if isinstance(func, ast.Lambda) else None)
            target, keys, assigns = loop_of.get(id(node), (None, [], {}))
            attr_parts = (_slot_parts(attr, target)
                          if not isinstance(attr, ast.Constant) else None)
            entry_parts = _slot_parts(entry, target)
            if entry_parts is None:
                unresolved.append(f'{module}.py:{getattr(node, "lineno", "?")} '
                                  'entry is not a literal or f-string')
                continue
            pairs = []

            def entries(parts, key):
                """[(entry, basis)] -- a template denoting one or more config keys."""
                got = _render_entries(parts, key, mod_choices)
                dyn = [p for p in parts if isinstance(p, tuple) and p[0] == 'dyn']
                basis = ''
                if dyn:
                    basis = ('self.%s expanded over that option\'s declared '
                             'choices' % (dyn[0][1] or '?'))
                return [(g, basis) for g in got]
            if isinstance(attr, ast.Constant):
                if isinstance(entry, ast.Constant):
                    pairs = [(str(attr.value), str(entry.value), '')]
                else:
                    # A literal attribute with a templated entry: pairing it with
                    # every loop key would be a guess, so it is reported instead of
                    # invented. No such call exists today -- if one appears, the
                    # `unresolved` line says which and the rule stays honest.
                    unresolved.append(
                        f'{module}.py:{getattr(node, "lineno", "?")} literal '
                        'attr with a templated entry (not guessed)')
            elif isinstance(attr, ast.Name) and attr.id == target:
                for k in keys:
                    for e, basis in entries(entry_parts, k):
                        pairs.append((k, e, basis))
            elif (isinstance(attr, ast.Name) and assigns.get(attr.id)
                  is not None):
                # `attr = f'tax_{key}'` one statement above the call: following the
                # assignment is what the census does, and a gate that stops short
                # of it silently audits fewer keys than the evidence it claims to
                # encode (measured: 56 keys here against 67 in the census, and the
                # missing ones were exactly the `tax_*` filters).
                aparts = assigns[attr.id]
                for k in keys:
                    a = _render(aparts, k)
                    if not a:
                        continue
                    for e, basis in entries(entry_parts, k):
                        pairs.append((a, e, basis))
            elif attr_parts is not None:
                # the attr is itself a template: `f'donor_{key}'` inline
                for k in keys:
                    a = _render(attr_parts, k)
                    if not a:
                        continue
                    for e, basis in entries(entry_parts, k):
                        pairs.append((a, e, basis))
            pairs = [(a, e, b) for a, e, b in pairs if a and e]
            if pairs:
                rows.extend((module, a, e, conv, b) for a, e, b in pairs)
            else:
                unresolved.append(f'{module}.py:{getattr(node, "lineno", "?")} '
                                  'attr/entry not resolvable from the AST')
    return rows, unresolved


def route_cli_types(code_dir):
    """{module: {dest: type-name}} read from each stage's `arguments` table.

    The table passes argparse keyword arguments as a *trailing dict literal*
    (`['--evalue', 'maximum E-value cutoff', {'type': float}]`), not as `ast`
    keywords; an earlier census read `item.keywords`, saw no type anywhere, and
    inflated the divergences from 11 to 36.
    """
    out = {}
    for module in ('database', 'search', 'analyze'):
        path = join(code_dir, 'neoHGT', f'{module}.py')
        if not isfile(path):
            continue
        with open(path, encoding='utf-8', errors='replace') as fh:
            tree = ast.parse(fh.read())
        types = {}
        for node in ast.walk(tree):
            if not (isinstance(node, (ast.List, ast.Tuple)) and node.elts):
                continue
            for item in node.elts:
                if not (isinstance(item, (ast.List, ast.Tuple)) and item.elts):
                    continue
                first = item.elts[0]
                if not (isinstance(first, ast.Constant)
                        and isinstance(first.value, str)
                        and first.value.startswith('-')):
                    continue
                flags = [f.strip() for f in first.value.split('|')]
                longs = [f for f in flags if f.startswith('--')]
                if not longs:
                    continue
                dest = max(longs, key=len)[2:].replace('-', '_')
                kw = {}
                for sub in item.elts[1:]:
                    if isinstance(sub, ast.Dict):
                        kw.update({str(k.value): v for k, v in
                                   zip(sub.keys, sub.values)
                                   if isinstance(k, ast.Constant)})
                typ = kw.get('type')
                if isinstance(typ, ast.Name):
                    types[dest] = typ.id
                elif isinstance(kw.get('action'), ast.Constant):
                    act = str(kw['action'].value)
                    types[dest] = ('bool' if act in ('store_true', 'store_false')
                                   else 'action:' + act)
                else:
                    types.setdefault(dest, None)
        out[module] = types
    return out


def route_ledger_emit(code_dir, out_path):
    """Write the committed route ledger: one record per config key the code reads.

    Emitting needs the config file's *values* (that is what decides whether a read
    can diverge at all), so this path requires PyYAML and is run where the
    scientific environment exists -- the round gates. Checking the ledger never
    needs it: a checkout without PyYAML still verifies that every divergence in the
    committed ledger is either converted or explained.
    """
    import yaml
    cfg_path = join(code_dir, 'neoHGT', 'config.yml')
    with open(cfg_path, encoding='utf-8') as fh:
        raw = yaml.safe_load(fh)

    def flat(d, prefix=''):
        out = {}
        for k, v in (d or {}).items():
            key = f'{prefix}{k}'
            if isinstance(v, dict):
                out.update(flat(v, key + '.'))
            else:
                out[key] = v
        return out

    leaves = flat(raw)
    types = route_cli_types(code_dir)
    rows, unresolved = route_reads(code_dir)
    records = {}
    for module, attr, entry, conv, basis in rows:
        val = leaves.get(entry, '<absent>')
        cli = types.get(module, {}).get(attr, NO_DEST)
        if cli == NO_DEST:
            # no CLI flag names this attribute: there is no second route to
            # compare, and inventing one would manufacture a finding
            divergent = False
        elif val == '<absent>' or val is None or isinstance(val, (dict, list)):
            divergent = False
        elif cli in ('int', 'float'):
            divergent = not isinstance(val, (int, float)) or isinstance(val, bool)
        elif cli == 'bool':
            divergent = not isinstance(val, bool)
        elif cli is None:
            # an option with no `type=`: argparse hands over the raw string, so a
            # non-string in the config file is a divergence. This is the shape
            # round 46's bug had (`--fetch-retries 3` reaching `range()` as '3'),
            # so leaving it out of the definition is what let the ledger call the
            # only two harmful keys harmless.
            divergent = not isinstance(val, str)
        else:
            divergent = False
        records.setdefault(entry, {
            'attr': attr, 'modules': [], 'converter': conv,
            'cli_type': (None if cli is None else cli),
            'yaml_type': type(val).__name__,
            'in_config': val != '<absent>',
            'divergent': bool(divergent),
        })
        if basis:
            records[entry].setdefault('expansion_basis', basis)
        if module not in records[entry]['modules']:
            records[entry]['modules'].append(module)
    ledger = {'generated_from': code_dir,
              'note': ('one record per config key the stage modules read; a key '
                       'whose two routes hand over different objects needs a '
                       'converter or a named neutralizer (register rule 19)'),
              'unresolved_call_sites': unresolved,
              'keys': records}
    with open(out_path, 'w', encoding='utf-8') as fh:
        json.dump(ledger, fh, ensure_ascii=False, indent=1, sort_keys=True)
    return len(records)


def _callee_names(fn, tree):
    """Attributes a called function normalises inside itself.

    Narrow on purpose: `self.<attr>`, `getattr/setattr(self, '<attr>')` literals,
    and the keys of module-level containers the body references by name. Taking
    every string constant instead makes the check answer to wording -- round 48's
    probe passed `fetch.queries` because the string `fetch_queries` appears in the
    function's own error message, and failed `fetch.retries` because no message
    names it. Neither verdict would have meant anything.
    """
    names = set()
    for sub in ast.walk(fn):
        if (isinstance(sub, ast.Attribute)
                and isinstance(sub.value, ast.Name) and sub.value.id == 'self'):
            names.add(sub.attr)
        elif (isinstance(sub, ast.Call) and isinstance(sub.func, ast.Name)
              and sub.func.id in ('getattr', 'setattr', 'hasattr')
              and len(sub.args) > 1 and isinstance(sub.args[1], ast.Constant)
              and isinstance(sub.args[1].value, str)):
            names.add(sub.args[1].value)
    referenced = {n.id for n in ast.walk(fn) if isinstance(n, ast.Name)}
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for tgt in node.targets:
            if not (isinstance(tgt, ast.Name) and tgt.id in referenced):
                continue
            val = node.value
            if isinstance(val, ast.Dict):
                names.update(k.value for k in val.keys
                             if isinstance(k, ast.Constant)
                             and isinstance(k.value, str))
            elif isinstance(val, (ast.List, ast.Tuple, ast.Set)):
                names.update(e.value for e in val.elts
                             if isinstance(e, ast.Constant)
                             and isinstance(e.value, str))
    return names


def _loops_of(fn):
    """id(node) -> (loop target, keys, {assigned name: template}) inside `fn`.

    Shared by the normaliser search and the read-line search so they cannot
    disagree about which loop a call sits in, and so a templated attribute
    (`attr = f'tax_{key}'`) resolves the same way in both.
    """
    loops = {}
    for node in ast.walk(fn):
        if not isinstance(node, ast.For):
            continue
        it = node.iter
        if isinstance(it, ast.Call) and it.args:
            it = it.args[0]
        if isinstance(it, (ast.List, ast.Tuple)):
            keys = [e.value for e in it.elts if isinstance(e, ast.Constant)]
            target = node.target.id if isinstance(node.target, ast.Name) else None
            assigns = {}
            for stmt in node.body:
                if (isinstance(stmt, ast.Assign) and len(stmt.targets) == 1
                        and isinstance(stmt.targets[0], ast.Name)):
                    name = stmt.targets[0].id
                    parts = _slot_parts(stmt.value, target)
                    assigns[name] = (parts if name not in assigns
                                     and parts is not None else None)
            for sub in ast.walk(node):
                loops[id(sub)] = (target, keys, assigns)
    return loops


def _normaliser_calls(fn, want, defs, tree):
    """[(line, names)] for each call of `want` in `fn`, loop key lists included."""
    loops = _loops_of(fn)
    hits = []
    for node in ast.walk(fn):
        if not (isinstance(node, ast.Call)
                and ((isinstance(node.func, ast.Name)
                      and node.func.id == want)
                     or (isinstance(node.func, ast.Attribute)
                         and node.func.attr == want))):
            continue
        names = set()
        for sub in ast.walk(node):
            if (isinstance(sub, ast.Attribute)
                    and isinstance(sub.value, ast.Name)
                    and sub.value.id == 'self'):
                names.add(sub.attr)
            elif isinstance(sub, ast.Constant) and isinstance(sub.value, str):
                names.add(sub.value)
        target, keys, _asg = loops.get(id(node), (None, [], {}))
        if keys and target:
            names.update(f'__loop__:{k}' for k in keys)
        if isinstance(node.func, ast.Name) and node.func.id in defs:
            names.update(_callee_names(defs[node.func.id], tree))
        hits.append((getattr(node, 'lineno', 0), names))
    return hits


def _read_line(fn, attr, tree):
    """Where this attribute is read from config, if a line can be named.

    Three shapes: a literal attr, the loop variable itself, and a name assigned one
    statement above the call (`attr = f'tax_{key}'`). The third is what the first
    version of rule 19b could not see, and it reported three sound exemptions as
    unverifiable -- a check that fails because it cannot read the code is worse than
    no check, because the failure gets attributed to the code.
    """
    loops = _loops_of(fn)
    for node in ast.walk(fn):
        if not (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == 'get_config' and len(node.args) > 1):
            continue
        a = node.args[1]
        target, keys, assigns = loops.get(id(node), (None, [], {}))
        if isinstance(a, ast.Constant) and str(a.value) == attr:
            return getattr(node, 'lineno', 0)
        # An inline template attribute is checked *before* the Name-only guard:
        # the first version put the guard above it, so the branch below never ran
        # and three sound exemptions kept reporting as unverifiable.
        inline = _slot_parts(a, target)
        if inline and any(_render(inline, k) == attr for k in keys):
            return getattr(node, 'lineno', 0)
        if not isinstance(a, ast.Name):
            continue
        if a.id == target and attr in keys:
            return getattr(node, 'lineno', 0)
        parts = assigns.get(a.id)
        if parts and any(_render(parts, k) == attr for k in keys):
            return getattr(node, 'lineno', 0)
    return None


def route_claim_problems(code_dir, claimed):
    """Rule 19b: a written reason must be findable, cover the key, and run after the read."""
    problems, detail = [], {}
    trees, defs_by_module = {}, {}
    for module in sorted(set(CLAIM_MODULES) | {c[0] for c in ROUTE_CLAIMS.values()}):
        path = join(code_dir, 'neoHGT', f'{module}.py')
        if isfile(path):
            with open(path, encoding='utf-8', errors='replace') as fh:
                tree = ast.parse(fh.read())
            trees[module] = tree
            defs_by_module[module] = {n.name: n for n in ast.walk(tree)
                                       if isinstance(n, ast.FunctionDef)}
    missing_modules = sorted(set(CLAIM_MODULES) - set(defs_by_module))
    unjudged, citations_seen = [], 0
    test_path = join(code_dir, 'neoHGT', 'tests', 'test_gtdb_pipeline.py')
    test_text = ''
    if isfile(test_path):
        with open(test_path, encoding='utf-8', errors='replace') as fh:
            test_text = fh.read()
    checked = 0
    for entry, (module, want, attr, note) in ROUTE_CLAIMS.items():
        if entry not in claimed:
            continue
        checked += 1
        bad, skip, seen = _citation_problems(note, defs_by_module, missing_modules)
        citations_seen += seen
        for cite in bad:
            problems.append(f'路由豁免 {entry} 的理由写着 "{cite}"，但那个位置在代码里'
                            '没有定义：理由本身已经不成立（它辩护的东西不在这儿）')
        unjudged += skip
        tree, defs = trees.get(module), defs_by_module.get(module, {})
        found = None
        for fname, fn in defs.items():
            read_line = _read_line(fn, attr, tree)
            for lineno, names in _normaliser_calls(fn, want, defs, tree):
                if attr in names or f'__loop__:{attr}' in names:
                    order = None
                    if read_line:
                        order = lineno > read_line
                    found = {'function': fname, 'line': lineno,
                             'read_line': read_line, 'after_read': order}
                    break
            if found:
                break
        if not found:
            problems.append(
                f'路由豁免说 {entry} 由 {want}() 归一化，但在 {module}.py 里找不到'
                f'覆盖属性 {attr} 的调用点：豁免写的是不存在（或不生效）的东西')
            continue
        if found['after_read'] is False:
            problems.append(
                f'路由豁免 {entry} 的归一化 {want}()（{found["function"]}:'
                f'{found["line"]}）发生在配置读取（行 {found["read_line"]}）之前，'
                '先归一化后取值等于没归一化')
        elif found['after_read'] is None:
            proof = ROUTE_ORDER_PROOF.get(entry)
            if not proof or f'def {proof}(' not in test_text:
                problems.append(
                    f'路由豁免 {entry} 的归一化顺序无法从行号判定（读的是循环形键），'
                    f'而申明指向的用例 {proof or "（无）"} 不在测试里：'
                    '顺序既没被读到也没被测到，不能算已核')
        detail[entry] = found
    if checked and len(detail) != checked:
        problems.append(f'路由豁免核对只判定成功 {len(detail)}/{checked} 条：'
                        '剩下的没有被看，绿灯不算数')
    if checked and citations_seen == 0:
        problems.append(f'{checked} 条路由豁免的理由里一个代码位置都没提到：'
                        '理由没有可核对的东西，等于没有理由')
    return {'claims_checked': checked, 'verdicts': detail,
            'citations_seen': citations_seen,
            'citations_unjudged': sorted(set(unjudged)),
            'citations_missing_modules': missing_modules}, problems


# ---------------------------------------------------------------------------
# Rule 20: a default the documentation promises must be the value the code gives.
# ---------------------------------------------------------------------------
#
# Round 50 measured the surface before writing this (临时文件/r50_default_claims_probe
# .json, r50_default_delivery_census_v3.json, r50_promise_sites.json): 156 CLI
# options, 145 of them with no `default=` in their own table, so the promise lives
# in a sentence while the delivery lives in a chain -- two places maintained
# apart. The first defect found this way was real: `--donor-label` had no table
# default, no config key and no profile entry, and its consumer only acted on the
# literal 'prefixed', so a default GTDB run wrote the synthetic taxId
# (`105000000`) into hgts/<sample>.txt while 计划任务对照表 4.4 promised prefixed
# labels. It is now pinned by a claim below and by an executed fixture test.
#
# Claims are pinned to (module, function, attribute), not to an option name,
# because the same attribute is read with different fallbacks in different
# places: `keep_placeholder_names` is established from the profile in
# set_params_gtdb and read defensively with False in taxfilt_wf_gtdb. Taking
# "the first getattr site" would freeze the wrong number.
DEFAULT_PROMISES = [
    {'key': 'seq_source', 'flag': '--seq-source', 'module': 'database',
     'fn': 'set_source_parameters', 'attr': 'seq_source', 'promised': 'local',
     'shape': 'or_chain', 'docs': ('docs/决策落点对照.md', 'docs/方案修订记录核对.md'),
     'tokens': ('local',)},
    {'key': 'taxid_base', 'flag': '--gtdb-taxid-base', 'module': 'database',
     'fn': 'set_source_parameters', 'attr': 'id_base', 'promised': '100000000',
     'shape': 'or_chain', 'docs': ('docs/GTDB-数据库.md',),
     'tokens': ('100000000',)},
    {'key': 'gtdb_taxid_base', 'flag': '--gtdb-taxid-base', 'module': 'gtdb',
     'fn': None, 'attr': 'gtdb_taxid_base', 'promised': '100000000',
     'shape': 'table', 'docs': ('docs/GTDB-数据库.md',),
     'tokens': ('100000000',)},
    {'key': 'keep_placeholder', 'flag': '--keep-placeholder-names',
     'module': 'search', 'fn': 'set_params_gtdb', 'attr': 'keep_placeholder_names',
     'promised': 'yes', 'shape': 'profile_get', 'cfg_name': 'keep_placeholder',
     'docs': ('docs/GTDB-数据库.md', 'docs/风险应对登记.md'), 'tokens': ('yes',)},
    {'key': 'require_staxids', 'flag': '--require-staxids', 'module': 'search',
     'fn': 'set_params_gtdb', 'attr': 'require_staxids', 'promised': 'yes',
     'shape': 'mode_default', 'docs': ('docs/GTDB-数据库.md',),
     'tokens': ('报错', 'GTDB')},
    {'key': 'donor_label', 'flag': '--donor-label', 'module': 'analyze',
     'fn': 'write_hgt_list_gtdb', 'attr': 'donor_label', 'promised': 'prefixed',
     'shape': 'none_guard', 'docs': ('docs/计划任务对照表.md',),
     'tokens': ('prefixed', '前缀')},
    {'key': 'allow_network', 'flag': '--allow-network', 'module': 'database',
     'fn': 'set_source_parameters', 'attr': 'allow_network', 'promised': 'no',
     'shape': 'is_true', 'docs': ('docs/验收标准对照.md', 'docs/风险应对登记.md'),
     'tokens': ('关',)},
    {'key': 'close_size', 'flag': '--close-size', 'module': 'analyze',
     'fn': 'apply_gtdb_grouping', 'attr': 'close_size',
     'promised': {'gtdb': '200', 'gtdb_archaea': '50'}, 'shape': 'setattr_profile',
     'cfg_name': 'closesize', 'cfg_path': ('grouping', 'profiles'),
     'docs': ('docs/GTDB-数据库.md',), 'tokens': ('200', '50')},
]
DEFAULT_RULE_MODULES = ('database', 'search', 'analyze', 'gtdb', 'util')
DEFAULT_WORD_RE = re.compile(r'默认|default', re.I)


def _norm_default(val):
    """One comparison space for `yes`/`no`, numbers and strings."""
    if isinstance(val, bool):
        return 'yes' if val else 'no'
    if val is None:
        return 'none'
    return str(val).strip().strip('\'"').lower()


def _safe_literal(node):
    try:
        return ast.literal_eval(node)
    except Exception:
        return None


def _fn_def(tree, name):
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    return None


def _cfg_deep(cfg, key):
    """Find `key` in the shipped config. Two different values make it unjudgeable."""
    found = []

    def walk(node, path):
        if isinstance(node, dict):
            for k, val in node.items():
                if k == key:
                    found.append(('.'.join(path + [str(k)]), val))
                walk(val, path + [str(k)])
        elif isinstance(node, list):
            for i, val in enumerate(node):
                walk(val, path + [str(i)])
    walk(cfg, [])
    if not found:
        return None, None
    values = {json.dumps(v, sort_keys=True, default=str) for _p, v in found}
    if len(values) > 1:
        return 'AMBIGUOUS', [p for p, _v in found]
    return found[0][1], found[0][0]


def _assign_to(fn, attr):
    for node in ast.walk(fn):
        if (isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Attribute) and t.attr == attr
                        and isinstance(t.value, ast.Name) and t.value.id == 'self'
                        for t in node.targets)):
            yield node


def _resolve_or_chain(module, fn, attr, trees, cfg):
    """`self.x = getattr(args,'x',None) or cfg.get('x') or 'literal'`.

    A falsy operand does not end the chain -- `False or 'local'` delivers 'local' --
    so the first *truthy* resolvable operand wins, and only if none is truthy does
    the last operand count (`self.x = a or False` really does deliver False).
    """
    for node in _assign_to(fn, attr):
        val = node.value
        operands = (val.values if isinstance(val, ast.BoolOp)
                    and isinstance(val.op, ast.Or) else [val])
        last = None
        for op in operands:
            lit = _safe_literal(op)
            if lit is not None and lit != '' and lit is not False:
                return _norm_default(lit), f'{module}.py:{node.lineno} 字面量 `{ast.unparse(op)[:30]}`'
            last = (op, lit)
            if (isinstance(op, ast.Call) and getattr(op.func, 'attr', '') == 'get'
                    and op.args and isinstance(op.args[0], ast.Constant)):
                key = op.args[0].value
                if cfg is None:
                    continue
                got, path = _cfg_deep(cfg, key)
                if got == 'AMBIGUOUS':
                    return None, f'config 里 {key} 有多个不同取值（{path}）'
                if got:
                    return _norm_default(got), f'config `{path}`（{module}.py:{node.lineno}）'
            if isinstance(op, ast.Attribute) and isinstance(op.value, ast.Name):
                mod, const = op.value.id, op.attr
                t2 = trees.get(mod)
                if t2 is not None:
                    for n2 in ast.walk(t2):
                        if (isinstance(n2, ast.Assign)
                                and any(getattr(t, 'id', None) == const
                                        for t in n2.targets)):
                            lv = _safe_literal(n2.value)
                            if lv:
                                return _norm_default(lv), f'{mod}.{const} = {lv}'
        if last and last[1] is not None:
            return _norm_default(last[1]), \
                f'{module}.py:{node.lineno} 末位回退 `{ast.unparse(last[0])[:30]}`'
        break
    return None, None


def _resolve_profile_get(module, fn, attr, cfg_name, cfg):
    for node in _assign_to(fn, attr):
        for call in ast.walk(node.value):
            if (isinstance(call, ast.Call) and getattr(call.func, 'attr', '') == 'get'
                    and call.args and isinstance(call.args[0], ast.Constant)
                    and call.args[0].value == cfg_name):
                if cfg is not None:
                    got, path = _cfg_deep(cfg, cfg_name)
                    if got == 'AMBIGUOUS':
                        return None, f'config 里 {cfg_name} 有多个不同取值（{path}）'
                    if got is not None:
                        return _norm_default(got), f'config `{path}`（{module}.py:{node.lineno}）'
                if len(call.args) > 1:
                    dflt = _safe_literal(call.args[1])
                    if dflt is not None:
                        return _norm_default(dflt), \
                            f'{module}.py:{node.lineno} 代码内回退 {dflt!r}'
        break
    return None, None


def _resolve_ternary_none(module, fn, attr, cfg):
    """`self.x = True if getattr(args,'x',None) is None else convert(...)`."""
    for node in _assign_to(fn, attr):
        ifexp = node.value
        if isinstance(ifexp, ast.IfExp) and isinstance(ifexp.test, ast.Compare):
            return _norm_default(_safe_literal(ifexp.body)), \
                f'{module}.py:{node.lineno} "没提这个选项时" 的那一支'
    return None, None


def _resolve_mode_default(module, fn, attr, cfg):
    """`if getattr(args,'x',None) is None: self.x = is_gtdb_mode(self)`.

    The default is a mode predicate: yes on the GTDB/hybrid routes, upstream
    behaviour (no) on the NCBI route. The register documents defaults per
    route, so the promise is judged on the GTDB route it describes.
    """
    for node in _assign_to(fn, attr):
        if isinstance(node.value, ast.Call) and \
                getattr(node.value.func, 'id', '') == 'is_gtdb_mode':
            return 'yes', f'{module}.py:{node.lineno} GTDB 路线默认'
    return None, None


def _resolve_none_guard(module, fn, attr, cfg):
    """`if something is None: self.x = 'literal'`."""
    guard = False
    for node in ast.walk(fn):
        if isinstance(node, ast.If) and any(
                isinstance(c, ast.Compare) and isinstance(c.ops[0], ast.Is)
                and isinstance(c.comparators[0], ast.Constant)
                and c.comparators[0].value is None
                for c in ast.walk(node.test)):
            guard = True
    if not guard:
        return None, None
    for node in _assign_to(fn, attr):
        lit = _safe_literal(node.value)
        if lit is not None:
            return _norm_default(lit), f'{module}.py:{node.lineno} 显式补默认值'
    return None, None


def _resolve_is_true(module, fn, attr, cfg):
    """`self.x = arg2bool(getattr(args,'x',None)) is True` delivers False."""
    for node in _assign_to(fn, attr):
        if isinstance(node.value, ast.Compare) and isinstance(node.value.ops[0], ast.Is):
            return 'no', f'{module}.py:{node.lineno} `... is True`：没提这个选项就是 no'
    return None, None


def _resolve_setattr_profile(module, fn, attr, promise, cfg):
    """A per-domain profile: the loop `for attr, name in ((...))` plus config.

    The value is not written at the assignment -- `setattr(self, attr, prof[name])`
    names a loop variable -- so the shape to look for is a For whose iterable
    mentions this attribute, containing a `setattr(self, <var>, ...)`.
    """
    found_loop = False
    for node in ast.walk(fn):
        if not isinstance(node, ast.For):
            continue
        mentioned = any(isinstance(c, ast.Constant) and c.value == attr
                        for c in ast.walk(node.iter))
        if not mentioned:
            continue
        if any(isinstance(c, ast.Call) and getattr(c.func, 'id', '') == 'setattr'
               and len(c.args) == 3 and isinstance(c.args[0], ast.Name)
               and c.args[0].id == 'self'
               for stmt in node.body for c in ast.walk(stmt)):
            found_loop = True
    if not found_loop:
        return None, None
    if cfg is None:
        return None, '没有 PyYAML，档位值读不出来'
    section = cfg
    for part in promise.get('cfg_path', ()):
        section = section.get(part) if isinstance(section, dict) else None
    if not isinstance(section, dict):
        return None, f'config 里找不到 {".".join(promise.get("cfg_path", ())) or "该段"}'
    return ({k: _norm_default((v or {}).get(promise['cfg_name']))
             for k, v in section.items() if isinstance(v, dict)},
            f'config `{".".join(promise.get("cfg_path", ()))}.*'
            f'.{promise["cfg_name"]}`（{module}.py 循环里的 setattr）')


def _resolve_table(promise, trees):
    """A default written in the option table: `default=5`, or `default=CONST`.

    The constant form is what `gtdb --gtdb-taxid-base` actually writes
    (`neoHGT/gtdb.py:4144`, with `DEFAULT_ID_BASE = 100_000_000` at gtdb.py:54), and
    it was the one delivery shape rule 20 could not resolve when the rule first
    shipped -- a gap noted in ACCEPTANCE §63.4 and closed here rather than left as a
    comment.
    """
    module, flag = promise['module'], promise['flag']
    tree = trees.get(module)
    if tree is None:
        return None, f'{module}.py 不在这棵树上'
    for lst in option_tables(tree):
        for item in lst.elts:
            if (not isinstance(item, (ast.List, ast.Tuple)) or len(item.elts) < 3
                    or not isinstance(item.elts[0], ast.Constant)
                    or not isinstance(item.elts[2], ast.Dict)):
                continue
            if flag not in re.findall(r'--?[a-z][\w-]*', str(item.elts[0].value)):
                continue
            for k, val in zip(item.elts[2].keys, item.elts[2].values):
                if getattr(k, 'value', None) != 'default':
                    continue
                lit = _safe_literal(val)
                if lit is not None:
                    return _norm_default(lit), \
                        f'{module}.py:{item.lineno} 表内 default={lit!r}'
                if isinstance(val, ast.Name):
                    for n in ast.walk(tree):
                        if (isinstance(n, ast.Assign)
                                and any(getattr(t, 'id', None) == val.id
                                        for t in n.targets)):
                            lv = _safe_literal(n.value)
                            if lv is not None:
                                return _norm_default(lv), \
                                    f'{module}.py:{n.lineno} default={val.id} = {lv!r}'
                    return None, f'{module}.py 里的常量 {val.id} 解析不出来'
    return None, None


SHAPES = {'table': lambda p, trees, cfg: _resolve_table(p, trees),
          'or_chain': lambda p, trees, cfg: _resolve_or_chain(
              p['module'], _fn_def(trees.get(p['module']), p['fn']), p['attr'],
              trees, cfg),
          'profile_get': lambda p, trees, cfg: _resolve_profile_get(
              p['module'], _fn_def(trees.get(p['module']), p['fn']), p['attr'],
              p['cfg_name'], cfg),
          'ternary_none': lambda p, trees, cfg: _resolve_ternary_none(
              p['module'], _fn_def(trees.get(p['module']), p['fn']), p['attr'], cfg),
          'mode_default': lambda p, trees, cfg: _resolve_mode_default(
              p['module'], _fn_def(trees.get(p['module']), p['fn']), p['attr'], cfg),
          'none_guard': lambda p, trees, cfg: _resolve_none_guard(
              p['module'], _fn_def(trees.get(p['module']), p['fn']), p['attr'], cfg),
          'is_true': lambda p, trees, cfg: _resolve_is_true(
              p['module'], _fn_def(trees.get(p['module']), p['fn']), p['attr'], cfg),
          'setattr_profile': lambda p, trees, cfg: _resolve_setattr_profile(
              p['module'], _fn_def(trees.get(p['module']), p['fn']), p['attr'],
              p, cfg)}


def _promise_written(code_dir, promise):
    """Is the promise still written where the claim says it is?

    The other half of the check: a code default that no document promises any more
    should not be certified as "consistent" by a rule that only reads code, because
    then the register keeps verifying a sentence nobody maintains. Every named file
    must exist; at least one of them must carry the flag, a 默认/default word and
    the promised value on the same line.
    """
    for rel in promise['docs']:
        if not isfile(join(code_dir, rel)):
            return f'文档不存在：{rel}', None
    for rel in promise['docs']:
        with open(join(code_dir, rel), encoding='utf-8', errors='replace') as fh:
            for lineno, line in enumerate(fh, 1):
                if promise['flag'] not in line or not DEFAULT_WORD_RE.search(line):
                    continue
                hits = [tok for tok in promise['tokens'] if tok in line]
                # a dict promise (per-domain defaults) requires every number: "默认 50"
                # alone used to read as a promise the archaeal profile kept while the
                # bacterial one said 200
                if len(hits) == len(promise['tokens']) or (
                        hits and not isinstance(promise['promised'], dict)):
                    return None, f'{rel}:{lineno}'
    return f'{promise["flag"]} 的默认句不再与取值 {promise["tokens"]} 写在同一行', None


def default_promise_problems(code_dir):
    """Rule 20: for every documented default, compare promise with delivery."""
    problems = []
    trees = {}
    for module in DEFAULT_RULE_MODULES:
        path = join(code_dir, 'neoHGT', f'{module}.py')
        if isfile(path):
            with open(path, encoding='utf-8', errors='replace') as fh:
                trees[module] = ast.parse(fh.read())
    cfg = None
    cfg_state = 'no PyYAML in this interpreter'
    cfg_path = join(code_dir, 'neoHGT', 'config.yml')
    if isfile(cfg_path):
        try:
            import yaml
        except ImportError:
            yaml = None
        if yaml is not None:
            with open(cfg_path, encoding='utf-8') as fh:
                cfg = yaml.safe_load(fh)
            cfg_state = 'read'
    judged, agreed, disagreed, unjudgeable = 0, [], [], []
    for promise in DEFAULT_PROMISES:
        judged += 1
        # a table default has no establishing function: the parser delivers it
        fn = (None if promise['shape'] == 'table'
              else _fn_def(trees.get(promise['module']), promise['fn']))
        if fn is None and promise['shape'] != 'table':
            problems.append(
                f"默认申明 {promise['key']} 指着 {promise['module']}.py 的 "
                f"{promise['fn']}()，那个函数不在了：交付值无从判定")
            unjudgeable.append(promise['key'])
            continue
        delivered, evidence = SHAPES[promise['shape']](promise, trees, cfg)
        why, written_at = _promise_written(code_dir, promise)
        if written_at is None:
            problems.append(
                f"默认申明 {promise['key']} 指着 {promise['docs']}，"
                f"但{why or '那里找不到这句'}：申明与写下的句子脱钩了")
        if delivered is None:
            unjudgeable.append(f"{promise['key']}（{evidence or '形态没认出来'}）")
            continue
        want = promise['promised']
        if isinstance(want, dict):
            ok = isinstance(delivered, dict) and all(
                delivered.get(k) == _norm_default(v) for k, v in want.items())
            got_repr = delivered
        else:
            ok = delivered == _norm_default(want)
            got_repr = delivered
        if ok:
            agreed.append(f"{promise['key']}={got_repr}@{evidence}")
        else:
            disagreed.append(promise['key'])
            problems.append(
                f"文档说 {promise['flag']} 默认 {want!r}，代码交付的是 {got_repr!r}"
                f"（依据 {evidence}）：承诺与实际交付分家了")
    if judged and not agreed and not disagreed:
        problems.append(f'{judged} 条默认申明一条都没判定成：'
                        '全部落在"看不见"，绿灯不算核对过')
    return ({'promises': judged, 'agree': agreed, 'disagree': disagreed,
             'unjudgeable': unjudgeable, 'config_state': cfg_state,
             'modules_read': sorted(trees)},
            problems)



def route_equivalence_problems(code_dir, ledger_path):
    """Rule 19: a divergent read needs a converter or a stated reason, in writing."""
    problems = []
    if not isfile(ledger_path):
        return ({'state': 'ledger missing'},
                [f'路由等价性清单不存在（{ledger_path}）：第 19 条规则没有可比对的账本，'
                 '这不是"没有问题"，这是没有核对'])
    try:
        with open(ledger_path, encoding='utf-8') as fh:
            ledger = json.load(fh)
    except (OSError, ValueError) as exc:
        return ({'state': 'ledger unreadable'},
                [f'路由等价性清单读不出来（{exc}）：核对退化为猜测'])
    keys = ledger.get('keys') or {}
    if not keys:
        return ({'state': 'ledger empty'},
                ['路由等价性清单里一个键都没有：账本空着不等于代码干净'])
    divergent = sorted(k for k, v in keys.items()
                       if v.get('divergent') and not v.get('converter'))
    explained = sorted(set(divergent) & set(ROUTE_NEUTRALIZERS))
    missing = sorted(set(divergent) - set(ROUTE_NEUTRALIZERS))
    for entry in missing:
        v = keys[entry]
        problems.append(
            f'配置键 {entry} 经配置读入时是 {v.get("yaml_type")}，而它的 CLI 选项声明 '
            f'type={v.get("cli_type")}：两条路交出不同的对象，既没有转换器也没有写明由谁归一化')
    stale = sorted(set(ROUTE_NEUTRALIZERS) - set(divergent))
    for entry in stale:
        problems.append(f'路由豁免清单里的 {entry} 已不再是分歧键（键改名、值类型对齐或已加'
                        '转换器）：豁免不能比它辩护的问题活得久')
    # the ledger's own population is re-derived from the code: a key the code now
    # reads but the ledger never heard of means someone added a knob and skipped
    # the audit.
    rows, unresolved = route_reads(code_dir)
    seen = {e for (_m, _a, e, _c, _b) in rows}
    unlisted = sorted(seen - set(keys))
    for entry in unlisted[:10]:
        problems.append(f'代码里新出现经配置读取的键 {entry} 不在路由账本里：'
                        '先跑 --emit-route 重新发账本，再判定它是否有后果')
    report = {'ledger_keys': len(keys), 'reads': len(rows),
              'divergent': len(divergent),
              'divergent_explained': explained,
              'divergent_undeclared': missing,
              'ledger_stale_exemptions': stale,
              'ledger_unlisted_keys': unlisted,
              'ledger_unresolved': ledger.get('unresolved_call_sites') or []}
    if len(rows) < 40:
        problems.append(f'路由核对只解析到 {len(rows)} 处读取：调用形态没认全，'
                        '绿灯不算数')
    claim_report, claim_problems = route_claim_problems(code_dir, set(explained))
    problems += claim_problems
    report['route_claims'] = claim_report
    report['route_claim_problems'] = claim_problems
    return report, problems


# ---------------------------------------------------------------------------
# Rule 21: a value that can arrive as None must never be used as a number.
# ---------------------------------------------------------------------------
#
# Three shipped defects had exactly this shape, and all three were found by
# executing rather than reading: `--bandwidth` reached `isinstance(bw, float)`
# (round 45), `--fetch-queries`/`--fetch-retries` reached `cids[i:i+cq]` and
# `range()` (round 46), `--donor-label` reached a comparison against 'prefixed'
# while the docs promised it (round 50). This rule watches the shape instead of
# the symptom.
#
# It deliberately does not try to classify every option's delivery: that classifier
# has been rewritten six times this round and disagreed with itself each time
# (临时文件/r51_delivery_taxonomy.json still reports 71 vs r87's 63 "delivers None"),
# and 20 of config.yml's 101 leaves are keys present-but-empty, so "the config
# supplies it" is not a stable predicate. What is stable is the negative direction:
# read a name in a numeric position, and demand a guard that names it. An
# unguarded case is a problem; the counts are printed so an empty examination
# cannot read as a pass.
NUMERIC_PARENTS = (ast.BinOp, ast.Slice, ast.Compare)
# Below this many options the "nothing was judgeable" guards stay silent: the
# shipped CLI has 156, and a one-option synthetic tree legitimately has nothing
# else to look at -- firing there would make the guard indistinguishable from
# noise, which is how guards get switched off. Same reasoning as rule 19's
# "fewer than 40 reads means the call shapes were not recognised".
VACUITY_MIN_OPTIONS = 20


def _option_index(code_dir):
    """(stage, flag, dest) plus the table default, for every CLI surface."""
    out = []
    for stage in CLI_STAGES:
        path = join(code_dir, 'neoHGT', f'{stage}.py')
        if not isfile(path):
            continue
        with open(path, encoding='utf-8', errors='replace') as fh:
            tree = ast.parse(fh.read())
        seen = set()
        for lst in option_tables(tree):
            for item in lst.elts:
                if (not isinstance(item, (ast.List, ast.Tuple)) or len(item.elts) < 2
                        or not isinstance(item.elts[0], ast.Constant)):
                    continue
                flags = re.findall(r'--?[a-z][\w-]*', str(item.elts[0].value))
                if not flags:
                    continue
                flag = max(flags, key=len)
                if flag in seen:
                    continue
                seen.add(flag)
                dest = flag.lstrip('-').replace('-', '_')
                default = 'ABSENT'
                if len(item.elts) > 2 and isinstance(item.elts[2], ast.Dict):
                    for k, val in zip(item.elts[2].keys, item.elts[2].values):
                        if getattr(k, 'value', None) != 'default':
                            continue
                        try:
                            default = ast.literal_eval(val)
                        except Exception:
                            default = val
                out.append({'stage': stage, 'flag': flag, 'dest': dest,
                            'default': default, 'tree': tree})
    return out


def _numeric_index(tree, module):
    """Every numeric-position read of `self.<attr>`/`args.<attr>` in one module.

    Built once per module and filtered per option afterwards: the first version
    walked all nine trees for every one of the 156 options, which pushed the
    register gate from ~25 s to ~115 s -- a gate that slow gets skipped, and a
    skipped gate checks nothing.

    Each record carries the attribute names mentioned by the tests it sits under
    (`guard_names`), so "is this read guarded for *these* names" is a set
    intersection rather than a re-traversal.
    """
    parent, func = {}, {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            parent[id(child)] = node
        if isinstance(node, ast.FunctionDef):
            start, end = node.lineno, node.end_lineno or node.lineno
            for child in ast.walk(node):
                func.setdefault(id(child), (node.name, start, end))
    out = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Attribute)
                and isinstance(node.value, ast.Name)
                and node.value.id in ('self', 'args')):
            continue
        p = parent.get(id(node))
        if isinstance(p, ast.Assign) and any(t is node for t in p.targets):
            continue
        kind = None
        if isinstance(p, ast.BinOp):
            # `'... %s' % self.x` and `'text: ' + str(self.x)` are messages: the
            # `%`/`+` is formatting, not arithmetic, and reading them as numeric
            # uses is how the first three reports of this rule were all false
            # alarms on `--taxdump`.
            others = [n for n in ast.walk(p) if n is not node]
            texty = any(isinstance(o, ast.Constant) and isinstance(o.value, str)
                        for o in others)
            kind = ('message' if texty and isinstance(p.op, (ast.Mod, ast.Add))
                    else 'arithmetic')
        elif isinstance(p, ast.Slice) or (
                isinstance(p, ast.Call) and getattr(p.func, 'id', '') == 'range'):
            # `len(x)` is deliberately not a bound: the attribute it counts is
            # usually a container loaded from disk, not the option's value.
            kind = 'bound'
        elif (isinstance(p, ast.Compare)
                and any(isinstance(c, ast.Constant)
                        and isinstance(c.value, (int, float))
                        and not isinstance(c.value, bool) for c in p.comparators)):
            kind = 'comparison'
        if kind is None or kind == 'message':
            continue
        fn = func.get(id(node), ('<module>', 0, 10 ** 9))
        guard_names, walker = set(), node
        while walker is not None and getattr(walker, 'lineno', 0) >= fn[1]:
            test = None
            if isinstance(walker, (ast.If, ast.IfExp)):
                test = walker.test
            elif isinstance(walker, ast.BoolOp):
                test = walker
            if test is not None:
                guard_names.update(re.findall(r'\b(?:self|args)\.(\w+)',
                                              ast.unparse(test)))
            walker = parent.get(id(walker))
        out.append({'module': module, 'line': node.lineno, 'kind': kind,
                    'function': fn[0], 'name': node.attr,
                    'guard_names': guard_names})
    return out


def _drop_reassigned(assign_lists, names, dest):
    """Keep only the names that carry the option's value, not one it is copied into.

    `self.taxdump` starts as `args.taxdump` and later holds a taxonomy dict loaded
    from disk; counting `len(self.taxdump)` as "the option reached a numeric
    position" is a false alarm of the same kind as the six r88 reported. A name is
    dropped as soon as some assignment gives it a value that did not come from this
    option -- and the number dropped is reported, because silently shrinking the
    population is how a rule starts agreeing with itself.
    """
    keep, dropped = set(), set()
    pat = rf"\bargs\.{re.escape(dest)}\b|['\"]{re.escape(dest)}['\"]"
    for name in names:
        touched, clean = False, True
        for rows in assign_lists.values():
            for attr, rhs in rows:
                if attr != name:
                    continue
                touched = True
                if not re.search(pat, rhs):
                    clean = False
        (keep if (clean or not touched) else dropped).add(name)
    # no fallback to `names`: an option whose attribute is overwritten before it is
    # read has nothing left to blame the option for, and keeping it anyway is how
    # `--taxdump` produced six findings that were all about a loaded dict
    return keep, sorted(dropped)


def numeric_none_problems(code_dir):
    """Rule 21: guard every numeric read of a value that can reach it as None."""
    problems = []
    options = _option_index(code_dir)
    trees = {}
    for stage in CLI_STAGES:
        path = join(code_dir, 'neoHGT', f'{stage}.py')
        if isfile(path):
            with open(path, encoding='utf-8', errors='replace') as fh:
                trees[stage] = ast.parse(fh.read())
    for module in ('util', 'gtdb'):
        path = join(code_dir, 'neoHGT', f'{module}.py')
        if isfile(path):
            with open(path, encoding='utf-8', errors='replace') as fh:
                trees[module] = ast.parse(fh.read())

    rows_all, _unresolved = route_reads(code_dir)
    entries_for_attr = {}
    for (_mod, attr, entry, _c, _b) in rows_all:
        entries_for_attr.setdefault(attr, set()).add(entry)
    cfg = None
    cfg_state = 'no PyYAML'
    cfg_path = join(code_dir, 'neoHGT', 'config.yml')
    if isfile(cfg_path):
        try:
            import yaml
        except ImportError:
            yaml = None
        if yaml is not None:
            with open(cfg_path, encoding='utf-8') as fh:
                cfg = yaml.safe_load(fh)
            cfg_state = 'read'

    can_be_none, examined, unguarded = 0, 0, []
    # Both indices are built once per module and then filtered per option; walking
    # nine syntax trees for each of the 156 options cost the gate 90 extra seconds,
    # and a slow gate is a gate that gets skipped.
    assign_lists = {
        stage: [(t.attr, ast.unparse(node.value))
                for node in ast.walk(tree) if isinstance(node, ast.Assign)
                for t in node.targets
                if isinstance(t, ast.Attribute)
                and isinstance(t.value, ast.Name) and t.value.id == 'self']
        for stage, tree in trees.items()}
    index = {stage: _numeric_index(tree, stage) for stage, tree in trees.items()}

    def storage_for(dest):
        names = {dest}
        pat = rf"\bargs\.{re.escape(dest)}\b|['\"]{re.escape(dest)}['\"]"
        for rows in assign_lists.values():
            for attr, rhs in rows:
                if re.search(pat, rhs):
                    names.add(attr)
        return names

    def numeric_sites(names):
        for stage, recs in index.items():
            for rec in recs:
                if rec['name'] in names:
                    yield dict(rec, guarded=bool(rec['guard_names'] & names))

    # A detector-level control, independent of how many options happen to be
    # None-capable today: count the numeric reads of *every* CLI-derived attribute.
    # If that is zero, the walk itself is blind -- which is the failure mode a
    # per-subset count cannot see, because after rounds 45/46/50 were fixed the
    # None-capable subset really may contain no numeric use at all.
    all_names = set()
    for opt in options:
        all_names |= storage_for(opt['dest'])
    total_sites = len({(s['module'], s['line'], s['kind'], s['name'])
                       for s in numeric_sites(all_names)})
    seen, dropped_names, reassigned_out = set(), set(), 0
    for opt in options:
        names = storage_for(opt['dest'])
        names, dropped = _drop_reassigned(assign_lists, names, opt['dest'])
        dropped_names.update(dropped)
        if not names:
            # every attribute this option feeds was overwritten from elsewhere
            # first, so the option is not where a numeric use gets its value
            reassigned_out += 1
            continue
        supplied = False
        if opt['default'] not in ('ABSENT', None):
            supplied = True
        entries = set(entries_for_attr.get(opt['dest'], ()))
        for n in names:
            entries.update(entries_for_attr.get(n, ()))
        if cfg is not None and not supplied:
            for entry in entries:
                cur = cfg
                for part in entry.split('.'):
                    if not isinstance(cur, dict) or part not in cur:
                        cur = None
                        break
                    cur = cur[part]
                if cur is not None:
                    supplied = True
                    break
        if supplied:
            continue
        can_be_none += 1
        for site in numeric_sites(names):
            key = (site['module'], site['line'], site['kind'], site['name'])
            if key in seen:
                # the same flag is declared on several CLI surfaces; reporting
                # one line four times reads as four findings
                continue
            seen.add(key)
            examined += 1
            if not site['guarded']:
                unguarded.append(site)
                problems.append(
                    f"{site['module']}.py:{site['line']} 在{site['kind']}位置用了 "
                    f"{opt['flag']} 的值（{site['function']}()），但这个选项可能以 None "
                    f"到达这里，而同一函数里没有以该属性名设的守卫：None 进数值位置"
                    f"就是 TypeError 或静默改变行为")
    if len(options) >= VACUITY_MIN_OPTIONS and can_be_none == 0:
        problems.append(f'{len(options)} 个选项里没有一条被判为"可能是 None"：'
                        '判定形态没认全，绿灯不算数')
    if (len(options) >= VACUITY_MIN_OPTIONS and can_be_none
            and examined == 0 and total_sites == 0):
        problems.append('全部选项属性里一个数值使用都找不到：这条规则在看空气，'
                        f'{can_be_none} 个可能交付 None 的选项没有被看过')
    return ({'options': len(options), 'can_be_none': can_be_none,
             'numeric_reads_examined': examined,
             'numeric_reads_in_all_cli_attrs': total_sites,
             'options_reassigned_before_any_use': reassigned_out,
             'names_dropped_as_reassigned': sorted(dropped_names)[:20],
             'names_dropped_count': len(dropped_names),
             'unguarded': [{'where': f"{s['module']}.py:{s['line']}",
                            'kind': s['kind'], 'function': s['function'],
                            'name': s['name']} for s in unguarded],
             'config_state': cfg_state}, problems)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', default=None,
                    help='the adaptation plan; when it is not available (e.g. '
                         'a repository checkout), the committed --plan-ids '
                         'digest is used instead')
    ap.add_argument('--plan-ids', default=None,
                    help='frozen id digest (default: <code>/docs/'
                         '方案编号清单.json), used when --plan is absent')
    ap.add_argument('--emit-ids', default=None,
                    help='write the plan id digest to this path and exit')
    ap.add_argument('--emit-route', default=None,
                    help='write the route-equivalence ledger (needs PyYAML, so '
                         'this runs in the round gates, not in a bare checkout) '
                         'and exit')
    ap.add_argument('--route-ledger', default=None,
                    help='committed route ledger (default: <code>/docs/'
                         '路由等价性清单.json); checking it needs no YAML parser')
    ap.add_argument('--code', default=os.environ.get('NEOHGT_CODE'))
    ap.add_argument('--upstream', default=None,
                    help='copy of the upstream project, used only to work out '
                         'which config.yml keys this adaptation introduced; when '
                         'it is absent the config-key rule judges the smaller '
                         'name-matching population and says so in its summary')
    ap.add_argument('--cache', default=os.environ.get('NEOHGT_CACHE'))
    ap.add_argument('--out', default=None,
                    help='JSON report path; required unless --emit-ids is used')
    ap.add_argument('--ledger', default=None,
                    help='plan §6 task ledger (default: <code>/docs/'
                         '计划任务对照表.md)')
    ap.add_argument('--defects', default=None,
                    help='appendix F defect register (default: <code>/docs/'
                         '既有缺陷登记.md)')
    ap.add_argument('--deliverables', default=None,
                    help='§12 deliverable register (default: <code>/docs/'
                         '交付物清单.md)')
    ap.add_argument('--risks', default=None,
                    help='§9 risk register (default: <code>/docs/'
                         '风险应对登记.md)')
    ap.add_argument('--acceptance', default=None,
                    help='§1.4/§10.7 DoD register (default: <code>/docs/'
                         '验收标准对照.md)')
    ap.add_argument('--decisions', default=None,
                    help='§8 ADR register (default: <code>/docs/'
                         '决策落点对照.md)')
    ap.add_argument('--tests', default=None,
                    help='§10.2/10.3/10.5 verification register (default: '
                         '<code>/docs/测试用例对照.md)')
    ap.add_argument('--assumptions', default=None,
                    help='§2.3 search-layer assumption register (default: '
                         '<code>/docs/搜索层假设核对.md)')
    ap.add_argument('--revisions', default=None,
                    help='§0.5 review-revision register (default: '
                         '<code>/docs/方案修订记录核对.md)')
    ap.add_argument('--milestones', default=None,
                    help='§13 milestone register (default: <code>/docs/'
                         '里程碑核对.md)')
    ap.add_argument('--discovery', default=None,
                    help='the running self-discovery list (default: <code>/'
                         'docs/GTDB-数据库.md)')
    args = ap.parse_args()
    if not args.upstream:
        # the layout this project has always been run in: two sibling folders
        # under one workspace root, the adapted one and the untouched upstream.
        # Absspath so `--code .` (how CI invokes this) resolves the same way.
        args.upstream = os.path.abspath(join(args.code, os.pardir, 'neoHGT-main'))
    if args.out:
        args.out = guarded_out(args.out, args.code)
    if args.emit_ids:
        args.emit_ids = guarded_out(args.emit_ids, args.code)
    if args.emit_route:
        args.emit_route = guarded_out(args.emit_route, args.code)
        n = route_ledger_emit(args.code, args.emit_route)
        print(f'{n} config keys read by the stage modules -> {args.emit_route}')
        return 0
    if not args.out and not args.emit_ids:
        ap.error('--out is required unless --emit-ids is used')

    digest = None
    if args.plan:
        if not isfile(args.plan):
            raise SystemExit(f'plan not found: {args.plan}')
        with open(args.plan, encoding='utf-8') as f:
            plan_text = f.read()
        inv = plan_inventory(plan_text)
        id_source = 'plan'
        if args.emit_ids:
            write_inventory(args.emit_ids, plan_text, args.plan)
            print(f'wrote plan id digest to {args.emit_ids}')
            return 0
    else:
        path = args.plan_ids or join(args.code, 'docs', '方案编号清单.json')
        if not isfile(path):
            raise SystemExit(
                f'neither --plan nor its digest is available ({path}); run '
                f'"--plan <plan> --emit-ids {path}" once to create it')
        inv, digest = load_inventory(path)
        plan_text, id_source = '', f'digest:{basename(path)}'
        if args.emit_ids:
            raise SystemExit('--emit-ids needs --plan')
    ids = inv['test_ids']
    cites = cite_map(join(args.code, 'neoHGT', 'tests'))

    rows, gaps = [], []
    # A repository checkout has the code, the registers and the digest, but not
    # the workspace evidence tree, so artifact-only ids cannot be verified here.
    # They are reported as unverifiable rather than silently passed or counted as
    # gaps, which would make the CI step fail for the wrong reason.
    cache_present = bool(args.cache) and isdir(args.cache)
    unverifiable = []
    concrete = [i for i in ids if not i.endswith(GROUP_SUFFIX)]
    for i in ids:
        if i.endswith(GROUP_SUFFIX):
            fam = i.split('.')[0]
            members = [m for m in concrete if m.startswith(fam + '.')]
            traced = bool(cites.get(i)) or any(
                m in cites or EVIDENCE.get(m) for m in members)
            rows.append({'id': i, 'kind': 'test' if cites.get(i) else 'group',
                         'members': members, 'tests': cites.get(i, []),
                         'traced': traced})
        else:
            test = cites.get(i)
            artifact = EVIDENCE.get(i)
            artifact_path = (join(args.cache, artifact)
                             if artifact and args.cache
                             else None)
            ok_artifact = bool(artifact_path and exists(artifact_path))
            # an id evidenced only by a workspace artifact is not "untraceable"
            # here, it is simply not checkable without that artifact
            here = bool(test) or ok_artifact or (bool(artifact)
                                                 and not cache_present)
            rows.append({'id': i, 'kind': 'test' if test else 'artifact',
                         'tests': test or [],
                         'artifact': artifact, 'artifact_exists': ok_artifact,
                         'traced': here,
                         'verifiable': bool(test) or ok_artifact})
        if not rows[-1]['traced']:
            gaps.append(i)
        elif not rows[-1].get('verifiable', True):
            unverifiable.append(i)

    report = {'plan': args.plan, 'id_source': id_source,
              'ids': len(ids), 'rows': rows, 'gaps': gaps,
              'ids_unverifiable_in_this_checkout': unverifiable}

    # --- the §7 option list: is every planned flag real, and does the code
    # actually read it? Only the §7 tables are the CLI contract; other sections
    # mention third-party flags (diamond --taxonmap, makeblastdb -taxid_map).
    plan_flags = inv['section7_flags']
    flags = cli_flags(args.code)
    real = set().union(*flags.values())
    missing_flags = [x for x in plan_flags
                     if x not in real and x.strip('-') not in
                     ('help', 'db', 'in', 'out', 'dir', 'action')]
    dead = dead_gtdb_options(args.code, flags)
    untested_flags = flags_without_tests(plan_flags, args.code)
    report['plan_flags'] = len(plan_flags)
    report['flags_missing_from_cli'] = missing_flags
    report['options_parsed_but_never_read'] = dead
    report['flags_without_tests'] = untested_flags

    # --- the two registers the plan asks for: the §6 numbered tasks (is every
    # planned task really implemented, or explicitly not?) and the appendix F
    # defect list (is every shared-root-cause defect fixed, or registered with
    # a reason and a reproducer?).
    root_dir = dirname(abspath(args.code))
    # a CI checkout has no workspace evidence tree; ledger citations of
    # `缓存测试/`... artifacts are then skipped (and counted) instead of failing
    tolerate_ws = not isdir(join(root_dir, WORKSPACE_DIRS['缓存测试']))
    ws_skips = []
    ids6, idsF = inv['tasks'], inv['defects']
    idsD, idsR = inv['deliverables'], inv['risks']
    idsDod, idsAdr = inv['dods'], inv['decisions']
    ledgers = [
        check_ledger('对照表', ids6,
                     args.ledger or join(args.code, 'docs',
                                         '计划任务对照表.md'),
                     args.code, root_dir, LEDGER_STATUSES,
                     tolerate_missing_workspace=tolerate_ws, ws_skips=ws_skips),
        check_ledger('登记册', idsF,
                     args.defects or join(args.code, 'docs',
                                          '既有缺陷登记.md'),
                     args.code, root_dir, DEFECT_STATUSES,
                     tolerate_missing_workspace=tolerate_ws, ws_skips=ws_skips),
        check_ledger('交付物表', idsD,
                     args.deliverables or join(args.code, 'docs',
                                               '交付物清单.md'),
                     args.code, root_dir, DELIVER_STATUSES,
                     tolerate_missing_workspace=tolerate_ws, ws_skips=ws_skips),
        check_ledger('风险表', idsR,
                     args.risks or join(args.code, 'docs', '风险应对登记.md'),
                     args.code, root_dir, RISK_STATUSES,
                     tolerate_missing_workspace=tolerate_ws, ws_skips=ws_skips),
        check_ledger('验收表', idsDod,
                     args.acceptance or join(args.code, 'docs',
                                             '验收标准对照.md'),
                     args.code, root_dir, DOD_STATUSES,
                     tolerate_missing_workspace=tolerate_ws, ws_skips=ws_skips),
        check_ledger('决策表', idsAdr,
                     args.decisions or join(args.code, 'docs',
                                            '决策落点对照.md'),
                     args.code, root_dir, DECISION_STATUSES,
                     tolerate_missing_workspace=tolerate_ws, ws_skips=ws_skips),
        check_ledger('用例表', inv['verification'],
                     args.tests or join(args.code, 'docs', '测试用例对照.md'),
                     args.code, root_dir, VERIF_STATUSES,
                     tolerate_missing_workspace=tolerate_ws, ws_skips=ws_skips),
        check_ledger('假设表', inv['assumptions'],
                     args.assumptions or join(args.code, 'docs',
                                              '搜索层假设核对.md'),
                     args.code, root_dir, ASSUMPTION_STATUSES,
                     tolerate_missing_workspace=tolerate_ws, ws_skips=ws_skips),
        check_ledger('修订表', inv['revisions'],
                     args.revisions or join(args.code, 'docs',
                                            '方案修订记录核对.md'),
                     args.code, root_dir, REVISION_STATUSES,
                     tolerate_missing_workspace=tolerate_ws, ws_skips=ws_skips),
        check_ledger('里程碑表', inv['milestones'],
                     args.milestones or join(args.code, 'docs', '里程碑核对.md'),
                     args.code, root_dir, MILESTONE_STATUSES,
                     tolerate_missing_workspace=tolerate_ws, ws_skips=ws_skips),
    ]
    task_problems = sum((l['problems'] for l in ledgers), [])
    for l, key in ((ledgers[0], 'tasks'), (ledgers[1], 'defects'),
                   (ledgers[2], 'deliverables'), (ledgers[3], 'risks'),
                   (ledgers[4], 'dods'),
                   (ledgers[5], 'decisions'),
                   (ledgers[6], 'verification'),
                   (ledgers[7], 'assumptions'),
                   (ledgers[8], 'revisions'),
                   (ledgers[9], 'milestones')):
        report[f'{key}_ledger'] = l['path']
        report[f'{key}_ids'] = l['ids']
        report[f'{key}_rows'] = l['rows']
        report[f'{key}_unlisted'] = l['unlisted']
        report[f'{key}_not_in_plan'] = l['extra']
        report[f'{key}_statuses'] = l['counts']
    # appendix C/D speak the *target* command line: every --flag they use must
    # be real, belong to a third-party tool, or be registered as a deliberate
    # renaming in the deliverables document.
    doc = args.deliverables or join(args.code, 'docs', '交付物清单.md')
    declared = deviation_flags(doc, args.code)
    unexplained = [x for x in inv['appendix_flags']
                   if x not in real and x not in THIRD_PARTY_FLAGS
                   and x not in declared]
    task_problems += [f'附录 C/D 用了 {x}，CLI 里没有、也没在偏差表里登记'
                      for x in unexplained]
    report['appendix_flags'] = len(inv['appendix_flags'])
    report['appendix_flags_unexplained'] = unexplained
    report['appendix_flags_declared_deviation'] = sorted(declared)
    # a planned option nobody hands to the parser is an option nobody verified:
    # wrong dest, wrong type, or parsed-then-ignored all hide here
    task_problems += [f'{x} 是方案 §7 承诺的选项，但没有任何用例把它交给解析器'
                      for x in untested_flags if x not in declared]
    report['flags_untested_problems'] = [x for x in untested_flags
                                         if x not in declared]
    # a tested flag is not a tested *value*: --mimag was handed to the parser by
    # a test while 'any', one of its four declared choices, reached no branch
    opts_choice = cli_choice_options(args.code)
    declared_values = declared_choice_values(doc, args.code)
    untested_values, backlog_values = choice_values_without_tests(
        args.code, declared_values)
    task_problems += [f'{x} 是 CLI 声明的取值，但没有任何用例把它交给该选项'
                      for x in untested_values]
    report['cli_choice_options'] = len(opts_choice)
    report['cli_choice_values'] = sum(len(e['choices']) for e in opts_choice)
    report['choice_values_without_tests'] = untested_values
    report['choice_values_backlog'] = backlog_values
    report['choice_values_declared_deviation'] = sorted(declared_values)
    # the documents we ship are copy-pasted: an example that names an option the
    # very sub-command in that example does not accept is a broken command line
    subs = cli_subcommands(args.code)
    used_doc = doc_command_flags(args.code, subs)
    declared_doc = declared_doc_usages(doc, subs, flags)
    doc_problems = doc_flag_problems(used_doc, flags, declared_doc)
    task_problems += doc_problems
    report['doc_subcommands'] = list(subs)
    report['doc_command_flags'] = len(used_doc)
    report['doc_command_usages'] = sum(len(v) for v in used_doc.values())
    report['doc_command_problems'] = doc_problems
    report['doc_command_declared_deviation'] = sorted(
        f'{w} {flag}' for w, stage, flag in declared_doc)
    symsB = inv['appendix_b_symbols']
    problemsB = appendix_b_problems(symsB, args.code, doc)
    task_problems += problemsB
    report['appendix_b_symbols'] = len(symsB)
    report['appendix_b_problems'] = problemsB
    # a running list whose heading counts the items by hand drifts the moment an
    # item is appended; that is the one register every round edits without looking
    disc = self_discovery_check(
        args.discovery or join(args.code, 'docs', 'GTDB-数据库.md'))
    task_problems += disc['problems']
    report['self_discovery'] = disc
    # a config key that nothing reads is a knob that only looks like a setting
    cfg_report, cfg_problems = config_key_problems(args.code, doc, args.upstream)
    task_problems += cfg_problems
    # ... and a config *section* chosen by an option value needs a written,
    # checked disposition before it can pass for a setting
    dyn_report, dyn_problems = dynamic_section_problems(args.code, doc)
    task_problems += dyn_problems
    report['dynamic_config_sections'] = dyn_report
    report['dynamic_config_problems'] = dyn_problems
    # ... and the sections the gate itself reads must still be there to be read
    head_report, head_problems = ledger_heading_problems(args.code, doc)
    task_problems += head_problems
    report['ledger_headings'] = head_report
    report['ledger_heading_problems'] = head_problems
    reuse_report, reuse_problems = reuse_inert_problems(args.code)
    task_problems += reuse_problems
    report['reuse_inert'] = reuse_report
    report['reuse_inert_problems'] = reuse_problems
    docval_report, docval_problems = doc_value_problems(args.code, subs,
                                               args.plan)
    task_problems += docval_problems
    report['doc_command_values'] = docval_report
    report['doc_value_problems'] = docval_problems
    route_report, route_problems = route_equivalence_problems(
        args.code, args.route_ledger or join(args.code, 'docs',
                                             '路由等价性清单.json'))
    task_problems += route_problems
    report['route_equivalence'] = route_report
    report['route_equivalence_problems'] = route_problems
    defaults_report, defaults_problems = default_promise_problems(args.code)
    task_problems += defaults_problems
    report['default_promises'] = defaults_report
    report['default_promise_problems'] = defaults_problems
    numeric_report, numeric_problems = numeric_none_problems(args.code)
    task_problems += numeric_problems
    report['numeric_none'] = numeric_report
    report['numeric_none_problems'] = numeric_problems
    report['config_keys'] = cfg_report
    checkers = checker_tools(args.code)
    checker_problems, checker_declared = checker_tool_problems(
        checkers, doc)
    task_problems += checker_problems
    report['checker_tools'] = {k: v for k, v in checkers.items()}
    report['checker_tools_problems'] = checker_problems
    report['checker_tools_declared'] = checker_declared
    # a test that asserts nothing is a green light with no information in it
    vacuous = assertion_free_tests(args.code)
    declared_vac = declared_vacuous_tests(doc)
    undeclared_vac = [x for x in vacuous if x not in declared_vac]
    task_problems += [
        f'{x} 这个用例里没有任何可观察的断言，代码什么都不做它也会通过'
        for x in undeclared_vac]
    report['assertion_free_tests'] = vacuous
    report['assertion_free_undeclared'] = undeclared_vac
    report['assertion_free_declared_legacy'] = sorted(declared_vac - set(vacuous))
    report['ledger_rows'] = ledgers[0]['rows']
    report['plan_tasks'] = ledgers[0]['ids']
    report['tasks_unlisted'] = ledgers[0]['unlisted']
    report['tasks_not_in_plan'] = ledgers[0]['extra']
    counts = ledgers[0]['counts']
    report['tasks_done'] = counts.get('已实现', 0)
    report['tasks_partial'] = counts.get('部分实现', 0)
    report['tasks_not_done'] = counts.get('未做', 0)
    report['verification_ids'] = ledgers[6]['ids']
    report['verification_rows'] = ledgers[6]['rows']
    report['verification_statuses'] = ledgers[6]['counts']
    report['digest_source'] = (digest or {}).get('source', '')
    report['digest_generated'] = (digest or {}).get('generated_utc', '')
    report['task_problems'] = task_problems
    report['verdict'] = ('PASS' if not gaps and not missing_flags and not dead
                         and not task_problems
                         else 'FAIL')
    out = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(report, f, indent=1, ensure_ascii=False)
    with open(os.path.splitext(out)[0] + '.md', 'w') as f:
        f.write('# 计划验证条目 ↔ 证据可追溯表\n\n'
                f'plan: `{args.plan}`  \nids: {len(ids)}  |  '
                f'gaps: {len(gaps)}  |  verdict: **{report["verdict"]}**\n\n'
                '| id | 依据 | 说明 |\n| --- | --- | --- |\n')
        for r in rows:
            if r['kind'] == 'test':
                why = ', '.join(r['tests'])
            elif r['kind'] == 'artifact':
                why = (f'`{r["artifact"]}` '
                       f'{"存在" if r["artifact_exists"] else "缺失"}')
            else:
                why = ', '.join(m for m in r['members'] if m in cites) or \
                    '整族由具体成员覆盖'
            f.write(f'| {r["id"]} | {"测试" if r["kind"] == "test" else "产物"}'
                    f' | {why} |\n')
        f.write('\n## §7 选项清单核对\n\n'
                f'plan §7 flags: {len(plan_flags)} | '
                f'CLI 中缺失: {missing_flags or "无"} | '
                f'解析后从未被读取: {dead or "无"}\n')
        for l in ledgers:
            caption = {'对照表': '§6 编号任务', '登记册': '附录 F 既有缺陷',
                       '交付物表': '§12 交付物',
                       '风险表': '§9 风险缓解',
                       '验收表': '§1.4/§10.7 验收标准',
                       '决策表': '§8 设计决策'}.get(l['label'], l['label'])
            f.write(f'\n## {l["label"]}核对（{caption}）\n\n'
                    f'登记文件: `{l["path"]}`  \n'
                    f'方案编号: {l["ids"]} | 表内行: {l["rows"]} | 状态分布: '
                    + ', '.join(f'{k} {v}' for k, v in
                               sorted(l['counts'].items())) + '\n\n'
                    f'未在表中出现的编号: {l["unlisted"] or "无"}  \n'
                    f'表里凭空多出的编号: {l["extra"] or "无"}  \n'
                    f'问题: {l["problems"] or "无"}\n')
        f.write('\n## 附录 C/D 命令形态核对\n\n'
                f'示例命令里的选项: {report["appendix_flags"]} | '
                f'既不在 CLI、也不是第三方工具、又没登记偏差: '
                f'{unexplained or "无"} | 已登记偏差: '
                f'{report["appendix_flags_declared_deviation"] or "无"}\n')
        f.write('\n## 附录 B 符号面核对\n\n'
                f'方案承诺的符号: {report["appendix_b_symbols"]} | '
                f'问题: {problemsB or "无"}\n')
        f.write('\n## §7 选项的用例覆盖\n\n'
                f'§7 选项: {report["plan_flags"]} | 没有任何用例交给解析器的: '
                f'{report["flags_without_tests"] or "无"} | '
                f'其中未登记偏差的: '
                f'{report["flags_untested_problems"] or "无"}\n')
        f.write('\n## 选项取值的用例覆盖\n\n'
                f'带 choices 的选项: {report["cli_choice_options"]} | '
                f'取值总数: {report["cli_choice_values"]} | '
                f'未交给任何用例的取值(判定阶段): '
                f'{report["choice_values_without_tests"] or "无"} | '
                f'未判定阶段 (应为空): '
                f'{len(report["choice_values_backlog"])} | '
                f'申明为偏差的: '
                f'{report["choice_values_declared_deviation"] or "无"}\n')
        f.write('\n## config.yml 的旋钮核对\n\n'
                f'本适配新增的键: {len(cfg_report["declared_own_keys"])} | '
                f'被代码读到的: {len(cfg_report["read_keys"])} | '
                f'没人读的: {cfg_report["never_read"] or "无"} | '
                f'代码读了但文件里没有的: '
                f'{list(cfg_report["read_but_absent"]) or "无"}\n'
                f'读取写法两种：下标 {cfg_report["subscript_paths"]} 条路径、'
                f'入口串 {cfg_report["entry_paths"]} 条（展开不了的 get_config '
                f'调用: {cfg_report["entry_unexpanded"] or "无"}）\n'
                f'按取值选择配置段: {dyn_report["candidates"] or "无"} | '
                f'登记表处置: {dyn_report["declarations"] or "无"} | '
                f'问题: {dyn_problems or "无"}\n'
                f'门禁要读的章节标题: {head_report["sections"]}/'
                f'{len(DECLARED_SECTIONS)} 仍在 —— '
                f'{head_report["read"]} | 问题: {head_problems or "无"}\n')
        f.write('\n## 文档示例里的选项取值\n\n'
                f'带取值的示例: {docval_report["pairs"]} 个 | 字面判定: '
                f'{docval_report["judged"]} 个 | 占位符（不判定，只计数）: '
                f'{len(docval_report["placeholders"])} 个 | 围栏不成对的文件: '
                f'{docval_report["unbalanced_fences"] or "无"} | '
                f'方案文件: {docval_report["plan_state"]}'
                f'{"（判定 " + str(docval_report["plan_judged"]) + " 个取值）" if docval_report["plan_judged"] else ""} | '
                f'问题: {docval_problems or "无"}\n')
        f.write('\n## 复用路径（--from-scores）不生效的选项\n\n'
                f'代码申明的 REUSE_INERT: {reuse_report["declared"]} | '
                f'调用图判定的: {reuse_report["derived"]} | '
                f'额外豁免: {reuse_report["exceptions"]} | '
                f'复用路径可达的函数: {reuse_report["live_functions"]} | '
                f'被跳过的步骤: {reuse_report["blocked_functions"]} | '
                f'问题: {reuse_problems or "无"}\n')
        f.write('\n## 核对工具自身的执行者\n\n'
                f'`tools/verify_*.py` 共 {len(checkers)} 个；没有任何执行者的: '
                f'{[k for k, v in checkers.items() if not v] or "无"} | '
                f'已申明的: {checker_declared or "无"} | '
                f'问题: {checker_problems or "无"}\n')
        f.write('\n## 断言缺失的用例核对\n\n'
                f'没有任何可观察断言的 test_* 函数: {vacuous} \n'
                f'未登记为既有缺陷的 (应为空): {undeclared_vac or "无"} | '
                f'登记表里已不再出现的旧条目: '
                f'{report["assertion_free_declared_legacy"] or "无"}\n')
        f.write('\n## 我们文档里的 neoHGT 命令示例\n\n'
                f'子命令: {", ".join(subs)} | 示例中用到的选项: '
                f'{report["doc_command_flags"]} 个（{report["doc_command_usages"]} '
                f'处）| 该子命令其实没有这个选项: '
                f'{doc_problems or "无"} | 已登记为偏差的: '
                f'{report["doc_command_declared_deviation"] or "无"}\n')
        f.write('\n## 自我发现清单核对\n\n'
                f'标题: `{disc["heading"]}`  \n'
                f'声明轮次区间: 第{disc["rounds"][0]}至{disc["rounds"][1]}轮 | '
                f'标题声明条数: {disc["declared_count"]} | 实际条数: '
                f'{disc["items"]} | 条目标注的轮次: '
                f'{sorted(set(disc["tagged_rounds"]))} \n'
                f'问题: {disc["problems"] or "无"}\n')
    claims = route_report.get('route_claims') or {}
    print(f'{len(ids)} plan ids; {len(gaps)} untraced; '
          f'{len(plan_flags)} §7 flags, {len(missing_flags)} missing from the '
          f'CLI, {len(dead)} parsed-but-unused; '
          f'{ledgers[0]["ids"]} §6 tasks + {ledgers[1]["ids"]} appendix-F '
          f'defects + {ledgers[2]["ids"]} §12 deliverables + '
          f'{ledgers[3]["ids"]} §9 risks + {ledgers[4]["ids"]} DoD criteria '
          f'+ {ledgers[5]["ids"]} ADRs + {ledgers[6]["ids"]} §10.2-10.5 '
          f'verification ids + {ledgers[7]["ids"]} §2.3 assumptions + '
          f'{ledgers[8]["ids"]} §0.5 revisions + '
          f'{ledgers[9]["ids"]} §13 milestones; '
          f'{len(inv["appendix_flags"])} appendix C/D flags '
          f'({len(unexplained)} unexplained), '
          f'{report["appendix_b_symbols"]} appendix-B symbols '
          f'({len(problemsB)} unexplained), '
          f'{disc["items"]} self-discovery entries counted by hand against a '
          f'title saying {disc["declared_count"]}, '
          f'{len(cfg_report["declared_own_keys"])} config keys this adaptation '
          f'introduced [{cfg_report["population"]}] '
          f'({len(cfg_report["never_read"])} unread), '
          f'{cfg_report["entry_paths"]} entry-string reads '
          f'({len(cfg_report["entry_unexpanded"])} unexpandable), '
          f'{len(dyn_report["candidates"])} runtime-picked config sections '
          f'({len(dyn_report["declarations"])} dispositioned in the ledger), '
          f'{head_report["sections"]}/{len(DECLARED_SECTIONS)} ledger sections '
          'the gate reads still present as headings, '
          f'{len(reuse_report["declared"])} options declared inert on the '
          f'--from-scores route ({len(reuse_problems)} disagreements with the '
          f'call graph), '
          f'{docval_report["judged"]}/{docval_report["pairs"]} documented '
          f'option values judged '
          f'({len(docval_problems)} rejected or blind), '
          f'{docval_report["plan_judged"]} plan-example values {docval_report["plan_state"]}, '
          f'{len(checkers)} checker tools '
          f'({sum(1 for v in checkers.values() if not v)} with no runner, '
          f'{len(checker_problems)} of them undeclared), '
          f'{len(vacuous)} assertion-free tests ({len(undeclared_vac)} '
          f'undeclared), '
          f'{report["doc_command_flags"]} doc example options '
          f'({len(doc_problems)} not accepted by the sub-command shown), '
          f'{route_report.get("ledger_keys", 0)} routed config keys '
          f'({route_report.get("divergent", 0)} whose two routes differ: '
          f'{len(route_report.get("divergent_explained", []))} explained, '
          f'{len(route_report.get("divergent_undeclared", []))} not), '
          f'{len(route_report.get("ledger_unlisted_keys", []))} keys missing '
          f'from the ledger, '
          f'{claims.get("claims_checked", 0)} written reasons verified in the '
          f'code ({claims.get("citations_seen", 0)} code citations checked, '
          f'{len(claims.get("citations_unjudged", []))} not judgeable here), '
          f'{defaults_report.get("promises", 0)} documented defaults compared with '
          f'the delivery ({len(defaults_report.get("agree", []))} agree, '
          f'{len(defaults_report.get("disagree", []))} disagree, '
          f'{len(defaults_report.get("unjudgeable", []))} not judgeable '
          f'[{defaults_report.get("config_state")}]), '
          f'{numeric_report.get("can_be_none", 0)}/'
          f'{numeric_report.get("options", 0)} options may deliver None, '
          f'{numeric_report.get("numeric_reads_examined", 0)}/'
          f'{numeric_report.get("numeric_reads_in_all_cli_attrs", 0)} numeric reads '
          f'of CLI attributes examined ({len(numeric_report.get("unguarded", []))} '
          f'unguarded), '
          f'{len(task_problems)} problem(s) -> {report["verdict"]}')
    for g in gaps:
        print('  GAP', g)
    if unverifiable:
        print(f'  NOTE: {len(unverifiable)} id(s) are evidenced only by files '
              f'under the workspace cache, which this checkout does not have: '
              f'{", ".join(unverifiable)}')
    if ws_skips:
        print(f'  NOTE: {len(ws_skips)} register citation(s) of workspace-tree '
              f'artifacts are not verifiable in this checkout (no 缓存测试/ '
              f'cache directory here); they are strictly checked when the '
              f'workspace is present')
    for m in missing_flags:
        print('  MISSING FLAG', m)
    for d in dead:
        print('  DEAD OPTION', d)
    for t in task_problems:
        print('  TASK PROBLEM', t)
    print(f'wrote {out}')
    return 0 if report['verdict'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
