#!/usr/bin/env python3
"""Check the user documentation: links, anchors, language and leftovers.

The checks are deliberately conservative:

* a relative Markdown link must resolve to a file in the repository, and an
  in-repository anchor must match a heading in its target file;
* an English document must not contain CJK text except for controlled proper
  names and Chinese file names that are part of the public interface;
* a shell variable used in a fenced example must be defined in the document
  (or be a well-known placeholder);
* public documents must not refer to internal QC registers, workspace-only
  evidence trees or development round numbers.

Exit code 0 means no problem was found.
"""

import os
import re
import sys
from urllib.parse import unquote
from os.path import abspath, dirname, isdir, isfile, join, relpath

CJK = re.compile(r'[\u3400-\u4dbf\u4e00-\u9fff]')
LINK = re.compile(r'\[([^\]]*)\]\(([^)\s]+)(?:\s+"[^"]*")?\)')
FENCE = re.compile(r'^\s*(```+|~~~+)')
VAR = re.compile(r'\$(?:\{)?([A-Za-z_][A-Za-z0-9_]*)(?:\})?')

# Variables that examples may use without a definition in the same file.
KNOWN_VARS = {
    'GTDB_DIR', 'NEOHGT_ROOT', 'NEOHGT_CODE', 'NEOHGT_CACHE', 'NEOHGT_TMP',
    'NEOHGT_PYTHON', 'NEOHGT_DIAMOND', 'NEOHGT_MAIN', 'DB', 'OUT', 'TMP',
    'HOME', 'PATH', 'PWD', 'CONDA_PREFIX', 'GITHUB_OUTPUT', 'GITHUB_WORKSPACE',
    'TMPDIR', 'MPLCONFIGDIR', 'XDG_CACHE_HOME', 'PIP_CACHE_DIR',
    'NUMEXPR_CACHE_DIR', 'PYTHONPATH', 'BASH_SOURCE', 'FOLDER_NAME',
    'DATABASE_DIRECTORY', 'RESULT', 'SAMPLE', 'INPUT',
}

# English documents may contain these proper names or public Chinese file names.
ALLOWED_CJK = (
    '曾子超', 'Zeng Zichao', '简体中文', 'GTDB-数据库.md', '新增函数清单.json',
    '方案编号清单.json', '路由等价性清单.json',
)

FORBIDDEN = (
    ('internal QC', re.compile(r'internal\s+QC', re.I)),
    ('内部QC', re.compile(r'内部\s*QC')),
    ('self-discovery', re.compile(r'self[- ]discovery', re.I)),
    ('自我发现', re.compile(r'自我发现')),
    ('交付物清单', re.compile(r'交付物清单')),
    ('既有缺陷登记', re.compile(r'既有缺陷登记')),
    ('workspace cache name', re.compile(r'缓存测试')),
    ('workspace temp name', re.compile(r'临时文件')),
    ('development round', re.compile(r'\bround\s+\d+\b', re.I)),
    ('plan section', re.compile(r'plan\s+§|方案\s*§')),
    ('DoD id', re.compile(r'\bDoD[-\s]?\d')),
    ('ADR id', re.compile(r'\bADR[-\s]?\d')),
    ('RK id', re.compile(r'\bRK\d')),
)

DOCS = (
    'README.md', 'README.zh-CN.md', 'CHANGELOG.md', 'CONTRIBUTING.md',
    'docs/README.md', 'docs/GTDB-数据库.md', 'docs/en/GTDB-database.md',
    'docs/cli-reference.md', 'docs/config-reference.md',
    'docs/troubleshooting.md',
)
ENGLISH = (
    'README.md', 'CHANGELOG.md', 'CONTRIBUTING.md', 'docs/README.md',
    'docs/en/GTDB-database.md', 'docs/cli-reference.md',
    'docs/config-reference.md', 'docs/troubleshooting.md',
)


def strip_code(text):
    """Remove fenced blocks so heading/link parsing does not see code."""
    out, fence = [], False
    for line in text.splitlines():
        if FENCE.match(line):
            fence = not fence
            out.append('')
            continue
        out.append('' if fence else line)
    return '\n'.join(out)


def slug(text):
    text = re.sub(r'`([^`]*)`', r'\1', text)
    text = re.sub(r'!?\[([^\]]*)\]\([^)]*\)', r'\1', text)
    text = re.sub(r'[*_~]', '', text)
    text = text.strip().lower()
    text = re.sub(r'[^\w\s\u3400-\u4dbf\u4e00-\u9fff-]', '', text)
    text = re.sub(r'\s+', '-', text)
    return text


def headings(text):
    seen, out = {}, set()
    for line in strip_code(text).splitlines():
        m = re.match(r'^(#{1,6})\s+(.*?)\s*$', line)
        if not m:
            continue
        base = slug(m.group(2))
        n = seen.get(base, 0)
        seen[base] = n + 1
        out.add(base if n == 0 else f'{base}-{n}')
    return out


def check_links(root, rel):
    problems = []
    path = join(root, rel)
    text = open(path, encoding='utf-8').read()
    text_nc = strip_code(text)
    target_heads = {}
    for m in LINK.finditer(text_nc):
        target = m.group(2)
        if target.startswith(('http://', 'https://', 'mailto:', '#')):
            if not target.startswith('#'):
                continue
            file_part, anchor = '', target[1:]
        else:
            file_part, _, anchor = target.partition('#')
        if file_part:
            dest = abspath(join(dirname(path), unquote(file_part)))
            if not os.path.exists(dest):
                problems.append(f'{rel}: link target does not exist: {target}')
                continue
            if anchor and isfile(dest):
                tpath = dest
            elif anchor:
                continue
            else:
                continue
        else:
            tpath = path
        if not anchor:
            continue
        if tpath not in target_heads:
            target_heads[tpath] = headings(open(tpath, encoding='utf-8').read())
        if anchor not in target_heads[tpath]:
            problems.append(f'{rel}: anchor does not exist: {target}')
    return problems


def check_english(root, rel):
    problems = []
    text = open(join(root, rel), encoding='utf-8').read()
    in_fence = False
    for n, line in enumerate(text.splitlines(), 1):
        if FENCE.match(line):
            in_fence = not in_fence
            continue
        if not CJK.search(line):
            continue
        if any(token in line for token in ALLOWED_CJK):
            continue
        problems.append(f'{rel}:{n}: CJK text outside the allow-list: {line.strip()[:100]}')
    return problems


def check_forbidden(root, rel):
    problems = []
    text = open(join(root, rel), encoding='utf-8').read()
    for n, line in enumerate(text.splitlines(), 1):
        for name, rx in FORBIDDEN:
            if rx.search(line):
                problems.append(f'{rel}:{n}: forbidden internal reference ({name}): {line.strip()[:100]}')
    return problems


def check_vars(root, rel):
    problems = []
    text = open(join(root, rel), encoding='utf-8').read()
    defined = set(KNOWN_VARS)
    for line in text.splitlines():
        for m in re.finditer(r'\b(?:export\s+)?([A-Za-z_][A-Za-z0-9_]*)=', line):
            defined.add(m.group(1))
    in_fence = False
    for n, line in enumerate(text.splitlines(), 1):
        if FENCE.match(line):
            in_fence = not in_fence
            continue
        if not in_fence:
            continue
        for name in VAR.findall(line):
            if name in defined:
                continue
            problems.append(f'{rel}:{n}: shell variable ${name} is not defined in this document')
    return problems


def main():
    root = abspath(sys.argv[1] if len(sys.argv) > 1 else '.')
    missing = [d for d in DOCS if not isfile(join(root, d))]
    problems = [f'missing document: {d}' for d in missing]
    for rel in DOCS:
        if rel in missing:
            continue
        problems += check_links(root, rel)
        problems += check_forbidden(root, rel)
        if rel in ENGLISH:
            problems += check_english(root, rel)
    for rel in ('README.md', 'README.zh-CN.md', 'docs/GTDB-数据库.md',
                'docs/en/GTDB-database.md', 'docs/cli-reference.md',
                'docs/config-reference.md', 'docs/troubleshooting.md',
                'CONTRIBUTING.md'):
        if isfile(join(root, rel)):
            problems += check_vars(root, rel)
    if problems:
        for p in problems:
            print('PROBLEM', p)
        print(f'{len(problems)} problem(s)')
        return 1
    print(f'{len(DOCS)} document(s) checked: 0 problem(s)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
