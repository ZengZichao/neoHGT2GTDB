#!/usr/bin/env python
"""Production reachability: can a command line actually arrive at this function?

Execution coverage answers "did any test call this?", which is not the question
a user cares about. `HybridProvider.resolve_shared_protein` was called by tests,
printed its setting into the build report, and was never invoked by any build
path (round 30). This tool answers the product-level question:

  1. take the functions this adaptation *added* (diffed against the pristine
     upstream tree, so upstream's own dead code is nobody's problem here);
  2. build a name-level call graph over the package, including the edges that
     Python resolves only at runtime -- `Class.method = func` rebinding, a
     function handed to another function as a value, and attribute calls;
  3. start from the entry points `scripts/neoHGT` dispatches to and report every
     added function that no entry point can reach.

A negative control (rewriting one real call site into a dummy) proves the
detector can see the failure it claims to catch.

    usage: python tools/verify_production_reachability.py \
               --code . --baseline ../neoHGT-main --out report.json
"""
import argparse
import ast
import json
import os
import re
import sys
from os.path import basename, exists, isdir, isfile, join

MODULES = ['neoHGT/gtdb.py', 'neoHGT/util.py', 'neoHGT/analyze.py',
           'neoHGT/database.py', 'neoHGT/search.py']

# what scripts/neoHGT can start
SEEDS = {
    'main', 'cli_main', 'call_with_source', '_ncbi_call', 'build_from_source',
    'compile_gtdb_database', 'Database.__call__', 'Search.__call__',
    'Analyze.__call__', 'predict_hgt_gtdb', 'taxfilt_wf_gtdb',
    'parse_def_table_checked', 'set_params_gtdb', 'taxid_wf_gtdb',
    'write_hgt_list_gtdb', 'cmd_probe', 'cmd_taxdump', 'cmd_gtdbtk',
    'cmd_map_report', 'run', 'write_gtdb_outputs', 'finalize_gtdb_faa',
    'plan_remote_fetches', 'extract_genomes_stream', 'extract_genomes',
    'build_taxonmap', 'compile_database', 'selfaln_wf', 'search_wf',
    'input_wf', 'args_wf', 'read_input', 'read_input_gtdb', 'make_db_from_files',
    'update_dmp_guard', 'bind_stream_extraction', 'extract_genomes_stream_ncbi',
    'download_parallel', 'download_genomes', 'genome_lineages',
    'genome_metadata', 'get_summary', 'sample_by_taxonomy', 'filter_genomes',
    'build_taxdump', 'get_categories', 'report_diff', 'retrieve_categories',
    'build_blast_db', 'build_diamond_db',
}


def _read(path):
    with open(path, encoding='utf-8', errors='replace') as fh:
        return fh.read()


def top_level_defs(text):
    """(qualified name, name) for module-level functions and methods only.

    Nested functions are left out on purpose: they are reached by calling their
    enclosing function, and a name-level graph has no way to resolve a closure,
    so counting them would only produce noise like `flush` or `progress`.
    """
    out = {}
    tree = ast.parse(text)
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out[node.name] = node.lineno
        elif isinstance(node, ast.ClassDef):
            for m in node.body:
                if isinstance(m, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    out[f'{node.name}.{m.name}'] = m.lineno
    return out


def edges(text):
    """caller -> set(call names), including values passed as arguments."""
    graph = {}
    tree = ast.parse(text)

    def note(caller, name):
        graph.setdefault(caller, set()).add(name)

    class Walk(ast.NodeVisitor):
        def __init__(self, caller):
            self.caller = caller

        def visit_Call(self, node):
            f = node.func
            if isinstance(f, ast.Name):
                note(self.caller, f.id)
            elif isinstance(f, ast.Attribute):
                note(self.caller, f.attr)
            self.generic_visit(node)

        def visit_Attribute(self, node):
            # obj.method used as a value (callback, rebinding) still names
            # `method`; without this edge the round-30 wiring would look dead
            if isinstance(node.ctx, ast.Load):
                note(self.caller, node.attr)
            # keep descending: in `f(x).get(y)` the receiver is itself a call
            self.generic_visit(node)

        def visit_Name(self, node):
            if isinstance(node.ctx, ast.Load):
                note(self.caller, node.id)

    for top in tree.body:
        if isinstance(top, ast.ClassDef):
            for m in top.body:
                if isinstance(m, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    Walk(f'{top.name}.{m.name}').visit(m)
            for stmt in top.body:
                if not isinstance(stmt, (ast.FunctionDef,
                                         ast.AsyncFunctionDef)):
                    Walk(top.name).visit(stmt)
        elif isinstance(top, (ast.FunctionDef, ast.AsyncFunctionDef)):
            Walk(top.name).visit(top)
        else:
            Walk('<module>').visit(top)
    return graph


def rebindings(text):
    """`Class.slot = func` and `slot = func` at module level."""
    dyn, alias = {}, {}
    for node in ast.parse(text).body:
        if not isinstance(node, ast.Assign) or not isinstance(node.value,
                                                              ast.Name):
            continue
        for t in node.targets:
            if isinstance(t, ast.Name):
                alias[t.id] = node.value.id
            elif isinstance(t, ast.Attribute) and isinstance(t.value,
                                                             ast.Name):
                key = f'{t.value.id}.{t.attr}'
                dyn.setdefault(key, set()).add(node.value.id)
                dyn.setdefault(t.attr, set()).add(node.value.id)
    return dyn, alias


def cli_seeds(code_dir):
    """Entry points derived from the shipped CLI, not from a hand-kept list.

    `scripts/neoHGT` instantiates one object per sub-command and calls it, so
    every capitalised name it references contributes a `<Class>.__call__` seed.
    Hard-coding this list is how `Gtdb.__call__` -- the `neoHGT gtdb`
    sub-command -- went missing and made its helpers look dead.
    """
    seeds = set(SEEDS)
    script = join(code_dir, 'scripts', 'neoHGT')
    if not isfile(script):
        return seeds
    tree = ast.parse(_read(script))
    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and node.id[:1].isupper():
            seeds.add(f'{node.id}.__call__')
            seeds.add(node.id)
        # the CLI builds its sub-command objects as
        # getattr(import_module('neoHGT.' + name), name.capitalize())(), so the
        # class names never appear as literals; recover them from the tuple of
        # module names instead, which is the same source of truth
        if (isinstance(node, ast.Assign)
                and any(getattr(t, 'id', '') == 'module_names'
                        for t in node.targets)):
            try:
                for nm in ast.literal_eval(node.value):
                    cls = nm.capitalize()
                    seeds |= {f'{cls}.__call__', f'{cls}.main', 'main'}
            except (ValueError, SyntaxError):
                pass
    return seeds


def import_aliases(text):
    """`from x import y as z` -> {z: y}; an aliased call still names the origin."""
    out = {}
    for node in ast.walk(ast.parse(text)):
        # `from util import resolve_taxon_ref as _resolve_ref` happens inside a
        # method in analyze.py, so this cannot be limited to module level
        if isinstance(node, ast.ImportFrom):
            for a in node.names:
                if a.asname:
                    out[a.asname] = a.name
    return out


def analyse(code_dir, baseline=None, added_list=None):
    """(unreachable, reachable_count, added_count).

    `baseline` diffs against a pristine upstream tree; `added_list` is the same
    information frozen into a committed file, for checkouts that have only the
    adapted code (CI). The frozen list is intersected with what exists now, so
    a deleted function cannot keep a slot open forever.
    """
    graph, defs, dyn, alias = {}, {}, {}, {}
    aliases = {}
    for rel in MODULES:
        path = join(code_dir, rel)
        if not exists(path):
            continue
        text = _read(path)
        defs[rel] = top_level_defs(text)
        for caller, names in edges(text).items():
            graph.setdefault(caller, set()).update(names)
        d, a = rebindings(text)
        for k, v in d.items():
            dyn.setdefault(k, set()).update(v)
        alias.update(a)
        aliases.update(import_aliases(text))

    added = defs
    if added_list:
        wanted = set(added_list)
        added = {rel: {n: ln for n, ln in names.items()
                       if n in wanted or n.split('.')[-1] in wanted}
                 for rel, names in defs.items()}
        added = {k: v for k, v in added.items() if v}
    elif baseline:
        # only what this adaptation introduced is this tool's business
        added = {}
        for rel, names in defs.items():
            bpath = join(baseline, rel)
            base = top_level_defs(_read(bpath)) if isfile(bpath) else {}
            added[rel] = {n: ln for n, ln in names.items()
                          if n not in base and n.split('.')[-1] not in base}

    seen, stack = set(), sorted(cli_seeds(code_dir))
    while stack:
        cur = stack.pop()
        if cur in seen:
            continue
        seen.add(cur)
        targets = {cur, alias.get(cur, '')}
        targets |= dyn.get(cur, set())
        # `from util import resolve_taxon_ref as _resolve_ref` : a call through
        # the local alias names the same function
        if cur in aliases:
            targets.add(aliases[cur])
        if '.' in cur:
            targets.add(cur.split('.')[-1])
        else:
            targets |= {k for k in graph if k.endswith('.' + cur)}
        # a slot the code dispatches through (`Analyze.set_parameters = ...`) is
        # a real entry point to the wrapper: walking only the callee lists left
        # such a function looking unreachable, because nothing ever calls it by
        # its own name -- the rebinding *is* the call
        for t in dyn.get(cur, set()) | {alias.get(cur, '')} \
                | ({aliases[cur]} if cur in aliases else set()):
            if t and t not in seen:
                stack.append(t)
        for t in targets:
            for callee in graph.get(t, ()):
                stack.append(callee)

    unreachable = {}
    for rel, names in added.items():
        miss = sorted(n for n in names
                      if n not in seen and n.split('.')[-1] not in seen
                      and not n.split('.')[-1].startswith('_'))
        if miss:
            unreachable[rel] = miss
    return unreachable, len(seen), sum(len(v) for v in added.values())


def control(code_dir):
    """Delete a real call site and prove the resolver drops out of reach."""
    text = _read(join(code_dir, 'neoHGT/database.py'))
    mutated = text.replace('resolve=(provider.resolve_shared_protein',
                           'resolve=(lambda _p: _p')
    if mutated == text:
        return {'applied': False,
                'note': 'the round-30 call site is gone; this control can no '
                        'longer prove anything and must be retargeted'}
    graph, _ = {}, {}
    # re-run the analysis with the mutated file substituted
    import tempfile
    with tempfile.TemporaryDirectory(dir=os.environ.get('NEOHGT_TMP', '/tmp')) \
            as tmp:
        for rel in MODULES:
            dst = join(tmp, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            src = _read(join(code_dir, rel))
            with open(dst, 'w', encoding='utf-8') as fh:
                fh.write(mutated if rel == 'neoHGT/database.py' else src)
        # copy the rest of the package so imports of the tree stay complete,
        # including the CLI script that supplies the seeds
        for sub in ('scripts',):
            src_dir = join(code_dir, sub)
            if isdir(src_dir):
                dst_dir = join(tmp, sub)
                os.makedirs(dst_dir, exist_ok=True)
                for fn in sorted(os.listdir(src_dir)):
                    # a __pycache__ left by any earlier import crashed this
                    # copy loop with IsADirectoryError, which read as a product
                    # regression instead of an environment artifact
                    if not os.path.isfile(join(src_dir, fn)):
                        continue
                    with open(join(src_dir, fn), 'rb') as fi, \
                            open(join(dst_dir, fn), 'wb') as fo:
                        fo.write(fi.read())
        un, _, _ = analyse(tmp)
    found = 'HybridProvider.resolve_shared_protein' in un.get('neoHGT/gtdb.py',
                                                             [])
    return {'applied': True,
            'what': 'rewrite build_from_source`s `resolve='
                    'provider.resolve_shared_protein` into a dummy lambda',
            'resolver_becomes_unreachable': found,
            'baseline_unreachable': un.get('neoHGT/gtdb.py', [])}


def _added_names(code_dir, baseline, added_list):
    """The {module: [names]} map the reachability check judges."""
    defs = {}
    for rel in MODULES:
        path = join(code_dir, rel)
        if not exists(path):
            continue
        names = top_level_defs(_read(path))
        if added_list:
            wanted = set(added_list)
            names = {n: ln for n, ln in names.items()
                     if n in wanted or n.split('.')[-1] in wanted}
        elif baseline and isdir(baseline):
            bpath = join(baseline, rel)
            base = top_level_defs(_read(bpath)) if isfile(bpath) else {}
            names = {n: ln for n, ln in names.items()
                     if n not in base and n.split('.')[-1] not in base}
        if names:
            defs[rel] = names
    return defs


def defs_by_module(code_dir):
    """{module: {name: lineno}} of every top-level definition that exists now."""
    out = {}
    for rel in MODULES:
        path = join(code_dir, rel)
        if exists(path):
            out[rel] = top_level_defs(_read(path))
    return out


def frozen_lists(path):
    """(added names, upstream names) from the committed list.

    The upstream half is what makes the list self-checking: with only the added
    names, a function introduced after the list was frozen belongs to no set and
    is never asked about -- the gate keeps reporting "0 problems" while silently
    judging a smaller and smaller share of the code.
    """
    with open(path) as fh:
        data = json.load(fh)
    return set(data.get('functions', [])), set(data.get('upstream', []))


def untracked_names(defs, added, upstream):
    """Definitions that are in neither the added list nor the upstream list."""
    out = {}
    for rel, names in defs.items():
        odd = [n for n in names
               if n not in added and n.split('.')[-1] not in added
               and n not in upstream and n.split('.')[-1] not in upstream]
        if odd:
            out[rel] = sorted(odd)
    return out


def list_drift(defs, baseline, added, upstream):
    """(missing, stale, upstream_missing) between the tree and the frozen lists.

    Local-only: recomputing the diff needs the pristine upstream tree, which CI
    does not have -- which is exactly why the upstream half is committed.
    """
    missing, stale, up_missing = [], [], []
    fresh_added, fresh_up = set(), set()
    for rel, names in defs.items():
        bpath = join(baseline, rel)
        base = top_level_defs(_read(bpath)) if isfile(bpath) else {}
        for n in names:
            if n in base or n.split('.')[-1] in base:
                fresh_up.add(n)
            else:
                fresh_added.add(n)
        for n in base:
            if n not in names and n not in upstream:
                up_missing.append(f'{rel}:{n}')
    missing = sorted(fresh_added - added)
    stale = sorted(n for n in added if n not in
                   {x for ns in defs.values() for x in ns})
    return missing, stale, sorted(set(up_missing))


def declared_exceptions(doc_path):
    """Names registered as intentionally unreachable in the deliverables doc."""
    if not isfile(doc_path):
        return set()
    text = _read(doc_path)
    sect = text.split('可达性偏差', 1)[1] if '可达性偏差' in text else ''
    return {m.group(1) for m in re.finditer(r'`([A-Za-z_][\w.]*)\(\)`', sect)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--code', default='.', help='adapted package root')
    ap.add_argument('--baseline', default=None,
                    help='pristine upstream tree (default: $NEOHGT_MAIN)')
    ap.add_argument('--deliverables', default=None)
    ap.add_argument('--out', default=None)
    ap.add_argument('--control', action='store_true',
                    help='also run the negative control (slower)')
    ap.add_argument('--emit-added', default=None,
                    help='write the list of added functions and exit (for CI)')
    ap.add_argument('--added', default=None,
                    help='committed list of added functions to judge')
    ap.add_argument('--check-list', action='store_true',
                    help='also compare the committed list with a fresh diff '
                         'against --baseline (needs the upstream tree)')
    args = ap.parse_args()
    code = os.path.abspath(args.code)
    base = args.baseline or os.environ.get('NEOHGT_MAIN') or \
        join(os.path.dirname(code), 'neoHGT-main')
    frozen, frozen_up = None, None
    if args.added and isfile(args.added):
        frozen, frozen_up = frozen_lists(args.added)
    unreachable, reach, added = analyse(
        code, base if isdir(base) else None, frozen)
    defs = defs_by_module(code)
    if args.emit_added:
        names = sorted({n for ns in
                        _added_names(code, base, None).values() for n in ns})
        ups = sorted({n for n in {x for ns in defs.values() for x in ns}
                      if n not in names})
        with open(args.emit_added, 'w') as fh:
            json.dump({'generated_from': os.path.abspath(code),
                       'baseline': base if isdir(base) else None,
                       'functions': names, 'upstream': ups,
                       'counts': {'added': len(names), 'upstream': len(ups)}},
                      fh, indent=1)
        print(f'{len(names)} added functions (+{len(ups)} upstream names) '
              f'-> {args.emit_added}')
        return 0
    declared = declared_exceptions(
        args.deliverables or join(code, 'docs', '交付物清单.md'))
    unjudged = {rel: [n for n in names if n not in declared
                      and n.split('.')[-1] not in declared]
                for rel, names in unreachable.items()}
    unjudged = {k: v for k, v in unjudged.items() if v}
    untracked = (untracked_names(defs, frozen, frozen_up)
                 if frozen and frozen_up else {})
    drift = None
    if args.check_list:
        if not isdir(base):
            raise SystemExit('--check-list needs the pristine upstream tree')
        drift = list_drift(defs, base, frozen or set(), frozen_up or set())
    report = {
        'purpose': 'functions this adaptation added that no CLI entry point can '
                   'reach (the class of defect test coverage cannot see)',
        'baseline': base if isdir(base) else None,
        'added_functions': added,
        'reachable_names': reach,
        'unreachable': unreachable,
        'declared_exceptions': sorted(declared),
        'unjudged': unjudged,
        'untracked_names': untracked,
        'list_drift': ({'missing_from_list': drift[0], 'stale_in_list': drift[1],
                        'upstream_changed': drift[2]} if drift else None),
    }
    if args.control:
        report['negative_control'] = control(code)
    out = args.out or join(code, '..', 'neoHGT2GTDB-缓存测试', 'reports',
                           'production_reachability.json')
    out = os.path.abspath(out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as fh:
        json.dump(report, fh, indent=1, ensure_ascii=False)
    print(f'added functions: {added}; reachable names: {reach}')
    for rel, names in unreachable.items():
        print(f'  {rel}: {", ".join(names)}')
    problems = sum(len(v) for v in unjudged.values())
    extra = sum(len(v) for v in untracked.values())
    print('problems: %d undeclared%s' % (
        problems,
        (' + %d untracked (in neither the added nor the upstream list)' % extra
         if extra else '')))
    for rel, names in untracked.items():
        print(f'  UNTRACKED {rel}: {", ".join(names)}')
    if drift:
        print(f'list drift: missing {len(drift[0])}, stale {len(drift[1])}, '
              f'upstream changed {len(drift[2])}')
        for kind, items in zip(('missing', 'stale', 'upstream'), drift):
            for x in items[:10]:
                print(f'  DRIFT {kind}: {x}')
        problems += sum(len(x) for x in drift)
    if 'negative_control' in report:
        nc = report['negative_control']
        print('negative control:', nc.get('applied'), '->',
              nc.get('resolver_becomes_unreachable'))
    print('wrote', out)
    return 1 if problems else 0


if __name__ == '__main__':
    raise SystemExit(main())
