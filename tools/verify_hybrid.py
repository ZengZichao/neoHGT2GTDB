#!/usr/bin/env python3
"""IT-5: merge the real GTDB r220 tree with the real NCBI taxdump.

GTDB has no eukaryotes or viruses, so reproducing the paper's Asgard scenario
needs a hybrid tree. This checks the merge against genuine data:

  * single root, so find_lca() stays usable;
  * GTDB domains grafted under the NCBI cellular root (131567);
  * cross-source LCA behaves as documented (plan RK7);
  * the NCBI subtree source and version are recorded.

Reads the GTDB dataset read-only. Writes only the report path given by --out.
"""

import argparse
import json
import os
import sys
from os import makedirs
from os.path import join

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from neoHGT import gtdb
from neoHGT.util import (find_lca, get_lineage, read_taxdump,
                         taxid_at_rank, timestamp)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gtdb-dir', default=os.environ.get('GTDB_DIR'))
    ap.add_argument('--ncbi-taxdump', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--no-selfcheck', action='store_true')
    args = ap.parse_args()

    result = {'started': str(timestamp()), 'gtdb_dir': args.gtdb_dir,
              'ncbi_taxdump_dir': args.ncbi_taxdump}

    ncbi = read_taxdump(args.ncbi_taxdump)
    result['ncbi_nodes'] = len(ncbi)
    top = {t['rank'] for t in ncbi.values() if t['rank'] in
           ('superkingdom', 'domain')}
    names = sorted({t['name'] for tid, t in ncbi.items()
                    if t['rank'] in ('superkingdom', 'domain')})
    result['ncbi_top_rank_spelling'] = sorted(top)
    result['ncbi_superkingdoms'] = names[:20]

    prov = gtdb.GtdbProvider(args.gtdb_dir, selfcheck=not args.no_selfcheck)
    prov.prepare()
    result['gtdb_release'] = prov.info['release']
    result['gtdb_layout'] = prov.info['layout']
    result['gtdb_selfcheck'] = (prov.report or {}).get('status')

    h = gtdb.HybridProvider(
        prov, args.ncbi_taxdump,
        ncbi_superkingdoms=('Eukaryota', 'Viruses'), lca_policy='root')
    td, tmap, stats = h.load_taxonomy()
    result['gtdb_nodes'] = len(prov.taxdump)
    result['merged_nodes'] = len(td)
    result['ncbi_nodes_grafted'] = stats.get('n_ncbi_nodes')

    # single root
    root_lca = find_lca(list(td), td)
    result['single_root'] = (root_lca == gtdb.ROOT_TID)
    bacteria = tmap['d__Bacteria']
    archaea = tmap['d__Archaea']
    result['gtdb_domain_parent'] = td[bacteria]['parent']
    result['domains_reattached_to_cellular_root'] = (
        td[bacteria]['parent'] == gtdb.CELLULAR_ROOT_TID
        and td[archaea]['parent'] == gtdb.CELLULAR_ROOT_TID)

    euk = next(tid for tid, t in td.items()
               if t['name'] == 'Eukaryota' and tid.isdigit()
               and int(tid) < gtdb.DEFAULT_ID_BASE)
    lca_cross = find_lca([bacteria, euk], td)
    result['lca_bacteria_vs_eukaryota'] = {
        'taxid': lca_cross, 'name': td[lca_cross]['name'],
        'rank': td[lca_cross]['rank']}
    # a bacterium and a eukaryote may only meet at/above the cellular root;
    # that is the RK7 semantic shift the plan warns about
    result['cross_source_lca_is_at_or_above_cellular_root'] = lca_cross in (
        gtdb.CELLULAR_ROOT_TID, gtdb.ROOT_TID)

    # GTDB ranks must still resolve inside the merged tree
    species = next(k for k in tmap if k.startswith('d__Bacteria')
                   and k.count(';') == 6)
    tid = tmap[species]
    result['gtdb_species'] = species.split(';')[-1]
    result['rank_lookups'] = {
        r: (taxid_at_rank(tid, r, td) or 'None') for r in gtdb.GTDB_RANKS}
    # in memory GTDB keeps its own rank names; once written to disk the domain
    # rank becomes 'superkingdom' for DIAMOND's sake, so verify both views
    import tempfile
    tmp = tempfile.mkdtemp(prefix='hyb_td_', dir=os.environ.get('TMPDIR'))
    gtdb.write_taxdump(td, tmp)
    td_disk = read_taxdump(tmp)
    result['disk_domain_rank'] = td_disk[bacteria]['rank']
    result['rank_lookup_domain_ok'] = bool(
        taxid_at_rank(tid, 'domain', td)) and bool(
        taxid_at_rank(tid, 'superkingdom', td_disk))
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
    result['lineage_depth_in_merged_tree'] = len(get_lineage(tid, td))
    # DIAMOND validates the rank column and aborts makedb on NCBI's newer rank
    # spellings ('cellular root', 'clade'), so the file a hybrid build hands to
    # diamond must be free of them while the in-memory tree keeps them
    result['ranks_flattened_for_diamond'] = \
        stats.get('ranks_flattened_for_diamond', {})
    result['nodes_dmp_is_diamond_legal'] = bool(
        {t['rank'] for t in td_disk.values()} <= set(gtdb.DIAMOND_RANKS))
    result['in_memory_ranks_preserved'] = bool(
        td[gtdb.CELLULAR_ROOT_TID]['rank'] == 'cellular root'
        or gtdb.CELLULAR_ROOT_TID not in td)
    # the cellular root is what the whole graft hangs from: its absence made
    # every later check vacuous while the verdict still read PASS
    result['cellular_root_present'] = gtdb.CELLULAR_ROOT_TID in td

    # lca policies
    pairs = [('s1', str(gtdb.DEFAULT_ID_BASE)), ('s2', euk)]
    for policy in ('root', 'domain', 'refuse'):
        hp = gtdb.HybridProvider(prov, args.ncbi_taxdump,
                                ncbi_superkingdoms=('Eukaryota', 'Viruses'),
                                lca_policy=policy)
        hp.taxdump, hp.taxid_map = td, tmap
        result[f'policy_{policy}_kept'] = len(hp.resolve_shared_protein(pairs))

    # the three policies must produce three *different* answers on a shared
    # protein pair; recording the counts without judging them is how the
    # round-30 "all three behave identically" defect kept printing PASS
    policies = [result.get(f'policy_{p}_kept', 0)
                for p in ('root', 'domain', 'refuse')]
    result['policies_distinguishable'] = len(set(policies)) == 3

    result['finished'] = str(timestamp())
    ok = (result['single_root']
          and result['cellular_root_present']
          and result['domains_reattached_to_cellular_root']
          and result['cross_source_lca_is_at_or_above_cellular_root']
          and result['rank_lookup_domain_ok']
          and result['nodes_dmp_is_diamond_legal']
          and result['in_memory_ranks_preserved']
          and result['policies_distinguishable'])
    result['verdict'] = 'PASS' if ok else 'FAIL'
    out = gtdb.guard_write(args.out)
    makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(result, f, indent=2, default=str)
    print(json.dumps(result, indent=2, default=str))
    print('IT-5 verdict:', result['verdict'])
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
