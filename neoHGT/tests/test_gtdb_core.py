"""Unit and property tests for neoHGT.gtdb (UT-0.x, UT-1.x, UT-3.1b).

The fixtures here are synthetic except where a test explicitly marks itself
as requiring the local GTDB dataset (which is read-only and never copied).
"""

import os
import shutil
import tempfile
import unittest
from os.path import isfile, join

from neoHGT import gtdb
from neoHGT import util
from neoHGT.util import read_taxdump, find_lca, taxid_at_rank

LIN_ENT = ('d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;'
           'o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;'
           's__Escherichia coli')
LIN_SALM = ('d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;'
            'o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;'
            's__Salmonella enterica')
LIN_ARCH = ('d__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Nanoarchaeales;'
            'f__Nanoarchaeaceae;g__Nanoarchaeum;s__Nanoarchaeum equitans')
LIN_CLUST = ('d__Bacteria;p__Bacillota_I;c__Bacilli_A;o__Bacillales;'
             'f__Bacillaceae_A;g__Bacillus;s__Bacillus_A cereus')
LIN_PLACE = ('d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;'
             'o__ML615J-28;f__CAG-698;g__336-3;s__336-3 bacterium BG-634')

ALL_LINS = [LIN_ENT, LIN_ENT, LIN_SALM, LIN_ARCH, LIN_CLUST, LIN_PLACE]

STEM_ENT = ('GB_GCA_000008885.1_d_Bacteria_p_Pseudomonadota_c_'
            'Gammaproteobacteria_o_Enterobacterales_f_Enterobacteriaceae_A'
            '_g_Wigglesworthia')
STEM_I1 = ('RS_GCF_000008885.2_d_Bacteria_p_Bacillota_I_c_Bacilli_A_o_'
           'Mycoplasmatales_f_Mycoplasmataceae_g_Mycoplasma')
STEM_I2 = ('GB_GCA_937958695.1_d_Bacteria_p_DUMJ01_c_DUMJ01_o_DUMJ01_f_'
           'DUMJ01_g_DUMJ01')
STEM_I3 = ('GB_GCA_000063505.1_d_Bacteria_p_Cyanobacteriota_c_Cyanobacteriia'
           '_o_PCC-6307_f_Cyanobiaceae_g_Synechococcus_C')
STEM_I4 = ('GB_GCA_000169075.1_d_Bacteria_p_Pseudomonadota_c_Gammapro;'
           'teobacteria_o_Pseudomonadales_f_DSM-21967_g_Oceanicoccus')


class TestParseStem(unittest.TestCase):
    """UT-0.2: local file-name stems stop at genus and carry cluster letters."""

    def test_basic(self):
        got = gtdb.parse_stem(STEM_ENT)
        self.assertEqual(got['accession'], 'GB_GCA_000008885.1')
        self.assertEqual(got['domain'], 'Bacteria')
        self.assertEqual(got['phylum'], 'Pseudomonadota')
        self.assertEqual(got['family'], 'Enterobacteriaceae_A')
        self.assertEqual(got['genus'], 'Wigglesworthia')

    def test_no_species_is_invented(self):
        for stem in (STEM_ENT, STEM_I1, STEM_I2, STEM_I3):
            self.assertNotIn('species', gtdb.parse_stem(stem))

    def test_cluster_letters_and_hyphens(self):
        got = gtdb.parse_stem(STEM_I1)
        self.assertEqual(got['phylum'], 'Bacillota_I')
        self.assertEqual(got['class'], 'Bacilli_A')
        self.assertEqual(got['genus'], 'Mycoplasma')
        got = gtdb.parse_stem(STEM_I3)
        self.assertEqual(got['order'], 'PCC-6307')
        self.assertEqual(got['genus'], 'Synechococcus_C')

    def test_placeholder_stems(self):
        got = gtdb.parse_stem(STEM_I2)
        self.assertEqual([got[r] for r in gtdb.GTDB_RANKS[:6]],
                         ['Bacteria'] + ['DUMJ01'] * 5)

    def test_bad_stem_raises(self):
        with self.assertRaises(ValueError):
            gtdb.parse_stem('GB_GCA_000008885.1_d_Bacteria_p_X')


class TestSeqIdCoverageProperty(unittest.TestCase):
    """AT-4 / DoD-3: 100% of protein ids must parse and map to a taxId.

    The invariant is checked over the committed fixture proteomes and over
    randomized id variants of every fixture stem, so a change in the id
    conventions cannot pass by only handling the handful of cases above.
    """

    FXDIR = join(os.path.dirname(os.path.abspath(__file__)),
                 'data', 'GTDB', 'protein')

    def _stems(self):
        from os.path import isfile, splitext
        if not os.path.isdir(self.FXDIR):
            self.skipTest(f'no fixture proteomes in {self.FXDIR}')
        return sorted(splitext(f)[0] for f in os.listdir(self.FXDIR)
                      if isfile(join(self.FXDIR, f)) and f.endswith('.faa'))

    def test_every_real_sequence_id_resolves(self):
        import json
        from os.path import basename
        fx = join(os.path.dirname(self.FXDIR), 'fixtures.json')
        expected = 0
        if isfile(fx):
            with open(fx) as f:
                manifest = json.load(f)
            expected = sum(g['kb_records'] for g in manifest['genomes'])
        stems = self._stems()
        s2t = {stem: str(100000000 + i) for i, stem in enumerate(stems)}
        n = 0
        for fname in sorted(os.listdir(self.FXDIR)):
            if not fname.endswith('.faa'):
                continue
            with open(join(self.FXDIR, fname)) as f:
                for line in f:
                    if not line.startswith('>'):
                        continue
                    sid = line[1:].split()[0]
                    got = gtdb.taxid_of_seqid(sid, stem2taxid=s2t)
                    self.assertIsNotNone(
                        got, f'{sid}: no taxId via route C')
                    n += 1
        self.assertGreater(n, 0)
        if expected:
            # every sequence counted in the manifest must resolve: the property
            # holds over the whole fixture set, not just the ones we happened to read
            self.assertEqual(n, expected)

    def test_randomised_variants_of_every_stem_resolve(self):
        import random
        stems = self._stems()
        s2t = {stem: str(100000000 + i) for i, stem in enumerate(stems)}
        rng = random.Random(42)
        for stem in stems:
            for _ in range(50):
                contig = f'{rng.choice(["NC", "BA", "CP"])}_' \
                         f'{rng.randrange(1, 10_000_000)}.{rng.randrange(1, 9)}'
                sid = f'{stem}{gtdb.SERIAL_ANCHOR}{contig}_{rng.randrange(1, 9999)}'
                self.assertEqual(gtdb.taxid_of_seqid(sid, stem2taxid=s2t),
                                 s2t[stem], sid)

    def test_unmappable_id_is_reported_not_guessed(self):
        stems = self._stems()
        s2t = {stem: str(100000000 + i) for i, stem in enumerate(stems)}
        ghost = ('GB_GCA_000000000.1_d_Bacteria_p_X_c_X_o_X_f_X_g_X'
                 + gtdb.SERIAL_ANCHOR + 'NC_000001.1_1')
        self.assertIsNone(gtdb.taxid_of_seqid(ghost, stem2taxid=s2t))


class TestParseSeqId(unittest.TestCase):
    """UT-3.1b: route C anchors for both Id conventions."""

    def test_local_anchor(self):
        sid = STEM_ENT + gtdb.SERIAL_ANCHOR + 'BA000021.3_1'
        got = gtdb.parse_seqid(sid)
        self.assertEqual(got['stem'], STEM_ENT)
        self.assertEqual(got['genome'], 'GB_GCA_000008885.1')

    def test_official_anchor(self):
        got = gtdb.parse_seqid('NC_009515.1_2', layout='official-plain')
        self.assertEqual(got['contig'], 'NC_009515.1')
        self.assertEqual(got['locus'], '2')

    def test_unparsable(self):
        self.assertIsNone(gtdb.parse_seqid('nonsense'))

    def test_route_c_lookup(self):
        s2t = {STEM_ENT: '1006000001'}
        self.assertEqual(
            gtdb.taxid_of_seqid(STEM_ENT + gtdb.SERIAL_ANCHOR + 'x.1_5',
                                stem2taxid=s2t), '1006000001')
        self.assertIsNone(
            gtdb.taxid_of_seqid(
                'RS_GCF_999999999.1_d_Bacteria_p_X_c_X_o_X_f_X_g_X'
                + gtdb.SERIAL_ANCHOR + 'y.1_1', stem2taxid=s2t))

    def test_route_c_fallback_via_prefix(self):
        p2g = {'NC_009515.1': 'RS_GCF_000008885.1'}
        g2t = {'RS_GCF_000008885.1': '1006000002'}
        self.assertEqual(
            gtdb.taxid_of_seqid('NC_009515.1_77', layout='official-plain',
                                prefix2genome=p2g, genome2taxid=g2t),
            '1006000002')


class TestSplitLineage(unittest.TestCase):
    """UT-0.x: the seven-rank contract is enforced, not assumed."""

    def test_ok(self):
        got = gtdb.split_lineage(LIN_ENT)
        self.assertEqual([r for r, _ in got], list(gtdb.GTDB_RANKS))
        self.assertEqual(got[-1][1], 's__Escherichia coli')

    def test_incomplete(self):
        with self.assertRaises(ValueError):
            gtdb.split_lineage('d__Bacteria;p__X')

    def test_wrong_prefix(self):
        with self.assertRaises(ValueError):
            gtdb.split_lineage('k__Bacteria;p__X;c__X;o__X;f__X;g__X;s__X')


class TestTaxdumpConversion(unittest.TestCase):
    """UT-1.1, UT-1.2, UT-1.3, UT-1.7."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(dir=os.environ.get('TMPDIR') or None)
        cls.td, cls.tmap, cls.stats = gtdb.build_taxdump_from_gtdb(ALL_LINS)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_node_count_matches_unique_prefix_paths(self):
        uniq = set()
        for lin in set(ALL_LINS):
            uniq.update(gtdb.lineage_paths(lin))
        self.assertEqual(self.stats['n_gtdb_nodes'], len(uniq))
        self.assertEqual(self.stats['n_nodes'], len(uniq) + 1)  # + root

    def test_every_parent_exists_and_root_self_references(self):
        for tid, rec in self.td.items():
            self.assertIn(rec['parent'], self.td)
        self.assertEqual(self.td[gtdb.ROOT_TID]['parent'], gtdb.ROOT_TID)

    def test_rank_counts(self):
        want = {'domain': 2, 'phylum': 3, 'class': 3, 'order': 4,
                'family': 4, 'genus': 5, 'species': 5}
        self.assertEqual(self.stats['by_rank'], want)

    def test_ids_are_in_reserved_int32_range(self):
        for tid in self.tmap.values():
            self.assertGreaterEqual(int(tid), gtdb.DEFAULT_ID_BASE)
            self.assertLess(int(tid), 2 ** 31)
        self.assertEqual(len(set(self.tmap.values())), len(self.tmap))

    def test_manifest_appends_without_renumbering(self):
        man = {k: {'taxid': v, 'rank': '', 'label': '', 'name': '',
                   'ncbi_taxid': ''} for k, v in self.tmap.items()}
        more = ALL_LINS + [
            'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;'
            'o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;'
            's__Escherichia fergusonii']
        _, tmap2, _ = gtdb.build_taxdump_from_gtdb(more, manifest=man)
        for k, v in self.tmap.items():
            self.assertEqual(tmap2[k], v, f'existing Id changed for {k}')
        self.assertGreater(len(tmap2), len(self.tmap))

    def test_release_swap_never_reissues_or_reuses_a_taxid(self):
        """What a real release upgrade does to the Id space, pinned offline.

        Measured between GTDB r214 and r220 (89,172 shared nodes kept their
        Ids in both directions, 24,145 / 57,207 new nodes collided with none of
        them, and no Id of a vanished node was reissued -- see
        缓存测试/reports/R24_cross_release_id_stability.json). Allocation is
        rank-blocked, so a new node legitimately lands *below* the old maximum;
        the invariant that matters is that it never lands on an Id the
        manifest already knows, including one it just lost.
        """
        man = {k: {'taxid': v, 'rank': '', 'label': '', 'name': '',
                   'ncbi_taxid': ''} for k, v in self.tmap.items()}
        old_ids = set(self.tmap.values())

        lin = sorted(set(ALL_LINS))
        dropped, renamed = lin[0], lin[1]
        genus = renamed.rsplit(';', 1)[0]
        next_release = [x for x in lin if x not in (dropped, renamed)]
        next_release += [f'{genus};s__Renamed sp900000001',
                         f'{genus};s__Added sp900000002']

        _, new, stats = gtdb.build_taxdump_from_gtdb(next_release,
                                                     manifest=man)
        shared = set(self.tmap) & set(new)
        self.assertTrue(shared, 'the two releases must overlap at all')
        for k in shared:
            self.assertEqual(new[k], self.tmap[k],
                             f'existing node {k} changed taxId')

        fresh = set(new) - set(self.tmap)
        self.assertTrue(fresh)
        self.assertFalse({new[k] for k in fresh} & old_ids,
                         'a node new to this release inherited a historical '
                         'taxId, so old scores.tsv rows would silently mean a '
                         'different organism')

        gone = set(self.tmap) - set(new)
        # node keys are the lineage prefix strings, so the dropped lineage is
        # itself the species node that disappeared
        self.assertIn(dropped, gone,
                      'the dropped lineage must show up as superseded')
        self.assertFalse({self.tmap[k] for k in gone} & set(new.values()),
                         'the taxId of a vanished node is back in circulation')
        self.assertEqual(stats['n_superseded'], len(gone))
        self.assertEqual(len(set(new.values())), len(new))

    def test_round_trip_through_files(self):
        # UT-1.3: write_taxdump output must read back identically, except that
        # the domain rank is stored as NCBI's 'superkingdom' (DIAMOND rejects
        # 'domain', spike S1).
        out = join(self.tmp, 'taxdump')
        gtdb.write_taxdump(self.td, out)
        back = read_taxdump(out)
        self.assertEqual(set(back), set(self.td))
        for tid in self.td:
            self.assertEqual(back[tid]['parent'], self.td[tid]['parent'])
            self.assertEqual(back[tid]['rank'],
                             gtdb.dmp_rank(self.td[tid]['rank']))
            self.assertEqual(back[tid]['name'], self.td[tid]['name'])

    def test_domain_rank_is_diamond_compatible(self):
        out = join(self.tmp, 'taxdump_alias')
        gtdb.write_taxdump(self.td, out)
        ranks = {v['rank'] for v in read_taxdump(out).values()}
        self.assertIn('superkingdom', ranks)
        self.assertNotIn('domain', ranks)
        # users may still ask for 'domain'
        self.assertEqual(util.rank_for_dump('domain', read_taxdump(out)),
                         'superkingdom')

    def test_names_are_unprefixed(self):
        # UT-1.4 contract: names.dmp must not carry rank prefixes.
        for rec in self.td.values():
            self.assertNotRegex(rec['name'], r'^[dpcofgs]__')
        self.assertIn('Escherichia coli',
                      {v['name'] for v in self.td.values()})
        self.assertIn('Bacillus cereus',
                      {v['name'] for v in self.td.values()})

    def test_lca_matches_lineage_prefix(self):
        # UT-1.7: cross-check find_lca against the longest common prefix.
        import random
        keys = [k for k in self.tmap if path_is_species(k)]
        random.seed(0)
        for a, b in [(x, y) for x in keys for y in keys]:
            want = common_prefix_path(a, b) or gtdb.ROOT_TID
            got_tid = find_lca([self.tmap[a], self.tmap[b]], self.td)
            got = key_of_tid(self.tmap, got_tid)
            self.assertEqual(got, want, f'LCA of {a} and {b}')

    def test_taxid_at_rank_works(self):
        tid = self.tmap[LIN_ENT]
        self.assertEqual(
            self.tmap['d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;'
                      'o__Enterobacterales;f__Enterobacteriaceae'],
            taxid_at_rank(tid, 'family', self.td))
        self.assertEqual(self.tmap['d__Bacteria'],
                         taxid_at_rank(tid, 'domain', self.td))
        self.assertIsNone(taxid_at_rank(tid, 'kingdom', self.td))

    def test_single_root_enforced(self):
        with self.assertRaises(AssertionError):
            gtdb.build_taxdump_from_gtdb(
                [LIN_ENT, LIN_ARCH], root_parent='42')  # parent absent


def path_is_species(k):
    return k.count(gtdb.LINEAGE_SEP) == gtdb.DEPTH - 1


def common_prefix_path(a, b):
    sa, sb = a.split(gtdb.LINEAGE_SEP), b.split(gtdb.LINEAGE_SEP)
    out = []
    for x, y in zip(sa, sb):
        if x != y:
            break
        out.append(x)
    return gtdb.LINEAGE_SEP.join(out)


def key_of_tid(tmap, tid):
    for k, v in tmap.items():
        if v == tid:
            return k
    return gtdb.ROOT_TID if tid == gtdb.ROOT_TID else None


class TestWorkspaceGuard(unittest.TestCase):
    """UT-5.x: outputs cannot escape the workspace or hit the read-only input."""

    def setUp(self):
        self.root = os.environ.get(
            'NEOHGT_ROOT',
            os.path.dirname(os.path.dirname(os.path.dirname(
                os.path.abspath(__file__)))))

    def test_inside_ok(self):
        out = gtdb.guard_write(join(self.root, 'neoHGT2GTDB-缓存测试/x/y'),
                               root=self.root)
        self.assertTrue(out.endswith('/x/y'))

    def test_dotdot_escape_rejected(self):
        with self.assertRaises(gtdb.GuardError):
            gtdb.guard_write(join(self.root, '..', '..', 'etc', 'passwd'),
                             root=self.root)

    def test_absolute_escape_rejected(self):
        with self.assertRaises(gtdb.GuardError):
            gtdb.guard_write('/tmp/evil', root=self.root)

    def test_read_only_input_rejected(self):
        inp = join(self.root, 'neoHGT2GTDB-缓存测试')
        with self.assertRaises(gtdb.GuardError):
            gtdb.guard_write(join(inp, 'db', 'meta.json'), root=self.root,
                             must_not_contain=inp)

    def test_expanduser_env_var(self):
        os.environ['_NEOHGT_TEST_SUB'] = 'sub'
        self.assertEqual(gtdb.expand('/a/$_NEOHGT_TEST_SUB/b'), '/a/sub/b')

    def test_without_a_workspace_writes_are_not_contained(self):
        """No configured root must mean "anywhere", not a TypeError.

        README promises that outputs may live anywhere unless NEOHGT_ROOT is
        set; guard_write used to feed None to realpath() in exactly that case.
        """
        saved_env = {k: os.environ.pop(k, None)
                     for k in ('NEOHGT_ROOT', 'NEOHGT_CODE')}
        saved = gtdb.load_config_file
        gtdb.load_config_file = lambda: {'workspace': {}}
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        try:
            out = gtdb.guard_write(join(d, 'output'))
            self.assertEqual(out, os.path.realpath(join(d, 'output')))
            with self.assertRaises(gtdb.GuardError):
                gtdb.guard_write(join(d, 'input', 'x'), must_not_contain=join(d, 'input'))
        finally:
            gtdb.load_config_file = saved
            for k, v in saved_env.items():
                if v is not None:
                    os.environ[k] = v


if __name__ == '__main__':
    unittest.main()
