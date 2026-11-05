"""Unit tests for the GTDB-aware additions to neoHGT.util.

Covers UT-1.4, UT-1.5 and UT-1.6 of the GTDB adaptation plan. Existing
functions are exercised unchanged as well, to lock requirement RT-3 (no
behavior change for the NCBI path).
"""

import unittest

from neoHGT import util


class TestRankVocabulary(unittest.TestCase):
    """UT-1.6: rank word tables."""

    def test_gtdb_ranks(self):
        self.assertEqual(
            util.rank_order('gtdb'),
            ('domain', 'phylum', 'class', 'order', 'family', 'genus',
             'species'))
        self.assertNotIn('superkingdom', util.rank_order('gtdb'))

    def test_ncbi_ranks_unchanged(self):
        self.assertIn('superkingdom', util.rank_order('ncbi'))
        with self.assertRaises(ValueError):
            util.rank_order('nonsense')

    def test_rank_normalize_aliases(self):
        self.assertEqual(util.rank_normalize('superkingdom', 'gtdb'), 'domain')
        self.assertEqual(util.rank_normalize('kingdom', 'gtdb'), 'domain')
        self.assertEqual(util.rank_normalize('domain', 'ncbi'),
                         'superkingdom')
        self.assertEqual(util.rank_normalize('species_latin', 'gtdb'),
                         'species')

    def test_rank_prefix(self):
        self.assertEqual(util.rank_prefix('species', 'gtdb'), 's__')
        self.assertEqual(util.rank_prefix('domain', 'gtdb'), 'd__')
        self.assertEqual(util.rank_prefix('superkingdom', 'gtdb'), 'd__')
        self.assertEqual(util.rank_prefix('species', 'ncbi'), 's__')

    def test_rank_plural_label(self):
        self.assertEqual(util.rank_plural_label('domain', 'gtdb'), 'domains')
        self.assertEqual(util.rank_plural_label('genus', 'gtdb'), 'genera')

    def test_rank_plural_legacy_semantics_kept(self):
        # RT-3: the original function must keep its original behavior.
        self.assertEqual(util.rank_plural('superkingdom'), 'superkingdoms')
        self.assertEqual(util.rank_plural('phylum'), 'phyla')
        self.assertEqual(util.rank_plural('class'), 'classes')
        self.assertEqual(util.rank_plural('family'), 'families')
        self.assertEqual(util.rank_plural('genus'), 'genera')
        self.assertEqual(util.rank_plural('species'), 'species')


class TestStripRankPrefix(unittest.TestCase):
    """UT-1.4 / UT-1.5: prefix and cluster-letter handling."""

    def test_official_labels(self):
        self.assertEqual(util.strip_rank_prefix('s__Escherichia coli'),
                         'Escherichia coli')
        self.assertEqual(util.strip_rank_prefix('d__Bacteria'), 'Bacteria')

    def test_local_labels(self):
        self.assertEqual(util.strip_rank_prefix('g_Nanoarchaeum'),
                         'Nanoarchaeum')
        self.assertEqual(util.strip_rank_prefix('s_Escherichia_coli'),
                         'Escherichia coli')
        self.assertEqual(util.strip_rank_prefix('d_Archaea'), 'Archaea')

    def test_cluster_letter_at_species(self):
        self.assertEqual(util.strip_rank_prefix('s__Bacillus_A cereus'),
                         'Bacillus cereus')
        self.assertEqual(
            util.strip_rank_prefix('s__Enterobacter hormaechei_C'),
            'Enterobacter hormaechei')
        self.assertEqual(
            util.strip_rank_prefix('s__Methanobrevibacter_A smithii'),
            'Methanobrevibacter smithii')

    def test_cluster_letter_at_other_ranks(self):
        self.assertEqual(util.strip_rank_prefix('p__Bacillota_I'), 'Bacillota')
        self.assertEqual(util.strip_rank_prefix('c__Bacilli_A'), 'Bacilli')
        self.assertEqual(
            util.strip_rank_prefix('f__Enterobacteriaceae_A'),
            'Enterobacteriaceae')
        self.assertEqual(util.strip_rank_prefix('g__Synechococcus_C'),
                         'Synechococcus')

    def test_placeholder_names(self):
        self.assertEqual(util.strip_rank_prefix('g__DUMJ01'), 'DUMJ01')
        self.assertEqual(util.strip_rank_prefix('g_336-3'), '336-3')
        self.assertEqual(util.strip_rank_prefix('o__ML615J-28'), 'ML615J-28')
        self.assertEqual(util.strip_rank_prefix('f__DSM-6578'), 'DSM-6578')

    def test_keep_cluster_letter_when_asked(self):
        self.assertEqual(
            util.strip_rank_prefix('s__Bacillus_A cereus',
                                   strip_cluster=False),
            'Bacillus_A cereus')

    def test_no_prefix(self):
        self.assertEqual(util.strip_rank_prefix('Escherichia coli'),
                         'Escherichia coli')
        self.assertEqual(util.strip_rank_prefix(''), '')

    def test_capital_guards(self):
        # UT-1.4: this trio is what makes the "silent hit loss" bug impossible.
        self.assertTrue(util.is_capital(
            util.strip_rank_prefix('s__Escherichia coli')))
        self.assertFalse(util.is_capital('s__Escherichia coli'))
        self.assertFalse(util.is_capital('d_Archaea'))

    def test_latin_guards(self):
        # UT-1.5
        self.assertTrue(util.is_latin_gtdb('s__Bacillus_A cereus'))
        self.assertFalse(util.is_latin('Bacillus_A cereus'))
        self.assertTrue(util.is_latin_gtdb('s__Escherichia coli'))


class TestRankOfLabel(unittest.TestCase):
    """Rank inference from prefixed labels."""

    def test_prefixes(self):
        self.assertEqual(util.rank_of_label('f__Enterobacteriaceae'), 'family')
        self.assertEqual(util.rank_of_label('g_Nanoarchaeum'), 'genus')
        self.assertEqual(util.rank_of_label('d__Bacteria'), 'domain')
        self.assertIsNone(util.rank_of_label('Escherichia coli'))


class TestResolveTaxonRef(unittest.TestCase):
    """UT-3.5: all four accepted reference forms plus ambiguity."""

    def setUp(self):
        self.taxdump = {
            '1': {'parent': '1', 'rank': 'no rank', 'name': 'root'},
            '1000000000': {'parent': '1', 'rank': 'domain',
                           'name': 'Bacteria'},
            '1005000001': {'parent': '1000000000', 'rank': 'genus',
                           'name': 'Escherichia'},
            '1006000002': {'parent': '1005000001', 'rank': 'species',
                           'name': 'Escherichia coli'},
            '1005000009': {'parent': '1000000000', 'rank': 'genus',
                           'name': 'Synechococcus'},
            '1005000010': {'parent': '1000000000', 'rank': 'genus',
                           'name': 'Synechococcus'}}
        self.labels = {
            'd__Bacteria': '1000000000',
            'g__Escherichia': '1005000001',
            's__Escherichia coli': '1006000002',
            'g__Synechococcus': '1005000009',
            'g__Synechococcus_C': '1005000010'}

    def test_numeric(self):
        self.assertEqual(
            util.resolve_taxon_ref('1006000002', self.taxdump, self.labels),
            '1006000002')
        with self.assertRaises(ValueError):
            util.resolve_taxon_ref('999999', self.taxdump, self.labels)

    def test_prefixed_label(self):
        self.assertEqual(
            util.resolve_taxon_ref('g__Escherichia', self.taxdump,
                                   self.labels), '1005000001')
        self.assertEqual(
            util.resolve_taxon_ref('g_Escherichia', self.taxdump, self.labels),
            '1005000001')

    def test_unique_bare_name(self):
        self.assertEqual(
            util.resolve_taxon_ref('Escherichia coli', self.taxdump,
                                   self.labels), '1006000002')

    def test_ambiguous_name_lists_candidates(self):
        with self.assertRaises(ValueError) as ctx:
            util.resolve_taxon_ref('Synechococcus', self.taxdump,
                                   self.labels)
        msg = str(ctx.exception)
        self.assertIn('matches 2 taxIds', msg)
        self.assertIn('1005000009', msg)
        self.assertIn('1005000010', msg)

    def test_unknown_name(self):
        with self.assertRaises(ValueError):
            util.resolve_taxon_ref('NoSuchTaxon', self.taxdump, self.labels)

    def test_synonyms(self):
        self.assertEqual(
            util.resolve_taxon_ref('Proteobacteria', self.taxdump,
                                   self.labels,
                                   {'Proteobacteria': 'd__Bacteria'}),
            '1000000000')


class TestGtdbProduct(unittest.TestCase):
    """GTDB deflines carry no product name."""

    def test_returns_empty(self):
        self.assertEqual(util.gtdb_product(
            '# 1 # 879 # -1 # ID=1_1;partial=10;start_type=ATG'), '')


class TestPlaceholderNames(unittest.TestCase):
    """UT-2.6: GTDB coded labels are official names, not filter fodder."""

    def test_coded_group_names_are_placeholders(self):
        for name in ('g__DUMJ01', 'g_336-3', 'o__ML615J-28', 'f__DSM-6578',
                     'p__2-12-FULL-45-22', 'g__UBAR193', 'c__CPR1'):
            self.assertTrue(util.is_placeholder_gtdb(name), name)

    def test_candidate_species_are_placeholders(self):
        self.assertTrue(util.is_placeholder_gtdb(
            's__MK2026 sp010841923'))
        self.assertTrue(util.is_placeholder_gtdb('Ambinsecta sp002725225'))

    def test_real_names_are_not_placeholders(self):
        for name in ('s__Escherichia coli', 'g__Halobacterium',
                     'f__Enterobacteriaceae', 'd__Bacteria',
                     's__Bacillus_A cereus', '[Gaiellea] virgaureae', ''):
            self.assertFalse(util.is_placeholder_gtdb(name), name)


if __name__ == '__main__':
    unittest.main()
