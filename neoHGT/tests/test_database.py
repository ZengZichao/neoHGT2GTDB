"""
////////////////////////////////////////////////////////////
│	
│	Filename:	test_database.py
│	Comment:	See documentation
│	
│	=======================================================
│	Authorship:	Zhihui Li (Qiyun Zhu and Katharina Dittmar)
│	Copyright:	Distributed under the terms of the Modified BSD License.
│	
////////////////////////////////////////////////////////////
"""
#!/usr/bin/env python3

from unittest import TestCase, main
from os import remove, makedirs, listdir
from os.path import join, dirname, realpath, isfile, isdir
from shutil import rmtree, copy
from tempfile import mkdtemp
import gzip

import pandas as pd

from neoHGT.database import Database
from neoHGT.util import taxdump_from_text


class DatabaseTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datadir = join(dirname(realpath(__file__)), 'data')

        # whether to test remote functions, which highly depend on the network
        # connection and the current status of the NCBI server
        self.test_remote = False

    def tearDown(self):
        rmtree(self.tmpdir)

    def test___call__(self):
        # TODO
        pass

    def test_set_parameters(self):
        # TODO
        pass

    def test_connect_server(self):
        # TODO
        pass

    def test_retrieve_taxdump(self):
        # TODO
        pass

    def test_retrieve_summary(self):
        # TODO
        pass

    def test_retrieve_categories(self):
        # TODO
        pass

    def test_filter_genomes(self):
        me = Database()
        header = ('# assembly_accession', 'assembly_level', 'ftp_path')
        data = (('GCF_000000001.1', 'Chromosome', ''),
                ('GCF_000000002.1', 'Complete Genome', ''),
                ('GCF_000000003.2', 'Scaffold', ''),
                ('GCF_000000004.1', 'Contig', ''),
                ('GCA_000000004.1', 'Contig', ''))
        df = pd.DataFrame(data, columns=header)
        me.complete = False
        me.genoids = None
        me.exclude = False

        # drop duplicates
        me.df = df.copy()
        me.filter_genomes()
        self.assertEqual(me.df.shape[0], 4)
        self.assertListEqual(me.df['genome'].tolist(), [
            'G000000001', 'G000000002', 'G000000003', 'G000000004'])
        self.assertEqual(me.df.query(
            'accession == "GCF_000000004.1"').shape[0], 1)

        # complete genomes only
        me.complete = True
        me.df = df.copy()
        me.filter_genomes()
        self.assertListEqual(me.df['accnov'].tolist(), [
            'GCF_000000001', 'GCF_000000002'])

        # include certain genomes
        me.complete = False
        me.genoids = 'G000000001,G000000003'
        me.df = df.copy()
        me.filter_genomes()
        self.assertListEqual(me.df['accession'].tolist(), [
            'GCF_000000001.1', 'GCF_000000003.2'])

        # exclude certain genomes
        me.genoids = ['GCF_000000002.1', 'GCF_000000004']
        me.exclude = True
        me.df = df.copy()
        me.filter_genomes()
        self.assertListEqual(me.df['accession'].tolist(), [
            'GCF_000000001.1', 'GCF_000000003.2'])

    def test_sort_genomes(self):
        me = Database()
        header = (
            'genome', 'taxid', 'refseq_category', 'assembly_level',
            'relation_to_type_material')
        data = (
            ('G1', '585056', '', 'Chromosome', None),  # E. coli UMN026
            ('G2', '1038927', 'representative genome', 'Chromosome', None),
            # E. coli O104:H4
            ('G3', '2580236', '', 'Contig', None),  # sync E. coli
            ('G4', '622', '', 'Scaffold', 'yes'),  # Shigella
            ('G5', '548', '', 'Scaffold', None),  # Klebsiella
            ('G6', '126792', 'reference genome', 'Contig', None))  # plasmid
        df = pd.DataFrame(data, columns=header)
        me.taxdump = taxdump_from_text(taxdump_proteo)
        me.df = df.copy()
        me.sort_genomes()
        self.assertListEqual(me.df['genome'].tolist(), [
            'G6', 'G2', 'G4', 'G1', 'G5', 'G3'])

    def test_filter_by_taxonomy(self):
        me = Database()
        header = ('organism_name', 'taxid', 'species', 'species_taxid')
        data = (('Escherichia coli UMN026', '585056', 'E. coli', '562'),
                ('Escherichia coli O104:H4', '1038927', 'E. coli', '562'),
                ('Klebsiella aerogenes', '548', 'Klebsiella aerogenes', '548'),
                ('unclassified Gammaproteobacteria', '118884', '', ''),
                ('Plasmid pPY113', '126792', '', ''))
        df = pd.DataFrame(data, columns=header)

        # organism names must be capital and Latinate
        me.capital = True
        me.block = None
        me.latin = True
        me.taxids = None
        me.exclude = False
        me.taxdump = taxdump_from_text(taxdump_proteo)
        me.df = df.copy()
        me.filter_by_taxonomy()
        self.assertNotIn('species_taxid', me.df.columns)
        self.assertListEqual(me.df.index.tolist(), [0, 1, 2])
        self.assertListEqual(me.df['species'].tolist(), ['562', '562', '548'])

        # block word
        me.block = 'plasmid'
        me.latin = False
        me.df = df.copy()
        me.filter_by_taxonomy()
        self.assertListEqual(me.df.index.tolist(), [0, 1, 2])

        # no Escherichia
        me.taxids = '561'
        me.exclude = True
        me.df = df.copy()
        me.filter_by_taxonomy()
        self.assertListEqual(me.df.index.tolist(), [2])

    def test_sample_by_taxonomy(self):
        me = Database()

        # do nothing
        me.rank = None
        me.sample = None
        me.above = None
        self.assertIsNone(me.sample_by_taxonomy())

        # regular
        header = ('genome', 'taxid')
        data = (('G1',  '585056'),  # E. coli UMN026
                ('G2', '1038927'),  # E. coli O104:H4
                ('G3', '2580236'),  # sync E. coli
                ('G4',     '622'),  # Shigella
                ('G5',     '548'),  # Klebsiella
                ('G6',  '126792'))  # plasmid
        df = pd.DataFrame(data, columns=header)
        me.taxdump = taxdump_from_text(taxdump_proteo)

        # up to one genome per genus
        me.rank = 'genus'
        me.sample = 1
        me.df = df.copy()
        me.sample_by_taxonomy()
        self.assertListEqual(me.df.columns.tolist(), list(header) + ['genus'])
        self.assertListEqual(sorted(me.selected), ['G1', 'G4', 'G5'])

        # up to two genomes for entire cellular life
        me.rank = 'superkingdom'
        me.sample = 2
        me.sample_by_taxonomy()
        self.assertListEqual(sorted(me.selected), ['G1', 'G2'])

        # up to one genome per species
        me.rank = 'species'
        me.sample = 1
        me.sample_by_taxonomy()
        self.assertListEqual(sorted(me.selected), [
            'G1', 'G3', 'G4', 'G5', 'G6'])

        # up to one genome per species, Latinate only
        me.rank = 'species_latin'
        me.sample = 1
        me.df['latin'] = [True, True, True, True, True, False]
        me.sample_by_taxonomy()
        self.assertListEqual(sorted(me.selected), [
            'G1', 'G3', 'G4', 'G5'])

        # species and above (no difference)
        me.rank = 'species_latin'
        me.above = True
        me.sample = 10
        me.sample_by_taxonomy()
        self.assertListEqual(sorted(me.selected), [
            'G1', 'G2', 'G3', 'G4', 'G5'])

    def test_download_genomes(self):
        # TODO
        pass

    def test_extract_genomes(self):
        # TODO
        pass

    def test_genome_lineages(self):
        me = Database()
        me.output = self.tmpdir
        me.taxdump = taxdump_from_text(taxdump_proteo)
        data = (('G1', '1224', ''),     # Proteobacteria
                ('G2', '562',  '562'),  # Escherichia coli
                ('G3', '622',  '622'),  # Shigella dysenteriae
                ('G4', '548',  '548'))  # Klebsiella aerogenes
        me.df = pd.DataFrame(data, columns=[
            'genome', 'taxid', 'species']).set_index('genome')
        for rank in ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                     'family', 'genus']:
            me.df[rank] = ''
        me.genome_lineages()
        with open(join(self.tmpdir, 'lineages.txt'), 'r') as f:
            obs = dict(x.split('\t') for x in f.read().splitlines())
        proteo = 'k__Bacteria; p__Proteobacteria;'
        self.assertEqual(obs['G1'], proteo + ' c__; o__; f__; g__; s__')
        entero = proteo + ' c__Gammaproteobacteria; o__Enterobacterales;' +\
            ' f__Enterobacteriaceae;'
        self.assertEqual(
            obs['G2'], entero + ' g__Escherichia; s__Escherichia coli')
        self.assertEqual(
            obs['G3'], entero + ' g__Shigella; s__Shigella dysenteriae')
        self.assertEqual(
            obs['G4'], entero + ' g__Klebsiella; s__Klebsiella aerogenes')
        remove(join(self.tmpdir, 'lineages.txt'))

    def test_genome_metadata(self):
        me = Database()
        me.output = self.tmpdir
        me.df = pd.Series({
            'genome': 'G1',
            'accession': 'GCF_000123456.1',
            'asm_name': 'ASM123v1',
            'bioproject': 'PRJNA123456',
            'biosample': 'SAMN00123456',
            'assembly_level': 'Chromosome',
            'organism_name': 'hypothetical organism',
            'infraspecific_name': '',
            'isolate': '',
            'taxid': '12345',
            'ftp_path': ('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/'
                         '123/456/GCF_000123456.1_ASM123v1'),
            'proteins': 100,
            'residues': 12500,
            'whatever': 'nonsense'}).to_frame().T
        me.genome_metadata()
        with open(join(self.tmpdir, 'genomes.tsv'), 'r') as f:
            obs = f.read().splitlines()
        exp = ('genome', 'proteins', 'residues', 'assembly_level', 'accession',
               'bioproject', 'biosample', 'asm_name', 'organism_name',
               'infraspecific_name', 'isolate', 'taxid', 'ftp_path')
        self.assertEqual(obs[0], '\t'.join(exp))
        exp = ('G1', '100', '12500', 'Chromosome', 'GCF_000123456.1',
               'PRJNA123456', 'SAMN00123456', 'ASM123v1',
               'hypothetical organism', '', '', '12345',
               ('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/123/456/'
                'GCF_000123456.1_ASM123v1'))
        self.assertEqual(obs[1], '\t'.join(exp))
        remove(join(self.tmpdir, 'genomes.tsv'))

    def test_build_taxdump(self):
        me = Database()
        me.output = self.tmpdir
        me.tmpdir = join(self.datadir, 'DnaK', 'taxdump')
        me.taxdump = taxdump_from_text(taxdump_proteo)
        data = (('G1', '1224'),     # Proteobacteria
                ('G2', '562'),      # Escherichia coli
                ('G3', '585056'),   # E. coli UMN026
                ('G4', '1038927'))  # E. coli O104:H4
        me.df = pd.DataFrame(data, columns=[
            'genome', 'taxid']).set_index('genome')
        me.build_taxdump()
        with open(join(self.tmpdir, 'taxdump', 'nodes.dmp'), 'r') as f:
            obs = set(x.split('\t')[0] for x in f.read().splitlines())
        exp = {'1', '131567', '2', '1224', '1236', '91347', '543', '561',
               '562', '585056', '1038927'}
        self.assertSetEqual(obs, exp)
        rmtree(join(self.tmpdir, 'taxdump'))

    def test_build_taxonmap(self):
        me = Database()
        me.output = self.tmpdir
        me.taxdump = taxdump_from_text(taxdump_proteo)
        me.p2tids = {'P1': {'766'},  # Rickettsiales
                     'P2': {'570', '548'},  # Klebsiella
                     'P3': {'620', '622'},  # Shigella
                     'P4': {'561', '562'},  # Escherichia
                     'P5': {'126792', '28211'}}  # root
        me.build_taxonmap()
        exp = {'P1': '766', 'P2': '570', 'P3': '620', 'P4': '561', 'P5': '1'}
        self.assertDictEqual(me.taxonmap, exp)
        with gzip.open(join(self.tmpdir, 'taxon.map.gz'), 'rt') as f:
            obs = dict(x.split('\t') for x in f.read().splitlines())
        self.assertDictEqual(obs, exp)
        remove(join(self.tmpdir, 'taxon.map.gz'))

    def test_compile_database(self):
        me = Database()
        me.output = self.tmpdir

        # don't compile
        me.compile = 'none'
        me.compile_database()
        self.assertListEqual(listdir(self.tmpdir), [])

        # get database files
        copy(join(self.datadir, 'DnaK', 'linear.faa'),
             join(self.tmpdir, 'db.faa'))
        makedirs(join(self.tmpdir, 'taxdump'))
        copy(join(self.datadir, 'DnaK', 'taxdump', 'nodes.dmp'),
             join(self.tmpdir, 'taxdump', 'nodes.dmp'))
        copy(join(self.datadir, 'DnaK', 'taxdump', 'names.dmp'),
             join(self.tmpdir, 'taxdump', 'names.dmp'))
        with open(join(self.datadir, 'DnaK', 'prot2tid.txt'), 'r') as f:
            me.taxonmap = dict(x.split('\t') for x in f.read().splitlines())

        # set parameters
        me.threads = 1
        me.tmpdir = self.tmpdir
        me.makeblastdb = 'makeblastdb'
        me.diamond = 'diamond'

        # compile blast database
        me.compile = 'blast'
        me.compile_database()
        self.assertTrue(isdir(join(self.tmpdir, 'blast')))
        for ext in ('phr', 'pin', 'pog', 'psd', 'psi', 'psq'):
            self.assertTrue(isfile(join(self.tmpdir, 'blast', f'db.{ext}')))
        rmtree(join(self.tmpdir, 'blast'))

        # compile diamond database
        me.compile = 'diamond'
        me.compile_database()
        self.assertTrue(isdir(join(self.tmpdir, 'diamond')))
        self.assertTrue(isfile(join(self.tmpdir, 'diamond', 'db.dmnd')))
        rmtree(join(self.tmpdir, 'diamond'))

        # compile both databases
        me.compile = 'both'
        me.compile_database()
        self.assertTrue(isdir(join(self.tmpdir, 'blast')))
        for ext in ('phr', 'pin', 'pog', 'psd', 'psi', 'psq'):
            self.assertTrue(isfile(join(self.tmpdir, 'blast', f'db.{ext}')))
        self.assertTrue(isdir(join(self.tmpdir, 'diamond')))
        self.assertTrue(isfile(join(self.tmpdir, 'diamond', 'db.dmnd')))
        rmtree(join(self.tmpdir, 'blast'))
        rmtree(join(self.tmpdir, 'diamond'))

        # clean up
        remove(join(self.tmpdir, 'db.faa'))
        rmtree(join(self.tmpdir, 'taxdump'))

    def test_build_blast_db(self):
        me = Database()
        me.output = self.tmpdir
        me.makeblastdb = 'makeblastdb'
        me.tmpdir = self.tmpdir
        copy(join(self.datadir, 'DnaK', 'linear.faa'),
             join(self.tmpdir, 'db.faa'))
        with open(join(self.datadir, 'DnaK', 'prot2tid.txt'), 'r') as f:
            me.taxonmap = dict(x.split('\t') for x in f.read().splitlines())
        me.build_blast_db()
        self.assertTrue(isdir(join(self.tmpdir, 'blast')))
        for ext in ('phr', 'pin', 'pog', 'psd', 'psi', 'psq'):
            self.assertTrue(isfile(join(self.tmpdir, 'blast', f'db.{ext}')))
        rmtree(join(self.tmpdir, 'blast'))
        remove(join(self.tmpdir, 'db.faa'))

    def test_build_diamond_db(self):
        me = Database()
        me.output = self.tmpdir
        me.diamond = 'diamond'
        me.threads = 1
        me.tmpdir = self.tmpdir
        copy(join(self.datadir, 'DnaK', 'linear.faa'),
             join(self.tmpdir, 'db.faa'))
        with open(join(self.datadir, 'DnaK', 'prot2tid.txt'), 'r') as f:
            me.taxonmap = dict(x.split('\t') for x in f.read().splitlines())
        makedirs(join(self.tmpdir, 'taxdump'))
        copy(join(self.datadir, 'DnaK', 'taxdump', 'nodes.dmp'),
             join(self.tmpdir, 'taxdump', 'nodes.dmp'))
        copy(join(self.datadir, 'DnaK', 'taxdump', 'names.dmp'),
             join(self.tmpdir, 'taxdump', 'names.dmp'))
        me.build_diamond_db()
        self.assertTrue(isdir(join(self.tmpdir, 'diamond')))
        self.assertTrue(isfile(join(self.tmpdir, 'diamond', 'db.dmnd')))
        rmtree(join(self.tmpdir, 'diamond'))
        remove(join(self.tmpdir, 'db.faa'))
        remove(join(self.tmpdir, 'taxdump', 'nodes.dmp'))
        remove(join(self.tmpdir, 'taxdump', 'names.dmp'))

    def test_check_local_file(self):
        me = Database()

        # file does not exist
        file = join(self.tmpdir, 'tmp.in')
        self.assertFalse(me.check_local_file(file))

        # empty file will be deleted
        open(file, 'w').close()
        self.assertFalse(me.check_local_file(file))
        self.assertFalse(isfile(file))

        # file exists and has content
        with open(file, 'w') as f:
            f.write('Hello world!')
        self.assertTrue(isfile(file))
        self.assertTrue(me.check_local_file(file))

        # overwrite existing file
        self.assertFalse(me.check_local_file(file, overwrite=True))
        self.assertFalse(isfile(file))


"""Constants"""

class TestStreamExtraction(TestCase):
    """P7-5: --stream-extraction must reproduce the upstream extraction."""

    # accession -> (sequence lines, definition suffix); the same accession
    # appearing in two genomes is a shared (duplicate) protein
    PROTS = {
        'WP_000000001.1': (['MMMMMMMMMM', 'MMMMMMMMMM'], 'heat shock protein'),
        'WP_000000002.1': (['AAAAAAAAAA'], 'ribosomal protein'),
        'WP_000000003.1': (['QQQQQ'], 'sigma factor'),
        'WP_000000004.1': (['WW', 'WW'], 'cell wall protein'),
    }
    GENOMES = (  # genome, taxid, accession layout
        ('G1', '562', ('WP_000000001.1', 'WP_000000002.1')),
        ('G2', '622', ('WP_000000001.1', 'WP_000000003.1')),
        ('G3', '766', ('WP_000000002.1', 'WP_000000004.1')),
    )

    def setUp(self):
        from neoHGT.database import Database
        self.me = Database
        self.taxdump = taxdump_from_text(taxdump_proteo)
        self.base = mkdtemp()
        self.addCleanup(rmtree, self.base, ignore_errors=True)

    def _instance(self, name):
        from types import SimpleNamespace
        import pandas as pd
        root = join(self.base, name)
        faa = join(root, 'download', 'faa')
        makedirs(faa)
        stems = {g: f'GCF_0000000{i}.1_ASM{i}'
                 for i, (g, _, _) in enumerate(self.GENOMES, start=1)}
        for g, tid, accs in self.GENOMES:
            with gzip.open(join(faa, f'{stems[g]}_protein.faa.gz'), 'wt') as f:
                for a in accs:
                    seqs, desc = self.PROTS[a]
                    f.write(f'>{a} {desc}\n')
                    for s in seqs:
                        f.write(s + '\n')
        me = SimpleNamespace()
        me.output = root
        me.down = join(root, 'download')
        me.tmpdir = join(root, 'tmp')
        makedirs(me.tmpdir)
        me.taxdump = self.taxdump
        me.df = pd.DataFrame(
            [(g, tid, f'https://x/y/{stems[g]}')
             for g, tid, _ in self.GENOMES],
            columns=['genome', 'taxid', 'ftp_path'])
        return me

    def _read(self, me, fname):
        path = join(me.output, fname)
        opener = gzip.open if fname.endswith('.gz') else open
        with opener(path, 'rt') as f:
            return f.read()

    def test_outputs_are_identical_to_upstream(self):
        a, b = self._instance('a'), self._instance('b')
        self.me.extract_genomes(a)
        from neoHGT.database import extract_genomes_stream_ncbi
        extract_genomes_stream_ncbi(b)

        self.assertEqual(self._read(a, 'db.faa'), self._read(b, 'db.faa'))
        self.assertEqual(self._read(a, 'genome.map.gz'),
                         self._read(b, 'genome.map.gz'))
        pd.testing.assert_series_equal(a.df['proteins'], b.df['proteins'])
        pd.testing.assert_series_equal(a.df['residues'], b.df['residues'])

        # shared proteins are emitted once, but their full length is counted
        # for every genome that carries them (G2: 20 shared + 5 own)
        self.assertEqual(a.df.set_index('genome')['residues'].to_dict(),
                         {'G1': 30, 'G2': 25, 'G3': 14})
        self.assertEqual(a.df.set_index('genome')['proteins'].to_dict(),
                         {'G1': 2, 'G2': 2, 'G3': 2})

    def test_taxon_map_matches_build_taxonmap(self):
        a, b = self._instance('a2'), self._instance('b2')
        self.me.extract_genomes(a)
        self.me.build_taxonmap(a)
        from neoHGT.database import extract_genomes_stream_ncbi
        extract_genomes_stream_ncbi(b)
        self.assertEqual(self._read(a, 'taxon.map.gz'),
                         self._read(b, 'taxon.map.gz'))
        # LCA of shared proteins, as computed by the upstream dict version
        taxon = dict(x.split('\t') for x in
                     self._read(b, 'taxon.map.gz').splitlines())
        self.assertEqual(taxon['WP_000000001.1'], '543')   # 562 vs 622
        self.assertEqual(taxon['WP_000000002.1'], '1224')  # 562 vs 766
        self.assertEqual(taxon['WP_000000003.1'], '622')

    def test_compile_maps_match_upstream_temp_files(self):
        a, b = self._instance('a3'), self._instance('b3')
        self.me.extract_genomes(a)
        self.me.build_taxonmap(a)
        from neoHGT.database import extract_genomes_stream_ncbi
        extract_genomes_stream_ncbi(b)
        with open(join(a.tmpdir, 'prot.accession2taxid'), 'w') as f:
            f.write('accession\taccession.version\ttaxid\tgi\n')
            for p, tid in sorted(a.taxonmap.items()):
                f.write(f'{p.rsplit(".", 1)[0]}\t{p}\t{tid}\tna\n')
            blast_a = join(a.tmpdir, 'taxon.map')
            with open(blast_a, 'w') as g:
                for p, tid in sorted(a.taxonmap.items()):
                    g.write(f'{p}\t{tid}\n')
        self.assertEqual(open(join(a.tmpdir, 'prot.accession2taxid')).read(),
                         open(join(b.tmpdir, 'prot.accession2taxid')).read())
        self.assertEqual(open(blast_a).read(),
                         open(join(b.tmpdir, 'taxon.map')).read())

    def test_binding_swaps_methods(self):
        from neoHGT.database import (bind_stream_extraction,
                                     extract_genomes_stream_ncbi, arg2bool)
        me = self._instance('c')
        self.assertTrue(arg2bool('yes'))
        self.assertFalse(arg2bool(None))
        bind_stream_extraction(me)
        self.assertIs(me.extract_genomes.__func__,
                      extract_genomes_stream_ncbi)


taxdump_proteo = (
    '1,root,1,no rank',
    '131567,cellular organisms,1,no rank',
    '2,Bacteria,131567,superkingdom',
    '1224,Proteobacteria,2,phylum',
    '28211,Alphaproteobacteria,1224,class',
    '766,Rickettsiales,28211,order',
    '1236,Gammaproteobacteria,1224,class',
    '91347,Enterobacterales,1236,order',
    '543,Enterobacteriaceae,91347,family',
    '561,Escherichia,543,genus',
    '562,Escherichia coli,561,species',
    '585056,Escherichia coli UMN026,562,no rank',
    '1038927,Escherichia coli O104:H4,562,no rank',
    '2580236,synthetic Escherichia coli Syn61,561,species',
    '620,Shigella,543,genus',
    '622,Shigella dysenteriae,620,species',
    '570,Klebsiella,543,genus',
    '548,Klebsiella aerogenes,570,species',
    '118884,unclassified Gammaproteobacteria,1236,no rank',
    '126792,Plasmid pPY113,1,species')


if __name__ == '__main__':
    main()
