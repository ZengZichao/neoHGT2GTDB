"""
////////////////////////////////////////////////////////////
│	
│	Filename:	database.py
│	Comment:	See documentation
│	
│	=======================================================
│	Authorship:	Zhihui Li (Qiyun Zhu and Katharina Dittmar)
│	Copyright:	Distributed under the terms of the Modified BSD License.
│	
////////////////////////////////////////////////////////////
"""
#!/usr/bin/env python3

import sys
from os import remove, makedirs, cpu_count
from os.path import join, isfile, isdir, getsize, basename, exists
from shutil import which, rmtree
from time import sleep
import gzip
import tarfile
from tempfile import mkdtemp
import pandas as pd
#//////////////////////////////////////////////////
from datetime import datetime
from multiprocessing import Pool, cpu_count
#//////////////////////////////////////////////////
#import objgraph # for debugging memory leaks

from neoHGT.util import (
	timestamp, load_configs, get_config, arg2bool, run_command,
	list_from_param, read_taxdump, contain_words, is_latin, is_capital,
	taxid_at_rank, taxids_at_ranks, is_ancestral, find_lca, rank_plural)


description = """build reference protein sequence and taxonomy database"""

arguments = [
	'default protocol',
	['--default',    'apply the default protocol', {'action': 'store_true'}],

	'basic',
	['-o|--output',  'output database directory', {'required': True}],
	['-c|--cats',    'enter one or more of the following RefSeq genome '
					 'categories, separated by comma: archaea, bacteria, '
					 'fungi, invertebrate, plant, protozoa, vertebrate_'
					 'mammalian, vertebrate_other, viral; "all" for all '
					 'categories, or "microbe" for archaea, bacteria, fungi, '
					 'and protozoa', {'default': 'microbe'}],

	'custom list',
	['-t|--taxids',  'provide a custom list of TaxIDs; only genomes under '
					 'these taxonomic groups will be included'],
	['-g|--genoids', 'provide a custom list of genome IDs, which are NCBI '
					 'assembly accessions (e.g., "GCF_000123456.1"), with or '
					 'without version number, or simply "G000123456"'],
	['-e|--exclude', 'exclude instead of include TaxIDs and genomes defined '
					 'by custom lists', {'action': 'store_true'}],

	'taxon sampling',
	['-s|--sample',  'sample up to this number of genomes per taxonomic '
					 'group at given rank (0 for none, omit for all)',
					 {'type': int}],
	['-r|--rank',    'taxonomic rank at which sampling will be performed',
					 {'default': 'species'}],
	['--above',      'sampling will also be performed on ranks from the '
					 'designated one to phylum',
					 {'action': 'store_true'}],

	'genome sampling',
	['--genbank',    'also search GenBank if a genome is not found in RefSeq; '
					 'otherwise only include RefSeq genomes',
					 {'action': 'store_true'}],
	['--complete',   'only include complete genomes or chromosomes',
					 {'action': 'store_true'}],
	['--reference',  'include NCBI-defined reference genomes',
					 {'action': 'store_true'}],
	['--represent',  'include NCBI-defined representative genomes',
					 {'action': 'store_true'}],
	['--typemater',  'include NCBI-defined type material genomes',
					 {'action': 'store_true'}],

	'taxonomic filter',
	['--capital',    'organism name must be capitalized',
					 {'choices': ['yes', 'no']}],
	['--block',      'ignore organism names that contain any of these words'],
	['--latin',      'genomes must have Latinate species names, e.g. "Vibrio '
					 'cholerae", but not "Vibrio sp. 123"',
					 {'choices': ['yes', 'no']}],

	'download',
	['--manual',     'export URLs of sampled genomes without downloading them'
					 ', so the user may download them manually, then resume '
					 'the database pipeline.', {'action': 'store_true'}],

	['--overwrite',  'overwrite existing files with newly downloaded ones; '
					 'othewise use existing files whenever available, i.e., '
					 'resume an interrupted run', {'action': 'store_true'}],
	['--retries',    'number of trials for downloading each data file',
					 {'type': int}],
	['--delay',      'seconds between two retries', {'type': int}],
	['--timeout',    'seconds seconds before program gives up waiting',
					 {'type': int}],

	'compile',
	['--compile',    'compile database using this program if available',
					 {'choices': ['diamond', 'blast', 'both', 'none']}],
	['--diamond',    'diamond executable'],
	['--makeblastdb', 'makeblastdb executable'],
	['--threads',    'number of threads for diamond makedb', {'type': int}],
	['--tmpdir',     'temporary directory for diamond makedb']
]

server = 'rsync://ftp.ncbi.nlm.nih.gov'

# Arguments for the GTDB data source. They are appended as a separate
# statement so that the table-driven block above (and the NCBI behavior it
# encodes) stays untouched.
arguments += [
    'GTDB source',
    ['--source', 'data source for building the database (default: the '
     'config file\'s database.source, else ncbi). "hybrid" grafts an NCBI '
     'classification skeleton onto the GTDB tree but ingests no NCBI '
     'protein: the searchable sequences are the GTDB ones',
     {'choices': ['ncbi', 'gtdb', 'hybrid'], 'default': None}],
    ['--gtdb-dir', 'root directory of a local GTDB dataset, which is opened '
     'read-only; required when --seq-source is "local"'],
    ['--local-layout', 'layout of the local GTDB dataset; "auto" detects it '
     'from file names, sequence ids and compression ("gtdb-r220-renamed" is '
     'accepted as an alias of "local-renamed")',
     {'choices': ['auto', 'local-renamed', 'gtdb-r220-renamed',
                  'official-plain', 'official-gz']}],
    ['--gtdb-release', 'GTDB release label written to meta.json (inferred '
     'from file names when omitted)'],
    ['--gtdb-metadata', 'explicit GTDB metadata TSV path(s), comma-delimited'],
    ['--gtdb-taxonomy', 'explicit GTDB taxonomy TSV path(s), comma-delimited'],
    ['--gtdb-taxid-base', 'first taxId of the reserved GTDB range (default '
     '100000000; DIAMOND sizes db.dmnd by the highest taxId, so lower is '
     'much smaller, but it must stay above real NCBI taxIds)',
     {'type': int}],
    ['--gtdb-id-manifest', 'previous gtdb/taxid.map.gz, so that taxIds stay '
     'stable across releases'],
    ['--gtdb-no-selfcheck', 'skip the input self-consistency checks (not '
     'recommended: the build would not be reproducible)',
     {'action': 'store_true'}],
    ['--seq-source', 'where protein sequences are read from',
     {'choices': ['local', 'ncbi', 'gtdb-tar', 'ena']}],
    ['--allow-network', 'acknowledge that the remote sequence channels '
     '(ncbi, gtdb-tar, ena) were requested; this version has no download '
     'executor, so both values only export gtdb/fetch_plan.tsv for manual '
     'use and then refuse',
     {'choices': ['yes', 'no']}],
    ['--gtdb-mirror', 'GTDB data host used by the gtdb-tar channel',
     {'choices': ['ecogenomic', 'aau', 'auto']}],
    ['--downloader', 'downloader named in the exported fetch plan (this '
     'version performs no transfers itself)',
     {'choices': ['curl', 'wget']}],
    ['--verify-md5', 'md5_check column of the exported fetch plan (MD5 '
     'verification itself is not implemented in this version)',
     {'choices': ['yes', 'no']}],
    ['--workspace-out', 'override workspace.cache_dir for the output tree'],
    ['--workspace-tmp', 'override workspace.tmp_dir for intermediate files'],

    'GTDB genome selection',
    ['--reps', 'restrict to GTDB species-cluster representative genomes '
     '(equivalent to -r species -s 1)', {'action': 'store_true'}],
    ['--genome-list', 'file(s) of GTDB accessions to include, one per line; '
     'accepts the bundled *_Sample_<rank>.IDlist files'],
    ['--rank-sample', 'keep one genome per taxonomic group at this rank '
     '(cheap tier for development and CI)',
     {'choices': ['domain', 'phylum', 'class', 'order', 'family', 'genus']}],
    ['--gtdb-taxa', 'include only genomes under these GTDB taxa; accepts '
     'prefixed labels (f__Enterobacteriaceae) or GTDB taxIds'],
    ['--min-completeness', 'minimum Checkm2 completeness (%%) required',
     {'type': float}],
    ['--max-contamination', 'maximum Checkm2 contamination (%%) allowed',
     {'type': float}],
    ['--mimag', 'require a MIMAG quality tier',
     {'choices': ['high', 'medium', 'low', 'any']}],
    ['--type-material', 'which type-material definition to use',
     {'choices': ['gtdb', 'ncbi']}],
    ['--domain', 'limit the build to one GTDB domain',
     {'choices': ['bacteria', 'archaea', 'both'], 'default': 'both'}],
    ['--expand-members', 'grow a representative selection to every member '
     'genome of the species clusters it stands for (needs non-representative '
     'sequence files, which the local r220 copy does not contain)',
     {'choices': ['yes', 'no']}],
    ['--keep-placeholder-names', 'keep GTDB placeholder group names such as '
     'g__DUMJ01 (recommended, they are valid taxon labels)',
     {'choices': ['yes', 'no']}],
    ['--hybrid-ncbi-taxdump', 'directory with NCBI nodes.dmp/names.dmp used '
     'by --source hybrid. A hybrid database still contains GTDB proteins '
     'only: the NCBI dump grafts a classification skeleton, it adds no '
     'searchable sequences'],
    ['--hybrid-lca-policy', 'how to treat proteins shared between GTDB and '
     'NCBI subtrees',
     {'choices': ['root', 'domain', 'refuse'], 'default': 'root'}],

    'GTDB output',
    ['--keep-faa', 'what to do with db.faa once the DIAMOND database is built',
     {'choices': ['keep', 'compress', 'drop'], 'default': 'keep'}],
    ['--input-fingerprints', 'also record sha256 of the input files used '
     '(recommended; the local GTDB copy ships no MD5SUM.txt)',
     {'choices': ['yes', 'no']}],
]

#//////////////////////////////////////////////////
def download_parallel(formatted_input):
	try:
		download_command, genome = formatted_input
		ec, _ = run_command(download_command)
		if ec == 0:
			print('  ' + genome, flush=True)
		if ec > 0:
			print(f'WARNING: Cannot retrieve {genome}.')
			return genome
	except Exception as e:
		print(e)
#//////////////////////////////////////////////////

class Database(object):

	def __init__(self):
		self.arguments = arguments
		self.description = description

	def __call__(self, args):
		print(f'Database building started at {timestamp()}.')

		# load configurations
		self.cfg = load_configs()

		# read and validate arguments
		self.set_parameters(args)

		# retrieve taxonomy database
		self.retrieve_taxdump()

		# retrieve genome assembly summary
		self.retrieve_summary()

		# retrieve genome categories
		self.retrieve_categories()

		# filter genomes
		self.filter_genomes()

		# sort genomes by quality
		self.sort_genomes()

		# identify taxonomy of genomes
		self.filter_by_taxonomy()

		# sample genomes by taxonomy
		self.sample_by_taxonomy()

		# sample genomes by quality
		self.sample_by_quality()

		# filter genomes to sampled one
		self.filter_to_sampled()

		# download genomes
		self.download_genomes()

		# extract genomes
		self.extract_genomes()

		# identify genome lineages
		self.genome_lineages()

		# write genome metadata
		self.genome_metadata()

		# build taxonomy database
		self.build_taxdump()

		# build protein-to-taxonomy map
		self.build_taxonmap()

		# build protein sequence database
		self.compile_database()

		# clean up
		if hasattr(self, 'mkdtemp'):
			rmtree(self.tmpdir)

		print(f'Database building finished at {timestamp()}.')

	def set_parameters(self, args):
		"""Workflow for validating and setting arguments.

		Parameters
		----------
		args : dict
			command-line arguments
		"""
		# load arguments
		for key, val in vars(args).items():
			setattr(self, key, val)

		# load configurations
		for key in 'capital', 'block', 'latin':
			get_config(self, key, f'taxonomy.{key}')
		for key in 'retries', 'delay', 'timeout':
			get_config(self, key, f'download.{key}')
		for key in 'diamond', 'makeblastdb':
			get_config(self, key, f'program.{key}')
		for key in 'threads', 'tmpdir':
			get_config(self, key, f'local.{key}')

		# convert boolean values
		for key in 'capital', 'latin':
			setattr(self, key, arg2bool(getattr(self, key, None)))

		# make temporary directory
		if not self.tmpdir:
			self.tmpdir = mkdtemp()
			setattr(self, 'mkdtemp', True)  # mark for cleanup
		if not isdir(self.tmpdir):
			raise ValueError(f'Invalid temporary directory: {self.tmpdir}.')

		# check local executables
		for key, exe in {'blast': 'makeblastdb', 'diamond': 'diamond'}.items():
			if self.compile in (key, 'both'):
				if getattr(self, exe) is None:
					setattr(self, exe, exe)
				if which(getattr(self, exe)) is None:
					raise ValueError(
						f'Invalid {exe} executable: {getattr(self, exe)}.')

		# default protocol
		if self.default:
			print('The default protocol is selected for database building.')
			print('The program will download all protein sequences of NCBI '
				  'RefSeq genomes of bacteria, archaea, fungi and protozoa, '
				  'keep one genome per species, plus all NCBI-defined '
				  'reference and representative genomes.')
			self.cats = 'microbe'
			self.sample = 1
			self.rank = 'species'
			self.reference = True
			self.represent = True

			if self.diamond is None:
				self.diamond = 'diamond'
			if which(self.diamond) is None:
				print('WARNING: Cannot find DIAMOND in this computer. '
					  'You will need to manually compile the database '
					  'after download is complete.')
				self.compile = 'none'
			else:
				self.compile = 'diamond'

		# determine number of CPUs to use
		if self.compile in ('diamond', 'both') and not self.threads:
			self.threads = cpu_count()
			if self.threads is None:
				self.threads = 1

		# create target directories
		makedirs(self.output, exist_ok=True)
		self.down = join(self.output, 'download')
		makedirs(self.down, exist_ok=True)

	def retrieve_taxdump(self):
		"""Retrieve NCBI taxdump.
		"""
		fname = 'taxdump.tar.gz'
		rfile = f'pub/taxonomy/{fname}'
		lfile = join(self.down, fname)

		# download taxdump
		if not self.check_local_file(lfile, self.overwrite):
			print('Downloading NCBI taxonomy database...', end='', flush=True)
			run_command(f'rsync -Ltq {server}/{rfile} {self.down}')
			print(' done.')

		# read taxdump
		print('Reading NCBI taxonomy database...', end='', flush=True)
		with tarfile.open(lfile, 'r:gz') as f:
			f.extract('names.dmp', self.tmpdir)
			f.extract('nodes.dmp', self.tmpdir)
		self.taxdump = read_taxdump(self.tmpdir)
		print(' done.')
		print(f'  Total number of TaxIDs: {len(self.taxdump)}.')

	def retrieve_summary(self, genbank=False):
		"""Retrieve genome assembly summary.
		"""

		def get_summary(target):
			key = target.lower()
			fname = f'assembly_summary_{key}.txt'
			rfile = f'genomes/{key}/{fname}'
			lfile = join(self.down, fname)

			# download summary
			if not self.check_local_file(lfile, self.overwrite):
				print(f'Downloading {target} assembly summary...', end='',
					  flush=True)
				run_command(f'rsync -Ltq {server}/{rfile} {self.down}')
				print(' done.')

			# read summary
			print(f'Reading {target} assembly summary...', end='', flush=True)
			df = pd.read_table(lfile, skiprows=1, low_memory=False)
			# ===== bug fix =====
			col = df.columns[0]
			if col.startswith('#'):
				df = df.rename(columns={col: col.lstrip('# ')})
			print(' done.')
			return df

		self.df = get_summary('RefSeq')
		if self.genbank:
			self.df = pd.concat([self.df, get_summary('GenBank')])
		print(f'  Total number of genomes: {self.df.shape[0]}.')

	def retrieve_categories(self):
		"""Retrieve genome categories.
		"""
		# parse categories
		cats = self.cats
		if cats == 'all':
			return
		elif cats == 'microbe':
			self.cats = ['archaea', 'bacteria', 'fungi', 'protozoa']
		else:
			self.cats = cats.split(',')
		print(f'Genome categories: {", ".join(self.cats)}')

		def get_categories(target):
			key = target.lower()

			# validate categories
			ec, out = run_command(
				f'rsync --list-only --no-motd {server}/genomes/{key}/')
			#//////////////////////////////////////////////////
			print(f'Running Command: rsync --list-only --no-motd {server}/genomes/{key}/')
			#//////////////////////////////////////////////////
			cats = [line.rsplit(None, 1)[-1] for line in out if
					line.startswith('d') and not line.endswith('.')]
			for cat in self.cats:
				if cat not in cats:
					raise ValueError(
						f'😌 No worries. This error ("{cat}" is not a valid {target} genome category) occurs either when you typed a category with a wrong name (for example virus (invalid) instead of viral) or when the NCBI server has temporarily shuted-down the connection from you (in which case a cup of coffee ☕️ will help)')

			# get genome list per category
			print(f'Downloading genome list per {target} category...')
			makedirs(join(self.down, 'cats'), exist_ok=True)
			ldir = join(self.down, 'cats')
			fname = 'assembly_summary.txt'
			file_ = join(self.tmpdir, fname)

			asms = []
			for cat in self.cats:
				lfile = join(ldir, f'{key}_{cat}.txt')

				# use local file
				if self.check_local_file(lfile):
					with open(lfile, 'r') as f:
						asms_ = f.read().splitlines()

				# download remote file
				else:
					rfile = f'genomes/{key}/{cat}/{fname}'
					run_command(f'rsync -Ltq {server}/{rfile} {self.tmpdir}')
					with open(file_, 'r') as f:
						asms_ = [x.split('\t', 1)[0] for x in f.read().
								 splitlines() if not x.startswith('#')]
					with open(lfile, 'w') as f:
						f.write(''.join([x + '\n' for x in asms_]))

				print(f'  {cat}: {len(asms_)}')
				asms += asms_

			if isfile(file_):
				remove(file_)
			print('Done.')
			return asms

		# filter genomes by category
		asmset = set(get_categories('RefSeq'))
		if self.genbank:
			asmset.update(get_categories('GenBank'))
		self.df = self.df[self.df['assembly_accession'].isin(asmset)]
		print(f'  Total number of genomes in categories: {self.df.shape[0]}.')

	def filter_genomes(self):
		"""Filter genomes based on genome metadata.
		"""
		print('Filtering genomes...')
		n = self.df.shape[0]

		def report_diff(msg):
			nonlocal n
			n_ = self.df.shape[0]
			if n_ < n:
				print('  ' + msg.format(n - n_))
			n = n_

		# complete genomes only
		if self.complete:
			self.df = self.df[self.df['assembly_level'].isin(
				{'Complete Genome', 'Chromosome'})]
			report_diff('Dropped {} non-complete genomes.')

		# non-redundant genome IDs
		# typically not necessary, just in case
		self.df = self.df.rename(columns={'assembly_accession': 'accession'})
		#//////////////////////////////////////////////////
		# bug: split('.', 1) is not in the standard of modern Python
		self.df['accnov'] = self.df['accession'].str.split('.').str[0]
		# # bug: split('_', 1) is not in the standard of modern Python
		self.df['genome'] = 'G' + self.df['accnov'].str.split('_').str[-1]
		#//////////////////////////////////////////////////
		self.df = self.df.drop_duplicates(subset=['genome'])

		# include/exclude genome Ids
		if self.genoids:
			self.genoids = set(list_from_param(self.genoids))
			print(f'{"Ex" if self.exclude else "In"}cluding '
				  f'{len(self.genoids)} custom genome IDs...')
			self.df = self.df[(self.df['accession'].isin(self.genoids) |
							   self.df['accnov'].isin(self.genoids) |
							   self.df['genome'].isin(self.genoids)) !=
							  self.exclude]
			report_diff('Dropped {} genomes.')

			# further sampling / filtering will not occur
			if not self.exclude:
				self.capital = False
				self.latin = False
				self.block = ''
				self.sample = None

		# genomes without download link
		self.df = self.df.query('ftp_path != "na"')
		report_diff('Dropped {} genomes without download link.')
		print('Done.')

	def sort_genomes(self):
		"""Sort genomes by quality as informed by metadata.
		"""
		# sort by reference > representative > type material > other
		self.df['rc_seq'] = self.df.apply(
			lambda x: 0 if x['refseq_category'] == 'reference genome'
			else (1 if x['refseq_category'] == 'representative genome'
				  else (2 if pd.notnull(x['relation_to_type_material'])
				  else 3)), axis=1)

		# sort by complete > scaffold > contig
		self.df['al_seq'] = self.df['assembly_level'].map(
			{'Chromosome': 0, 'Complete Genome': 0, 'Scaffold': 1,
			 'Contig': 2})

		# sort genomes by three criteria
		self.df = self.df.sort_values(by=['rc_seq', 'al_seq', 'genome'])

		# clean up temporary columns
		self.df = self.df.drop(columns=['al_seq', 'rc_seq'])

	def filter_by_taxonomy(self):
		"""Filter genomes by taxonomy.
		"""
		print('Filtering genomes by taxonomy...')
		n = self.df.shape[0]

		def report_diff(msg):
			nonlocal n
			n_ = self.df.shape[0]
			if n_ < n:
				print('  ' + msg.format(n - n_))
			n = n_

		# remove non-capitalized organism names
		if self.capital:
			self.df = self.df[self.df['organism_name'].apply(is_capital)]
			report_diff(
				'Dropped {} genomes without capitalized organism name.')

		# block certain words in organism names
		if self.block:
			self.block = list_from_param(self.block)
			self.df = self.df[~self.df['organism_name'].apply(
				contain_words, args=(self.block,))]
			report_diff('Dropped {} genomes with one or more blocked words in '
						'organism name.')

		# remove original species information
		self.df = self.df.drop(columns=['species_taxid'])

		# drop genomes whose taxIds are not in taxdump
		self.df = self.df.dropna(subset=['taxid'])
		self.df['taxid'] = self.df['taxid'].astype(str)
		self.df = self.df[self.df['taxid'].isin(self.taxdump)]
		report_diff('Dropped {} genomes without valid taxId.')

		# assign genomes to species (represented by taxId not name)
		self.df['species'] = self.df['taxid'].apply(
			taxid_at_rank, rank='species', taxdump=self.taxdump)

		# drop genomes without species taxId
		self.df = self.df.dropna(subset=['species'])
		report_diff('Dropped {} genomes without valid species taxId.')

		# drop genomes without Latinate species name
		self.df['latin'] = self.df['species'].apply(
			lambda x: is_latin(self.taxdump[x]['name']))
		if self.latin:
			self.df = self.df.query('latin')
			report_diff('Dropped {} genomes without Latinate species name.')
		print('Done.')

		# include/exclude taxIds
		if self.taxids:
			self.taxids = set(list_from_param(self.taxids))
			print(f'{"Ex" if self.exclude else "In"}cluding '
				  f'{len(self.taxids)} custom TaxIDs...')

			self.df = self.df[self.df['taxid'].apply(
				lambda x: is_ancestral(x, self.taxids, self.taxdump)
				!= self.exclude)]
			report_diff('Dropped {} genomes.')

	def sample_by_taxonomy(self):
		"""Sample genomes at designated taxonomic rank.
		"""
		# don't sample; keep all
		if self.sample is None:
			self.selected = set(self.df['genome'])
			return

		print('Sampling genomes based on taxonomy...')
		self.selected, n = set(), 0
		rank, sample, latin = self.rank, self.sample, False

		# Latinate species names
		if rank == 'species_latin':
			rank, latin = 'species', True
		print(f'Up to {sample} genome(s) will be sampled per {rank}' +
			  (' (Latinate names only) ' if latin else '') + '.')

		# assign genomes to given rank
		if rank not in self.df.columns:
			self.df[rank] = self.df['taxid'].apply(
				taxid_at_rank, rank=rank, taxdump=self.taxdump)
		df_ = self.df.dropna(subset=[rank])

		# keep Latinate species names only
		if latin:
			df_ = df_.query('latin')

		# select up to given number of genomes of each taxonomic group
		self.selected = set(df_.groupby(rank).head(sample)['genome'])
		n = df_[rank].unique().shape[0]
		print(f'  Sampled {len(self.selected)} genomes from {n} '
			  f'{rank_plural(rank)}.')

		# sample at ranks above the current one
		if not self.above:
			return

		# determine ranks to sample at
		ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
		if rank not in ranks:
			raise ValueError(
				f'Cannot determine taxonomic ranks above "{rank}", because it '
				'is not among the six standard ranks: {", ".join(ranks)}.')
		ranks = ranks[:ranks.index(rank)][::-1]
		print(f'Sampling will also be performed on: {", ".join(ranks)}.')

		for r in ranks:
			self.df['tmp_rank'] = self.df['taxid'].apply(
				taxid_at_rank, rank=r, taxdump=self.taxdump)
			df_ = self.df.dropna(subset=['tmp_rank'])

			# calculate number of genomes yet to be sampled per group
			goal = pd.Series(sample, df_['tmp_rank'].unique())
			done = df_.query('genome in @self.selected')[
				'tmp_rank'].value_counts()
			left = goal.subtract(done, fill_value=0)
			left = left[left > 0].astype(int)

			# sample specific number of genomes per group
			s = []
			for tid, num in left.items():
				s.extend(df_.query(f'tmp_rank == "{tid}"').head(num)['genome'])
			print(f'  Sampled {len(s)} more genomes from {left.shape[0]} '
				  f'{rank_plural(r)}.')
			self.selected.update(s)

		print(f'  Sampled a total of {len(self.selected)} genomes at '
			  f'{rank} and above.')
		self.df = self.df.drop(columns=['tmp_rank'])

	def sample_by_quality(self):
		"""Sample genomes according to quality categories.
		"""
		# add reference genomes
		if self.reference:
			df_ = self.df.query('refseq_category == "reference genome"')
			print(f'Included {df_.shape[0]} reference genomes.')
			self.selected.update(df_['genome'])

		# add representative genomes
		if self.represent:
			df_ = self.df.query('refseq_category == "representative genome"')
			print(f'Included {df_.shape[0]} representative genomes.')
			self.selected.update(df_['genome'])

		# add type material genomes
		if self.typemater:
			df_ = self.df[self.df['relation_to_type_material'].notna()]
			print(f'Included {df_.shape[0]} type material genomes.')
			self.selected.update(df_['genome'])

	def filter_to_sampled(self):
		"""Filter genomes to sampled ones.
		"""
		n = len(self.selected)
		if n == 0:
			raise ValueError('No genome is retained after sampling.')
		self.df = self.df.query('genome in @self.selected')
		self.df = self.df.sort_values('genome')
		print(f'Total number of sampled genomes: {n}.')

	def download_genomes(self):
		"""Download genomes from NCBI.
		"""
		if self.manual:
			fname = 'urls.txt'
			self.df['ftp_path'].to_csv(
				join(self.output, 'urls.txt'), header=False, index=False)
			print(f'URLs of sampled genomes written to {fname}. You may '
				  'manually download their protein sequence data to faa/, '
				  'then restart the database building pipeline.')
			sys.exit(0)

		print('Downloading non-redundant genomic data from NCBI...',
			  flush=True)
		ldir = join(self.down, 'faa')
		makedirs(ldir, exist_ok=True)
		failed = []
		#//////////////////////////////////////////////////
		list_of_urls = []
		for row in self.df.itertuples():
			g = row.genome
			rdir = row.ftp_path.split('/', 3)[-1]
			stem = rdir.rsplit('/', 1)[-1]
			fname = f'{stem}_protein.faa.gz'
			lfile = join(ldir, fname)
			#//////////////////////////////////////////////////
			# downloaded files will be ignored
			for i in range(self.retries):
				download_command = f'rsync -Ltq {server}/{rdir}/{fname} {ldir}'
				temp_tuple = (download_command, g)
				if not self.check_local_file(lfile, overwrite=False):
					list_of_urls.append(temp_tuple)
		#//////////////////////////////////////////////////
		pool = Pool(cpu_count())
		print(f"There are {cpu_count()} CPU cores. Commencing parallel download sequence. Downloading {len(list_of_urls)} files.")
		results = pool.map(download_parallel, list_of_urls)
		pool.close()
		pool.join()
		with open(f'{ldir}/download_result_{datetime.now().strftime("%Y-%m-%d%H:%M:%S")}.txt', 'w') as fb:
			for line in results:
				if line is not None:
					fb.write(line + '\n')
		print('✅ Done.')

		# drop genomes that cannot be retrieved
		if len(failed):
			print('Failed to retrieve the following genomes:')
			print('  ' + ', '.join(failed))
			failed = set(failed)
			self.df = self.df.query('genome not in @failed')

	def extract_genomes(self):
		"""Extract proteins from genomes.

		Notes
		-----
		Write protein sequences of all genomes into db.faa.
		Write protein to genomes(s) map to genome.map.gz.
		"""
		print('Extracting downloaded genomic data...\n', end='', flush=True)
		ldir = join(self.down, 'faa')
		#/*====================*/
		prots = {} # memory leak warning
		#/*====================*/
		cp = None  # current protein accession
		fout = open(join(self.output, 'db.faa'), 'w')

		def write_prot():
			# some protein accessions may be duplicated
			try:
				fout.write(f'>{cp} {prots[cp]["name"]}\n{prots[cp]["seq"]}\n')
			except KeyError:
				return
			else:
				#/*====================*/
				tempseq = len(prots[cp]['seq'])
				tempgs = prots[cp]['gs']
				temptids = prots[cp]['tids']
				del prots[cp]
				#/*====================*/
				prots[cp] = {}
				prots[cp]['aa'] = tempseq
				prots[cp]['gs'] = tempgs
				prots[cp]['tids'] = temptids

		g2n, g2aa = {}, {}
		#//////////////////////////////////////////////////
		counter = 0
		#//////////////////////////////////////////////////
		for row in self.df.itertuples():
			g, tid = row.genome, row.taxid
			g2n[g], g2aa[g] = 0, 0
			stem = row.ftp_path.rsplit('/', 1)[-1]
			lfile = join(ldir, f'{stem}_protein.faa.gz')
			if exists(lfile):
				with gzip.open(lfile, 'rb') as f:
					try:
						content = f.read().decode().splitlines()
					except (OSError, EOFError, TypeError):
						print(f' skipping corrupted file {stem}.', end='',
							flush=True)
						continue
				cp = None
				for line in content:
					if line.startswith('>'):
						write_prot()
						p, name = line[1:].split(None, 1)
						g2n[g] += 1
						if p in prots:
							cp = None
							prots[p]['gs'].add(g)
							prots[p]['tids'].add(tid)
							g2aa[g] += prots[p]['aa']
						else:
							cp = p
							prots[p] = {
								'name': name, 'gs': {g}, 'tids': {tid}, 'seq': ''}
					elif cp:
						line = line.rstrip('*')
						prots[cp]['seq'] += line
						g2aa[g] += len(line)
				write_prot()
				#//////////////////////////////////////////////////
				counter += 1
				print(f'{counter}/{self.df.shape[0]}: Extracted {stem}.')
				#//////////////////////////////////////////////////
				#debugging memory leaks
				# objgraph.show_most_common_types(limit=5)
				#//////////////////////////////////////////////////
			else:
				print(f'❌ {lfile} does not exist.')
		fout.close()
		print('✅ Done.')
		print('Combined protein sequences written to db.faa.')

		print('☕️ Do not worry if you encounter a bug afterwards. It is likely that the program will run out of memory and invoke a kill command by the OS.')
		self.df['proteins'] = self.df['genome'].map(g2n)
		self.df['residues'] = self.df['genome'].map(g2aa)
		self.p2tids = {k: v['tids'] for k, v in prots.items()}

		# report summary
		print(f'  Total number of genomes extracted: {len(g2n)}.')
		print(f'  Total number of unique proteins extracted: {len(prots)}.')
		print('  Number of proteins shared by multiple genomes: {}.'.format(
			sum([1 for k, v in prots.items() if len(v['gs']) > 1])))
		print('  Number of proteins shared by multiple TaxIDs: {}.'.format(
			sum([1 for k, v in prots.items() if len(v['tids']) > 1])))

		# write protein to genomes(s) map
		fname = 'genome.map.gz'
		with gzip.open(join(self.output, fname), 'wb') as f:
			for k, v in sorted(prots.items()):
				f.write('{}\t{}\n'.format(k, ','.join(
					sorted(v['gs']))).encode())
		print(f'Protein-to-genome(s) map written to {fname}.')

	def genome_lineages(self):
		"""Generate lineage information for genomes.

		Notes
		-----
		Write genome lineages to lineages.txt.
		"""
		# identify taxa at standard ranks
		ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
				 'family', 'genus', 'species']
		self.df[ranks[:-1]] = self.df['taxid'].apply(lambda x: pd.Series(
			taxids_at_ranks(x, ranks[:-1], self.taxdump)))

		# report number of taxa represented at each rank
		print('Number of taxonomic groups represented:')
		for rank in ranks:
			print(f'  {rank}: {self.df[rank].nunique()}.')

		# missing ranks are absent, not a value; normalise them before building
		# the lineage so None/NaN and '' behave the same on pandas 2 and 3
		self.df[ranks] = self.df[ranks].astype(object).fillna('')

		# merge superkingdom and kingdom
		self.df['kingdom'] = self.df[['superkingdom', 'kingdom']].apply(
			lambda x: x.iloc[1] if x.iloc[1] else x.iloc[0], axis=1)

		# generate lineage string
		self.df['lineage'] = self.df[ranks[1:]].fillna('').apply(
			lambda col: col.apply(lambda val, name: '{}__{}'.format(
				name[0], self.taxdump[val]['name'] if val else ''),
				args=(col.name,))).apply('; '.join, axis=1)

		# write table
		fname = 'lineages.txt'
		self.df['lineage'].to_csv(
			join(self.output, fname), sep='\t', header=False)
		print(f'Genome lineages written to {fname}.')

	def genome_metadata(self):
		"""Write and report genome metadata.

		Notes
		-----
		Write genome metadata to genomes.tsv.
		"""
		self.df = self.df.set_index('genome')
		self.df = self.df[[
			'proteins', 'residues', 'assembly_level',
			'accession', 'bioproject', 'biosample', 'asm_name',
			'organism_name', 'infraspecific_name', 'isolate', 'taxid',
			'ftp_path']]
		fname = 'genomes.tsv'
		self.df.to_csv(join(self.output, fname), sep='\t')
		print(f'Genome metadata written to {fname}.')

	def build_taxdump(self):
		"""Build taxonomy database.

		Notes
		-----
		Write refined taxonomy database to taxdump/nodes.dmp and names.dmp.
		"""
		tids = self.df['taxid'].unique().tolist()

		# identify all ancestral TaxIDs
		ancs = set(tids)
		for tid in tids:
			cid = tid
			while True:
				pid = self.taxdump[cid]['parent']
				if cid == pid or pid in ancs:
					break
				ancs.add(pid)
				cid = pid

		# shrink taxdump files
		dir_ = join(self.output, 'taxdump')
		makedirs(dir_, exist_ok=True)
		for fname in ('nodes.dmp', 'names.dmp'):
			fo = open(join(dir_, fname), 'w')
			fi = open(join(self.tmpdir, fname), 'r')
			for line in fi:
				row = line.rstrip('\r\n').replace('\t|', '').split('\t')
				if row[0] in ancs:
					if fname == 'nodes.dmp' or 'scientific name' in row:
						fo.write(line)
			fi.close()
			fo.close()
		print('Taxonomy database written to taxdump/.')

	def build_taxonmap(self):
		"""Build protein-to-taxonomy map.
		"""
		# assign shared protein to lowest common ancestor (LCA)
		self.taxonmap = {p: max(tids) if len(tids) == 1 else find_lca(
			tids, self.taxdump) for p, tids in self.p2tids.items()}

		# write taxonomy map
		fname = 'taxon.map.gz'
		with gzip.open(join(self.output, fname), 'wb') as f:
			for p, tid in sorted(self.taxonmap.items()):
				f.write(f'{p}\t{tid}\n'.encode())
		print(f'Protein-to-taxonomy map written to {fname}.')

	def compile_database(self):
		"""Compile database using external programs.
		"""
		if self.compile == 'none':
			return
		if self.compile in ('blast', 'both'):
			self.build_blast_db()
		if self.compile in ('diamond', 'both'):
			self.build_diamond_db()

	def build_blast_db(self):
		"""Build BLAST database.
		"""
		# create temporary taxon map
		taxonmap = join(self.tmpdir, 'taxon.map')
		with open(taxonmap, 'w') as f:
			for p, tid in sorted(self.taxonmap.items()):
				f.write(f'{p}\t{tid}\n')

		# build BLAST database
		makedirs(join(self.output, 'blast'), exist_ok=True)
		cmd = ' '.join((
			self.makeblastdb,
			'-dbtype', 'prot',
			'-in', join(self.output, 'db.faa'),
			'-out', join(self.output, 'blast', 'db'),
			'-title', 'db',
			'-parse_seqids',
			'-taxid_map', taxonmap))
		print('Build BLAST database...', flush=True)
		ec, out = run_command(cmd)
		if ec:
			raise ValueError(f'makeblastdb failed with error code {ec}.')

		# clean up
		remove(taxonmap)
		print('Done.')

	def build_diamond_db(self):
		"""Build DIAMOND database.
		"""
		# create temporary taxon map
		taxonmap = join(self.tmpdir, 'prot.accession2taxid')
		with open(taxonmap, 'w') as f:
			f.write('accession\taccession.version\ttaxid\tgi\n')
			for p, tid in sorted(self.taxonmap.items()):
				f.write(f'{p.rsplit(".", 1)[0]}\t{p}\t{tid}\tna\n')

		# build DIAMOND database
		makedirs(join(self.output, 'diamond'), exist_ok=True)
		cmd = ' '.join((
			self.diamond, 'makedb',
			'--threads', str(self.threads),
			'--in', join(self.output, 'db.faa'),
			'--taxonmap', taxonmap,
			'--taxonnodes', join(self.output, 'taxdump', 'nodes.dmp'),
			'--taxonnames', join(self.output, 'taxdump', 'names.dmp'),
			'--db', join(self.output, 'diamond', 'db'),
			'--tmpdir', self.tmpdir))
		print('Build DIAMOND database...', flush=True)
		ec, out = run_command(cmd)
		if ec:
			raise ValueError(f'diamond failed with error code {ec}.')

		# clean up
		remove(taxonmap)
		print('Done.')

	@staticmethod
	def check_local_file(file, overwrite=False):
		"""Check existing local file.

		Parameters
		----------
		file : str
			expected local file path and name
		overwrite : bool, optional
			whether overwrite existing file (default: False)

		Returns
		-------
		bool
			whether file exists
		"""
		if isfile(file):
			if getsize(file) == 0:
				remove(file)
			elif overwrite:
				print(f'  WARNING: existing local file {basename(file)} will be overwritten.')
				remove(file)
			else:
				#print(f'  Using local file {basename(file)}.')
				return True
		return False


# ===========================================================================
# GTDB / hybrid data source support
#
# Everything below is appended at module level and reaches the class only
# through the two assignments at the bottom of this file. The NCBI pipeline
# above is left byte-for-byte intact, which is what makes the zero-regression
# requirement (plan G4 / DoD-7 / RT-1) verifiable. Scope note: "byte-for-byte"
# covers the NCBI *logic*; the shared argparse table below it gains the GTDB
# options, so `--help` lists options the NCBI branch does not read (see the
# --force help).
# ===========================================================================

import os as _os
import re as _re
import json as _json
import shutil as _shutil
import subprocess as _subprocess
from os.path import realpath as _realpath, dirname as _dirname
from neoHGT import gtdb as _gtdb
from neoHGT.util import taxids_at_ranks as _taxids_at_ranks

gtdb = _gtdb  # readable name for the helpers below
dirname = _dirname  # the NCBI code above never used it; the GTDB helpers do

# Bytes of db.faa per protein, measured on three real builds of the local r220
# dataset (14,373 / 905,829 / 3,247,172 proteins -> 455 / 417 / 466 B per
# protein). The plan's 406 GB figure covers the whole per-genome tree
# (373k genomes, not only representatives), which is why it overshoots by ~2.8x.
GTDB_BYTES_PER_PROTEIN = 466.0
# Peak multiple of db.faa size needed during 'diamond makedb'
# (db.faa + .dmnd + --tmpdir + taxonmap), per plan §11.2 item 2.
GTDB_PEAK_FACTOR = 3.0

GTDB_EXTRA_ARGS = [
    'GTDB safety',
    ['--force', 'proceed even when the disk-space pre-flight fails (read by '
     'the GTDB route: it also unlocks the experimental BLAST compile and '
     'publishes builds with missing genomes; the NCBI route ignores it)',
     {'action': 'store_true'}],
    ['--max-genomes', 'build at most this many genomes (protects the shared '
     'disk during development)', {'type': int}],
    ['--limit', 'alias of --max-genomes', {'type': int}],
    ['--taxdump-trim', 'write only the taxa ancestral to the selected genomes '
     'instead of the full GTDB tree (smaller files, but taxIds are then '
     'specific to this selection)', {'action': 'store_true'}],
]
arguments += GTDB_EXTRA_ARGS


def set_source_parameters(self, args):
    """Validate and normalize the options of a non-NCBI build.

    Parameters
    ----------
    self : Database
        instance whose attributes are set
    args : argparse.Namespace
        command-line arguments
    """
    self.source = (getattr(args, 'source', None)
                   or ((self.cfg or {}).get('database') or {}).get('source')
                   or 'ncbi')
    # config defaults for the GTDB section (plan §6.1); the command line wins
    gcfg = ((self.cfg or {}).get('database') or {}).get('gtdb') or {}
    self.gtdb_dir = (getattr(args, 'gtdb_dir', None) or gcfg.get('dir')
                     or _os.environ.get('GTDB_DIR'))
    self.seq_source = (getattr(args, 'seq_source', None)
                       or gcfg.get('seq_source') or 'local')
    self.local_layout = (getattr(args, 'local_layout', None)
                         or gcfg.get('local_layout') or 'auto')
    release = getattr(args, 'gtdb_release', None) or gcfg.get('release')
    self.gtdb_release = None if release in (None, 'auto') else str(release)
    # explicit table paths (plan §7.1); a plain string is also accepted so that
    # the CLI can take "a.tsv,b.tsv"
    self.gtdb_metadata = (list_from_param(getattr(args, 'gtdb_metadata', None))
                          if getattr(args, 'gtdb_metadata', None) else None)
    self.gtdb_taxonomy = (list_from_param(getattr(args, 'gtdb_taxonomy', None))
                          if getattr(args, 'gtdb_taxonomy', None) else None)
    self.id_base = (getattr(args, 'gtdb_taxid_base', None)
                    or gcfg.get('taxid_base') or _gtdb.DEFAULT_ID_BASE)
    self.gtdb_mirror = (getattr(args, 'gtdb_mirror', None)
                        or gcfg.get('mirror') or 'ecogenomic')
    # remote-channel knobs: they only shape the exported fetch plan (the actual
    # transfer needs --allow-network yes), but they must be honoured, not
    # silently ignored (plan §7.1)
    self.downloader = getattr(args, 'downloader', None) or 'curl'
    self.verify_md5 = arg2bool(getattr(args, 'verify_md5', None)) is not False
    _gtdb.check_id_base(self.id_base)
    self.allow_network = arg2bool(getattr(args, 'allow_network', None)) is True
    self.keep_faa = getattr(args, 'keep_faa', None) or 'keep'
    self.max_genomes = (getattr(args, 'max_genomes', None)
                        or getattr(args, 'limit', None))
    if getattr(args, 'max_genomes', None) and getattr(args, 'limit', None) \
            and str(args.limit) != str(args.max_genomes):
        print(f'  NOTE: --limit {args.limit} is ignored; --max-genomes '
              f'{args.max_genomes} wins.')
    self.force = bool(getattr(args, 'force', False))
    self.fingerprints = arg2bool(
        getattr(args, 'input_fingerprints', None)
        if getattr(args, 'input_fingerprints', None) is not None
        else gcfg.get('verify_input_fingerprints'))
    # quality thresholds stay opt-in: applying them silently would make
    # "--reps" counts disagree with the metadata representative totals (DoD-2)
    self.no_selfcheck = arg2bool(getattr(args, 'gtdb_no_selfcheck', None)) or \
        gcfg.get('selfcheck') == 'skip'

    # workspace defaults and boundaries (plan §5.8)
    cfg = self.cfg or {}
    ws = cfg.get('workspace') or {}
    self.cache_dir = (getattr(args, 'workspace_out', None)
                      or ws.get('cache_dir')
                      or _os.environ.get('NEOHGT_CACHE'))
    self.ws_tmp = (getattr(args, 'workspace_tmp', None) or ws.get('tmp_dir')
                   or _os.environ.get('NEOHGT_TMP'))
    self.ws_root = _os.environ.get('NEOHGT_ROOT') or (
        dirname(_realpath(self.cache_dir)) if self.cache_dir else None)

    # Plan §5.8 / UT-5.x: the scratch directory must live inside the workspace.
    # Upstream defaults an omitted --tmpdir to mkdtemp() in the *system* temp
    # dir, which only lands inside the workspace by accident when tools/env.sh
    # exports TMPDIR; anyone running the CLI without that wrapper would spread
    # multi-GB intermediates over /tmp. Relocate it here -- this function runs
    # only on the gtdb/hybrid branch, so the NCBI route keeps its upstream
    # default byte-for-byte.
    if not getattr(args, 'tmpdir', None) and self.ws_tmp:
        _gtdb.guard_write(self.ws_tmp, root=self.ws_root)
        _os.makedirs(self.ws_tmp, exist_ok=True)
        old = getattr(self, 'tmpdir', None)
        self.tmpdir = mkdtemp(prefix='gtdb_', dir=self.ws_tmp)
        if getattr(self, 'mkdtemp', False) and old and not _os.listdir(old):
            _os.rmdir(old)      # only ever removes the empty dir we just made
        print(f'  Scratch directory inside the workspace: {self.tmpdir}')

    if self.source in ('gtdb', 'hybrid'):
        if self.seq_source == 'local' and not self.gtdb_dir:
            raise ValueError(
                '--source gtdb with --seq-source local requires --gtdb-dir '
                '(the read-only GTDB dataset root).')
        if not self.gtdb_dir:
            raise ValueError('--gtdb-dir is required for --source gtdb/hybrid.')
        if not isdir(self.gtdb_dir):
            raise ValueError(
                f'GTDB dataset directory does not exist: {self.gtdb_dir}.')
        if _re.search(r'\b(RS_GCF_|GB_GCA_)\S+', str(self.gtdb_taxa or '')):
            raise ValueError(
                '--gtdb-taxa takes taxon labels or taxIds, not genome '
                f'accessions; use --genome-list for "{self.gtdb_taxa}".')

    # hybrid mode needs an NCBI taxdump of its own: it is not part of the local
    # GTDB copy (plan §3.4), so ask for it before any work happens (§5.7)
    self.hybrid_lca_policy = (getattr(args, 'hybrid_lca_policy', None)
                              or 'root')
    self.hybrid_ncbi_taxdump = (getattr(args, 'hybrid_ncbi_taxdump', None)
                                or (self.cfg or {}).get(
                                    'database', {}).get(
                                        'hybrid', {}).get('ncbi_taxdump'))
    if self.source == 'hybrid':
        if not self.hybrid_ncbi_taxdump:
            raise ValueError(
                '--source hybrid needs --hybrid-ncbi-taxdump <dir with '
                'nodes.dmp/names.dmp>: the local GTDB copy ships no NCBI '
                'taxonomy (plan §3.4). Use --source gtdb for a pure GTDB '
                'database.')
        if not isdir(self.hybrid_ncbi_taxdump):
            raise ValueError(
                f'--hybrid-ncbi-taxdump is not a directory: '
                f'{self.hybrid_ncbi_taxdump}.')
        if self.hybrid_lca_policy not in _gtdb.HYBRID_LCA_POLICIES:
            raise ValueError(
                f'Unknown --hybrid-lca-policy "{self.hybrid_lca_policy}": '
                + ', '.join(_gtdb.HYBRID_LCA_POLICIES) + '.')

    # output and temporary directories must stay inside the workspace and
    # must never land on the read-only input dataset
    if self.gtdb_dir:
        if _gtdb.is_within(_realpath(self.output), _realpath(self.gtdb_dir)):
            raise ValueError(
                f'-o {self.output} is inside the read-only GTDB dataset '
                f'{self.gtdb_dir}; choose an output under the workspace.')
    _gtdb.guard_write(self.output, root=self.ws_root,
                      must_not_contain=self.gtdb_dir)
    if self.ws_tmp:
        _gtdb.guard_write(self.ws_tmp, root=self.ws_root)

    # GTDB profile: names carry rank prefixes and placeholder labels, so the
    # capitalization / Latinity / block-word filters would silently discard
    # hits (plan R4, RK4). They stay off unless the user asks otherwise.
    if self.source in ('gtdb', 'hybrid'):
        if getattr(args, 'capital', None) is None:
            self.capital = False
        if getattr(args, 'latin', None) is None:
            self.latin = False
        if getattr(args, 'block', None) in (None, ''):
            self.block = []
        if getattr(args, 'keep_placeholder_names', None) is None:
            self.keep_placeholder_names = True
        else:
            self.keep_placeholder_names = arg2bool(
                args.keep_placeholder_names)
            # the build keeps every label it is given; this filter is applied
            # while searching, so saying nothing here would leave the user
            # believing the database had been pruned
            print('NOTE: --keep-placeholder-names has no effect on the '
                  'database; placeholder labels are filtered at search time '
                  '(--tax-latin / --tax-capital, or '
                  'taxonomy.profiles.gtdb.keep_placeholder in config.yml).')
    return self


def make_provider(self, args):
    """Build the source provider for this run (GTDB, or GTDB + NCBI).

    Parameters
    ----------
    self : Database
        instance holding the validated source options
    args : argparse.Namespace
        parsed arguments

    Returns
    -------
    SourceProvider
        a prepared provider; for ``--source hybrid`` the GTDB subtree is
        grafted under the NCBI cellular root (plan §5.7)
    """
    manifest = getattr(args, 'gtdb_id_manifest', None) or join(
        self.output, 'gtdb', 'taxid.map.gz')
    # remembered for write_gtdb_outputs(): the reserved-ids sidecar next to
    # this manifest is what keeps historical taxIds reserved across rebuilds
    self.gtdb_manifest = manifest
    provider = _gtdb.GtdbProvider(
        self.gtdb_dir, release=self.gtdb_release, layout=self.local_layout,
        id_base=self.id_base, manifest_path=manifest,
        selfcheck=not self.no_selfcheck,
        metadata_files=self.gtdb_metadata, taxonomy_files=self.gtdb_taxonomy)
    provider.prepare()
    if self.source == 'hybrid':
        # without this branch --source hybrid would silently build a pure GTDB
        # database and ignore --hybrid-ncbi-taxdump / --hybrid-lca-policy
        provider = _gtdb.HybridProvider(
            provider, self.hybrid_ncbi_taxdump,
            lca_policy=self.hybrid_lca_policy)
        provider.prepare()
        print(f'  hybrid: NCBI taxdump {self.hybrid_ncbi_taxdump}, '
              f'lca policy "{self.hybrid_lca_policy}"')
    return provider


def build_from_source(self, args):
    """Build a reference database from GTDB (or a GTDB + NCBI hybrid).

    Parameters
    ----------
    self : Database
        instance to populate
    args : argparse.Namespace
        command-line arguments

    Notes
    ----
    The steps mirror the NCBI pipeline so that search/analyze need no changes:
    db.faa, genomes.tsv, lineages.txt, taxon.map.gz, genome.map.gz, taxdump/,
    diamond/ and -- on top of the NCBI layout -- gtdb/ and meta.json.
    """
    started = timestamp()
    print(f'Database building started at {started} (source: {self.source}).')

    # 1. probe the dataset and prove it is self-consistent (plan §3.6, ADR-9)
    print(f'Probing GTDB dataset at {self.gtdb_dir}...')
    provider = make_provider(self, args)
    print(f'  layout: {provider.info["layout"]}, release: '
          f'{provider.info["release"]}, sequence Ids: '
          f'{provider.info["id_convention"]}')

    # 2. full taxonomy (needed to resolve labels and taxIds)
    print('Converting GTDB lineages into a taxdump...')
    taxdump_full, taxid_map, stats = provider.load_taxonomy()
    print(f'  {stats["n_nodes"]:,} nodes (GTDB '
          f'{stats["n_gtdb_nodes"]:,} + root), taxId range '
          f'{stats["id_min"]}-{stats["id_max"]}')

    # 3. select and sample genomes
    # -e/--exclude is rejected before any metadata/taxdump work: failing after
    # the full-table load made the refusal cost a minute on real releases.
    if getattr(args, 'exclude', False):
        raise ValueError(
            '-e/--exclude is not supported for --source gtdb yet; select the '
            'genomes to keep with --gtdb-taxa/--genome-list instead.')
    provider.opts = gtdb.selection_opts(args)
    df = provider.select(provider.load_genomes())
    n_before = df.shape[0]
    df = provider.sample(df)
    print(f'Selected {df.shape[0]:,} of {n_before:,} genomes'
          + (f' ({", ".join(provider.selection_labels)})'
             if provider.selection_labels else ''))
    if getattr(provider, 'sample_note', None):
        print(f'  sampling: {provider.sample_note}')
    if df.empty:
        raise ValueError(
            'The selection criteria matched no GTDB genome; check --gtdb-taxa, '
            '--rank-sample and the quality thresholds.')
    if self.max_genomes and df.shape[0] > self.max_genomes:
        df = df.head(int(self.max_genomes))
        print(f'  LIMITED to the first {self.max_genomes} selected genomes '
              '(--max-genomes).')
    self.df = df

    # 3b. optional expansion to the member genomes of the selected clusters
    if arg2bool(getattr(args, 'expand_members', None)):
        self.df = df = expand_selection(self, provider, df)

    # 4. disk pre-flight (plan §11.2, RK18)
    self.capacity = check_capacity(self, df)

    # 5. optionally trim the tree to the selected genomes. The default keeps
    # the full GTDB tree so that taxIds stay meaningful across subset builds
    # (plan §6.4 task 2.10).
    if getattr(args, 'taxdump_trim', False):
        keep = set()
        for lin in df['gtdb_taxonomy']:
            keep.update(gtdb.lineage_paths(lin))
        taxdump, taxid_map, stats = provider.load_taxonomy(trim=keep)
        self.taxdump, self.taxid_map, self.stats = taxdump, taxid_map, stats
        print(f'  taxdump trimmed to {len(taxdump):,} nodes for this build')
    else:
        self.taxdump, self.taxid_map, self.stats = (
            taxdump_full, taxid_map, stats)

    # 6. resolve the sequence sources
    fetches = list(provider.protein_sources(df))
    if self.seq_source != 'local':
        fetches = plan_remote_fetches(self, fetches)

    # 7. streaming extraction into db.faa + DIAMOND taxonmap (fixes RK2)
    print('Extracting proteins (streaming, one genome at a time)...')
    tmp_taxonmap = join(self.tmpdir, 'gtdb_prot.accession2taxid')
    tmp_genomemap = join(self.tmpdir, 'gtdb_genome.map')
    out_faa = join(self.output, 'db.faa')
    _gtdb.guard_write(out_faa, root=self.ws_root,
                      must_not_contain=self.gtdb_dir)
    if self.overwrite:
        # a rebuild that only deleted db.faa left the old extraction state and
        # taxon maps behind, so a later "resume" could pair a fresh db.faa with
        # stale maps (or a killed rebuild with stale ones)
        for stale in (extraction_state_path(self.output),
                      join(self.output, 'taxon.map.gz'),
                      join(self.output, 'genome.map.gz'),
                      out_faa + '.partial'):
            if _os.path.isfile(stale):
                _os.remove(stale)
                print(f'  --overwrite: removed stale {basename(stale)}')
    ext = None
    reuse = self.check_local_file(out_faa, self.overwrite)
    if reuse and self.seq_source == 'local':
        # prints its own reason when it declines, so the caller stays quiet
        ext = resume_extraction_state(self, out_faa, tmp_taxonmap, fetches)
    if ext is None:
        def progress(i, total):
            if i % 50 == 0 or i == total:
                print(f'  {i}/{total} genomes', end='\r', flush=True)

        ext = _gtdb.extract_genomes_stream(
            fetches, out_faa, tmp_taxonmap, provider.info['layout'],
            blast_taxidmap=join(self.tmpdir, 'gtdb_seq2taxid.txt'),
            genome_map=tmp_genomemap, on_progress=progress,
            # --hybrid-lca-policy only means something if the build actually
            # asks the provider when a protein has more than one candidate
            resolve=(provider.resolve_shared_protein
                     if self.source == 'hybrid' else None))
        print('')
        write_extraction_state(self.output, ext, fetches,
                               taxonmap=tmp_taxonmap,
                               fingerprint=selection_fingerprint(self))

    # 7b. empty-extraction protection: every broken-input accident used to
    # flow into an empty db.faa and die later inside diamond makedb with a
    # misleading error (B-04); refusing here names the real cause.
    unreadable = int(ext.get('genomes_unreadable', 0))
    if unreadable:
        print(f'  WARNING: {unreadable} genome files were unreadable and are '
              'counted as missing.')
    if int(ext.get('genomes_done', 0)) == 0 or int(ext.get('proteins', 0)) == 0:
        raise ValueError(
            f'Nothing was extracted ({ext.get("genomes_done", 0)} genomes, '
            f'{ext.get("proteins", 0)} proteins; '
            f'{ext.get("genomes_missing", 0)} genomes missing their sequence '
            'file). Check --gtdb-dir, the local layout and --seq-source; '
            'note the remote channels plan downloads but do not perform '
            'them in this version.')
    missing = int(ext.get('genomes_missing', 0)) + unreadable
    if missing:
        # completeness is a gate, not a footnote: a build that silently lost
        # genomes publishes a reference set that turns missing biology into
        # HGT false negatives (C-3.5). --force publishes anyway.
        if not self.force:
            raise ValueError(
                f'{missing} of {ext.get("genomes_planned", 0)} genomes had no '
                'usable sequence file; refusing to publish an incomplete '
                'reference set. Inspect the missing list in the build report '
                'and re-run with --force to publish anyway.')
        print(f'  WARNING: publishing with {missing} missing genomes '
              '(--force).')

    # 8. per-rank lineage table and genome metadata
    write_gtdb_outputs(self, provider, df, ext)

    # 9. compile
    self.gtdb_taxonmap = tmp_taxonmap
    self.gtdb_max_id_len = ext.get('max_id_len', 0)
    compile_gtdb_database(self, args)
    finalize_gtdb_faa(self, args)

    # 10. provenance, fingerprints and report
    write_gtdb_meta(self, provider, args, ext, started, df)
    write_build_report(self, provider, df, ext, started)

    if hasattr(self, 'mkdtemp'):
        _shutil.rmtree(self.tmpdir, ignore_errors=True)
    print(f'Database building finished at {timestamp()}.')
    return self


def plan_remote_fetches(self, fetches):
    """Turn non-local sequence channels into a written fetch plan.

    Parameters
    ----------
    self : Database
        instance holding seq_source / allow_network
    fetches : list
        GenomeFetch items

    Returns
    -------
    list
        the same fetches, with url set

    Notes
    ----
    Network access is disabled unless --allow-network yes is given, so these
    channels export a plan (like the NCBI --manual path) instead of fetching.
    """
    plan = join(self.output, 'gtdb', 'fetch_plan.tsv')
    makedirs(dirname(plan), exist_ok=True)
    if self.seq_source == 'ncbi':
        for i, f in enumerate(fetches):
            acc = _re.sub(r'^(RS|GB)_', '', f.accession)
            kind = 'refseq' if f.accession.startswith('RS_') else 'genbank'
            fetches[i] = _gtdb.GenomeFetch(
                f.accession, f.genome, f.taxid, path=None,
                url=f'https://ftp.ncbi.nlm.nih.gov/genomes/all/{kind}/'
                    f'{acc[:4]}/{acc[:7]}/{acc[:10]}/{acc}_protein.faa.gz',
                kind='remote-ncbi')
    elif self.seq_source == 'gtdb-tar':
        # plan §5.4: this channel is the official *tar bundle*, streamed once
        # and filtered by whitelist -- not one HTTP request per genome, which
        # is what used to be written here while still calling itself gtdb-tar
        if not self.gtdb_release:
            # a missing release used to produce a literal releases/releaseNone/
            # URL in the exported plan
            raise ValueError(
                '--seq-source gtdb-tar needs a release label: pass '
                '--gtdb-release r220 (or similar) so the exported plan does '
                'not point at releases/releaseNone/.')
        host = _gtdb.gtdb_host(getattr(self, 'gtdb_mirror', None))
        tar = (f'https://{host}/releases/release{self.gtdb_release}/'
               'gtdb_proteins_aa_reps.tar.gz')
        for i, f in enumerate(fetches):
            fetches[i] = _gtdb.GenomeFetch(
                f.accession, f.genome, f.taxid, path=None, url=tar,
                kind='remote-gtdb-tar')
    elif self.seq_source == 'ena':
        for i, f in enumerate(fetches):
            acc = _re.sub(r'^(RS|GB)_', '', f.accession)
            fetches[i] = _gtdb.GenomeFetch(
                f.accession, f.genome, f.taxid, path=None,
                url=f'https://www.ebi.ac.uk/ena/browser/api/fasta/{acc}',
                kind='remote-ena')
    with open(plan, 'w') as f:
        f.write(_gtdb.provenance(
            release=getattr(self, 'gtdb_release', None),
            extra={'seq_source': self.seq_source}))
        f.write('# fetch plan written by neoHGT. This version performs NO '
                'network transfer at all (no download executor is '
                'implemented); run the commands below by hand, then re-run '
                'the build with --seq-source local.\n')
        f.write('# the command column honours --downloader; the md5_check column '
                'honours --verify-md5 (for gtdb-tar, verify against MD5SUM.txt '
                'before ingesting)\n')
        if self.seq_source == 'gtdb-tar':
            f.write('# each gtdb-tar command re-reads the release bundle; to '
                    'ingest many genomes, run ONE pass instead: '
                    'curl -L <url> | tar -xz --files-from members.txt, where '
                    'members.txt lists the genomic_files_reps/... members of '
                    'this build (plan §5.4 "stream once, filter by whitelist")\n')
        f.write('accession\turl\tkind\tdownloader\tmd5_check\tcommand\n')
        for x in fetches:
            url = x.url or ''
            name = basename(url) or f'{x.accession}.faa.gz'
            dl = getattr(self, 'downloader', 'curl')
            if x.kind == 'remote-gtdb-tar':
                member = f'genomic_files_reps/{x.accession}_protein.faa.gz'
                fetch = (f'curl -L {url}' if dl == 'curl'
                         else f'wget -q -O - {url}')
                cmd = (f"{fetch} | tar -xzO --wildcards '{member}' "
                       f'> {x.accession}_protein.faa.gz')
            else:
                cmd = (f'curl -L -C - -o {name} {url}' if dl == 'curl'
                       else f'wget -c -O {name} {url}')
            f.write(f'{x.accession}\t{url}\t{x.kind}\t{dl}\t'
                    f'{"yes" if getattr(self, "verify_md5", True) else "no"}'
                    f'\t{cmd}\n')
    print(f'  remote channel "{self.seq_source}": wrote fetch plan to '
          f'{plan}')
    if not self.allow_network:
        raise ValueError(
            f'--seq-source {self.seq_source} needs network access, which this '
            'run has not authorized. Fetch the listed files (or use '
            '--seq-source local, which reads the on-disk GTDB copy), then '
            're-run with --allow-network yes.')
    # No channel has a download executor in this version: with the flag set,
    # the old behaviour planned the fetches, skipped every genome and built an
    # empty database. The URLs this code can generate for the NCBI layout and
    # for the r220-era GTDB bundle were measured 404, and the ENA endpoint
    # returns *nucleotide* sequence, which would poison a protein database.
    # Refusing is the honest behaviour until an executor with alphabet and
    # checksum validation exists.
    raise ValueError(
        f'--seq-source {self.seq_source} has no download executor in this '
        'version: neoHGT would fetch nothing and build an empty database. '
        f'Run the commands in {plan} manually, place the files in a local '
        'directory and rebuild with --seq-source local.')


def write_gtdb_outputs(self, provider, df, ext):
    """Write taxdump/, lineages.txt and genomes.tsv of a GTDB build."""
    out = self.output
    # taxdump
    gtdb.write_taxdump(self.taxdump, join(out, 'taxdump'))
    print('Taxonomy database written to taxdump/.')

    # taxId plus one column per GTDB rank
    df = df.copy()
    df['taxid'] = df['gtdb_taxonomy'].map(self.taxid_map)
    segs = df['gtdb_taxonomy'].str.split(';', expand=True)
    for i, rank in enumerate(gtdb.GTDB_RANKS):
        df[rank] = segs[i]
    df['seq_path'] = [
        (x.path or '') for x in provider.protein_sources(df)]
    per_genome = ext.get('per_genome', {})
    df['proteins_extracted'] = [per_genome.get(a, (0, 0))[0]
                                for a in df['accession']]
    df['residues'] = [per_genome.get(a, (0, 0))[1] for a in df['accession']]
    self.df = df

    # lineages.txt in official GTDB label form (plan §5.6)
    with open(join(out, 'lineages.txt'), 'w') as f:
        f.write(_gtdb.provenance(provider.info))
        f.write('# GTDB lineages (official rank prefixes preserved)\n')
        for lin, tid in zip(df['gtdb_taxonomy'], df['taxid']):
            f.write(f'{tid}\t{lin}\n')
    print('Genome lineages written to lineages.txt.')

    # genome table, indexed by the internal genome id
    cols = ['accession', 'accnov', 'genome', 'domain', 'taxid', 'tax_source',
            'species', 'organism_name', 'gtdb_taxonomy', 'proteins',
            'proteins_extracted', 'residues', 'completeness', 'contamination',
            'quality_score', 'category', 'assembly_level', 'ncbi_taxid',
            'seq_path']
    tab = df.copy()
    for col in cols:
        if col not in tab:
            tab[col] = ''
    tab = tab[cols].set_index('genome')
    tab.to_csv(join(out, 'genomes.tsv'), sep='\t')
    print(f'Genome metadata written to genomes.tsv ({len(tab):,} rows).')

    # gtdb/ side files, including route C tables
    gtdb.write_side_files(
        out, provider.info, provider.df_all, self.taxid_map, self.taxdump,
        self.stats, stem2taxid=ext.get('stem2taxid'),
        prefix2genome=ext.get('prefix2genome'),
        previous_manifest=getattr(self, 'gtdb_manifest', None))
    print('GTDB derived tables written to gtdb/.')
    return self


def check_capacity(self, df):
    """Disk pre-flight for a GTDB build (plan §11.2, risk RK18).

    Parameters
    ----------
    self : Database
        instance holding ``output`` and ``force``
    df : pandas.DataFrame
        the selected genomes (needs the ``proteins`` column)

    Returns
    -------
    dict
        planned proteins, estimated bytes and the free space that was seen

    Raises
    ------
    ValueError
        When the output volume cannot hold the estimated peak and ``--force``
        was not given.
    """
    plan_proteins = int(df['proteins'].sum())
    est_bytes = int(plan_proteins * GTDB_BYTES_PER_PROTEIN)
    free = _shutil.disk_usage(self.output).free
    print(f'Estimated db.faa size: {est_bytes / 2 ** 30:.1f} GiB '
          f'({plan_proteins:,} proteins); peak ~'
          f'{est_bytes * GTDB_PEAK_FACTOR / 2 ** 30:.1f} GiB; '
          f'free on output volume: {free / 2 ** 30:.1f} GiB')
    if free < est_bytes * GTDB_PEAK_FACTOR:
        if not self.force:
            raise ValueError(
                f'Insufficient free space on {self.output}: '
                f'{free / 2 ** 30:.1f} GiB available, ~'
                f'{est_bytes * GTDB_PEAK_FACTOR / 2 ** 30:.1f} GiB needed at '
                'peak. Reduce the selection (--rank-sample/--gtdb-taxa), use '
                '--keep-faa compress or --keep-faa drop, or pass --force to '
                'ignore this check.')
        print('  WARNING: continuing despite the disk pre-flight (--force).')
    return {'plan_proteins': plan_proteins, 'estimated_bytes': est_bytes,
            'free_bytes': free, 'peak_factor': GTDB_PEAK_FACTOR}


def expand_selection(self, provider, df):
    """Grow a representative selection to the member genomes it stands for.

    Parameters
    ----------
    self : Database
        instance holding output / workspace boundaries
    provider : GtdbProvider
        prepared provider with df_all and layout information
    df : pandas.DataFrame
        the current genome selection

    Returns
    -------
    pandas.DataFrame
        selection extended with the member genomes that have sequence

    Notes
    ----
    The member -> representative relation comes from the metadata column
    'gtdb_genome_representative' (the official sp_clusters.tsv is absent from
    local copies, plan §3.4). Every member is listed in gtdb/members.tsv.gz,
    whether or not its sequence is available; if none are, the build stops with
    an actionable message instead of quietly shipping a smaller database.
    """
    print('Indexing local sequence files for member expansion...')
    index = gtdb.stem_index(provider.info)
    members, rep = gtdb.expand_to_members(provider.df_all, df, seq_index=index)
    path = gtdb.guard_write(join(self.output, 'gtdb', 'members.tsv.gz'),
                            root=self.ws_root, must_not_contain=self.gtdb_dir)
    makedirs(dirname(path), exist_ok=True)
    cols = [c for c in ('accession', 'genome', 'organism_name',
                       'gtdb_taxonomy', 'checkm2_completeness',
                       'checkm2_contamination', 'representative',
                       'has_local_fasta') if c in members.columns]
    with gzip.open(path, 'wt', encoding='utf-8') as f:
        f.write(gtdb.provenance(provider.info,
                                extra={'selection_genomes': int(df.shape[0])}))
        f.write('\t'.join(cols) + '\n')
        for r in members[cols].itertuples(index=False):
            f.write('\t'.join('' if v is None or (isinstance(v, float) and
                                                  v != v) else str(v)
                              for v in r) + '\n')
    have = (int(members['has_local_fasta'].sum())
            if 'has_local_fasta' in members else int(members.shape[0]))
    print(f'Member genomes: {rep["members_in_metadata"]:,} in metadata, '
          f'{have:,} with local sequence; table written to '
          f'gtdb/members.tsv.gz.')
    if members.empty:
        raise ValueError(
            'The selected genomes carry no member genomes in this metadata '
            'table, so --expand-members cannot grow the selection.')
    usable = (members[members['has_local_fasta']]
              if 'has_local_fasta' in members else members)
    if usable.empty:
        raise ValueError(
            f'None of the {rep["members_in_metadata"]:,} member genomes has a '
            'local sequence file. This dataset only stores species-cluster '
            'representatives; the full member list is in '
            f'{path}. Obtain the non-representative protein files (official '
            'GTDB *_genomic_files_all tarballs via --seq-source gtdb-tar, or '
            '--seq-source ena) and point --gtdb-dir at them, or drop '
            '--expand-members.')
    return pd.concat([df, usable], ignore_index=True)


def extraction_state_path(output):
    """Path of the per-build extraction state used to resume reliably."""
    return join(output, 'gtdb', 'extraction_state.json')


def selection_sha(fetches, fingerprint=None):
    """Fingerprint the genome selection of a build (resume safety check).

    The accession list alone used to be the whole fingerprint: two runs with
    the same genomes but a different ``--gtdb-taxid-base``, a different local
    layout, or a different dataset root share it, so a resume would keep the
    old db.faa (and its old taxIds) while meta.json claimed the new settings.
    The extra context makes those resumes refuse.
    """
    parts = sorted(x.accession for x in fetches)
    if fingerprint:
        parts.append(fingerprint)
    return gtdb.sha256_text('\n'.join(parts))


def selection_fingerprint(self):
    """Everything besides the accession list that changes the extracted bytes.

    Includes the taxId base (taxIds live in taxon.map.gz, not in db.faa, but a
    resume that keeps db.faa also keeps the old taxId map), the local layout
    (it decides the seqid parsing) and the realpath of the input dataset (two
    mirrors of one release are not interchangeable once provenance matters).
    """
    return '\n'.join((
        f'id_base={getattr(self, "id_base", "")}',
        f'layout={getattr(self, "local_layout", "") or ""}',
        f'gtdb_dir={(_os.path.realpath(self.gtdb_dir) if getattr(self, "gtdb_dir", None) else "")}',
    ))


EXTRACTION_STATE_SCHEMA = 2


def write_extraction_state(output, ext, fetches, taxonmap=None,
                           fingerprint=None):
    """Persist extraction counters so an interrupted build can be resumed.

    Parameters
    ----------
    output : str
        database directory
    ext : dict
        counters from gtdb.extract_genomes_stream()
    fetches : list
        genome selection this state belongs to
    taxonmap : str, optional
        path of the streamed DIAMOND taxonmap, recorded for provenance
    fingerprint : str, optional
        selection_fingerprint() of this run (id base, layout, dataset root)

    Returns
    -------
    str
        the state file path
    """
    path = gtdb.guard_write(extraction_state_path(output))
    makedirs(dirname(path), exist_ok=True)
    payload = {k: v for k, v in ext.items() if k != 'unresolved'}
    # the taxonmap sits in the scratch directory, which the workspace rules
    # require to be disposable; record the name, not a path that may vanish
    payload['taxonmap_scratch_name'] = basename(taxonmap) if taxonmap else ''
    payload['selection_sha256'] = selection_sha(fetches, fingerprint)
    payload['schema_version'] = EXTRACTION_STATE_SCHEMA
    payload['unresolved'] = list(ext.get('unresolved', []))[:2000]
    with open(path, 'w') as f:
        _json.dump(payload, f, default=str)
    return path


def resume_extraction_state(self, out_faa, tmp_taxonmap, fetches):
    """Rebuild the extraction artefacts of a completed pass from disk.

    Returns
    -------
    dict or None
        the stored counters, or None when a consistent resume is impossible

    Notes
    ----
    Reusing db.faa is only safe together with the database's taxon.map.gz and
    genome.map.gz, from which the DIAMOND taxonmap and the per-genome counters
    are regenerated. Without the state file the counters cannot be recovered
    without re-reading hundreds of GB, so the build refuses rather than
    publishing a meta.json that disagrees with db.faa. A state written for a
    different genome selection is also refused: db.faa would no longer match
    the genomes.tsv / meta.json of this run.
    """
    state_p = extraction_state_path(self.output)
    if not (isfile(out_faa) and isfile(state_p)):
        print('  no extraction_state.json next to db.faa; extracting again.')
        return None
    with open(state_p) as f:
        ext = _json.load(f)
    fp = selection_fingerprint(self)
    if ext.get('selection_sha256') != selection_sha(fetches, fp):
        print('  selection, taxId base, layout or dataset changed since '
              'db.faa was built; re-extracting instead of resuming.')
        return None
    version = ext.get('schema_version')
    if version is not None and int(version) != EXTRACTION_STATE_SCHEMA:
        print(f'  extraction_state.json schema {version} != '
              f'{EXTRACTION_STATE_SCHEMA}; re-extracting.')
        return None
    # "exists" is not "complete": compare db.faa against the byte count the
    # state itself recorded (C-3.7: a killed re-run used to leave a truncated
    # db.faa that a third run resumed into a silently short database)
    want_bytes = ext.get('bytes')
    got_bytes = _os.path.getsize(out_faa)
    if want_bytes and int(want_bytes) != got_bytes:
        print(f'  db.faa is {got_bytes:,} B but the state recorded '
              f'{int(want_bytes):,} B; the file is truncated or stale, '
              're-extracting.')
        return None
    taxon_src = join(self.output, 'taxon.map.gz')
    gmap_src = join(self.output, 'genome.map.gz')
    if not isfile(taxon_src):
        raise ValueError(
            f'{out_faa} exists but taxon.map.gz does not: cannot resume '
            'safely. Re-run with --overwrite to rebuild from scratch.')
    taxon_rows = 0
    with open(tmp_taxonmap, 'w') as fo:
        fo.write(gtdb.DIAMOND_TAXONMAP_HEADER)
        with gzip.open(taxon_src, 'rt') as fi:
            for line in fi:
                x = line.rstrip('\r\n').split('\t')
                if len(x) >= 2:
                    fo.write(f'{x[0]}\t{x[0]}\t{x[1]}\tna\n')
                    taxon_rows += 1
    if ext.get('proteins') and taxon_rows != int(ext['proteins']):
        raise ValueError(
            f'taxon.map.gz holds {taxon_rows:,} proteins but the extraction '
            f'state recorded {int(ext["proteins"]):,}: the artifacts '
            f'disagree. Re-run with --overwrite to rebuild from scratch.')
    if isfile(gmap_src):
        with open(join(self.tmpdir, 'gtdb_genome.map'), 'w') as fo, \
                gzip.open(gmap_src, 'rt') as fi:
            for line in fi:
                fo.write(line)
    with open(join(self.tmpdir, 'gtdb_seq2taxid.txt'), 'w') as fo:
        fo.write(gtdb.BLAST_TAXIDMAP_HEADER)
        with gzip.open(taxon_src, 'rt') as fi:
            for line in fi:
                x = line.rstrip('\r\n').split('\t')
                if len(x) >= 2:
                    fo.write(f'{x[0]}\t{x[1]}\n')
    ext['per_genome'] = {k: tuple(v) for k, v in
                         (ext.get('per_genome') or {}).items()}
    ext['resumed'] = True
    print(f'  Resumed from existing db.faa + {basename(state_p)} '
          f'({ext.get("genomes_done", 0):,} genomes, '
          f'{ext.get("proteins", 0):,} proteins).')
    return ext


def compile_gtdb_database(self, args):
    """Compile the DIAMOND (optionally BLAST) database of a GTDB build."""
    compile_ = getattr(args, 'compile', None) or 'diamond'
    if compile_ == 'none':
        return
    if compile_ in ('blast', 'both'):
        if not getattr(self, 'force', False):
            raise ValueError(
                'The BLAST route is unavailable for GTDB databases (plan '
                'ADR-8): makeblastdb caps sequence ids at 50 characters, '
                f'while this build produced ids up to '
                f'{getattr(self, "gtdb_max_id_len", 270)} characters. Use '
                '--compile diamond, or pass --force to attempt it anyway.')
        print('  WARNING: attempting the experimental BLAST route (--force).')
        makedirs(join(self.output, 'blast'), exist_ok=True)
        # argv list, not a shell string: workspace paths with spaces or shell
        # metacharacters used to be split or executed
        cmd = [
            self.makeblastdb, '-dbtype', 'prot',
            '-in', join(self.output, 'db.faa'),
            '-out', join(self.output, 'blast', 'db'),
            '-title', 'db', '-parse_seqids',
            '-taxid_map', join(self.tmpdir, 'gtdb_seq2taxid.txt')]
        ec, _ = run_command(cmd)
        if ec:
            raise ValueError(f'makeblastdb failed with error code {ec}.')
        print('BLAST database written to blast/.')
    if compile_ in ('diamond', 'both'):
        makedirs(join(self.output, 'diamond'), exist_ok=True)
        # DIAMOND indexes nodes.dmp/names.dmp by taxId and is slow with the
        # full GTDB tree, so makedb gets only the taxa this map references;
        # the complete tree remains in <db>/taxdump for search and analyze.
        used = gtdb.taxids_in_taxonmap(self.gtdb_taxonmap)
        compact_dir = join(self.tmpdir, 'taxdump_compact')
        gtdb.write_compact_taxdump(self.taxdump, used, compact_dir)
        print(f'  DIAMOND sees {len(used):,} taxa (compact taxdump), while '
              f'<db>/taxdump keeps all {len(self.taxdump):,} nodes.')
        # The taxonmap was streamed during extraction, so 333M rows never live
        # in memory. Both key columns hold the identical full sequence Id:
        # stripping ".version" (as the NCBI branch does) would break matching
        # for ids like '..._SerialNumber_AE017199.1_1'.
        # argv list, not a shell string (same reason as makeblastdb above)
        cmd = [
            self.diamond, 'makedb',
            '--threads', str(self.threads or 1),
            '--no-parse-seqids',
            '--in', join(self.output, 'db.faa'),
            '--taxonmap', self.gtdb_taxonmap,
            '--taxonnodes', join(compact_dir, 'nodes.dmp'),
            '--taxonnames', join(compact_dir, 'names.dmp'),
            '--db', join(self.output, 'diamond', 'db'),
            '--tmpdir', self.tmpdir]
        print('Build DIAMOND database...', flush=True)
        ec, out = run_command(cmd)
        if ec:
            raise ValueError(
                f'diamond makedb failed with error code {ec}:\n'
                + '\n'.join(out or []))
        print('DIAMOND database written to diamond/db.dmnd.')

    # protein-to-taxon map kept as a compressed artifact (search --taxmap)
    src = self.gtdb_taxonmap
    dst = join(self.output, 'taxon.map.gz')
    if isfile(src):
        with open(src) as fi, gzip.open(dst + '.partial', 'wt') as fo:
            next(fi, None)
            for line in fi:
                x = line.rstrip('\r\n').split('\t')
                if len(x) >= 3:
                    fo.write(f'{x[1]}\t{x[2]}\n')
        _os.replace(dst + '.partial', dst)
        print(f'Protein-to-taxonomy map written to {basename(dst)}.')

    gm = join(self.tmpdir, 'gtdb_genome.map')
    if isfile(gm):
        with open(gm) as fi, gzip.open(join(self.output, 'genome.map.gz')
                                       + '.partial', 'wt') as fo:
            for line in fi:
                fo.write(line)
        _os.replace(join(self.output, 'genome.map.gz') + '.partial',
                    join(self.output, 'genome.map.gz'))
        print('Protein-to-genome map written to genome.map.gz.')

    # db.faa disposition is handled by finalize_gtdb_faa: it used to sit behind
    # the `--compile none` early return, so --keep-faa did nothing unless a
    # DIAMOND database was built
    return self


def finalize_gtdb_faa(self, args):
    """Honour `--keep-faa` for a GTDB build, whatever was compiled.

    Parameters
    ----------
    self : Database
        instance holding keep_faa / output / tmpdir
    args : argparse.Namespace
        parsed arguments (the compile choice is reported back)

    Returns
    -------
    Database
        the same instance
    """
    faa = join(self.output, 'db.faa')
    if self.keep_faa == 'compress' and isfile(faa):
        with open(faa) as fi, gzip.open(faa + '.gz', 'wt') as fo:
            for line in fi:
                fo.write(line)
        remove(faa)
        print('  db.faa compressed (--keep-faa compress).')
    elif self.keep_faa == 'drop' and isfile(faa):
        # the choice is the user's explicit request: keeping the file while
        # claiming to honour "drop" was a no-op tier
        remove(faa)
        print('  db.faa removed (--keep-faa drop); the compiled database '
              'under diamond/ (and blast/) is the surviving artifact.')
    return self


def write_gtdb_meta(self, provider, args, ext, started, df=None):
    """Write meta.json plus input fingerprints.

    `df` is the selected genome table; the manifest's protein reconciliation
    has to be computed against it rather than against the provider's loaded
    table, which holds every genome of the release.
    """
    stats = dict(ext)
    if self.fingerprints:
        files = list(provider.info['taxonomy']) + list(
            provider.info['metadata'])
        fp = join(self.output, 'gtdb', 'input_fingerprints.tsv')
        makedirs(dirname(fp), exist_ok=True)
        stats['input_sha256'] = gtdb.fingerprint(files, fp,
                                                 info=provider.info)
    tm_sha = None
    if isfile(self.gtdb_taxonmap):
        tm_sha = gtdb.sha256_file(self.gtdb_taxonmap)
    dv = None
    # --compile none never locates diamond, and shutil.which(None) raises
    if self.diamond and which(self.diamond):
        try:
            ec, out = run_command(f'{self.diamond} --version')
            dv = '\n'.join(out or []).strip()
        except Exception:
            dv = None
    opts = {k: v for k, v in vars(args).items()}
    gtdb.write_meta_json(
        join(self.output, 'meta.json'), provider, opts, stats,
        taxonmap_sha=tm_sha, diamond_version=dv, started=str(started),
        finished=str(timestamp()), config_snapshot=self.cfg,
        selection=df)
    print('Build manifest written to meta.json.')
    return self


def write_build_report(self, provider, df, ext, started):
    """Write gtdb/build_report.md (plan §5.6)."""
    path = join(self.output, 'gtdb', 'build_report.md')
    makedirs(dirname(path), exist_ok=True)
    rep = provider.report or {}
    lines = [
        '# GTDB build report',
        '',
        f'- dataset: `{provider.info["gtdb_dir"]}` (opened read-only)',
        f'- release: `{provider.info["release"]}`, layout: '
        f'`{provider.info["layout"]}`, id convention: '
        f'`{provider.info["id_convention"]}`',
        f'- sequence-Id anchor coverage (`{gtdb.SERIAL_ANCHOR}`): '
        f'{provider.info["serial_anchor_coverage"]:.2%}',
        f'- longest defline observed: {provider.info["max_first_line_len"]} '
        'chars',
        f'- genomes selected: {df.shape[0]:,}; proteins emitted: '
        f'{ext.get("proteins", 0):,}; duplicates skipped: '
        f'{ext.get("proteins_dup", 0):,}; without taxId: '
        f'{ext.get("proteins_no_taxid", 0):,}',
        f'- missing sequence files: {ext.get("genomes_missing", 0)}',
        f'- db.faa: {ext.get("bytes", 0) / 2 ** 30:.1f} GiB',
        f'- quality classes in the dataset: ' + (
            ', '.join(f'{k}={v:,}' for k, v in
                      sorted(getattr(provider, 'quality_classes', {}).items()))
            or 'not reported by this source'),
        '',
        '## Input self-consistency checks',
        '',
        '| check | expected | observed | status |',
        '| --- | --- | --- | --- |',
    ]
    for c in rep.get('checks', []):
        lines.append(f'| {c["name"]} | {c["expected"]} | {c["observed"]} | '
                     f'{c["status"]} |')
    lines += ['', '## Provenance of the local copy', '']
    for p in provider.info['provenance_scripts']:
        lines.append(f'- `{p}`')
    lines += ['', '## Official files absent from this local copy', '']
    for p in provider.info['absent_official']:
        lines.append(f'- `{p}`')
    if ext.get('unresolved'):
        lines += ['', '## Unresolved genomes', '']
        lines += [f'- {a}' for a in ext['unresolved'][:200]]
    if provider.name == 'hybrid':
        # plan §5.7: a hybrid database must say which NCBI subtree was grafted
        # in and how, because that graft is what makes eukaryote/virus donors
        # reachable at all
        st = provider.stats or {}
        lines += ['', '## Hybrid graft (plan §5.7)', '',
                  f'- NCBI taxdump: `{getattr(provider, "ncbi_dir", "?")}`',
                  f'- kept ranks: {", ".join(provider.keep_ranks)}',
                  f'- LCA policy for shared proteins: `{st.get("lca_policy")}`',
                  f'- proteins with more than one taxId candidate: '
                  f'{ext.get("proteins_shared", 0):,} '
                  f'(reassigned {ext.get("proteins_shared_reassigned", 0):,}, '
                  f'left unassigned {ext.get("proteins_shared_dropped", 0):,})'
                  ' — the policy is applied to exactly these, so a count of 0 '
                  'means it had nothing to resolve in this build',
                  f'- NCBI nodes merged: {st.get("n_ncbi_nodes", 0):,}',
                  f'- nodes in merged tree: {st.get("n_nodes_merged", 0):,}',
                  '- ranks flattened for DIAMOND in `nodes.dmp` (the true rank '
                  'stays in the tree): '
                  f'{st.get("ranks_flattened_for_diamond") or "none"}']
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f'Build report written to {path}.')
    return self


_ncbi_call = Database.__call__


def call_with_source(self, args):
    """Run the NCBI pipeline unchanged, or the GTDB/hybrid one.

    Parameters
    ----------
    self : Database
        instance
    args : argparse.Namespace
        parsed arguments

    Notes
    ----
    The NCBI branch delegates to the original __call__ without any change, so
    its behavior stays bit-identical (plan G4, DoD-7, RT-1). The GTDB branch
    reuses set_parameters for argument handling, then follows its own
    provider-driven pipeline.
    """
    # one source resolution for the whole stage: set_source_parameters()
    # reads args.source, then config, then 'ncbi' -- the router used to read
    # args.source alone, which made `database.source: gtdb` in config.yml a
    # dead key that the shipped comment promised would work
    source = (getattr(args, 'source', None)
              or ((self.cfg or {}).get('database') or {}).get('source')
              or 'ncbi')
    if source == 'ncbi':
        if arg2bool(getattr(args, 'stream_extraction', None)):
            bind_stream_extraction(self)
        return _ncbi_call(self, args)
    print(f'Database building started at {timestamp()}.')
    se = getattr(args, 'stream_extraction', None)
    if se is not None and arg2bool(se) is False:
        # the GTDB route has no non-streaming extractor any more (RK2), so say
        # out loud that the request cannot be honoured rather than letting the
        # user believe the in-memory path ran
        print('NOTE: --stream-extraction no is an NCBI-route switch; a GTDB '
              'build always streams, because holding a 333M-protein database '
              'in memory is not possible (plan §11.2 item 2). Continuing with '
              'streaming extraction.')
    self.cfg = load_configs()
    self.set_parameters(args)
    set_source_parameters(self, args)
    return build_from_source(self, args)


Database.__call__ = call_with_source


# ===========================================================================
# P7-5: memory-bounded extraction for the NCBI route (plan §6.7). Opt-in via
# --stream-extraction yes; with the default 'no' the NCBI pipeline is the
# untouched upstream code.
# ===========================================================================

arguments += [
    'NCBI extraction memory fix (P7)',
    ['--stream-extraction', 'extract proteins of the NCBI route without '
     'keeping the protein-to-genome / protein-to-taxid maps in RAM: dedup '
     'through a disk-backed index and external sorting in --tmpdir. Produces '
     'byte-identical db.faa / genome.map.gz / taxon.map.gz; needs extra '
     'temporary disk (~30 bytes per protein occurrence).',
     {'choices': ['yes', 'no']}],
]


def stream_paths(self):
    """Temporary files used by the streaming NCBI extraction."""
    return (join(self.tmpdir, 'stream_prot2genome.tsv'),
            join(self.tmpdir, 'stream_prot2taxid.tsv'))


def _external_sort(src, dst):
    """Sort a TSV by its first then second field, on disk.

    Notes
    ----
    LC_ALL=C forces byte ordering, which is what Python's sorted() does on
    str keys; without it a UTF-8 locale could order the map differently from
    the upstream implementation.
    """
    env = dict(_os.environ, LC_ALL='C')
    cmd = ['sort', '-t', '\t', '-k1,1', '-k2,2', '-o', dst, src]
    p = _subprocess.run(cmd, env=env)
    if p.returncode:
        raise ValueError(f'sort failed with exit code {p.returncode}.')
    return dst


def _grouped(src):
    """Yield (key, [values]) from a key-sorted two-column TSV, streaming."""
    key, vals = None, []
    with open(src) as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            if len(x) < 2:
                continue
            if x[0] != key:
                if key is not None:
                    yield key, vals
                key, vals = x[0], [x[1]]
            else:
                vals.append(x[1])
    if key is not None:
        yield key, vals


def extract_genomes_stream_ncbi(self):
    """Streaming equivalent of Database.extract_genomes.

    Notes
    -----
    The upstream method keeps one dict entry per unique protein; a full NCBI
    build therefore runs out of memory (the plan's RK2, and the '☕️' warning
    printed by the original). Here proteins are written out as soon as they
    are complete, duplicate detection uses a sqlite index in --tmpdir, and the
    protein-to-genome / protein-to-taxid maps are built by external sort plus
    a streaming group-by. The emitted files match the original byte for byte.
    """
    import sqlite3
    print('Extracting downloaded genomic data... (streaming)')
    ldir = join(self.down, 'faa')
    pg, pt = stream_paths(self)
    seen_db = join(self.tmpdir, 'stream_seen.sqlite3')
    for p in (pg, pt, seen_db):
        if isfile(p):
            remove(p)
    con = sqlite3.connect(seen_db)
    con.execute('PRAGMA journal_mode=OFF')
    con.execute('PRAGMA synchronous=OFF')
    con.execute('CREATE TABLE seen(acc TEXT PRIMARY KEY, aa INTEGER)')
    cur = con.cursor()
    fout = open(join(self.output, 'db.faa'), 'w')
    fg, ft = open(pg, 'w'), open(pt, 'w')

    g2n, g2aa = {}, {}
    counter = 0
    for row in self.df.itertuples():
        g, tid = row.genome, row.taxid
        g2n[g], g2aa[g] = 0, 0
        stem = row.ftp_path.rsplit('/', 1)[-1]
        lfile = join(ldir, f'{stem}_protein.faa.gz')
        if not exists(lfile):
            print(f'❌ {lfile} does not exist.')
            continue
        try:
            with gzip.open(lfile, 'rb') as f:
                content = f.read().decode().splitlines()
        except (OSError, EOFError, TypeError):
            print(f' skipping corrupted file {stem}.', end='', flush=True)
            continue

        cp, cname, buf = None, '', []

        def flush():
            """Emit the pending protein and remember its length."""
            nonlocal cp, buf
            if cp is None:
                return
            seq = ''.join(buf)
            cur.execute('INSERT OR IGNORE INTO seen VALUES(?,?)',
                        (cp, len(seq)))
            fout.write(f'>{cp} {cname}\n{seq}\n')
            cp, buf = None, []

        for line in content:
            if line.startswith('>'):
                flush()
                p, name = line[1:].split(None, 1)
                g2n[g] += 1
                hit = cur.execute(
                    'SELECT aa FROM seen WHERE acc=?', (p,)).fetchone()
                if hit is not None:
                    # shared protein: not written again, counted only
                    g2aa[g] += hit[0]
                    cp = None
                else:
                    cp, cname, buf = p, name, []
                fg.write(f'{p}\t{g}\n')
                ft.write(f'{p}\t{tid}\n')
            elif cp:
                line = line.rstrip('*')
                buf.append(line)
                g2aa[g] += len(line)
        flush()
        counter += 1
        print(f'{counter}/{self.df.shape[0]}: Extracted {stem}.')

    fout.close()
    fg.close()
    ft.close()
    con.close()
    remove(seen_db)
    print('✅ Done.')
    print('Combined protein sequences written to db.faa.')

    # protein-to-genome(s) map, from the external sort
    srt_g = _external_sort(pg, pg + '.sorted')
    fname = 'genome.map.gz'
    n_multi_genome = 0
    with gzip.open(join(self.output, fname), 'wb') as out:
        for p, gs in _grouped(srt_g):
            if len(gs) > 1:
                n_multi_genome += 1
            out.write('{}\t{}\n'.format(p, ','.join(sorted(set(gs)))).encode())
    print(f'Protein-to-genome(s) map written to {fname}.')

    # protein-to-taxonomy map: LCA of the taxIds, computed per protein only
    srt_t = _external_sort(pt, pt + '.sorted')
    n_multi_taxid = 0
    blast_map = join(self.tmpdir, 'taxon.map')
    diamond_map = join(self.tmpdir, 'prot.accession2taxid')
    fb, fd = open(blast_map, 'w'), open(diamond_map, 'w')
    fd.write('accession\taccession.version\ttaxid\tgi\n')
    with gzip.open(join(self.output, 'taxon.map.gz'), 'wb') as out:
        for p, tids in _grouped(srt_t):
            uniq = set(tids)
            if len(uniq) > 1:
                n_multi_taxid += 1
                tid = find_lca(sorted(uniq), self.taxdump)
            else:
                tid = max(uniq)
            out.write(f'{p}\t{tid}\n'.encode())
            fb.write(f'{p}\t{tid}\n')
            fd.write(f'{p.rsplit(".", 1)[0]}\t{p}\t{tid}\tna\n')
    fb.close()
    fd.close()
    print('Protein-to-taxonomy map written to taxon.map.gz.')
    remove(pg)
    remove(pt)

    self.df['proteins'] = self.df['genome'].map(g2n)
    self.df['residues'] = self.df['genome'].map(g2aa)
    self.p2tids = None
    self.taxonmap = None
    self.stream_maps = True
    n_unique = sum(1 for _ in _grouped(srt_g))
    print(f'  Total number of genomes extracted: {len(g2n)}.')
    print(f'  Total number of unique proteins extracted: {n_unique:,}.')
    print(f'  Number of proteins shared by multiple genomes: '
          f'{n_multi_genome:,}.')
    print(f'  Number of proteins shared by multiple TaxIDs: {n_multi_taxid:,}.')
    return self


def build_taxonmap_stream(self):
    """No-op: the streaming extractor already wrote taxon.map.gz."""
    print('Protein-to-taxonomy map already written (streaming extraction).')
    return self


def build_blast_db_stream(self):
    """BLAST build using the taxid map written during extraction."""
    makedirs(join(self.output, 'blast'), exist_ok=True)
    cmd = [
        self.makeblastdb,
        '-dbtype', 'prot',
        '-in', join(self.output, 'db.faa'),
        '-out', join(self.output, 'blast', 'db'),
        '-title', 'db',
        '-parse_seqids',
        '-taxid_map', join(self.tmpdir, 'taxon.map')]
    print('Build BLAST database...', flush=True)
    ec, out = run_command(cmd)
    if ec:
        raise ValueError(f'makeblastdb failed with error code {ec}.')
    remove(join(self.tmpdir, 'taxon.map'))
    print('Done.')
    return self


def build_diamond_db_stream(self):
    """DIAMOND build using the accession map written during extraction."""
    makedirs(join(self.output, 'diamond'), exist_ok=True)
    taxonmap = join(self.tmpdir, 'prot.accession2taxid')
    cmd = [
        self.diamond, 'makedb',
        '--threads', str(self.threads),
        '--in', join(self.output, 'db.faa'),
        '--taxonmap', taxonmap,
        '--taxonnodes', join(self.output, 'taxdump', 'nodes.dmp'),
        '--taxonnames', join(self.output, 'taxdump', 'names.dmp'),
        '--db', join(self.output, 'diamond', 'db'),
        '--tmpdir', self.tmpdir]
    print('Build DIAMOND database...', flush=True)
    ec, out = run_command(cmd)
    if ec:
        raise ValueError(f'diamond failed with error code {ec}.')
    remove(taxonmap)
    print('Done.')
    return self


def bind_stream_extraction(self):
    """Swap the NCBI methods for their memory-bounded equivalents."""
    from types import MethodType
    self.extract_genomes = MethodType(extract_genomes_stream_ncbi, self)
    self.build_taxonmap = MethodType(build_taxonmap_stream, self)
    self.build_blast_db = MethodType(build_blast_db_stream, self)
    self.build_diamond_db = MethodType(build_diamond_db_stream, self)
    print('  Streaming extraction enabled (--stream-extraction yes).')
    return self

