"""
////////////////////////////////////////////////////////////
│	
│	Filename:	search.py
│	Comment:	See documentation
│	
│	=======================================================
│	Authorship:	Zhihui Li (Qiyun Zhu and Katharina Dittmar)
│	Copyright:	Distributed under the terms of the Modified BSD License.
│	
////////////////////////////////////////////////////////////
"""
#!/usr/bin/env python3

import re
from os import remove, makedirs, cpu_count
from os.path import join, isdir, isfile
from shutil import which, rmtree
from tempfile import mkdtemp
from time import time, sleep
from math import log
from urllib.parse import quote, urlencode
from urllib.request import urlopen, Request, HTTPError, URLError

from neoHGT.util import (
	timestamp, load_configs, get_config, arg2bool, list_from_param, file2id,
	id2file_map, read_taxdump, read_input_prots, read_prot2taxid,
	get_product, seqid2accver, write_fasta, run_command, contain_words,
	is_latin, is_capital, is_ancestral, taxid_at_rank)


description = """batch sequence homology searching and filtering"""

arguments = [
	'basic',
	['-i|--input',   'input protein set file, or directory where one or more '
					 'input files are located', {'required': True}],
	['-o|--output',  'directory where search results are to be saved',
					 {'required': True}],
	['-m|--method',  'search method',
					 {'choices': ['auto', 'diamond', 'blast', 'remote',
					  'precomp'], 'default': 'auto'}],
	['-s|--precomp', 'file or directory of precomputed search results (when '
					 'method = precomp)'],

	'database',
	['-d|--db',      'reference protein sequence database'],
	['-t|--taxdump', 'directory of taxonomy database files (nodes.dmp and '
					 'names.dmp)'],
	['--taxmap',     'sequence Id to taxId mapping file (not necessary if '
					 'protein database already contains taxonomy)'],

	'search behaviors',
	['-k|--maxhits', 'maximum number of hits to preserve per query (0 for '
					 'unlimited)', {'type': int}],
	['--minsize',    'minimum length of query sequence (aa)', {'type': int}],
	['--queries',    'number of queries per run (0 for whole sample)',
					 {'type': int}],
	['--maxchars',   'maximum number of characters per run (0 for unlimited)',
					 {'type': int}],

	'search cutoffs',
	['--maxseqs',    'maximum number of sequences to return', {'type': int}],
	['--evalue',     'maximum E-value cutoff', {'type': float}],
	['--identity',   'minimum percent identity cutoff', {'type': int}],
	['--coverage',   'minimum percent query coverage cutoff', {'type': int}],
	['--extrargs',   'extra arguments for choice of search method'],

	'taxonomic filters',
	['--tax-include', 'include taxa under those taxIds'],
	['--tax-exclude', 'exclude taxa under those taxIds'],
	['--tax-unique',  'ignore more than one hit with same taxID',
					  {'choices': ['yes', 'no']}],
	['--tax-unirank', 'ignore more than one hit under same taxon at this '
					  'rank'],
	['--tax-capital', 'ignore taxon names that are not capitalized',
					  {'choices': ['yes', 'no']}],
	['--tax-latin',   'ignore species names that are not Latinate',
					  {'choices': ['yes', 'no']}],
	['--tax-block',   'ignore taxon names containing those words'],

	'local search behaviors',
	['-p|--threads', 'number of threads (0 for all CPU cores)',
					 {'type': int}],
	['--tmpdir',     'temporary directory'],
	['--diamond',    'diamond executable'],
	['--blastp',     'blastp executable'],
	['--blastdbcmd', 'blastdbcmd executable'],

	'remote search behaviors',
	['--algorithm',  'remote search algorithm'],
	['--retries',    'maximum number of retries per search',
					 {'type': int}],
	['--delay',      'seconds between two search requests',
					 {'type': int}],
	['--timeout',    'seconds before program gives up waiting',
					 {'type': int}],
	['--entrez',     'entrez query text'],
	['--server',     'remote search server URL'],

	'self-alignment options',
	['--aln-method',  'self-alignment method',
					  {'choices': ['auto', 'native', 'fast', 'lookup',
					   'precomp']}],
	['--aln-precomp', 'file or directory of precomputed sequence Id to score '
					  'maps (when self-alignment method = precomp)'],
	['--aln-server',  'remote align server URL'],

	'remote fetch options',
	['--fetch-enable',  'whether to enable remote fetch',
						{'choices': ['auto', 'yes', 'no']}],
	['--fetch-queries', 'maximum number of query entries per search'],
	['--fetch-retries', 'maximum number of retries per search'],
	['--fetch-delay',   'seconds between two fetch requests',
						{'type': int}],
	['--fetch-timeout', 'seconds before program gives up waiting',
						{'type': int}],
	['--fetch-server',  'remote fetch server URL'],
]


class Search(object):

	def __init__(self):
		self.arguments = arguments
		self.description = description

	def __call__(self, args):
		print(f'Homology search started at {timestamp()}.')

		# load configurations
		self.cfg = load_configs()

		# read and validate arguments
		self.args_wf(args)

		# read and validate input data
		self.input_wf()

		# perform homology search for each sample
		for sid, sample in sorted(self.data.items()):
			if 'done' in sample:
				continue
			prots = sample['prots']

			print(f'Batch homology search of {sid} started at {timestamp()}.')

			# collect sequences to search
			id2idx = {}
			seqs = []
			for i, prot in enumerate(prots):
				if 'hits' in prot:
					continue
				id_ = prot['id']
				seqs.append((id_, prot['seq']))
				id2idx[id_] = i
			print(f'Number of queries: {len(seqs)}.')

			# divide sequences into batches
			batches = ([seqs] if self.method == 'precomp' else
					   self.subset_seqs(seqs, self.queries, self.maxchars))

			# run search for each batch
			n = 0
			for batch in batches:

				# batch homology search
				print(getattr(self, 'pcmap', '') or '')  # was `print(self.pcmap)`: pcmap only exists on -m precomp, the bare attribute crashed every native-engine batch
				# todo: general fix:
				# temporary fix by replacing self.pcmap[sid] with self.pcmap.values()[0]
				res = self.search_wf(
					batch, list(self.pcmap.values())[0] if self.method == 'precomp'
					else None)

				# update taxIds of hits
				self.taxid_wf(res)

				# update taxonomic information of new taxIds
				self.taxinfo_wf(res)

				# perform taxonomy-based filtering
				self.taxfilt_wf(res)

				# update samples with search results
				indices = [id2idx[x[0]] for x in batch]
				self.update_search_results(prots, res, set(indices))

				# perform self-alignment
				seqs2a = [x for x in batch if 'score' not in prots[
					id2idx[x[0]]]]
				if seqs2a:
					for id_, score in self.selfaln_wf(seqs2a, res).items():
						prots[id2idx[id_]]['score'] = score

				# write search results to file
				with open(join(self.output, f'{sid}.tsv'), 'a') as f:
					self.write_search_results(f, prots, indices)

				n += len(batch)
				print(f'  {n} queries completed.')

			print(f'Batch homology search of {sid} ({len(prots)} proteins) '
				  f'ended at {timestamp()}.')

		# clean up
		if hasattr(self, 'mkdtemp'):
			rmtree(self.tmpdir)

		print(f'Batch homology search finished at {timestamp()}.')

	"""master workflows"""

	def args_wf(self, args):
		"""Workflow for validating and setting arguments.

		Parameters
		----------
		args : dict
			command-line arguments

		Notes
		-----
		Workflow:
		1. Load command-line arguments.
		2. Update arguments from configurations.
		3. Validate input and output directories.
		4. Determine search method and parameters.
		5. Determine taxonomy database and map.
		6. Determine self-alignment method and parameters.
		7. Determine remote fetch method and parameters.
		8. Print major settings if desired.
		"""
		# load arguments
		for key, val in vars(args).items():
			setattr(self, key, val)

		# check input directory and data
		if isfile(self.input):
			self.input_map = {file2id(self.input): self.input}
		elif isdir(self.input):
			self.input_map = {k: join(self.input, v) for k, v in id2file_map(
				self.input).items()}
			if len(self.input_map) == 0:
				raise ValueError(
					f'No input data are found under: {self.input}.')
		else:
			raise ValueError(
				f'Invalid input data file or directory: {self.input}.')

		# check / create output directory
		makedirs(self.output, exist_ok=True)
		self.prev_map = id2file_map(self.output, 'tsv')

		"""determine search strategy"""

		# load search parameters
		get_config(self, 'evalue', 'search.evalue', float)
		for key in 'method', 'minsize', 'maxseqs', 'identity', 'coverage':
			get_config(self, key, f'search.{key}')
		for key in 'diamond', 'blastp', 'blastdbcmd':
			get_config(self, key, f'program.{key}')

		if self.method not in {'auto', 'diamond', 'blast', 'remote',
							   'precomp'}:
			raise ValueError(f'Invalid search method: {self.method}.')

		# look for precomputed search results
		if self.method == 'precomp' and not self.precomp:
			raise ValueError('Must specify location of pre-computed search '
							 'results.')
		if self.precomp:
			if isfile(self.precomp):
				if len(self.input_map) > 1:
					raise ValueError('A directory of multiple pre-computed '
									 'search result files is needed.')
				self.pcmap = {file2id(self.precomp): self.precomp}
			elif isdir(self.precomp):
				self.pcmap = {k: join(self.precomp, v) for k, v in id2file_map(
					self.precomp).items()}
				if len(self.pcmap) == 0:
					raise ValueError('Cannot locate any pre-computed search '
									 f'results at: {self.precomp}.')
			else:
				raise ValueError('Invalid pre-computed search result file or '
								 f'directory: {self.precomp}.')
			if self.method == 'auto':
				self.method = 'precomp'

		# check local search executables and databases
		diamond_db = self.check_diamond()
		blast_db = self.check_blast()

		# choose a local search method if available, or do remote search
		if self.method == 'auto':
			if diamond_db:
				self.method = 'diamond'
				self.db = diamond_db
			elif blast_db:
				self.method = 'blast'
				self.db = blast_db
			else:
				self.method = 'remote'

		# load method-specific arguments
		for key in 'queries', 'maxchars', 'extrargs':
			get_config(self, key, f'{self.method}.{key}')

		# load remote search settings
		if self.method == 'remote':
			for key in 'db', 'algorithm', 'delay', 'timeout', 'entrez':
				get_config(self, key, f'remote.{key}')
			get_config(self, 'server', 'server.search')

		# determine number of threads
		if self.method in ('diamond', 'blast') and not self.threads:

			# use all available CPUs if threads is set to zero or left empty
			self.threads = cpu_count()

			# do single-threading if CPU count not working
			if self.threads is None:
				print('WARNING: Cannot determine number of CPUs. Will do '
					  'single-threading instead.')
				self.threads = 1

			# apply BLAST CPU number cap
			if self.method == 'blast' and self.threads > 8:
				print('WARNING: BLAST can only use a maximum of 8 CPUs.')
				self.threads = 8

		# check / create temporary directory
		if self.method in ('diamond', 'blast'):
			dir_ = self.tmpdir
			if not dir_:
				self.tmpdir = mkdtemp()
				setattr(self, 'mkdtemp', True)  # mark for cleanup
			elif not isdir(dir_):
				raise ValueError(f'Invalid temporary directory: {dir_}.')

		"""determine taxonomy database and filters"""

		# initialize protein-to-taxId map
		self.prot2tid = {}

		# assign taxonomy database
		for key in 'taxdump', 'taxmap':
			get_config(self, key, f'database.{key}')

		if self.method != 'remote':

			# check local taxonomy database
			if not self.taxdump:
				print('WARNING: Local taxonomy database is not specified. '
					  'Will try to retrieve taxonomic information from remote '
					  'server.')
			elif not isdir(self.taxdump):
				raise ValueError(
					f'Invalid taxonomy database directory: {self.taxdump}.')
			else:
				for fname in 'names.dmp', 'nodes.dmp':
					if not isfile(join(self.taxdump, fname)):
						raise ValueError(
							f'Taxonomy database file {fname} is not found.')

			# check local taxonomy map
			if self.taxmap and not isfile(self.taxmap):
				raise ValueError(
					f'Invalid protein-to-taxId map: {self.taxmap}.')

		# load taxonomic filters and convert to lists
		for key in 'include', 'exclude', 'block':
			attr = f'tax_{key}'
			get_config(self, attr, f'taxonomy.{key}')
			setattr(self, attr, list_from_param(getattr(self, attr)))

		# load taxonomy switches
		for key in 'unique', 'unirank', 'capital', 'latin':
			get_config(self, f'tax_{key}', f'taxonomy.{key}')

		"""determine self-alignment strategy"""

		# load configurations
		get_config(self, 'aln_method', 'search.selfaln')
		get_config(self, 'aln_server', 'server.selfaln')
		met_ = self.aln_method
		if met_ not in ('auto', 'native', 'fast', 'lookup', 'precomp'):
			raise ValueError(f'Invalid self-alignment method: {met_}.')

		# look for precomputed self-alignment results
		if met_ == 'precomp' and not self.aln_precomp:
			raise ValueError('Must specify location of pre-computed self-'
							 'alignment scores.')
		pre_ = self.aln_precomp
		if pre_:
			if isfile(pre_):
				if len(self.input_map) > 1:
					raise ValueError('A directory of multiple pre-computed '
									 'self-alignment result files is needed.')
				self.aln_pcmap = {file2id(pre_): pre_}
			elif isdir(pre_):
				self.aln_pcmap = {k: join(pre_, v) for k, v in
								  id2file_map(pre_).items()}
				if len(self.aln_pcmap) == 0:
					raise ValueError('Cannot locate any pre-computed self-'
									 f'alignment results at: {pre_}.')
			else:
				raise ValueError('Invalid pre-computed self-alignment result '
								 f'file or directory: {pre_}.')
			if met_ == 'auto':
				self.aln_method = 'precomp'

		# use the same search method for self-alignment, otherwise use fast
		if met_ in ('auto', 'native'):
			self.aln_method = 'fast' if self.method == 'precomp' else 'native'

		"""determine fetch strategy"""

		# load configurations
		get_config(self, 'fetch_server', 'server.fetch')
		for key in 'enable', 'queries', 'retries', 'delay', 'timeout':
			get_config(self, f'fetch_{key}', f'fetch.{key}')

		# determine remote or local fetching
		if self.fetch_enable == 'auto':
			self.fetch_enable = 'yes' if (
				self.method == 'remote' or not self.taxdump) else 'no'

		"""final steps"""

		# convert boolean values
		for key in 'tax_unique', 'tax_capital', 'tax_latin':
			setattr(self, key, arg2bool(getattr(self, key, None)))

		# convert fractions to percentages
		for metric in 'identity', 'coverage':
			val = getattr(self, metric)
			if val and val < 1:
				setattr(self, metric, val * 100)

		# print major settings
		print('Settings:')
		print(f'  Search method: {self.method}.')
		print(f'  Self-alignment method: {self.aln_method}.')
		print(f'  Remote fetch enabled: {self.fetch_enable}.')

	def check_diamond(self):
		"""Check if DIAMOND is available.

		Returns
		-------
		str or None
			valid path to DIAMOND database, or None if not available

		Raises
		------
		ValueError
			If settings conflict.
		"""
		if self.method in ('diamond', 'auto'):
			if not self.diamond:
				self.diamond = 'diamond'
			if which(self.diamond):
				try:
					db = self.db or self.cfg['database']['diamond']
				except KeyError:
					pass
				if db:
					if isfile(db) or isfile(f'{db}.dmnd'):
						return db
					elif self.method == 'diamond':
						raise ValueError(f'Invalid DIAMOND database: {db}.')
				elif self.method == 'diamond':
					raise ValueError(
						'A protein database is required for DIAMOND search.')
			elif self.method == 'diamond':
				raise ValueError(
					f'Invalid diamond executable: {self.diamond}.')

	def check_blast(self):
		"""Check if BLAST is available.

		Returns
		-------
		str or None
			valid path to BLAST database, or None if not available

		Raises
		------
		ValueError
			If settings conflict.
		"""
		if self.method in ('blast', 'auto'):
			if not self.blastp:
				self.blastp = 'blastp'
			if which(self.blastp):
				try:
					db = self.db or self.cfg['database']['blast']
				except KeyError:
					pass
				if db:
					if all(isfile(f'{db}.{x}') for x in ('phr', 'pin', 'psq')):
						return db
					elif self.method == 'blast':
						raise ValueError(f'Invalid BLAST database: {db}.')
				elif self.method == 'blast':
					raise ValueError(
						'A protein database is required for BLAST search.')
			elif self.method == 'blast':
				raise ValueError(f'Invalid blastp executable: {self.blastp}.')

	def input_wf(self):
		"""Master workflow for processing input data.

		Notes
		-----
		Workflow:
		1. Read proteins of each protein set.
		2. Process search results from previous runs.
		3. Import precomputed self-alignment scores.
		4. Fetch sequences for proteins to be searched.
		5. Drop sequences shorter than threshold.
		6. Read or initiate taxonomy database.
		7. Read protein-to-taxId map.
		"""
		# initiate data structure
		self.data = {}

		# read proteins of each sample
		print('Reading input proteins...')
		nprot = 0
		for id_, fname in self.input_map.items():
			prots = read_input_prots(fname)
			n = len(prots)
			if n == 0:
				raise ValueError(f'No protein entries are found for {id_}.')
			print(f'  {id_}: {n} proteins.')
			self.data[id_] = {'prots': prots}
			nprot += n
		print(f'Done. Read {nprot} proteins from {len(self.data)} samples.')

		# process search results from previous runs
		ndone = 0
		if self.prev_map:
			print('Processing search results from previous runs...')
			for id_ in self.data:
				if id_ not in self.prev_map:
					continue
				prots_ = self.data[id_]['prots']
				n = len(prots_)
				m = len(self.parse_prev_results(join(
					self.output, self.prev_map[id_]), prots_))
				if m == n:
					self.data[id_]['done'] = True
				ndone += m
			print(f'Done. Found results for {ndone} proteins, remaining '
				  f'{nprot - ndone} to search.')

		# check if search is already completed
		if (ndone == nprot):
			return

		# import precomputed self-alignment scores
		if self.aln_method == 'precomp' and hasattr(self, 'aln_pcmap'):
			n, m = 0, 0
			print('Importing precomputed self-alignment scores...', end='')
			for sid, file in self.aln_pcmap.items():

				# read scores
				scores = {}
				with open(file, 'r') as f:
					for line in f:
						line = line.rstrip('\r\n')
						if not line or line.startswith('#'):
							continue
						id_, score = line.split('\t')
						scores[id_] = float(score)

				# assign scores if not already
				for prot in self.data[sid]['prots']:
					if 'score' in prot:
						continue
					n += 1
					id_ = prot['id']
					try:
						prot['score'] = scores[id_]
						m += 1
					except KeyError:
						pass
			print(' done.')
			print(f'  Imported scores for {n} proteins.')
			dif = n - m
			if dif > 0:
				raise ValueError(f'Missing scores for {dif} proteins.')

		# fetch sequences for unsearched proteins
		seqs2q = self.check_missing_seqs(self.data)
		n = len(seqs2q)
		if n > 0:
			print(f'Sequences of {n} proteins are to be retrieved.')
			if self.method == 'blast':
				print('Fetching sequences from local BLAST database...',
					  end='')
				seqs = self.blast_seqinfo(seqs2q)
				n = self.update_prot_seqs(seqs)
				print(' done.')
				print(f'  Obtained sequences of {n} proteins.')
				seqs2q = self.check_missing_seqs(self.data)
				n = len(seqs2q)
				if n > 0:
					print(f'  Remaining {n} proteins.')
			if n > 0 and self.fetch_enable == 'yes':
				print(f'Fetching {n} sequences from remote server...',
					  flush=True)
				seqs = self.remote_seqinfo(seqs2q)
				n = self.update_prot_seqs(seqs)
				print(f'Done. Obtained sequences of {n} proteins.')
				seqs2q = self.check_missing_seqs(self.data)
				n = len(seqs2q)
			if n > 0:
				raise ValueError(f'  Cannot obtain sequences of {n} proteins.')

		# drop short sequences
		if self.minsize:
			print(f'Dropping sequences shorter than {self.minsize} aa...',
				  end='')
			for sid, sample in self.data.items():
				for i in reversed(range(len(sample['prots']))):
					if len(sample['prots'][i]['seq']) < self.minsize:
						del sample['prots'][i]
			print(' done.')

		# read or initiate taxonomy database
		# read external taxdump
		if self.taxdump is not None:
			print('Reading local taxonomy database...', end='')
			self.taxdump = read_taxdump(self.taxdump)
			print(' done.')
			print(f'  Read {len(self.taxdump)} taxa.')

		# read taxdump generated by previous runs
		elif (isfile(join(self.output, 'names.dmp')) and
				isfile(join(self.output, 'nodes.dmp'))):
			print('Reading custom taxonomy database...', end='')
			self.taxdump = read_taxdump(self.output)
			print(' done.')
			print(f'  Read {len(self.taxdump)} taxa.')

		# build taxdump from scratch
		else:
			print('Initiating custom taxonomy database...', end='')
			self.taxdump = {'1': {
				'name': 'root', 'parent': '1', 'rank': 'no rank'}}
			self.update_dmp_files(['1'])
			print(' done.')

		# record invalid taxIds to save compute
		self.badtaxids = set()

	def search_wf(self, seqs, file=None):
		"""Master workflow for batch homology search.

		Parameters
		----------
		seqs : list of tuple
			query sequences (Id, sequence)
		file : str, optional
			file of precomputed results

		Returns
		-------
		dict
			sequence Id to hit table

		Notes
		-----
		Workflow:
		1. Generate an Id-to-length map.
		2. Import precomputed search results (if method = precomp).
		3. Call choice of method (remote, blast or diamond) to search.
		"""
		# generate an Id-to-length map
		lenmap = {x[0]: len(x[1]) for x in seqs}

		# import pre-computed search results
		if self.method == 'precomp':
			print('Importing pre-computed search results...', end='')
			res = self.parse_hit_table(file, lenmap)
			print(' done.')
			print(f'  Found hits for {len(res)} proteins.')

		# run de novo search
		elif self.method == 'remote':
			res = self.remote_search(seqs)
			sleep(self.delay)
		elif self.method == 'blast':
			res = self.blast_search(seqs)
		elif self.method == 'diamond':
			res = self.diamond_search(seqs)

		return res

	def taxid_wf(self, prots):
		"""Master workflow for associating hits with taxIds.

		Parameters
		----------
		prots : dict of list
			proteins (search results)

		Notes
		-----
		Workflow:
		1. Update taxmap with taxIds directly available from hit tables.
		2. Get taxIds for sequence Ids without them:
		2.1. Query external taxon map, if available.
		2.2. Do local fetch (from a BLAST database) if available.
		2.3. Do remote fetch (from NCBI server) if available.
		3. Delete hits without associated taxIds.
		"""
		added_taxids = set()

		# update taxon map with taxIds already in hit tables
		ids2q, added = self.update_hit_taxids(prots)
		added_taxids.update(added)

		# attempt to look up taxIds from external taxon map
		if ids2q and self.method != 'remote' and self.taxmap is not None:

			# load taxon map on first use (slow and memory-hungry)
			if isinstance(self.taxmap, str):
				print('Reading protein-to-TaxID map (WARNING: may be slow and '
					  'memory-hungry)...', end='', flush=True)
				self.taxmap = read_prot2taxid(self.taxmap)
				print(' done.')
				print(f'  Read {len(self.taxmap)} records.')

			ids2q, added = self.update_hit_taxids(prots, self.taxmap)
			added_taxids.update(added)

		# attempt to look up taxIds from local BLAST database
		if (ids2q and self.method in ('blast', 'precomp') and self.db
				and self.blastdbcmd):
			newmap = {x[0]: x[1] for x in self.blast_seqinfo(ids2q)}
			ids2q, added = self.update_hit_taxids(prots, newmap)
			added_taxids.update(added)

		# attempt to look up taxIds from remote server
		if ids2q and self.fetch_enable == 'yes':
			print(f'Fetching taxIds of {len(ids2q)} sequences from remote '
				  'server...', flush=True)
			newmap = {x[0]: x[1] for x in self.remote_seqinfo(ids2q)}
			print(f'Done. Obtained taxIds of {len(newmap)} sequences.')
			ids2q, added = self.update_hit_taxids(prots, newmap)
			added_taxids.update(added)

		# drop hits whose taxIds cannot be identified
		n = len(ids2q)
		if n > 0:
			print(f'WARNING: Cannot obtain taxIds for {n} sequences. These '
				  'hits will be dropped.')
			for hits in prots.values():
				for i in reversed(range(len(hits))):
					if hits[i]['id'] in ids2q:
						del hits[i]

	def taxinfo_wf(self, prots):
		"""Master workflow for associating hits with taxonomic information.

		Parameters
		----------
		prots : dict of list
			proteins (search results)

		Notes
		-----
		Workflow:
		1. Obtain a list of taxIds represented by hits.
		2. List taxIds that are missing in current taxonomy database.
		3. Get taxonomic information for taxIds by remote fetch, if available.
		4. Append new taxonomic information to dump files, if available.
		5. Delete hits whose taxIds are not associated with information.
		"""
		# list taxIds without information
		tids2q = set()
		for prot, hits in prots.items():
			for hit in hits:
				tid = hit['taxid']
				if tid not in self.taxdump:
					tids2q.add(tid)
		tids2q = sorted(tids2q)

		# retrieve information for these taxIds
		if tids2q and self.fetch_enable == 'yes':
			print(f'Fetching {len(tids2q)} taxIds and ancestors from remote '
				  'server...', flush=True)
			xml = self.remote_taxinfo(tids2q)
			added = self.parse_taxonomy_xml(xml)
			print(f'Done. Obtained taxonomy of {len(tids2q)} taxIds.')
			self.update_dmp_files(added)
			tids2q = [x for x in tids2q if x not in added]

		# drop taxIds whose information cannot be obtained
		if tids2q:
			print(f'WARNING: Cannot obtain information of {len(tids2q)} '
				  'taxIds. These hits will be dropped.')
			tids2q = set(tids2q)
			for hits in prots.values():
				for i in reversed(range(len(hits))):
					if hits[i]['taxid'] in tids2q:
						del hits[i]
			self.badtaxids.update(tid)

	def taxfilt_wf(self, prots):
		"""Workflow for filtering hits by taxonomy.

		Parameters
		----------
		prots : dict of list of dict
			proteins (search results)

		Notes
		-----
		Workflow:
		1. Bottom-up filtering (delete in place)
		1.1. Delete taxIds already marked as bad.
		1.2. Delete taxIds whose ancestors are not in the "include" list.
		1.3. Delete taxIds, any of whose ancestors is in the "exclude" list.
		1.4. Delete empty taxon names.
		1.5. Delete taxon names that are not capitalized.
		1.6. Delete taxon names in which any word is in the "block" list.
		2. Bottom-up filtering (mark and batch delete afterwards)
		2.1. Mark taxIds that already appeared in hit table for deletion.
		2.2. Mark taxIds whose ancestor at given rank already appeared in hit
			table for deletion.
		"""
		# filtering of taxIds and taxon names
		for id_, hits in prots.items():

			# bottom-up filtering by independent criteria
			for i in reversed(range(len(hits))):
				todel = False
				tid = hits[i]['taxid']
				taxon = self.taxdump[tid]['name']
				if tid in self.badtaxids:
					todel = True
				elif (self.tax_include and not
					  is_ancestral(tid, self.tax_include, self.taxdump)):
					todel = True
				elif (self.tax_exclude and
					  is_ancestral(tid, self.tax_exclude, self.taxdump)):
					todel = True
				elif taxon == '':
					todel = True
				elif self.tax_capital and not is_capital(taxon):
					todel = True
				elif self.tax_block and contain_words(taxon, self.tax_block):
					todel = True
				elif self.tax_latin:
					tid_ = taxid_at_rank(tid, 'species', self.taxdump)
					if not tid_ or not is_latin(tid_):
						todel = True
				if todel:
					del hits[i]
					self.badtaxids.add(tid)

			# top-down filtering by sorting-based criteria
			todels = []
			used = set()
			used_at_rank = set()
			for i in range(len(hits)):
				tid = hits[i]['taxid']
				if self.tax_unique:
					if tid in used:
						todels.append(i)
						continue
					else:
						used.add(tid)
				if self.tax_unirank:
					tid_ = taxid_at_rank(tid, self.tax_unirank, self.taxdump)
					if tid_ and tid_ in used_at_rank:
						todels.append(i)
					else:
						used_at_rank.add(tid_)
			for i in reversed(todels):
				del hits[i]

	def selfaln_wf(self, seqs, prots=None):
		"""Master workflow for protein sequence self-alignment.

		Parameters
		----------
		seqs : list of tuple
			query sequences (Id, sequence)
		prots : dict of list of dict, optional
			hit tables, only relevant when amet = lookup

		Returns
		-------
		dict
			Id-to-score map

		Notes
		-----
		Workflow:
		1. If amet = lookup, just look up, and raise if any sequences don't
		have self-hits.
		2. If amet = fast, run fast built-in algorithm on each sequence.
		3. If amet = native, find the corresponding search method.
		4. Run self-alignment in batches only when method = precomp, otherwise
		the query sequences are already subsetted into batches.
		5. If some sequences don't have self hits in batch self-alignments,
		try to get them via single self-alignments.
		6. If some sequences still don't have self hits, do built-in algorithm,
		but note that the output may be slightly different from others.
		7. If some sequences still don't have self hits, raise.
		"""
		res = []

		# just look up (will fail if some are not found)
		if self.aln_method == 'lookup':
			res = self.lookup_selfaln(seqs, prots)

		# use built-in algorithm
		elif self.aln_method == 'fast':
			for id_, seq in seqs:
				bitscore, evalue = self.fast_selfaln(seq)
				res.append((id_, bitscore, evalue))

		# use the same search method for self-alignment
		elif self.aln_method == 'native':

			# divide sequences into batches, when search results are
			# precomputed
			batches = ([seqs] if self.method != 'precomp' else
					   self.subset_seqs(seqs, self.queries, self.maxchars))

			# do self-alignments in batches to save compute
			for batch in batches:
				res_ = []

				# call search method
				if self.method == 'remote':
					res_ = self.remote_selfaln(batch)
					sleep(self.delay)
				elif self.method == 'blast':
					res_ = self.blast_selfaln(batch)
				elif self.method == 'diamond':
					res_ = self.diamond_selfaln(batch)

				# merge results
				res += res_

		# if some hits are not found, do single alignments
		left = set([x[0] for x in seqs]) - set([x[0] for x in res])
		if left:
			print('WARNING: The following sequences cannot be self-aligned '
				  'in a batch. Do individual alignments instead.')
			print('  ' + ', '.join(left))
			for id_, seq in seqs:
				if id_ not in left:
					continue
				res_ = None

				# call search method
				if self.method == 'remote':
					res_ = self.remote_selfaln([(id_, seq)])
					sleep(self.delay)
				elif self.method == 'blast':
					res_ = self.blast_selfaln([(id_, seq)])
				elif self.method == 'diamond':
					res_ = self.diamond_selfaln([(id_, seq)])

				# if failed, do built-in alignment
				if not res_:
					print(f'WARNING: Sequence {id_} cannot be self-aligned '
						  'using the native method. Do fast alignment '
						  'instead.')
					bitscore, evalue = self.fast_selfaln(seq)
					res_ = [(id_, bitscore, evalue)]

				# merge results
				res += res_

			# check if all sequences have results
			left = set([x[0] for x in seqs]) - set([x[0] for x in res])
			if left:
				raise ValueError('Cannot calculate self-alignment metrics for '
								 'the following sequences:\n  ' + ', '.join(
									 sorted(left)))
		return {x[0]: x[1] for x in res}

	"""input/output functions"""

	@staticmethod
	def subset_seqs(seqs, queries=None, maxchars=None):
		"""Generate subsets of sequences based on cutoffs.

		Parameters
		----------
		seqs : list of tuple
			sequences to subset (id, sequence)
		queries : int, optional
			number of query sequences per subset
		maxchars : int, optional
			maximum total length of query sequences per subset

		Returns
		-------
		list of list of tuple
			subsets

		Raises
		------
		ValueError
			If any sequence exceeds maxchars.
		"""
		if not maxchars:

			# no subsetting
			if not queries:
				return [seqs]

			# subsetting only by queries
			subsets = []
			for i in range(0, len(seqs), queries):
				subsets.append(seqs[i:i + queries])
			return subsets

		# subsetting by maxchars, and by queries if applicable
		subsets = [[]]
		cquery, cchars = 0, 0
		for id_, seq in seqs:
			chars = len(seq)
			if chars > maxchars:
				raise ValueError(f'Sequence {id_} exceeds maximum allowed '
								 f'length {maxchars} for search.')
			if cchars + chars > maxchars or queries == cquery > 0:
				subsets.append([])
				cquery, cchars = 0, 0
			subsets[-1].append((id_, seq))
			cquery += 1
			cchars += chars
		return subsets

	def update_search_results(self, prots, res, indices=set()):
		"""Update proteins with new search results.

		Parameters
		----------
		prots : list of dict
			proteins to update
		res : dict
			search results
		indices : set of int, optional
			indices of proteins to be updated
			if omitted, only proteins with hits will be updated
		"""
		for i, prot in enumerate(prots):
			if 'hits' in prot:
				continue
			if indices and i not in indices:
				continue
			id_ = prot['id']
			if id_ in res:
				prot['hits'] = []
				n = 0
				for hit in res[id_]:
					prot['hits'].append(hit)
					n += 1
					if self.maxhits and n == self.maxhits:
						break
			elif indices and i in indices:
				prot['hits'] = []

	@staticmethod
	def write_search_results(f, prots, indices=None):
		"""Write search results to a file.

		Parameters
		----------
		f : file handle
			file to write to (in append mode)
		prots : array of hash
			protein set
		indices : list of int, optional
			limit to these proteins
		"""
		for i in indices if indices else range(len(prots)):
			prot = prots[i]
			f.write(f'# ID: {prot["id"]}\n')
			f.write(f'# Length: {len(prot["seq"])}\n')
			f.write(f'# Product: {prot["product"]}\n')
			f.write(f'# Score: {prot["score"]}\n')
			f.write(f'# Hits: {len(prot["hits"])}\n')
			for hit in prot['hits']:
				f.write('\t'.join([hit[x] for x in (
					'id', 'identity', 'evalue', 'score', 'coverage',
					'taxid')]) + '\n')

	@staticmethod
	def parse_prev_results(fp, prots):
		"""Parse previous search results.

		Parameters
		----------
		file : str
			file containing search results
		prots : list
			protein records

		Returns
		-------
		list of str
			completed protein Ids
		"""
		done = []
		with open(fp, 'r') as f:
			for line in f:
				if line.startswith('# ID: '):
					done.append(line[6:].rstrip('\r\n'))
		doneset = set(done)
		for prot in prots:
			if prot['id'] in doneset:
				prot['score'] = 0
				prot['hits'] = []
		return done

	@staticmethod
	def check_missing_seqs(data):
		"""Get a list of proteins whose sequences remain to be retrieved.

		Parameters
		----------
		data : dict
			protein sets

		Returns
		-------
		list of str
			Ids of proteins without sequences
		"""
		res = set()
		for sid, sample in data.items():
			if 'done' in sample:
				continue
			for prot in sample['prots']:
				if not prot['seq'] and 'hits' not in prot:
					res.add(prot['id'])
		return sorted(res)

	def update_dmp_files(self, ids):
		"""Write added taxonomic information to custom taxdump files.

		Parameters
		----------
		ids : list of str
			added taxIds

		Notes
		-----
		Taxonomic information will be appended to nodes.dmp and names.dmp in
		the working directory.
		"""
		fo = open(join(self.output, 'nodes.dmp'), 'a')
		fa = open(join(self.output, 'names.dmp'), 'a')
		for id_ in sorted(ids, key=int):
			fo.write('\t|\t'.join((
				id_, self.taxdump[id_]['parent'],
				self.taxdump[id_]['rank'])) + '\t|\n')
			fa.write('\t|\t'.join((
				id_, self.taxdump[id_]['name'])) + '\t|\n')
		fo.close()
		fa.close()

	"""sequence query functions"""

	def blast_seqinfo(self, ids):
		"""Retrieve information of given sequence Ids from local BLAST database.

		Parameters
		----------
		ids : list of str
			query sequence Ids

		Returns
		-------
		list of tuple
			(id, taxid, product, sequence)

		Notes
		-----
		When making database (using makeblastdb), one should do -parse_seqids
		to enable search by name (instead of sequence) and -taxid_map with a
		seqId-to-taxId map to enable taxId query.
		"""
		# run blastdbcmd
		# fields: accession, taxid, sequence, title
		cmd = ' '.join((
			self.blastdbcmd,
			'-db', self.db,
			'-entry', ','.join(ids),
			'-outfmt', '"%a %T %s %t"'))
		out = run_command(cmd)[1]

		# parse output
		res = []
		header = True
		for line in out:
			# catch invalid database error
			if header:
				# letter case is dependent on BLAST version
				if line.lower().startswith('blast database error'):
					raise ValueError(f'Invalid BLAST database: {self.db}.')
				header = False

			# if one sequence Id is not found, program will print:
			#   Error: [blastdbcmd] Entry not found: NP_123456.1
			# if none of sequence Ids are found, program will print:
			#   Error: [blastdbcmd] Entry or entries not found in BLAST
			#   database
			if (line.startswith('Error') or 'not found' in line or
					line.startswith('Please refer to')):
				continue

			# limit to 4 partitions because title contains spaces
			x = line.split(None, 3)

			# if database was not compiled with -taxid_map, taxIds will be 0
			if x[1] in ('0', 'N/A'):
				x[1] = ''

			# title will be empty if -parse_seqids was not triggered
			if len(x) == 3:
				x.append('')

			# parse title to get product
			else:
				x[3] = get_product(x[3])

			res.append((x[0], x[1], x[3], x[2]))
		return res

	def remote_seqinfo(self, ids):
		"""Retrieve information of given sequence Ids from remote server.

		Parameters
		----------
		ids : list of str
			query sequence Ids (e.g., accessions)

		Returns
		-------
		list of tuple
			(id, taxid, product, sequence)

		Raises
		------
		ValueError
			All sequence Ids are invalid.
			Failed to retrieve info from server.
		"""
		return self.parse_fasta_xml(self.remote_fetches(
			ids, 'db=protein&rettype=fasta&retmode=xml&id={}'))

	def remote_fetches(self, ids, urlapi):
		"""Fetch information from remote server in batch mode

		Parameters
		----------
		ids : list of str
			query entries (e.g., accessions)
		urlapi : str
			URL API, with placeholder for query entries

		Returns
		-------
		str
			fetched information

		Notes
		-----
		The function dynamically determines the batch size, starting from a
		large number and reducing by half on every other retry. This is because
		the NCBI server is typically busy and frequently runs into the "502 Bad
		Gateway" issue. To resolve, one may subset queries and retry.
		"""
		cq = self.fetch_queries  # current number of queries
		cids = ids  # current list of query Ids
		res = ''
		while True:
			batches = [cids[i:i + cq] for i in range(0, len(cids), cq)]
			failed = []

			# batch fetch sequence information
			for batch in batches:
				try:
					res += self.remote_fetch(urlapi.format(','.join(batch)))
					print(f'  Fetched information of {len(batch)} entries.',
						  flush=True)
					sleep(self.fetch_delay)
				except ValueError:
					failed.extend(batch)

			# reduce batch size by half on each trial
			if failed and cq > 1:
				cids = failed
				cq = int(cq / 2)
				print('Retrying with smaller batch size...', flush=True)
			else:
				cids = []
				break
		if cids:
			print(f'WARNING: Cannot retrieve information of {len(cids)} '
				  'entries.')
		return res

	def remote_fetch(self, urlapi):
		"""Fetch information from remote server.

		Parameters
		----------
		urlapi : str
			URL API

		Returns
		-------
		str
			fetched information

		Raises
		------
		ValueError
			Fetch failed.
		"""
		url = f'{self.fetch_server}?{urlapi}'
		for i in range(self.fetch_retries):
			if i:
				print('Retrying...', end=' ', flush=True)
				sleep(self.fetch_delay)
			try:
				with urlopen(url, timeout=self.fetch_timeout) as response:
					return response.read().decode('utf-8')
			except (HTTPError, URLError) as e:
				print(f'{e.code} {e.reason}.', end=' ', flush=True)
		print('', flush=True)
		raise ValueError('Failed to fetch information from remote server.')

	def update_prot_seqs(self, seqs):
		"""Update protein sets with retrieved sequences.

		Parameters
		----------
		seqs : list of tuple
			protein sequences (id, taxid, product, sequence)

		Returns
		-------
		int
			number of proteins with sequence added

		Notes
		------
		Different protein sets may contain identical protein Ids.
		"""
		# hash proteins by id
		prots = {x[0]: (x[2], x[3]) for x in seqs}

		# if queries are accessions without version, NCBI will add version
		acc2ver = {}
		for id_, info in prots.items():
			acc_ = re.sub(r'\.\d+$', '', id_)
			acc2ver[acc_] = id_

		n = 0
		for sid, sample in self.data.items():
			for prot in sample['prots']:
				if prot['seq']:
					continue

				# try to match protein Id (considering de-versioned accession)
				id_ = prot['id']
				if id_ not in prots:
					try:
						id_ = acc2ver[id_]
					except KeyError:
						continue
				if id_ not in prots:
					continue

				# update protein information
				for i, key in enumerate(['product', 'seq']):
					if not prot[key]:
						prot[key] = prots[id_][i]
				n += 1
		return n

	@staticmethod
	def parse_fasta_xml(xml):
		"""Parse sequence information in FASTA/XML format retrieved from NCBI
		server.

		Parameters
		----------
		xml : str
			sequence information in XML format

		Returns
		-------
		list of str
			[id, taxid, product, sequence]

		Notes
		-----
		NCBI EFectch record type = TinySeq XML

		.. _NCBI RESTful API:
			https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_
			values_of__retmode_and/?report=objectonly
		.. _NCBI RESTful API example:
			https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein
			&rettype=fasta&retmode=xml&id=NP_454622.1,NP_230502.1,NP_384288.1
		"""
		seqs = []
		for m in re.finditer(r'<TSeq>(.+?)<\/TSeq>', xml, re.DOTALL):
			s_ = m.group(1)
			seq = []
			for key in 'accver', 'taxid', 'defline', 'sequence':
				m_ = re.search(r'<TSeq_%s>(.+)<\/TSeq_%s>' % (key, key), s_)
				seq.append(m_.group(1) if m_ else '')
			seq[2] = get_product(seq[2])
			seqs.append(seq)
		return seqs

	"""homology search functions"""

	def blast_search(self, seqs):
		"""Run BLAST to search a sequence against a database.

		Parameters
		----------
		seqs : list of tuple
			query sequences (id, sequence)

		Returns
		-------
		dict of list of dict
			hit table per query sequence

		Raises
		------
		ValueError
			If BLAST run fails.

		Notes
		-----
		- In ncbi-blast+ 2.7.1, the standard tabular format (-outfmt 6) is:
			qaccver saccver pident length mismatch gapopen qstart qend sstart
			send evalue bitscore
		- In older versions, fields 1 and 2 are qseqid and sseqid. The
			difference is that an sseqid may read "ref|NP_123456.1|" instead of
			"NP_123456.1".
		- staxids are ;-delimited, will be "N/A" if not found or the database
			does not contain taxIds.
		- The BLAST database should ideally be prepared as:
			makeblastdb -in seqs.faa -dbtype prot -out db -parse_seqids \
				-taxid map seq2taxid.txt
		- Unlike blastn, blastp does not have -perc_identity.
		"""
		tmpin = join(self.tmpdir, 'tmp.in')
		with open(tmpin, 'w') as f:
			write_fasta(seqs, f)
		cmd = [self.blastp,
			   '-query', tmpin,
			   '-db', self.db]
		args = {x: getattr(self, x, None) for x in (
			'evalue', 'coverage', 'maxseqs', 'threads', 'extrargs')}
		if args['evalue']:
			cmd.extend(['-evalue', str(args['evalue'])])
		if args['coverage']:
			cmd.extend(['-qcov_hsp_perc', str(args['coverage'])])
		if args['maxseqs']:
			cmd.extend(['-max_target_seqs', str(args['maxseqs'])])
		if args['threads'] is not None:
			cmd.extend(['-num_threads', str(args['threads'])])
		if args['extrargs']:
			cmd.append(args['extrargs'])
		cmd.append('-outfmt "6 qaccver saccver pident evalue bitscore qcovhsp'
				   ' staxids"')
		ec, out = run_command(' '.join(cmd))
		remove(tmpin)
		if ec:
			raise ValueError(f'blastp failed with error code {ec}.')
		return self.parse_def_table(out)

	def diamond_search(self, seqs):
		"""Run DIAMOND to search a sequence against a database.

		Parameters
		----------
		seqs : list of tuple
			query sequences (id, sequence)

		Returns
		-------
		dict of list of dict
			hit table per query sequence

		Raises
		------
		ValueError
			If DIAMOND run fails.

		Notes
		-----
		The DIAMOND database should ideally be prepared as:
			diamond makedb --in seqs.faa --db db \
				--taxonmap prot.accession2taxid
		"""
		tmpin = join(self.tmpdir, 'tmp.in')
		with open(tmpin, 'w') as f:
			write_fasta(seqs, f)
		cmd = [self.diamond, 'blastp',
			   '--query', tmpin,
			   '--db', self.db,
			   '--threads', str(self.threads),
			   '--tmpdir', self.tmpdir]
		args = {x: getattr(self, x, None) for x in (
			'evalue', 'identity', 'coverage', 'maxseqs', 'extrargs')}
		if args['evalue']:
			cmd.extend(['--evalue', str(args['evalue'])])
		if args['identity']:
			cmd.extend(['--id', str(args['identity'])])
		if args['coverage']:
			cmd.extend([' --query-cover', str(args['coverage'])])
		if args['maxseqs']:
			cmd.extend(['--max-target-seqs', str(args['maxseqs'])])
		if args['extrargs']:
			cmd.append(args['extrargs'])
		cmd.extend(['--outfmt',
					'6 qseqid sseqid pident evalue bitscore qcovhsp staxids'])
		ec, out = run_command(' '.join(cmd), merge=False)
		remove(tmpin)
		if ec:
			raise ValueError(f'diamond failed with error code {ec}.')
		return self.parse_def_table(out)

	def remote_search(self, seqs):
		"""Perform BLAST search through a remote server.

		Parameters
		----------
		seqs : list of tuple
			query sequences (id, sequence)

		Returns
		-------
		dict of list of dict
			hit table per query sequence

		Notes
		-----
		As of 2021, use of RESTful APIs to conduct BLAST searches is
		deprecated by NCBI.

		.. _NCBI's official reference of RESTful APIs:
			https://ncbi.github.io/blast-cloud/dev/using-url-api.html
		.. _NCBI's official sample Perl script:
			https://blast.ncbi.nlm.nih.gov/docs/web_blast.pl
		.. _NCBI has restrictions on the frequency and bandwidth of remote
			BLAST searches. See this page:
			https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=
			BlastDocs&DOC_TYPE=DeveloperInfo
		.. _Instead, NCBI recommends setting up custom BLAST servers. See:
			https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&
			DOC_TYPE=CloudBlast
		"""
		print(f'Submitting {len(seqs)} queries for search.', end='',
			  flush=True)
		data = [('CMD', 'Put'),
				('PROGRAM', 'blastp'),
				('DATABASE', self.db),
				('WORD_SIZE', 6)]
		if self.algorithm:
			data.append(('BLAST_PROGRAMS', self.algorithm))
		if self.evalue:
			data.append(('EXPECT', self.evalue))
		if self.maxseqs:
			data.append(('MAX_NUM_SEQ', self.maxseqs))
		if self.entrez:
			data.append(('EQ_TEXT', self.entrez))
		query = ''.join([f'>{id_}\n{seq}\n' for id_, seq in seqs])
		data.append(('QUERY', query))
		data = urlencode(data)
		if self.extrargs:
			data += '&' + self.extrargs.lstrip('&')
		data = data.encode()

		trial = 0
		while True:
			if trial:
				if trial == (self.retries or 0) + 1:
					raise ValueError(
						f'Remote search failed after {trial} trials.')
				print(f'Retry {trial} times.', end='', flush=True)
				sleep(self.delay)
			trial += 1

			# get request Id
			req = Request(self.server, data=data)
			with urlopen(req) as response:
				res = response.read().decode('utf-8')
			m = re.search(r'^    RID = (.*$)', res, re.MULTILINE)
			if not m:
				print('WARNING: Failed to obtain RID.')
				continue
			rid = m.group(1)
			print(f' RID: {rid}.', end='', flush=True)
			sleep(10)

			# check status
			data_ = urlencode([('CMD', 'Get'),
							   ('FORMAT_OBJECT', 'SearchInfo'),
							   ('RID', rid)]).encode()
			req_ = Request(self.server, data=data_)
			starttime = time()
			success = False
			while True:
				with urlopen(req_) as response:
					res = response.read().decode('utf-8')
				m = re.search(r'\s+Status=(.+)', res, re.MULTILINE)
				if not m:
					print('WARNING: Failed to retrieve remote search status.')
					break
				status = m.group(1)
				if status == 'WAITING':
					if time() - starttime > self.timeout:
						print('WARNING: Remote search timeout.')
						break
					print('.', end='', flush=True)
					sleep(self.delay)
					continue
				elif status in ('FAILED', 'UNKNOWN'):
					print('WARNING: Remote search failed.')
					break
				elif status == 'READY':
					if 'ThereAreHits=yes' not in res:
						print('WARNING: Remote search returned no result.')
						break
					print(' Success.', end='', flush=True)
					success = True
					break
				else:
					print(f'WARNING: Unknown remote search status: {status}.')
					break
			if not success:
				continue
			sleep(10)

			# retrieve result
			data_ = [('CMD', 'Get'),
					 ('ALIGNMENT_VIEW', 'Tabular'),
					 ('FORMAT_TYPE', 'Text'),
					 ('RID', rid)]
			if self.maxseqs:
				data_.extend([('MAX_NUM_SEQ',  self.maxseqs),
							  ('DESCRIPTIONS', self.maxseqs)])
			data_ = urlencode(data_).encode()
			req_ = Request(self.server, data=data_)
			with urlopen(req_) as response:
				res = response.read().decode('utf-8')
			if '# blastp' not in res or '# Query: ' not in res:
				print('WARNING: Invalid format of remote search results.')
				continue
			print(' Results retrieved.')
			break

		# fields (as of 2018): query acc.ver, subject acc.ver, % identity,
		# alignment length, mismatches, gap opens, q. start, q. end, s. start,
		# s. end, evalue, bit score, % positives
		m = re.search(r'<PRE>(.+?)<\/PRE>', res, re.DOTALL)
		out = m.group(1).splitlines()
		lenmap = {id_: len(seq) for id_, seq in seqs}
		return self.parse_m8_table(out, lenmap)

	def parse_hit_table(self, file, lenmap=None):
		"""Determine hit table type and call corresponding parser.

		Parameters
		----------
		file : str
			hit table file
		lenmap : dict, optional
			map of sequence Ids to lengths (only needed for m8)

		Returns
		-------
		list of dict
			hit table
		"""
		ism8 = None
		lines = []
		with open(file, 'r') as f:
			for line in f:
				line = line.rstrip('\r\n')
				if line and not line.startswith('#'):
					lines.append(line)
					if ism8 is None:
						x = line.split('\t')
						ism8 = len(x) > 8
		return (self.parse_m8_table(lines, lenmap) if ism8 else
				self.parse_def_table(lines))

	def parse_def_table(self, lines):
		"""Parse search results in default tabular format.

		Parameters
		----------
		lines : list of str
			search result in default tabular format
			fields: qseqid sseqid pident evalue bitscore qcovhsp staxids

		Returns
		-------
		dict of list of dict
			hits per query
		"""
		res = {}
		ths = {x: getattr(self, x, 0) for x in (
			'evalue', 'identity', 'coverage', 'maxhits')}
		for line in lines:
			line = line.rstrip('\r\n')
			if not line or line.startswith('#'):
				continue
			x = line.split('\t')

			try:
				# filter by thresholds
				if ths['evalue']:
					if x[3] != '*' and ths['evalue'] < float(x[3]):
						continue
				if ths['identity']:
					if x[2] != '*' and ths['identity'] > float(x[2]):
						continue
				if ths['coverage']:
					if x[5] != '*' and ths['coverage'] > float(x[5]):
						continue
				
				# pass if maximum targets reached
				if ths['maxhits']:
					if x[0] in res and ths['maxhits'] == len(res[x[0]]):
						continue

				# add hit to list
				res.setdefault(x[0], []).append({
					'id': seqid2accver(x[1]), 'identity': x[2], 'evalue': x[3],
					'score': x[4], 'coverage': x[5], 'taxid': x[6] if x[6] not
					in {'', 'N/A', '0'} else ''})

			except Exception as e:
				print(e)
		return res

	def parse_m8_table(self, lines, lenmap):
		"""Parse search results in BLAST's standard tabular format (m8).

		Parameters
		----------
		lines : list of str
			search result in BLAST m8 tabular format
			fields: qseqid sseqid pident length mismatch gapopen qstart qend
					sstart send evalue bitscore
		lenmap : dict
			map of sequence Ids to lengths (needed for calculating coverage)

		Returns
		-------
		list of dict
			hit table

		Raises
		------
		ValueError
			Query Id not found in length map.
		"""
		res = {}
		ths = {x: getattr(self, x, 0) for x in (
			'evalue', 'identity', 'coverage', 'maxhits')}
		for line in lines:
			line = line.rstrip('\r\n')
			if not line or line.startswith('#'):
				continue
			x = line.split('\t')

			# calculate coverage
			if x[0] not in lenmap:
				raise ValueError(f'Invalid query sequence Id: {x[0]}.')
			try:
				cov = (int(x[7]) - int(x[6]) + 1) / lenmap[x[0]] * 100
			except ValueError:
				cov = 0

			# filter by thresholds
			if ths['evalue']:
				if x[10] != '*' and ths['evalue'] < float(x[10]):
					continue
			if ths['identity']:
				if x[2] != '*' and ths['identity'] > float(x[2]):
					continue
			if ths['coverage']:
				if cov and ths['coverage'] > cov:
					continue

			# pass if maximum targets reached
			if ths['maxhits']:
				if x[0] in res and ths['maxhits'] == len(res[x[0]]):
					continue

			# add hit to list
			res.setdefault(x[0], []).append({
				'id': seqid2accver(x[1]), 'identity': x[2], 'evalue': x[10],
				'score': x[11], 'coverage': f'{cov:.2f}', 'taxid': ''})
		return res

	"""taxonomy query functions"""

	def update_hit_taxids(self, prots, taxmap={}):
		"""Update hits with taxIds, and update master sequence Id to taxId map.

		Parameters
		----------
		prots : dict of list of dict
			proteins (e.g., search results)
		taxmap : dict, optional
			reference sequence Id to taxId map

		Returns
		-------
		list, list
			sequence Ids still without taxIds
			taxIds added to master map
		"""
		idsotid = set()  # proteins without taxIds
		newtids = set()  # newly added taxIds

		for prot, hits in prots.items():
			for hit in hits:
				id_, tid = hit['id'], hit['taxid']

				# taxId already in hit table
				if tid:
					if id_ not in self.prot2tid:
						self.prot2tid[id_] = tid
						newtids.add(tid)
					continue

				# taxId already in taxon map
				try:
					hit['taxid'] = self.prot2tid[id_]
					continue
				except KeyError:
					pass

				# taxId in reference taxon map:
				try:
					tid = taxmap[id_]
					hit['taxid'] = tid
					self.prot2tid[id_] = tid
					newtids.add(tid)
					continue
				except KeyError:
					pass

				# not found
				idsotid.add(id_)
		return sorted(idsotid), sorted(newtids)

	def remote_taxinfo(self, ids):
		"""Retrieve complete taxonomy information of given taxIds from remote
		server.

		Parameters
		----------
		ids : list of str
			query taxIds

		Returns
		-------
		str
			taxonomy information in XML format

		Raises
		------
		ValueError
			TaxID list is invalid.
		ValueError
			Failed to retrieve info from server.
		"""
		res = self.remote_fetches(ids, 'db=taxonomy&id={}')

		# this error occurs when taxIds are not numeric
		if '<ERROR>ID list is empty' in res:
			raise ValueError('Invalid taxId list.')
		return res

	def parse_taxonomy_xml(self, xml):
		"""Parse taxonomy information in XML format retrieved from NCBI server.

		Parameters
		----------
		xml : str
			taxonomy information in XML format

		Returns
		-------
		list of str
			taxIds added to taxonomy database

		Notes
		-----
		The function will update taxonomy database.
		"""
		added = []

		# get result for each query
		p = re.compile(r'<Taxon>\n'
					   r'\s+<TaxId>(\d+)<\/TaxId>\n.+?'
					   r'\s+<ScientificName>([^<>]+)<\/ScientificName>.+?'
					   r'\s+<ParentTaxId>(\d+)<\/ParentTaxId>.+?'
					   r'\s+<Rank>([^<>]+?)<\/Rank>(.+?)\n'
					   r'<\/Taxon>',
					   re.DOTALL | re.VERBOSE)
		for m in p.finditer(xml):
			tid = m.group(1)
			if tid in self.taxdump:
				continue

			# add query taxId to taxdump
			self.taxdump[tid] = {
				'name': m.group(2), 'parent': m.group(3), 'rank': m.group(4)}
			added.append(tid)

			# get lineage
			m1 = re.search(r'<LineageEx>(.+?)<\/LineageEx>', m.group(5),
						   re.DOTALL)
			if not m1:
				continue

			# move up through lineage
			p1 = re.compile(r'\s+<Taxon>\n'
							r'\s+<TaxId>(\d+)<\/TaxId>\n'
							r'\s+<ScientificName>([^<>]+)<\/ScientificName>\n'
							r'\s+<Rank>([^<>]+)<\/Rank>\n'
							r'\s+<\/Taxon>\n',
							re.DOTALL | re.VERBOSE)

			for m2 in reversed(list(p1.finditer(m1.group(1)))):
				tid_ = m2.group(1)
				pid = self.taxdump[tid]['parent']
				if pid == '':
					self.taxdump[tid]['parent'] = tid_
				elif pid != tid_:
					raise ValueError(
						f'Broken lineage for {tid}: {pid} <=> {tid_}.')
				tid = tid_
				if tid in self.taxdump:
					continue
				self.taxdump[tid] = {
					'name': m2.group(2), 'parent': '', 'rank': m2.group(3)}
				added.append(tid)

			# stop at root
			if self.taxdump[tid]['parent'] == '':
				self.taxdump[tid]['parent'] = '1'
		return added

	"""self-alignment functions"""

	@staticmethod
	def lookup_selfaln(seqs, hits):
		"""Look up self-alignment metrics of sequences from their hit tables.

		Parameters
		----------
		seqs : list of tuple
			query sequences (id, sequence)
		hits : dict of list of dict
			hit tables

		Returns
		-------
		list of tuple
			(id, bitscore, evalue)
		"""
		res = []
		for id_, seq in seqs:
			msg = (f'Cannot find a self-hit for sequence {id_}. Consider '
				   'setting self-alignment method to other than "lookup".')
			if id_ not in hits:
				raise ValueError(msg)
			found = False
			for hit in hits[id_]:
				if hit['id'] == id_:
					res.append((id_, hit['score'], hit['evalue']))
					found = True
					break
			if not found:
				raise ValueError(msg)
		return res

	@staticmethod
	def fast_selfaln(seq):
		"""Calculate self-alignment statistics using built-in algorithm.

		Parameters
		----------
		seq : str
			query sequence

		Returns
		-------
		tuple of (str, str)
			bitscore and evalue

		Notes
		-----
		Statistics are calculated following:

		.. _Official BLAST documentation:
			https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html

		Default BLASTp parameters are assumed (matrix = BLOSUM62, gapopen
		= 11, gapextend = 1), except for that the composition based statistics
		is switched off (comp-based-stats = 0).

		Result should be identical to that by DIAMOND, but will be slightly
		different from that by BLAST.
		"""
		# BLOSUM62 is the default aa substitution matrix for BLAST / DIAMOND
		blosum62 = {'A': 4, 'R': 5, 'N': 6,  'D': 6, 'C': 9,
					'Q': 5, 'E': 5, 'G': 6,  'H': 8, 'I': 4,
					'L': 4, 'K': 5, 'M': 5,  'F': 6, 'P': 7,
					'S': 4, 'T': 5, 'W': 11, 'Y': 7, 'V': 4}

		# calculate raw score (S)
		n, raw = 0, 0
		for c in seq.upper():
			try:
				n += 1
				raw += blosum62[c]

			# in case there are non-basic amino acids
			except KeyError:
				pass

		# BLAST's empirical values when gapopen = 11, gapextend = 1. See:
		# ncbi-blast-2.7.1+-src/c++/src/algo/blast/core/blast_stat.c, line #268
		lambda_, K = 0.267, 0.041

		# calculate bit score (S')
		bit = (lambda_ * raw - log(K)) / log(2)

		# calculate e-value (E)
		e = n ** 2 * 2 ** -bit

		return f'{bit:.1f}', f'{e:.3g}'

	def blast_selfaln(self, seqs):
		"""Run BLAST to align sequences to themselves.

		Parameters
		----------
		seqs : list of tuple
			query sequences (id, sequence)

		Returns
		-------
		list of tuple
			(id, bitscore, evalue)
		"""
		tmpin = join(self.tmpdir, 'tmp.in')
		with open(tmpin, 'w') as f:
			write_fasta(seqs, f)
		cmd = ' '.join((
			self.blastp,
			'-query', tmpin,
			'-subject', tmpin,
			'-num_threads', str(self.threads),
			'-outfmt', '6'))
		extrargs = getattr(self, 'extrargs', None)
		if extrargs:
			cmd += ' ' + extrargs
		ec, out = run_command(cmd)
		if ec:
			raise ValueError(f'blastp failed with error code {ec}.')
		remove(tmpin)
		return(self.parse_self_m8(out))

	def diamond_selfaln(self, seqs):
		"""Run DIAMOND to align sequences to themselves.

		Parameters
		----------
		seqs : list of tuple
			query sequences (id, sequence)

		Returns
		-------
		list of tuple
			(id, bitscore, evalue)
		"""
		# generate temporary query file
		tmpin = join(self.tmpdir, 'tmp.in')
		with open(tmpin, 'w') as f:
			write_fasta(seqs, f)

		# generate temporary database
		tmpdb = join(self.tmpdir, 'tmp.dmnd')
		cmd = ' '.join((
			self.diamond, 'makedb',
			'--in', tmpin,
			'--db', tmpdb,
			'--threads', str(self.threads),
			'--tmpdir', self.tmpdir))
		ec, out = run_command(cmd, merge=False)
		if ec:
			raise ValueError(f'diamond failed with error code {ec}.')

		# perform search
		cmd = ' '.join((
			self.diamond, 'blastp',
			'--query', tmpin,
			'--db', tmpdb,
			'--threads', str(self.threads),
			'--tmpdir', self.tmpdir))
		extrargs = getattr(self, 'extrargs', None)
		if extrargs:
			cmd += ' ' + extrargs
		ec, out = run_command(cmd, merge=False)
		if ec: # RC=0 indicates successful completion
			raise ValueError(f'diamond failed with error code {ec}.')

		remove(tmpin)
		remove(tmpdb)
		return(self.parse_self_m8(out))

	def remote_selfaln(self, seqs):
		"""Perform BLAST search through a remote server.

		Parameters
		----------
		seqs : list of tuple
			query sequences (id, sequence)

		Returns
		-------
		list of tuple
			(id, bitscore, evalue)
		"""
		# further split sequences into halves (to comply with URI length limit)
		batches = self.subset_seqs(seqs, maxchars=int(
			(self.maxchars + 1) / 2)) if self.maxchars else [seqs]
		result = []
		for batch in batches:

			# generate query URL
			query = ''.join([f'>{id_}\n{seq}\n' for id_, seq in batch])
			query = quote(query)
			url = (f'{self.aln_server}?CMD=Put&PROGRAM=blastp&'
				   f'DATABASE={self.db}&QUERY={query}&SUBJECTS={query}')
			if self.extrargs:
				url += '&' + self.extrargs.lstrip('&')
			print(f'Submitting {len(batch)} queries for self-alignment.',
				  end='', flush=True)

			trial = 0
			while True:
				if trial:
					if trial == (self.retries or 0) + 1:
						raise ValueError('Remote self-alignment failed after '
										 f'{trial} trials.')
					print(f'Retry {trial} times.', end='', flush=True)
					sleep(self.delay)
				trial += 1

				# get request Id
				with urlopen(url) as response:
					res = response.read().decode('utf-8')
				m = re.search(r'^    RID = (.*$)', res, re.MULTILINE)
				if not m:
					print('WARNING: Failed to obtain RID.')
					continue
				rid = m.group(1)
				print(f' RID: {rid}.', end='', flush=True)
				sleep(1)

				# check status
				url_ = (f'{self.aln_server}?CMD=Get&FORMAT_OBJECT=SearchInfo&'
						f'RID={rid}')
				starttime = time()
				success = False
				while True:
					with urlopen(url_) as response:
						res = response.read().decode('utf-8')
					m = re.search(r'\s+Status=(.+)', res, re.MULTILINE)
					if not m:
						print('WARNING: Failed to retrieve remote self-'
							  'alignment status.')
						break
					status = m.group(1)
					if status == 'WAITING':
						if time() - starttime > self.timeout:
							print('WARNING: Remote self-alignment timeout.')
							break
						print('.', end='', flush=True)
						sleep(self.delay)
						continue
					elif status in ('FAILED', 'UNKNOWN'):
						print('WARNING: Remote self-alignment failed.')
						break
					elif status == 'READY':
						if 'ThereAreHits=yes' not in res:
							print('WARNING: Remote self-alignment returned no '
								  'result.')
							break
						success = True
						break
					else:
						print('WARNING: Unknown remote self-alignment status: '
							  f'{status}.')
						break
				if not success:
					continue
				sleep(1)

				# retrieve result
				url_ = (f'{self.aln_server}?CMD=Get&ALIGNMENT_VIEW=Tabular&'
						f'FORMAT_TYPE=Text&RID={rid}')
				with urlopen(url_) as response:
					res = response.read().decode('utf-8')
				if '# blastp' not in res or '# Query: ' not in res:
					print('WARNING: Invalid format of remote self-alignment '
						  'results.')
					continue
				print(' Results retrieved.')
				break

			m = re.search(r'<PRE>(.+?)<\/PRE>', res, re.DOTALL)
			out = m.group(1).splitlines()
			result += self.parse_self_m8(out)
		return result

	@staticmethod
	def parse_self_m8(lines):
		"""Extract self-alignment results from m8 format table.

		Parameters
		----------
		lines : list of str
			hit table in BLAST m8 format
			fields: qseqid sseqid pident length mismatch gapopen qstart qend
			sstart send evalue bitscore

		Returns
		-------
		list of tuple
			hit table (id, bitscore, evalue)
		"""
		res = []
		used = set()
		for line in lines:
			x = line.rstrip('\r\n').split('\t')
			if x[0].startswith('#'):
				continue
			if len(x) < 12:
				continue
			if x[0] != x[1]:
				continue
			if x[0] in used:
				continue
			res.append((x[1], x[11], x[10]))
			used.add(x[0])
		return res


# ===========================================================================
# GTDB adaptation (plan P3), appended at module level so that the NCBI
# behavior above stays exactly as it was.
# ===========================================================================

from neoHGT import gtdb as _gtdb
from os.path import basename as _basename

OUTFMT_HINT = ('diamond blastp ... --outfmt 6 qseqid sseqid pident evalue '
			   'bitscore qcovhsp staxids')

arguments += [
	'GTDB mode',
	['--source', 'data source of the reference database; "auto" detects a '
	 'GTDB database from its gtdb/ subdirectory',
	 {'choices': ['auto', 'ncbi', 'gtdb', 'hybrid'], 'default': 'auto'}],
	['--gtdb-stem2taxid', 'explicit gtdb/stem2taxid.tsv.gz used as the '
	 'fallback (route C) when the database was built without embedded taxIds'],
	['--keep-unresolved', 'keep hits whose taxId cannot be resolved, instead '
	 'of dropping them', {'choices': ['yes', 'no']}],
	['--require-staxids', 'fail immediately when the hit table lacks the '
	 'staxids column', {'choices': ['yes', 'no']}],
	['--keep-placeholder-names', 'exempt GTDB placeholder labels such as '
	 'g__DUMJ01 or "MK2026 sp010841923" from --tax-capital/--tax-latin; they '
	 'are official labels, not mis-capitalized or non-Latinate names',
	 {'choices': ['yes', 'no']}],
]


def gtdb_db_dir(self):
	"""Find the GTDB database root this run is working with.

	Returns
	-------
	str or None
		database root that contains a gtdb/ directory, else None
	"""
	if getattr(self, 'source', 'auto') == 'ncbi':
		# the explicit escape hatch stays, but it must be loud: forcing NCBI
		# semantics onto a GTDB database is exactly the silent mismatch the
		# mode detection exists to prevent
		for attr in ('gtdb_stem2taxid', 'taxdump', 'taxmap', 'input',
					 'output', 'db'):
			hint = getattr(self, attr, None)
			if isinstance(hint, str) and hint and _gtdb.find_db_root(hint):
				print(f'WARNING: --source ncbi was given, but {attr} points '
					  f'inside the GTDB database at '
					  f'{_gtdb.find_db_root(hint)}; continuing with NCBI '
					  'semantics as requested. Drop --source (or use '
					  '--source auto) to run the GTDB route.')
				break
		return None
	for attr in ('gtdb_stem2taxid', 'taxdump', 'taxmap', 'input', 'output',
				 'db'):
		hint = getattr(self, attr, None)
		if isinstance(hint, str) and hint:
			root = _gtdb.find_db_root(hint)
			if root:
				return root
	return None


def is_gtdb_mode(self):
	"""Whether this search run targets a GTDB database.

	The value is the one `set_params_gtdb` resolved once from --db; re-probing
	here would redo filesystem checks on every batch and could disagree with the
	decisions already baked into `self.gtdb_labels` and the profile defaults.
	"""
	return bool(getattr(self, 'gtdb_root', None))


def guard_precomp_multisample(self):
	"""Refuse pre-compacted tables from a directory when several are loaded.

	Notes
	-----
	The upstream search loop hands `list(self.pcmap.values())[0]` to every batch
	(search.py:170-173, its own "# todo: general fix"), so with `-m precomp`
	pointing at a directory of per-sample tables, every sample after the first
	would be answered from the wrong file and come back with no hits -- in a
	GTDB run, where taxIds are synthetic, that looks like "no homologs found"
	rather than "wrong input file". Plan appendix F6 asks for this to be
	registered; on the GTDB surface it is refused instead of silently mis-searched.
	"""
	if getattr(self, 'method', None) != 'precomp':
		return
	if len(getattr(self, 'pcmap', None) or {}) < 2:
		return
	raise ValueError(
		f'GTDB search cannot use {len(self.pcmap)} pre-computed tables at '
		'once: neoHGT reads only the first one for every query batch, so the '
		'other samples would silently return zero hits. Run one sample per '
		'invocation (-s <that sample>.tsv), or put all queries in a single '
		'table. See plan appendix F6.')


def guard_precomp_knobs(self):
	"""Refuse the batch and program-argument options under -m precomp.

	Notes
	-----
	`args_wf` asks config.yml for `queries`, `maxchars` and `extrargs` under the
	section named by the search method (search.py:305-306), so with
	`-m precomp` it consults `precomp.queries` and friends -- names the shipped
	config cannot supply, because there is no `precomp:` section. That gap is not
	a missing feature: on this route nothing reads those three values. The query
	set is searched as one batch instead of being subsetted (search.py:160-161),
	and self-alignment is forced to the built-in method (search.py:412-414), so
	the only code that would honour them -- the batched native self-alignment
	that calls `subset_seqs` and passes `extrargs` to BLAST/DIAMOND
	(search.py:930-948) -- is unreachable; its `self.method != 'precomp'` test is
	always true there. Refusing therefore also refuses a `precomp:` section that
	does not exist, which is the honest outcome: a knob that cannot do anything
	must not be settable. Census: 缓存测试/reports/R54_entry_expansion.json.
	"""
	if getattr(self, 'method', None) != 'precomp':
		return
	given = [flag for flag, val in (('--queries', getattr(self, 'queries', None)),
	                               ('--maxchars',
	                                getattr(self, 'maxchars', None)),
	                               ('--extrargs',
	                                getattr(self, 'extrargs', None)))
	         if val not in (None, '', 0, False)]
	if not given:
		return
	raise ValueError(
		f'{", ".join(given)} cannot take effect with -m precomp: a '
		'pre-computed search table is read whole, so --queries / --maxchars '
		'have no batch to subset and --extrargs has no search program to '
		'extend (self-alignment falls back to the built-in method on this '
		'route). Drop ' + ' / '.join(given) + ', or choose -m diamond, '
		'-m blast or -m remote.')


_orig_taxid_wf = Search.taxid_wf
_orig_parse_def_table = Search.parse_def_table
_orig_set_params = Search.args_wf
_orig_update_dmp = Search.update_dmp_files
_orig_taxfilt_wf = Search.taxfilt_wf


FETCH_COUNTS = {'fetch_queries': '--fetch-queries',
		'fetch_retries': '--fetch-retries'}


def coerce_fetch_counts(self):
	"""Make the two routes agree on the fetch batch options.

	These are the only numeric flags in their group that declare no `type=` in
	the `arguments` table, so argparse hands the code a string while the shipped
	config.yml hands it an int. The consumers slice and count with them --
	`cids[i:i + cq]` (search.py:1290) and `range(self.fetch_retries)`
	(search.py:1336) -- so the CLI route dies with
	`TypeError: 'str' object cannot be interpreted as an integer` before it ever
	reaches the network, while the config route runs. Measured in
	缓存测试/reports/R46_fetch_route_types.json.
	"""
	for attr, flag in FETCH_COUNTS.items():
		val = getattr(self, attr, None)
		if val is None or (isinstance(val, int) and not isinstance(val, bool)):
			continue
		try:
			num = int(str(val).strip())
		except ValueError:
			raise ValueError(
				f'Invalid {flag}: {val!r} is not a whole number; it counts '
				f'{"queries per batch" if attr == "fetch_queries" else "retries per search"}.')
		if num < 1:
			raise ValueError(
				f'Invalid {flag}: {num} -- at least 1 is needed or the fetch '
				'loop would never run.')
		setattr(self, attr, num)


def set_params_gtdb(self, args):
	"""Apply GTDB defaults after the regular argument handling."""
	res = _orig_set_params(self, args)
	coerce_fetch_counts(self)
	self.source = getattr(args, 'source', 'auto') or 'auto'
	self.gtdb_stem2taxid = getattr(args, 'gtdb_stem2taxid', None)
	self.keep_unresolved = arg2bool(getattr(args, 'keep_unresolved', None))
	self.gtdb_root = gtdb_db_dir(self)
	# --require-staxids defaults to GTDB modes only: as an unconditional True
	# it rejected the project's own 6-column result tables on plain NCBI runs
	# (an undeclared NCBI behaviour change). The mode must be resolved first.
	if getattr(args, 'require_staxids', None) is None:
		self.require_staxids = is_gtdb_mode(self)
	else:
		self.require_staxids = arg2bool(args.require_staxids)
	if self.source in ('gtdb', 'hybrid') and not self.gtdb_root:
		# analyse.py raises for the same situation; one option must not mean
		# "strict" in one stage and "silently fall back" in the other
		raise ValueError(
			f'--source {self.source}: no GTDB database markers (gtdb/'
			'labels.tsv or gtdb/taxid.map.gz) found for this run. Check the '
			'--db/--taxdump paths; refusing to continue with NCBI semantics '
			'because the hit table and the reference database would '
			'disagree.')
	if self.source == 'auto' and self.gtdb_root:
		print(f'MODE SWITCH: GTDB semantics inferred from {self.gtdb_root} '
		      '(pass --source ncbi to force NCBI behaviour).')
	if self.gtdb_root:
		# GTDB profile defaults (config: taxonomy.profiles.gtdb). Names in a
		# GTDB taxdump are already un-prefixed, but 1,658 of the r220 taxa are
		# numeric placeholder labels ('2-12-FULL-45-22') which a capital or
		# Latin filter would silently discard.
		prof = (((self.cfg or {}).get('taxonomy') or {})
				.get('profiles') or {}).get('gtdb') or {}
		if getattr(args, 'tax_capital', None) is None:
			self.tax_capital = bool(prof.get('capital', False))
		if getattr(args, 'tax_latin', None) is None:
			self.tax_latin = bool(prof.get('latin', False))
		if getattr(args, 'tax_block', None) is None and 'block' in prof:
			self.tax_block = prof['block'] or []
		# placeholder labels are official GTDB names; a style filter would
		# otherwise discard 1,658 r220 taxa' hits in silence
		if getattr(args, 'keep_placeholder_names', None) is None:
			self.keep_placeholder_names = bool(prof.get('keep_placeholder',
														True))
		else:
			self.keep_placeholder_names = arg2bool(args.keep_placeholder_names)
		# GTDB taxIds are synthetic; querying NCBI EFetch for them would
		# return wrong or no taxonomy, so remote fetching is impossible here.
		# 'auto' has already been resolved to 'no' by the time this runs when a
		# taxdump is present, so the user's own request is what must be echoed
		if (getattr(args, 'fetch_enable', None) in ('auto', 'yes')
			or self.fetch_enable == 'yes'):
			print('GTDB mode: forcing --fetch-enable no (NCBI has no record '
				  'of synthetic GTDB taxIds).')
		self.fetch_enable = 'no'
		guard_precomp_multisample(self)
		guard_precomp_knobs(self)
		meta = _gtdb.read_meta(self.gtdb_root)
		self.gtdb_layout = (meta.get('local_layout')
							or _gtdb.LEGACY_RENAMED_LAYOUT)
		self.gtdb_labels = _gtdb.load_labels(self.gtdb_root)
		self.gtdb_stems = None
		# an explicit taxmap is memory-hungry for 333M proteins
		if isinstance(self.taxmap, str) and 'stem2taxid' in self.taxmap:
			pass
		elif self.taxmap:
			print('WARNING: a full protein-to-taxId map can need tens of GB '
				  'for the GTDB representative set; prefer the database '
				  'embedded taxIds or gtdb/stem2taxid.tsv.gz.')
	return res


def _load_route_c(self):
	"""Load the route C tables on first use.

	A missing/empty table used to print "0 genome stems" and continue, which
	routed *every* subject through the unresolved branch and deleted all of
	its hits (RK5 broken in the other direction). An empty table is fatal.
	"""
	if getattr(self, 'gtdb_stems', None) is None:
		root = self.gtdb_root
		stems = _gtdb.load_stem2taxid(root)
		if self.gtdb_stem2taxid and isfile(self.gtdb_stem2taxid):
			stems = _gtdb.load_stem2taxid(
				dirname(dirname(self.gtdb_stem2taxid))) or stems
		self.gtdb_stems = stems
		self.gtdb_g2t = _gtdb.load_genome2taxid(root)
		self.gtdb_p2g = _gtdb.load_prefix2genome(root)
		print(f'  Route C tables: {len(self.gtdb_stems):,} genome stems, '
			  f'{len(self.gtdb_p2g):,} nucleotide prefixes.')
		if not self.gtdb_stems and not self.gtdb_stem2taxid:
			raise ValueError(
				'Route C lookup table is empty: no gtdb/stem2taxid.tsv.gz '
				'next to the database and no --gtdb-stem2taxid given. Every '
				'subject would resolve to "unresolved" and every hit be '
				'deleted; rebuild the database or pass --taxmap instead.')
	return self.gtdb_stems


def open_taxmap_index(self):
	"""Open --taxmap through a disk-backed index (plan P7, route B).

	Parameters
	----------
	self : Search
		instance whose `taxmap` is the path of a sorted protein-to-taxId TSV

	Returns
	-------
	gtdb.TaxidIndex or None
		the cached index, or None when no map was given
	"""
	if not isinstance(getattr(self, 'taxmap', None), str):
		return None
	idx = getattr(self, 'gtdb_taxmap', None)
	if idx is None:
		try:
			idx = _gtdb.TaxidIndex(self.taxmap)
		except ValueError as e:
			raise ValueError(
				f'--taxmap could not be indexed: {e}')
		self.gtdb_taxmap = idx
		print(f'  Route B: indexed {len(idx):,} rows of '
			  f'{_basename(self.taxmap)} on disk (not loaded into memory).')
	return idx


def taxid_wf_gtdb(self, prots):
	"""Resolve hit taxIds for a GTDB database: A, then C, then report.

	Parameters
	----------
	prots : dict of list
		proteins (search results)

	Notes
	----
	Route A is the taxId embedded in the search engine output (DIAMOND
	staxids); route C derives the genome from the sequence Id and looks it up
	in gtdb/stem2taxid.tsv.gz. Unresolvable hits are written to
	unresolved_taxids.txt and unresolved_hits.tsv instead of vanishing
	silently, which was plan risk RK5.
	"""
	if not is_gtdb_mode(self):
		return _orig_taxid_wf(self, prots)

	# route A: taxIds already present in the hit table
	ids2q, added = self.update_hit_taxids(prots)
	n_a = len(ids2q)

	# route C: parse the sequence Id
	if ids2q:
		_load_route_c(self)
		newmap = {}
		for pid in ids2q:
			tid = _gtdb.taxid_of_seqid(
				pid, self.gtdb_layout, self.gtdb_stems, self.gtdb_p2g,
				self.gtdb_g2t)
			if tid:
				newmap[pid] = tid
		if newmap:
			ids2q, added2 = self.update_hit_taxids(prots, newmap)
			print(f'  Route C resolved {len(newmap):,} sequence Ids.')
		else:
			ids2q = sorted(ids2q)

	# route B: the user-supplied --taxmap, accessed through a disk-backed
	# index so that a release-scale map does not have to fit in RAM
	if ids2q:
		idx = open_taxmap_index(self)
		if idx is not None:
			ids2q, added = self.update_hit_taxids(prots, idx)
			if added:
				print(f'  Route B: --taxmap resolved '
					  f'{len(added):,} new sequence Ids.')

	# never contact NCBI in GTDB mode
	if ids2q:
		n = len(ids2q)
		bad = set(ids2q)
		n_hits = sum(len(h) for h in prots.values())
		pct = (100.0 * sum(1 for hits in prots.values()
		                   for hit in hits if hit['id'] in bad) / n_hits
		       if n_hits else 0.0)
		print(f'WARNING: {n} sequence Ids could not be assigned a GTDB taxId '
			  f'({pct:.1f}% of the hits in this batch). Unless '
			  '--keep-unresolved yes was given, their hits are DELETED now '
			  '-- an unresolved donor silently shifts the LCA of every group '
			  'it belongs to. The ids are listed in unresolved_taxids.txt.')
		dn = self.tmpdir if getattr(self, 'tmpdir', None) and isdir(
			self.tmpdir) else self.output
		if dn:
			with open(join(dn, 'unresolved_taxids.txt'), 'w') as f:
				f.write('\n'.join(ids2q) + '\n')
			with open(join(dn, 'unresolved_hits.tsv'), 'w') as f:
				f.write('query\tsubject\tidentity\tevalue\tscore\tcoverage\n')
				for prot, hits in prots.items():
					for hit in hits:
						if hit['id'] in bad:
							f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
								prot, hit['id'], hit['identity'],
								hit['evalue'], hit['score'],
								hit['coverage']))
		if not self.keep_unresolved:
			for hits in prots.values():
				for i in reversed(range(len(hits))):
					if hits[i]['id'] in bad:
						del hits[i]
		self.n_unresolved = n
		return
	print(f'  TaxIds: {len(self.prot2tid):,} subjects resolved '
		  f'(route A covered {len(added)}).')
	return


def parse_def_table_checked(self, lines):
	"""Parse a default tabular hit table, validating the staxids column.

	Parameters
	----------
	lines : list of str
		search result in default tabular format

	Raises
	------
	ValueError
		When the table has fewer than 7 columns, i.e. no staxids column. The
		previous code swallowed the resulting IndexError and printed it, which
		showed up as "0 hits" with no explanation (plan risk RK5).
	"""
	head = None
	for line in lines:
		line = line.rstrip('\r\n')
		if not line or line.startswith('#'):
			continue
		head = line.split('\t')
		break
	if head is not None and len(head) < 7:
		if getattr(self, 'require_staxids', True):
			raise ValueError(
				f'The hit table has {len(head)} columns but 7 are required '
				'(qseqid sseqid pidentity evalue bitscore qcovhsp staxids). '
				'Without the staxids column no hit can be assigned a taxon. '
				f'Re-run the search as: {OUTFMT_HINT}')
		# --require-staxids no: the guard is the user's to waive, but waiving
		# it must be visible -- taxa then come from the subject Id alone
		print(f'WARNING: the hit table has {len(head)} columns and no staxids; '
		      'continuing because --require-staxids no was given. Taxa will be '
		      'resolved from the subject Id only (route C), so any subject '
		      'whose Id does not carry a GTDB stem will be reported without a '
		      'taxon.')
	if head is not None and len(head) == 8 and is_gtdb_mode(self):
		# the router treats >8 columns as m8; exactly 8 lands here, where the
		# def-table reader would take qstart (a coordinate) for the taxId and
		# produce a plausible-looking but wrong HGT table
		raise ValueError(
			'The hit table has 8 columns: this looks like a standard BLAST '
			'outfmt 6 table (qseqid sseqid pident length mismatch gapopen '
			'qstart qend), not the 7-column table this pipeline needs. '
			f'Re-run the search as: {OUTFMT_HINT}')
	if head is not None and len(head) > 8 and is_gtdb_mode(self):
		print(f'NOTE: the hit table has {len(head)} columns; the first 7 are '
			  'used (qseqid sseqid pident evalue bitscore qcovhsp staxids).')
	return _orig_parse_def_table(self, lines)


def update_dmp_guard(self, ids):
	"""Refuse to grow a GTDB taxdump at run time.

	Parameters
	----------
	ids : list of str
		taxIds the caller wants to append

	Notes
	----
	The "add missing taxa while searching" path only makes sense when every
	taxId can be looked up at NCBI. GTDB taxIds are synthetic, so the tree has
	to be complete before the search starts (plan P3 task 3.2).
	"""
	if is_gtdb_mode(self):
		raise ValueError(
			f'Cannot add {len(ids)} taxId(s) to a GTDB taxonomy at run time: '
			'GTDB taxIds are synthetic and have no NCBI records, so the tree '
			'must be complete beforehand. Pass -t <database>/taxdump (built '
			'by "neoHGT database --source gtdb").')
	return _orig_update_dmp(self, ids)


_orig_is_capital, _orig_is_latin = is_capital, is_latin


def taxfilt_wf_gtdb(self, prots):
	"""Filter hits with GTDB's naming reality taken into account.

	Parameters
	----------
	prots : dict of list of dict
		hits per query

	Notes
	-----
	Two departures from the NCBI filter (plan §6.4 tasks 2.6 and 3.5): a rank
	alias that a GTDB dump spells differently is resolved before use, and
	placeholder labels ('g__DUMJ01', 'MK2026 sp010841923' -- 1,658 taxa in r220)
	are exempted from --tax-capital/--tax-latin unless the user asks for the
	strict reading with --keep-placeholder-names no. The upstream filter looks
	both predicates up in this module, so they are rebound for the duration of
	the call and restored afterwards (this process is single-threaded here, and
	the NCBI route never enters this branch).
	"""
	if is_gtdb_mode(self) and self.tax_unirank:
		from neoHGT.util import rank_for_dump
		self.tax_unirank = rank_for_dump(self.tax_unirank, self.taxdump)
	if not is_gtdb_mode(self):
		return _orig_taxfilt_wf(self, prots)
	from neoHGT.util import (is_capital_gtdb, is_latin_gtdb,
	                         is_placeholder_gtdb)
	strict = not getattr(self, 'keep_placeholder_names', False)

	# upstream's --tax-latin resolves the species rank and then tests the taxId
	# it got back *as a name*, so against a real taxdump no hit is ever
	# Latinate; here the taxon name is looked up before judging it
	def _latin(value):
		rec = (getattr(self, 'taxdump', None) or {}).get(str(value))
		name = (rec or {}).get('name') or str(value)
		return is_latin_gtdb(name) or (
			(not strict) and is_placeholder_gtdb(name))

	# the GTDB-aware predicates are the prefix-tolerant versions of the upstream
	# ones, so a name that still carries 's__'/'g__' is judged on its bare name
	# (plan DoD-5), and coded placeholder labels survive unless the strict
	# reading was asked for with --keep-placeholder-names no
	globals()['is_capital'] = (
		lambda n: is_capital_gtdb(n) or (
			(not strict) and is_placeholder_gtdb(n)))
	globals()['is_latin'] = _latin
	try:
		return _orig_taxfilt_wf(self, prots)
	finally:
		globals()['is_capital'] = _orig_is_capital
		globals()['is_latin'] = _orig_is_latin


Search.args_wf = set_params_gtdb
Search.taxid_wf = taxid_wf_gtdb
Search.parse_def_table = parse_def_table_checked
Search.update_dmp_files = update_dmp_guard
Search.taxfilt_wf = taxfilt_wf_gtdb
