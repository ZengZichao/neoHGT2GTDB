"""
////////////////////////////////////////////////////////////
│	
│	Filename:	util.py
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
from os import listdir
from os.path import (
    join, isfile, basename, splitext, dirname, realpath, expanduser)
import subprocess
import gzip
import bz2
import lzma
import datetime

import yaml


zipdict = {'.gz': gzip, '.bz2': bz2, '.xz': lzma, '.lz': lzma}


def timestamp():
    """Get current time.
    """
    return datetime.datetime.now()


def load_configs():
    """Load configurations from file.

    Returns
    -------
    dict or None
        configurations, or None if not found
    """
    file = find_config_file()
    if file is not None:
        with open(file, 'r') as f:
            return yaml.load(f, Loader=yaml.SafeLoader)


def find_config_file():
    """Find configuration file path.

    Returns
    -------
    str or None
        configuration file path, or None if not found
    """
    fname = 'config.yml'
    # current directory
    fp = fname
    # home directory
    if not isfile(fp):
        fp = join(expanduser('~'), '.neoHGT', fname)
    # program directory
    if not isfile(fp):
        fp = join(dirname(realpath(__file__)), fname)
    # final check
    if isfile(fp):
        return fp


def get_config(obj, attr, entry, func=None):
    """Load a single configuration into an attribute.

    Parameters
    ----------
    obj : class
        target class to edit attribute
    attr : str
        target attribute
    entry : str
        configuration entry (period-delimited)
    func : function, optional
        function to manipulate variable
    """
    if obj.cfg is None:
        return
    if getattr(obj, attr, None) is not None:
        return
    keys = entry.split('.')

    # go down hierarchies till last level
    d_ = obj.cfg
    for key in keys[:-1]:
        try:
            d_ = d_[key]
        except KeyError:
            return

    # get value of last level
    try:
        val = d_[keys[-1]]
    except KeyError:
        return

    if val is None:
        return
    if func:
        val = func(val)
    setattr(obj, attr, val)


def arg2bool(val):
    """Convert an argument to boolean.

    Parameters
    ----------
    val : str or bool or None
        original value

    Returns
    -------
    bool
        converted value
    """
    if val is None:
        return False
    if isinstance(val, bool):
        return val
    val = val.lower()
    if val in ('yes', 'true', 't', 'y', '1'):
        return True
    elif val.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ValueError('Boolean value expected.')


def read_file(fp):
    """Read a regular or compressed file.

    Parameters
    ----------
    fp : str
        input filepath

    Returns
    -------
    file handle
        file content
    """
    ext = splitext(fp)[1]
    zipfunc = getattr(zipdict[ext], 'open') if ext in zipdict else open
    return zipfunc(fp, 'rt')


def list_from_param(param):
    """Read list of entries from file or string.

    Parameters
    ----------
    param : str
        parameter which may be a file path or a comma-delimited string

    Returns
    -------
    list
        list of entries
    """
    if not param:
        return []
    elif isinstance(param, list):
        return param
    elif isinstance(param, str):
        if isfile(param):
            with read_file(param) as f:
                return f.read().splitlines()
        else:
            return param.split(',')


def dict_from_param(param):
    """Read dict from file or string.

    Parameters
    ----------
    param : str
        parameter which may be a file path or a string in format of
        key1:value1,key2:value2,...

    Returns
    -------
    dict
        dictionary of keys to values

    Raises
    ------
    ValueError
        If any section in string is not in "key:value" format.
    ValueError
        In any line in file is not in "key<tab>value" format.
    """
    if not param:
        return {}
    elif isinstance(param, dict):
        return param
    elif isinstance(param, str):
        if isfile(param):
            with read_file(param) as f:
                try:
                    return dict(x.split('\t') for x in f.read().splitlines())
                except ValueError:
                    raise ValueError(f'Invalid dictionary file: "{param}".')
        else:
            try:
                return dict(x.split(':') for x in param.split(','))
            except ValueError:
                raise ValueError(f'Invalid dictionary string: "{param}".')


def run_command(cmd, capture=True, merge=True):
    """Run an external command and retrieve screen output.

    Parameters
    ----------
    cmd : str
        command to execute
    capture : bool, optional
        whether capture screen output (stdout) (default: True)
    merge : bool, optional
        whether merge stderr into stdout (default: True)

    Returns
    -------
    int
        exit code
    list or None
        screen output split by line, or None if not capture
    """
    res = subprocess.run(
        cmd,
        # an argv list must not meet a shell: quoting is then the caller's
        # problem, and space-bearing paths stop being command syntax
        shell=not isinstance(cmd, (list, tuple)),
        stdout=(subprocess.PIPE if capture else None),
        stderr=(subprocess.STDOUT if merge else subprocess.DEVNULL))
    return (res.returncode,
            res.stdout.decode('utf-8').splitlines() if capture else None)


def file2id(fname):
    """Extract sample Id from filename.

    Parameters
    ----------
    fname : str
        target filename

    Returns
    -------
    str
        sample Id
    """
    id_, ext = splitext(basename(fname))
    if ext in zipdict:
        id_ = splitext(id_)[0]
    return id_


def id2file_map(dir_, ext=None, ids=None):
    """Generate a map of Ids to files.

    Parameters
    ----------
    dir_ : str
        directory to search for matching files
    ext : str, optional
        filename extension
    ids : iterable, optional
        Id list

    Returns
    -------
    dict
        Id to file map
    """
    res = {}
    if ext is not None:
        if not ext.startswith('.'):
            ext = '.' + ext
        n = len(ext)
    for fname in listdir(dir_):
        id_ = None
        if ext is not None:
            if fname.endswith(ext):
                id_ = fname[:-n]
        else:
            id_, ext_ = splitext(fname)
            if ext_ in zipdict:
                id_ = splitext(id_)[0]
        if id_ is None:
            continue
        if ids is not None and id_ not in ids:
            continue
        if id_ in res:
            raise ValueError(f'Ambiguous files for Id: {id_}.')
        res[id_] = fname
    return res


def read_input_prots(fp):
    """Read input proteins.

    Parameters
    ----------
    fp : str
        input filepath

    Returns
    -------
    list of dict
        information of proteins

    Notes
    -----
    A proteome may contain duplicated protein Ids, in which case the copy
    number will be preserved.
    """
    prots = []
    isfasta = None
    used = {}
    lines = []
    with read_file(fp) as f:
        for line in f:
            line = line.rstrip('\r\n')
            if not line or line.startswith('#'):
                continue

            # determine FASTA vs. plain list by the first character ">"
            if isfasta is None:
                isfasta = line.startswith('>')
            if isfasta:
                lines.append(line)

            # parse plain list
            else:
                # merge duplicated Ids
                if line in used:
                    prot_ = prots[used[line]]
                    prot_['dups'] = prot_.get('dups', 0) + 1
                # add new entry
                else:
                    used[line] = len(prots)
                    prots.append({'id': line, 'product': '', 'seq': ''})

    # parse FASTA
    if isfasta:
        for id_, product, seq in read_fasta(lines):
            if id_ in used:
                prot_ = prots[used[id_]]
                prot_['dups'] = prot_.get('dups', 0) + 1
            else:
                used[id_] = len(prots)
                prots.append({'id': id_, 'product': product, 'seq': seq})

    return prots


def read_fasta(lines):
    """Read sequences from FASTA-format string.

    Parameters
    ----------
    lines : iterable of str
        FASTA-format lines

    Returns
    -------
    list of list
        [id, product, sequence]
    """
    seqs = []
    for line in lines:
        line = line.rstrip('\r\n')
        if line.startswith('>'):
            x = line[1:].split(None, 1)
            seqs.append([x[0], get_product(x[1]) if len(x) > 1 else '', ''])
        else:
            seqs[-1][2] += line.rstrip('*')
    return seqs


def write_fasta(seqs, f):
    """Write sequences to a multi-FASTA file.

    Parameters
    ----------
    seqs : list of iterable
        [id, sequence]
    f : file handle
        file to write
    """
    for id_, seq in seqs:
        f.write(f'>{id_}\n{seq}\n')


def _get_taxon(tid, taxdump):
    """Get information of a given taxId from taxonomy database.

    Parameters
    ----------
    tid : str
        taxId to query
    taxdump : dict
        taxonomy database

    Returns
    -------
    dict
        information of taxon

    Raises
    ------
    ValueError
        If taxId is not found in taxonomy database.
    """
    try:
        return taxdump[tid]
    except KeyError:
        raise ValueError(f'TaxID {tid} is not found in taxonomy database.')


def describe_taxon(tid, taxdump):
    """Generate a string that best describes a taxon.

    Parameters
    ----------
    tid : str
        taxId to describe
    taxdump : dict
        taxonomy database

    Returns
    -------
    str
        description of taxon

    Raises
    ------
    ValueError
        If taxId is not found in taxonomy database.
    """
    taxon = _get_taxon(tid, taxdump)
    name, rank = taxon['name'], taxon['rank']
    return (f'{name} (no rank)' if not rank or rank == 'no rank'
            else f'{rank} {name}')


def is_capital(name):
    """Check if a taxon name is capitalized.

    Parameters
    ----------
    name : str
        taxon name to check

    Returns
    -------
    bool
        whether taxon name is capitalized

    Notes
    -----
    A valid taxon name may start with "[".
    """
    try:
        return name.lstrip('[')[0].isupper()
    except IndexError:
        return False


def is_latin(name):
    """Check if a species name is Latinate.

    Parameters
    ----------
    name : str
        species name to check

    Returns
    -------
    bool
        whether species name is Latinate
    """
    if name == '':
        return False
    elif name.count(' ') != 1:
        return False
    if name[0] == '[':
        i = name.find(']')
        if i == -1:
            return False
        name = name[1:i] + name[i + 1:]
    name = name.replace(' ', '')
    if not name.istitle():
        return False
    elif not name.isalpha():
        return False
    return True


def contain_words(text, words):
    """Check if a string contains any of given words

    Parameters
    ----------
    text : str
        query string
    words : list of str
        target words

    Returns
    -------
    bool
        whether string contains at least one of given words
    """
    return re.search(r'\b{}\b'.format('|'.join(words)), text,
                     re.IGNORECASE) is not None


def get_product(title):
    """Extract protein product from sequence title.

    Parameters
    ----------
    title : str
        original title which reads "product [organism]"

    Returns
    -------
    str
        product

    Notes
    -----
    NCBI protein product xyz [organism]
    """
    title = re.sub(r'^\s+|\s+$', '', title)
    title = re.sub(r'\s+\[.+\]$', '', title)
    return title


def seqid2accver(seqid):
    """Convert NCBI-style sequence Id to accession.version.

    Parameters
    ----------
    seqid : str
        original sequence Id (may be seqid or accver)

    Returns
    -------
    str
        sequence Id in accver format

    Notes
    -----
    Newer NCBI BLAST programs use accver instead of seqid in the standard
    tabular format, producing the desired sequence Id like "NP_123456.1".
    However older programs may produce "ref|NP_123456.1|", which needs to
    be converted to the actual sequence Id.
    """
    m = re.match(r'^[^|]+\|([^|]+)\|$', seqid)
    return m.group(1) if m else seqid


def read_taxdump(dir_):
    """Read NCBI taxdump or compatible taxonomy systems.

    Parameters
    ----------
    dir_ : str
        directory containing nodes.dmp and names.dmp

    Returns
    -------
    dict
        taxonomy database
    """
    taxdump = {}
    with open(join(dir_, 'nodes.dmp'), 'r') as f:
        for line in f:
            x = line.rstrip('\r\n').replace('\t|', '').split('\t')
            taxdump[x[0]] = {'parent': x[1], 'rank': x[2]}
    with open(join(dir_, 'names.dmp'), 'r') as f:
        for line in f:
            x = line.rstrip('\r\n').replace('\t|', '').split('\t')
            if len(x) < 4 or x[3] == 'scientific name':
                try:
                    taxdump[x[0]]['name'] = x[1]
                except KeyError:
                    pass
    return taxdump


def read_prot2taxid(file):
    """Read protein-to-taxId map.

    Parameters
    ----------
    file : str
        protein-to-taxId mapping file

    Returns
    -------
    dict
        protein-to-taxId map

    Notes
    -----
    Two formats are supported:
    - "plain": name <tab> taxId
    - "ncbi": accn <tab> accn.ver <tab> taxId ...
    """
    isncbi = None
    header = ['accession', 'accession.version', 'taxid']
    res = {}
    with read_file(file) as f:
        for line in f:
            x = line.rstrip('\r\n').split('\t')
            if len(x) == 1:
                continue
            if isncbi is None:
                if len(x) >= 3 and x[:3] == header:
                    isncbi = True
                    continue
                else:
                    isncbi = False
            if isncbi:
                res[x[1]] = x[2]
            else:
                res[x[0]] = x[1]
    return res


def get_lineage(qid, taxdump):
    """Get taxIds of self and all ancestral hierarchies in order.

    Parameters
    ----------
    qid : str
        query taxId
    taxdump : dict
        taxonomy database

    Returns
    -------
    list of str
        taxIds in order (low-to-hight)
    """
    cid = qid
    pid = ''
    lineage = [qid]
    while True:
        taxon = _get_taxon(cid, taxdump)
        pid = taxon['parent']
        if pid == cid or pid == '0':
            break
        lineage.append(pid)
        cid = pid
    return lineage


def is_ancestral(qid, tids, taxdump):
    """Loop up the taxdump hierarchies for a particular taxId.

    Parameters
    ----------
    qid : str
        query taxId
    tids : set of str
        target taxIds
    taxdump : dict
        taxonomy database

    Returns
    -------
    bool
        whether any target taxId is reached
    """
    cid = qid
    pid = ''
    while True:
        if cid in tids:
            return True
        pid = _get_taxon(cid, taxdump)['parent']
        if pid == cid or pid == '0':
            break
        cid = pid
    return False


def taxid_at_rank(qid, rank, taxdump):
    """Find taxId at certain rank for a query taxId.

    Parameters
    ----------
    qid : str
        query taxId
    rank : str
        target rank
    taxdump : dict
        taxonomy database

    Returns
    -------
    str or None
        taxId at target rank, or None if not found
    """
    cid = qid
    pid = ''
    while True:
        taxon = _get_taxon(cid, taxdump)
        if taxon['rank'] == rank:
            return cid
        pid = taxon['parent']
        if pid == cid or pid == '0':
            break
        cid = pid
    return None


def taxids_at_ranks(qid, ranks, taxdump):
    """Find taxId at certain rank for a query taxId.

    Parameters
    ----------
    qid : str
        query taxId
    rank : list of str
        target ranks
    taxdump : dict
        taxonomy database

    Returns
    -------
    dict of str or None
        taxIds at target ranks, or None if not found
    """
    cid = qid
    pid = ''
    res = {x: None for x in ranks}
    rankset = set(ranks)
    while True:
        taxon = _get_taxon(cid, taxdump)
        rank = taxon['rank']
        if rank in rankset:
            res[rank] = cid
        pid = taxon['parent']
        if pid == cid or pid == '0':
            break
        cid = pid
    return res


def find_lca(tids, taxdump):
    """Find the lowest common ancestor (LCA) of given taxIds.

    Parameters
    ----------
    tids : iterable of str
        query taxIds
    taxdump : dict
        taxonomy database

    Returns
    -------
    str
        taxId of LCA
    """
    cal = None  # common ancestral lineage (CAL)
    for tid in tids:
        lineage = get_lineage(tid, taxdump)[::-1]
        n = len(lineage)

        # let current lineage be CAL
        if cal is None:
            cal = lineage
            continue

        # attempt to find LCA between current CAL and current lineage
        idx = None
        for i, tid_ in enumerate(cal):
            if i >= n or tid_ != lineage[i]:
                break
            idx = i + 1

        # LCA not found (i.e., even roots are different)
        if idx is None:
            raise ValueError('Cannot find LCA of taxIds in database.')

        # reduce CAL
        if idx < len(cal):
            cal = cal[:idx]

    # LCA is the lowest taxId of CAL
    return cal[-1]


def sort_by_hierarchy(tids, taxdump):
    """Sort a sequence of taxIds by hierarchy from low to high.

    Parameters
    ----------
    tids : list of str
        taxIds to sort
    taxdump : dict
        taxonomy database

    Returns
    -------
    list of str
        sorted taxIds
    """
    # start with any taxId from pool
    seq = [tids[0]]
    pool = [x for x in tids[1:]]

    # loop until pool is drained
    while pool:
        found = False
        for i, tid in enumerate(pool):

            # add to end of sequence
            if taxdump[seq[-1]]['parent'] == tid:
                seq.append(tid)
                found = True

            # add to beginning of sequence
            elif taxdump[tid]['parent'] == seq[0]:
                seq.insert(0, tid)
                found = True

            # remove from pool
            if found:
                pool.pop(i)
                break

        # if none can be added, then they are not sortable (i.e., not a
        # sequence in the taxonomic hierarchy)
        if not found:
            raise ValueError('Cannot sort taxIds by hierarchy.')
    return seq


def refine_taxdump(tids, taxdump):
    """Refine taxonomy database to given taxIds and their ancestors.

    Parameters
    ----------
    tids : iterable of str
        taxIds to include
    taxdump : dict
        taxonomy database

    Returns
    -------
    dict
        refined taxonomy database

    Notes
    -----
    `children` is a list of immediate child taxIds of current taxId.
    """
    ancs = set().union(*[get_lineage(x, taxdump) for x in tids])
    return {k: v for k, v in taxdump.items() if k in ancs}


def add_children(taxdump):
    """Add a `children` property to each taxId in taxonomy database.

    Parameters
    ----------
    taxdump : dict
        taxonomy database

    Notes
    -----
    `children` is a list of immediate child taxIds of current taxId.
    """
    children = {}
    for tid, taxon in taxdump.items():
        pid = taxon['parent']
        if pid != tid and pid != '0':
            children.setdefault(pid, []).append(tid)
    for tid, taxon in taxdump.items():
        taxon['children'] = children.get(tid, [])


def get_descendants(tid, taxdump):
    """Get all descendant taxIds (including self) under a given taxId.

    Parameters
    ----------
    tid : str
        query taxId
    taxdump : dict
        taxonomy database

    Returns
    -------
    list of str
        taxIds of descendants
    """
    res = []
    for cid in _get_taxon(tid, taxdump)['children']:
        res.extend([cid] + get_descendants(cid, taxdump))
    return res


def save_figure(fig, file):
    """Save a figure to a local image file.

    Parameters
    ----------
    fig : matplotlib.pyplot.figure
        figure to save
    file : str
        filename
    """
    fig.tight_layout()
    fig.savefig(file, bbox_inches='tight')


def taxdump_from_text(text):
    """Read taxdump from comma-delimited text.

    Parameters
    ----------
    text : list of str
        multi-line, comma-delimited text
        columns: taxId, name, parent taxId, rank

    Returns
    -------
    dict of dict
        taxonomy database
    """
    res = {}
    for line in text:
        x = line.split(',')
        res[x[0]] = {'name': x[1], 'parent': x[2], 'rank': x[3]}
    return res


def rank_plural(rank):
    """Convert taxonomic rank to plural form.

    Parameters
    ----------
    rank : str
        taxonomic rank

    Returns
    -------
    str
        taxonomic rank in plural form
    """
    if rank.endswith('phylum'):
        return rank[:-2] + 'a'
    elif rank.endswith('class'):
        return rank + 'es'
    elif rank.endswith('family'):
        return rank[:-1] + 'ies'
    elif rank.endswith('genus'):
        return rank[:-2] + 'era'
    elif rank.endswith('species'):
        return rank
    else:
        return rank + 's'


# ---------------------------------------------------------------------------
# Source-aware taxonomy abstractions (GTDB adaptation, additive only).
#
# The functions above keep their original NCBI semantics; everything below is
# new and is dispatched by "source" so that the NCBI code path is untouched.
# ---------------------------------------------------------------------------

NCBI_RANKS = (
    'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family',
    'genus', 'species')
GTDB_RANKS = (
    'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
GTDB_RANK_PREFIXES = {
    'domain': 'd__', 'phylum': 'p__', 'class': 'c__', 'order': 'o__',
    'family': 'f__', 'genus': 'g__', 'species': 's__'}
NCBI_RANK_PREFIXES = {
    'superkingdom': 'k__', 'kingdom': 'k__', 'phylum': 'p__', 'class': 'c__',
    'order': 'o__', 'family': 'f__', 'genus': 'g__', 'species': 's__'}
RANK_ALIASES = {
    'superkingdom': 'domain', 'domain': 'domain', 'regio': 'domain',
    'kingdom': 'domain', 'phylum': 'phylum', 'division': 'phylum',
    'class': 'class', 'order': 'order', 'family': 'family',
    'subspecies': 'species', 'species': 'species', 'species_latin': 'species',
    'no rank': 'no rank'}
RANK_PREFIX_RE = re.compile(r'^([dkpcofgs])_(?:_\s*)?(.*)$')
CLUSTER_SUFFIX_RE = re.compile(r'^(?P<base>.+)_(?P<letter>[A-Z])(?P<rest>\s.*|)$')


def rank_order(source='ncbi'):
    """Get the ordered list of ranks used by a data source.

    Parameters
    ----------
    source : str, optional
        data source name ('ncbi', 'gtdb' or 'hybrid')

    Returns
    -------
    tuple of str
        ranks from high to low

    Raises
    ------
    ValueError
        If the source is unknown.
    """
    if source in ('gtdb', 'hybrid'):
        return GTDB_RANKS
    if source == 'ncbi':
        return NCBI_RANKS
    raise ValueError(f'Unknown taxonomy source: {source}.')


def rank_prefix(rank, source='ncbi'):
    """Get the name prefix of a rank (e.g. 'species' -> 's__').

    Parameters
    ----------
    rank : str
        taxonomic rank
    source : str, optional
        data source name

    Returns
    -------
    str
        rank prefix, or empty string if the source has no prefix convention
    """
    rank = rank_normalize(rank, source)
    if source in ('gtdb', 'hybrid'):
        return GTDB_RANK_PREFIXES.get(rank, '')
    return NCBI_RANK_PREFIXES.get(rank, '')


def rank_normalize(rank, source='ncbi'):
    """Normalize a rank name into the vocabulary of a data source.

    Parameters
    ----------
    rank : str
        taxonomic rank as provided by the user
    source : str, optional
        target data source name

    Returns
    -------
    str
        normalized rank name

    Notes
    ----
    'superkingdom' and 'kingdom' both become 'domain' for GTDB, and 'domain'
    stays 'domain'. Unknown ranks are returned as-is (lowercased) so that the
    NCBI behavior of long-tail ranks ('no rank', 'clade', ...) is preserved.
    """
    if rank is None:
        return rank
    raw = str(rank).strip()
    low = raw.lower()
    alias = RANK_ALIASES.get(low)
    if alias is None:
        return low
    if source in ('gtdb', 'hybrid'):
        return alias
    if alias == 'domain' and source == 'ncbi':
        return 'superkingdom'
    return alias


def rank_plural_label(rank, source='ncbi'):
    """Convert a rank of a given source to its plural display form.

    Parameters
    ----------
    rank : str
        taxonomic rank
    source : str, optional
        data source name

    Returns
    -------
    str
        plural form
    """
    rank = rank_normalize(rank, source)
    if rank == 'domain':
        return 'domains'
    if rank == 'superkingdom':
        return 'superkingdoms'
    return rank_plural(rank)


def strip_rank_prefix(name, strip_cluster=True):
    """Remove a rank prefix and, optionally, a species-cluster suffix.

    Parameters
    ----------
    name : str
        taxon name, possibly like 's__Escherichia coli', 'd__Bacteria',
        or the local file-name form 'g_Nanoarchaeum'
    strip_cluster : bool, optional
        whether to remove a trailing cluster letter ('_A'..'_Z')

    Returns
    -------
    str
        bare taxon name

    Notes
    ----
    GTDB circumscribes species and higher taxa with alphabetic cluster
    suffixes which may occur at any rank ('p__Bacillota_I',
    's__Bacillus_A cereus'). NCBI's is_latin() rejects such names, so they
    must be stripped before name-style checks are applied.
    """
    if name is None:
        return ''
    name = str(name)
    m = RANK_PREFIX_RE.match(name)
    body = name
    local = False
    if m:
        body = m.group(2)
        # 'g_Name' (local file-name style) uses underscores as word
        # separators, while 'g__Name' (official) uses them as the prefix.
        local = not name[len(m.group(1)):].startswith('__')
    if strip_cluster:
        m = CLUSTER_SUFFIX_RE.match(body)
        if m:
            body = m.group('base') + m.group('rest')
    if local:
        body = body.replace('_', ' ')
    return body.strip()


def rank_for_dump(rank, taxdump):
    """Pick the rank spelling that is actually used inside a taxdump.

    Parameters
    ----------
    rank : str
        rank name requested by the user (any alias)
    taxdump : dict
        taxonomy database

    Returns
    -------
    str
        the alias present in the database, preferring the normalized GTDB
        name, then its NCBI counterpart

    Notes
    ----
    A GTDB taxdump stores the domain rank as 'superkingdom', because DIAMOND
    validates the rank column against the NCBI vocabulary (spike S1). Users
    may still say 'domain', so rank lookups resolve through this function.
    """
    want = rank_normalize(rank, 'gtdb')
    present = {t['rank'] for t in taxdump.values()}
    if want in present:
        return want
    for alias in ('superkingdom', 'domain', 'kingdom'):
        if alias in present and RANK_ALIASES.get(alias) == want:
            return alias
    return want


def rank_of_label(label):
    """Infer the GTDB/NCBI rank from a prefixed taxon label.

    Parameters
    ----------
    label : str
        taxon label such as 'f__Enterobacteriaceae' or 'f_Enterobacteriaceae'

    Returns
    -------
    str or None
        normalized rank name, or None if the label carries no prefix
    """
    m = RANK_PREFIX_RE.match(str(label))
    if not m:
        return None
    table = {'d': 'domain', 'p': 'phylum', 'c': 'class', 'o': 'order',
             'f': 'family', 'g': 'genus', 's': 'species', 'k': 'domain'}
    return table.get(m.group(1))


def resolve_taxon_ref(ref, taxdump, labels=None, synonyms=None, source=None):
    """Resolve a user-provided taxon reference into a taxId.

    Parameters
    ----------
    ref : str
        numeric taxId, prefixed label ('g__Escherichia', local form
        'g_Escherichia'), or bare scientific name
    taxdump : dict
        taxonomy database
    labels : dict, optional
        label -> taxId map (GTDB official labels are unique, names may not be)
    synonyms : dict, optional
        alias -> canonical label map
    source : str, optional
        data source name, used for rank alias handling

    Returns
    -------
    str
        resolved taxId

    Raises
    ------
    ValueError
        If the reference is unknown or ambiguous.
    """
    ref = str(ref).strip()
    if not ref:
        raise ValueError('Empty taxon reference.')

    # 1. numeric taxId
    if ref.isdigit():
        if ref in taxdump:
            return ref
        raise ValueError(f'TaxID {ref} is not found in taxonomy database.')

    # 2. user alias table
    if synonyms and ref in synonyms:
        ref = synonyms[ref]

    # 3. exact label match ('g__Escherichia' / 'g_Escherichia')
    if labels and ref in labels:
        return labels[ref]
    if labels and rank_of_label(ref):
        bare = strip_rank_prefix(ref, strip_cluster=False)
        cand = {v for k, v in labels.items()
                if strip_rank_prefix(k, strip_cluster=False) == bare}
        if len(cand) == 1:
            return cand.pop()
        if len(cand) > 1:
            raise ValueError(
                f'Taxon reference "{ref}" is ambiguous: taxIds '
                f'{", ".join(sorted(cand))}. Use a taxId or a full lineage.')

    # 4. bare name match against the taxdump
    want = strip_rank_prefix(ref)
    cand = [tid for tid, taxon in taxdump.items()
            if taxon.get('name', '').lower() == want.lower()]
    if not cand:
        cand = [tid for tid, taxon in taxdump.items()
                if strip_rank_prefix(taxon.get('name', '')).lower()
                == want.lower()]
    if len(cand) == 1:
        return cand[0]
    if len(cand) > 1:
        preview = sorted(cand)[:10]
        raise ValueError(
            f'Taxon name "{ref}" matches {len(cand)} taxIds '
            f'({", ".join(preview)}{"..." if len(cand) > 10 else ""}). '
            'Use a taxId, a prefixed label, or a lineage string.')
    raise ValueError(
        f'Taxon "{ref}" is not found in taxonomy database. If it is a rank '
        'name, use --tax-unirank or a prefixed label (e.g. g__Name).')


def is_latin_gtdb(name):
    """Check whether a GTDB taxon name is Latinate.

    Parameters
    ----------
    name : str
        taxon name with rank prefix and/or cluster suffix

    Returns
    -------
    bool
        whether the bare name is a two-word Latinate binomial
    """
    return is_latin(strip_rank_prefix(name))


def is_capital_gtdb(name):
    """Check whether a GTDB taxon name is capitalized after unprefixing.

    Parameters
    ----------
    name : str
        taxon name, possibly prefixed ('s__Escherichia coli')

    Returns
    -------
    bool
        whether the bare name starts with an uppercase letter
    """
    return is_capital(strip_rank_prefix(name))


def is_placeholder_gtdb(name):
    """Check whether a GTDB taxon label is a coded placeholder name.

    Parameters
    ----------
    name : str
        taxon name, possibly prefixed ('g__DUMJ01')

    Returns
    -------
    bool
        whether the label is an accession-derived code or a bare 'sp012345678'
        strain designation instead of a Latinised name

    Notes
    ----
    These are official GTDB labels: r220 has 1,658 taxa whose name is a code
    ('DUMJ01', '336-3', 'ML615J-28', '2-12-FULL-45-22') plus candidate species
    like 'MK2026 sp010841923'. A --tax-capital/--tax-latin filter would discard
    their hits without saying so, which is what --keep-placeholder-names
    prevents (plan §6.4 task 2.6).
    """
    words = strip_rank_prefix(name or '').replace('_', ' ').split()
    if not words or not any(c.isdigit() for w in words for c in w):
        return False
    if any(re.fullmatch(r'sp\d+', w.lower()) for w in words):
        return True
    return all(not any(c.islower() for c in w) for w in words)


def gtdb_product(title):
    """Extract protein product from a GTDB (Prodigal) sequence title.

    Parameters
    ----------
    title : str
        original title, which carries no product annotation

    Returns
    -------
    str
        always an empty string

    Notes
    ----
    GTDB proteins are annotated by Prodigal, whose defline reads
    "# start # end # strand # ID=...;partial=..." and contains no product
    name. Returning "" avoids the garbage strings that get_product() would
    otherwise produce.
    """
    return ''
