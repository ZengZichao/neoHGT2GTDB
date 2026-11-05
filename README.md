English (primary) | [简体中文](README.zh-CN.md)


# neoHGT2GTDB: 🧬 Horizontal Gene Transfer on the GTDB reference frame

> **This is NOT the official neoHGT / neoHGT2.** neoHGT2GTDB is an independent, unofficial adaptation of [neoHGT2](https://github.com/cgneo/neoHGT), maintained in this repository. It is not endorsed by, affiliated with, or supported by the upstream author(s) or [Qiyun lab](https://github.com/qiyunlab). The official project remains at [cgneo/neoHGT](https://github.com/cgneo/neoHGT); no official logo or branding is used here.

**neoHGT2GTDB** is a deep, GTDB-focused adaptation of [neoHGT2](https://github.com/cgneo/neoHGT). It keeps the upstream project intact — the NCBI route, the multi-process database downloader, the search and analysis pipeline — and adds a complete GTDB route on top of it: reference databases built from a local, read-only [GTDB](https://gtdb.ecogenomic.org/) release (verified on R09-RS220), synthetic taxIds recorded in an append-only manifest, GTDB-flavoured ranks and donor labels, an optional hybrid NCBI+GTDB tree, and a machine-checked documentation/test gate for everything above.

The command-line entry point stays `neoHGT`, inherited from upstream, so pre-existing neoHGT/neoHGT2 scripts keep working; the GTDB additions live under the `gtdb` subcommand utilities and the new `--source gtdb|--hybrid` database option (see [GTDB database (experimental)](#gtdb-database-experimental)).

**Naming, in one place:** the *repository* and the *pip distribution* are named **neoHGT2GTDB**; the *import package* and the *command* remain `neoHGT`. That is a deliberate compatibility choice: every documented example, existing neoHGT/neoHGT2 script, and the machine-checked documentation gate are pinned to the upstream spelling, and renaming the executable would break all of them without adding anything. `pip install .` therefore installs a distribution called `neoHGT2GTDB` that provides the `neoHGT` command and the `neoHGT` package.

Project lineage: [HGTector](https://github.com/qiyunlab/HGTector) → neoHGT (v1) → [neoHGT2](https://github.com/cgneo/neoHGT) → **neoHGT2GTDB** (this repository). If you are interested in reproducing the results from the [publication](https://www.nature.com/articles/s41467-023-36487-z), the previous version neoHGT1 may still be found [here](https://github.com/cgneo/neoHGT-v1). Bugs in the NCBI route inherited from neoHGT2 are best reported in [upstream's issue tracker](https://github.com/cgneo/neoHGT/issues); bugs in the GTDB adaptation belong to this repository's issue tracker, and collaborators are welcome.

## Table of Contents
1. [About](#about)
2. [Installation](#installation)
3. [Database](#database)
4. [Compiling](#compiling)
5. [GTDB database (experimental)](#gtdb-database-experimental)
6. [Search](#search)
7. [Analyze](#analyze)
8. [Uninstall](#lastly)
9. [License](#license)

## About
**neoHGT2GTDB** predicts the horizontal/lateral gene transfer (HGT) for species of interest. Like the upstream [neoHGT](https://github.com/cgneo/neoHGT) project it descends from, it is based on [HGTector](https://github.com/qiyunlab/HGTector) developed by Qiyun lab — upstream fixed known HGTector bugs, parallelized the download process with multi-processing, and improved overall stability. What neoHGT2GTDB itself adds is the GTDB route: it rebuilds the reference frame on GTDB's standardized, rank-complete prokaryotic taxonomy (with an optional hybrid mode that grafts the top of the NCBI tree onto the GTDB tree), which measurably improves the statistical behaviour of the self/close/distal grouping, and it documents and machine-checks that behaviour end to end.

  

## Installation
The project is written entirely in Python3.<br>
Recommended installation method:

Step 1/3:

```
conda create -n neoHGT -c conda-forge python=3 pyyaml pandas matplotlib scikit-learn
```

Step 2/3:

Next, you will need to manually install diamond by visiting [the official conda webpage](https://anaconda.org/bioconda/diamond/files?version=0.9.26).

You will see a list of files, please look at the beginning of the filename, which describes the operating system. After downloading the one that corresponds to your system, you can run

```
conda activate neoHGT
```
```
conda install <diamond___.tar.bz2>
```

followed by the file you just downloaded (dragging the file to your terminal) to install diamond.

Step 3/3:

Run the following command:

```
conda activate neoHGT
```
```
pip install .
```

Run this from a clone of this (neoHGT2GTDB) repository. To install the upstream project instead, use `pip install git+https://github.com/cgneo/neoHGT.git`.


## Database

Compared to HGTector, the neoHGT2 downloader that this project inherits not only fixes some known bugs of HGTector but also uses multi-processes to enable parallel downloading, saving you time on the database downloading process.

neoHGT2GTDB keeps the upstream ```server.py```, which downloads and installs the database without a Graphical User Interface.

```server.py``` is **NCBI-only**: it fetches the prebuilt `hgtdb_20230102`
database and knows nothing about GTDB. A GTDB (or `--source hybrid`) database is
built from the command line — see [GTDB database (experimental)](#gtdb-database-experimental)
below.

But the recommended way is to download manually:

## Option1 (Manual)
### (Most comprehensive and most up-to-date but time-consuming):

First, we need ```taxonmap```, where you have two options to obtain the file(s):

This following link will produce a smaller, more compact ```taxdump``` (around 0.4GB after decompression). You can paste the link on your browser to start downloading:
```
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.Z
```

This following link will produce the latest, more comprehensive ```taxdump``` (around 1.8GB after decompression):

```
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.Z
```

___

Next, we will need ```prot.accession2taxid```, and there are two options:

This following link will produce a smaller, more compact ```prot.accession2taxid``` (around 7GB). You can paste the link on your browser to start downloading:

```
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
```

This following link will produce the latest, more comprehensive ```prot.accession2taxid``` (around 12GB):

```
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
```

## Option2 (Manual):
### (Time-saving but less comprehensive and less up-to-date)

There exists a pre-built database on [this OneDrive link](https://arizonastateu-my.sharepoint.com/personal/qzhu44_asurite_asu_edu/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fqzhu44%5Fasurite%5Fasu%5Fedu%2FDocuments%2FPublic%2FHGTector) where you can just download the latest ```hgtdb___.tar.xz``` file, and you can skip the database compiling steps and jump to the [Compiling](#compiling) step.

___

## Database downloading

After obtaining ```taxonmap``` and ```prot.accession2taxid```, we will need to run ```neoHGT database``` command to download around 291,000 genome files to build the diamond database (~90GB) depending on how many categories you have chosen.

The command to build the database is:

```
neoHGT database -c bacteria,archaea,viral -o <FOLDER_NAME>
```

where ```-c``` or equivalently ```-cats``` means category.

This following table explains the command and options available for you to build a customized database.

| Command | Options |
|---------|---------|
| -c (-cats)      | archaea, bacteria, fungi, invertebrate, metagenomes, mitochondrion, plant, plasmid, plastid, protozoa, unknown, vertebrate_mammalian, vertebrate_other, and viral (derived from the [folder names of NCBI RefSeq database](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)) |
| -o      | Output database directory |

___

Known issue:

Don't panic if you see this error message:
```
File "/Users/___/miniconda3/envs/neoHGT/lib/python3.11/site-packages/neoHGT/database.py", line 351, in get_categories
    raise ValueError(
ValueError: "___" is not a valid RefSeq genome category
```

This error occurs either when you typed a category with a wrong name (for example ```virus``` folder does not exist on NCBI server, you should type ```viral``` instead) or when the NCBI server has temporarily shut down the connection from you (refused your connection). There is nothing wrong with the program itself.

The best thing you can do when having this error is to test your connection by entering the following command to the terminal:

```
rsync --list-only --no-motd rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/
```

If you see something like ```rsync: failed to connect to ftp.ncbi.nlm.nih.gov: Connection refused (61)
rsync error: error in socket IO (code 10)```, it confirms that, indeed, there is a network issue (likely from the NCBI server side).

Otherwise, if the previous command returns you a list that starts with ```drwxr-sr-x```, followed by something like ```lrwxrwxrwx```, it means the connection is restored, and you can safely execute the program again.

For more advanced options such as excluding specific taxids, please refer to HGTector's [Database](https://github.com/qiyunlab/HGTector/blob/master/doc/database.md) documentation.

After seeing the program finished with
```
✅ Done
Combined protein sequences written to db.faa.
```
it means that you can proceed with manual compiling.

___

## Compiling

Only DIAMOND is used for compiling, because it is optimized and faster.

### If you followed the instructions of <b>Option1</b>:

```
cd <DATABASE_DIRECTORY>
```

```
mkdir diamond
```

And then carefully replace the entry with ```<‼️>``` with the files you downloaded and extracted. Make sure the files are all under the root folder of your ```<DATABASE_DIRECTORY>```. Below, the command is split over several lines to make things clear, but in practice you should stick to one line without pressing the ```return``` key.

```
diamond makedb --threads 64 --in db.faa
--taxonmap <‼️prot.accession2taxid.FULL>
--taxonnodes <‼️taxdump>/nodes.dmp
--taxonnames <‼️taxdump>/names.dmp
--db diamond/db
```

### If you followed the instructions of <b>Option2</b>:

Step 1/2:

```
cd <DATABASE_DIRECTORY>
```

```
echo $'accession.version\ttaxid' | cat - <(zcat taxon.map.gz) > prot.accession2taxid.FULL
```
In this step, you may encounter a strange error like ```zcat: can't stat: taxon.map.gz (taxon.map.gz.Z): No such file or directory```. If you see this, just simply rename the ```taxon.map.gz``` file into ```taxon.map.gz.Z```.

Step 2/2:

Make sure you are still in the ```<DATABASE_DIRECTORY>```, then enter the following command:

```
mkdir diamond
```
```
diamond makedb --threads 64 --in db.faa --taxonmap prot.accession2taxid.FULL --taxonnodes taxdump/nodes.dmp --taxonnames taxdump/names.dmp --db diamond/db
```

## GTDB database (experimental)

This is the centrepiece of neoHGT2GTDB, and the deep improvement over neoHGT2:
besides the NCBI route inherited from upstream, it can build its reference database
directly from a **local, read-only GTDB release** (verified on R09-RS220). Nothing is downloaded: the
taxonomy is converted from GTDB's own `*_taxonomy_r220.tsv` plus
`*_metadata_r220.tsv`, protein sequences are streamed from the per-genome
`.faa` files, and every synthetic taxId is recorded in an append-only manifest
so that ids stay stable across rebuilds (both `database` and
`neoHGT gtdb taxdump` reuse the manifest of the output directory, or of
`--gtdb-id-manifest`, when one exists).

No workspace setup is required: outputs may live anywhere unless
`$NEOHGT_ROOT` is set, in which case writes must stay inside it. The one hard
rule in every configuration is that the read-only GTDB input dataset is never
written to. The remote sequence channels (`--seq-source ncbi|gtdb-tar|ena`)
have **no download executor in this version**: they write
`gtdb/fetch_plan.tsv` for manual execution and then refuse, instead of
silently building an empty database.

The utilities live under `neoHGT gtdb <subcommand>`
(`probe, taxdump, map-report, gtdbtk, reconcile, tree-dist`), which is the same
interface as `python -m neoHGT.gtdb <subcommand>`; the older flat form
`neoHGT gtdb --action probe|taxdump|map-report` still parses.

```sh
# 0) read-only pre-flight: layout detection + 19 consistency checks + sha256 fingerprints
neoHGT gtdb probe --dir /path/to/GTDB_r220 --out "$NEOHGT_CACHE/logs/layout_probe.json"
# 0b) same checks, but write the accession-level difference lists to disk
neoHGT gtdb reconcile --dir /path/to/GTDB_r220 --out "$NEOHGT_CACHE/reconcile" \
  --report "$NEOHGT_CACHE/reports/reconcile.json"

# 1) development tier: one representative genome per family (564 archaeal genomes, ~2.5 min)
neoHGT database --source gtdb --gtdb-dir /path/to/GTDB_r220 \
  --rank-sample family --gtdb-taxa d__Archaea --seq-source local \
  --compile diamond --threads 16 -o "$NEOHGT_CACHE/db/ar_family"

# 2) search (the DIAMOND db carries the GTDB taxIds, so staxids come out directly)
diamond blastp --query sample.faa --db "$NEOHGT_CACHE/db/ar_family/diamond/db" \
  --threads 16 --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids \
  --max-target-seqs 25 --sensitive -o hits.tsv
neoHGT search -i sample.faa -m precomp -s hits.tsv \
  -t "$NEOHGT_CACHE/db/ar_family/taxdump" -o "$NEOHGT_CACHE/search/sample"

# 3) analyse with GTDB-flavoured ranks and donors
neoHGT analyze -i "$NEOHGT_CACHE/search/sample/sample.tsv" \
  -t "$NEOHGT_CACHE/db/ar_family/taxdump" -o "$NEOHGT_CACHE/analyze/sample" \
  --self-rank genus --donor-rank genus --donor-label prefixed \
  --donor-ncbi-map yes
```

Incompatible CLI shorthand (argparse prefix matching): the new GTDB options
share prefixes with upstream long options, so the abbreviations
`--a --g --ge --gen --l --rep --ty --typ --type` on `database`,
`--r --re --s` on `search` and `--donor-n` on `analyze` became ambiguous and
now fail with "ambiguous option". Spell the long options out; this is the one
breaking change for pre-existing NCBI scripts.

Notes on behaviour that differs from the NCBI route:

- Remote taxonomy fetch is **forced off** (`--fetch-enable no`), because
  synthetic GTDB taxIds have no NCBI record. On the GTDB routes a hit table
  without a `staxids` column is rejected loudly instead of silently dropping
  hits (`--require-staxids` defaults to yes here); on the NCBI route the
  upstream behaviour is preserved. An 8-column table (a standard BLAST
  `outfmt 6`) is refused in GTDB mode: its `qstart` column would be read as
  a taxId.
- `nodes.dmp` is written with the NCBI rank vocabulary (`domain` is stored as
  `superkingdom`) and 13 columns, which `diamond makedb` requires; on the CLI
  both `domain` and `superkingdom` are accepted and mapped automatically.
- Donors are reported as official GTDB labels (`g__Halorubrum`), and
  `--donor-ncbi-map yes` additionally writes `hgts/<sample>.donor_ncbi.tsv`
  mapping each donor back to its NCBI taxId, offline, from the `ncbi_taxid`
  column of the GTDB metadata.
- `neoHGT gtdb reconcile --dir <release> --out <dir>` diffs metadata, GTDB
  taxonomy and the local `.faa` index and writes the five accession-level
  difference lists (`missing_faa`, `extra_faa`, `nonrep_faa`,
  `metadata_not_in_taxonomy`, `taxonomy_not_in_metadata`) even when they are
  empty, so an audit trail exists rather than only a pass/fail count; it exits
  non-zero on any difference. On this r220 copy all five lists are empty
  (596,859 genomes / 113,104 representatives / 113,104 FASTA, 24 s, 1.7 GB).
- Every derived table under `gtdb/` (and `lineages.txt`, the reconcile lists,
  the fetch plan, member expansion) starts with a provenance line
  `# neoHGT version=…; gtdb_release=…; layout=…; generated=…`, so any artifact
  can be traced to the code, release and moment that produced it.
  `taxon.map.gz` and `genome.map.gz` stay comment-free on purpose: DIAMOND,
  BLAST and `read_prot2taxid` consume them verbatim.
- GTDB-Tk output can be used as an input taxonomy without making gtdbtk a
  dependency: `neoHGT gtdb gtdbtk --in BINs.classify.max_summary --db <db>
  --out input_tax.tsv` converts GTDB-Tk lineages (v2 `classification`, v1
  `classify_info`, de novo `gtdbtk.lineage`, or seven rank columns) into the
  synthetic taxIds of that database, matching the deepest assigned rank, and
  the result feeds straight into `neoHGT analyze --input-tax`.
- GTDB renames higher taxa between releases (`p__Proteobacteria` →
  `p__Pseudomonadota`) and the official `synonyms.*.tsv` is not part of a release
  copy, so `neoHGT analyze --taxon-synonyms aliases.tsv` (`alias<TAB>current`)
  maps the old names you kept using onto the current ones for
  `--self-tax/--close-tax/--tax-include/--tax-exclude/--input-tax`, and prints
  which deprecated name was rewritten (plan risk RK11).
- `--keep-placeholder-names` (default `yes`, from
  `taxonomy.profiles.gtdb.keep_placeholder`) exists because 1,658 r220 taxa use
  accession-derived codes as official labels (`2-12-FULL-45-22`, `g__DUMJ01`,
  `s__MK2026 sp010841923`). The GTDB profile turns `--tax-capital`/`--tax-latin` off
  by default for exactly that reason; if you turn one of them on, those codes are
  still kept unless you pass `--keep-placeholder-names no`.
- `-m precomp` pointing at a **directory** of per-sample tables is refused under
  GTDB: the upstream loop feeds every batch the first table (its own
  `# todo: general fix`), so every later sample would silently return zero hits —
  which with synthetic GTDB taxIds looks like "no homologs" rather than "wrong
  input file" (plan appendix F6, registered in `docs/既有缺陷登记.md`).
- `--expand-members yes` grows a representative selection to the member genomes
  of the species clusters it stands for (relation taken from the metadata column
  `gtdb_genome_representative`, since `sp_clusters.tsv` is missing locally) and
  writes the full list to `gtdb/members.tsv.gz`. On this r220 copy the 793
  Enterobacteriaceae representatives expand to 94,393 members, none of which
  ship protein files, so the build stops with an actionable message rather than
  silently shipping a smaller database.
- `neoHGT gtdb tree-dist --gtdb-dir <release> --rank genus --focus g__X`
  reports patristic distances from the release's per-rank trees, which measure
  objectively how many taxa are "close". Archaeal genus trees put only 27 of
  1,846 genera within 0.5 substitutions of the query genus (bacterial: 158 of
  23,111 for Escherichia) - the reason `gtdb_archaea` needs a much smaller
  `closesize` than the bacterial or full-set tiers. Adding `--db <db>
  --close-tax-out close.txt` resolves those labels to the database's synthetic
  taxIds so they can be passed straight to `neoHGT analyze --close-tax`, i.e. a
  phylogeny-defined close group; on the archaeal tier this raised predictions
  from 371 (9.2% of proteins) to 849 (21.1%) while still scoring 0 of the 53
  conserved ar53 marker genes as HGTs.
- `--stream-extraction yes` (NCBI route only; default off, behaviour unchanged)
  replaces the upstream protein extraction with a disk-backed equivalent:
  measured 882 -> 138 MB per million proteins, i.e. ~287 GB -> ~45 GB projected
  for 333M proteins, at 2.3x extraction time and ~100 MiB per million proteins
  of temp disk, with byte-identical `db.faa` / `genome.map.gz` / `taxon.map.gz`.
- Re-running the identical command reuses `db.faa` and the recorded
  `gtdb/extraction_state.json`, so an interrupted multi-hour extraction is not
  repeated; the regenerated database is byte-identical to a single successful
  build except for timestamps.
- `neoHGT search` **appends** to `<output>/<sample>.tsv` (upstream behaviour,
  unchanged here): running it twice into the same `-o` directory writes every
  protein twice and `analyze` will score both copies. Use a fresh `-o`
  directory per run, or move the previous table aside first.
- `db.dmnd` grows with the numeric *range* of taxIds, so the id base defaults
  to `100000000` (`--gtdb-taxid-base`, minimum `10000000`).
- `--source hybrid` grafts the top of the NCBI tree onto the GTDB tree so that
  eukaryotic and viral taxa can share one tree with GTDB prokaryotes: both GTDB
  domains are re-parented under the cellular root (taxId 131567),
  `--hybrid-ncbi-taxdump` points at the NCBI `nodes.dmp`/`names.dmp` (required)
  and `--hybrid-lca-policy root|domain|refuse` decides what happens to proteins
  that both sources claim. Current NCBI dumps spell the top ranks
  `cellular root`, `acellular root`, `realm` and `clade`, and DIAMOND rejects
  every one of them at `makedb`, so `nodes.dmp` collapses any rank outside its
  vocabulary to `no rank` while the merged tree, `gtdb/lineages.tsv` and
  `gtdb/labels.tsv` keep the true rank; what was collapsed is recorded in
  `meta.json` -> `taxonomy.hybrid` and in the "Hybrid graft" section of
  `gtdb/build_report.md`. Verified on real data: the merge itself
  (`tools/verify_hybrid.py`: 25 NCBI nodes grafted onto 146,380 GTDB nodes,
  single root, cross-source LCA at the cellular root) and a 6-genome hybrid
  database whose DIAMOND hits, search table and `scores.tsv` are byte-identical
  to the same 6 genomes built with `--source gtdb`. Searchable sequences are
  still GTDB proteins only, because downloading NCBI *protein* sequences was not
  authorized in this environment.
- The BLAST route is **not usable** with GTDB ids: `makeblastdb` caps local ids
  at 50 characters while GTDB sequence ids are 151-270 characters long. Only
  DIAMOND is supported.

Measured budgets (R09-RS220, 113,104 representative genomes / 333,384,550
proteins in the full set):

| tier | genomes | proteins | db.faa | db.dmnd | wall clock |
| --- | --- | --- | --- | --- | --- |
| 6 genomes (smoke; `缓存测试/db/it6`-class builds) | 6 | 14,373 | 6.2 MB | 614 MB | ≈1 min (81 s; resume 99 s) |
| family sample, Archaea | 564 | 905,829 | 360 MB | 985 MB | 2 min 37 s |
| one family, `--reps` | 793 | 3,247,172 | 1.41 GiB | 2.07 GiB | 4 min 50 s |
| genus sample, Bacteria | 23,112 | ~60-80 M | 28-37 GB | 42-110 GB | extrapolated |
| full representative set | 113,104 | 333,384,550 | ~145 GiB | ~1.5-3x db.faa | ~2-3 h extraction + 6-15 h makedb (extrapolated) |

Measured: `db.faa` takes **466 bytes per protein** (three tiers: 455 / 417 / 466),
which is what the disk pre-flight assumes, with a 3x peak factor for `makedb`.

Every document exists in both languages, English first: this README has a
Chinese mirror at [`README.zh-CN.md`](README.zh-CN.md); the full GTDB guide is
[`docs/GTDB-数据库.md`](docs/GTDB-数据库.md) with its English counterpart
[`docs/en/GTDB-database.md`](docs/en/GTDB-database.md); and each Chinese
register under `docs/` has its English original in
[`docs/en/`](docs/en/).

## Search

The following approach is considered "safer" as you will manually execute the search using ```diamond``` or ```blastp```, which can reduce the chances of bugs due to third-library support (on diamond for example).

Step 1/3:

```
conda activate neoHGT
```

Step 2/3:

(Below, the command is split over several lines to make things clear, but in practice you should stick to one line without pressing the ```return``` key)

```
diamond blastp
--query <input_file.faa>
--db <diamond_db>
--outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids
> <output_directory>/result.txt
```

Step 3/3:

```
neoHGT search -i <input_file.faa>
-m precomp -s <output_directory>/result.txt
-t <taxdump_dir>
-o <output_directory>
```

Checkout the ```.tsv``` file in your ```<output_directory>```, and if it's non-empty, you can proceed to the [Analyze](#analyze) step.

___

If for some reasons ```diamond``` doesn't work, the alternative will be ```blastp```, which can be found [here](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html). But first, you will need to compile blastp database by following this [article](https://www.ncbi.nlm.nih.gov/books/NBK569841/).

The command to make blast database is more or less the following:

```
makeblastdb -parse_seqids -in db.faa
-out <blast_db_path>/<name_of_db>
-title <name_of_db>
-taxid_map <taxon.map>
```

Here the taxon.map file can be obtained by running:

```
gunzip -k taxon.map.gz
```

Then, enter the following command:

```
blastp -query <input_file.faa>
-db <blast_db_path>/<name_of_db>
-outfmt "6 qaccver saccver pident evalue bitscore qcovhsp staxids"
> <output_directory>/result_blast.txt
```

And lastly, 

```
neoHGT search -i <input_file.faa>
-m precomp -s <output_directory>/result_blast.txt
-t <taxdump_dir>
-o <output_directory>
```

Checkout the ```.tsv``` file in your ```<output_directory>```, and if it's non-empty, you can proceed to the [Analyze](#analyze) step.

___

Apart from the manual approach, this following command will perform the search process all in one place, but it is not recommended because a known compatibility bug was discovered and confirmed due to the dependency upgrade of ```diamond```:

```
neoHGT search -i <input_file.faa.gz> -o . -m diamond -p 32 -d <diamond_db> -t <taxdump_dir>
```

To see more details, please refer to [HGTector's issue page](https://github.com/qiyunlab/HGTector/issues/114).

## Analyze

Following upstream neoHGT, the output format of the graph is `.pdf` rather than `.png`, for better resolution and publication purposes.

The command for analyzing the obtained ```.tsv``` file is:

```
neoHGT analyze -i <output.tsv>
-t <taxdump_dir>
-o <output_dir>
```

Which will generate several files, and the most important one is ```scatter.pdf```, which is the main output of the pipeline.

___

For the rigor of research publication purposes, you might be interested in having complete control of the variable such as ```--self-tax``` and ```--close-tax``` because otherwise the program will auto-infer these values.

```
neoHGT analyze -i <output.tsv>
-t <taxdump_dir>
-o <output_dir>
--self-tax <NCBI Taxonomy ID>
--close-tax <NCBI Taxonomy ID>
```
___

Above, the commands were split over several lines to make things clear, but in practice you should stick to one line without pressing the ```return``` key.


## Lastly

To uninstall, you may run:
```
conda remove -n neoHGT --all
```

Note: `bioconda` is a conda *channel*, not a package, so there is nothing to uninstall there; if you also want to remove DIAMOND, either delete the whole environment as above or run `conda remove -n neoHGT diamond`.


## Publication

The upstream project neoHGT (v1) was a result of the author's internship at LOB - Laboratoire d'Optique et Biosciences, where the team published a paper on [Nature Communications](https://www.nature.com/articles/s41467-023-36487-z); neoHGT2GTDB keeps that work as its scientific baseline and citation.

> Filée, J., Becker, H. F., Mellottee, L., Eddine, R. Z., Li, Z., Yin, W., Lambry, J.-C., Liebl, U., & Myllykallio, H. (2023, February 15). Bacterial origins of thymidylate metabolism in Asgard Archaea and Eukarya. Nature Communications. Retrieved February 18, 2023, from https://www.nature.com/articles/s41467-023-36487-z



## Works cited

> Zhu Q, Kosoy M, Dittmar K. HGTector: [an automated method facilitating genome-wide discovery of putative horizontal gene transfers](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-717). *BMC Genomics*. 2014. 15:717.


## License

This repository — the unofficial neoHGT2GTDB adaptation — and its upstream projects are all distributed under the [BSD 3-clause](http://opensource.org/licenses/BSD-3-Clause) license; see the full [statement](LICENSE), which carries every upstream copyright notice.

Copyright (c) 2026-, 曾子超 (Zeng Zichao, [ORCID 0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)) and the neoHGT2GTDB contributors.

Copyright (c) 2023-, Zhihui LI.

Copyright (c) 2013-, [Qiyun Zhu](mailto:qiyunzhu@gmail.com) and [Katharina Dittmar](mailto:katharinad@gmail.com) — authors of upstream HGTector and neoHGT, from which this project derives.
