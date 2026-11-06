English (primary) | [简体中文](README.zh-CN.md)


# neoHGT2GTDB: 🧬 Horizontal Gene Transfer on the GTDB reference frame

> **This is NOT the official neoHGT / neoHGT2.** neoHGT2GTDB is an independent, unofficial adaptation of [neoHGT2](https://github.com/cgneo/neoHGT), maintained in this repository. It is not endorsed by, affiliated with, or supported by the upstream author(s) or [Qiyun lab](https://github.com/qiyunlab). The official project remains at [cgneo/neoHGT](https://github.com/cgneo/neoHGT); no official logo or branding is used here.

**neoHGT2GTDB** is a GTDB-focused adaptation of [neoHGT2](https://github.com/cgneo/neoHGT). It keeps the upstream NCBI route, database downloader, search and analysis pipeline, and adds a GTDB route: reference databases built from a local, read-only [GTDB](https://gtdb.ecogenomic.org/) release, synthetic taxIds recorded in an append-only manifest, GTDB ranks and donor labels, and an optional hybrid NCBI+GTDB tree. The GTDB route has been exercised end to end on a local R09-RS220 copy; the official flat layout has not yet been tested through a complete build (see [GTDB database (experimental)](#gtdb-database-experimental)).

The command-line entry point stays `neoHGT`, inherited from upstream, so pre-existing neoHGT/neoHGT2 scripts keep working. GTDB additions live under the `gtdb` utility subcommands and the new `--source gtdb|hybrid` database option.

**Naming, in one place:** the *repository* and the *pip distribution* are named **neoHGT2GTDB**; the *import package* and the *command* remain `neoHGT`. This is a compatibility choice: documented examples and existing scripts use the upstream spelling. `pip install .` therefore installs a distribution called `neoHGT2GTDB` that provides the `neoHGT` command and the `neoHGT` package.

Project lineage: [HGTector](https://github.com/qiyunlab/HGTector) → neoHGT (v1) → [neoHGT2](https://github.com/cgneo/neoHGT) → **neoHGT2GTDB** (this repository). If you are interested in reproducing the results from the [publication](https://www.nature.com/articles/s41467-023-36487-z), the previous version neoHGT1 may still be found [here](https://github.com/cgneo/neoHGT-v1). Bugs in the NCBI route inherited from neoHGT2 are best reported in [upstream's issue tracker](https://github.com/cgneo/neoHGT/issues); bugs in the GTDB adaptation belong to this repository's issue tracker, and collaborators are welcome.

## Table of Contents

1. [About](#about)
2. [Installation](#installation)
3. [GTDB database (experimental)](#gtdb-database-experimental)
4. [NCBI database](#ncbi-database)
5. [Compiling](#compiling)
6. [Search](#search)
7. [Analyze](#analyze)
8. [Documentation](#documentation)
9. [Lastly](#lastly)
10. [Publication](#publication)
11. [Works cited](#works-cited)
12. [License](#license)

## About

**neoHGT2GTDB** predicts horizontal/lateral gene transfer (HGT) for species of interest. Like the upstream [neoHGT](https://github.com/cgneo/neoHGT) project it descends from, it is based on [HGTector](https://github.com/qiyunlab/HGTector) developed by Qiyun lab. The GTDB route rebuilds the reference frame on GTDB's standardized, rank-complete prokaryotic taxonomy and optionally grafts the top of the NCBI tree onto it. A head-to-head comparison of the statistical behaviour of the GTDB and NCBI routes has not yet been completed; the GTDB guide reports the measurements that are available and the limits of each dataset.

## Installation

The project is written entirely in Python 3.

Step 1/3:

```sh
conda create -n neoHGT -c conda-forge python=3 pyyaml pandas matplotlib scikit-learn
conda activate neoHGT
```

Step 2/3: install DIAMOND from bioconda. The GTDB route has been tested with DIAMOND 2.1.10.

```sh
conda install -c bioconda diamond
diamond version
```

Step 3/3: install this repository.

```sh
conda activate neoHGT
pip install .
```

Run this from a clone of this repository. To install the upstream project instead, use `pip install git+https://github.com/cgneo/neoHGT.git`.

## GTDB database (experimental)

The GTDB route builds a reference database directly from a **local, read-only GTDB release**. Nothing is downloaded: the taxonomy is converted from GTDB's `*_taxonomy_r220.tsv` and `*_metadata_r220.tsv`, protein sequences are read from the per-genome `.faa` files, and every synthetic taxId is recorded in an append-only manifest so IDs stay stable across rebuilds.

### Environment

Outputs may live anywhere unless a workspace is configured. Setting `NEOHGT_ROOT` makes the tool refuse writes outside that root; the read-only GTDB input is never written to either way. The examples below use these variables:

```sh
export GTDB_DIR=/path/to/GTDB_r220
export NEOHGT_ROOT=/path/to/workspace
export NEOHGT_CACHE="$NEOHGT_ROOT/cache"
export NEOHGT_TMP="$NEOHGT_ROOT/tmp"
```

Alternatively, source the bundled wrapper, which derives the code directory from its own location:

```sh
export GTDB_DIR=/path/to/GTDB_r220
. /path/to/neoHGT2GTDB/tools/env.sh
```

The wrapper also honours `tools/env.local` (git-ignored) for machine-specific paths. See [docs/GTDB-数据库.md](docs/GTDB-数据库.md) for the full environment and workspace reference.

### Quick start

```sh
# 0) read-only preflight: layout detection and dataset self-checks.
#    n_checks is dataset-dependent; on the r220 reference copy it is 19.
neoHGT gtdb probe --dir "$GTDB_DIR" \
  --out "$NEOHGT_CACHE/logs/layout_probe.json"

# 1) build: one representative genome per family, archaea (~564 genomes)
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --rank-sample family --domain archaea --seq-source local \
  --compile diamond --threads 16 \
  -o "$NEOHGT_CACHE/db/ar_family"

# 2) search: the DIAMOND database carries GTDB taxIds, so staxids comes out directly
diamond blastp --query sample.faa \
  --db "$NEOHGT_CACHE/db/ar_family/diamond/db" --threads 16 \
  --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids \
  --max-target-seqs 25 --sensitive -o hits.tsv
neoHGT search -i sample.faa -m precomp -s hits.tsv \
  -t "$NEOHGT_CACHE/db/ar_family/taxdump" \
  -o "$NEOHGT_CACHE/search/sample"

# 3) analyze with GTDB ranks and official donor labels
neoHGT analyze -i "$NEOHGT_CACHE/search/sample/sample.tsv" \
  -t "$NEOHGT_CACHE/db/ar_family/taxdump" \
  -o "$NEOHGT_CACHE/analyze/sample" \
  --self-rank genus --donor-rank genus --donor-label prefixed \
  --donor-ncbi-map yes
```

The full guide covers [data preparation](docs/GTDB-数据库.md#数据准备), [probe](docs/GTDB-数据库.md#只读预检-probe), [database](docs/GTDB-数据库.md#建库-database), [search](docs/GTDB-数据库.md#搜索-search), [analyze](docs/GTDB-数据库.md#分析-analyze), [hybrid mode](docs/GTDB-数据库.md#混合源-hybrid), [resource budgets](docs/GTDB-数据库.md#资源预算), and [troubleshooting](docs/GTDB-数据库.md#故障排查). An English version is at [docs/en/GTDB-database.md](docs/en/GTDB-database.md).

### Behaviour that differs from the NCBI route

- Remote taxonomy fetch is **forced off** (`--fetch-enable no`), because synthetic GTDB taxIds have no NCBI record. A hit table without a `staxids` column is rejected (`--require-staxids` defaults to `yes`); an 8-column standard BLAST `outfmt 6` table is refused in GTDB mode because its `qstart` column would be read as a taxId.
- `nodes.dmp` is written with the NCBI rank vocabulary (`domain` is stored as `superkingdom`) and 13 columns, which `diamond makedb` requires. Both `domain` and `superkingdom` are accepted on the CLI.
- Donors are reported as official GTDB labels (`g__Halorubrum`); `--donor-ncbi-map yes` writes `hgts/<sample>.donor_ncbi.tsv`, mapping each donor back to its NCBI taxId offline from the metadata `ncbi_taxid` column.
- `neoHGT gtdb reconcile` diffs metadata, GTDB taxonomy and the local `.faa` index and writes five accession-level difference lists even when they are empty. It exits non-zero on any difference.
- Every derived table under `gtdb/` starts with a provenance line `# neoHGT version=…; gtdb_release=…; layout=…; generated=…`. `taxon.map.gz` and `genome.map.gz` stay comment-free because DIAMOND, BLAST and `read_prot2taxid` consume them verbatim.
- GTDB-Tk output can be converted into the database's synthetic taxIds with `neoHGT gtdb gtdbtk`; the result feeds `neoHGT analyze --input-tax`.
- The analyze stage accepts `--taxon-synonyms aliases.tsv` (`alias<TAB>current`) to map old GTDB names to current ones for `--self-tax`, `--close-tax`, `--tax-include`, `--tax-exclude` and `--input-tax`.
- `--keep-placeholder-names` (default `yes`) keeps official accession-derived labels such as `g__DUMJ01`; 1,658 r220 taxa use such labels.
- `-m precomp` pointing at a directory of per-sample tables is refused under GTDB, because the upstream loop would feed every batch the first table.
- `--expand-members yes` grows a representative selection to the member genomes of its species clusters and writes `gtdb/members.tsv.gz`.
- `neoHGT gtdb tree-dist` reports patristic distances from the release's per-rank trees and can resolve them to synthetic taxIds with `--db` and `--close-tax-out`.
- `--stream-extraction yes` (NCBI route only; default off) replaces the upstream protein extraction with a disk-backed implementation: measured 882 → 138 MB per million proteins, at 2.3× extraction time and about 100 MiB of temp disk per million proteins.
- Re-running the identical build command reuses `db.faa` and `gtdb/extraction_state.json`, so an interrupted extraction is not repeated.
- `neoHGT search` **appends** to `<output>/<sample>.tsv`; use a fresh `-o` directory per run, or move the previous table aside first.
- `db.dmnd` grows with the numeric *range* of taxIds, so the ID base defaults to `100000000` (`--gtdb-taxid-base`, minimum `10000000`).
- `--source hybrid` grafts the top of the NCBI tree onto the GTDB tree. `--hybrid-ncbi-taxdump` points at the NCBI `nodes.dmp`/`names.dmp` (required), and `--hybrid-lca-policy root|domain|refuse` decides what happens to proteins that both sources claim. The merged tree keeps the true ranks in memory and in `gtdb/lineages.tsv`, while `nodes.dmp` collapses ranks outside DIAMOND's vocabulary to `no rank`. Searchable sequences are still GTDB proteins only.
- The BLAST route is **not usable** with GTDB IDs: `makeblastdb` caps local IDs at 50 characters while GTDB sequence IDs are 151–270 characters long. Only DIAMOND is supported.
- `--keep-faa drop` **deletes `db.faa`** after a successful build. Use `--keep-faa compress` if you want to keep the sequences in compressed form.

### Incompatible CLI shorthand

argparse prefix matching means that new long options can make an old abbreviation ambiguous. These abbreviations now fail with "ambiguous option" and must be spelled out:

| Sub-command | Abbreviations that now fail | Previously meant |
| --- | --- | --- |
| `database` | `--a`, `--e`, `--geno`, `--l`, `--ra`/`--ran`, `--rep`, `--s`, `--ta`/`--tax`, `--ty`/`--typ`/`--type` | `--above`/`--latin`/`--genoids`/`--rank`/`--represent`/`--sample`/`--taxids`/`--typemater` and others |
| `search` | `--r`/`--re`, `--s` | `--retries`, `--server` |
| `analyze` | `--t`/`--ta`/`--tax`, `--donor-n` | `--taxdump`, `--donor-name` |

In the `database` sub-command, `--h`, `--k` and `--w` are also newly ambiguous (`--k` is ambiguous in `search` too), but they had no previous long-option meaning. `--g`/`--ge`/`--gen` were already ambiguous upstream between `--genbank` and `--genoids`; they are not caused by this adaptation.

### Measured budgets

R09-RS220, 113,104 representative genomes / 333,384,550 proteins in the full set:

| tier | genomes | proteins | db.faa | db.dmnd | wall clock |
| --- | --- | --- | --- | --- | --- |
| 6 genomes (smoke) | 6 | 14,373 | 6.2 MB | 614 MB | ≈1 min |
| family sample, Archaea | 564 | 905,829 | 360 MB | 985 MB | 2 min 37 s |
| one family, `--reps` | 793 | 3,247,172 | 1.41 GiB | 2.07 GiB | 4 min 50 s |
| genus sample, Bacteria | 23,112 | ~60-80 M | 28-37 GB | 42-110 GB | extrapolated |
| full representative set | 113,104 | 333,384,550 | ~145 GiB | ~1.5-3x db.faa | ~2-3 h extraction + 6-15 h makedb (extrapolated) |

`db.faa` takes **466 bytes per protein** (three tiers: 455 / 417 / 466), which is what the disk pre-flight assumes, with a 3× peak factor for `makedb`.

## NCBI database

Compared to HGTector, the neoHGT2 downloader that this project inherits fixes some known bugs of HGTector and uses multi-processing for parallel downloading.

neoHGT2GTDB keeps the upstream `server.py`, which downloads and installs the database without a Graphical User Interface.

`server.py` is **NCBI-only**: it fetches the prebuilt `hgtdb_20230102` database and knows nothing about GTDB. A GTDB (or `--source hybrid`) database is built from the command line; see [GTDB database (experimental)](#gtdb-database-experimental) above.

The recommended way is to download manually.

### Option1 (Manual)

#### (Most comprehensive and most up-to-date but time-consuming)

First, we need `taxonmap`. You have two options:

The following link produces a smaller, more compact `taxdump` (around 0.4 GB after decompression). Paste it in your browser to start downloading:

```text
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.Z
```

The following link produces the latest, more comprehensive `taxdump` (around 1.8 GB after decompression):

```text
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.Z
```

___

Next, we need `prot.accession2taxid`. There are two options:

The following link produces a smaller, more compact `prot.accession2taxid` (around 7 GB):

```text
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
```

The following link produces the latest, more comprehensive `prot.accession2taxid` (around 12 GB):

```text
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
```

### Option2 (Manual)

#### (Time-saving but less comprehensive and less up-to-date)

There is a pre-built database on [this OneDrive link](https://arizonastateu-my.sharepoint.com/personal/qzhu44_asurite_asu_edu/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fqzhu44%5Fasurite%5Fasu%5Fedu%2FDocuments%2FPublic%2FHGTector). Download the latest `hgtdb___.tar.xz` file; the database download and construction steps can then be skipped, and you can go directly to [Compiling](#compiling).

___

### Database downloading

After obtaining `taxonmap` and `prot.accession2taxid`, run `neoHGT database` to download around 291,000 genome files and build the DIAMOND database (~90 GB), depending on the categories you choose.

```sh
neoHGT database -c bacteria,archaea,viral -o <FOLDER_NAME>
```

`-c` or `--cats` means category. This table explains the options for a custom database:

| Command | Options |
|---------|---------|
| -c (-cats) | archaea, bacteria, fungi, invertebrate, metagenomes, mitochondrion, plant, plasmid, plastid, protozoa, unknown, vertebrate_mammalian, vertebrate_other, and viral (derived from the [folder names of the NCBI RefSeq database](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)) |
| -o | Output database directory |

___

Known issue: do not panic if you see an error like this:

```text
ValueError: "___" is not a valid refseq genome category
```

This occurs when the category name is wrong (for example `virus` does not exist on the NCBI server; use `viral` instead) or when the NCBI server has temporarily refused the connection. The program itself is fine.

Test your connection with:

```sh
rsync --list-only --no-motd rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/
```

If you see `rsync: failed to connect to ftp.ncbi.nlm.nih.gov: Connection refused (61)` and `rsync error: error in socket IO (code 10)`, there is a network problem, likely on the NCBI side. If the command returns a list starting with `drwxr-sr-x` followed by entries such as `lrwxrwxrwx`, the connection is restored and you can rerun the program.

For advanced options such as excluding specific taxIds, refer to HGTector's [Database](https://github.com/qiyunlab/HGTector/blob/master/doc/database.md) documentation.

When the program finishes with:

```text
✅ Done
Combined protein sequences written to db.faa.
```

you can proceed to manual compiling.

___

## Compiling

Only DIAMOND is used for compiling, because it is optimized and faster.

### If you followed Option1

```sh
cd <DATABASE_DIRECTORY>
mkdir diamond
```

Replace the entries marked `<‼️>` with the files you downloaded and extracted. Keep all files under the root of `<DATABASE_DIRECTORY>`. The command is split over lines here for clarity; in practice keep it on one line.

```sh
diamond makedb --threads 64 --in db.faa --taxonmap <‼️prot.accession2taxid.FULL> --taxonnodes <‼️taxdump>/nodes.dmp --taxonnames <‼️taxdump>/names.dmp --db diamond/db
```

### If you followed Option2

Step 1/2:

```sh
cd <DATABASE_DIRECTORY>
```

```sh
echo $'accession.version\ttaxid' | cat - <(zcat taxon.map.gz) > prot.accession2taxid.FULL
```

If you see `zcat: can't stat: taxon.map.gz (taxon.map.gz.Z): No such file or directory`, rename `taxon.map.gz` to `taxon.map.gz.Z`.

Step 2/2:

```sh
mkdir diamond
```

```sh
diamond makedb --threads 64 --in db.faa --taxonmap prot.accession2taxid.FULL --taxonnodes taxdump/nodes.dmp --taxonnames taxdump/names.dmp --db diamond/db
```

## Search

The manual approach is considered safer because you run `diamond` or `blastp` yourself.

Step 1/3:

```sh
conda activate neoHGT
```

Step 2/3:

```sh
diamond blastp --query <input_file.faa> --db <diamond_db> --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids > <output_directory>/result.txt
```

Step 3/3:

```sh
neoHGT search -i <input_file.faa> -m precomp -s <output_directory>/result.txt -t <taxdump_dir> -o <output_directory>
```

Check the `.tsv` file in `<output_directory>`; if it is non-empty, proceed to [Analyze](#analyze).

___

If `diamond` does not work, the alternative is `blastp`, available [here](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html). Compile a BLAST database following [this article](https://www.ncbi.nlm.nih.gov/books/NBK569841/).

```sh
makeblastdb -parse_seqids -in db.faa -out <blast_db_path>/<name_of_db> -title <name_of_db> -taxid_map <taxon.map>
```

The `taxon.map` file can be obtained with:

```sh
gunzip -k taxon.map.gz
```

Then run:

```sh
blastp -query <input_file.faa> -db <blast_db_path>/<name_of_db> -outfmt "6 qaccver saccver pident evalue bitscore qcovhsp staxids" > <output_directory>/result_blast.txt
```

And finally:

```sh
neoHGT search -i <input_file.faa> -m precomp -s <output_directory>/result_blast.txt -t <taxdump_dir> -o <output_directory>
```

Check the `.tsv` file in `<output_directory>`; if it is non-empty, proceed to [Analyze](#analyze).

___

The one-step search is not recommended because of a known compatibility bug caused by a DIAMOND dependency upgrade:

```sh
neoHGT search -i <input_file.faa.gz> -o . -m diamond -p 32 -d <diamond_db> -t <taxdump_dir>
```

For details, see [HGTector's issue page](https://github.com/qiyunlab/HGTector/issues/114).

## Analyze

Following upstream neoHGT, the graph output format is `.pdf` rather than `.png`, for better resolution and publication.

```sh
neoHGT analyze -i <output.tsv> -t <taxdump_dir> -o <output_dir>
```

This generates several files; the most important is `scatter.pdf`, the main output of the pipeline.

___

For publication-grade rigour you may want full control of `--self-tax` and `--close-tax`; otherwise the program infers them.

```sh
neoHGT analyze -i <output.tsv> -t <taxdump_dir> -o <output_dir> --self-tax <NCBI Taxonomy ID> --close-tax <NCBI Taxonomy ID>
```

___

## Documentation

The user documentation is available in both English and Chinese:

| Document | Language | Contents |
| --- | --- | --- |
| [README.md](README.md) | English | This page: installation, GTDB quick start, NCBI route |
| [README.zh-CN.md](README.zh-CN.md) | Chinese | Chinese mirror of this page |
| [docs/GTDB-数据库.md](docs/GTDB-数据库.md) | Chinese | Complete GTDB guide |
| [docs/en/GTDB-database.md](docs/en/GTDB-database.md) | English | Complete GTDB guide |
| [docs/README.md](docs/README.md) | English | Documentation index |
| [docs/cli-reference.md](docs/cli-reference.md) | English | Complete CLI option reference |
| [docs/config-reference.md](docs/config-reference.md) | English | `config.yml` reference |
| [docs/troubleshooting.md](docs/troubleshooting.md) | English | Common problems and fixes |
| [CHANGELOG.md](CHANGELOG.md) | English | Release notes |
| [CONTRIBUTING.md](CONTRIBUTING.md) | English | Development, tests and documentation maintenance |
| [CITATION.cff](CITATION.cff) | — | Machine-readable citation metadata |

The three JSON files under `docs/` ([`新增函数清单.json`](docs/新增函数清单.json), [`方案编号清单.json`](docs/方案编号清单.json), [`路由等价性清单.json`](docs/路由等价性清单.json)) are maintainer-facing machine manifests. They keep their Chinese file names and do not have English counterparts; their purpose, schema and regeneration commands are described in [docs/README.md](docs/README.md).

## Lastly

To uninstall:

```sh
conda remove -n neoHGT --all
```

Note: `bioconda` is a conda *channel*, not a package. To remove DIAMOND separately, run `conda remove -n neoHGT diamond` or delete the whole environment as above.

## Publication

The upstream project neoHGT (v1) was a result of the author's internship at LOB - Laboratoire d'Optique et Biosciences, where the team published a paper in [Nature Communications](https://www.nature.com/articles/s41467-023-36487-z); neoHGT2GTDB keeps that work as its scientific baseline and citation.

> Filée, J., Becker, H. F., Mellottee, L., Eddine, R. Z., Li, Z., Yin, W., Lambry, J.-C., Liebl, U., & Myllykallio, H. (2023, February 15). Bacterial origins of thymidylate metabolism in Asgard Archaea and Eukarya. Nature Communications. Retrieved February 18, 2023, from https://www.nature.com/articles/s41467-023-36487-z

## Works cited

> Zhu Q, Kosoy M, Dittmar K. HGTector: [an automated method facilitating genome-wide discovery of putative horizontal gene transfers](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-717). *BMC Genomics*. 2014. 15:717.

## License

This repository — the unofficial neoHGT2GTDB adaptation — and its upstream projects are distributed under the [BSD 3-clause](https://opensource.org/licenses/BSD-3-Clause) license; see the full [statement](LICENSE), which carries every upstream copyright notice.

Copyright (c) 2026-, 曾子超 (Zeng Zichao, [ORCID 0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)) and the neoHGT2GTDB contributors.

Copyright (c) 2023-, Zhihui LI.

Copyright (c) 2013-, [Qiyun Zhu](mailto:qiyunzhu@gmail.com) and [Katharina Dittmar](mailto:katharinad@gmail.com) — authors of upstream HGTector and neoHGT, from which this project derives.
