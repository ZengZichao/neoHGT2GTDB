> English | [Simplified Chinese](../GTDB-%E6%95%B0%E6%8D%AE%E5%BA%93.md)
>
> Condensed version: [GTDB database (experimental)](../../README.md#gtdb-database-experimental)

# GTDB database reference

neoHGT2GTDB builds reference databases from a local, read-only [GTDB](https://gtdb.ecogenomic.org/) release and runs homology search and HGT analysis against them. GTDB provides a standardized, hierarchy-complete classification of bacteria and archaea, with exactly one representative genome per species cluster. This tool converts GTDB lineages into an NCBI-style taxdump with stable synthetic taxIds, so DIAMOND, `neoHGT search` and `neoHGT analyze` keep working through their usual interfaces.

GTDB contains no eukaryotes or viruses. To place the top of the NCBI tree (eukaryotes, viruses, and so on) on the same tree as the GTDB prokaryotic tree, use `--source hybrid`. That mode needs a separate NCBI taxdump, and this version can search GTDB protein sequences only (see [Hybrid mode](#hybrid-mode)).

All performance numbers on this page were measured on a local GTDB R09-RS220 copy. Entries marked "extrapolated" are estimates derived from measured throughput and size, not executed runs. Different hardware, filesystems and datasets will give different results.

## Table of contents

1. [Overview](#overview)
2. [Data preparation](#data-preparation)
3. [Environment and workspace](#environment-and-workspace)
4. [Quick start](#quick-start)
5. [Read-only preflight (probe)](#read-only-preflight-probe)
6. [Mismatch lists (reconcile)](#mismatch-lists-reconcile)
7. [Building a taxdump](#building-a-taxdump)
8. [Building a database](#building-a-database)
9. [Search](#search)
10. [Analyze](#analyze)
11. [Other utilities](#other-utilities)
12. [Hybrid mode](#hybrid-mode)
13. [Output files and meta.json](#output-files-and-metajson)
14. [Resuming and release upgrades](#resuming-and-release-upgrades)
15. [Resource budgets](#resource-budgets)
16. [Known limitations](#known-limitations)
17. [Troubleshooting](#troubleshooting)
18. [Appendix: parameter calibration and measurements](#appendix-parameter-calibration-and-measurements)

## Overview

The GTDB route provides:

- Reference-library builds from a local, read-only GTDB release, with no network download.
- Conversion of GTDB lineages into NCBI-style `nodes.dmp`/`names.dmp`, with a stable synthetic taxId for every taxon.
- GTDB ranks (`domain`/`phylum`/`class`/`order`/`family`/`genus`/`species`) and official labels such as `g__Escherichia` in search and analysis.
- Genome selection by representative status, taxon subset, genome list, or sampling.
- Expansion of representatives into species-cluster members (`--expand-members`).
- Conversion of GTDB-Tk classifications into library taxIds for `analyze --input-tax`.
- Patristic distances between taxa from the release's per-rank trees, usable as `analyze --close-tax` input.
- `--source hybrid`, which grafts the top of the NCBI tree onto the GTDB tree.

Building, searching and analysing are separate steps. A database is built once; search and analysis can be repeated for any number of samples.

## Data preparation

### Release contents

A database build needs:

| Item | Description |
| --- | --- |
| Taxonomy files | `*_taxonomy_rNNN.tsv` or `.tsv.gz`, one GTDB lineage per accession |
| Metadata files | `*_metadata_rNNN.tsv`, with representative status, quality classes, NCBI taxId and related columns |
| Protein sequences | One `.faa` file per representative genome |
| Optional: IDlists | `*_Sample_<rank>.IDlist` for fast taxon subsets |
| Optional: trees | `tree_bac120/`, `tree_ar53/` collapsed trees per rank |

Three directory layouts are recognized:

| Layout | Characteristics |
| --- | --- |
| `local-renamed` | A locally curated copy: sequence IDs rewritten, files under `metadata_bac120/`, `metadata_ar53/`, and so on |
| `official-plain` | An unpacked official release with taxonomy files directly in the dataset root |
| `official-gz` | Like `official-plain`, with `.gz` taxonomy/metadata files |

`probe` and `taxdump` detect the layout automatically; `--local-layout` overrides the detection. `gtdb taxdump` needs only taxonomy files, while `neoHGT database --source gtdb` also needs metadata.

### Read-only policy

The dataset at `--gtdb-dir` is read-only throughout. Before writing, the tool rejects any output path inside the input dataset; after a build it also audits input mtimes. Do not place output directories inside the GTDB release.

### Explicit file paths

When automatic discovery does not fit your directory structure, name the files explicitly:

```sh
neoHGT gtdb probe --dir "$GTDB_DIR" \
  --gtdb-taxonomy "$GTDB_DIR/bac120_taxonomy_r220.tsv" \
  --gtdb-metadata "$GTDB_DIR/bac120_metadata_r220.tsv" \
  --out "$NEOHGT_CACHE/logs/layout_probe.json"
```

Separate multiple files with commas. `--gtdb-release` overrides the release label inferred from file names.

## Environment and workspace

### Two modes

**No workspace (default).** Outputs may be written anywhere the user has permission; the tool protects only the read-only GTDB input dataset. This suits users who already manage their output directories.

**Workspace (recommended).** When `NEOHGT_ROOT` is set, every output must resolve inside it. This simplifies cleanup, migration and auditing, and prevents accidental writes to the dataset or system directories.

### Environment variables

| Variable | Purpose |
| --- | --- |
| `GTDB_DIR` | Root of the local GTDB release (read-only) |
| `NEOHGT_ROOT` | Workspace root; when set, every output must stay inside it |
| `NEOHGT_CODE` | Source directory; derived by `tools/env.sh` |
| `NEOHGT_CACHE` | Databases, search results, analysis results and logs |
| `NEOHGT_TMP` | Scratch space for builds and searches; safe to empty |
| `NEOHGT_PYTHON` | Python interpreter |
| `NEOHGT_DIAMOND` | DIAMOND executable |

`tools/env.sh` is a convenience wrapper. It derives `NEOHGT_CODE` from its own location, defaults the workspace root to the parent of the code directory, and creates the cache and scratch directories. It never writes outside the workspace.

```sh
# Option 1: explicit settings (suitable for scripts and cluster jobs)
export GTDB_DIR=/path/to/GTDB_r220
export NEOHGT_ROOT=/path/to/workspace
export NEOHGT_CACHE="$NEOHGT_ROOT/cache"
export NEOHGT_TMP="$NEOHGT_ROOT/tmp"
export DB="$NEOHGT_CACHE/db/ar_family"
export OUT="$NEOHGT_CACHE/analyze/sample"

# Option 2: use the bundled wrapper
export GTDB_DIR=/path/to/GTDB_r220
. /path/to/neoHGT2GTDB/tools/env.sh
```

To avoid repeating the settings in every shell, put `GTDB_DIR` and interpreter paths in the repository's `tools/env.local`; that file is git-ignored and never committed.

The examples below use `$GTDB_DIR`, `$NEOHGT_CACHE` and `$NEOHGT_TMP`. Replace them with your own paths; every command can be copied as is.

## Quick start

The commands below form a minimal GTDB build, search and analysis pipeline.

```sh
# 0) Read-only preflight: detect the layout and run dataset self-checks
neoHGT gtdb probe --dir "$GTDB_DIR" \
  --out "$NEOHGT_CACHE/logs/layout_probe.json"

# 1) Build: one representative per family for archaea, about 564 genomes
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --rank-sample family --domain archaea --seq-source local \
  --compile diamond --threads 16 \
  -o "$NEOHGT_CACHE/db/ar_family"

# 2) Search: the DIAMOND database carries taxIds, so staxids is available
diamond blastp --query sample.faa \
  --db "$NEOHGT_CACHE/db/ar_family/diamond/db" --threads 16 \
  --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids \
  --max-target-seqs 25 --sensitive -o hits.tsv
neoHGT search -i sample.faa -m precomp -s hits.tsv \
  -t "$NEOHGT_CACHE/db/ar_family/taxdump" \
  -o "$NEOHGT_CACHE/search/sample"

# 3) Analyze with GTDB ranks and official donor labels
neoHGT analyze -i "$NEOHGT_CACHE/search/sample/sample.tsv" \
  -t "$NEOHGT_CACHE/db/ar_family/taxdump" \
  -o "$NEOHGT_CACHE/analyze/sample" \
  --self-rank genus --donor-rank genus --donor-label prefixed \
  --donor-ncbi-map yes
```

Inputs, outputs and options for each step are described below. A reference library should span multiple families or higher ranks; a single-family library lacks distant taxa and usually cannot form a distal group.

## Read-only preflight (probe)

```sh
neoHGT gtdb probe --dir "$GTDB_DIR" \
  --out "$NEOHGT_CACHE/logs/layout_probe.json"
```

`probe` reads the dataset without writing to it and produces a JSON report containing:

- the layout, release, sequence-ID convention and `_SerialNumber_` anchor coverage;
- genome-file counts, domain sets and taxonomy/metadata file lists;
- available IDlist tiers and whether representative genomes have local sequences;
- the self-check result `selfcheck`, whose `n_checks` is the number of checks actually run on this dataset.

The number of self-checks is dataset-dependent: each metadata table, taxonomy file and IDlist tier contributes checks. The r220 reference copy reports 19, but other datasets may differ. Use `selfcheck.n_checks` and `selfcheck.status` from the output; do not treat 19 as a fixed contract.

All checks run by default. Skipping them weakens the reproducibility guarantee:

| Sub-command | Option that skips self-checks |
| --- | --- |
| `neoHGT gtdb probe` | `--no-selfcheck` |
| `neoHGT gtdb taxdump` | `--no-selfcheck` |
| `neoHGT database --source gtdb` | `--gtdb-no-selfcheck` |

`probe` also accepts `--no-deep`, which skips the per-accession alignment checks for a fast look at large datasets. Invalid input exits non-zero; `probe` never writes a database.

## Mismatch lists (reconcile)

`probe` reports whether counts agree. To learn *which* accessions disagree, use `reconcile`:

```sh
neoHGT gtdb reconcile --dir "$GTDB_DIR" \
  --out "$NEOHGT_CACHE/reconcile" \
  --report "$NEOHGT_CACHE/reports/reconcile.json"
```

It computes set differences among the full metadata table, the GTDB taxonomy and the local `.faa` index, and writes five TSVs:

| File | Meaning |
| --- | --- |
| `missing_faa.tsv` | Accessions marked as representatives in metadata but with no local sequence file |
| `extra_faa.tsv` | Accessions with a local `.faa` that metadata does not recognize |
| `nonrep_faa.tsv` | Local files of non-representative (member) genomes |
| `metadata_not_in_taxonomy.tsv` | In metadata but not in taxonomy |
| `taxonomy_not_in_metadata.tsv` | In taxonomy but not in metadata |

Each table has a comment header with source, release and generation time, so it leaves a trace even when empty. Any difference makes the command exit non-zero. On the r220 reference copy all five tables are empty; run this once for every new dataset before building.

## Building a taxdump

```sh
neoHGT gtdb taxdump --dir "$GTDB_DIR" \
  --out "$NEOHGT_CACHE/db/r220_taxdump"
```

The output is `taxdump/nodes.dmp`, `taxdump/names.dmp`, and the mapping tables under `gtdb/`. On the r220 reference copy:

- 146,380 nodes (GTDB 146,379 + root);
- rank counts: domain 2, phylum 217, class 607, order 2,030, family 5,460, genus 24,959, species 113,104;
- about 42.9 s, peak memory about 1.5 GB, artifacts about 192 MB.

### Datasets with lineage files only

A taxdump is a function of the lineage strings and does not need metadata. An official release may place taxonomy files directly in the dataset root; `gtdb taxdump` recognizes that layout and enters `lineage_only` mode:

```sh
neoHGT gtdb taxdump --dir /path/to/GTDB_r214/parent_directory \
  --gtdb-id-manifest /path/to/old_db/gtdb/taxid.map.gz \
  --out "$NEOHGT_CACHE/db/r214_taxdump"
```

This mode runs only the checks that do not need metadata (for example, lineage depth always 7; one accession must not have two lineages), and prints the skipped checks. The r214 reference copy measured about 113,318 nodes in 7.9 s. A full database build still requires metadata.

Malformed lineages are rejected before any write, with the file and line reported. For example, a lineage depth other than 7 or two lineages for one accession exits with status 1 and does not create the output directory.

### Reusing an old ID manifest

When upgrading a release, point `--gtdb-id-manifest` at the old library's `gtdb/taxid.map.gz` to keep existing taxIds unchanged. Taxa that disappear from the new release are written to `gtdb/delnodes.txt` and their taxIds are kept in `gtdb/reserved_ids.txt`, so later taxa cannot reuse them. `--gtdb-taxid-base` must match the old library (default `100000000`).

## Building a database

### Basic command

```sh
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --reps --seq-source local --compile diamond --threads 16 \
  -o "$NEOHGT_CACHE/db/r220_reps"
```

At the `database` stage `--source` accepts `ncbi`, `gtdb` and `hybrid`; there is no `auto`. When omitted, the value comes from `database.source` in `config.yml` (shipped value `ncbi`), and the command line always wins.

### Common options

| Option | Purpose |
| --- | --- |
| `--gtdb-dir` | GTDB release root; falls back to `$GTDB_DIR` or `config.yml` |
| `--gtdb-taxa` | Keep only descendants of these GTDB taxa; accepts `f__Enterobacteriaceae` or GTDB taxIds |
| `--reps` | Use species-cluster representatives only |
| `--genome-list` | Read accessions to include from a file, one per line |
| `--rank-sample <rank>` | Take the highest-quality genome per taxon |
| `--sample N --rank <rank>` | Take at most N genomes per `<rank>` taxon; can be combined with `--above` |
| `--max-genomes N` | Truncate the selection to the first N genomes for smoke tests |
| `--expand-members yes` | Expand representatives into species-cluster members |
| `--min-completeness` / `--max-contamination` | Filter by genome completeness and contamination |
| `--mimag high\|medium\|low\|any` | Filter by MIMAG quality class |
| `--type-material gtdb\|ncbi` | Filter by type-material definition; `--typemater` is an alias for `ncbi` |
| `--domain bacteria\|archaea\|both` | Build one domain or both (default `both`) |
| `--seq-source local\|ncbi\|gtdb-tar\|ena` | Sequence source; only `local` has a reader in this version |
| `--compile none\|diamond\|blast\|both` | Search database to compile; the GTDB route supports DIAMOND only |
| `--threads N` | DIAMOND compilation threads |
| `--tmpdir DIR` | Scratch directory; defaults to `workspace.tmp_dir` or the system temporary directory |
| `--keep-faa keep\|compress\|drop` | What to do with `db.faa` after the build; default `keep` |
| `--gtdb-taxid-base N` | First synthetic taxId; default `100000000`, minimum `10000000` |
| `--gtdb-id-manifest PATH` | Reuse an old library's `gtdb/taxid.map.gz` |
| `--input-fingerprints yes` | Record sha256 and size of input files |
| `--workspace-out DIR` / `--workspace-tmp DIR` | Override the output/scratch directory for this run |
| `--overwrite` | Overwrite existing artifacts and re-extract |
| `--force` | Skip the capacity precheck or force an unsupported combination |

`--keep-faa drop` **deletes** `db.faa` after a successful build. Use `keep` or `compress` if the file is still needed. `compress` writes `db.faa.gz` and removes the original.

Defaults: `--seq-source` is `local`; `--allow-network` defaults to `no` (network access stays off); `--gtdb-taxid-base` defaults to `100000000`.

### Sampling and `--above`

`--sample N --rank <rank>` first takes at most N genomes per `<rank>` taxon. Genomes whose lineage does not reach that rank do not participate in this round; with `--above`, the program starts one level above the sampling rank and tops up each coarser taxon to N (up to phylum). `--rank-sample <rank>` is a lower-cost alternative: it takes the highest-quality genome per taxon and ignores `--above`.

On the GTDB representative set, lineages normally extend to species, so `--above` often leaves the result unchanged. It matters mainly for shallower inputs (GTDB-Tk results, member-genome expansion).

### Quality-class options

The r220 metadata table has 596,859 rows (bacteria 584,382 / archaea 12,477), of which 113,104 are representatives. The counts below come from the production selection code over the full table:

| Option | Rule | Full table | Representatives |
| --- | --- | ---: | ---: |
| `--type-material gtdb` | `gtdb_type_species_of_genus == 't'` | 4,545 | 3,347 |
| `--type-material ncbi` | Non-empty `ncbi_type_material_designation` other than `na` | 22,256 | 16,371 |
| `--mimag high` | `mimag_high_quality == 't'` | — | 27,292 |
| `--mimag medium` | `mimag_medium_quality == 't'` | — | 85,752 |
| `--mimag low` | `mimag_low_quality == 't'` | — | 8 |
| `--mimag any` | At least one MIMAG tier is `t` | 596,497 | 113,052 |

`--complete` is rejected on the GTDB route: upstream checks `ncbi_assembly_level`, which describes how a sequence was submitted, not how complete the genome is. Use `--min-completeness` and `--max-contamination` instead. On the r220 representative set, completeness ≥ 50 keeps 82,662 genomes, contamination ≤ 10 keeps 82,464, and both together keep 82,161.

`--mimag` errors out when the required columns are absent instead of silently returning the full table.

### Expanding to member genomes

GTDB publishes one representative genome per species cluster. `--expand-members yes` uses the metadata column `gtdb_genome_representative` to expand representatives into members:

```sh
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --gtdb-taxa f__Enterobacteriaceae --reps --expand-members yes \
  --compile none -o "$NEOHGT_CACHE/db/entero_members"
```

On the r220 reference copy, 793 Enterobacteriaceae representatives expand to 94,393 member genomes, written to `gtdb/members.tsv.gz`. If member genomes have no local protein files, the program stops before extraction and points to the official `*_genomic_files_all` tar or a remote sequence channel.

### Compilation and capacity

`--compile diamond` builds `diamond/db.dmnd` with DIAMOND. GTDB sequence IDs are 151–270 characters long, above the 50-character limit of `makeblastdb`, so `--compile blast` or `both` is rejected by default; `--force` attempts it, but the BLAST route is not usable.

A disk-capacity precheck runs before the build: free space must be at least 3× the estimated `db.faa` size, otherwise the build stops. `--force` skips the precheck. `--compile none` writes `db.faa`, the taxdump and side tables without a DIAMOND database, and records `diamond_version` as `null` in `meta.json`.

## Search

The GTDB DIAMOND database already carries taxIds, so the `staxids` column of a hit table is directly usable. The recommended path is manual search plus precomputed results:

```sh
diamond blastp --query q.faa --db "$DB/diamond/db" --threads 8 \
  --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids \
  --max-target-seqs 25 --sensitive > hits.tsv

neoHGT search -i q.faa -m precomp -s hits.tsv -t "$DB/taxdump" \
  -o "$NEOHGT_CACHE/search/sample" --tmpdir "$NEOHGT_TMP/search/$$"
```

Behaviour in GTDB mode:

- `--fetch-enable no` is forced, because synthetic taxIds have no NCBI records.
- `--require-staxids` defaults to `yes`; a hit table without a `staxids` column is an error. An explicit `no` downgrades it to a warning.
- An 8-column standard BLAST `outfmt 6` table is rejected: it has no `staxids`, and its `qstart` column would be read as a taxId.
- Hits whose taxId cannot be determined are written to `unresolved_taxids.txt` and `unresolved_hits.tsv`.
- `search` appends to `<output>/<sample>.tsv`. Reusing the same `-o` directory writes every protein twice; use a fresh directory per run, or move the old table aside first.

`--keep-placeholder-names` defaults to `yes`, from `taxonomy.profiles.gtdb.keep_placeholder` in `config.yml`. GTDB has 1,658 taxa whose official labels are accession-derived codes (for example `2-12-FULL-45-22`, `g__DUMJ01`, `s__MK2026 sp010841923`). When `--tax-capital` or `--tax-latin` is enabled, these labels are exempted; `no` removes the exemption.

`--taxon-synonyms` accepts a two-column `alias<TAB>current` TSV mapping old GTDB names (such as `p__Proteobacteria`) to current ones (`p__Pseudomonadota`). The alias table is applied before resolving `--self-tax`, `--close-tax`, `--tax-include`, `--tax-exclude` and `--input-tax`.

`--taxmap` is the fallback when a hit table has no `staxids` and sequence IDs cannot be resolved otherwise. It needs a plain-text, byte-sorted mapping file (usually an unpacked `taxon.map.gz`) and uses a disk-anchored index instead of loading the whole mapping into memory. Passing `.gz` is rejected.

## Analyze

```sh
neoHGT analyze -i "$NEOHGT_CACHE/search/sample/sample.tsv" \
  -t "$DB/taxdump" -o "$NEOHGT_CACHE/analyze/sample" \
  --self-rank domain --donor-rank genus --donor-label prefixed \
  --donor-ncbi-map yes
```

### Ranks and labels

- `--self-rank`, `--donor-rank` and `--tax-unirank` accept `domain` and map it to the `superkingdom` stored in the taxdump; `superkingdom` is also accepted.
- `--donor-label prefixed` is the default and outputs official GTDB labels such as `g__Escherichia`.
- `--donor-ncbi-map yes` writes `hgts/<sample>.donor_ncbi.tsv`, mapping donors back to NCBI taxIds through the metadata `ncbi_taxid` column, entirely offline.

### Grouping parameters

At the `search` and `analyze` stages `--source` accepts `auto`, `ncbi`, `gtdb` and `hybrid`, defaulting to `auto`. `auto` detects the source from the library directory; `gtdb` errors out when the GTDB markers are missing; `hybrid` prints a warning and runs the upstream NCBI path when they are missing.

GTDB grouping uses a default `closesize` of 200 (50 for the `gtdb_archaea` profile), because the GTDB tree is denser than the NCBI tree. The domain is determined by a majority vote of the domain-level ancestors of each protein's best hit (a tie counts as undetermined). The log lines `GTDB profile: input domain = ...` and `Close group has N taxa.` report the values actually used.

`--close-size` defaults to 200 for bacteria and 50 for archaea. It is a lower bound, not a quota. When the tree cannot supply enough taxa the group is truncated; cross-sample comparisons must use the number in the log. `--bandwidth` accepts `auto`, `grid`, `silverman`, or a number between 0.1 and 1.0. The direction of `--noise` is not guaranteed to be monotonic; tune it against your own sample's curve.

`--from-scores` reuses an existing `scores.tsv`, skips the grouping stage, and notes in the log which grouping parameters do not affect this run. Both paths should agree on the protein set and ranking; floating-point columns may differ in the last digit because `scores.tsv` stores six significant digits.

### Measured close_size

On a library of 564 archaeal representative genomes with 4,025 query proteins:

| close_size | 10 | 50 | 100 | 200 | 400 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Predicted HGTs | 371 | 371 | 132 | 0 | 0 |

At 200 and 400 the close group swallows every best hit, the distal group is empty, and the log prints `WARNING: No hit is assigned to distal group.` The archaeal profile therefore defaults to 50; larger values need a library close to the full representative set to be meaningful.

Bacterial samples do not share one curve. On a 4,896-genome "one representative per family" library, *Escherichia* stays at 771 for 10/50/100/200, Alicyclobacillaceae at 753, and f__ATN3 at 714; Arcobacteraceae gives 436 at 10/50 and drops to zero at 100/200 because its self group has only two taxa, and a larger close_size jumps to the whole phylum and swallows every hit. "The bacterial tier 200 works" is therefore a per-sample property, not a property of the library.

### donor=0

`donor=0` means the LCA of the distal hits is coarser than the rank requested by `--donor-rank`; it does not mean "no donor". In an E. coli query against a 4,896-genome cross-family library, genus-level requests annotate 10.4% of predictions, family-level 13.1%, and domain-level 100%. Among the `0` LCA nodes, 591 are at superkingdom, 65 at class, 21 at phylum and 14 at order. Loosen `--donor-rank` for higher resolution.

### Explicit close groups and input taxonomy

`--close-tax` accepts a taxId list produced by `tree-dist --close-tax-out`, giving a phylogenetically defined close group. `analyze` first prunes the taxonomy tree to the range covered by the hits; taxIds that are completely unreachable are dropped with a warning. If none are reachable, the command errors out and explains the usage.

`--input-tax` accepts a `genome_stem<TAB>taxId` table produced by `gtdb gtdbtk` to assign an input taxonomy to the query genome.

## Other utilities

### GTDB-Tk conversion

GTDB-Tk is not a runtime dependency; `gtdb gtdbtk` only reads its result table:

```sh
neoHGT gtdb gtdbtk --in BINs.classify.max_summary \
  --db "$NEOHGT_CACHE/db/ar_family" \
  --out "$NEOHGT_CACHE/gtdbtk/input_tax.tsv" \
  --report "$NEOHGT_CACHE/gtdbtk/gtdbtk_report.json"
```

It recognizes v2 `classification`, v1 `classify_info`, de novo `gtdbtk.lineage`/`lineage`, and seven-rank-column tables. The taxId is the deepest matchable prefix; lineages absent from the library are reported as unmapped and no ID is invented. The report JSON includes the `matched_rank` distribution and the reason for each unmatched genome.

### Mapping donor lists (map-report)

```sh
neoHGT gtdb map-report --db "$DB" --in hgt_list.txt \
  --out donor_ncbi.tsv
```

The input is one GTDB label or synthetic taxId per line. The output TSV gives the NCBI taxId mapping; unmapped lines are recorded as unmapped, not as failures. `--to` defaults to `ncbi`.

### Phylogenetic distances (tree-dist)

```sh
neoHGT gtdb tree-dist --gtdb-dir "$GTDB_DIR" --rank genus \
  --focus g__QDMS2 --top-k 20 --out dist.tsv --report dist.json
```

`--tree` can name a tree file directly. The program computes patristic distances and writes a distance-sorted taxon table plus the distance distribution. With `--gtdb-dir` it tries the two domain trees for that rank and selects the one containing the focal taxon.

Measured on r220 genus-level trees:

| Focus | Taxa in tree | Median distance | ≤0.1 | ≤0.3 | ≤0.5 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `g__QDMS2` (archaea) | 1,846 | 3.00 | 0 | 1 | 27 |
| `g__Halorubrum` (archaea) | 1,846 | 3.03 | 1 | 10 | 37 |
| `g__Escherichia` (bacteria) | 23,111 | 2.36 | 33 | 111 | 158 |

Using tree distances as a close group:

```sh
neoHGT gtdb tree-dist --gtdb-dir "$GTDB_DIR" --rank genus \
  --focus g__QDMS2 --max-distance 0.5 --db "$DB" \
  --out close.tsv --close-tax-out close.txt
neoHGT analyze -i search/sample.tsv -t "$DB/taxdump" -o "$OUT" \
  --close-tax "$(cat close.txt)" --self-rank genus --donor-rank genus
```

`--db` resolves labels to synthetic taxIds through the library's `gtdb/labels.tsv` and keeps only taxa actually present in the library. On the 564-genome archaeal library, 21 of the 27 genera at ≤0.5 had no hits and were dropped with warnings; using the remaining taxa as the close group raised predictions from 371 (9.2% of proteins) to 849 (21.1%), while all 53 ar53 conserved marker genes still scored 0.

## Hybrid mode

`--source hybrid` grafts the top of the NCBI tree onto the GTDB tree:

- GTDB's `d__Bacteria` and `d__Archaea` are re-parented under the NCBI cellular root (taxId 131567).
- Top-level eukaryote and virus nodes are kept; other NCBI nodes do not enter the tree.
- `--hybrid-ncbi-taxdump` names a directory containing `nodes.dmp`/`names.dmp` and is required.
- `--hybrid-lca-policy root|domain|refuse` decides how proteins shared across sources are handled:
  - `root`: keep every candidate and place the LCA at the cellular root;
  - `domain`: prefer the candidate's own domain node; multi-candidate or cross-domain proteins are recorded in `x_source_proteins`;
  - `refuse`: discard candidates whose taxId is outside the GTDB range (`>= --gtdb-taxid-base`); when none remain, write `na`.

```sh
neoHGT database --source hybrid --gtdb-dir "$GTDB_DIR" \
  --hybrid-ncbi-taxdump /path/to/ncbi_taxdump \
  --hybrid-lca-policy domain --reps --seq-source local \
  --compile diamond -o "$NEOHGT_CACHE/db/hybrid"
```

The merged tree keeps a single root. On the r220 reference copy, the pure GTDB tree has 146,380 nodes; grafting 25 NCBI nodes produces 146,404 nodes, a net increase of 24 because one grafted node is the root itself. Both domains have parent 131567, and the bacteria-versus-eukaryote LCA falls on the cellular root.

DIAMOND validates the rank column of `nodes.dmp`. Current NCBI taxdumps use `cellular root`, `acellular root`, `realm`, `clade` and other modern ranks that DIAMOND rejects. `write_taxdump()` writes those ranks as `no rank` in `nodes.dmp`, while the in-memory tree, `gtdb/lineages.tsv` and `gtdb/labels.tsv` keep the real ranks. Flattened nodes are listed in `meta.json` under `taxonomy.hybrid.ranks_flattened_for_diamond` and in `gtdb/build_report.md`.

A hybrid library currently searches GTDB proteins only (`--seq-source local`). Grafted NCBI nodes provide tree semantics and contribute no hits; searching eukaryotic or viral proteins requires an NCBI reference protein library.

## Output files and meta.json

### Library directory layout

```text
<db>/db.faa
<db>/genomes.tsv
<db>/lineages.txt
<db>/taxon.map.gz
<db>/genome.map.gz
<db>/meta.json
<db>/taxdump/{nodes.dmp,names.dmp}
<db>/diamond/db.dmnd
<db>/gtdb/taxid.map.gz
<db>/gtdb/{labels.tsv,lineages.tsv,genome2taxid.tsv.gz,stem2taxid.tsv.gz,
           prefix2genome.tsv.gz,ncbi_map.tsv.gz,build_report.md,
           input_fingerprints.tsv,extraction_state.json,fetch_plan.tsv,
           members.tsv.gz,delnodes.txt,reserved_ids.txt}
```

`fetch_plan.tsv` appears only with `--seq-source ncbi|gtdb-tar|ena`; `members.tsv.gz` only with `--expand-members yes`; `delnodes.txt` and `reserved_ids.txt` only when an old ID manifest is reused and nodes changed.

### meta.json

| Field | Contents |
| --- | --- |
| `tool` / `source` / `seq_source` | Version, data source, sequence channel |
| `gtdb_release` / `local_layout` | Release label and layout |
| `id_convention` / `taxonomy` | Synthetic taxId rules, starting base, rank counts, superseded nodes |
| `extraction_resumed` | Whether this run reused an existing `db.faa` |
| `input_provenance` | Input root, read-only audit, taxonomy/metadata paths, sequence directories, `input_sha256` |
| `counts` | Genomes, proteins, duplicates, hits without taxId, taxa, `db.faa` size, longest sequence ID |
| `diamond_version` / `taxid_map_sha256` | DIAMOND version and taxonmap digest |
| `selection` | Selection method and actual counts |
| `selfcheck` / `selfcheck_items` | Self-check status and per-item results |
| `config_snapshot` / `build_started` / `build_finished` | Effective configuration and timestamps |

### Provenance header

Derived tables under `gtdb/`, `lineages.txt`, `input_fingerprints.tsv`, `fetch_plan.tsv`, `members.tsv.gz`, `delnodes.txt`, `reserved_ids.txt` and the reconcile lists all start with:

```text
# neoHGT version=2.0b3+gtdb1; gtdb_release=220; layout=local-renamed; generated=<timestamp>
```

The version comes from `neoHGT.__version__` and changes with the installed package.

`taxon.map.gz` and `genome.map.gz` deliberately carry no comments: DIAMOND, BLAST and `read_prot2taxid` consume them verbatim. The `--input-tax` file produced by `gtdbtk` likewise cannot carry comment lines. Readers of internal `gtdb/` tables skip `#` lines.

## Resuming and release upgrades

### Resuming

The extraction stage records counters in `<db>/gtdb/extraction_state.json`. Running the same command again into the same `-o` directory:

- If `db.faa`, the state file and the selection fingerprint all exist and agree, extraction is skipped and the derived artifacts are rebuilt from `taxon.map.gz`/`genome.map.gz`.
- If the selection fingerprint changed, the old `db.faa` is not reused and extraction runs again.
- If `db.faa` exists but `taxon.map.gz` is missing, the command errors out and suggests `--overwrite`.
- Use `--overwrite` to force a fresh extraction.

### Upgrading to a new release

1. Run `probe` and `reconcile` on the new dataset. `n_checks` depends on the dataset structure; if `idlist_vs_taxa[*]` or `rep_genomes_have_local_fasta` fails, the copy is incomplete and you should stop.
2. Build the new taxdump with the old library's `gtdb/taxid.map.gz`:
   ```sh
   neoHGT gtdb taxdump --dir /path/to/new_release \
     --gtdb-release 232 \
     --gtdb-id-manifest /path/to/old_db/gtdb/taxid.map.gz \
     --out "$NEOHGT_CACHE/db/r232_taxdump"
   ```
3. Pass `--gtdb-id-manifest` to the real build as well, and keep `--gtdb-taxid-base` unchanged.
4. After building, check `taxonomy.rank_counts` and `n_superseded` in `meta.json`, and spot-check that old nodes kept their taxIds.
5. Labels may change between releases. Reinterpret `g__` labels in historical `hgts/*.txt` against the new `gtdb/labels.tsv`, and regenerate `gtdb/ncbi_map.tsv.gz`.
6. Library size affects the usable `close_size` range; recalibrate after switching libraries (see [Appendix](#appendix-parameter-calibration-and-measurements)).

TaxIds of deleted taxa are kept in `reserved_ids.txt` and are never reused. In a bidirectional r214/r220 check, shared node IDs changed zero times, with zero collisions and zero reuse.

## Resource budgets

| Tier | Genomes | Proteins | `db.faa` | `db.dmnd` | Wall clock |
| --- | ---: | ---: | ---: | ---: | --- |
| 6-genome smoke | 6 | 14,373 | 6.2 MB | 614 MB | about 1 min |
| Family sample, archaea | 564 | 905,829 | 360 MB | 985 MB | 2 min 37 s |
| One family, `--reps` | 793 | 3,247,172 | 1.41 GiB | 2.07 GiB | 4 min 50 s |
| Genus sample, bacteria (extrapolated) | 23,112 | 60–80 million | 28–37 GB | 42–110 GB | extrapolated |
| Full representative set (extrapolated) | 113,104 | 333,384,550 | about 145 GiB | about 1.5–3× `db.faa` | extraction about 2–3 h + makedb 6–15 h |

The `db.faa` size baseline is **466 bytes per protein** (three measured tiers: 455 / 417 / 466). The capacity precheck estimates from this value and applies a 3× peak factor for `makedb`.

Two DIAMOND constraints:

1. `db.dmnd` grows with the numeric range of taxIds, not with the node count. The same 13 MB library gives 5.7 GB at `--gtdb-taxid-base 1e9` and 620 MB at `1e8`. The default is therefore `100000000`, and values below `10000000` are rejected.
2. The rank column of `nodes.dmp` must use the NCBI vocabulary and all 13 columns must be present. This implementation maps GTDB `domain → superkingdom` and handles modern NCBI ranks automatically.

Measured cross-volume read bandwidth: 60 MB/s single stream and 127 MB/s with 8 concurrent streams. Tiered builds are recommended (archaeal family first, then genus, then the full set), and `--keep-faa compress` frees space after a build.

## Known limitations

- GTDB contains bacteria and archaea only, with no eukaryotes or viruses. Use `--source hybrid` when those are needed.
- A hybrid library currently searches GTDB proteins only; NCBI nodes provide tree semantics only.
- Remote sequence channels (`--seq-source ncbi|gtdb-tar|ena`) only export `gtdb/fetch_plan.tsv` for manual execution; this version has no download executor. `--allow-network yes` is a prerequisite confirmation for those channels, not a download switch.
- Lineages embedded in file names stop at genus level; species-level assignments come from metadata/taxonomy.
- A local GTDB copy may have been renamed or have had sequence IDs rewritten. The tool records input sha256 and runs self-checks, but no official MD5SUM is available for comparison.
- `--keep-faa drop` deletes `db.faa`. Use `compress` if the sequences must be kept.
- GTDB taxIds are synthetic and must not be mixed with NCBI taxIds in the same historical result tree. Cross-release stability depends on `gtdb/taxid.map.gz`.
- The official flat layout can be used for `taxdump`, but a full build still requires metadata. Complete builds have been run only on a local copy that carries metadata.
- The BLAST route is not usable with GTDB IDs: `makeblastdb` caps local IDs at 50 characters, while GTDB sequence IDs are 151–270 characters.
- A head-to-head comparison against the NCBI route has not been completed; the performance and grouping numbers here come from GTDB-internal measurements.
- Genus-sample and full-set sizes and times are extrapolated, not executed.

## Troubleshooting

| Symptom | Cause | Action |
| --- | --- | --- |
| `no metadata table ...` | The dataset has taxonomy files only | A build needs `--gtdb-metadata`; `taxdump` alone can continue |
| `does not have 7 ranks` | Lineage strings truncated or using a non-standard separator | Use the official taxonomy files |
| `representative ... has no local FASTA` | Representative `.faa` missing | Add the sequence files, or use `--seq-source gtdb-tar\|ena` to generate a fetch plan |
| `Refusing to write ... outside workspace` | `NEOHGT_ROOT` is set and the output is outside it | Move the output inside the workspace or adjust `NEOHGT_ROOT` |
| `Refusing to write ... input dataset is read-only` | The output path is inside `--gtdb-dir` | Choose another output directory |
| `ambiguous option` | A long-option abbreviation matches several options | Spell out the full long option |
| `No hit is assigned to distal group` | The close group is too large and swallowed every hit | Reduce `--close-size` or enlarge the reference library |
| `Invalid taxonomic rank` | The NCBI taxdump uses modern ranks | Use a taxdump generated by this tool; it flattens unsupported ranks for DIAMOND |
| `makeblastdb` fails or IDs are truncated | GTDB IDs exceed 50 characters | Use DIAMOND |
| `taxon.map.gz` missing | The build did not finish or artifacts were moved | Rebuild with `--overwrite` |
| `db.faa` missing | `--keep-faa drop` was used | Rebuild, or unpack `db.faa.gz` from a `compress` build |
| Out of memory | Full-set build or NCBI extraction is memory-hungry | Build in tiers; on the NCBI route use `--stream-extraction yes` |
| Out of disk | Capacity precheck failed | Free space, or pass `--force` explicitly (not recommended) |

## Appendix: parameter calibration and measurements

### Streaming extraction (NCBI route)

`--stream-extraction yes` replaces the upstream protein extraction with a disk-backed implementation. It applies to the NCBI route only and is off by default. Measured at 200k and 1M proteins:

| Implementation | 200k proteins | 1M proteins | Marginal cost | Extrapolated to 333.4M proteins | Time for 1M | Temp disk |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Upstream | 365 MB | 1,071 MB | about 882 MB per million | about 287 GB | 20.3 s | 0 |
| Streaming | 193 MB | 304 MB | about 138 MB per million | about 45 GB | 45.9 s | 104 MiB |

Memory drops about 6.4×, at the cost of about 2.3× extraction time and about 100 MiB of temporary disk per million proteins. `db.faa`, `genome.map.gz` and `taxon.map.gz` are byte-identical to the upstream implementation.

### `--taxmap` disk index

`--taxmap` uses a disk index with one anchor per 2048 rows (about 0.2% of the data) and resolves queries by binary search plus a local scan. Measured at 10,000 rows with spacing 512, it creates 20 anchors and a 21-line index file; a second open reuses the cache. The mapping file must be plain text and byte-sorted by its first column; `.gz` is rejected.

### Parameter sensitivity

- `--bandwidth` is one of the strongest knobs and its effect is sample-dependent. On the halo sample, `auto`/`grid`/`silverman` give 371/145/222 predictions; on bathy they give 367/346/366.
- The direction of `--noise` is not guaranteed to be monotonic. halo gives 374/371/325 at 25/50/75, while bathy gives 285/421 at 25/75; on bimodal or long-tailed distributions the direction can flip.
- `--input-cov` mainly changes the input classification, not the number of predictions. When the hit classification is clean, 50–95 is usually flat; near the coverage threshold, raising it can demote the input classification from species to genus.
- `--from-scores` agrees with a full analysis on the prediction set and ranking; floating-point columns may differ in the last digit because of `scores.tsv` precision.

### Test fixtures

`neoHGT/tests/data/GTDB/` is a KB-scale fixture in the repository (about 118 KiB) with protein excerpts from six genomes, 30 `db.faa` records, 32 taxdump nodes and 168 hit rows. It can be used directly as a read-only searchable library; the offline end-to-end tests are described in `CONTRIBUTING.md`.
