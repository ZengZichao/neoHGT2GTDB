# Command-line reference

`neoHGT` has four sub-commands: `database`, `search`, `analyze`, and `gtdb`. This page lists every option accepted by each sub-command. The page is maintained by hand against the shipped `arguments` tables and is checked against them by the traceability gate (`tools/verify_plan_traceability.py`); when this page and `neoHGT <sub-command> --help` disagree, the code is authoritative.

Run `neoHGT <sub-command> --help` for the same information in a terminal.

## Conventions

- `-o/--output` is required where shown; relative paths are resolved against the current directory.
- Values in the `Default` column are argparse defaults. A dash means the option has no argparse default; the effective value then comes from `config.yml` or from an auto-detection rule. See [config-reference.md](config-reference.md).
- Long options may be abbreviated only when the abbreviation is unambiguous. New GTDB options make several short abbreviations ambiguous; spell out the full option in scripts (see the README's "Incompatible CLI shorthand").
- `--source` has different value sets at different stages: `database` accepts `ncbi|gtdb|hybrid` (no `auto`), while `search` and `analyze` accept `auto|ncbi|gtdb|hybrid`.

### `neoHGT database`

| Option | Default | Description |
| --- | --- | --- |
| `--default` | — | apply the default protocol |
| `-o`, `--output` | — | output database directory |
| `-c`, `--cats` | microbe | enter one or more of the following RefSeq genome categories, separated by comma: archaea, bacteria, fungi, invertebrate, plant, protozoa, vertebrate_mammalian, vertebrate_other, viral; "all" for all categories, or "microbe" for archaea, bacteria, fungi, and protozoa |
| `-t`, `--taxids` | — | provide a custom list of TaxIDs; only genomes under these taxonomic groups will be included |
| `-g`, `--genoids` | — | provide a custom list of genome IDs, which are NCBI assembly accessions (e.g., "GCF_000123456.1"), with or without version number, or simply "G000123456" |
| `-e`, `--exclude` | — | exclude instead of include TaxIDs and genomes defined by custom lists |
| `-s`, `--sample` | — | sample up to this number of genomes per taxonomic group at given rank (0 for none, omit for all) |
| `-r`, `--rank` | species | taxonomic rank at which sampling will be performed |
| `--above` | — | sampling will also be performed on ranks from the designated one to phylum |
| `--genbank` | — | also search GenBank if a genome is not found in RefSeq; otherwise only include RefSeq genomes |
| `--complete` | — | only include complete genomes or chromosomes |
| `--reference` | — | include NCBI-defined reference genomes |
| `--represent` | — | include NCBI-defined representative genomes |
| `--typemater` | — | include NCBI-defined type material genomes |
| `--capital` | — | organism name must be capitalized (choices: yes, no) |
| `--block` | — | ignore organism names that contain any of these words |
| `--latin` | — | genomes must have Latinate species names, e.g. "Vibrio cholerae", but not "Vibrio sp. 123" (choices: yes, no) |
| `--manual` | — | export URLs of sampled genomes without downloading them, so the user may download them manually, then resume the database pipeline. |
| `--overwrite` | — | overwrite existing files with newly downloaded ones; othewise use existing files whenever available, i.e., resume an interrupted run |
| `--retries` | — | number of trials for downloading each data file |
| `--delay` | — | seconds between two retries |
| `--timeout` | — | seconds seconds before program gives up waiting |
| `--compile` | — | compile database using this program if available (choices: diamond, blast, both, none) |
| `--diamond` | — | diamond executable |
| `--makeblastdb` | — | makeblastdb executable |
| `--threads` | — | number of threads for diamond makedb |
| `--tmpdir` | — | temporary directory for diamond makedb |
| `--source` | — | data source for building the database (default: the config file's database.source, else ncbi). "hybrid" grafts an NCBI classification skeleton onto the GTDB tree but ingests no NCBI protein: the searchable sequences are the GTDB ones (choices: ncbi, gtdb, hybrid) |
| `--gtdb-dir` | — | root directory of a local GTDB dataset, which is opened read-only; required when --seq-source is "local" |
| `--local-layout` | — | layout of the local GTDB dataset; "auto" detects it from file names, sequence ids and compression ("gtdb-r220-renamed" is accepted as an alias of "local-renamed") (choices: auto, local-renamed, gtdb-r220-renamed, official-plain, official-gz) |
| `--gtdb-release` | — | GTDB release label written to meta.json (inferred from file names when omitted) |
| `--gtdb-metadata` | — | explicit GTDB metadata TSV path(s), comma-delimited |
| `--gtdb-taxonomy` | — | explicit GTDB taxonomy TSV path(s), comma-delimited |
| `--gtdb-taxid-base` | — | first taxId of the reserved GTDB range (default 100000000; DIAMOND sizes db.dmnd by the highest taxId, so lower is much smaller, but it must stay above real NCBI taxIds) |
| `--gtdb-id-manifest` | — | previous gtdb/taxid.map.gz, so that taxIds stay stable across releases |
| `--gtdb-no-selfcheck` | — | skip the input self-consistency checks (not recommended: the build would not be reproducible) |
| `--seq-source` | — | where protein sequences are read from (choices: local, ncbi, gtdb-tar, ena) |
| `--allow-network` | — | acknowledge that the remote sequence channels (ncbi, gtdb-tar, ena) were requested; this version has no download executor, so both values only export gtdb/fetch_plan.tsv for manual use and then refuse (choices: yes, no) |
| `--gtdb-mirror` | — | GTDB data host used by the gtdb-tar channel (choices: ecogenomic, aau, auto) |
| `--downloader` | — | downloader named in the exported fetch plan (this version performs no transfers itself) (choices: curl, wget) |
| `--verify-md5` | — | md5_check column of the exported fetch plan (MD5 verification itself is not implemented in this version) (choices: yes, no) |
| `--workspace-out` | — | override workspace.cache_dir for the output tree |
| `--workspace-tmp` | — | override workspace.tmp_dir for intermediate files |
| `--reps` | — | restrict to GTDB species-cluster representative genomes (equivalent to -r species -s 1) |
| `--genome-list` | — | file(s) of GTDB accessions to include, one per line; accepts the bundled *_Sample_<rank>.IDlist files |
| `--rank-sample` | — | keep one genome per taxonomic group at this rank (cheap tier for development and CI) (choices: domain, phylum, class, order, family, genus) |
| `--gtdb-taxa` | — | include only genomes under these GTDB taxa; accepts prefixed labels (f__Enterobacteriaceae) or GTDB taxIds |
| `--min-completeness` | — | minimum Checkm2 completeness (%%) required |
| `--max-contamination` | — | maximum Checkm2 contamination (%%) allowed |
| `--mimag` | — | require a MIMAG quality tier (choices: high, medium, low, any) |
| `--type-material` | — | which type-material definition to use (choices: gtdb, ncbi) |
| `--domain` | both | limit the build to one GTDB domain (choices: bacteria, archaea, both) |
| `--expand-members` | — | grow a representative selection to every member genome of the species clusters it stands for (needs non-representative sequence files, which the local r220 copy does not contain) (choices: yes, no) |
| `--keep-placeholder-names` | — | keep GTDB placeholder group names such as g__DUMJ01 (recommended, they are valid taxon labels) (choices: yes, no) |
| `--hybrid-ncbi-taxdump` | — | directory with NCBI nodes.dmp/names.dmp used by --source hybrid. A hybrid database still contains GTDB proteins only: the NCBI dump grafts a classification skeleton, it adds no searchable sequences |
| `--hybrid-lca-policy` | root | how to treat proteins shared between GTDB and NCBI subtrees (choices: root, domain, refuse) |
| `--keep-faa` | keep | what to do with db.faa once the DIAMOND database is built (choices: keep, compress, drop) |
| `--input-fingerprints` | — | also record sha256 of the input files used (recommended; the local GTDB copy ships no MD5SUM.txt) (choices: yes, no) |
| `--force` | — | proceed even when the disk-space pre-flight fails (read by the GTDB route: it also unlocks the experimental BLAST compile and publishes builds with missing genomes; the NCBI route ignores it) |
| `--max-genomes` | — | build at most this many genomes (protects the shared disk during development) |
| `--limit` | — | alias of --max-genomes |
| `--taxdump-trim` | — | write only the taxa ancestral to the selected genomes instead of the full GTDB tree (smaller files, but taxIds are then specific to this selection) |
| `--stream-extraction` | — | extract proteins of the NCBI route without keeping the protein-to-genome / protein-to-taxid maps in RAM: dedup through a disk-backed index and external sorting in --tmpdir. Produces byte-identical db.faa / genome.map.gz / taxon.map.gz; needs extra temporary disk (~30 bytes per protein occurrence). (choices: yes, no) |

### `neoHGT search`

| Option | Default | Description |
| --- | --- | --- |
| `-i`, `--input` | — | input protein set file, or directory where one or more input files are located |
| `-o`, `--output` | — | directory where search results are to be saved |
| `-m`, `--method` | auto | search method (choices: auto, diamond, blast, remote, precomp) |
| `-s`, `--precomp` | — | file or directory of precomputed search results (when method = precomp) |
| `-d`, `--db` | — | reference protein sequence database |
| `-t`, `--taxdump` | — | directory of taxonomy database files (nodes.dmp and names.dmp) |
| `--taxmap` | — | sequence Id to taxId mapping file (not necessary if protein database already contains taxonomy) |
| `-k`, `--maxhits` | — | maximum number of hits to preserve per query (0 for unlimited) |
| `--minsize` | — | minimum length of query sequence (aa) |
| `--queries` | — | number of queries per run (0 for whole sample) |
| `--maxchars` | — | maximum number of characters per run (0 for unlimited) |
| `--maxseqs` | — | maximum number of sequences to return |
| `--evalue` | — | maximum E-value cutoff |
| `--identity` | — | minimum percent identity cutoff |
| `--coverage` | — | minimum percent query coverage cutoff |
| `--extrargs` | — | extra arguments for choice of search method |
| `--tax-include` | — | include taxa under those taxIds |
| `--tax-exclude` | — | exclude taxa under those taxIds |
| `--tax-unique` | — | ignore more than one hit with same taxID (choices: yes, no) |
| `--tax-unirank` | — | ignore more than one hit under same taxon at this rank |
| `--tax-capital` | — | ignore taxon names that are not capitalized (choices: yes, no) |
| `--tax-latin` | — | ignore species names that are not Latinate (choices: yes, no) |
| `--tax-block` | — | ignore taxon names containing those words |
| `-p`, `--threads` | — | number of threads (0 for all CPU cores) |
| `--tmpdir` | — | temporary directory |
| `--diamond` | — | diamond executable |
| `--blastp` | — | blastp executable |
| `--blastdbcmd` | — | blastdbcmd executable |
| `--algorithm` | — | remote search algorithm |
| `--retries` | — | maximum number of retries per search |
| `--delay` | — | seconds between two search requests |
| `--timeout` | — | seconds before program gives up waiting |
| `--entrez` | — | entrez query text |
| `--server` | — | remote search server URL |
| `--aln-method` | — | self-alignment method (choices: auto, native, fast, lookup, precomp) |
| `--aln-precomp` | — | file or directory of precomputed sequence Id to score maps (when self-alignment method = precomp) |
| `--aln-server` | — | remote align server URL |
| `--fetch-enable` | — | whether to enable remote fetch (choices: auto, yes, no) |
| `--fetch-queries` | — | maximum number of query entries per search |
| `--fetch-retries` | — | maximum number of retries per search |
| `--fetch-delay` | — | seconds between two fetch requests |
| `--fetch-timeout` | — | seconds before program gives up waiting |
| `--fetch-server` | — | remote fetch server URL |
| `--source` | auto | data source of the reference database; "auto" detects a GTDB database from its gtdb/ subdirectory (choices: auto, ncbi, gtdb, hybrid) |
| `--gtdb-stem2taxid` | — | explicit gtdb/stem2taxid.tsv.gz used as the fallback (route C) when the database was built without embedded taxIds |
| `--keep-unresolved` | — | keep hits whose taxId cannot be resolved, instead of dropping them (choices: yes, no) |
| `--require-staxids` | — | fail immediately when the hit table lacks the staxids column (choices: yes, no) |
| `--keep-placeholder-names` | — | exempt GTDB placeholder labels such as g__DUMJ01 or "MK2026 sp010841923" from --tax-capital/--tax-latin; they are official labels, not mis-capitalized or non-Latinate names (choices: yes, no) |

### `neoHGT analyze`

| Option | Default | Description |
| --- | --- | --- |
| `-i`, `--input` | — | input search result file, or directory where one or more input files are located |
| `-o`, `--output` | — | directory where analysis results will be saved |
| `-t`, `--taxdump` | — | directory of taxonomy database files (nodes.dmp and names.dmp); required if they are not found in input directory |
| `-k`, `--maxhits` | — | maximum number of sequence similarity search hits per gene (protein) to preserve |
| `--evalue` | — | maximum E-value cutoff |
| `--identity` | — | minimum percent identity cutoff |
| `--coverage` | — | minimum percent query coverage cutoff |
| `--input-tax` | — | taxIds of input samples, format: sample1:taxId1,sample2:taxId2..., will auto-infer if omitted |
| `--input-cov` | — | for auto-inference: an input sample will be assigned to a taxon if it represents this percentage or more best hits |
| `--self-tax` | — | taxIds of "self" group, delimited by comma, will auto-infer if omitted |
| `--close-tax` | — | taxIds of "close" group, delimited by comma, will auto-infer if omitted |
| `--self-rank` | — | for auto-inference: "self" group must be at or above this rank |
| `--close-size` | — | for auto-inference: "close" group must have at least this number of taxa |
| `--weighted` | — | score is sum of weighted bit scores; otherwise simple counts (choices: yes, no) |
| `--outliers` | — | detect and remove outliers using selected method (choices: none, zscore, boxplot) |
| `--orphans` | — | keep orphans (proteins without non-self hits) in statistical analysis (choices: yes, no) |
| `--bandwidth` | — | bandwidth for Gaussian KDE (auto, grid, silverman, or a number between 0.1 and 1.0) |
| `--bw-steps` | — | steps for auto and grid kernel bandwidth optimization |
| `--low-part` | — | maximum percentage below threshold for automatic bandwidth optimization |
| `--noise` | — | percent valley-to-peak distance to exclude from cluster |
| `--fixed` | — | use this percentage as threshold if KDE clustering fails |
| `--silhouette` | — | silhouette score threshold for cluster refinement |
| `--self-low` | — | HGT has low "self" score (an optional criterion) (choices: yes, no) |
| `--distal-top` | — | find match in "distal" group which is LCA of hits with bit score at most this percentage lower than best hit |
| `--donor-name` | — | report taxon name instead of taxId of donor. |
| `--donor-rank` | — | report donor at this rank (genus, phylum, etc.). |
| `--from-scores` | — | if score table already exists, use it and skip search result parsing and taxonomy inference; otherwise overwrite it |
| `--source` | auto | data source of the reference database; "auto" detects a GTDB database from its gtdb/ subdirectory (choices: auto, ncbi, gtdb, hybrid) |
| `--donor-label` | — | how to express the donor taxon in hgts/*.txt: taxId, bare name, or the official prefixed GTDB label (g__Escherichia) (choices: taxid, name, prefixed) |
| `--donor-ncbi-map` | — | also write <sample>.donor_ncbi.tsv mapping donors to NCBI taxIds (offline, from gtdb/ncbi_map.tsv.gz) (choices: yes, no) |
| `--taxon-synonyms` | — | TSV of deprecated -> current GTDB names (alias<TAB>current); GTDB renames higher taxa between releases and the official synonyms.*.tsv is not in a local release copy |

### `neoHGT gtdb probe`

| Option | Default | Description |
| --- | --- | --- |
| `--dir` | — | GTDB dataset root (read-only) |
| `--out` | — | Destination `layout_probe.json` |
| `--gtdb-release` | — | Expected release label; inferred from file names when omitted |
| `--local-layout` | `auto` | Layout override: `auto`, `local-renamed`, `official-plain`, `official-gz` |
| `--gtdb-metadata` | — | Explicit metadata TSV path(s), comma-delimited |
| `--gtdb-taxonomy` | — | Explicit taxonomy TSV path(s), comma-delimited |
| `--no-selfcheck` | — | Skip the dataset self-checks (not recommended) |
| `--no-deep` | — | Skip the per-accession alignment checks |

### `neoHGT gtdb reconcile`

| Option | Default | Description |
| --- | --- | --- |
| `--dir` | — | GTDB dataset root (read-only) |
| `--out` | — | Destination directory for the five mismatch lists |
| `--report` | — | Destination JSON summary of the counts |
| `--gtdb-release` | — | Expected release label |
| `--local-layout` | `auto` | Layout override |

### `neoHGT gtdb taxdump`

| Option | Default | Description |
| --- | --- | --- |
| `--dir` | — | GTDB dataset root (read-only) |
| `--out` | — | Destination directory |
| `--gtdb-release` | — | Expected release label |
| `--local-layout` | `auto` | Layout override |
| `--gtdb-taxid-base` | `100000000` | First synthetic taxId; minimum `10000000` |
| `--gtdb-id-manifest` | — | Previous `gtdb/taxid.map.gz`, to keep taxIds stable |
| `--gtdb-metadata` | — | Explicit metadata TSV path(s) |
| `--gtdb-taxonomy` | — | Explicit taxonomy TSV path(s) |
| `--no-selfcheck` | — | Skip the dataset self-checks |

### `neoHGT gtdb map-report`

| Option | Default | Description |
| --- | --- | --- |
| `--db` | — | Database directory built with `neoHGT --source gtdb` |
| `--in` | — | Input list with one donor label or synthetic taxId per line |
| `--out` | — | Destination TSV |
| `--to` | `ncbi` | Target taxonomy; v1 implements `ncbi` only |

### `neoHGT gtdb gtdbtk`

| Option | Default | Description |
| --- | --- | --- |
| `--in` | — | GTDB-Tk classify or de novo summary table |
| `--db` | — | Database built with `neoHGT --source gtdb` |
| `--out` | — | Destination `genome<TAB>taxId` file |
| `--report` | — | Destination JSON report |
| `--min-rank` | — | Drop genomes assigned below this rank: `domain`, `phylum`, `class`, `order`, `family`, `genus`, `species` |

### `neoHGT gtdb tree-dist`

| Option | Default | Description |
| --- | --- | --- |
| `--gtdb-dir` | — | GTDB release root; uses the trees under `tree_*/` |
| `--tree` | — | A specific Newick tree file |
| `--rank` | genus | Rank tree to use with `--gtdb-dir`: `phylum`, `class`, `order`, `family`, `genus` |
| `--focus` | — | Label of the focal taxon, e.g. `g__Halorubrum` |
| `--max-distance` | — | Maximum patristic distance to report |
| `--top-k` | — | Report the top K closest taxa |
| `--out` | — | Destination TSV of `label<TAB>distance` |
| `--report` | — | Destination JSON summary |
| `--db` | — | Also resolve labels to synthetic taxIds using `gtdb/labels.tsv` |
| `--close-tax-out` | — | Write resolved taxIds as a comma-separated list for `analyze --close-tax` |

## Legacy flat form

The older form `neoHGT gtdb --action probe|taxdump|map-report` is still accepted. New commands should use `neoHGT gtdb <sub-command>`.
