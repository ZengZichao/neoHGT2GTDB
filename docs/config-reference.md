# `config.yml` reference

neoHGT2GTDB reads defaults from `config.yml`. The program looks for the file in this order:

1. the current directory;
2. `~/.neoHGT/config.yml`;
3. the directory of the installed package.

Command-line options override configuration values. For the GTDB-specific options, the command line also overrides `database.gtdb`.

## Workspace

| Key | Type | Meaning |
| --- | --- | --- |
| `workspace.code_dir` | path or empty | Source directory used to derive the workspace root (its parent) |
| `workspace.tmp_dir` | path or empty | Scratch directory for builds and searches |
| `workspace.cache_dir` | path or empty | Cache/reports directory |

The workspace is optional. When no root is configured, outputs may be written anywhere and only the read-only GTDB input is protected. When a root is configured, every output must resolve inside it. The `NEOHGT_ROOT`, `NEOHGT_TMP` and `NEOHGT_CACHE` environment variables take precedence over these keys. `--workspace-out` and `--workspace-tmp` override the output and scratch directories for one `database` run.

## `database`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `source` | `ncbi`/`gtdb`/`hybrid` | `ncbi` | Data source when `database` is run without `--source` |
| `hybrid.ncbi_taxdump` | path or empty | empty | Directory with NCBI `nodes.dmp`/`names.dmp` for `--source hybrid` |
| `gtdb.dir` | path or empty | empty | Read-only GTDB release root; `$GTDB_DIR` and `--gtdb-dir` are alternatives |
| `gtdb.release` | string or `auto` | `auto` | Release label; inferred from file names when `auto` |
| `gtdb.local_layout` | `auto`/`local-renamed`/`official-plain`/`official-gz` | `auto` | Layout detection override; `gtdb-r220-renamed` is accepted as an alias of `local-renamed` |
| `gtdb.seq_source` | `local`/`ncbi`/`gtdb-tar`/`ena` | `local` | Sequence channel |
| `gtdb.taxid_base` | integer | `100000000` | First synthetic taxId; must be at least `10000000` |
| `gtdb.selfcheck` | `required` | `required` | Run self-checks and abort on failure |
| `gtdb.verify_input_fingerprints` | `yes`/`no` | `yes` | Record input sha256, size and line counts |
| `gtdb.mirror` | `ecogenomic`/`aau`/`auto` | `ecogenomic` | Host used by the `gtdb-tar` channel |
| `diamond` | path or empty | empty | DIAMOND database path for the NCBI route |
| `blast` | path or empty | empty | BLAST database path for the NCBI route |
| `taxdump` | path or empty | empty | NCBI `nodes.dmp`/`names.dmp` directory |
| `taxmap` | path or empty | empty | Sequence ID to taxId mapping for the NCBI route |

Quality thresholds are deliberately not read from this section. Pass `--min-completeness`, `--max-contamination`, `--mimag` or `--type-material` per build, so the selection in `meta.json` matches the requested filters.

## `program`

| Key | Meaning |
| --- | --- |
| `diamond` | DIAMOND executable |
| `blastp` | `blastp` executable |
| `blastdbcmd` | `blastdbcmd` executable |
| `makeblastdb` | `makeblastdb` executable |

An empty value means "find the program on `PATH`".

## `server`

| Key | Meaning |
| --- | --- |
| `search` | NCBI BLAST search endpoint |
| `selfaln` | NCBI self-alignment endpoint |
| `fetch` | NCBI Entrez fetch endpoint |

## `search`

| Key | Default | Meaning |
| --- | --- | --- |
| `method` | `auto` | `auto`, `diamond`, `blast`, `remote` or `precomp`; `auto` prefers `precomp` > `diamond` > `blast` > `remote` |
| `minsize` | `30` | Minimum query protein length in amino acids |
| `maxseqs` | `500` | Maximum number of sequences returned |
| `evalue` | `1e-5` | Maximum E-value |
| `identity` | `0` | Minimum percent identity |
| `coverage` | `0` | Minimum percent query coverage |
| `maxhits` | `0` | Maximum hits preserved (`0` for unlimited) |
| `selfaln` | `auto` | Self-alignment method: `auto`, `native`, `fast`, `lookup`, `precomp` |

## `local`, `diamond`, `blast`

| Key | Default | Meaning |
| --- | --- | --- |
| `local.threads` | `0` | Threads (`0` for all cores) |
| `local.tmpdir` | empty | Temporary directory |
| `diamond.queries` | `0` | Queries per DIAMOND run (`0` for the whole sample) |
| `diamond.maxchars` | `0` | Maximum characters per run (`0` for unlimited) |
| `diamond.extrargs` | empty | Extra DIAMOND arguments |
| `blast.queries` | `100` | Queries per BLAST run |
| `blast.maxchars` | `100000` | Maximum characters per BLAST run |
| `blast.extrargs` | empty | Extra `blastp` arguments |

## `remote`

| Key | Default | Meaning |
| --- | --- | --- |
| `db` | `refseq_select_prot` | Remote database |
| `algorithm` | `kmerBlastp` | Remote BLAST algorithm |
| `queries` | `0` | Queries per run |
| `maxchars` | `7000` | Maximum characters per run |
| `retries` | `5` | Retries per search |
| `delay` | `60` | Seconds between requests |
| `timeout` | `1800` | Seconds before giving up |
| `entrez` | see file | Entrez query text |
| `extrargs` | see file | Extra URL arguments |

Remote search is not used by the GTDB route.

## `fetch` and `download`

| Key | Default | Meaning |
| --- | --- | --- |
| `fetch.enable` | `auto` | Enable remote taxonomy fetch: `auto`, `yes`, `no`; the GTDB route forces `no` |
| `fetch.queries` | `100` | Maximum query entries per fetch |
| `fetch.retries` | `3` | Retries per fetch |
| `fetch.delay` | `5` | Seconds between fetches |
| `fetch.timeout` | `60` | Seconds before giving up |
| `download.retries` | `3` | Retries per downloaded file |
| `download.delay` | `10` | Seconds between download retries |
| `download.timeout` | `60` | Seconds before giving up |

## `taxonomy`

| Key | Default | Meaning |
| --- | --- | --- |
| `include` | empty | Include only these taxIds and their descendants |
| `exclude` | empty | Exclude these taxIds and their descendants |
| `unique` | `yes` | Ignore more than one hit with the same taxId |
| `unirank` | empty | Ignore more than one hit under the same taxon at this rank |
| `capital` | `yes` | Ignore taxon names that are not capitalized |
| `latin` | `no` | Ignore species names that are not Latinate |
| `block` | list | Ignore taxon names containing any of these words |
| `profiles.gtdb.capital` | `no` | GTDB names carry rank prefixes; style filters are off by default |
| `profiles.gtdb.latin` | `no` | See above |
| `profiles.gtdb.block` | empty list | No blocked words for GTDB |
| `profiles.gtdb.keep_placeholder` | `yes` | Exempt official placeholder labels from style filters |
| `profiles.ncbi.capital` | `yes` | Upstream NCBI behaviour |
| `profiles.ncbi.latin` | `no` | Upstream NCBI behaviour |

## `analyze`

| Key | Meaning |
| --- | --- |
| `maxhits` | Maximum hits per gene for analysis |
| `evalue` | Analysis E-value cutoff |
| `identity` | Analysis identity cutoff |
| `coverage` | Analysis coverage cutoff |

These keys provide a second, stricter filter without repeating the search.

## `grouping`

| Key | Default | Meaning |
| --- | --- | --- |
| `inputcov` | `75` | Minimum percentage of best hits needed to assign the input genome to a taxon |
| `selfrank` | empty | Minimum rank of the self group |
| `closesize` | `10` | Minimum number of taxa in the close group |
| `profiles.gtdb.closesize` | `200` | Close-group size for GTDB bacteria |
| `profiles.gtdb.inputcov` | `75` | Input coverage for GTDB bacteria |
| `profiles.gtdb_archaea.closesize` | `50` | Close-group size for GTDB archaea |
| `profiles.gtdb_archaea.inputcov` | `75` | Input coverage for GTDB archaea |

The close group is a lower bound: when the tree cannot supply that many taxa it is truncated. The effective size is printed as `Close group has N taxa.`

## `predict`

| Key | Default | Meaning |
| --- | --- | --- |
| `weighted` | `yes` | Weight genes by bit score instead of counting occurrences |
| `outliers` | `zscore` | High-end outlier detection: `none`, `zscore`, `boxplot` |
| `orphans` | `no` | Keep proteins without non-self hits |
| `bandwidth` | `auto` | KDE bandwidth: `auto`, `grid`, `silverman`, or a number between 0.1 and 1.0 |
| `bwsteps` | `20` | Steps for bandwidth optimization |
| `lowpart` | `75` | Maximum low-end quantile for automatic bandwidth selection |
| `noise` | `50` | Percentage of valley-to-peak distance excluded from the atypical cluster |
| `fixed` | `25` | Fallback atypical quantile when KDE clustering fails |
| `silhouette` | `0.5` | Silhouette threshold for cluster refinement |
| `selflow` | `no` | Require the self score to be atypically low |

`noise` is not guaranteed to be monotonic across samples: on bimodal or long-tailed distributions, increasing it can loosen rather than tighten the prediction. Tune it against your own sample's curve.

## `donor`

| Key | Default | Meaning |
| --- | --- | --- |
| `distaltop` | `10` | LCA is taken over distal hits within this percentage of the best distal bit score |
| `name` | `no` | Report taxon names instead of taxIds |
| `rank` | empty | Report the donor at this rank |

On GTDB databases, `analyze --donor-label prefixed` writes the official prefixed GTDB label, and `--donor-ncbi-map yes` additionally writes a donor-to-NCBI-taxId table.

## Examples

Use a larger close group for a dense bacterial GTDB library:

```yaml
grouping:
  profiles:
    gtdb:
      closesize: 300
      inputcov: 75
```

Keep GTDB style filters off but block a custom word:

```yaml
taxonomy:
  profiles:
    gtdb:
      capital: no
      latin: no
      block:
        - synthetic
```

Set a workspace and a local GTDB root:

```yaml
workspace:
  code_dir: /path/to/neoHGT2GTDB
  tmp_dir: /path/to/workspace/tmp
  cache_dir: /path/to/workspace/cache

database:
  gtdb:
    dir: /path/to/GTDB_r220
```
