> English (primary) | Chinese mirror: [GTDB-数据库.md](../GTDB-数据库.md) — keep the two in lockstep. Chinese mirror linked above.

# GTDB database (experimental)

neoHGT2GTDB supports building the reference library under the [GTDB](https://gtdb.ecogenomic.org/) taxonomy — this is where this adaptation deeply improves on upstream
[neoHGT2](https://github.com/cgneo/neoHGT) (whose command-line entry is likewise `neoHGT`). GTDB provides standardized,
hierarchy-complete prokaryotic classification (exactly one representative genome per species cluster), which markedly improves the
statistical behavior of self/close/distal grouping. GTDB **contains no eukaryotes or viruses**; when needed, use `--source hybrid` (that mode requires obtaining the NCBI
taxdump separately, see the "Known limitations" section).

All numbers on this page are measured on this machine (GTDB R09-RS220 local dataset, 2026-11-01/02).
For the condensed English version, see ["GTDB database (experimental)" in `README.md`](../README.md).

## Workspace conventions (must be set up first)

All artifacts may only be written to three directories, and the dataset pointed to by `--gtdb-dir` is **read-only throughout** (the tool self-checks write
permission at startup and rejects any output path pointing into it):

| Directory | Purpose |
| --- | --- |
| `neoHGT2GTDB-项目代码/` | Source code, tests, documentation |
| `neoHGT2GTDB-临时文件/` (the workspace temp directory, kept under its Chinese name) | `--tmpdir`, DIAMOND intermediate state, disposable |
| `neoHGT2GTDB-缓存测试/` (the workspace cache directory, kept under its Chinese name) | `db/`, `search/`, `analyze/`, `fixtures/`, `logs/`, `benchmarks/`, `reports/` |

First load the environment (it redirects matplotlib/XDG/TMPDIR to inside the project, avoiding writes outside the workspace):

```sh
. tools/env.sh
```

However, GTDB/hybrid database builds **no longer depend on this wrapper**: when `--tmpdir` is not given, upstream opens the scratch directory at
`mkdtemp()` (i.e. wherever the environment TMPDIR points), so the past state of "intermediate artifacts inside the workspace" was really a side effect of env.sh.
Now the `--source gtdb|hybrid` branch moves it under `workspace.tmp_dir` (`临时文件/gtdb_XXXX`),
and, in passing, `rmdir`s the empty directory upstream just created (`rmdir` can only delete empty directories, so it cannot delete anything else); at the end of the build the
whole directory is cleaned up as before. Measured: pointing `TMPDIR` somewhere else and running a 1-genome build without `--tmpdir`, the log prints
`Scratch directory inside the workspace: …/neoHGT2GTDB-临时文件/gtdb_7sh52ard`, the directory used as the
decoy **stays empty throughout**, and no residue is left after the build ends (`缓存测试/db/nomprefail`). The NCBI route's defaults
stay exactly as upstream left them (RT-2: do not change the default behavior of existing options), so an NCBI-route build should still source env.sh first.

The three conventions can be checked mechanically with `tools/verify_workspace.py` (IT-7 / DoD-9):

```sh
python tools/verify_workspace.py --out "$NEOHGT_CACHE/reports/IT7_workspace.json" \
    --since "2026-11-01 17:25:00" --probe-report \
    "$NEOHGT_CACHE/reports/IT7_tmp_disposable_probe.json"
```

It checks four things: (1) outside the three project directories there are no files **produced by this project** (pre-existing files older than the start time are
listed separately under `predates_project` and do not count as violations); (2) there are no symlinks inside the project directories pointing outside the workspace;
(3) the cache directory keeps the subdirectory structure of the table above; (4) an mtime audit of the dataset (reusing `verify_read_only.py`).
The claim "**`临时文件/` can be safely emptied**" is not guessed by reading the manifest: the scratch paths appearing in the manifest are listed only as a **manifest**;
the real verdict comes from a behavioral probe — resume the build with a **brand-new empty `--tmpdir`**, then have
`tools/verify_resume.py` compare against the first build file by file (PASS proves nothing depends on the old scratch contents).

## Read-only pre-flight (recommended before every build)

```sh
neoHGT gtdb probe --dir "$GTDB_DIR" --out "$NEOHGT_CACHE/logs/layout_probe.json"
# 等价入口（不经 CLI 包装）：
python -m neoHGT.gtdb probe --dir "$GTDB_DIR" --out "$NEOHGT_CACHE/logs/layout_probe.json"
```

The two entry points are equivalent: `neoHGT gtdb <子命令>` merely forwards the arguments as-is to `neoHGT.gtdb.main()`
(exit codes included; the command-line entry `neoHGT/gtdb.py::cli_main` prints `ValueError`/`SelfCheckError` such as "invalid input"
as a single `ERROR: …` line and returns exit code 2, instead of dumping a
traceback — `main()` itself still raises as before, so library callers are unaffected); the subcommands are
`probe / taxdump / map-report / gtdbtk / reconcile / tree-dist`; the earlier flat form
`neoHGT gtdb --action probe|taxdump|map-report` still works.

It outputs the layout label, release, sequence ID convention and anchor coverage, and runs **19** self-consistency checks
(faa file count vs the metadata representative-set count, taxonomy line count vs metadata line count, lineage depth always 7,
the rank-prefix vocabulary, the line counts of the six `*_Sample_<rank>.IDlist` tiers vs the number of unique taxa, the domain sets,
accession alignment, and whether representative genomes have local sequences). **Any mismatch aborts the build** (it can be skipped explicitly with
`--gtdb-no-selfcheck`, but reproducibility is then no longer guaranteed).

Measured on local r220: layout `gtdb-r220-renamed`, `_SerialNumber_` anchor coverage 1.00,
longest defline 270 characters, 19/19 checks passed, took 28 s.

## Writing the mismatch lists to disk (DoD-2)

`probe` only tells you whether the **counts** match; to get "which accessions mismatch" onto disk, use
`reconcile`:

```sh
python -m neoHGT.gtdb reconcile --dir "$GTDB_DIR" \
    --out "$NEOHGT_CACHE/reconcile" --report "$NEOHGT_CACHE/reports/DoD2_reconcile.json"
```

It computes set differences among three sides — the full metadata table, the GTDB taxonomy, and the local `.faa` files — and writes out five
TSVs (each one has a `#` comment header line stating the source / release / generation time plus an `accession`
header line, so **it leaves a trace even when empty**); the exit code is 1 on inconsistency:

| File | Meaning |
| --- | --- |
| `missing_faa.tsv` | accessions marked as representatives in metadata but with no local sequence file |
| `extra_faa.tsv` | accessions with a local `.faa` that GTDB metadata does not recognize at all |
| `nonrep_faa.tsv` | local files of **non-representative** (member) genomes; should be empty under the default layout |
| `metadata_not_in_taxonomy.tsv` | in metadata but not in `*_taxonomy.tsv` |
| `taxonomy_not_in_metadata.tsv` | in taxonomy but not in metadata |

Measured on local r220 (`reports/DoD2_reconcile.json`): metadata 596,859 genomes
(of which 113,104 representatives) / 113,104 local FASTA / 596,859 taxonomy lines, **all five lists
are 0**, `consistent: true`; took 24.3 s, peak memory 1.7 GB.

## Building a standalone taxdump

```sh
python -m neoHGT.gtdb taxdump --dir "$GTDB_DIR" --out "$NEOHGT_CACHE/db/r220_taxdump"
```

Measured: **146,380 nodes** (GTDB 146,379 + root); rank counts = domain 2 / phylum 217 /
class 607 / order 2,030 / family 5,460 / genus 24,959 / species 113,104;
conversion 42.9 s, peak memory 1.5 GB; artifacts 192 MB.

**Datasets with only lineage files can also produce a taxdump**. The taxonomy tree is a function of the lineage strings alone and does not need the
metadata table; after an official release is unpacked, `bac120/ar53_taxonomy_rNNN.tsv(.gz)` sits directly in the dataset root,
but the original implementation's discovery logic only looked inside `tree_bac120/` and `tree_ar53/`, so the "official layout" was ironically not recognized. Now:

```sh
# 真实数据实测：本地 r214 副本（无 metadata 表、谱系文件在根目录）
neoHGT gtdb taxdump --dir /path/to/GTDB_r214/parent_directory \
    --gtdb-id-manifest <r220 库>/gtdb/taxid.map.gz --out <新目录>
```

It prints `NOTICE: no metadata table ...` and `lineage_only: true`, runs only the checks that work without metadata
(lineage depth always 7; one accession must not have two lineages), and explicitly lists the 6 skipped §3.6 checks; measured
113,318 nodes / 6.3 seconds / `official-plain`, and it also writes `gtdb/{taxid.map.gz,delnodes.txt,reserved_ids.txt}`
(`缓存测试/gtdb_r214_taxdump_v2/`; in the first run's `缓存测试/gtdb_r214_taxdump/`, that "skipped" list had **made-up names** —
the node count and mappings remain valid; see item 10 of "Self-discovery" in this file).
The database build (`neoHGT database --source gtdb`) still **requires** the metadata table, and the error will name `--gtdb-metadata`.

Malformed lineage strings are rejected **before any writes** — the local r226 copy is exactly such input (a third party glued the ranks together with `_`,
truncated them to genus level, and left the `.faa` suffix on the accession column):

```
[FAIL] taxonomy_files_readable: .../ar53_taxonomy_r226.tsv:1: lineage of RS_GCF_945873965.1.faa
       does not have 7 ranks (provenance check failed).
```
Exit code 1, and the output directory is not created. For the measured numbers of the checks and of both directions (old mapping + new lineages, old mapping + even older lineages), see
`缓存测试/reports/R24_lineage_only_and_release_upgrade_SOP.json`.

## Building databases

```sh
# 极小档实测（r220 古菌每科 1 个代表，564 基因组 / 905,829 蛋白，墙钟 2 分 37 秒）
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" --rank-sample family \
  --domain archaea --compile diamond -o "$NEOHGT_CACHE/db/r220_ar_family"

# 谱系子集（示例：肠杆菌科代表集）
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --gtdb-taxa f__Enterobacteriaceae --reps --seq-source local \
  --compile diamond --threads 4 -o "$NEOHGT_CACHE/db/r220_entero" \
  --tmpdir "$NEOHGT_TMP/build/$$"

# 使用外部基因组清单（可直接吃本地现成 *_Sample_<rank>.IDlist）
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --genome-list "$GTDB_DIR/metadata_bac120/bac120_metadata_r220.tsv_Sample_Genus.IDlist" \
  --compile diamond -o "$NEOHGT_CACHE/db/r220_bac_genus"

# 全量代表集（113,104 基因组 / 333,384,550 蛋白，零下载）
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" --reps --seq-source local \
  --min-completeness 50 --max-contamination 10 \
  --compile diamond --keep-faa compress \
  -o "$NEOHGT_CACHE/db/r220_reps" --tmpdir "$NEOHGT_TMP/build/$$"
```

`--max-genomes N` truncates the selection to the first N genomes, for controlled smoke tests.
Quality-class filtering goes through `--min-completeness / --max-contamination / --mimag / --type-material gtdb|ncbi`
(upstream's `--typemater` is just an alias for `--type-material ncbi` on the GTDB route; `--complete` is **rejected** rather than
silently ignored — see the section "Semantics of the quality-class options under GTDB" below for the real counts of each column).
When `--gtdb-dir / --seq-source / --local-layout / --gtdb-taxid-base` are not given explicitly, they are taken from
`database.gtdb` in `config.yml` (`dir` being the local read-only dataset root); the command line always wins.

The sampling trio `-s/--sample`, `-r/--rank`, `--above` has the same semantics on the GTDB route as on the NCBI route
(wired in this update: previously `--above` never made it into the provider's opts and was silently ignored under GTDB):

- First take at most n genomes per `<rank>` taxon;
- Genomes whose lineage string **does not reach that rank** (e.g. identified only to genus `g__G2` with no `s__`) sit on this library
  out; the program prints "how many entries were not identified to that rank" and hints at `--above`;
- With `--above` given, starting from one level above the sampling rank and going upward level by level (stopping at phylum; consistent with upstream, the domain is not used as a sampling tier),
  each coarser taxon is topped up to n — this is exactly the channel by which the "shallow lineage" genomes above enter the library.
- `--rank-sample <rank>` is another low-cost tier: take only the single highest-quality-scoring genome per taxon, ignoring `--above`.

One honest measured conclusion: in r220 the **representative lineages always have 7 segments down to species** (the
depth histogram of `gtdb_taxonomy` in `genomes.tsv` is `{7: N}`), and "n genomes per fine taxon" necessarily makes every coarser taxon
satisfy n as well (unless the candidates are exhausted), so on full-depth data `--above` is a **provably no-op**:
the same family with the same set of parameters measured 235 vs 235 genomes. Where it does act is on inputs with shallower lineages
(GTDB-Tk results, member-genome expansion), and this path is locked by a unit test
(`test_above_tops_up_coarser_groups`: a taxon identified only to genus gets in only with `--above`).

On the `search` side, `--keep-placeholder-names` is wired up (default `yes`, from
`taxonomy.profiles.gtdb.keep_placeholder` in `config.yml`): 1,658 taxa in GTDB have
numbered names (`2-12-FULL-45-22`, `g__DUMJ01`, `s__MK2026 sp010841923`), and these **are** official
labels; if the user explicitly turns on `--tax-capital` / `--tax-latin`, `util.is_placeholder_gtdb()`
exempts such labels, otherwise tens of thousands of hits would silently vanish. `--keep-placeholder-names no` returns to the strict semantics.
The exemption takes effect only on GTDB databases, and it swaps in the predicate only for that one filtering call, restoring it when the call ends.

> This piece of semantics was fixed in two places (the original implementation was broken in both situations):
> - **`--tax-latin yes` under GTDB used to necessarily drop every hit**. Upstream `taxfilt_wf()` first parses out the species-level
>   **taxId**, then hands that taxId to `is_latin()` as a "name" — a digit string is never a binomial Latin name, so
>   every hit was judged "not Latin". The GTDB route now swaps that taxId for its own name before the test
>   (`neoHGT/search.py::taxfilt_wf_gtdb`); the NCBI route is untouched word for word (RT-1 line-by-line comparison passes),
>   and this upstream defect is recorded on row 3.5 of internal QC ledger 《计划任务对照表》 (workspace-side, not part of the public release).
> - **Strict semantics does not mean rejecting prefixed official names**. `--keep-placeholder-names no` now only cancels the numbered-label exemption;
>   the predicate still uses the prefix-tolerant version (`s__Halobacterium` still counts as capitalized), because DoD-5 requires judging "by the bare name with the prefix stripped".
>   Test cases: `neoHGT/tests/test_gtdb_pipeline.py::TestSearchChoiceValues` and
>   `TestPlaceholderExemption`.
>
> The same round also fixed `--require-staxids`: it was parsed, converted by `arg2bool`, and stored on the object, **but no code ever read it**,
> so `no` had no effect at all; worse, because `arg2bool(None)` is `False`, "the user never mentioned this option" was treated as "the user does not want it".
> Now the GTDB/hybrid routes default to `--require-staxids yes`: a hit table missing the staxids column errors out immediately (RK5's protection);
> only an explicit `--require-staxids no` downgrades this to a WARNING and continues (classification then rests on the subject Id alone);
> the NCBI route restores upstream semantics (a missing column follows upstream behavior), and this default no longer applies to NCBI runs.

`--compile none` also got a fix: after writing `db.faa`/taxdump/side tables, `write_gtdb_meta()`
handed the not-yet-located `self.diamond=None` to `shutil.which()`, raising `TypeError` and failing the entire build at
the last step — the library "looked built" but had no `meta.json`. It now checks for None first, records `diamond_version` in
`meta.json` as `null`, and `TestMetaWithoutDiamond` locks this in (including the other half, "the version is still
recorded faithfully when DIAMOND exists").

The item-by-item correspondence between the plan §6's numbered tasks and the implementation is written in
internal QC ledger 《计划任务对照表》 (workspace-side, not part of the public release) and checked by the traceability gate (IDs complete,
statuses drawn from a fixed vocabulary, referenced paths actually existing, and "partially implemented / not done" entries required to state clearly what is missing).

The 14 "pre-existing defects shared with this adaptation" listed in the plan's Appendix F are registered in
internal QC ledger 《既有缺陷登记》 (workspace-side, not part of the public release): each entry is either **fixed** (with evidence attached) or
**registered, not changed**, stating clearly why it was not changed (mostly "fixing it would require deleting existing upstream lines", which conflicts with RT-1's append-only constraint)
and how to reproduce it. The same gate also checks this table, and line-number references are verified as well (a form like `search.py:1166`
pointing at a nonexistent line gets reported).

The plan's §9 risk register (19 entries) is checked in
internal QC ledger 《风险应对登记》 (workspace-side, not part of the public release): each entry records "where the mitigation planned in the plan finally landed", with status drawn from
`已实现缓解 / 部分实现 / 登记接受`; today's score is 16 with mitigation implemented and 3 partially implemented
(RK7's known-HGT sensitivity validation, RK9's remote resume execution, RK19's "NFS→NFS output" fallback,
stuck respectively on NCBI reference-library authorization, network-access authorization, and the "must not write outside the workspace" safety rule).

The eighth and ninth tables: the 9 NCBI assumptions that the plan §2.3 wrote down line by line for `search.py` are checked in
internal QC ledger 《搜索层假设核对》 (workspace-side, not part of the public release) (5 adapted for GTDB, 3 confirmed to need no change,
1 **contrary to what the plan concluded at the time** — the plan considered the BLAST route "feasible", but measurement showed `makeblastdb` cannot build indexes for
GTDB sequence IDs of 151–270 characters, so `--compile blast/both` was changed to refuse by default and only continue with `--force`, see ADR-8); the 19 review revisions of the plan §0.5 are checked in
internal QC ledger 《方案修订记录核对》 (workspace-side, not part of the public release) (18 adopted, 1 partially adopted: REV-12's fixture source and
§10.1 coverage points ⑥⑦ were not copied over, written out item by item in the deliverables document's checklist). These two families are likewise
"numbered, truth-checkable" plan content that previously had no register at all.

The tenth table is the §13 milestone check (internal QC ledger 《里程碑核对》, workspace-side, not part of the public release).
Beyond these tables, the gate itself has now also entered CI: `docs/方案编号清单.json` is a digest freezing every numbered list in the plan
(276 items: 28 tasks / 14 defects / 12 deliverables / 19 risks / 10 DoD / 11 decisions / 36 verification-matrix entries /
9 search-layer assumptions / 19 review revisions / 37 test IDs / 27 §7 options / 27 Appendix C/D options / 27 Appendix B symbols),
so **the checks can run even without the plan document in the repository**. CI's "Register gate" step runs the ten tables against this digest;
recomputing the digest is only possible locally with the real plan in hand (`--plan <方案> --emit-ids docs/方案编号清单.json`),
and the test `test_frozen_plan_digest_equals_the_plan_it_came_from` asserts locally that the digest matches the plan field by field,
going red the moment the digest goes stale. Additionally, this tool itself refuses to treat paths outside the workspace as output targets
(`guarded_out()`, rule 11: `--out /tmp/x` is not a redirect and is not intercepted by `tools/run.sh`, so it has to do its own intercepting).

The plan's §10.2 / §10.3 / §10.4 / §10.5 verification matrix (36 numbered entries: 20 UT, 9 IT, 4 scientific-validity items, 3 RT) is checked in
internal QC ledger 《测试用例对照》 (workspace-side, not part of the public release): each entry states which test case (with `文件:行号`) or which evidence artifact
covers it, with status drawn from `已覆盖 / 部分覆盖 / 未覆盖`. Currently 33 covered, 0 partially covered, 3 uncovered (IT-3's data-supplementation route, and the 10.4-1 and 10.4-4 head-to-head NCBI reference-library comparisons, all stuck on network access not being authorized). §10.4-3's parameter sensitivity has measured curves: all four knobs — close_size, input_cov, bandwidth, noise — have measured curves, see `缓存测试/reports/S10_4_3_sensitivity.md`. This table grew out of an inventory: back then it was discovered that
the derived fixtures required by §10.1 did not exist, and that the GTDB branches of `search`/`analyze` had never had an **offline** end-to-end test —
laying the verification matrix out as a machine-checkable table exists precisely so that "written in the plan but actually absent from the code" never silently happens again.

Among these, RK11 was completed in this update: GTDB renames higher taxa (`p__Proteobacteria` →
`p__Pseudomonadota`), and the official `synonyms.*.tsv` is not in the local copy, so
`analyze` gained `--taxon-synonyms <alias<TAB>current 的 TSV>`; all five references —
`--self-tax`,
`--close-tax`,
`--tax-include`,
`--tax-exclude`,
`--input-tax` — pass through the alias table before resolution, and when an old name is used, a "this is a deprecated name, X is now used instead" hint is printed.

The 10 acceptance criteria of the plan §1.4/§10.7 (DoD-0…DoD-9) are cross-checked in
internal QC ledger 《验收标准对照》 (workspace-side, not part of the public release): status drawn from `达标 / 部分达标 / 被外部条件阻塞`,
each entry with an evidence path. Current state: 9 met + 1 partially met (DoD-2: the representative set's **genome count** is reconciled by a build-time enforcement item,
while the **protein entry count** is a metadata sum plus extrapolation from measured tiers, because the full extraction was left to the user's decision and was not run).
The same round folded §3.6's three rows of whole-library file counts (`.faa`/`.fna`/`.tsv`) into the read-only audit:
`tools/verify_read_only.py` also emits `files_by_extension`, measured
`faa=146686 fna=113360 tsv=112972`, equal item by item to the numbers in the §3.6 table
(`缓存测试/logs/read_only_check_r12.txt`). These three rows depend on the copy's specific contents, so they are treated as
"leave a trace on every audit" rather than a "build-time equality assertion"; `--no-inventory` turns this inventory off.

The 11 design decisions of the plan §8 (ADR-1…ADR-10 including 4b) land in
internal QC ledger 《决策落点对照》 (workspace-side, not part of the public release): status drawn from
`已落实 / 部分落实 / 偏离并说明 / 未落实`, each entry with its code location. Current state: 10 implemented + 1 partially implemented
(ADR-4's "v1 only promises the representative set" holds in code, but its run-through evidence only goes up to the 3.25M-protein tier; the full tier was not run).

The plan's §12 deliverables (D1–D12) and the **interface shapes of Appendices B/C/D** are checked in
internal QC ledger 《交付物清单》 (workspace-side, not part of the public release): each of the 12 deliverables gets a status and an artifact location (D11's full-tier benchmark
and D12's NCBI head-to-head comparison are honestly marked as partially delivered); the 27 functions/classes/constants promised by Appendix B must be findable
in the code (being located in another module of the same package is allowed), and any that cannot be found must be registered with an explanation in that file's "Appendix B symbol deviations" table;
the 27 `--flag`s in the Appendix C/D example commands must each be either a real option of this tool's CLI, an option of a third-party tool,
or recorded in the "Appendix C/D command-shape deviations" table.

> The second command above (the Enterobacteriaceae one) is a measured one: `f__Enterobacteriaceae` has exactly **793** representative
> genomes in r220 (153 genera), `genomes.tsv` likewise has 793 rows, matching metadata one by one; no quality thresholds were layered on
> (`--min-completeness/--max-contamination` must be given explicitly, see above).

> **Semantics of the quality-class options under GTDB**: the local r220 metadata full table has **596,859** rows
> (bacteria 584,382 / archaea 12,477), of which **113,104** rows are the representative set. All counts below were measured by the **production selection code**
> `neoHGT/gtdb.py::GtdbProvider.select()` over this full table (`缓存测试/reports/R29_quality_class_counts.json`),
> not sampling estimates — the first time these were recorded, only the first 400 thousand rows of the table were read, yielding "only 14,745 entries / 3.7% carry a designation";
> that was **sampling bias** (type material concentrates toward the end of the file), and it has been replaced with the full-table numbers.
>
> - `--type-material gtdb` → `gtdb_type_species_of_genus == 't'`: **4,545** rows in the full table, **3,347** rows among representatives.
> - `--type-material ncbi` (upstream's `--typemater` flag is exactly its alias in GTDB mode) → `ncbi_type_material_designation`
>   with `na`/empty strings necessarily excluded: **22,256** rows in the full table (`na` 574,603), **16,371** rows among representatives = 14.5%, in 5 kinds of
>   designation: type material 16,313 / synonym type material 28 / reftype 19 / neotype 9 / pathotype 2.
> - `--mimag high|medium|low|any` → `mimag_*_quality == 't'`: among representatives 27,292 / 85,752 / 8, the three tiers mutually exclusive;
>   `any` (the fourth choice in §7, meaning "must carry some MIMAG tier") keeps **113,052 / 113,104** among representatives and
>   596,497 / 596,859 in the full table, i.e. only 52 / 362 genomes have "no tier at all". A missing column raises an error instead of silently returning the full table:
>   GTDB-Tk input tables and releases that do not publish MIMAG classes simply lack these columns, and `--mimag` on such tables was a
>   **silent no-op**.
> - `--complete` under GTDB **errors out directly** instead of silently doing nothing: upstream looks at `ncbi_assembly_level`, which on the r220 representative set
>   keeps only 8,593 (Chromosome 1,070 + Complete Genome 7,523), while the 104,511 discarded genomes still have a median completeness of
>   94–96% — assembly level describes "how the sequence was submitted", not "how complete the genome is" (`ncbi_refseq_category`
>   likewise cannot be used to judge quality: 87.4% of the representative set is `na`). Plan §7 R5 explicitly replaces it with `--min-completeness/--max-contamination`;
>   measured on the representative set, completeness≥50 keeps 82,662, contamination≤10 keeps 82,464, and both together 82,161.
> - The build report prints these classes' **actual counts** on a separate line (`quality classes in the dataset: …`), because plan §7 requires
>   the report to give real numbers rather than making the user guess; the counts are tallied on the **unfiltered** table, so it describes the dataset, not this run's filtering.

### Database directory layout

```
<db>/db.faa  genomes.tsv  lineages.txt  taxon.map.gz  genome.map.gz  meta.json
<db>/taxdump/{nodes.dmp,names.dmp}     # 全量 GTDB 树（rank 用 NCBI 词表）
<db>/diamond/db.dmnd
<db>/gtdb/taxid.map.gz                 # 跨版本 taxId 稳定性凭据（只追加不重排）
<db>/gtdb/{labels.tsv,lineages.tsv,genome2taxid.tsv.gz,stem2taxid.tsv.gz,
          prefix2genome.tsv.gz,ncbi_map.tsv.gz,build_report.md,
          input_fingerprints.tsv,extraction_state.json,fetch_plan.tsv,
          members.tsv.gz,delnodes.txt,reserved_ids.txt}
```

`fetch_plan.tsv` appears only with `--seq-source ncbi|gtdb-tar|ena`; `members.tsv.gz` appears only with
`--expand-members yes`; `delnodes.txt` and `reserved_ids.txt` appear only when "the previous version's
manifest is reused and nodes really did disappear".

`meta.json` records the GTDB release, `local_layout`, `id_convention`, input provenance (including
the list of locally adapted third-party scripts and sha256 fingerprints), the DIAMOND version, the taxonmap digest, the results of the 19 self-consistency checks,
and statistics (genomes/proteins/nodes/db.faa size/longest sequence ID).

**Provenance header of derived files**: every table under `gtdb/`, plus `lineages.txt`,
`gtdb/input_fingerprints.tsv`, `gtdb/fetch_plan.tsv`, `gtdb/members.tsv.gz`,
`gtdb/delnodes.txt`, `gtdb/reserved_ids.txt`, and the `reconcile` mismatch lists, all have as their first line

```
# neoHGT version=2.0b3; gtdb_release=220; layout=gtdb-r220-renamed; generated=<本地时间>[; 附加字段]
```

Thus any artifact can be traced back to "which version of the code, which release, generated when".
The two exceptions are **intentional**: `taxon.map.gz` and `genome.map.gz` are consumed as-is by DIAMOND, BLAST and
`read_prot2taxid`, and the `--input-tax` file (the output of `gtdb gtdbtk`) is parsed as-is by
`dict_from_param`, so they cannot carry comment lines. The functions that read the internal `gtdb/` tables
(`load_labels`, `load_manifest`, `load_genome2taxid`, `load_stem2taxid`,
`load_prefix2genome`, `load_ncbi_map`) all skip `#` lines; `load_genome2taxid`
now explicitly locates the header row and reports a clear error when a column is missing.

### Resuming an interrupted build (IT-6)

After each successful step, the extraction phase writes its counters into `<db>/gtdb/extraction_state.json`. When run again with the **same command**
into the same `-o` directory:

- If `db.faa` exists and the state file exists **and the selection fingerprint matches** → **skip extraction**, and rebuild the DIAMOND taxonmap (two columns with identical values), the BLAST seq2taxid, and the genome map from `taxon.map.gz` /
  `genome.map.gz`;
  the protein count, dedup count and no-taxId count in `genomes.tsv` / `meta.json` are exactly those of a single successful build.
- The state file records the sha256 fingerprint of this run's genome selection; if **the selection changed, the old `db.faa` will not be reused**
  (otherwise it would produce a library inconsistent with `genomes.tsv`), and the program automatically switches to re-extracting.
- If `db.faa` exists but `taxon.map.gz` is missing → error out immediately and hint at rebuilding with `--overwrite`;
  metadata inconsistent with `db.faa` is never published.
- Use `--overwrite` when re-extraction is needed (it deletes the old `db.faa` and re-streams the extraction).

To compare a "single successful build" with a "resumed build", use `tools/verify_resume.py`: files with identical sha256 are recorded
`identical`; gz artifacts like `taxon.map.gz`/`genome.map.gz` whose lines match after stripping comments are recorded
`text_equal`; text tables differing only in provenance timestamp lines are recorded `timestamp_only`;
`meta.json` is allowed to differ in only four items — `build_started/build_finished/extraction_resumed` and the
`tmpdir` echoed from the command line (a location, not build content). The comparison is **streaming** — `gtdb/lineages.tsv`
is already 96 MB at the 6-genome tier, and larger at full size, so it cannot be read whole into memory. Any other difference is
`mismatched`, and the verdict is FAIL.

## meta.json field reference and the GTDB release-upgrade SOP

### meta.json fields

| Field | Contents |
| --- | --- |
| `tool` / `source` / `seq_source` | version strings, data source (`gtdb`/`hybrid`), sequence channel (`local`/`ncbi`/`gtdb-tar`/`ena`) |
| `gtdb_release` / `local_layout` | the release version number and the layout determination result (`gtdb-r220-renamed` / `official-plain` / `official-gz`) |
| `id_convention` / `taxonomy.id_base` / `taxonomy.rank_counts` / `taxonomy.n_superseded` | synthetic taxId rules, the starting base, node counts per rank, number of superseded nodes |
| `extraction_resumed` | whether this run reused an existing `db.faa` (IT-6 traceability) |
| `input_provenance` | the read-only input root, `source_writable_by_user` (the result of `os.access(W_OK)`) and `read_only_policy` (stating that enforcement is `guard_write` + a post-build mtime audit, not file permissions), absolute paths of the taxonomy/metadata files, sequence directories, the list of locally adapted third-party scripts, missing official auxiliary files, `_SerialNumber_` anchor coverage, `input_sha256` (file → sha256) |
| `counts` | `genomes{,_planned,_missing}`, `proteins`, `proteins_planned` (the protein count the genomes selected this run should have per metadata), `proteins_extracted_vs_planned` (`ok` / `N short` / `not extracted yet`, obtained by `neoHGT/gtdb.py::proteins_reconciliation` wrapping `neoHGT/gtdb.py::reconcile_protein_counts`; the two fields derive from the same dict — an early version unpacked it by position and wrote the extracted count into the "expected" field; the baseline is **the genomes selected this run** (passed in explicitly as `selection=df` by `neoHGT/database.py::write_gtdb_meta`), not the whole release's table), `proteins_duplicated`, `proteins_without_taxid`, `taxa`, `db_faa_bytes`, `max_seqid_len`, `by_domain_selected`, `representatives_in_dataset` (whole-dataset scope, used for the DoD-2 three-way counts) |
| `diamond_version` / `taxid_map_sha256` | the DIAMOND version used for compilation, the taxonmap digest |
| `selection` | how this selection was made (`--reps`/`--genome-list`/`--rank-sample`/`--gtdb-taxa`/quality thresholds) and the actual hit counts |
| `selfcheck` / `selfcheck_items` | the overall status of the 19 checks and each item's `expected/observed/status` |
| `config_snapshot` / `build_started` / `build_finished` | the snapshot of the effective configuration and the timestamps |

Apart from `build_started/build_finished/extraction_resumed`, the `meta.json` of two identical builds
must be identical field by field (enforced by `tools/verify_resume.py`).

### Upgrading to a new release (e.g. r220 → r232)

1. **Probe first**: `neoHGT gtdb probe --dir <新数据集> --gtdb-release 232 --out ...`.
   The 19 checks do not hard-code r220 numbers (the IDlist counts are computed live from the current dataset), so a new release should come out all green;
   if `idlist_vs_taxa[*]` or `rep_genomes_have_local_fasta` fails, the copy is incomplete — **stop**.
   Then run `neoHGT gtdb reconcile --dir <新数据集> --out ... --report ...` once
   (see the "Writing the mismatch lists to disk" section) to leave on-disk evidence for each mismatched accession, and only then decide whether to build.
2. **Reuse the old ID manifest** to keep taxIds stable across versions (append-only, never reshuffled):
   ```sh
   neoHGT gtdb taxdump --dir <新数据集> --gtdb-release 232 \
     --gtdb-id-manifest <旧库>/gtdb/taxid.map.gz --out <新库目录>
   ```
   That command only produces the taxonomy-side artifacts (`<新库>/taxdump/`, `<新库>/gtdb/` (including the new
   `taxid.map.gz`), `<新库>/meta.json`) and does not touch sequences; measured on r220 it took 42.9 s and 1.5 GB of memory.
   `--gtdb-taxid-base` must match the old library (default `100000000`), otherwise the old and new ID spaces would overlap.

   **Deleted taxa keep their taxIds permanently**: the
   `gtdb/reserved_ids.txt` sitting next to the previous version's `gtdb/taxid.map.gz` is read back automatically, and taxa newly added in the new release will **not** take over IDs
   used historically (otherwise the old `scores.tsv`/`hgts/*.txt` would point at a different organism without a single character changing —
   plan RK6). When rebuilding the same release (nothing added or removed), the file is not created and the mapping is identical entry by entry —
   measured on r220: `146,379` nodes and the two builds' `lineage→taxId` fully equal.

   **This SOP has been run through on a real cross-release scenario** (local r214/r220 real data, no network):
   taking r220's `taxid.map.gz` against r214's lineages, the 89,172 shared node IDs changed **zero** times, the 24,145
   r214-only nodes had **zero collisions** with known IDs in the old table and did not take over the ID of any of the 57,207 vanished nodes,
   all 57,207 of which landed in `delnodes.txt`; the reverse direction (r214 mapping + r220 lineages) likewise zero changes / zero collisions /
   zero reuse, with `delnodes.txt` at 24,145 entries. Watch the criterion: allocation is blocked by rank
   (`slot = id_base + rank_ordinal * 1e6`), so a new node landing inside its own block is **correct**;
   "a new ID must be greater than the old maximum" is the wrong criterion (the first version of the check misused it, see
   `缓存测试/reports/ACCEPTANCE.md` §37.2).
   The numbers come from `缓存测试/reports/R24_lineage_only_and_release_upgrade_SOP.json`,
   and the offline test `neoHGT/tests/test_gtdb_core.py::test_release_swap_never_reissues_or_reuses_a_taxid`
   pins the same invariant in CI.
3. **The real database build** likewise carries `--gtdb-id-manifest <旧库>/gtdb/taxid.map.gz`; after building, check
   `meta.json`'s `taxonomy.rank_counts` and `n_superseded`, and spot-check that several old nodes' taxIds are unchanged.
4. **Rebuild the interpretation of search results**: lineage labels change (GTDB renames are very common), so the
   `g__` labels in historical `hgts/*.txt` need to be reinterpreted against the new `gtdb/labels.tsv`; `gtdb/ncbi_map.tsv.gz` must also be regenerated.
5. **Recalibrate the grouping parameters**: a change of library size shifts the usable range of `closesize` (see the measured numbers in the "Analyze" section:
   the 564-genome tier was only calibrated up to close_size=100 — 10/50/100 → 371/371/132 hits; 200/400 were not measured separately at that tier and, per that section's conclusion, only make sense at library density close to the full representative set). After switching libraries, rerun
   `tools/run_hgt_demo.py --close-sizes 10,50,100,200,400` before settling on defaults.
6. **Keep the old library directory** (do not delete it); the new results directory records the `db 路径 + meta.json sha256` for reproducibility.

## Search

The DIAMOND library has taxIds embedded, so `staxids` is directly available; the recommended path remains "manual search + precomputed results":

```sh
diamond blastp --query q.faa --db "$DB/diamond/db" --threads 8 \
  --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids \
  --max-target-seqs 25 --sensitive > hits.tsv
neoHGT search -i q.faa -m precomp -s hits.tsv -t "$DB/taxdump" \
  -o "$NEOHGT_CACHE/search/sample" --tmpdir "$NEOHGT_TMP/search/$$"
```

In GTDB mode the program **forces `--fetch-enable no`** (synthetic taxIds have no corresponding NCBI records),
no longer "patching the taxdump while searching", and when the hit table lacks the `staxids` column it **errors out directly** instead of silently dropping hits.
Hits whose taxId cannot be determined are written to `unresolved_taxids.txt` + `unresolved_hits.tsv`.
In addition, the 8-column standard BLAST `outfmt 6` table is **rejected outright** in GTDB mode: such a table has no
`staxids` column, and its `qstart` column would be misread as a taxId.

> **Note (pre-existing upstream behavior, not introduced by the GTDB branch)**: `search` writes
> `<output>/<sample>.tsv` in **append** mode, so repeated runs into the same `-o` directory write each protein twice,
> and `analyze` counts all of those duplicate rows. For a rerun, use a fresh `-o` directory, or rename the old table first to preserve it.

## Analyze

```sh
neoHGT analyze -i "$NEOHGT_CACHE/search/sample/sample.tsv" -t "$DB/taxdump" \
  -o "$NEOHGT_CACHE/analyze/sample" \
  --self-rank domain --donor-rank genus --donor-label prefixed \
  --donor-ncbi-map yes
```

- `--self-rank/--donor-rank/--tax-unirank` accept `domain` (automatically mapped to the
  `superkingdom` actually stored in the taxdump) as well as the `superkingdom` alias.
- `--donor-label prefixed` outputs the official GTDB labels (`g__Escherichia`).
- `--donor-ncbi-map yes` additionally saves `hgts/<sample>.donor_ncbi.tsv`, mapping donors back to NCBI taxIds
  (source: the metadata's `ncbi_taxid` column; measured coverage on r220 is 100%, fully offline).
- `--bandwidth` / `predict.bandwidth` can be `auto`, `grid`, `silverman`, or a number between 0.1–1.0;
  writing `bandwidth: 1` (a YAML integer) and `--bandwidth 1` now produce the same answer,
  and out-of-range values are rejected with "source + usable range" stated, instead of just tossing out `Invalid bandwidth`.
- The allowed values of `--source` **differ by stage** (the gate checks the documentation examples pair by pair):

  | Subcommand | Allowed values | Default |
  | --- | --- | --- |
  | `database` | `ncbi\|gtdb\|hybrid` | `ncbi` |
  | `search` | `auto\|ncbi\|gtdb\|hybrid` | `auto` |
  | `analyze` | `auto\|ncbi\|gtdb\|hybrid` | `auto` |

  There is no `auto` at the build stage: which kind of library to build is **specified** by the user, with no room to "guess one from the directory";
  search/analyze face an existing library, and only then can detection happen. When writing documentation, do not put `auto` into build examples.
- `--source` takes effect on both analyze routes by the same criteria: `ncbi` never enters the GTDB branch;
  `gtdb` errors out directly when the GTDB marker is not found (naming `gtdb/labels.tsv` and the way out); `hybrid`, when not found, prints a
  WARNING saying this run actually ran the upstream NCBI face; `auto` is the "use it if detected, don't if not" one.
- The GTDB grouping's `closesize` default is raised to 200 (50 for the archaea tier), because the GTDB tree is far denser than NCBI's;
  the domain is determined automatically by a **vote of the domain-level ancestors of each protein's best hit** (requiring a strict majority; a tie counts as undetermined),
  and the log prints `GTDB profile: input domain = archaea -> using "gtdb_archaea"`.
  When `--from-scores` reuses a scores table there is no hit table, so the vote falls back to the `scores.tsv` donor column itself,
  and the two routes therefore give the same tier for the same sample; a NOTE line is also added saying the grouping stage was skipped and
  these two values do not affect this run's results. When the domain cannot be determined the log writes
  `input domain = not determined ... falling back to the default profile "gtdb"`,
  never conflating "could not detect" and "it really is bacteria" into the same sentence. An explicit `--close-size` always wins.
  Adjustable via `grouping.profiles` in `config.yml`.
- **The reference library must span at least multiple families**: measured, using the "Enterobacteriaceae representative set" (793 genomes, a single family) with an E. coli query,
  4,311 proteins all predicted 0 HGTs at `close_size=10/50/100/200`, with no distal group — because the library
  contains no distant taxa at all. Likewise, if the query genome is itself contained in the library, its best hits get absorbed into the self group.
  Only cross-family/cross-order sampling (e.g. `--rank-sample family`) is statistically meaningful.
- Measured (`db/r220_ar_family`, 564 archaeal representative genomes, 4,025-protein query):
  `close_size` 10/50/100 predicted 371/371/132 HGTs respectively, showing that at this tier the archaea tier 50 is usable,
  while 200/400 need library density close to the full representative set to be meaningful.

## Using GTDB-Tk classifications as input (P7)

GTDB-Tk is **not a runtime dependency** (plan N2): this tool only reads its results table, translating "user genomes → GTDB lineages"
into the GTDB library's synthetic taxIds, for use by `analyze --input-tax`.

```sh
neoHGT gtdb gtdbtk --in BINs.classify.max_summary \
  --db "$NEOHGT_CACHE/db/r220_ar_family" \
  --out "$NEOHGT_CACHE/gtdbtk/input_tax.tsv" \
  --report "$NEOHGT_CACHE/gtdbtk/gtdbtk_report.json" [--min-rank genus]

neoHGT analyze -i search/sample.tsv -t "$DB/taxdump" -o "$NEOHGT_CACHE/analyze/sample" \
  --input-tax "$NEOHGT_CACHE/gtdbtk/input_tax.tsv" \
  --self-rank genus --donor-rank genus --donor-label prefixed
```

- Four table shapes are recognized automatically: v2 `classification`, v1 `classify_info`, de novo's
  `gtdbtk.lineage` / `lineage`, and the "7 rank columns" style; v2's `#####` banner is skipped,
  and the `#user_genome` header is kept.
- **Per-row fallback**: when some ranks are unassigned, v2 leaves `classification` empty; the program automatically switches to
  `classify_info` and truncates the segments after NA.
- The taxId is the **deepest matchable prefix** (e.g. if assigned only to genus, the genus node is used); the report JSON contains the
  `matched_rank` distribution and the **reason** for each unmatched genome.
- The taxId source is always the target library's `gtdb/taxid.map.gz`; when the library lacks that lineage it is judged unmapped,
  and IDs are never invented out of thin air.
- Measured (`db/r220_ar_family`, a `taxid.map.gz` with 146,381 nodes): all 6/6 rows of the synthetic gtdbtk table
  mapped successfully; the query genome assigned at genus level → `105000360`, and `analyze --input-tax`
  then ran end-to-end without trouble (`缓存测试/gtdbtk/`).

## Expanding lineage subsets to member genomes (P7, optional)

GTDB publishes only one representative genome per species cluster; the local image stores protein files only for representative genomes.
`--expand-members yes` uses the metadata's `gtdb_genome_representative` column to expand "representatives" into
"all the members they represent" (the official `sp_clusters.tsv` is missing locally; plan §3.4 is exactly what requires constructing the relation this way):

```sh
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --gtdb-taxa f__Enterobacteriaceae --reps --expand-members yes --compile none \
  -o "$NEOHGT_CACHE/db/entero_members"
```

Measured (r220): 793 representatives → **94,393 member genomes**, all written into
`gtdb/members.tsv.gz` (with `representative` and `has_local_fasta` columns).
This copy has **no local protein files for any member genomes** (0/94,393), so the program aborts before extraction
with an actionable hint (pointing out that the official `*_genomic_files_all` tar or `--seq-source gtdb-tar|ena` is needed), rather than quietly
producing a smaller library. Log at `logs/expand_members.log`.

## Phylogenetic-distance-assisted grouping (P7)

The local release bundle also carries the **collapsed trees per rank** (`tree_{bac120,ar53}/{x}_{rank}.tree`);
they provide an objective yardstick for "exactly how many taxa count as close", which is precisely what `closesize` tries to proxy:

```sh
neoHGT gtdb tree-dist --gtdb-dir "$GTDB_DIR" --rank genus \
  --focus g__QDMS2 --top-k 20 --out dist.tsv --report dist.json
# 或直接指定树文件：--tree "$GTDB_DIR/tree_ar53/ar53_r220_genus.tree"
```

The program uses an in-house lightweight Newick parser (no third-party dependencies) to compute **patristic distances**,
outputting a distance-sorted taxon table and the distance distribution (the counts of taxa at `<= 0.01/0.05/0.1/0.2/0.3/0.5`).
In `--gtdb-dir` mode it tries the two domain trees for that rank in turn, automatically selecting the one containing the focal taxon.

Measured on r220 (genus-level tree):

| Focus | Total taxa in that tree | Median distance | ≤0.1 | ≤0.3 | ≤0.5 |
| --- | --- | --- | --- | --- | --- |
| `g__QDMS2` (archaea, the genus of this run's query genome) | 1,846 | 3.00 | 0 | 1 | 27 |
| `g__Halorubrum` (archaea) | 1,846 | 3.03 | 1 | 10 | 37 |
| `g__Escherichia` (bacteria) | 23,111 | 2.36 | 33 | 111 | 158 |

How to read it: at genus level, archaea have only single digits to a few dozen "truly close" taxa (only 27 genera at ≤0.5), while the bacteria side has 158 genera;
this agrees completely with §8's `close_size` sweep — the archaea tier must be much smaller than the bacteria/full tiers,
otherwise the close group eats all the best hits and the distal group is empty.

**Feeding the tree distances straight into `analyze`** (phylogeny-assisted close group):

```sh
neoHGT gtdb tree-dist --gtdb-dir "$GTDB_DIR" --rank genus --focus g__QDMS2 \
  --max-distance 0.5 --db "$DB" --out close.tsv --close-tax-out close.txt
neoHGT analyze -i search/sample.tsv -t "$DB/taxdump" -o "$OUT" \
  --close-tax "$(cat close.txt)" --self-rank genus --donor-rank genus
```

`--db` resolves the labels into synthetic taxIds via the library's `gtdb/labels.tsv`, keeping only taxa that **actually exist among this library's
genomes**. Another layer of protection lives in `analyze`: a taxId given via `--close-tax` that is completely unreachable in the hit table
(`analyze` first prunes the taxonomy tree down to what the hits involve) is **dropped with a warning** rather than
crashing the program while reading the tree; if none are reachable, an error explains the usage. Measured (`db/r220_ar_family`):
of the 27 genera at ≤0.5, 21 had no hits in that 564-genome archive → dropped with warnings, and the close group continued with the remaining taxa.
Rerunning at this tier with the resolved close group: predictions rose from 371 (9.2% of proteins) to 849 (21.1%),
while all 53 ar53 conserved marker genes still scored 0 — i.e. the phylogenetically defined close group did not misreport conserved genes
as HGTs (the same-tier numbers are in the README's GTDB section).

### Mapping donor lists back to NCBI taxIds (`map-report`)

The HGT lists produced by `analyze` are text with one donor taxId or GTDB label per line; to batch-translate them back to NCBI taxIds offline, use
`map-report`: it reads the library directory's `gtdb/ncbi_map.tsv.gz`
(built from the metadata's `ncbi_taxid` column) and `gtdb/labels.tsv`, parsing the input line by line
(label or synthetic taxId both work):

```sh
neoHGT gtdb map-report --db <库目录> --in <HGT清单> -o <输出TSV>
```

`--to` defaults to `ncbi` (v1 implements only this one target system). The output TSV gives the mapping result line by line,
and it prints the three counts rows / mapped / unmapped; lines that cannot be mapped are recorded as unmapped, not counted as failures.

## Hybrid source `--source hybrid` (measured)

`--source hybrid` grafts the top of the NCBI taxonomy tree onto the GTDB tree: GTDB's `d__Bacteria` and
`d__Archaea` are re-hung under NCBI's cellular root (taxId 131567), and, by top-level rank
(constructor argument `keep_ranks=('superkingdom',)`; the alias table also treats the `domain`/
`kingdom` written by the current NCBI dump as top-level), a number of eukaryote and virus top-level nodes are kept; the remaining NCBI nodes do not enter the tree;
how proteins shared across sources are handled is decided by `--hybrid-lca-policy` (`root`/`domain`/`refuse`, plan §5.7 item 3).
The NCBI tree's location is given by `--hybrid-ncbi-taxdump`; a missing file, or one without 131567, both produce a clear error
(suggesting to switch to `--source gtdb`).

> **These three settings were originally only "validated and printed"; the build flow never called the adjudication function (now wired in)**: when the extraction phase
> met duplicate protein IDs, it wrote `taxon.map` with "first taxId encountered wins", so `root` / `domain` / `refuse`
> made no difference whatsoever to a real build, and that line in the build report was merely restating the user's input (the IT-5 row in the table below, "the three lca_policy
> behaviors match expectations", verifies the merged tree and does not cover the build-time adjudication — read it together with this paragraph). The semantics now:
>
> - `root`: every candidate is kept, and the LCA of a cross-source shared protein lands on the cellular root (the shift recorded by RK7);
> - `domain`: walk up the merged tree to find that candidate's own **domain node** (`d__Bacteria` /
>   `d__Archaea` hanging under 131567, or an NCBI-side top-level node grafted in; both the `domain`/`superkingdom`
>   spellings of old and new dumps are recognized) and label accordingly; proteins with multiple candidates or crossing domains are recorded in `x_source_proteins`, and only when no domain node is found does it fall to root;
> - `refuse`: discard candidates whose taxId is outside the GTDB interval (`>= --gtdb-taxid-base`); if no candidates remain to keep as a result,
>   that `taxon.map` line is written as `na` rather than guessing one.
>
> The adjudication happens in a single **streaming replacement** pass after extraction (over `taxon.map.gz`'s source file and the BLAST taxidmap),
> with memory growing only with "the number of shared proteins". The hybrid section of the build report now gives the real counts:
> `proteins with more than one taxId candidate: N (reassigned M, left unassigned K)` —
> in a hybrid library of purely GTDB sequences N is usually 0, and this line then honestly writes 0, i.e. "this build had nothing to adjudicate".

Three pieces of evidence run through with real data (both the GTDB and NCBI dumps read-only):

| Gate | Method | Result |
| --- | --- | --- |
| IT-5 tree merge | `python tools/verify_hybrid.py --gtdb-dir $GTDB_DIR --ncbi-taxdump 缓存测试/input/ncbi_taxdump --out 缓存测试/reports/IT5_hybrid_r3.json` | PASS: pure GTDB tree 146,380 nodes, 146,404 after the merge (25 NCBI nodes grafted); single root; the two domains' parent node is 131567; the LCA of bacteria vs eukaryote lands on the cellular root (the semantic shift recorded by RK7); the three `lca_policy` behaviors match expectations |
| IT-5b full pipeline | 6-genome hybrid library → DIAMOND blastp → `neoHGT search` → `neoHGT analyze` (script `临时文件/verify_hybrid6.sh`) | **Byte-for-byte identical** to the `--source gtdb` library of the same 6 genomes: identical `hits.tsv`, identical search-table data rows, identical `scores.tsv` (42/42 hits carry staxid, HGT calls likewise 0) — grafting on top does not change the results for prokaryotic queries |
| IT-6e hybrid resume build | `tools/verify_resume.py --fresh 临时文件/hybrid6_step1 --resumed 缓存测试/db/hybrid6` (the second run with a brand-new empty `--tmpdir`) | PASS: 7 files byte-identical, 6 gz files identical after stripping comments, 5 differing only in timestamps, 0 mismatches, 0 missing |

The DIAMOND rank-vocabulary pit (found by measurement): the rank column of `nodes.dmp` is validated by DIAMOND, and the
**current** NCBI taxdump uses modern ranks — `cellular root`(131567), `acellular root`(10239
viruses), `realm`, `clade` — so `diamond makedb` exits straight away with `Error: Invalid taxonomic rank:
cellular root`. Therefore `gtdb.write_taxdump()` writes ranks not in the `gtdb.DIAMOND_RANKS` vocabulary as `no rank` (exactly what NCBI used for these nodes before changing the vocabulary in 2024). **Only the file is changed**:
the in-memory tree, `gtdb/lineages.tsv` and `gtdb/labels.tsv` keep the real ranks, and the tree read back from
`taxdump/nodes.dmp` also still works by taxId and prefixed labels. The flattened nodes are listed one by one in `meta.json`'s
`taxonomy.hybrid.ranks_flattened_for_diamond` and in the
"Hybrid graft" section of `gtdb/build_report.md` (on real data: `acellular root`×1, `cellular root`×1,
`clade`×1, `realm`×6). GTDB's own 7 ranks are all inside the vocabulary after the `domain→superkingdom` alias,
and `gtdb.illegal_dmp_ranks()` is always empty on a pure GTDB tree, so this fallback cannot mask a GTDB-side problem
(`test_gtdb_ranks_are_never_flattened` locks that; `test_modern_ncbi_ranks_survive_diamond`
locks that the flattening happens only in the file).

Another implementation constraint: the merge is done on a **copy**. The early implementation rewrote the tree cached by the wrapped GTDB provider directly,
so NCBI nodes leaked even into the "pure GTDB tree" (whichever was asked first decided what it looked like); now
`HybridProvider.load_taxonomy()` copies first and grafts second, and `test_merge_does_not_mutate_the_gtdb_tree`
locks this in.

The part still not working end to end: in a hybrid library **only GTDB proteins are searchable** (`--seq-source local`), because
downloading NCBI protein sequences has not been authorized; the grafted-in NCBI nodes currently serve only tree semantics. To make eukaryotic/viral donors
genuinely hittable (as the original paper's Asgard scenario requires), the NCBI reference protein library must additionally be obtained — see the D12 §10.4
open items.
## Measured resource budgets

| Tier | Genomes | Proteins | db.faa | diamond db | Notes |
| --- | --- | --- | --- | --- | --- |
| 6-genome smoke (IT-6 tier, `fixtures/lists/e2e6.IDlist`) | 6 | 14,373 | 6.2 MB | 614 MB (base `1e8`) | one-shot build 81 s; resumed run 99 s (skips extraction) |
| 6-genome smoke (IT-2 tier) | 6 | 27,667 | 12.3 MB | **620 MB** (base `1e8`) / **5.6 GB** (base `1e9`) | makedb peak RSS ≈ 3.8 GB; proves `db.dmnd` inflates with the taxId **numeric range** |
| `--rank-sample family --domain archaea` | **564** | **905,829** | **360 MB** | **985 MB** | wall clock **2 min 37 s**, peak RSS 6.5 GB, can run in CI |
| `--gtdb-taxa f__Enterobacteriaceae --reps` | **793** | **3,247,172** | **1.41 GiB** | **2.07 GiB** | wall clock **4 min 50 s**; `genomes.tsv` 793 rows = exact metadata count |
| `--rank-sample genus` (bacteria) | 23,112 | 60–80 million | 28–37 GB | 42–110 GB | extrapolated at the measured 466 B/protein |
| Full representative set | 113,104 | 333,384,550 | **≈145 GiB** (extrapolated) | ≈1.5–3× db.faa → 220–450 GB | extraction extrapolated at the 793-tier throughput ≈2–3 h + makedb 6–15 h; not executed |

**Size baseline**: measured `db.faa` = **466 B/protein** (three tiers: 455 / 417 / 466); the capacity precheck estimates from exactly this value (peak factor 3×). The plan's "406 GB / 333M proteins" corresponds to **all 373,000 per-genome files** (including the non-representative set), so extrapolating directly from it would overestimate the representative set by about 2.8×.

**Measured cross-volume read bandwidth**: single stream 60 MB/s, 8 concurrent streams 127 MB/s (NFS `/data2` → local `/data`). Below the plan's 200 MB/s threshold → **tiered builds are advised** (`archaea_family` first, then `genus`, finally the full set), and after building, free space with `--keep-faa compress`. A capacity precheck runs before building (free space must be ≥ 3× the estimated `db.faa`); if insufficient, the build aborts, unless `--force` is given explicitly.

**Two measured hard constraints of DIAMOND** (must be known):

1. `db.dmnd` size grows with the taxId **numeric range** (not the node count): the same 13 MB library gives 5.7 GB with `--gtdb-taxid-base 1e9` (measured in the S1 spike; re-measured with the 12.3 MB library at the IT-2 tier it gives 5.6 GB, see the table above), and 620 MB with `1e8`. **This implementation therefore defaults to `100000000`** (still about 30× higher than the largest NCBI taxId in production, ≈3.4e6), and rejects values below `10000000`.
2. The rank column of `nodes.dmp` must use the NCBI vocabulary (GTDB's `domain` gets written as `superkingdom`), and all 13 columns must be filled in; otherwise `diamond makedb` fails outright. This implementation already handles it automatically.

**BLAST route**: `makeblastdb` has a **50-character limit** on local ids, while GTDB sequence IDs are 151–270 characters → in GTDB mode `--compile blast|both` **exits with an error** (plan ADR-8); only an explicit `--force` makes it attempt the run. v1 only promises DIAMOND.

## Making `--taxmap` usable at full scale (P7, route B revived)

The plan originally **vetoed** route B (reading the 333 million "protein→taxId" pairs into a Python dict would need tens of GB of memory). It is now a **disk-anchored index**: `gtdb.TaxidIndex` keeps only one (key, byte-offset) anchor per 2048 rows (about 0.2% of the data volume) and caches the anchor table in `<map>.idx`; a query binary-searches to the anchor + a local linear scan, **no longer loading the mapping into memory**.

```sh
# 需要明文且按第一列字节序排序（即建库写出的顺序）
gunzip -k "$DB/taxon.map.gz"
neoHGT search -i q.faa -m precomp -s hits.tsv -t "$DB/taxdump" \
  --taxmap "$DB/taxon.map" -o "$OUT" --tmpdir "$TMP"
```

- Route B is only the fallback when the hit table **has no** `staxids` column and the sequence IDs cannot be resolved via route C either; the resolution order is still A → C → B → report `unresolved_taxids.txt`.
- Passing in a `.gz` is explicitly rejected with a prompt to decompress first (gzip cannot be randomly accessed by byte offset).
- Measured: 10,000 rows / spacing 512 → 20 anchors, 21-line index file; the second open reuses the cache directly (<0.001 s). Tests: 5 cases in `TestTaxmapIndex` + the workflow-level `test_route_b_after_route_c` (IDs without `staxids` that route C does not recognize are filled in by `--taxmap`).

## Memory fix on the NCBI route (P7, optional)

Upstream `database.py`'s `extract_genomes()` keeps a record in memory for **every unique protein** and directly prints the warning that it "may get killed by the OS"; the GTDB branch was streaming from the start, and this change feeds the same capability back into the NCBI route as an **explicit switch** (default `no`; the code and behavior of the default path are completely unchanged):

```sh
neoHGT database --source ncbi -c bacteria --stream-extraction yes -o "$DB" --tmpdir "$TMP"
```

- Deduplication now uses a sqlite index inside `--tmpdir`; `genome.map.gz` / `taxon.map.gz` / BLAST `-taxid_map` / DIAMOND `prot.accession2taxid` are generated by **external sorting + streaming grouping**, so peak memory no longer grows linearly with the total protein count; sorting fixes `LC_ALL=C` to keep the same byte order as upstream's `sorted()`.
- The cost is extra space in `--tmpdir` (about 30 bytes per protein occurrence).
- Equivalence is asserted item by item by `neoHGT/tests/test_database.py::TestStreamExtraction`: for the same input, `db.faa`, `genome.map.gz`, `taxon.map.gz`, and the per-genome protein/residue counts are **exactly identical** to the upstream implementation; measured peak RSS is in `benchmarks/stream_extraction_slope.txt` (two scales, 200k / 1M proteins):

  | Implementation | 200k proteins | 1M proteins | Marginal cost | Extrapolated to 333.4M proteins | Extraction time (1M) | tmpdir(1M) |
  | --- | --- | --- | --- | --- | --- | --- |
  | Upstream | 365 MB | 1,071 MB | ≈882 MB per million proteins | **≈287 GB** | 20.3 s | 0 |
  | Streaming | 193 MB | 304 MB | ≈138 MB per million proteins | **≈45 GB** | 45.9 s | 104 MiB |

  That is, memory drops **about 6.4×**, at the cost of about 2.3× the extraction time and about 100 MiB of temp disk per million proteins.
- Shared proteins that need an LCA are computed **group by group** in the streaming implementation (only a few taxIds per group), so it both saves memory and stays consistent with upstream's `max(tids) / find_lca()` rule.

## KB-scale test fixtures and offline end-to-end (§10.1 / IT-4)

`neoHGT/tests/data/GTDB/` is a KB-scale fixture (~119 KiB) that **is by itself a searchable library directory**, generated by `tools/make_fixtures.py` from read-only local r220 data:

```
GTDB/
├── protein/<6 个基因组的 .faa 摘录>     # 每个只留前 5 条蛋白（复制输入，不是软链）
├── metadata_bac120_sample.tsv / metadata_ar53_sample.tsv / taxonomy_sample.tsv
├── layout.txt / fixtures.json / meta.json
├── taxdump/{nodes.dmp,names.dmp}        # 转换器产物（32 个节点）
├── gtdb/{taxid.map.gz,labels.tsv,lineages.tsv,stem2taxid.tsv.gz,
│        genome2taxid.tsv.gz,ncbi_map.tsv.gz}
├── db.faa                                # 摘录合并（30 条蛋白）
└── search/sample.tsv                     # 7 列含 staxids 的命中样例（168 行）
```

> This inventory is not "write it once and forget it": `tools/verify_fixture_provenance.py` reopens the read-only sources,
> checking entry by entry that every fixture sequence really exists in its source with matching residues, that the source sha256 and byte count are unchanged,
> that `protein/` contains no unregistered orphan files, that the actual record count of `db.faa` equals the value published by the inventory, and that the genome count is ≤ 10.
> The report is at `缓存测试/reports/R35_fixture_provenance.json` (workspace cache/temp directory, kept under its Chinese name) (measured: all 6 sources read, 0 problems).

Rebuild (the source data is read-only; at most 10 genomes copied per run, each copy recorded with its sha256 into `缓存测试/fixtures/manifest.tsv`; run `--dry-run` first to preview the copy list):

```bash
. tools/env.sh
python tools/make_fixtures.py --dry-run
python tools/make_fixtures.py --limit 6
```

The hit sample is **constructed by lineage distance** and does not depend on the local DIAMOND version, so it reproduces byte for byte across environments:

- Most proteins have themselves as their best hit, followed by descending scores in the order "same family → same domain → cross-domain", exercising the self/close grouping;
- The last 2 proteins of each genome **have no self hit**; their best hits come from the most lineage-distant genome, i.e. an artificial HGT.

`neoHGT/tests/test_gtdb_pipeline.py::TestFixtureEndToEnd` runs the whole chain offline with **the real CLI parser** (`parse_args` of `scripts/neoHGT`), equivalent to:

```bash
neoHGT search  -i GTDB/protein/RS_GCF_000026225.1_…_g_Escherichia.faa \
               -o <out> -m precomp -s GTDB/search/sample.tsv \
               -t GTDB/taxdump --source gtdb --threads 1
neoHGT analyze -i <out> -o <out2> -t GTDB/taxdump --source gtdb \
               --input-tax <genome stem>:<species taxId> \
               --self-rank genus --donor-rank genus --donor-label prefixed \
               --close-size 1 --silhouette 0
```

What is asserted is not "it ran", but that **the 2 planted hits are recovered**, that the donor is `g__Nanoarchaeum`, that `scores.tsv` contains the `distal` column, that `scatter.pdf` and the various KDE/hist PDFs are all present, that every hit has a numeric taxId (DoD-3), and that `urlopen`/`run_command` are not called even once during the whole run (DoD-4). `--silhouette 0` is because a 5-protein sample cannot support interpreting the silhouette coefficient; it should not be turned off on real data. In these two tests the fixture is only a **read-only source**: the tests first copy it into a temp directory and then run, and at the end a whole-tree sha256 comparison proves that not one byte of the fixture committed to the repository was written.

## Zero-regression gate (RT-1, required reading for contributors)

The GTDB adaptation must **not touch any existing behavior of the NCBI path**, so besides the unit tests there is a mechanically checkable gate:

```sh
python tools/verify_regression.py \
    --baseline /path/to/pristine/neoHGT-main \
    --code     "$NEOHGT_CODE" \
    --out      "$NEOHGT_CACHE/regression/RT1_failure_set_diff.txt"
```

It does two things:

1. **Failure-set comparison entry by entry** (not comparing counts): run `neoHGT/tests` in both trees and diff by test id; `newly_failing` must be empty. Upstream already has 13 failures of its own (plot-format assertions + external binaries missing from PATH); that list is frozen in `ci/rt1_baseline_failures.txt`, and the CI step compares against it directly.
2. **Deleted-line audit**: `neoHGT/{database,search,analyze,util}.py` and `config.yml` must have `deleted_lines = 0` relative to the baseline (pure additions + method rebinding); the only permitted change is the 1 line in `scripts/neoHGT` that adds `gtdb` to `module_names`.

The companion `tools/verify_ncbi_fixtures.py` covers the other half, "byte-for-byte comparison of fixture outputs": with the same interpreter and the same random seed, run the NCBI analysis path once in each of the two trees (`neoHGT/tests/data/Ecoli/search` → silverman, then `from_scores` → grid), compare sha256 file by file, and first run each tree twice to prove in-tree reproducibility. matplotlib writes `/CreationDate` and `/ModDate` into PDFs, so the raw bytes of two runs of the same tree differ anyway; PDFs are therefore compared with those two fields removed, everything else byte for byte. Note: **the `Ecoli/analyze/scores.tsv` already committed in the repository cannot be reproduced even by the original upstream code under the current scipy/scikit-learn** (the two trees produce identical output, but different from the fixture, and upstream's `test___call__` fails the same way in the baseline for that reason), so this gate compares "are the two trees consistent", not "does it equal that old fixture".

There is also a "name matches reality" gate: `tools/verify_plan_traceability.py` extracts all `UT-*/IT-*/AT-*/RT-*` verification entry numbers from the plan document and requires every number to be either cited by name in some test file or mapped to an evidence artifact that **really exists**; otherwise the exit code is 1. The product is the cross-reference table `缓存测试/reports/plan_traceability.md`. This gate works: UT-0.3, UT-2.1/UT-2.x, UT-3.3, UT-3.4, UT-4.x, and UT-5.x were all once "written in the plan with no corresponding test in reality", found by its scan and then filled in (AT-4 was found earlier by manually reading the plan and was included in the same table); change any number in the plan to a nonexistent value and you will see it report FAIL.

If pytest produced no summary line (missing interpreter, collection errors, etc.), the tool **refuses to compare** instead of reporting a pass — empty sets being "equal" to each other is exactly the easiest fake green to disguise. Same in CI: assert the baseline file is non-empty, assert the summary line exists, assert the parsed `FAILED` lines are non-empty; if any of these is not satisfied, the check fails.

## Known limitations

### Measured cautions for the grouping parameters (three samples)

`缓存测试/reports/S10_4_3_sensitivity.md` is the result of re-scanning `close_size` / `input_cov` / `bandwidth` / `noise` on three samples: halo (archaea), bathy (an archaeon from a different phylum), entero (bacteria, a single-family library). Three cautions that would change the conclusions:

1. **`--close-size` is a lower bound, not a quota.** The default picks the tier by domain (bacteria `grouping.profiles.gtdb` = 200, archaea `gtdb_archaea` = 50); when the tree cannot yield that many taxa it is truncated (bathy has only 20 close taxa at its default tier), so "used the default" does not mean "the close group has 200/50". The `Close group has N taxa.` line in the log is the actually effective value; cross-sample comparisons must look at it.
2. **The direction of `--noise` is not guaranteed to be monotonic.** The config comment says "larger is stricter", but two samples measured in opposite directions (halo 25→75: 374→325 stricter; bathy 25→75: 285→421 looser). What it changes is the excluded interval endpoint in the valley-to-peak distance, and on bimodal or long-tailed distributions the direction flips. To tune this parameter you must look at your own sample's curve; do not quote the comment.
3. **`--input-cov` affects the input classification, not the number of predictions.** When the hit classification is clean (best hits 100% within the same taxon) the whole 50–95 range is flat; only when the coverage sits near the threshold (measured 89.47% on the bacterial sample) does cov=95 demote the input classification from species to genus. It is therefore a switch for "whether the sample's own classification is trustworthy", and `--self-rank` decides whether this demotion affects the results.
4. **`donor=0` does not mean "no donor".** It means the LCA of the distal hit is coarser than the rank requested via `--donor-rank`: measured on E. coli against the 4,896-genome cross-family library, 10.4% are annotatable when genus level is requested, 13.1% at family level, 100% at domain level, and among those `0` LCA nodes themselves, 591 fall at superkingdom, 65 at class, 21 at phylum. `--distal-top` is completely flat across 10–100 on this dataset and is not a usable knob; for higher resolution, loosen `--donor-rank`.
5. **`--from-scores` and a full run must give the same result, but once they did not.** GTDB's profile defaults and the `domain → superkingdom` rank translation were originally done only inside `read_input()`, which `--from-scores` skips, so `--donor-rank domain` could not label a single prediction on the reuse path. It has now been extracted into `apply_gtdb_grouping()` shared by both routes (measured after the fix the two routes agree exactly: 1235 predictions, 100% labeled), and is pinned by `neoHGT/tests/test_gtdb_pipeline.py::TestFixtureEndToEnd`.
6. **The same library with the same parameter does not give consistent `close_size` curves across four bacterial samples (measured).** On the 4,896-genome "one representative per family" library: *Escherichia* is 771 constant across 10/50/100/200, Alicyclobacillaceae is 753 constant, f__ATN3 is 714 constant, while **Arcobacteraceae gives 436 at 10/50 and drops straight to zero at 100/200** — its self group has only 2 taxa (that genus has only one other representative in this library), and at `close_size≥100` the inferred close group jumps in one step to the entire phylum (13,208 taxa), swallowing all the hits; the distal group is empty, and the log gives `WARNING: No hit is assigned to distal group. Cannot predict HGTs.`. So "the bacterial tier 200 is usable" is a **per-sample** property, not a property of the library; data taken from `缓存测试/reports/R26_bac_close_size_curves.json` (prediction-only tier, DIAMOND not re-run).

Also confirmed along the way: `--from-scores` reusing the already-written `scores.tsv` is **equal to a full run on the prediction set** (both bathy tiers give 724→367), so scans over threshold-like parameters (`bandwidth`/`noise`/`fixed`) can safely take the reuse path, saving a full scoring round each time. But "equal" only holds for the set of proteins called and the ranking; the floating-point columns can differ in the last digit because `scores.tsv` is persisted with 6 significant digits (see item 5 of this subsection), so do not assert byte-for-byte equality when comparing.

- No eukaryotes/viruses: GTDB has only `d__Bacteria` and `d__Archaea` (measured).
- hybrid mode needs the NCBI taxdump: the local GTDB copy does not contain that file; this time `input/taxdump.tar.gz` was downloaded
  only after separate authorization from the consumer side. **The hybrid tree has passed joint testing with a real NCBI taxdump + a real GTDB r220**:
  IT-5 (tree merging and LCA semantics), IT-5b (6-genome library build → DIAMOND search → full analyze chain),
  IT-6e (hybrid resumable build). Still unauthorized is hybrid's **NCBI protein sequence** channel
  (`--seq-source ncbi`), so what is searchable in the hybrid library remains GTDB proteins only — the grafted-in NCBI
  taxonomy nodes currently serve only the tree's semantics and contribute no hits; to genuinely deliver the original paper's Asgard scenario in which
  "eukaryotic/viral donors can be hit", the NCBI reference protein library must be obtained separately (see open item D12).
- The remote sequence channels (`--seq-source ncbi|gtdb-tar|ena`) only export a fetch plan,
  `gtdb/fetch_plan.tsv`, for manual execution — this version has no download executor and will not initiate any networked download;
  `--allow-network yes` is only a prerequisite confirmation for using these channels, not a networking switch.
- Lineage embedded in file names stops at genus level; at species level, metadata/taxonomy is always authoritative.
- The local copy was renamed and had its sequence IDs rewritten by third-party scripts, and bears traces of a merge with r207; there is no official MD5SUM,
  so this tool builds its own input sha256 fingerprints and enforces a self-consistency check.
- `--keep-faa drop` currently **only warns and does not delete** (deletion requires backing up first, in keeping with the project's security rules).
- GTDB taxIds are synthetic IDs and must not be mixed with NCBI taxIds on the same historical result tree; cross-release ID stability depends on
  `gtdb/taxid.map.gz` — verified bidirectionally in measurements between local r214/r220 (0 changes / 0 conflicts /
  0 reuse; see the end of the "building the taxdump separately" section and `缓存测试/reports/R24_lineage_only_and_release_upgrade_SOP.json`).
- The official release's **flat layout** (lineage table placed directly in the data root) is now recognized, and `gtdb taxdump` can work from
  the lineage table alone; but a full library build still needs the metadata table, and the two other local copies (r214, r226) carry no metadata,
  so "one run straight through from the official dataset" has still never been measured. The column mapping itself is calibrated against the official schema: the 11
  `REQUIRED_COLS` columns hit 11/11 in the real `bac120_metadata_r220.tsv` header.
- The lineage strings in the r226 local copy were glued together with `_` and truncated to genus level by a third party — **unqualified input**: the deep check will
  reject it before any write (exit code 1); the official original `*_taxonomy_r226.tsv` must be used instead.
