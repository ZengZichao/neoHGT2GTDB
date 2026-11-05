> English (primary) | Chinese mirror: [测试用例对照.md](../测试用例对照.md) — keep the two in lockstep. Chinese mirror linked above.

# Plan §10 verification-matrix check (UT / IT / scientific validity / RT)

The unit-test criteria of §10.2, the integration-test scenarios of §10.3, the four pieces of evidence for "is it actually more trustworthy" in §10.4,
and the regression and compatibility items of §10.5 are the **verification plan itself** laid down by the plan. This file states, entry by entry, which test
or which piece of evidence covers each ID, and it accepts the machine check of `tools/verify_plan_traceability.py` (the same set of
rules as the other nine registers):

1. Every ID in the four tables must appear exactly once, and the register must not contain IDs absent from the plan;
2. The status may only take `已覆盖` / `部分覆盖` / `未覆盖`, and the latter two must state clearly which half is missing;
3. Referenced backtick paths (optionally with `:行号`) must really exist, and the line number must not exceed the file length; references to test files are always written in the `路径::符号` form, and the gate will confirm the file really contains that `def` / `class` — the measured lesson of Round 20: inserting one test into a test file silently misaligns every line-number reference, and such a misalignment happens to be invisible to the "path exists" check.

The motivation for writing this table was an inventory taken in Round 15: the derived fixtures §10.1 requires did not exist at the time, and the GTDB
branches of `search`/`analyze` had never had any **offline** test that actually ran. Spreading the verification matrix out into a checkable table exists
precisely so that the silent gap of "written in the plan but actually absent from the code" never appears again.

In the plan, §10.4 is an ordered list rather than a table, so its IDs are written here as `10.4-<n>` (in one-to-one correspondence
with "item n" in §10.4); without stable IDs there can be no registration, and without registration no gating.


## §10.2 Unit / property tests

| ID | Status | Location | Notes |
| --- | --- | --- | --- |
| UT-0.1 | Covered | `neoHGT/tests/test_gtdb_pipeline.py::test_renamed_layout_recognized`, `neoHGT/tests/test_gtdb_pipeline.py::test_official_flat_layout_is_discovered` | The local renamed layout is recognized as gtdb-r220-renamed; official-plain / official-gz **are recognized both on synthetic trees and on the real official layout** (an unpacked release puts `bac120_taxonomy_rNNN.tsv(.gz)` directly in the data root rather than inside `tree_bac120/` — in Round 24, using the local r214 `parent_directory`, it was discovered that the original implementation only recognized wrapper directories and could not see root-level tables; it has now been fixed as "wrapper directories first, root level as the fallback", with the precedence locked by `neoHGT/tests/test_gtdb_pipeline.py::test_wrapper_directories_take_precedence_over_root_files`); when undecidable it reports a clear error and lists candidates (`neoHGT/tests/test_gtdb_pipeline.py::test_undecidable_layout_lists_candidates`), and forcing a wrong layout is likewise rejected; `--gtdb-taxonomy/--gtdb-metadata` now take effect **before discovery fails** (`neoHGT/tests/test_gtdb_pipeline.py::test_table_overrides_are_applied_before_discovery_fails`). |
| UT-0.2 | Covered | `neoHGT/tests/test_gtdb_core.py::test_cluster_letters_and_hyphens` | p_Bacillota_I, c_Bacilli_A, f_Enterobacteriaceae_A, o_ML615J-28, f_DSM-6578, g_336-3, g_DUMJ01 all parse correctly; it also asserts that when the filename carries no species-rank segment **no species may be invented out of thin air** (`neoHGT/tests/test_gtdb_core.py::test_no_species_is_invented`). |
| UT-0.3 | Covered | `neoHGT/tests/test_gtdb_pipeline.py::test_selfcheck_fails_loudly_on_one_inconsistent_file_count` | The check items of §3.6 are implemented one by one; artificially change one domain's faa file count by 1: with strict=False it pinpoints exactly which item failed, and with strict=True it raises SelfCheckError and refuses to build the database. |
| UT-1.1 | Covered | `neoHGT/tests/test_gtdb_core.py::test_node_count_matches_unique_prefix_paths` | Node count = number of unique lineage prefixes, every parent exists, the root points to itself, and the rank distribution matching the input are all asserted; the claim that the r220 full set should have 146,379 nodes belongs to the full-set check, with evidence in IT-1. |
| UT-1.2 | Covered | `neoHGT/tests/test_gtdb_core.py::test_manifest_appends_without_renumbering` | The bijection and the ID range (≥ taxid_base and < 2^31) live in `neoHGT/tests/test_gtdb_core.py::test_ids_are_in_reserved_int32_range`; appending new lineages without changing existing IDs is asserted by this item, and in chained rebuilds retired IDs continue to be kept in reserved_ids.txt (see RK6). Round 24 added the real scenario of "a release version swap" as the test `neoHGT/tests/test_gtdb_core.py::test_release_swap_never_reissues_or_reuses_a_taxid`: performing "delete one lineage + rename one + add two" at once, it asserts that shared-node IDs are unchanged, **new nodes must not receive any ID known to the old table** (including the IDs of just-disappeared nodes), and disappeared nodes are faithfully registered as superseded. Note the criterion is not "new IDs must be greater than the old maximum" — allocation is chunked by rank (`slot = id_base + rank_ordinal * 1e6`), and new nodes legitimately fall inside their own segment; the real-data check is in `缓存测试/reports/R24_lineage_only_and_release_upgrade_SOP.json` (r220 mapping onto r214 lineages: 89,172 shared nodes with 0 changes, 24,145 new nodes with 0 conflicts, 57,207 disappeared nodes all entering `delnodes.txt`; the reverse direction likewise). |
| UT-1.3 | Covered | `neoHGT/tests/test_gtdb_core.py::test_round_trip_through_files` | After the conversion artifacts are written to disk they are read back with util.read_taxdump and equal the in-memory structure field by field; the upstream reading function itself also has a regression test (`neoHGT/tests/test_util.py::test_read_taxdump`). |
| UT-1.4 | Covered | `neoHGT/tests/test_util_gtdb.py::test_capital_guards` | All three boundaries are locked down: s__Escherichia coli with the prefix stripped is true, the original prefixed string is false, and d_Archaea is false (i.e., "the prefix must be stripped first"). The actual exemption behavior on the search side is additionally covered by `neoHGT/tests/test_gtdb_pipeline.py::test_gtdb_profile_defaults_the_filters_off`. |
| UT-1.5 | Covered | `neoHGT/tests/test_util_gtdb.py::test_latin_guards` | s__Bacillus_A cereus is judged a Latin name once the cluster suffix is stripped; the cluster suffix is also stripped correctly at genus and phylum rank (p_Bacillota_I). |
| UT-1.6 | Covered | `neoHGT/tests/test_util_gtdb.py::test_gtdb_ranks` | GTDB returns 7 ranks and contains no superkingdom; the rank-prefix table in `neoHGT/tests/test_util_gtdb.py::test_rank_prefix`; alias normalization in `neoHGT/tests/test_util_gtdb.py::test_rank_normalize_aliases`, and resolve_taxon_ref("superkingdom") landing on domain is asserted by the tests starting from `neoHGT/tests/test_util_gtdb.py::test_numeric`. |
| UT-1.7 | Covered | `neoHGT/tests/test_gtdb_core.py::test_lca_matches_lineage_prefix` | On the converted GTDB tree, find_lca agrees with the "longest common prefix of the lineage strings"; there is also a test that the LCA in the cross-source (hybrid) case sits above both domains, `neoHGT/tests/test_gtdb_pipeline.py::test_cross_source_lca_sits_above_domains`. |
| UT-2.x | Covered | `neoHGT/tests/test_gtdb_pipeline.py::test_protein_total_survives_column_normalization`, `neoHGT/gtdb.py::reconcile_protein_counts` | Parsing by column name (shuffling the column order leaves the result unchanged), an error naming which column is missing when a column is missing (`neoHGT/tests/test_gtdb_pipeline.py::test_missing_column_names_the_column`), genome/accession uniqueness, **row counts must not be lost in the column mapping**, and "the total protein count is a comparable quantity" are all asserted: reading the raw TSV independently and summing `protein_count` must equal the sum of `proteins` after `load_metadata()` normalization, and silent NaN must not appear. Round 18 added another layer of reconciliation: `meta.json` now records both `proteins_planned` (the number of proteins the genomes selected this time should have) and `proteins_extracted_vs_planned` (how many extraction actually obtained, and by how much they differ), so accidents like "the sequences of an entire genome missing" show up in the manifest instead of surfacing weeks later as unexplainable HGT counts. (On a real 4,896-genome build this item caught and fixed two reconciliation defects: `write_meta_json()` originally took `provider.df` — the table of the **entire release** — as the baseline, and so reported a library that had extracted all 12,958,575 proteins as "1,944,720,987 short"; when changing the baseline it also unpacked the helper's return value by position, writing the **number of extracted entries** into the "should-have count" field, the same manifest claiming to be fully extracted and 1.94 billion entries short at the same time. Reconciliation now goes uniformly through `neoHGT/gtdb.py::proteins_reconciliation` (returns a dict, fields bound to values), and the baseline is passed in explicitly by the build path as `selection=df`; three tests lock the baseline (`neoHGT/tests/test_gtdb_pipeline.py::test_protein_reconciliation_uses_the_selection_not_the_release`), fields not swapping places (`neoHGT/tests/test_gtdb_pipeline.py::test_proteins_reconciliation_labels_both_sides_correctly`), and the publishing chain being complete (`neoHGT/tests/test_gtdb_pipeline.py::test_meta_publishes_the_reconciliation`)) The number 333,384,550 for the r220 full set is itself the value of this field on the full tier (see `缓存测试/reports/P0_spike_S3_cost.md` and the `counts.representatives_in_dataset.proteins` of `缓存测试/db/r220_ar_family/meta.json`). Round 27 added the "option → parse → take effect" chain: the 15 options of plan §7 previously had **no test handing them to the real parser** (the gate could only prove "the flag exists in the CLI and is read somewhere", blind to negative failures such as a wrong dest name, a wrong type, or the value never reaching the selector); now `neoHGT/tests/test_gtdb_pipeline.py::TestPlanOptionWiring` asserts entry by entry that the parse result lands on the expected dest, and asserts the **filtering consequence** for `--reps/--genome-list/--min-completeness/--max-contamination/--mimag/--type-material/--gtdb-taxa`; this rule also entered the gate (`tools/verify_plan_traceability.py::flags_without_tests` — a §7 option with neither a test nor a registered deviation is judged a problem; negative control in `neoHGT/tests/test_gtdb_pipeline.py::test_every_plan_option_is_handed_to_the_parser`). Along the way one real defect was fixed: for `--genome-list`, the bare `GCA_000000003` spelling the documentation promises to `_expand_acc` actually matched no genome (expanded to `GB_GCA_000000003` it equals neither the versioned `accession` nor the prefix-less `accnov`); matching is now done after normalizing both sides to "the prefixed form without version". Round 29 went one level deeper: **"handed to the parser" does not mean "the branch reads the right column"**. Hitting the same table with the three quality definitions one at a time revealed that `--type-material ncbi` and `gtdb` shared one and the same `gtdb_type_species_of_genus` branch, the upstream `--typemater` was consumed by nobody in GTDB mode, `--complete` was discarded after parsing, and `--mimag`/`--type-material` silently let things through `if col in sel`-style when a column was missing. Now `neoHGT/tests/test_gtdb_pipeline.py::TestQualitySelectionOptions` (12 tests) locks down: each definition reads its own column and `na`/empty strings do not count as designation, `--typemater` is an alias of `ncbi`, unknown definition values and missing columns both ValueError, the `--complete` error must name the replacement **and give the reason**, `quality_classes` counts correctly, and the filter label must carry the real counts (`GTDB type species (3 -> 1 genomes)`); the counts on the real full r220 table are in `缓存测试/reports/R29_quality_class_counts.json`, and the build report gained a `quality classes in the dataset:` line (its existence asserted by `neoHGT/tests/test_gtdb_pipeline.py::TestFixtureDatabaseBuild`). The same round added one more layer: the build path originally assembled the provider opts inline, so `--complete`/`--typemater` parsed successfully yet were not in this table (equivalently, out of the command line's reach); opts are now generated by `neoHGT/gtdb.py::selection_opts`, and test_parsed_options_reach_the_build_opts of `neoHGT/tests/test_gtdb_pipeline.py::TestPlanOptionWiring` asserts that every flag's parse result really enters the opts and that database.py installs exactly this table, and requires `--complete` to error starting from the CLI too. Round 30 dug this same layer deeper: **"the option is validated and printed into the report" does not mean "the build flow ever called it"**. New gate rule `tools/verify_plan_traceability.py::choice_values_without_tests` (the match scope is a single `test_` function body; a table-driven `for tier in (…)` counts as coverage, "the word mentioned elsewhere in the file" does not; the database stage is judged, while for the search/analyze stages this round only reports counts and registers them in the "option value deviations" section of `docs/交付物清单.md`); plus three classes of tests: `neoHGT/tests/test_gtdb_pipeline.py::TestSharedProteinResolution` (candidates for duplicate proteins must be handed to `--hybrid-lca-policy` for adjudication, `refuse` lands on `na` instead of guessing, behavior unchanged when there is no resolver), `neoHGT/tests/test_gtdb_pipeline.py::TestBuildArtifactChoices` (`--keep-faa keep/compress/drop`, `--input-fingerprints yes`, `--stream-extraction yes`, `--local-layout gtdb-r220-renamed`, `--gtdb-mirror auto`, `--seq-source gtdb-tar`, `--expand-members yes/no`, each running the real value one by one, where compress must leave `db.faa.gz` behind and no plaintext anymore, and drop must print the reason for not deleting), and `TestDatabaseChoiceValues` (`--domain` three values, `--mimag` four values, `--rank-sample` six tiers verifying the group counts one by one; yes/no-type options must really become booleans — the string `'no'` is itself truthy, an entire class of silent negative failures). The test file also had an **empty duplicate class** `TestFixtureDatabaseBuild(FakeInput)` deleted (an empty shell left by an earlier editing accident, waved through by the name check and deceiving readers). Round 32 added two new machines, each with its own negative control: `tools/verify_production_reachability.py` asks whether **the command line can reach it** (a different thing from "did the tests run"), judging each of the 196 functions this adaptation added, and the 5 unreachable ones must have their reasons stated in the "reachability deviations" of `docs/交付物清单.md`; test `neoHGT/tests/test_gtdb_pipeline.py::TestProductionReachability` (including a negative control: replacing the build-side `resolve=provider.resolve_shared_protein` with a dummy lambda must make that function unreachable again); `tools/verify_plan_traceability.py::self_discovery_check` checks whether the entry count and round range hand-written in the heading of the "self-discovery" list equal the list itself, and test `neoHGT/tests/test_gtdb_pipeline.py::TestSelfDiscoveryRegister` tests each of the four drifts once (count one behind, duplicated numbering, the newest entry missing its round, the range stuck at the previous round). The same round also changed the criterion of `neoHGT/tests/test_gtdb_pipeline.py::TestSearchGtdbModeDetection` from **string matching to AST matching**: the original criterion searched for the literal `"getattr(self, 'gtdb_root', None)"`, and the fifth decision point happened to be folded across two lines, so it reported green while genuinely missing one spot; the comparison after injecting the folded spelling back into a scratch copy is in `缓存测试/reports/R32_mode_negative_control.json` (the AST rule catches it, the literal rule does not). Round 33 added two: `neoHGT/tests/test_gtdb_pipeline.py::TestCliSurfaceAndDocExamples` (4 tests) pins that `tools/verify_plan_traceability.py::cli_flags` now reads the CLI surface from the **option table** (the old implementation scraped `gtdb`'s source with an `add_argument` regex, blind to the five options whose names occupy a line of their own in the table), and pins the new rule `tools/verify_plan_traceability.py::doc_command_flags`, which judges pairwise "whether that tier of an option in an example really exists", with exemptions effective per line only; `缓存测试/reports/R33_doc_flag_controls.json` is its negative control. Round 34 wired the value rule to the fourth CLI surface: the option table of `neoHGT gtdb` itself had never been asked before, and joining it immediately exposed `--action taxdump`/`--action map-report`, two values no test had ever handed over; `neoHGT/tests/test_gtdb_pipeline.py::TestGtdbActionChoices` (7 tests) makes all three values pass through the shipped parser, land in their respective functions, checks where the parameters land (`--gtdb-taxid-base` lands in id base rather than the manifest-parameter tier, `--no-selfcheck` arrives as "do not check"), and asserts that the three values land in **three different** destinations; the negative control uses a scratch tree (one new option with two values and no tests → two problems must be reported), and simultaneously pins "judged scope = CLI_STAGES, the backlog must be empty". Total value count 104 → 107. Round 40 turned the lens on configuration and tiers: the new rule `tools/verify_plan_traceability.py::config_key_problems`
  requires that "config.yml keys introduced by this adaptation must be read", with `neoHGT/tests/test_gtdb_pipeline.py::TestConfigKnobs` (2 tests)
  pinning both directions; `neoHGT/tests/test_gtdb_pipeline.py::TestGtdbGroupingProfile` (5 tests → 12 tests) pins that the grouping tiers really land
  (archaea 50, bacteria 200, explicit options take precedence, the generic `grouping.closesize=10` must not be taken for a tier,
  `set_parameters_gtdb` records the values actually passed on the command line, because `--from-scores` does not go through `read_input()`);
  Round 41 added to this class the layer of "both routes must give the same answer" (7 tests): the same batch of taxa walked via
  the full route (`data` hits the table) and the reuse route (`df['match']`) must land in the same tier; the donor column's
  `int64`/`float64` (a single null turns it to float) must both be readable; the reuse route must state plainly that grouping tiers are not in effect this run;
  when the domain cannot be determined it prints fallback instead of staying silent; a half-archaea half-bacteria tie must not luckily be assigned to a tier,
  plus a negative control (removing the scores table makes it "undeterminable" again); four mutants each had half of the fix deleted,
  and must be caught by the tests written for that half, with every mutant first asserting that `neoHGT.analyze.__file__` lands inside the copy
  (the first version did not assert this; `import` walked out of the scratch tree and loaded the shipped code, and all four mutants "passed"),
  see `缓存测试/reports/R41_falsification.json` and `缓存测试/reports/R41_domain_probe.json`.
  Round 42 moved the lens from "whether the tiers are right" to "whether the flag counts": `--source` had previously been read only once, inside `read_input()`,
  and `--from-scores` skips that step, so the same `--source ncbi` was the NCBI face in a full run
  yet entered the GTDB branch in a reuse run (`临时文件/r49_src_ncbi_*.log` is the pre-fix measurement,
  `r50_src_ncbi_reuse.log` the post-fix one: not a single GTDB line appears anymore, and `--source gtdb` meeting a library without the marker
  refuses with exit code 1, naming the marker and the way out, and writes no artifact at all, `r50_src_gtdb_nomarker.log`).
  The new class `neoHGT/tests/test_gtdb_pipeline.py::TestSourceHonouredOnBothRoutes` (7 tests) locks down
  the four values of `neoHGT/analyze.py::detect_gtdb_root` (ncbi does not go looking for the marker, auto gives the same answer on both routes,
  gtdb errors, hybrid prints WARNING), plus RT-1's boundary: without `--source`, the reuse route's NCBI face
  is word-for-word identical to upstream (even upstream's `--from-scores` pitfall is preserved as-is), and only an explicit pass adds back the skipped step.
  `neoHGT/tests/test_gtdb_pipeline.py::TestReuseRouteInertOptions` (7 tests) locks the three directions of rule 17:
  delete a real declaration, add a fake declaration, leave a stale exemption — each mutates on a **copy tree** and then calls the rule body itself,
  without touching shipped code; the comparison of the five mutants is in `缓存测试/reports/R42_falsification.json`
  (outside `neoHGT/` you must copy `tools/` along too — a package-only copy makes the tests that read the gate tool fail for the wrong reason).
  Round 43 added rule 18: the **values** of documentation examples must be accepted by that tier (`--source` has no `auto` in the database tier;
  writing `neoHGT database --source auto` is rejected by argparse).
  `neoHGT/tests/test_gtdb_pipeline.py::TestDocExampleValues` (7 tests) pins: the shipped documents have 0 problems,
  planting `--source auto` in a copy must be named; `probe|nonsense` names only nonsense without implicating probe;
  placeholders are not judged; **prose mentions outside fences are not judged** (the first version scanned whole lines and "caught" six false violations in a clean document);
  deliberately left unpaired fences must be reported as "not seen"; an empty document tree refuses to pass vacuously;
  measured in `缓存测试/reports/R43_doc_values_probe.json` and in `plan_traceability.json`'s
  `doc_command_values` section (97 examples with values / 66 literal judgments / 31 placeholders).
  Round 44 covered "a log gets written and nobody reads it": the new class
  `neoHGT/tests/test_gtdb_pipeline.py::TestExternalWriteGate` (5 tests) pins that IT-7 now reads back
  the `EXTERNAL` lines of `commands.log` — paths outside the workspace must be named; a missing log or not a single command both count
  as "nothing was looked at" rather than a pass; new files in this toolchain's namespace are judged failures while someone else's application is only reported as unattributable,
  and the real workspace must be clean under the current window; `TestDocExampleValues` grew from 7 tests to 10, the three new ones bringing the plan file
  under the same value reading (planting `--source automatic` in a copy must be named, `local` is not implicated,
  and with `--plan` missing it truthfully reports not checked).
  Round 45 fixed "the same value giving two answers on two routes": `bandwidth: 1` (a YAML integer, the maximum the help text allows)
  crashed on the configuration route at `perform_kde()`'s `isinstance(bw, float)`, while `--bandwidth 1`
  (the parser has no `type=`, so it arrives as the string ``'1'`) ran to completion — the end-to-end comparison is in
  `缓存测试/reports/R45_bandwidth_routes.txt` (before the fix cfg crashed / cli 222 entries; after the fix both routes give 222 entries).
  The new class `neoHGT/tests/test_gtdb_pipeline.py::TestBandwidthValueTypes` (4 tests) pins: both routes the same answer,
  `auto|grid|silverman` passed through as-is, out-of-range values must report the **source and the usable range** (the upstream line
  `Invalid bandwidth: 1.` does not say what counts as legal — exactly why this bug lived so long), and that `True`
  does not count as a number (in python bool is an int). `TestExternalWriteGate` also grew to 6 tests: one was added saying
  "non-judged prefixes only report the conclusion and no longer walk the tree just to print" (Round 44's IT-7 walked all of `/tmp`'s
  200,000 entries every time, with the number even fluctuating with other people's processes).
  Round 46 was the **mirror image** of the same shape: `--fetch-queries`/`--fetch-retries` in the `arguments` table
  did not declare `{'type': int}` (the same group's `--fetch-delay`/`--fetch-timeout` did declare it), so the command line
  handed `'50'`/`'3'` as strings to `cids[i:i + cq]` and `range(self.fetch_retries)`,
  measured as `TypeError: 'str' object cannot be interpreted as an integer`, happening before any request was made,
  while the configuration-file route (YAML's int) ran through (`缓存测试/reports/R46_fetch_route_types.json`,
  with the post-fix comparison in `R46_fetch_after_fix.json`). New class
  `neoHGT/tests/test_gtdb_pipeline.py::TestFetchCountCoercion` (7 tests): both routes the same answer,
  non-numbers refused by name; `<1` refused (`range(0)` is "wanted to fetch and did nothing"); bool does not count as a count;
  unset stays unset; and the last test uses a `urlopen` accounting stub to prove the **consumer** really no longer crashes halfway —
  testing the converter alone does not count as evidence.


  `TestProductionReachability` grew to 6 tests, and `test_a_wrapper_reached_only_through_a_rebinding_is_reachable`
  tests both directions. The measured profile effect is in `缓存测试/reports/R40_profile_effect.json`. Round 35 turned the fixture "excerpt" declaration into a checkable gate: `tools/verify_fixture_provenance.py`
  re-reads the read-only sources, comparing item by item sequence existence and residues, source sha256/byte counts, copy hashes, orphan fixture files,
  the published entry count of `db.faa` and the 10-genome cap; `neoHGT/tests/test_gtdb_pipeline.py::TestFixtureProvenance`
  (5 tests) builds a scratch tree for each of seven kinds of corruption to prove "every class gets reported", 0 problems on a clean tree, and verifies `--gtdb-dir` re-rooting;
  measured in `缓存测试/reports/R35_fixture_provenance.json`. Round 36 brought the reachability gate's **input** into the check as well: `neoHGT/tests/test_gtdb_pipeline.py::TestProductionReachability`
  grew from 2 tests to 5, requiring the real tree's frozen list to have mutually exclusive halves covering all existing definitions (0 `UNTRACKED`),
  a synthetic "function smuggled in" must be named; the old single-half list cannot answer the question; the three drift classes of `list_drift` —
  missing/stale/upstream changed — each tested once. Round 37 moved the lens onto the tests themselves: the new rule `tools/verify_plan_traceability.py::assertion_free_tests`
  scanned out 12 `test_*` with "no observable assertion at all", of which **two were our own**
  (`neoHGT/tests/test_gtdb_pipeline.py::TestPrecompMultiSampleGuard::test_single_table_is_allowed`
  and `::test_other_methods_are_untouched`, calling the guard but never looking at the result); they were changed to
  `assertIsNone` and `assertRaises` appearing in pairs, with a loop proving the guard distinguishes by method;
  the rule first scanned out 14, and after fixing our own 2, 12 remained — all 12 in upstream files,
  registered in the "pre-existing tests lacking assertions" section of `docs/交付物清单.md` without modification (RT-1 is purely additive);
  the rule's own test is `neoHGT/tests/test_gtdb_pipeline.py::TestAssertionFreeRule` (2 tests). Round 38 audited the "auditor" itself: the new rule `tools/verify_plan_traceability.py::checker_tools` checks whether every `tools/verify_*.py` has an executor (a CI step / nightly / another tool / a test), and `neoHGT/tests/test_gtdb_pipeline.py::TestCheckerToolUpkeep` (4 tests) requires "the set without executors to equal exactly the declared set", the declared tools' `--help` to still come up, synthetic orphans to be named, and nightly to reuse `tools/verify_read_only.py` instead of carrying its own weakened copy; the re-run evidence of the four declared tools is `缓存测试/reports/R38_ncbi_defaults.json`, `R38_hybrid.json`, `R38_resume_it6b.json`, `R38_resume_it6c.json`, `R38_resume_hybrid6.json`. |
| UT-3.1 | Covered | `neoHGT/tests/test_gtdb_pipeline.py::test_taxonmap_key_columns_are_identical` | The two taxonmap columns hold identical values and the row count = the db.faa sequence count (`neoHGT/tests/test_gtdb_pipeline.py::test_every_protein_gets_a_taxid` asserts every protein has a taxId); ID deformation of the kind rsplit would cause is locked down by a character-by-character equality assertion. |
| UT-3.1b | Covered | `neoHGT/tests/test_gtdb_core.py::test_route_c_lookup` | The local form hits via a stem-based reverse lookup; the official form hits via the prefix2genome degraded path, `neoHGT/tests/test_gtdb_core.py::test_route_c_fallback_via_prefix`; an unresolvable ID is reported explicitly rather than guessed, `neoHGT/tests/test_gtdb_core.py::test_unmappable_id_is_reported_not_guessed`. |
| UT-3.2 | Covered | `neoHGT/tests/test_gtdb_pipeline.py::test_six_columns_rejected` | 6 columns (staxids missing) → a clear ValueError; 7 columns pass normally, `neoHGT/tests/test_gtdb_pipeline.py::test_seven_columns_ok`. Round 31 added both values of this switch: `--require-staxids yes`/the default must both error immediately (RK5's protection), and only an explicit `no` degrades to a WARNING and continues — previously this option **was parsed, was converted, was stored on the object, but no code ever read it**, so `no` had no effect at all (test `neoHGT/tests/test_gtdb_pipeline.py::TestSearchChoiceValues`, which also covers, in the same class, the yes/no of `--tax-capital/--tax-latin/--tax-unique/--keep-placeholder-names` and every value of `--aln-method`, `--fetch-enable`). |
| UT-3.3 | Covered | `neoHGT/tests/test_gtdb_pipeline.py::test_route_a_then_c_then_unresolved_keeps_hits` | A fails → C hits → when both fail, unresolved_taxids.txt and unresolved_hits.tsv are written; --keep-unresolved keeps the hit records (`neoHGT/tests/test_gtdb_pipeline.py::test_unresolved_hit_dropped_when_not_kept` asserts that without keeping they are dropped and counted). |
| UT-3.4 | Covered | `neoHGT/tests/test_gtdb_pipeline.py::test_gtdb_mode_never_touches_update_dmp_files` | In GTDB mode, swapping update_dmp_files for a stub that throws and running the whole fallback chain never triggers it; the semantics of that function on the upstream NCBI path are still covered by `neoHGT/tests/test_search.py::test_update_dmp_files`. |
| UT-3.5 | Covered | `neoHGT/tests/test_util_gtdb.py::test_numeric` | All four kinds of input hit — numeric ID / prefixed label / filename-style label / unique prefix-stripped name (starting at `neoHGT/tests/test_util_gtdb.py::test_prefixed_label`); an ambiguous name errors and lists candidates (`neoHGT/tests/test_util_gtdb.py::test_ambiguous_name_lists_candidates`). |
| UT-3.6 | Covered | `neoHGT/tests/test_gtdb_pipeline.py::TestNoNetwork` | urlopen and util.run_command are replaced with assertion-raising stubs, and local ingestion + conversion makes zero calls throughout; the search and analyze sides likewise make zero calls (the end-to-end tests under `neoHGT/tests/test_gtdb_pipeline.py::setUp` assert on every run that calls is empty). |
| UT-5.x | Covered | `neoHGT/tests/test_gtdb_core.py::TestWorkspaceGuard` | Out-of-bounds paths (.. / absolute paths / $HOME / environment variables / symlinks) are all rejected, and writes toward the read-only input dataset are rejected; default outputs and scratch directories stay inside the workspace (`neoHGT/tests/test_gtdb_pipeline.py::TestWorkspaceGuards`, `neoHGT/tests/test_gtdb_pipeline.py::TestScratchDirectoryStaysInWorkspace`), and insufficient disk capacity errors out early in `neoHGT/tests/test_gtdb_pipeline.py::TestCapacityPreflight`. |
| UT-4.x | Covered | `neoHGT/tests/test_gtdb_pipeline.py::TestDonorReporting` | The donor is output in g__Escherichia form and additionally has an offline NCBI mapping (`neoHGT/tests/test_gtdb_pipeline.py::test_prefixed_labels_and_ncbi_map`); the profile is chosen per domain with archaea at closesize=50 (`neoHGT/tests/test_gtdb_pipeline.py::TestArchaeaProfile`); the close group size satisfies close_size (`neoHGT/tests/test_gtdb_pipeline.py::TestCloseGroupSizing`); an unreachable --close-tax is dropped with a warning (`neoHGT/tests/test_gtdb_pipeline.py::TestCloseGroupReachability`), and when the whole neighbor group is absent from the hits it no longer crashes inside sklearn but delivers a human-readable message before predicting (`neoHGT/tests/test_gtdb_pipeline.py::TestDegenerateGroupGuard`; the real trigger path is in `缓存测试/reports/R25_phylogeny_close_tax.json`: of 158 genus-level neighbors, 152 had no hits in this "one representative per family" library, and after filtering all 4,311 proteins were wiped out); the real-CLI end-to-end recovery instrumentation is in `neoHGT/tests/test_gtdb_pipeline.py::TestFixtureEndToEnd`. Round 31 actually ran every `choices` value on the analyze side and asserted the observed difference: `--weighted no` collapses to "nothing left to cluster" on the 5-protein fixture and is rejected by the guard; `--self-low yes` changes the predictions; `--outliers boxplot` predicts one fewer while `none/zscore` match the default; `--donor-label name` changes the donor column; `--donor-ncbi-map yes` additionally writes a `donor_ncbi.tsv` (`no` does not); only `--orphans no` prints "Removed N ORFans." (tests in `TestFixtureEndToEnd` and `TestSearchChoiceValues`; probe artifacts in `缓存测试/reports/R31_analyze_value_probe.json`). |

## §10.3 Integration tests

| ID | Status | Location | Notes |
| --- | --- | --- | --- |
| IT-0 | Covered | Probe and eight-item check in `缓存测试/reports/P0_spike_S0_layout.md` | The raw machine-readable artifact is `缓存测试/logs/layout_probe.json`; on mismatch the exit code is non-0 (the strict semantics of gtdb probe are covered by UT-0.3). In Round 24 the "mismatch" branch was genuinely run on two other local releases: `GTDB_r226` (its lineages glued together with `_` by a third party and truncated to genus level) and `GTDB_r214/parent_directory` (the official flat layout, no metadata table) were both refused with exit code 2, each printing one `ERROR: …` line naming what is missing / how to pass it, **creating no output at all** (`neoHGT2GTDB-临时文件/r226_clean.log`, `neoHGT2GTDB-临时文件/r214_probe.log`; test `neoHGT/tests/test_gtdb_pipeline.py::test_refused_dataset_reports_one_line_instead_of_a_traceback`). |
| IT-1 | Covered | `缓存测试/db/r220_taxdump/taxdump/nodes.dmp` | The full r220 metadata/taxonomy converted; the measured node count and rank counts are recorded in `缓存测试/reports/ACCEPTANCE.md`; the artifacts can be read back by util.read_taxdump. |
| IT-2 | Covered | `缓存测试/db/r220_entero/meta.json` | A small library of the f__Enterobacteriaceae representative set: db.dmnd is valid, the genomes.tsv row count matches the family's representative count, and meta.json contains provenance and selfcheck results (same directory tree). |
| IT-2b | Covered | `缓存测试/db/r220_ar_family/meta.json` | The fast tier --rank-sample family (564 archaeal genomes) measured 2 min 37 s / 360 MB / peak 985 MB, already usable as a CI tier; the budget is recorded in `缓存测试/reports/P0_spike_S3_cost.md`. |
| IT-3 | Not covered | The GenomeFetch and --seq-source dispatch of `neoHGT/gtdb.py` | The missing half: **the --seq-source ena / gtdb-tar / ncbi supplementary-fetch routes have never actually been executed**; downloading and network access are not authorized. What is covered is the code path itself — the URL and checksum logic have unit tests (`neoHGT/tests/test_gtdb_pipeline.py::test_ncbi_and_ena_channels_keep_their_urls`), and the diff sets are written to disk by the reconcile manifest (DoD-2). Once the user allows network access these can be executed directly, with no code changes needed. Round 24 additionally covered the **ingestion side** of the "official layout" half: `bac120/ar53_taxonomy_rNNN.tsv(.gz)` placed in the root after unpacking is now discovered, and `gtdb taxdump` can produce a taxdump on an official dataset with **only lineage files** (measured on local r214 at 7.9 s / 113,318 nodes, with `缓存测试/gtdb_r214_taxdump/meta.json` recording `lineage_only: true`); still missing are **running one full database build on an official dataset** (there is no official metadata table to try locally: the r214/r226 copies carry only lineage files; though the column mapping itself was calibrated against the official schema — all 11 of the `REQUIRED_COLS` columns hit the header of the real `bac120_metadata_r220.tsv` 11/11), and the download itself. |
| IT-4 | Covered | `neoHGT/tests/test_gtdb_pipeline.py::TestFixtureEndToEnd` | Two pieces of evidence coexist: the end-to-end on real data (`缓存测试/search/halo/hgt_demo.json`, 371 predictions + donor labels + scatter.pdf), and the **offline reproducible** version added this round — a KB-scale fixture, the real CLI parser, asserting that the 2 buried hits are recovered and the donor is g__Nanoarchaeum. |
| IT-5 | Covered | `缓存测试/reports/IT5_hybrid_r3.json` | On the real 6-genome set, hybrid is byte-for-byte identical to pure GTDB in DIAMOND hits, search tables, and scores.tsv; the graft happens on a copy (146,380 vs 146,404 nodes), and the NCBI-side viral/mammalian donors still await reference-library authorization, see D12. |
| IT-6 | Covered | `缓存测试/reports/IT6e_resume_hybrid.json` | Rerun after interruption: artifact manifests identical, resume decisions and first-run results deterministic; when library files are missing it refuses to "resume" and states what is missing. |
| IT-7 | Covered | `缓存测试/reports/IT7_workspace.json` | Mechanical checks of all project files, of escaping symlinks, and of scratch-directory disposability (both probes PASS); the read-only dataset audit is in `缓存测试/logs/read_only_check_r15.txt`. |

## §10.4 Scientific validity verification

| ID | Status | Location | Notes |
| --- | --- | --- | --- |
| 10.4-1 | Not covered | `tools/run_hgt_demo.py`, `缓存测试/search/halo/hgt_demo.json` | The missing half: **the head-to-head comparison of the same batch of queries on the GTDB library and on the NCBI reference library has never been done** — the NCBI reference protein library requires an online download, which is not authorized. The GTDB-side artifacts (HGT set, donor labels, scores.tsv) have been written to disk by run_hgt_demo, and the donor→NCBI offline mapping needed for the comparison has also been implemented (the ncbi_map of `neoHGT/tests/test_gtdb_core.py` and analyze's donor_ncbi.tsv); once the reference library is obtained this can be executed directly, with no code changes needed. The "partial delivery" registered under the same ID in §12 D12 is the same matter. |
| 10.4-2 | Covered | `缓存测试/reports/D12_markers.json` | False-positive check on conserved marker genes: of the ar53 representative-set markers in the local copy (53 single-copy marker genes, stricter than the 30 ribosomal proteins the plan used as an example), those judged as HGT on the same GTDB library number **0**; the command and its interpretation are in `tools/run_marker_fp.py`, and it remains 0 after switching to a phylogenetic-distance definition of the close group (`缓存测试/reports/D12_markers_phyloclose.json`). |
| 10.4-3 | Covered | `缓存测试/reports/S10_4_3_sensitivity.md`, `缓存测试/reports/S10_4_3_all_samples.json`, `缓存测试/benchmarks/` | **Three samples, four curves**: halo (archaeal f__QDMS2, 4,025 proteins / 75,742 hits), bathy (Thermoproteota of a **different phylum**, 2,303 proteins / 36,015 hits, `缓存测试/search/bathy/hgt_demo.json`), entero (a bacterial boundary sample; a single-family library has no distant donors so it is always 0, used only for input_cov). close_size is monotonic on both samples that have distal hits; bandwidth is the knob with the largest effect but **its magnitude depends on the sample** (halo differs by up to 2.6×, bathy 1.3×); input_cov is completely flat at 50–95 on both archaeal samples, while on the bacterial sample with 89.47% coverage cov=95 degrades the input classification from species to genus — i.e., this option's point of action is the input classification itself; the direction of noise is **opposite** between the two samples (contradicting the config comment "the larger the stricter"), and this has been written into the usage documentation. In addition, bathy reusing the scores table via `--from-scores` gives the same values as a full run (724→367), proving that the reuse path fixed in Round 18 also holds on real data. Note: `search/halo_p2/query.tsv` has the same sha256 as halo's query.tsv and was therefore not used as a second sample. Round 25 expanded the bacterial side from 1 sample to 4 (the same 4,896-genome library, `close_size=50`: 771 / 753 / 714 / 436 predictions, with 100% of hits carrying staxids), see `缓存测试/reports/R25_bacterial_multi_sample.json`; Round 26 then added four-tier `close_size` curves for three of those samples (predictions only, without rerunning DIAMOND): three samples are completely flat across 10/50/100/200 (771/753/714), while **the fourth (Arcobacteraceae) drops to zero from 100 on**, because its self group has only 2 taxa and at `close_size≥100` the neighbor group swells to the entire phylum and swallows all hits, see `缓存测试/reports/R26_bac_close_size_curves.json` — so "bacteria tier 200 usable" is a per-sample property, not a library property; the same round also discovered and fixed three fabricated names in the "skipped check items" list under lineage-only mode (`neoHGT/tests/test_gtdb_pipeline.py::test_reported_skips_are_real_check_names`); this round also measured the incompatibility of "phylogenetic neighbor group + sparsely sampled library": of 158 genus-level neighbors, 152 had no hits in that library, and once the neighbor group empties, upstream dies inside sklearn; this was changed to give a diagnosis before predicting (`缓存测试/reports/R25_phylogeny_close_tax.json`). |
| 10.4-4 | Not covered | The hgtdb_20230102 download point of `neoHGT/server.py` | The missing half: **the comparison against HGTector's official prebuilt library has never been executed**; that library requires an online download and is not authorized. The code side has already separated the two variables "change only the taxonomy" and "change the reference sequences too": the three paths --source gtdb / hybrid / ncbi do not shadow one another, and hybrid and gtdb have been proven byte-for-byte equivalent on the same sequences (IT-5), so this item is short exactly one reference library. |

## §10.5 Regression and compatibility

| ID | Status | Location | Notes |
| --- | --- | --- | --- |
| RT-1 | Covered | `缓存测试/regression/RT1_failure_set_diff.txt` | The failure sets of the baseline and after the adaptation are compared item by item with comm (this round 13/13, 0 new failures), and 5 upstream files have deleted_lines=0; the byte-for-byte comparison of Ecoli fixture artifacts is in `缓存测试/reports/RT1_ncbi_fixtures.json`; in CI it is a real gate (`.github/workflows/ci.yml`). |
| RT-2 | Covered | `缓存测试/reports/RT2_defaults.json` | Option by option, the argparse defaults before and after the adaptation are checked to be identical; the GTDB profile's relaxations happen only inside the --source gtdb/hybrid branches. |
| RT-3 | Covered | `neoHGT/tests/test_util.py::test_is_latin` | The signatures and behavior of util.py's existing functions are unchanged (this round's deleted_lines=0 on upstream files is its mechanical evidence); the GTDB semantics of rank_plural / is_latin / is_capital are implemented uniformly as new functions + provider dispatch, and the old assertions pass as-is (`neoHGT/tests/test_util_gtdb.py::test_rank_plural_legacy_semantics_kept` locks rank_plural's legacy semantics). |

## Rounds 47 / 48: register rule 19 (configuration route vs command-line route)

Round 47 turned the Round 46 census into rule 19, with the new class
`neoHGT/tests/test_gtdb_pipeline.py::TestRouteEquivalenceRule` (initial version 5 tests):
a planted divergence with no explanation must be named; a stale exemption left behind must be named;
a missing ledger file must be `ledger missing` + problem (not "an empty problem list" read as a pass);
a walkthrough that parses too few readings must refuse to vouch for itself; the real tree has 0 problems and "explanation count == divergence count".

Round 48 found the rule **looking at 16% less than the evidence it claims** (58 keys vs the census's 67 keys):
the section name of `f'{self.method}.{key}'` is actually written in the `choices` of `-m|--method`,
and the rule rejected it outright as "only knowable at runtime". After adding expansion by the declared choices the ledger went 58 → **73 keys**,
`unresolved_call_sites` went to zero, and **the divergences are still the same 8** (the wider net hid no defects — this sentence was written after measuring, not before).
The class therefore grew to **7 tests**:

| Test | What it pins down |
| --- | --- |
| `test_a_declared_slot_expands_and_an_undeclared_one_is_reported` | Pins both ends: if `choices` is declared it must be expanded; remove `choices` and it must fall back to `unresolved_call_sites` **stated explicitly**, not quietly audit less — the latter is how the gap grows back |
| `test_the_shipped_ledger_refuses_nothing_and_marks_its_expansions` | Today no call site is refused (`unresolved_call_sites == []`), and every key brought in by expansion carries an `expansion_basis` stating why it is in scope |

The ledger itself is `docs/路由等价性清单.json`, committed into the repository (73 keys / 8 divergences / 8 with explanations).
**Checking it needs no PyYAML** (the value types were baked in when the ledger was emitted), so this rule runs on a bare checkout too;
only re-emitting the ledger with `--emit-route` needs PyYAML; that step runs in every round's gate,
and the gate requires "the re-emitted copy to be byte-for-byte identical to the one in the repository"; if they differ it fails with `LEDGER STALE`.

## Round 49: rule 19b (the "reason" sentence written in the exemption must itself be checked)

Rule 19 only looks at "whether there is an explanation"; an explanation is prose, so the prose can go stale while the gate stays green. Round 49 turned the 8 prose sentences
into decidable objects, adding the class `neoHGT/tests/test_gtdb_pipeline.py::TestRouteClaimVerification`
(**7 tests**: the first 5 establish the contract, the last 2 joined this round together with the citation check):

| Test | What it pins down |
| --- | --- |
| `test_a_sound_claim_passes_and_says_which_lines_it_read` | A good reason must pass, and report the line numbers it read (`function` / `read_line` / normalized line); behind a green light there must be readable evidence |
| `test_a_claim_naming_nothing_real_is_a_problem` | When the normalized function named by the reason does not exist in the module, it must be a problem, and it must not enter `verdicts` ("not judged" and "judged pass" must not be mixed into one number) |
| `test_normalising_before_the_read_is_a_problem` | The same call, the same key table, just one line earlier: normalizing before reading the value equals not normalizing at all — line order is the only difference, and the check must see it |
| `test_an_unorderable_claim_needs_a_test_that_exists` | When a loop-form key cannot be ordered, the reason must point to a test that really exists; pointing to a nonexistent test is a problem |
| `test_the_shipped_claims_are_all_checked` | The real tree has 0 problems, and `claims_checked == 分歧且写明理由的键数` (not a single one skipped); if the ledger's shape changes it must error, not pass vacuously |
| `test_a_reason_naming_a_function_that_is_not_there_is_a_problem` | The **location** in a reason is also an assertion: `mod.py name` must be defined in the module; this test also proves the three normalization ones pass as usual, so the failure can only come from the citation |
| `test_the_shipped_reasons_cite_only_things_that_exist` | On the real tree all 14 citations resolve, `citations_unjudged == []`; `citations_seen >= 键数` guards against the idle state of "reasons no longer name locations, so there is nothing to check" |

Mute mutant and mutants: `缓存测试/reports/R49_test_falsification.json` (muting `route_claim_problems()`
entirely → **7/7 of the class fails**, including the two tests asserting "no problems") and
`R49_citation_falsification.json` (3 mutants of the citation check: muted entirely / ignoring only `mod.py name` /
always reporting bad; the expected tests catch each one, and the same unmutated copy has 0 failures).

## Round 50: rule 20 (the defaults the documentation promises vs the defaults the code delivers)

This round both added the rule and fixed one real product defect. The class has **7 tests** in total (the first 6 landed with the rule; the 7th was added after gate pass A
to cover the delivery form `default=<常量引用>`).

| Test | What it pins down |
| --- | --- |
| `neoHGT/tests/test_gtdb_pipeline.py::TestFixtureEndToEnd::test_the_donor_label_default_is_the_one_the_docs_promise` | Run one end-to-end with `--donor-label` removed entirely; the donor column must be the documented `g__Nanoarchaeum`. Before the fix it wrote the synthetic taxId `105000000` — this sentence was run out, not read out |
| `TestDefaultPromises::test_every_shipped_promise_agrees_and_none_is_unjudgeable` | On the real tree all 8 declarations are judged "agree", `unjudgeable == []`, `config_state == 'read'`; "pass" must be judged into existence, not be an absence of looking |
| `TestDefaultPromises::test_a_delivery_that_moved_is_reported_as_disagreement` | A delivery moved in code (`or 'local'` → `or 'remote'`) must be named, and the remaining declarations pass as usual — one change must not mute the whole rule |
| `TestDefaultPromises::test_a_promise_deleted_from_the_docs_is_reported` | Conversely: deleting the "default…" sentence from the docs is also a problem; the gate does not vouch for sentences nobody maintains |
| `TestDefaultPromises::test_an_establishing_function_that_is_gone_is_not_a_pass` | The declaration points at `write_hgt_list_gtdb()`; after the function is renamed it must error, not say "all good" |
| `TestDefaultPromises::test_a_profile_default_without_config_is_unjudgeable_not_agreed` | With the per-domain tiers removed, `close_size` lands in `unjudgeable` and must **not** be counted as agreed |
| `TestDefaultPromises::test_a_table_default_written_as_a_constant_is_resolved` | `default=DEFAULT_ID_BASE` is also a delivery: unmutated, `gtdb_taxid_base=100000000` lands in `agree`; after changing the fixture constant to 900000000 it must land in `disagree` — proving this is not "read it and move on" |
| `TestDefaultPromises::test_a_tree_where_nothing_can_be_judged_says_so` | Breaking all four kinds of establishing point (including in-table `default=`) → the rule must report of itself that "not a single one was judged" |

Mute mutant: `缓存测试/reports/R50_rule20_falsification.json` (replacing `default_promise_problems()` with
"report a clean table" → **7/7 of `TestDefaultPromises` fails**, re-verified by the pass-B gate).
The intermediate artifacts of the scope measurement and the two probe errors are recorded in
`缓存测试/reports/R50_default_scope.json` and entry 35 of `docs/GTDB-数据库.md`.

## Round 51: rule 21 (a value that can be None must not appear in an arithmetic position)

`neoHGT/tests/test_gtdb_pipeline.py::TestNumericNoneGuards` (**7 tests**). This rule watches the common shape of
three real defects in this project's history (Round 45 `isinstance(bw, float)`, Round 46
`cids[i:i+cq]` / `range()`, Round 50 the donor column written with a synthetic taxId).

| Test | What it pins down |
| --- | --- |
| `test_an_unguarded_range_argument_is_reported` | `range(self.narrow)` with no guard → must be named (this is the shape of the Round 46 crash) |
| `test_a_truthiness_guard_makes_the_same_code_pass` | The same code wrapped in `if self.narrow:` → must not be reported. The only differences are line order and the guard, and the rule must see that difference |
| `test_a_message_is_not_arithmetic` | `'… %s' % self.x` is formatting, not arithmetic — without this distinction `--taxdump` would report 6 false defects at once |
| `test_a_config_value_counts_as_delivery_and_stops_the_report` | `get_config(self,'narrow','search.narrow')` + the config supplying a value = delivered, no longer counts as "possibly None" |
| `test_an_attribute_reassigned_from_disk_is_not_the_options_value` | The attribute is overwritten by another source before use (`self.x = load(self.x)`) → this numeric use is not charged to the option, and the excluded count must be written into the report |
| `test_a_tree_where_every_option_has_a_default_is_refused_as_blind` | All 24 options have `default=` → the rule reports of itself "not all judgeable forms recognized" (on a small tree this rule deliberately stays silent, see `VACUITY_MIN_OPTIONS`) |
| `test_the_shipped_tree_is_clean_but_not_because_it_was_unexamined` | The real tree has 0 problems, and `options >= 150`, `can_be_none >= 50`, `numeric_reads_in_all_cli_attrs >= 1`, `config_state == 'read'` — clean must be the clean that comes after having looked |

Mute mutant: `缓存测试/reports/R51_rule21_falsification.json`. The criterion is an **explicit expected set**:
the 4 tests asserting "findings/counts exist" must fail, while the 3 tests asserting only "nothing was reported" pass naturally under a deaf rule;
they are named explicitly in `cannot_notice_under_a_mute` and do not count as coverage.
(The previous version reused r69/r80's "failure count == the N failed the suite reports"; those two numbers come from the same summary line and are **identical by construction**,
so it never checked anything — this is written in entry 36 of `docs/GTDB-数据库.md`.)
The discrepancies among the six readings of the delivery census and their root cause (20 config leaves of "key present, value empty") are also in entry 36.



## Round 52: the security rules' enforcement points checked in reverse

`TestRunnerChangeList` (6 tests, all running the **released** `tools/run.sh`, inside a synthetic workspace):

| Test | Why "always list everything" cannot paper over it |
|---|---|
| `test_created_file_is_reported_with_its_category` | A created file must appear with the `CREATED` category |
| `test_deleted_file_is_reported_as_gone` | A deletion must appear as `GONE`; "(none)" is no witness that "nothing changed" |
| `test_rename_inside_the_tree_is_reported` | `mv` preserves mtime and the `-newer` filter is completely blind to it; now both endpoints must appear |
| `test_quoted_argument_reaches_the_program_intact` | `-- python -c "print('x')"` must arrive intact (the old `eval "$*"` died on a syntax error) |
| `test_a_command_that_changes_nothing_lists_nothing` | A no-op must report nothing and the scratch directory's file set must be unchanged — without this, the first four can be satisfied by "list everything every time" |
| `test_the_deletion_reaches_the_command_log` | What the terminal saw must land in `commands.log` (the `GONE` line) |

`TestIntroducedConfigKeys` (8 tests) watches the config-key rule's **judged population**: the population must be
"the diff against the upstream config.yml", not "the key name contains gtdb"; when upstream is unavailable it must fall back to the name method and say so in the readings,
and "did not look" must not blend into "looked and found no problem" (`test_an_unreachable_upstream_dir_is_not_a_vacuous_pass`:
`upstream_dir=''` resolves to the released file being compared itself, and diffing it against itself yields an empty set — the greenest not-looking of all);
and the same rule must **bite** when a key nobody reads is planted into it (`test_the_rule_population_is_the_upstream_diff`
only asserting a clean tree of 0 would be satisfied by a deaf rule, so a planted-key paragraph was added). Three items on the product side: `database.source` is read and the command line still wins,
it is still ncbi under the released configuration, and `workspace_root()` really takes the path its docstring promises (env first).

Two mute mutants, both judged by **expected sets named per test**, and both requiring `skipped == 0`
(the copy workspace initially lacked the upstream copy, and two testimonies silently skipped — that kind of "pass" is worthless):
`R52_runner_falsification.json` (reverting to `find -newer` → tests 1-4 fail; reverting to `eval` → only test 6 fails;
unmutated 6/6 pass) and `R52_config_rule_falsification.json` (deaf rule → 2 tests; removing the `database.source` wiring → 3 tests;
removing the workspace lineage branch → 3 tests; population reverting to the name method → 1 test; unmutated 11/11 pass).
The readings of the six cases of the self-check script `临时文件/r100_runner_inventory_selftest.sh` are archived together with
`r99`'s "before/after" probe (including a positive control: writing a file is seen, rename and deletion are not → falsifying the old mechanism).

This round's process ledger (also kept as files): `r107`'s first version asked "did the key name ever appear as a string" and got "all 26 appeared" —
**a check that always answers yes is not a check**; only after switching to requiring a read shape of `.get('k')`/`['k']` did it surface 5 candidates,
from which 3 real defects were pinned down manually; the expected sets of `r101`/`r109` each missed once (the former assumed the old runner was blind only to rename/deletion
when in fact the category name was lost too; the latter's anchor appeared in two functions, and a blind edit turned the other one into a NameError).

### IT-7's external-write exemptions must be bitten too (4 tests added the same round)

`TestExternalWriteGate`'s namespace attribution was judged red once in Round 52: a new file appeared in
`~/.cache/pip/selfcheck/`, while measurement (`R52_pip_attribution.json`) proved that our four gate-shape commands wrote
0 files there and `commands.log` never contained a pip call. The attribution direction was right (commands bypassing the runner leave no trace),
but it books **what others did** to our account, and a gate like that sooner or later gets treated as noise. So exemptions may only be "used with evidence":
4 new tests assert respectively: evidence present → we are not judged; evidence missing → judged; the evidence's conclusion changed → judged;
an undeclared subpath (`pip/http`) → judged as before; declaring an entire namespace is rejected outright in the code.
`tools/env.sh` exports `PIP_CACHE_DIR`/`PIP_DISABLE_PIP_VERSION_CHECK` so the precondition is enforced rather than assumed.

## Round 53: an exemption's evidence must be newer than the code it measures

`TestExternalWriteGate` grew from 10 tests to 16; the 6 new ones specifically watch "can the evidence lie":

| Test | The failure mode it rules out |
|---|---|
| `test_evidence_older_than_the_code_it_measures_is_refused` | Evidence older than the newest code in scope → the exemption is refused, and it **reports which probe should be rerun** (a red gate must be actionable) |
| `test_touching_the_evidence_file_does_not_refresh_it` | `touch` cannot renew it: the evidence's own stated `newest_measured_code_mtime` is what counts |
| `test_a_stated_measurement_survives_an_old_file_mtime` | The reverse: an mtime made older by copying while the statement is still fresh → valid (both directions tested, otherwise the rule is "mtime worship") |
| `test_an_exemption_over_an_empty_scope_is_refused` | When no executed file is in scope it is **refused** rather than agreed (what was not measured cannot be vouched for) |
| `test_the_pip_detector_reads_the_command_field_not_the_prose` | The detector reads only the command field; in Round 42 an intent whose prose said "pip redirection" was once counted as a real pip call |
| `test_the_freshness_scope_is_the_code_the_gates_execute` | scope = `neoHGT/*.py`, `neoHGT/tests/*.py`, `tools/*.py|sh`; `docs/*.md` is not in it (documentation is not executed code), and a missing scope must be reported as `None` |

Mute mutants: four machines in `R53_freshness_falsification.json`, with distinct fingerprints
(deafening the comparison → 2 tests, an empty scope pretending to agree → 1 test, the detector reverting to the whole line → 1 test, taking the mtime as authoritative → 3 tests),
the unmutated copy passes 16/16 with `skipped_tests: []`. This machine itself had to be fixed twice before it was trustworthy:
once it forgot `re.M` (`^` anchors only the start of the whole string → parsing always empty), once it added `-rs`
(**pytest's `-r` replaces the default set rather than appending**, so the `FAILED` lines disappeared);
both times it produced `effective: false` instead of green. The third pass also corrected my mental model of my own code's reach
(deafening the "old-vs-new comparison" does not disable the "empty scope" guard, so that test must still pass).


## Round 54: rules 23 and 24 (config sections selected by value / section headings the gate reads) and entry-string reads

The register gate previously only recognized **subscript**-style config reads (`cfg['database']['gtdb']`) and could not see the upstream-idiomatic
`get_config(self, key, f'search.{key}')`; measurement showed the path sets of the two read styles are **completely disjoint**
(26 vs 73, intersection 0, `缓存测试/reports/R54_config_reads_precision.json`).
After the scan was added, the first version merged the literals of the five `for key in (…)` loops across the whole function, conjuring 409
"reads" out of thin air — after switching to expansion scoped by call site it converged to 73, of which only 6 point to sections that do not exist in the file
(`auto.*`, `precomp.*`, from the same call site `neoHGT/search.py:306`).

| Test | Criterion |
| --- | --- |
| `TestConfigEntryExpansion::test_two_loops_sharing_a_name_do_not_cross_multiply` | Loop variables sharing a name must not leak into each other: the whole set must equal exactly 4 entries, and even one extra entry is a failure (this round's real defect was precisely "too many") |
| `…::test_an_inner_loop_shadows_the_outer_one` | The inner loop shadows the outer one; the outer loop's values are not candidates |
| `…::test_a_renamed_leaf_still_expands` | Renamed leaves like `f'grouping.{key.replace("_","")}'` are still inside the expansion (missing it would create a false "nobody reads this") |
| `…::test_a_call_it_cannot_expand_is_reported_not_invented` | A `for key in 'e':` iterating a **bare string** is not guessed at: it goes into `unexpanded` with the line number named |
| `…::test_the_subscript_resolver_is_call_scoped_too` | The same call-scoping rule was also fixed into `config_reads()` (it originally merged by "the last loop overrides the previous one") |
| `…::test_the_shipped_tree_forms_no_impossible_paths` | Real tree: 73 reads, 0 spots that cannot be expanded, and `grouping.name`/`diamond.threads`/`remote.method` no longer appear |
| `TestDynamicConfigSections` (11 tests) | Six directions of rule 23: an undeclared dynamic section goes red; a two-line disposition with evidence turns green; a refusal naming one option too few → red; the evidence line not being a `raise` → red; a declaration omitting one covered option → red; the normalizing `if` leaving a gap → red; the attribute assigned back to its original value after normalization → red; declaring a section the code already forms → red; a register that cannot be read ≠ pass; a section name that really exists in config.yml needs no declaration; **change one side of the shipped refusal and the rule must go red again** |
| `TestLedgerRowDetection` (3 tests) | `ledger_rows()` was changed to recognize ledgers by their header: a wide prose table is not taken for a ledger, ledger rows spanning prose are still read, and the shipped deliverables table reads out exactly D1–D12 |
| `TestPrecompRouteRefusesBatchKnobs` (5 tests) | Each of the three options is refused by name; when given together all three appear; a non-precomp route is left alone; `--queries 0`/`--maxchars 0` (same semantics as this route) are not refused; `guard_precomp_knobs` really hangs in the GTDB parameter route and hangs exactly once |
| `TestIntroducedConfigKeys::test_a_key_read_only_through_an_entry_string_is_still_a_read` | After entry strings join the judged population, the rule still bites: a newly added "nobody reads it" knob is still reported, and dynamic sections already handled by rule 23 are not re-reported by this rule |

| `TestLedgerSectionHeadings` (5 tests) | Rule 24: the 7 section headings the gate itself must read must still be headings — a missing heading must name **which reader would go mute**, a prose mention of the word does not count as a heading, the shipped documents are 7/7 complete, an unreadable file may not pretend to be clean, and this round's incident is replayed once (deleting a heading line must go red) |
| `TestConfigKnobs::test_a_dynamic_section_is_reported_only_until_rule_23_owns_it` | The handoff between rules 22 and 23 is asserted **in both directions**: with the register table, `read_but_absent == {}`; without it, the three `precomp.extrargs/maxchars/queries` entries must be reported back as-is with `problems` non-empty (a missing table is not cleanliness), after which rule 23 takes over and at that point has no problem of its own |
| `TestFixtureEndToEnd::test_the_shipped_cli_route_refuses_the_inert_knobs` | The refusal must hold beyond being called directly: along the full route of **the shipped parser + `Search()(args)`**, both `--queries 5` and `--extrargs=--seg` must raise `ValueError` naming the flag (the other 5 tests only call the guard as a function, which cannot prove the command line reaches it) |

Two other tests of `TestIntroducedConfigKeys` were changed to call the rule **with the register table**: `config_key_problems` now needs
to read the disposition table to know which dynamic sections are not its business; passing `ledger=None` would re-report `precomp.*` as "read but absent" —
this is not a defect but the handoff, so the tests pass the register table in, and separately assert the shape "it reports when no table is passed".

The same round also fixed `tools/run.sh`'s move pairing, and `TestRunnerChangeList` grew from 6 tests to 8:

| Test | Criterion |
| --- | --- |
| `test_a_bulk_move_pairs_by_name_not_by_hash_order` | When bulk-`mv`-ing files of the same name, the two ends of every `MOVED` **must match by filename**, and `moved(ambiguous)=0` |
| `test_an_unpairable_move_is_declared_instead_of_guessed` | When two files share size+mtime but the names do not match up, it must **not** claim a pairing: `moved(ambiguous)=2`, no `MOVED` lines, yet the files on both ends still appear in the manifest |

Negative control `临时文件/r130_move_pairing_falsification.py`: installing the name-blind old loop back into two copy trees,
the second test (the deterministic witness) **fails on both passes**, and the remaining tests must not move on either pass; the first test passes one pass and fails the other —
its readings depend on awk's traversal order, so the report records it as `order_dependent_lock` rather than coverage
(the instability itself is the defect that was fixed).

### The six tests added after pass C: how the gate reads its own documents (`TestDocSectionLocator`)

The 6 tests that went red in pass C were not a regression but the gate exposing one of its own reading habits: sections were located with
`text.split(token, 1)[1]`, and I had quoted that section's heading verbatim in the D5 line, so the slice landed on
`[23, 31]` instead of `[363, 374]`. The fix is a shared `doc_section(text, start, stop=None)`:
the start must be the line containing the heading, and the end is the next heading of the same or higher level (or a `stop` named by the caller).

| Test | Criterion |
| --- | --- |
| `…::test_prose_mention_does_not_move_the_slice` | Mention the section name once in prose first, and the slice still starts at the heading and still ends before the next section |
| `…::test_a_slice_ends_at_its_own_heading_level` | A `###` subsection stays inside the section; only the next `##` closes it |
| `…::test_a_token_that_never_names_a_heading_has_no_section` | A name mentioned only in prose returns an empty string — "cannot be read" and "that section is fine" must be distinguishable |
| `…::test_the_stop_token_wins_on_the_plan_layout` | An appendix like the plan's, with **same-level subdivisions**, needs `stop` to close; without capping it would swallow the next appendix (this test was originally a weak assertion; it was strengthened only after the mute mutant of the `r136` survey pass exposed that it could not catch the mutation) |
| `…::test_the_shipped_ledger_hazard_is_still_there_and_still_harmless` | First assert that the D5 line **really does** still quote that section's heading, then require rule 23 to read the two-line disposition: deleting that quote deletes the precondition of this regression, and that must be seen |
| `…::test_no_reader_locates_our_ledger_by_first_mention` | Pins at the source level that all six readers go through `doc_section` and that first-mention spellings are always judged red (this machine recognizes only text shapes, so only mutations changing the call-site spelling trigger it — recorded honestly in the report as a boundary) |

Negative control `临时文件/r136_docloc_falsification.py`: four machines respectively revert rule 23's reader to first mention,
revert the appendix C/D readers to "to the end of the file", degrade `doc_section` wholesale into not recognizing headings, and make it read nothing at all;
the named lists come from the **survey pass** (`R54_docloc_survey.json`) rather than judgment, requiring 8/1/4/17 test failures,
`skipped: 0`, and no unplanned failures. The copy tree needs to see the upstream configuration for those two tests to run, so the machines byte-copy the upstream
`neoHGT/config.yml` to their own relative location and prove the match with sha256 — otherwise the control group skips two tests at the start
and reports coverage it never measured.

## Round 55: the round-level change witness `tools/round_witness.sh` (the gap in security rule 16) and the three documents I broke myself

| Test | Criterion |
| --- | --- |
| `TestRoundWitness::test_writes_between_commands_are_named` | A file written outside a command must appear in this round's manifest (the same test also proves that writes inside commands appear too); deletions must appear as GONE — an mtime scan cannot answer "who is no longer there" |
| `…::test_the_witness_does_not_report_its_own_snapshots` | The witness's own snapshot names must not appear in **change lines**; the first version wrote the assertion over the whole text and was blocked red by the legitimate "snapshot provenance" in the header — provenance may name them, change lines may not |
| `…::test_stop_without_start_refuses_and_writes_no_report` | stop without start must exit non-zero **and write no report**: an empty witness must not treat "no changes" as a result |
| `…::test_restarting_a_round_is_refused_and_the_window_survives` | Restarting within a round must be refused (otherwise it would erase every write made before the restart), and after the refusal the earlier write must still be in the final report |
| `…::test_a_round_that_changed_says_it_in_words_not_a_blank_list` | A round in which nothing changed must write "NOTHING CHANGED -- 若非属实，说明见证启动晚了"; it may not hand in a blank table |
| `…::test_a_round_number_that_is_not_one_is_refused` | The round number must be 1–3 digits; `0`/`abc`/`55x`/the empty string are all refused (a tool that swallows a wrong round number would write the manifest into the wrong round) |

The self-test `临时文件/r140_round_witness_selftest.sh` is five bash cases; besides the four above,
**the second case attacks the premise**: it manufactures "one write between two commands" in the same throwaway tree,
confirms that the change list of `tools/run.sh`'s subsequent commands really does not have it while the round-level witness does —
if it held the other way around, this tool would be redundant work, and the script must say so honestly rather than keep being repaired.
Reading: `verdict: PASS`; the things that went red in the first two passes were defects of this test of mine
(the nested `run.sh` failed to produce the target file; case 3's judging condition was written backwards), not of the tool under test.

This round also has three mistakes of my own to record (their shapes will recur):
an unquoted here-doc swallowed all backtick-wrapped command names in the document as shell command substitutions
(`docs/交付物清单.md` was at one point written as "新工具 ：把…"), and the damaged version is archived at
`backup/*.damaged_r55_heredoc`; after restoring from the backup, insertion was changed to go through a **file** with a read-back assertion;
one assertion was written too wide (the second test in the table above); and `rm -rf` was used once on a scratch tree I had built myself —
rule 14 does not care who the author is; from then on everything is archived, never deleted.
