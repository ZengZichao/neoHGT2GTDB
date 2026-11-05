> English (primary) | Chinese mirror: [搜索层假设核对.md](../搜索层假设核对.md) — keep the two in lockstep. Chinese mirror linked above.

# Plan §2.3 search-layer assumptions check (S1–S9)

Plan §2.3 lists item by item the **code-level assumptions** `search.py` makes under NCBI data shapes, and for each one writes out
the "GTDB consequence". These conclusions are decidable as true or false: either they still hold in today's code, or they have already been changed, or the plan's
inference at the time was itself wrong. This register checks all 9 items, and accepts
the machine check of `tools/verify_plan_traceability.py` (the same rule set as the other registers):

1. Every number in the §2.3 table must appear exactly once;
2. Status may only take `已按GTDB改造` / `确认无需改动` / `与方案预期不同并已调整` / `未处理`;
   the latter two kinds (including "differs from the plan") must state their basis clearly;
3. The referenced paths and `::符号` must really exist.

## S1–S9

| ID | Status | Location | Description |
| --- | --- | --- | --- |
| S1 | Adapted for GTDB | `neoHGT/search.py` (`set_params_gtdb()`, `open_taxmap_index()`), `neoHGT/tests/test_gtdb_pipeline.py::TestTaxmapIndex` | Option names are reused as before, but `-t` must point at the **in-database** `taxdump/`, and when `--taxmap` is switched to a GTDB database it becomes a disk index (`TaxidIndex`) so that 330 million mappings are not read into memory; the GTDB face additionally gains `--source`/`--gtdb-stem2taxid`/`--keep-unresolved`/`--require-staxids` |
| S2 | Confirmed no change needed | `neoHGT/gtdb.py` (`write_taxdump()`), `neoHGT/tests/test_gtdb_core.py::TestTaxdumpConversion` | The converter produces output exactly in the `nodes.dmp`/`names.dmp` same-directory, NCBI-delimiter, 13-column shape (the rank vocabulary was also specifically run through DIAMOND's whitelist), so `search.py`'s check that "the two dmp files must be in the same directory" holds as-is and was not changed |
| S3 | Adapted for GTDB | `neoHGT/search.py` (`set_params_gtdb()` forces `fetch_enable = 'no'`), `neoHGT/tests/test_gtdb_pipeline.py::TestNoNetwork` | The plan's worry — "assuming any taxId can be looked up at NCBI EFetch" — is a real risk: querying a synthetic taxId at NCBI either returns nothing or **returns the wrong real organism**. In GTDB mode automatic supplemental fetching is force-disabled and the reason is printed; UT-3.6/DoD-4 assert zero network calls with mocks |
| S4 | Adapted for GTDB | `neoHGT/search.py` (`update_dmp_guard()`), `neoHGT/tests/test_gtdb_pipeline.py::TestPlanVerificationIds` | The runtime supplemental-fetch path of "appending to the taxdump while searching" is rejected in GTDB mode (`update_dmp_files` must not be called), because the taxonomy must stay consistent with the `taxid.map.gz` written to disk at build time, or `scores.tsv` cannot be explained afterwards; the corresponding test replaces that function with a stub that raises |
| S5 | Adapted for GTDB | `neoHGT/search.py` (`parse_def_table_checked()`), `neoHGT/tests/test_gtdb_pipeline.py::TestHitTableValidation` | Upstream swallows the IndexError of `x[6]` into a one-line print, so hits **silently disappear**; the GTDB face validates the column count explicitly (6 columns → a clear error plus the `--outfmt` string that should be used), and adds `--require-staxids` |
| S6 | Adjusted, differs from plan | `neoHGT/database.py` (`compile_gtdb_database()`), `neoHGT/tests/test_gtdb_pipeline.py::TestBlastRoute`, `缓存测试/reports/P0_spike_S2_blast.md` | The plan said "feasible, as long as `makeblastdb -taxid_map`". **Measurement says otherwise**: `makeblastdb` cannot index GTDB sequence IDs of 151–270 characters (the S2 spike), so `--compile blast/both` is refused outright on the GTDB face with an explanation of the 50-character cap, and only an explicit `--force` proceeds; the supported route is DIAMOND (ADR-8). The BLAST route itself remains fully usable (`--compile blast --force` + `-taxid_map` with GTDB-derived taxIds) |
| S7 | Adapted for GTDB | `neoHGT/gtdb.py` (`parse_seqid()`'s `contig`/`contig_nover` two forms), `neoHGT/tests/test_gtdb_core.py::TestSeqIdCoverageProperty` | Upstream's "strip the `.version`" regex is wrong for IDs like `…_SerialNumber_AE017199.1_1` (the last segment is not a version), so route C's parsing keeps both the with-version and version-stripped forms, and the official form additionally has `prefix2genome` as a fallback; a property test asserts parsability for all real IDs |
| S8 | Confirmed no change needed | `neoHGT/search.py` (`OUTFMT_HINT`), `缓存测试/reports/RT1_ncbi_fixtures.json` | The search core is source-agnostic and unchanged except for the hit-table column conventions: RT-1 records that `search.py` to this day has `deleted_lines=0` (only 372 lines appended), and the NCBI-side products are byte-for-byte unchanged |
| S9 | Confirmed no change needed | `neoHGT/tests/test_gtdb_pipeline.py::TestFixtureEndToEnd` | Self-alignment is source-agnostic; the newly added offline end-to-end test asserts that `# Score:` is computed from a real self-alignment (not passed in from precomp), showing that this path works normally on a GTDB database |
