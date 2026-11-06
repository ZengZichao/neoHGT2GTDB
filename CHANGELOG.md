# Changelog

All notable changes to neoHGT2GTDB are recorded here. The project follows a `base+local` version scheme: the upstream base version is preserved and a `+gtdbN` local segment marks this adaptation.

## 2.0b3+gtdb1

Initial public release of the GTDB adaptation.

### Added

- GTDB database build from a local, read-only GTDB release, with automatic layout detection (`local-renamed`, `official-plain`, `official-gz`) and dataset self-checks.
- NCBI-style taxdump conversion with synthetic, append-only taxIds and cross-release ID stability through `gtdb/taxid.map.gz`, `gtdb/delnodes.txt` and `gtdb/reserved_ids.txt`.
- `neoHGT gtdb probe`, `reconcile`, `taxdump`, `map-report`, `gtdbtk` and `tree-dist` sub-commands.
- `neoHGT database --source gtdb|hybrid`, including representative selection, taxon subsets, genome lists, rank sampling, quality filters, member expansion and resumable extraction.
- GTDB ranks, official donor labels, placeholder-name exemptions, `--donor-ncbi-map` and `--taxon-synonyms` in search and analysis.
- `--source hybrid`: grafting the top of the NCBI tree onto the GTDB tree with a configurable LCA policy.
- Offline GTDB-Tk classification conversion and phylogeny-assisted close groups from per-rank trees.
- `--stream-extraction` on the NCBI route: disk-backed protein extraction with byte-identical outputs.
- `--taxmap` disk-anchored index for protein-to-taxId maps at full scale.
- English and Chinese user guides, a documentation index, a complete CLI reference, a configuration reference and a troubleshooting guide.
- Machine-readable manifests for added functions, plan IDs and configuration-route equivalence.

### Known limitations

- A full build from the official flat GTDB layout has not been tested end to end; a metadata table is still required for database builds.
- Hybrid libraries search GTDB proteins only; NCBI protein sequences are not ingested.
- The remote sequence channels only export a fetch plan; there is no download executor.
- The BLAST route is not usable with GTDB IDs because of the 50-character local-ID limit.
- Genus-sample and full-set resource numbers are extrapolated, not measured.
- A head-to-head comparison with the NCBI route has not been completed.
