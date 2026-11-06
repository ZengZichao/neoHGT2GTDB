# Contributing

Contributions to the GTDB adaptation are welcome. This page describes the development environment, the test suites, the documentation gates and the machine-readable manifests.

## Development environment

```sh
conda create -n neoHGT-dev -c conda-forge python=3.11 pyyaml numpy scipy pandas scikit-learn matplotlib pytest
conda activate neoHGT-dev
conda install -c bioconda diamond blast
pip install -e .
```

The GTDB offline tests use only the KB-scale fixtures under `neoHGT/tests/data/GTDB/`; a full GTDB release is not required. Tests that need a release skip themselves when `GTDB_DIR` is not set.

`tools/env.sh` is a convenience wrapper for a workspace layout. It reads optional machine paths from `tools/env.local` (git-ignored). Do not commit absolute paths.

## Tests

Run the fast unit tests:

```sh
python -m pytest neoHGT/tests/test_util.py neoHGT/tests/test_util_gtdb.py neoHGT/tests/test_gtdb_core.py -q
```

Run the GTDB integration tests on the committed fixtures:

```sh
python -m pytest neoHGT/tests/test_gtdb_pipeline.py -q
```

Run the NCBI regression suite:

```sh
python -m pytest neoHGT/tests/test_search.py neoHGT/tests/test_analyze.py neoHGT/tests/test_database.py -q
```

When `GTDB_DIR` points at a full release, the fixture-provenance and read-only guards run against it:

```sh
python tools/verify_fixture_provenance.py --code . --gtdb-dir "$GTDB_DIR"
```

## Documentation gates

The traceability gate checks documented options, option values, default promises, configuration routes and the frozen plan inventory:

```sh
python tools/verify_plan_traceability.py --code . --out /path/to/workspace/gate.json
```

The reachability tool checks that added functions can be reached from a CLI path and that no function is missing from the frozen list:

```sh
python tools/verify_production_reachability.py --code . \
  --added docs/新增函数清单.json --out /path/to/workspace/reachability.json
```

Both tools write their reports inside the configured workspace. They accept a `--code` directory, so they can run from outside the repository.

## Machine-readable manifests

The three manifests under `docs/` are part of the CI contract:

| File | Regenerate with |
| --- | --- |
| `docs/新增函数清单.json` | `python tools/verify_production_reachability.py --code . --emit-added docs/新增函数清单.json --baseline /path/to/upstream` |
| `docs/方案编号清单.json` | `python tools/verify_plan_traceability.py --plan /path/to/plan.md --emit-ids docs/方案编号清单.json` |
| `docs/路由等价性清单.json` | `python tools/verify_plan_traceability.py --code . --emit-route docs/路由等价性清单.json` |

Update a manifest whenever the corresponding code or plan changes. A stale manifest is a gate failure.

## Adding a CLI option

1. Add the option to the module's `arguments` table, with a clear help string, a `choices` list where applicable, and the correct destination attribute.
2. Add a test that passes every documented value to the real parser. The traceability gate rejects a `choices` value that no test has ever parsed.
3. Document the option in `docs/cli-reference.md` and, when it affects users, in the GTDB guide and README.
4. If the option reads a configuration key, add the key to `config.yml` and regenerate `docs/路由等价性清单.json`.
5. Run the traceability gate and the relevant tests.

## Adding a GTDB feature

- Keep the NCBI route unchanged. The regression suite and the append-only audit are there to catch accidental changes.
- Record every new function in `docs/新增函数清单.json` by regenerating it.
- For a new file format or side table, add a provenance line and a reader that skips `#` lines.
- For a new behavior that users must know about, update both language versions of the GTDB guide and the README.
- Add an offline test using the KB-scale fixtures whenever possible.

## Documentation style

- User guides describe what the software does and how to use it. They do not record development rounds, debugging sessions or internal task IDs.
- Commands must be copy-pasteable. Define every environment variable before its first use, or use explicit paths.
- English and Chinese versions of the GTDB guide keep the same section structure and the same numbers.
- `docs/cli-reference.md` and `docs/config-reference.md` are generated from the shipped code and configuration; update the generator or the source table rather than editing the output by hand when possible.
- Relative links must resolve inside the repository. Run `python tools/verify_docs.py .` before opening a pull request.
- External links (NCBI, OneDrive, Dropbox, conda pages) are not checked by CI; run `curl -I` on the URLs in the README and the GTDB guide before a release and replace dead links.

## Submitting changes

- Keep each change focused.
- Include the tests and documentation updates in the same change.
- Describe the user-visible behavior in the pull request.
- Do not commit local paths, dataset files or generated caches.
