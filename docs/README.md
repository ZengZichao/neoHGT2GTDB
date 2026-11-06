# Documentation index

This directory contains the user manuals and the machine-readable maintenance manifests for neoHGT2GTDB.

## User documentation

| File | Language | Audience | Contents |
| --- | --- | --- | --- |
| [`../README.md`](../README.md) | English | First-time users | Installation, GTDB quick start, NCBI route, search and analysis |
| [`../README.zh-CN.md`](../README.zh-CN.md) | Chinese | First-time users | Chinese mirror of the README |
| [`GTDB-数据库.md`](GTDB-数据库.md) | Chinese | GTDB users | Complete GTDB guide: data preparation, probe, build, search, analyze, hybrid mode, budgets, limitations, troubleshooting |
| [`en/GTDB-database.md`](en/GTDB-database.md) | English | GTDB users | English translation of the complete GTDB guide |
| [`cli-reference.md`](cli-reference.md) | English | All users | Every command-line option, grouped by sub-command |
| [`config-reference.md`](config-reference.md) | English | All users | `config.yml` keys and defaults |
| [`troubleshooting.md`](troubleshooting.md) | English | All users | Common errors and fixes |
| [`../CHANGELOG.md`](../CHANGELOG.md) | English | All users | Release notes |
| [`../CONTRIBUTING.md`](../CONTRIBUTING.md) | English | Contributors | Development setup, tests, documentation gates and manifest maintenance |
| [`../CITATION.cff`](../CITATION.cff) | — | All users | Machine-readable citation metadata |

The GTDB guide and the README are the primary sources for users. The CLI and configuration references are maintained by hand against the shipped code and configuration; the traceability gate (`tools/verify_plan_traceability.py`) checks their documented options, option values and default promises.

## Machine-readable manifests

The three JSON files below are used by the project's CI and maintenance tools. They keep Chinese file names and have no English counterparts. They are not required to run neoHGT2GTDB, but changes to the CLI, configuration or plan inventory must keep them in sync.

### `新增函数清单.json`

- **Purpose**: freeze the functions/methods that this adaptation added, plus the upstream names they are compared against. CI uses it to detect functions that no CLI path can reach and functions that are neither declared new nor known upstream.
- **Consumer**: `tools/verify_production_reachability.py --added docs/新增函数清单.json` (CI step "Production reachability").
- **Schema**:

| Field | Type | Meaning |
| --- | --- | --- |
| `schema_version` | integer | Manifest schema version |
| `generated_utc` | string | Generation timestamp (UTC) |
| `generated_by` | string | Tool that produced the file |
| `generated_from` | string | Source tree the scan ran on |
| `baseline` | string | Upstream baseline used for the "added" comparison |
| `functions` | array of strings | Added function/method names |
| `upstream` | array of strings | Names inherited from upstream |
| `counts` | object | `added` and `upstream` array lengths |

- **Regenerate**:
  ```sh
  python tools/verify_production_reachability.py --code . \
    --emit-added docs/新增函数清单.json \
    --baseline /path/to/pristine/neoHGT-main
  ```
  The baseline is optional for the frozen list; without it the upstream half is not recomputed.

### `方案编号清单.json`

- **Purpose**: a frozen digest of every numbered list in the design plan (tasks, defects, deliverables, risks, acceptance criteria, decisions, verification IDs, options and symbols). It lets the traceability gate check the plan-derived IDs without the plan document in the repository.
- **Consumer**: `tools/verify_plan_traceability.py --code .` (CI step "Register gate").
- **Schema**:

| Field | Type | Meaning |
| --- | --- | --- |
| `schema_version` | integer | Manifest schema version |
| `generated_utc` | string | Generation timestamp (UTC) |
| `generated_by` | string | Tool that produced the file |
| `source` | string | Name of the source plan document |
| `counts` | object | Number of entries per inventory category |
| `inventory` | object | The ID lists themselves, keyed by category |

- **Regenerate** (requires the plan document):
  ```sh
  python tools/verify_plan_traceability.py --plan /path/to/plan.md \
    --emit-ids docs/方案编号清单.json
  ```
  The plan is not part of this repository, so the committed digest is the public reference. The gate checks that every ID it consumes is present and that the counts match the arrays.

### `路由等价性清单.json`

- **Purpose**: record every `config.yml` key read by the stage modules, whether the CLI route and the configuration route hand over the same object, and — for divergent keys — the function that normalises the value and the reason it is safe.
- **Consumer**: `tools/verify_plan_traceability.py --code .` (route-equivalence rule).
- **Schema**:

| Field | Type | Meaning |
| --- | --- | --- |
| `schema_version` | integer | Manifest schema version |
| `generated_utc` | string | Generation timestamp (UTC) |
| `generated_by` | string | Tool that produced the file |
| `generated_from` | string | Source tree the scan ran on |
| `keys` | object | One record per configuration key |
| `note` | string | Explanation of the divergence rule |
| `unresolved_call_sites` | array | Call sites the scanner could not resolve (must be empty) |

Each record under `keys` contains `attr`, `modules`, `converter`, `cli_type`, `yaml_type`, `in_config` and `divergent`. Divergent records also contain `neutralizer` and `reason`.

- **Regenerate**:
  ```sh
  python tools/verify_plan_traceability.py --code . \
    --emit-route docs/路由等价性清单.json
  ```

### Updating the manifests

Update the manifests when:

- a function or method is added or removed (`新增函数清单.json`);
- a config key is added, removed or changes type (`路由等价性清单.json`);
- a plan-derived ID list changes (`方案编号清单.json`).

Run the corresponding `--emit-*` command, then run the full traceability gate and the reachability tool. A manifest that no longer matches the code is reported as a gate failure rather than silently ignored.

## Public exception lists

Two small JSON files make deliberate exceptions auditable:

### `test-exceptions.json`

Lists pre-existing upstream tests that exercise code but contain no assertion. They are kept unchanged under the append-only compatibility policy. New tests must contain an observable assertion, and the traceability gate reports any assertion-free test that is not listed here.

### `reachability-exceptions.json`

Lists public helper functions that are kept for library users and tests but are deliberately not reachable from a CLI path. The production-reachability tool reports these as declared exceptions instead of undeclared unreachable functions.

Both files are checked by CI. Adding an entry is a deliberate statement that the exception is part of the public interface, not a way to silence a new defect.

## What the public CI verifies

| Check | Where |
| --- | --- |
| Relative links, in-repository anchors, English-document language, shell variables in examples, internal-reference leftovers | `tools/verify_docs.py` |
| Frozen plan IDs, CLI option surface, documented example options and values | `tools/verify_plan_traceability.py` |
| Documented default values against the delivered code | `tools/verify_plan_traceability.py` |
| Configuration-route equivalence (73 keys) | `tools/verify_plan_traceability.py` |
| Added-function reachability and public exception list | `tools/verify_production_reachability.py` |
| GTDB fixture provenance | `tools/verify_fixture_provenance.py` |
| NCBI regression against the recorded baseline | `pytest` + `ci/rt1_baseline_failures.txt` |
| External URL status (advisory; hosts may rate-limit) | `tools/check_external_links.py` |
| Read-only GTDB input guard | CI workspace check |

Some checks require private maintenance records and the author's full datasets, which are not part of this checkout. The gate reports those checks as skipped instead of counting them as passed. The public documentation therefore states what is checked here and what is not; it does not claim that every maintenance record is verified by public CI.
