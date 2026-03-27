[PLANS]
- 2026-03-27T21:40:00+01:00 [USER] Create a new repository to compare Ensembl VEP results against biodatageeks/vepyr using Python DataFrame-based diffs, with container-first execution and staged support for cache/flag combinations.

[DECISIONS]
- 2026-03-27T21:40:00+01:00 [USER] Initialize a new git repository, but perform no other git mutations.
- 2026-03-27T21:40:00+01:00 [USER] v1 scope is run-and-compare, not compare-only.
- 2026-03-27T21:40:00+01:00 [USER] Canonical comparison is two-tier: variant/allele level and consequence-level drill-down.
- 2026-03-27T21:40:00+01:00 [USER] Default execution mode is container-first.
- 2026-03-27T21:40:00+01:00 [USER] v1 preset support starts with `--everything` plus Ensembl cache; architecture must stage future `merged` and `refseq`.
- 2026-03-27T21:40:00+01:00 [USER] Canonical equality is parsed semantic CSQ comparison rather than raw-string equality.
- 2026-03-27T21:40:00+01:00 [USER] Default `vepyr` source is a local checkout path.

[PROGRESS]
- 2026-03-27T21:40:00+01:00 [TOOL] Verified the workspace started empty and initialized git in `/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly`.
- 2026-03-27T22:05:00+01:00 [CODE] Added initial repository scaffold: project metadata, local `AGENTS.md`, staged preset manifest, container definitions, Python package, fixture VCFs, and unit/integration tests.
- 2026-03-27T22:05:00+01:00 [CODE] Implemented CLI commands `run`, `list-presets`, and `inspect-run`, plus semantic VCF normalization, tiered comparison, reporting, and runtime artifact handling.
- 2026-03-27T22:20:00+01:00 [CODE] Added root `requirements.txt`, documented local environment setup in `README.md`, created `.venv`, and expanded `.gitignore` to cover virtualenv, coverage, wheel metadata, IDE files, and local override artifacts.
- 2026-03-27T22:50:00+01:00 [CODE] Added repo-local `.env` configuration support with `.env.example`, automatic CLI loading, and README documentation for path-based runtime settings.
- 2026-03-27T23:00:00+01:00 [CODE] Cloned `biodatageeks/vepyr` into repo root under `vepyr/`, added it to `.gitignore`, and pointed local `.env` plus `.env.example` at that checkout path.
- 2026-03-27T22:49:52+0100 [CODE] Fixed local `vepyr` execution by preserving the configured `.../.vepyr/bin/python` symlink path instead of resolving it to the system interpreter, then aligned normalized count columns to `Int64` for `dataframely` schema validation.
- 2026-03-27T22:49:52+0100 [CODE] Added regression tests for preserved `vepyr_python` interpreter paths and normalized integer count dtypes; updated `README.md` with the verified local execution workflow and smoke-run result.

[DISCOVERIES]
- 2026-03-27T21:40:00+01:00 [TOOL] Existing `annotator_testing` runner compares VCF `INFO.CSQ` semantically and already contains `vepyr` adapter logic worth mirroring conceptually.
- 2026-03-27T21:40:00+01:00 [TOOL] `diffly` requires Python >=3.11 and compares Polars DataFrames using unique primary keys.
- 2026-03-27T21:40:00+01:00 [TOOL] `polars-bio` exposes `scan_vcf` / `read_vcf` and supports lazy VCF ingestion.
- 2026-03-27T21:40:00+01:00 [TOOL] Golden input exists at `~/Downloads/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf` and is about 2.6G.
- 2026-03-27T22:05:00+01:00 [TOOL] The local host Python lacked `pytest`, `polars`, `rich`, `diffly`, `dataframely`, and `polars_bio`, so runtime validation was limited to static compilation and CLI smoke tests that avoid those imports.
- 2026-03-27T22:35:00+01:00 [TOOL] `polars-bio` exposes VCF positions as `start`/`end` rather than `pos`; the adapter now maps `start + 1` back to VCF-style 1-based `POS`.
- 2026-03-27T22:35:00+01:00 [TOOL] Minimal fixture VCFs can fail in `polars-bio` / DataFusion with `VCF read error: unexpected EOL`, so normalization now falls back to a simple text parser for those cases.
- 2026-03-27T22:35:00+01:00 [TOOL] Docker is available locally, but no local `vepyr` checkout was found under `/Users/lukaszjezapkowicz/Desktop/magisterka/praca`, and no concrete VEP cache path has been resolved yet.
- 2026-03-27T22:50:00+01:00 [TOOL] `terraform.tfvars` in `annotator_testing` exposes `vep_cache_dir=/Users/lukaszjezapkowicz/.vep`, the GRCh38 FASTA under that cache, and the local `vep_bin` path; these values were mirrored into `.env.example` and local `.env`.
- 2026-03-27T23:00:00+01:00 [TOOL] Root-local `vepyr/` checkout now exists and contains `pyproject.toml`; `git status` confirms the directory is ignored.
- 2026-03-27T22:49:52+0100 [TOOL] Resolving the configured `VEPYR_DIFFLY_VEPYR_PYTHON` symlink path breaks the `vepyr` runtime on macOS because `/.../.vepyr/bin/python` resolves to `/Library/Frameworks/.../python3.12`, which loses the virtualenv site-packages and caused `ModuleNotFoundError: No module named 'vepyr'`.
- 2026-03-27T22:49:52+0100 [TOOL] Real smoke-run normalization produced unsigned aggregate counts from Polars (`UInt32`), while the contract schema expected signed `Int64`; explicit casting in normalization resolves that mismatch cleanly.
- 2026-03-27T22:49:52+0100 [TOOL] Existing external assets are sufficient for local end-to-end execution: VEP binary and cache under `/Users/lukaszjezapkowicz/.vep`, `vepyr` interpreter under `annotator_testing/runner/.vepyr/bin/python`, and parquet cache root under `annotator_testing/.cache/cache_testing/vepyr_cache`.

[OUTCOMES]
- 2026-03-27T21:40:00+01:00 [ASSUMPTION] UNCONFIRMED until verification: the first implementation will favor deterministic local unit tests and provide runtime containers/configuration without executing real VEP inside this session.
- 2026-03-27T22:05:00+01:00 [TOOL] Verified `PYTHONPATH=src python3 -m compileall src tests`, `PYTHONPATH=src python3 -m vepyr_diffly.cli list-presets`, and `PYTHONPATH=src python3 -m vepyr_diffly.cli inspect-run --run-dir /tmp/vepyr-diffly-smoke`.
- 2026-03-27T22:05:00+01:00 [TOOL] UNCONFIRMED: `pytest` suite and end-to-end containerized `run` command were not executed in this session because required Python/runtime dependencies are not installed locally.
- 2026-03-27T22:20:00+01:00 [TOOL] Installed local Python dependencies into `.venv` from `requirements.txt` and verified `.venv/bin/python -m vepyr_diffly.cli list-presets`.
- 2026-03-27T22:35:00+01:00 [TOOL] Verified `.venv/bin/pytest -q` with `7 passed`.
- 2026-03-27T22:35:00+01:00 [TOOL] UNCONFIRMED: full `vepyr-diff run ...` against real VEP remains blocked on user-provided or discoverable `vepyr` checkout path and VEP cache location.
- 2026-03-27T22:50:00+01:00 [TOOL] Verified `.env` loading via `vepyr_diffly.settings.load_repo_env()`, `.venv/bin/python -m vepyr_diffly.cli list-presets`, and `.venv/bin/pytest -q` after adding `python-dotenv`.
- 2026-03-27T23:00:00+01:00 [TOOL] Verified local `.env` resolves `VEPYR_DIFFLY_VEPYR_PATH=/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/vepyr` and `PYTHONPATH=src .venv/bin/python -m vepyr_diffly.cli list-presets` still works.
- 2026-03-27T22:49:52+0100 [TOOL] Verified `PYTHONPATH=src .venv/bin/pytest -q` with `8 passed` after adding regression coverage for the local runtime fixes.
- 2026-03-27T22:49:52+0100 [TOOL] Verified end-to-end local smoke run with `PYTHONPATH=src .venv/bin/python -m vepyr_diffly.cli run --output-dir runs/local-smoke-fix2`: variant tier `1000` joined equal, consequence tier `34741` joined equal, and zero mismatches in both tiers.
- 2026-03-27T22:49:52+0100 [TOOL] Verified `PYTHONPATH=src .venv/bin/python -m vepyr_diffly.cli inspect-run --run-dir runs/local-smoke-fix2` and confirmed `summary.json`, parquet diffs, TSVs, and runtime logs were written to the run directory.
