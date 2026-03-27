# vepyr-diffly

`vepyr-diffly` compares Ensembl VEP output against `vepyr` output using semantic, DataFrame-based diffs.

It does four things:

- run Ensembl VEP and `vepyr` on the same input VCF,
- normalize both annotated VCF outputs into comparable tabular forms,
- compare them semantically with Polars DataFrames,
- print a clear console summary and write detailed diff artifacts for fixing `vepyr`.

## Current Scope

The current implementation is no longer only scaffolded. It has been verified locally in `execution_mode=local` on the first `1000` variants from the golden HG002 VCF.

Current verified scope:

- one active preset: `ensembl_everything`
- semantic comparison of annotated VCF `INFO/CSQ`
- two comparison tiers:
  - variant / allele summary
  - consequence-level normalized rows
- local runtime execution against a real local VEP installation
- local runtime execution against an existing `vepyr` Python environment and parquet cache
- clear console and file-based reporting

The architecture is prepared for later `merged`, `refseq`, and additional flag combinations.

## Python Environment

The root [requirements.txt](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/requirements.txt) is the simplest way to prepare a local Python environment.

It installs:

- the local package itself in editable mode,
- runtime dependencies declared in [pyproject.toml](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/pyproject.toml),
- developer tools currently used by this repo: `pytest` and `ruff`.

Recommended setup:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
cp .env.example .env
```

Quick verification after installation:

```bash
python -m vepyr_diffly.cli list-presets
pytest
```

Notes:

- Python `>=3.11` is required.
- The Python environment alone is enough for local normalization, reporting, and tests.
- The verified runtime path today is local execution, not Docker.

## Repo Configuration

Runtime paths should live in repo-local `.env`, not be hardcoded in commands.

Files:

- [`.env.example`](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/.env.example)
  - committed template with the expected variables
- `.env`
  - local working copy, ignored by git

The CLI loads `.env` automatically on startup.

Currently supported variables:

- `VEPYR_DIFFLY_PRESET`
- `VEPYR_DIFFLY_EXECUTION_MODE`
- `VEPYR_DIFFLY_INPUT_VCF`
- `VEPYR_DIFFLY_OUTPUT_DIR`
- `VEPYR_DIFFLY_SAMPLE_FIRST_N`
- `VEPYR_DIFFLY_VEP_CACHE_DIR`
- `VEPYR_DIFFLY_VEP_CACHE_VERSION`
- `VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR`
- `VEPYR_DIFFLY_REFERENCE_FASTA`
- `VEPYR_DIFFLY_VEPYR_PATH`
- `VEPYR_DIFFLY_VEPYR_PYTHON`
- `VEPYR_DIFFLY_VEP_BIN`
- `VEPYR_DIFFLY_VEP_PERL5LIB`
- `VEPYR_DIFFLY_VEP_MERGED`

The example file already contains the VEP cache and FASTA paths discovered from:

- `/Users/lukaszjezapkowicz/Desktop/magisterka/praca/gdl-annotations-infra/modules/python/annotator_testing/terraform.tfvars`

For the working local setup used in this repo, the important variables are:

- `VEPYR_DIFFLY_EXECUTION_MODE=local`
- `VEPYR_DIFFLY_VEP_CACHE_DIR`
- `VEPYR_DIFFLY_VEP_CACHE_VERSION`
- `VEPYR_DIFFLY_REFERENCE_FASTA`
- `VEPYR_DIFFLY_VEP_BIN`
- `VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR`
- `VEPYR_DIFFLY_VEPYR_PYTHON`

If you want to reuse the already prepared `annotator_testing` `vepyr` environment, point:

- `VEPYR_DIFFLY_VEPYR_PYTHON` at `.../runner/.vepyr/bin/python`
- `VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR` at the parent directory that contains `parquet/115_GRCh38_vep`

Then you can run:

```bash
python -m vepyr_diffly.cli run
```

Verified smoke command:

```bash
python -m vepyr_diffly.cli run --output-dir runs/local-smoke
```

Verified result on `2026-03-27T22:49:00+01:00`:

- sample size: `1000`
- variant tier: `1000` joined equal, `0` mismatches
- consequence tier: `34741` joined equal, `0` mismatches

## CLI

```bash
python -m vepyr_diffly.cli list-presets
python -m vepyr_diffly.cli run --preset ensembl_everything --input-vcf /path/to/input.vcf --output-dir runs/demo --sample-first-n 1000 --vepyr-path /path/to/vepyr --vep-cache-dir /path/to/vep/cache
python -m vepyr_diffly.cli inspect-run --run-dir runs/demo
```

## How To Run

Typical local session:

```bash
cd /Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli list-presets
PYTHONPATH=src pytest -q
PYTHONPATH=src python -m vepyr_diffly.cli run
PYTHONPATH=src python -m vepyr_diffly.cli inspect-run --run-dir runs/local-smoke
```

If `.env` is configured, `run` uses it automatically.

Useful commands:

- list presets:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli list-presets
```

- run tests:

```bash
source .venv/bin/activate
PYTHONPATH=src pytest -q
```

- smoke run on the configured sample:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli run --output-dir runs/local-smoke
```

- smoke run with a custom sample size:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli run \
  --output-dir runs/local-smoke-500 \
  --sample-first-n 500
```

- run against another VCF:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli run \
  --input-vcf /absolute/path/to/input.vcf \
  --output-dir runs/other-input \
  --sample-first-n 1000
```

- run on the whole input VCF while `.env` still contains sampling:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli run \
  --output-dir runs/full-golden \
  --sample-first-n 1000
```

This is still a sampled run. For a true full run without sampling, remove `VEPYR_DIFFLY_SAMPLE_FIRST_N` from `.env` first, then run:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli run --output-dir runs/full-golden
```

This matters because if `VEPYR_DIFFLY_SAMPLE_FIRST_N` is still present in `.env`, it remains active unless you edit the config.

## Expected Output

The console output is intentionally short and should be readable during iterative fixing.

During a real run you should expect:

- a progress bar for VCF row ingestion,
- a short run header with preset, input path, and sample size,
- one summary table with one row for `variant` and one for `consequence`.

Typical successful smoke output:

```text
Preset: ensembl_everything
Input: /Users/lukaszjezapkowicz/Downloads/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf
Sample first N: 1000

Comparison Summary

Tier         Equal  Left only  Right only  Unequal  Joined equal
variant      yes    0          0           0        1000
consequence  yes    0          0           0        34741
```

Meaning of the summary columns:

- `Equal`
  - whether the whole tier is equal after normalization and diffing
- `Left only`
  - rows present only in Ensembl VEP output
- `Right only`
  - rows present only in `vepyr` output
- `Unequal`
  - rows matched by primary key but with different payload values
- `Joined equal`
  - rows matched by primary key and identical after normalization

What to expect in the run directory:

- `summary.json`
  - machine-readable summary of both tiers
- `summary.md`
  - human-readable summary for quick inspection
- `variant_mismatches.tsv`
  - practical file to inspect variant-level mismatches
- `consequence_mismatches.tsv`
  - practical file to inspect consequence-level mismatches
- `variant_diff.parquet`
  - diff rows in columnar form for further analysis
- `consequence_diff.parquet`
  - consequence diff rows in columnar form
- `runtime/effective_config.json`
  - exact resolved configuration used for the run
- `runtime/vep.log`
  - full VEP invocation and stderr/stdout
- `runtime/vepyr.log`
  - full `vepyr` invocation and stderr/stdout

If the run is equal:

- both mismatch TSVs should be empty or contain no diff rows,
- both parquet diff files should contain no mismatch rows,
- `summary.json` should report `equal: true` for both tiers.

If the run is not equal:

- start with `summary.md`,
- then inspect `variant_mismatches.tsv` and `consequence_mismatches.tsv`,
- use the runtime logs only if the failure happened before diff generation.

## Troubleshooting

If `pytest` fails:

- first run it again with verbose output:

```bash
source .venv/bin/activate
PYTHONPATH=src pytest -q -vv
```

- most failures in this repo should be deterministic and local to normalization, comparison, or report generation

If `run` fails before producing annotated VCFs:

- inspect `runtime/vep.log`
- inspect `runtime/vepyr.log`
- inspect `runtime/effective_config.json`

Typical causes:

- bad `VEPYR_DIFFLY_VEP_BIN`
- bad `VEPYR_DIFFLY_VEP_CACHE_DIR`
- bad `VEPYR_DIFFLY_REFERENCE_FASTA`
- bad `VEPYR_DIFFLY_VEP_PERL5LIB`
- `VEPYR_DIFFLY_VEPYR_PYTHON` does not have `vepyr` installed
- `VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR` points to the wrong directory

If `vepyr` fails to import:

- check that `VEPYR_DIFFLY_VEPYR_PYTHON` points to a virtualenv interpreter like `.../.vepyr/bin/python`
- do not replace it with a resolved system interpreter path

If the run completes but reports mismatches:

- `Left only > 0` means VEP produced rows that `vepyr` did not
- `Right only > 0` means `vepyr` produced rows that VEP did not
- `Unequal > 0` means both sides produced the same key but different values

For fix work, the most important files are:

- `summary.md`
- `variant_mismatches.tsv`
- `consequence_mismatches.tsv`

If you want to rerun cleanly:

- choose a new `--output-dir`
- or delete the old run directory manually before rerunning

The simplest safe pattern is:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli run --output-dir runs/local-smoke-next
```

## What The Repo Does

The repo currently exposes three CLI commands:

- `list-presets`
  - shows the staged preset matrix
- `run`
  - resolves one preset,
  - optionally samples the first `N` variants from the source VCF,
  - runs VEP and `vepyr`,
  - normalizes both annotated outputs,
  - compares them at two levels:
    - variant / allele level,
    - consequence / CSQ-entry level,
  - prints a summary and writes artifacts
- `inspect-run`
  - reads an existing `summary.json` and prints it back in a quick inspection form

The canonical comparison is semantic:

- row order does not matter,
- CSQ entry order does not matter after normalization,
- the comparison is based on parsed fields, not raw VCF line equality.

## Outputs

Each run directory contains:

- `summary.md`
- `summary.json`
- `variant_diff.parquet`
- `consequence_diff.parquet`
- `variant_mismatches.tsv`
- `consequence_mismatches.tsv`
- `runtime/effective_config.json`
- `runtime/*.log`

## Runtime Notes

- Canonical equality is semantic comparison of parsed `CSQ`.
- Raw VCF equivalence is out of scope.
- The current verified execution path is local runtime via `VEPYR_DIFFLY_EXECUTION_MODE=local`.
- `VEPYR_DIFFLY_VEPYR_PYTHON` should point at a Python interpreter that already has `vepyr` importable.
- `VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR` should point at the cache root, not directly at the leaf parquet feature directory.
- The default golden input is expected at `~/Downloads/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf`.
- The local smoke run now works end-to-end with the checked-in code and external assets already present on this machine.
- Docker/container adapters remain in the repo, but they are not the currently verified execution path.
