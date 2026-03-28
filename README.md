# vepyr-diffly

`vepyr-diffly` compares Ensembl VEP output against `vepyr` output using semantic, DataFrame-based diffs.

It does four things:

- run Ensembl VEP and `vepyr` on the same input VCF,
- normalize both annotated VCF outputs into comparable tabular forms,
- compare them semantically with Polars DataFrames,
- print a clear console summary and write detailed diff artifacts for fixing `vepyr`.

## Local Quickstart

This is the shortest practical path to test VEP vs `vepyr` locally on this machine.

### 1. Create the Python environment

```bash
cd /Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### 2. Create `.env`

```bash
cp .env.example .env
```

The CLI loads `.env` automatically. The minimum important variables are:

- `VEPYR_DIFFLY_EXECUTION_MODE=local`
- `VEPYR_DIFFLY_INPUT_VCF`
- `VEPYR_DIFFLY_OUTPUT_DIR`
- `VEPYR_DIFFLY_VEP_CACHE_DIR`
- `VEPYR_DIFFLY_VEP_CACHE_VERSION`
- `VEPYR_DIFFLY_REFERENCE_FASTA`
- `VEPYR_DIFFLY_VEP_BIN`
- `VEPYR_DIFFLY_VEPYR_PYTHON`
- `VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR`

The current [.env.example](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/.env.example) already points at the local paths that were verified in this workspace:

- input VCF: `/Users/lukaszjezapkowicz/Downloads/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf`
- VEP cache: `/Users/lukaszjezapkowicz/.vep`
- FASTA: `/Users/lukaszjezapkowicz/.vep/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa`
- VEP binary: `/Users/lukaszjezapkowicz/VepAnnotations/ensembl/ensembl-vep/vep`
- `vepyr` Python env: `/Users/lukaszjezapkowicz/Desktop/magisterka/praca/gdl-annotations-infra/modules/python/annotator_testing/runner/.vepyr/bin/python`
- `vepyr` cache root: `/Users/lukaszjezapkowicz/Desktop/magisterka/praca/gdl-annotations-infra/modules/python/annotator_testing/.cache/cache_testing/vepyr_cache`

### 3. Sanity-check the repo

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli list-presets
PYTHONPATH=src pytest -q
```

### 4. Run a local smoke test

For a full end-to-end smoke test, including annotation by both engines:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli run \
  --output-dir runs/local-smoke-1000 \
  --sample-first-n 1000
```

To change how many input variants are used in the smoke test, change `--sample-first-n`:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli run \
  --output-dir runs/local-smoke-5000 \
  --sample-first-n 5000
```

### 5. Reuse existing annotated VCFs without re-annotating

This is the faster iteration path when VEP output already exists and you only want to compare:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli compare-existing \
  --preset ensembl_everything \
  --left-vcf runs/local-smoke-fix2/runtime/vep.annotated.vcf \
  --right-vcf runs/local-smoke-fix2/runtime/vepyr.annotated.vcf \
  --output-dir runs/compare-only
```

### 6. What to expect on stdout

`run` has two expensive phases:

- annotation
- comparison

The compare phase prints explicit progress and ends with a summary table like:

```text
Comparison Summary
variant      yes   0   0   0   1000
consequence  yes   0   0   0   34741
```

The output directory then contains:

- `summary.json`
- `summary.md`
- `variant_diff.parquet`
- `consequence_diff.parquet`
- `variant_mismatches.tsv`
- `consequence_mismatches.tsv`
- `runtime/compare.progress.log`
- `runtime/vep.log`
- `runtime/vepyr.log`

## How The Comparison Works

The comparison is semantic, not byte-for-byte.

### Short Technical Summary

In short, the repo does this:

1. prepares one canonical input VCF for both annotators,
2. runs VEP and `vepyr` on exactly that same prepared input,
3. reads both annotated VCFs into tabular form,
4. normalizes `INFO/CSQ` into stable variant and consequence tables,
5. compares those tables with `diffly`,
6. writes a compact console summary plus detailed diff artifacts.

The main libraries used by the compare path are:

- `polars-bio` to scan annotated VCF files,
- `polars` to normalize and aggregate rows,
- `diffly` to compare normalized DataFrames by primary key,
- `dataframely` to validate the normalized schema,
- `rich` to print readable summaries and progress.

```text
input VCF
   |
   v
sample first N (optional)
   |
   v
split multi-allelic -> runtime/prepared_input.vcf
   |
   +--------------------+
   |                    |
   v                    v
 VEP                vepyr
   |                    |
   v                    v
vep.annotated.vcf   vepyr.annotated.vcf
   |                    |
   +----------+---------+
              |
              v
   normalize to variant tier + consequence tier
              |
              v
      diffly compare by primary key
              |
              v
 summary.json / summary.md / parquet diffs / mismatch TSV / progress logs
```

For small files the consequence compare can run directly on normalized tables. For large files it switches to a bucketized path: consequence rows are partitioned into hash buckets, each bucket is compared independently with `diffly`, and the bucket results are merged at the end. This keeps memory bounded and lets the compare use multiple CPU cores.

The pipeline is:

1. prepare the input VCF
2. annotate the same prepared input with VEP and `vepyr`
3. normalize both annotated VCFs into comparable tables
4. compare those tables with `diffly`
5. write console summary plus file artifacts

### Input preparation

Before annotation, `run` prepares a canonical input:

- if sampling is enabled, it writes `runtime/sampled_input.vcf`
- it then decomposes multi-allelic rows into single-alt rows
- it writes the final annotation input to `runtime/prepared_input.vcf`
- it records preparation stats in `runtime/input_preparation.json`

This step is required because multi-allelic input created misleading consequence-level diffs, especially for indels. After splitting to single-alt before annotation, the smoke `5000` run became fully equal again.

### Two comparison tiers

The repo compares two normalized tables.

Variant tier:

- primary key: `chrom`, `pos`, `ref`, `alt`
- tracks counts and basic summary columns per variant/allele
- useful to detect gross record-level divergence

Consequence tier:

- one normalized row per semantic `CSQ` consequence
- compared after parsing and normalizing `INFO/CSQ`
- this is the tier used to diagnose real annotation drift

### Why `diffly`

`diffly` is the final diff engine. We do not compare raw VCF text directly. Instead:

- VCF is normalized into tabular rows
- rows get stable primary keys
- duplicates are counted explicitly
- then `diffly` compares the normalized DataFrames

### Why bucketization exists

The full golden annotated VCFs are too large for a naive eager compare.

For large runs, consequence comparison works in buckets:

- normalize consequence rows into hash buckets
- compare bucket-by-bucket with `diffly`
- merge bucket diff artifacts at the end

This keeps memory bounded and gives progress reporting for large compares.

## Results

### A. Smoke Tests

The following smoke matrix was executed locally on `2026-03-28` with fresh per-size runs under [runs/smoke-matrix-full](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/smoke-matrix-full).

Each row means:

- sample `N` source variants from the golden VCF,
- split multi-allelic rows to canonical single-alt prepared input,
- annotate that prepared input with VEP and `vepyr`,
- compare the two annotated outputs with `compare-existing`.

| Sample size | Prepared rows | Split source rows | VEP s | `vepyr` s | Compare s | Variant equal | Consequence equal | Consequence joined equal | Run dir |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | --- |
| `10` | `10` | `0` | `0` | `47` | `0` | `true` | `true` | `432` | [run-10](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/smoke-matrix-full/run-10) |
| `100` | `100` | `0` | `3` | `46` | `1` | `true` | `true` | `8448` | [run-100](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/smoke-matrix-full/run-100) |
| `1000` | `1000` | `0` | `17` | `41` | `3` | `true` | `true` | `34741` | [run-1000](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/smoke-matrix-full/run-1000) |
| `5000` | `5021` | `21` | `50` | `49` | `6` | `true` | `true` | `93694` | [run-5000](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/smoke-matrix-full/run-5000) |
| `10000` | `10066` | `66` | `73` | `53` | `7` | `true` | `true` | `125755` | [run-10000](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/smoke-matrix-full/run-10000) |
| `25000` | `25235` | `235` | `163` | `56` | `14` | `true` | `true` | `257509` | [run-25000](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/smoke-matrix-full/run-25000) |
| `50000` | `50616` | `616` | `351` | `68` | `28` | `true` | `true` | `543416` | [run-50000](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/smoke-matrix-full/run-50000) |
| `100000` | `101350` | `1350` | `720` | `86` | `60` | `true` | `true` | `1099319` | [run-100000](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/smoke-matrix-full/run-100000) |

Observations:

- all smoke runs are fully equal at both tiers
- for larger samples, VEP dominates runtime
- `vepyr` stays materially faster than VEP on every verified smoke point
- multi-allelic decomposition becomes visible from `5000+` and grows with sample size

The raw matrix is also saved as [results.tsv](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/smoke-matrix-full/results.tsv).

### B. Golden Test

This section is intentionally not filled with a final trusted outcome yet.

Current status:

- the old full golden compare-only run proved that the full pipeline scales
- the old full golden consequence result is not yet a trustworthy biological verdict, because those older annotated VCFs were generated before the canonical single-alt preparation flow and the old `vepyr` output showed malformed multi-allelic `CSQ` rows
- the next trustworthy golden result should come from a fresh full run using the same canonical prepared-input path as the smoke matrix above

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

## Full Golden Run

Full golden execution was attempted on `2026-03-27/2026-03-28` with:

```bash
source .venv/bin/activate
VEPYR_DIFFLY_SAMPLE_FIRST_N= PYTHONPATH=src python -m vepyr_diffly.cli run \
  --output-dir runs/full-golden-20260327
```

Observed outcome:

- input VCF: `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf`
- output directory: `runs/full-golden-20260327`
- `vep.annotated.vcf`: about `16G`
- `vepyr.annotated.vcf`: about `15G`
- VEP annotated record count: `4048342`
- `vepyr` annotated record count: `4048342`

What succeeded:

- full Ensembl VEP annotation completed successfully
- full `vepyr` annotation completed successfully
- runtime logs and annotated outputs were preserved

What did not complete:

- the Python comparison stage did not produce `summary.json`, `summary.md`, parquet diff files, or mismatch TSVs for the full golden run
- the visible terminal ending was:

```text
/opt/homebrew/Cellar/python@3.14/3.14.2/Frameworks/Python.framework/Versions/3.14/lib/python3.14/multiprocessing/resource_tracker.py:396: UserWarning: resource_tracker: There appear to be 1 leaked semaphore objects to clean up at shutdown: {'/mp-g_hezd24'}
```

Interpretation at that point:

- the bottleneck was no longer annotation
- the old implementation still had a global consequence merge / compare bottleneck
- the most likely failure point was the previous full-frame consequence normalization strategy

What changed since then:

- `compare-existing` was refactored to use bucketized streaming consequence comparison
- consequence rows are now partitioned into hash buckets and compared bucket-by-bucket with `diffly`
- the smoke workflow for this new path is verified locally
- the exact full golden compare should now be re-run with the new implementation rather than relying on the old failed run

## CLI

```bash
python -m vepyr_diffly.cli list-presets
python -m vepyr_diffly.cli run --preset ensembl_everything --input-vcf /path/to/input.vcf --output-dir runs/demo --sample-first-n 1000 --vepyr-path /path/to/vepyr --vep-cache-dir /path/to/vep/cache
python -m vepyr_diffly.cli compare-existing --preset ensembl_everything --left-vcf /path/to/vep.annotated.vcf --right-vcf /path/to/vepyr.annotated.vcf --output-dir runs/compare-only
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

`run` now prepares a canonical annotation input before invoking VEP and `vepyr`:

- if sampling is enabled, it first writes `runtime/sampled_input.vcf`,
- then it decomposes multi-allelic rows into single-alt rows and writes `runtime/prepared_input.vcf`,
- it also records the preparation stats in `runtime/input_preparation.json`.

This matches the earlier `annotator_testing` behavior for meaningful VEP vs `vepyr` comparison on multi-allelic input.

- compare only on already existing annotated VCFs:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli compare-existing \
  --preset ensembl_everything \
  --left-vcf runs/local-smoke-fix2/runtime/vep.annotated.vcf \
  --right-vcf runs/local-smoke-fix2/runtime/vepyr.annotated.vcf \
  --output-dir runs/compare-only
```

- compare only on the full golden annotated outputs without re-annotating:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli compare-existing \
  --preset ensembl_everything \
  --left-vcf runs/full-golden-20260327/runtime/vep.annotated.vcf \
  --right-vcf runs/full-golden-20260327/runtime/vepyr.annotated.vcf \
  --output-dir runs/full-golden-compare
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
- prepared-input metadata under `runtime/input_preparation.json`,
- timestamped stage logs for normalization and comparison,
- for larger annotated VCFs: explicit consequence bucketization logs and per-bucket compare completion logs,
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
- mismatch TSV files are intentionally sampled on large runs
  - the parquet diff remains complete
  - the TSV is capped to the first `2000` mismatch rows to avoid heavy extra I/O

## Performance Notes

The current compare pipeline is optimized for large annotated VCFs more aggressively than the first implementation:

- left and right consequence bucketization run in parallel,
- left and right variant summaries are materialized in parallel,
- consequence comparison uses dedicated subprocess workers instead of a thread fallback path in the main process,
- consequence rows are hash-partitioned into buckets and compared bucket-by-bucket with `diffly`,
- full mismatch parquet artifacts are preserved, while mismatch TSV output is sampled to reduce avoidable write cost.

For smoke-sized inputs this mainly reduces overhead. For large annotated VCFs the main expected gains are bounded memory use and much better CPU utilization.
- `variant_diff.parquet`
  - diff rows in columnar form for further analysis
- `consequence_diff.parquet`
  - consequence diff rows in columnar form
- `runtime/effective_config.json`
  - exact resolved configuration used for the run
- `runtime/compare.progress.log`
  - timestamped stage log showing how far normalization and diffing have progressed
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
- `runtime/compare.progress.log`

If you want to know how far a long compare-only or full compare has already gone:

- read `runtime/compare.progress.log`
- the most important milestones are:
  - `parsed CSQ header`
  - `materializing variant summary`
  - `variant summary rows=...`
  - `bucketizing consequence rows`
  - `consequence progress ... variants (...%), chunk_parts=..., buckets=...`
  - `comparing ... buckets with diffly`
  - `completed bucket 00xx (x/y)`
  - `merging bucket diff artifacts`
  - `running variant tier diff`
  - `writing summaries`
  - `completed successfully`

For small inputs:

- the compare path chooses fewer buckets automatically to reduce overhead
- the verified smoke compare on `1000` variants currently uses `8` buckets

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

For larger annotated VCFs, consequence comparison now works like this:

- stream base VCF records in chunks,
- vectorize CSQ explode and normalization in Polars,
- hash-partition normalized consequence rows into bucket parquet directories,
- compare corresponding buckets with `diffly`,
- merge per-bucket diff artifacts into final parquet and TSV outputs.

## Outputs

Each run directory contains:

- `summary.md`
- `summary.json`
- `variant_diff.parquet`
- `consequence_diff.parquet`
- `variant_mismatches.tsv`
- `consequence_mismatches.tsv`
- `normalized/left.consequence_buckets/`
- `normalized/right.consequence_buckets/`
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

## Known Limitations

- Only the `ensembl_everything` preset is currently verified end-to-end.
- `merged` and `refseq` cache flavors are part of the intended architecture, but are not implemented and verified yet.
- Docker/container execution is scaffolded in the repo, but local execution is the only verified runtime path today.
- The verified smoke result is based on the first `1000` variants from the golden HG002 VCF.
- A full unsampled golden annotation run completed, and compare-only mode now exists for reusing those annotated VCFs without repeating annotation.
- The new bucketized compare path is verified on smoke data, but the full golden compare should still be re-run and recorded with this refactored implementation.
- The repo currently compares semantic normalized outputs, not byte-for-byte annotated VCF files.
- The current documentation and examples assume the local external assets already exist on this machine: VEP binary, VEP cache, FASTA, `vepyr` Python environment, and `vepyr` parquet cache.
