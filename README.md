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
- `VEPYR_DIFFLY_VEPYR_USE_FJALL`
- `VEPYR_DIFFLY_COMPARE_MODE`
- `VEPYR_DIFFLY_BUCKET_COUNT`
- `VEPYR_DIFFLY_COMPARE_WORKERS`
- `VEPYR_DIFFLY_MEMORY_BUDGET_MB`
- `VEPYR_DIFFLY_FINGERPRINT_ONLY`

The current [.env.example](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/.env.example) already points at the local paths that were verified in this workspace:

- input VCF: `/Users/lukaszjezapkowicz/Downloads/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf`
- VEP cache: `/Users/lukaszjezapkowicz/.vep`
- FASTA: `/Users/lukaszjezapkowicz/.vep/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa`
- VEP binary: `/Users/lukaszjezapkowicz/VepAnnotations/ensembl/ensembl-vep/vep`
- `vepyr` Python env: `/Users/lukaszjezapkowicz/Desktop/magisterka/praca/gdl-annotations-infra/modules/python/annotator_testing/runner/.vepyr/bin/python`
- `vepyr` cache root: `/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/.cache/vepyr_cache`

If the cache directory also contains `*.fjall` stores, you can enable the faster `vepyr` backend with:

- `VEPYR_DIFFLY_VEPYR_USE_FJALL=true`

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

To force the `vepyr` fjall backend for the run:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli run \
  --output-dir runs/local-smoke-fjall-5000 \
  --sample-first-n 5000 \
  --use-fjall
```

### 5. Reuse existing annotated VCFs without re-annotating

This is the faster iteration path when VEP output already exists and you only want to compare:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli compare-existing \
  --preset ensembl_everything \
  --left-vcf runs/local-smoke-1000/runtime/vep.annotated.vcf \
  --right-vcf runs/local-smoke-1000/runtime/vepyr.annotated.vcf \
  --output-dir runs/compare-only \
  --compare-mode fast \
  --memory-budget-mb 1024
```

Useful compare tuning flags:

- `--compare-mode fast`: exact precheck per bucket, then `diffly` only on buckets that differ
- `--compare-mode debug`: force `diffly` for every bucket and keep full debug artifacts
- `--bucket-count <N>`: override automatic consequence bucket count
- `--compare-workers <N>`: override compare worker process count
- `--memory-budget-mb <N>`: hard memory budget used by the bounded-memory planner for chunk sizes, bucket fanout, and worker count
- `--fingerprint-only`: run the bucket precheck only and skip exact diff artifact generation

### 5a. Run only `vepyr` on an existing prepared input

This is the safest way to produce only the missing `vepyr` output when `prepared_input.vcf` already exists and you do not want to rerun VEP:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli annotate-vepyr \
  --input-vcf runs/full-golden-fresh/runtime/prepared_input.vcf \
  --output-vcf runs/full-golden-fresh/runtime/vepyr.annotated.vcf \
  --cache-dir /Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/.cache/vepyr_cache \
  --reference-fasta /Users/lukaszjezapkowicz/.vep/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
```

Optional:

- add `--use-fjall` when the cache root also contains `*.fjall` stores
- add `--log-path /path/to/vepyr.log` to override the default log destination

### 5b. Benchmark compare-only on an existing annotated pair

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli benchmark-compare \
  --left-vcf runs/full-golden-fresh/runtime/vep.annotated.vcf \
  --right-vcf runs/full-golden-fresh/runtime/vepyr.annotated.vcf \
  --output-json /tmp/vepyr-diffly-benchmark.json
```

This command runs `compare-existing` in temporary directories, writes only the compact benchmark summary JSON, and cleans the heavy temp artifacts automatically.

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

`runtime/vep.log` and `runtime/vepyr.log` are written incrementally while annotation is running, so they can be tailed live during long runs.
Each log also records `STARTED`, `ENDED`, and `EXIT` markers to make elapsed-time checks easier.

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

For small files the compare can still run directly on normalized tables. For larger files the pipeline now uses a bounded-memory path for both tiers: variant rows and consequence rows are streamed in chunks, spilled into hash buckets on disk, compacted per bucket, and then compared bucket-by-bucket. In `fast` mode each bucket is prechecked first and `diffly` runs only on buckets that actually differ. This keeps working memory tied to the configured memory budget instead of full input size.

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
- run an exact precheck per bucket
- in `fast` mode, skip `diffly` for buckets already proven equal
- in `debug` mode, compare every bucket with `diffly`
- merge bucket diff artifacts at the end

This keeps memory bounded and gives progress reporting for large compares.

## Results

### A. Smoke Tests

The most recent post-fix smoke matrix that was re-verified locally is:

- `10`
- `100`
- `1000`
- `10000`

For each size the workflow was:

- `run --compare-mode fast`
- then `compare-existing --compare-mode debug`
- with the same annotated outputs
- under a bounded-memory budget

Observed outcome:

| Sample size | Prepared rows | `fast` variant equal | `fast` consequence equal | `debug` variant equal | `debug` consequence equal |
| --- | ---: | --- | --- | --- | --- |
| `10` | `10` | `true` | `true` | `true` | `true` |
| `100` | `100` | `true` | `true` | `true` | `true` |
| `1000` | `1000` | `true` | `true` | `true` | `true` |
| `10000` | `10066` | `true` | `true` | `true` | `true` |

Observations:

- all currently re-verified smoke points are fully equal at both tiers
- the `vepyr` import path issue is fixed; the runner now prefers the working `.vepyr` environment before falling back to the local checkout
- multi-allelic decomposition is still visible at `10000`, where `10000` source variants became `10066` prepared single-alt rows
- the compare path is green in both `fast` and `debug` on small and medium-sized fixtures after the RAM-safe refactor

### B. Golden Test

The current retained full compare-only result is [runs/full-golden-compare-fast](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast), run in `fast` mode against:

- [vep.annotated.vcf](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-fresh/runtime/vep.annotated.vcf)
- [vepyr.annotated.vcf](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-fresh/runtime/vepyr.annotated.vcf)

Effective plan from [summary.json](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast/summary.json):

- `memory_budget_mb=1024`
- `bucket_count=256`
- `variant_chunk_rows=10000`
- `consequence_chunk_rows=8192`
- `compare_workers=2`
- `parallelize_sides=false`

Observed compare-only timings:

| Stage | Seconds | Approx. |
| --- | ---: | ---: |
| Variant summary | `511.181` | `8m 31s` |
| Consequence bucketization | `1665.249` | `27m 45s` |
| Schema validation | `0.027` | `<1s` |
| Variant diff | `2.357` | `2s` |
| Consequence diff | `761.191` | `12m 41s` |
| Total compare-only | `2940.005` | `49m 00s` |

Observed result:

- variant tier: `equal=true`
- variant joined equal rows: `4096123`
- variant `left_only=0`, `right_only=0`, `unequal=0`
- consequence tier: `equal=false`
- consequence joined equal rows: `36914657`
- consequence `left_only=392`, `right_only=392`, `unequal=0`

What this means:

- at the variant level VEP and `vepyr` are fully equal
- at the consequence level there are a small number of rows that exist only on one side
- there are no cases where the same normalized consequence key exists on both sides with different payload values

Why this run still takes a long time:

- the dominant cost is not the final diff itself, but preprocessing the annotated VCFs into bounded-memory bucket artifacts
- `consequence_bucketization` alone is about `28` minutes
- even in `fast` mode, only `9/256` consequence buckets were proven equal by precheck, so `247/256` buckets still needed exact diff
- the inputs are very large annotated VCFs, so parsing `CSQ`, normalizing transcript consequences, spilling to disk, and compacting buckets dominate wall-clock time
- the run deliberately used a conservative `1024 MB` memory budget, which reduces RAM pressure but also reduces concurrency

Current strongest suspicions about the remaining differences:

- the mismatch pattern looks like transcript-selection or transcript-membership drift, not a gross variant-level bug
- the biggest clusters are in `inframe_insertion`, `downstream_gene_variant`, `5_prime_UTR_variant`, and transcript-specific insertion consequences
- many mismatches are paired `left_only/right_only` rows around the same locus but for different `ENST...` transcripts
- representative hotspots include `HADHB`, `LINC03025`, `ACIN1`, `CCDC66`, and `ABCF1`
- this suggests that VEP and `vepyr` often agree on the variant and broad effect family, but disagree on the exact transcript consequence set emitted for some loci

Most useful artifacts for debugging the golden mismatch:

- [summary.json](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast/summary.json)
- [summary.md](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast/summary.md)
- [consequence_mismatches.tsv](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast/consequence_mismatches.tsv)
- [consequence_diff.parquet](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast/consequence_diff.parquet)
- [runtime/compare.progress.log](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast/runtime/compare.progress.log)

## Current Scope

Current verified scope:

- one active preset: `ensembl_everything`
- semantic comparison of annotated VCF `INFO/CSQ`
- two comparison tiers: variant and consequence
- bounded-memory compare path for large annotated VCFs
- local runtime execution against a real local VEP installation
- local runtime execution against an existing `vepyr` Python environment and parquet cache
- compare-only reuse of existing annotated VCFs via `compare-existing`

## Python Environment

The root [requirements.txt](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/requirements.txt) installs:

- the local package in editable mode
- runtime dependencies from [pyproject.toml](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/pyproject.toml)
- `pytest`
- `ruff`

Recommended setup:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
cp .env.example .env
```

Quick verification:

```bash
source .venv/bin/activate
PYTHONPATH=src python -m vepyr_diffly.cli list-presets
PYTHONPATH=src pytest -q
```

## Repo Configuration

Runtime paths should live in repo-local `.env`.

Important variables:

- `VEPYR_DIFFLY_EXECUTION_MODE`
- `VEPYR_DIFFLY_INPUT_VCF`
- `VEPYR_DIFFLY_OUTPUT_DIR`
- `VEPYR_DIFFLY_VEP_CACHE_DIR`
- `VEPYR_DIFFLY_VEP_CACHE_VERSION`
- `VEPYR_DIFFLY_REFERENCE_FASTA`
- `VEPYR_DIFFLY_VEP_BIN`
- `VEPYR_DIFFLY_VEPYR_PYTHON`
- `VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR`
- `VEPYR_DIFFLY_COMPARE_MODE`
- `VEPYR_DIFFLY_MEMORY_BUDGET_MB`

Notes:

- `VEPYR_DIFFLY_VEPYR_PYTHON` should point to a working virtualenv interpreter such as `.../.vepyr/bin/python`
- `VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR` should point to the cache root, not directly to the deepest parquet leaf
- the CLI loads `.env` automatically

## CLI

Main commands:

```bash
python -m vepyr_diffly.cli list-presets
python -m vepyr_diffly.cli run --preset ensembl_everything --input-vcf /path/to/input.vcf --output-dir runs/demo --sample-first-n 1000
python -m vepyr_diffly.cli compare-existing --preset ensembl_everything --left-vcf /path/to/vep.annotated.vcf --right-vcf /path/to/vepyr.annotated.vcf --output-dir runs/compare-only --compare-mode fast --memory-budget-mb 1024
python -m vepyr_diffly.cli annotate-vepyr --input-vcf /path/to/prepared_input.vcf --output-vcf /path/to/vepyr.annotated.vcf
python -m vepyr_diffly.cli benchmark-compare --left-vcf /path/to/vep.annotated.vcf --right-vcf /path/to/vepyr.annotated.vcf --output-json /tmp/benchmark.json
python -m vepyr_diffly.cli inspect-run --run-dir runs/demo
```

## Expected Output

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

Main run artifacts:

- `summary.json`
- `summary.md`
- `variant_diff.parquet`
- `consequence_diff.parquet`
- `variant_mismatches.tsv`
- `consequence_mismatches.tsv`
- `runtime/effective_config.json`
- `runtime/compare.progress.log`
- `runtime/vep.log`
- `runtime/vepyr.log`

On large runs the TSV mismatch files are intentionally sampled, while the parquet diff artifacts remain complete.

## Troubleshooting

If tests fail:

```bash
source .venv/bin/activate
PYTHONPATH=src pytest -q -vv
```

If annotation fails before compare:

- inspect `runtime/vep.log`
- inspect `runtime/vepyr.log`
- inspect `runtime/effective_config.json`

Typical causes:

- wrong `VEPYR_DIFFLY_VEP_BIN`
- wrong `VEPYR_DIFFLY_VEP_CACHE_DIR`
- wrong `VEPYR_DIFFLY_REFERENCE_FASTA`
- wrong `VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR`
- `VEPYR_DIFFLY_VEPYR_PYTHON` does not have a working `vepyr` install

If `vepyr` fails to import:

- check that `VEPYR_DIFFLY_VEPYR_PYTHON` points to a real virtualenv interpreter like `.../.vepyr/bin/python`
- avoid resolving that path to the system interpreter

If a long compare is still running, tail:

```bash
tail -f runs/<run>/runtime/compare.progress.log
```

The most important milestones are:

- `parsed CSQ header`
- `materializing variant summary`
- `bucketizing consequence rows`
- `comparing ... buckets`
- `writing summaries`
- `completed successfully`

## Runtime Notes

- canonical equality is semantic comparison of parsed `CSQ`
- raw VCF text equality is out of scope
- the currently verified execution path is local runtime
- the compare path is now bounded-memory and driven by the configured memory budget
- for large runs, choosing a smaller memory budget trades throughput for lower RAM pressure

## Known Limitations

- only the `ensembl_everything` preset is currently verified end-to-end
- `merged` and `refseq` cache flavors are not implemented and verified yet
- Docker/container execution remains scaffolded in the repo, but local execution is the currently verified path
- the retained full golden result currently documents `compare-existing --compare-mode fast`; a fresh full end-to-end golden `run` with the same current codepath is still worth recording separately
- the repo compares semantic normalized outputs, not byte-for-byte annotated VCF files
- the current docs assume the local external assets already exist on this machine: VEP binary, VEP cache, FASTA, `vepyr` Python environment, and `vepyr` parquet cache
