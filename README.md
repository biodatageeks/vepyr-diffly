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

For local plugin cache builds, configure plugin source files in `.env` as needed:

- `VEPYR_DIFFLY_PLUGIN_CLINVAR_SOURCE`
- `VEPYR_DIFFLY_PLUGIN_SPLICEAI_SOURCE`
- `VEPYR_DIFFLY_PLUGIN_CADD_SNV_SOURCE`
- `VEPYR_DIFFLY_PLUGIN_CADD_INDEL_SOURCE`
- `VEPYR_DIFFLY_PLUGIN_ALPHAMISSENSE_SOURCE`
- `VEPYR_DIFFLY_PLUGIN_DBNSFP_SOURCE`

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
  --chromosomes 1,2,3 \
  --compare-mode fast \
  --memory-budget-mb 1024
```

Useful compare tuning flags:

- `--chromosomes 1` or `--chromosomes 1,2,3`: process only the selected chromosomes; standard aliases like `1` and `chr1` are treated as the same chromosome
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

### 5c. Build a local `vepyr` cache for selected chromosomes

Use the helper script in `scripts/` when you want to materialize the local cache outside a compare run.

Build only chromosome `Y` plus the default plugin set:

```bash
source .venv/bin/activate
python scripts/build_chr_cache.py --chromosomes Y
```

Build chromosomes `1` and `Y`:

```bash
source .venv/bin/activate
python scripts/build_chr_cache.py --chromosomes 1,Y
```

Build the full cache:

```bash
source .venv/bin/activate
python scripts/build_chr_cache.py
```

By default the script:

- installs the local checkout from `VEPYR_DIFFLY_VEPYR_PATH` with `pip install -e --no-build-isolation`
- builds the core cache under `VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR`
- builds core `variation.fjall` and `translation_sift.fjall` from the generated parquet cache
- derives the Ensembl source cache from `VEPYR_DIFFLY_VEP_CACHE_DIR`
- builds plugin caches for `clinvar`, `spliceai`, `cadd`, `alphamissense`, and `dbnsfp`
- requires each requested plugin to be provided as a local source file via `.env` or CLI
- does not download plugin sources from the Internet in the standard build path
- materializes `cadd` from two local source files into one shared cache root layout: `cadd/` and `cadd.fjall/`

Example local plugin source configuration:

```bash
export VEPYR_DIFFLY_PLUGIN_CLINVAR_SOURCE=/data/plugins/clinvar.vcf.gz
export VEPYR_DIFFLY_PLUGIN_CADD_SNV_SOURCE=/data/plugins/whole_genome_SNVs.tsv.gz
export VEPYR_DIFFLY_PLUGIN_CADD_INDEL_SOURCE=/data/plugins/gnomad.genomes.r4.0.indel.tsv.gz
export VEPYR_DIFFLY_PLUGIN_DBNSFP_SOURCE=/data/plugins/dbNSFP4.9c_grch38.gz
```

Useful flags:

- `--no-plugins` to skip plugin cache generation
- `--no-core-fjall` to skip `variation.fjall` / `translation_sift.fjall`
- `--only-plugins` to skip the core cache and build only plugin caches
- `--chromosomes 1,4,10` to restrict both core and plugin generation to the selected chromosomes
- `--plugins clinvar,spliceai` to limit the plugin set
- `--preview-rows 1000` to build both core and plugin caches from preview-sized source slices; for plugins this means first `N` data rows, while for the binary Ensembl core cache it means copying only the first preview-sized set of whole source shard files
- when `--chromosomes` is used for plugins, preview slicing is chromosome-aware and uses `tabix` automatically when a matching `.tbi` index exists
- `--clean-plugin-output` to delete only the current-layout outputs of the requested plugins before rebuilding (`<version>/<plugin>/` and `<version>/<plugin>.fjall/`), which is useful for repeatable size and timing comparisons
- `--remove-old-layout-cache` to delete recognized stale artifacts from the previous cache layout before building (`parquet/<version>`, root-level `*.fjall`, and root-level plugin directories)
- `--assume-sorted-plugin-input` to skip SQL `ORDER BY` for single-source plugin builds when the raw input is already sorted by `chrom,pos,ref,alt`; this is a safe opt-in and is intentionally ignored for `cadd`
- `--clinvar-source`, `--spliceai-source`, `--cadd-snv-source`, `--cadd-indel-source`, `--alphamissense-source`, `--dbnsfp-source` to override `.env` per run
- `--force-plugin-source` is retained for CLI compatibility but ignored for local-source builds
- `--skip-install` if `vepyr` is already installed in the active interpreter

Build only plugin caches without the core cache:

```bash
source .venv/bin/activate
python scripts/build_chr_cache.py --only-plugins --plugins clinvar,spliceai,cadd,alphamissense
```

Build only plugin preview caches from the first `1000` source rows:

```bash
source .venv/bin/activate
python scripts/build_chr_cache.py --only-plugins --plugins clinvar,spliceai,cadd,alphamissense --preview-rows 1000
```

Build repeatable chromosome-scoped preview caches and clean only the selected plugin outputs first:

```bash
source .venv/bin/activate
python scripts/build_chr_cache.py \
  --only-plugins \
  --plugins clinvar,spliceai,cadd,alphamissense,dbnsfp \
  --chromosomes 1 \
  --preview-rows 1000 \
  --clean-plugin-output \
  --skip-install
```

### 5d. Run a plugin round-trip validation check

Use the dedicated helper when you want an explicit `convert -> read back -> verify row counts and schema` check for plugin caches.

Example:

```bash
source .venv/bin/activate
python scripts/plugin_round_trip_test.py --plugins clinvar,spliceai,cadd,alphamissense --preview-rows 1000
```

What this script does for each plugin:

- creates a preview source if `--preview-rows` is set
- runs plugin conversion through the local `vepyr` build path
- reads back the generated parquet files
- verifies row counts: expected rows from the preview source vs actual parquet rows
- verifies schema: expected plugin columns and Arrow types vs the written parquet schema
- prints explicit per-plugin stages: `convert`, `read-back`, `verify rows`, `verify schema`, and final `round-trip`

Useful flags:

- `--plugins clinvar,spliceai` to limit the plugin set
- `--preview-rows 1000` to keep the check fast on very large sources
- `--skip-install` if the active environment already has the correct editable `vepyr`
- `--cache-dir /path/to/output` to keep outputs in a known location
- `--keep-cache` to keep the generated parquet files instead of deleting the temporary cache directory after the check

### 5e. Create `.tbi` indexes for plugin source files

Use the dedicated helper when you want to prepare tabix indexes for plugin
inputs in `plugins/`.

Dry-run all known plugin sources:

```bash
source .venv/bin/activate
python scripts/create_plugin_indexes.py --dry-run
```

Index only `dbnsfp`, `spliceai`, and `cadd`:

```bash
source .venv/bin/activate
python scripts/create_plugin_indexes.py --plugins dbnsfp,spliceai,cadd
```

For plain gzip inputs such as `clinvar.vcf.gz` and `AlphaMissense_hg38.tsv.gz`,
create sibling BGZF files first and index those:

```bash
source .venv/bin/activate
python scripts/create_plugin_indexes.py --plugins clinvar,alphamissense --recompress-plain-gzip
```

Notes:

- BGZF inputs are indexed in place.
- Plain gzip inputs cannot be indexed directly by `tabix`; they must be converted to BGZF first.
- The recompress mode writes sibling `.bgz` files and leaves the original `.gz` untouched.
- The helper prints `[i/n]` progress and periodic elapsed-time updates while `bgzip` or `tabix` is still running.
- These `.tbi` files can help future chromosome-scoped plugin-cache builds; they do not speed up annotation from already-built parquet/fjall caches.

### 5f. Run a plugin annotation smoke test from built cache

Use this helper when you want to verify that `vepyr.annotate()` really picks up
plugin values from the already-built plugin cache.

It:

- reads one variant from each selected plugin cache parquet directory
- writes a temporary mini VCF from those cache-backed variants
- runs `vepyr.annotate()` against the partitioned cache
- checks that at least one plugin-specific field is populated for each plugin

Example:

```bash
source .venv/bin/activate
python scripts/plugin_annotation_smoke_test.py --plugins clinvar,spliceai,cadd,alphamissense,dbnsfp --skip-install
```

Notes:

- The script uses `use_fjall=False` deliberately, so it exercises the parquet-backed annotation path.
- It expects plugin caches under `.cache/vepyr_cache/115_GRCh38_vep/<plugin>/`.
- It expects `VEPYR_DIFFLY_REFERENCE_FASTA` in `.env`, or you can pass `--reference-fasta`.
- If the existing local core cache contains corrupted parquet files from an interrupted older build, the smoke test can fail before plugin lookup. Rebuild the core cache into a clean cache root in that case.

### 5g. Issue Validation Artifact

For a step-by-step local proof that the current plugin-cache work satisfies the original plugin issue using reduced-scope evidence, see:

- [`github_issue_check.MD`](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/github_issue_check.MD)

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

For small files the compare can still run directly on normalized tables. For larger files the pipeline now uses a bounded-memory path for both tiers: variant rows and consequence rows are streamed in chunks, spilled into hash buckets on disk, compacted per bucket, and then compared bucket-by-bucket. A resource planner derives chunk sizes, bucket count, worker count, and side parallelism from `--memory-budget-mb`, so working memory scales with the configured budget instead of full input size.

The pipeline can now also be scoped to selected chromosomes. The filter is available on both `run` and `compare-existing`, is applied before annotation in `run`, and during streaming normalization in `compare-existing`. Standard naming aliases are normalized, so `1` matches both `1` and `chr1`, and `MT` matches `MT`, `M`, `chrMT`, and `chrM`.

The pipeline is:

1. prepare the input VCF
2. annotate the same prepared input with VEP and `vepyr`
3. normalize both annotated VCFs into comparable tables
4. compare those tables with `diffly`
5. write console summary plus file artifacts

### Input preparation

Before annotation, `run` prepares a canonical input:

- if sampling is enabled, it writes `runtime/sampled_input.vcf`
- if chromosome filtering is enabled, it writes `runtime/filtered_input.vcf`
- it then decomposes multi-allelic rows into single-alt rows
- it writes the final annotation input to `runtime/prepared_input.vcf`
- it records overall and per-chromosome preparation stats in `runtime/input_preparation.json`

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

For large runs, variant comparison now follows the same model:

- stream variant rows into hash buckets on disk
- compact each bucket into deterministic `bucket.parquet`
- compare the variant tier bucket-by-bucket instead of one global in-memory diff
- keep the old eager path only for safely small files

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

The newest retained full compare-only result is [runs/full-golden-compare-fast-per-chrom-20260330](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast-per-chrom-20260330), run in `fast` mode against:

- [vep.annotated.vcf](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-fresh/runtime/vep.annotated.vcf)
- [vepyr.annotated.vcf](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-fresh/runtime/vepyr.annotated.vcf)

Effective plan from [summary.json](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast-per-chrom-20260330/summary.json):

- `memory_budget_mb=1024`
- `bucket_count=256`
- `variant_chunk_rows=10000`
- `consequence_chunk_rows=8192`
- `compare_workers=2`
- `parallelize_sides=false`

Observed compare-only timings:

| Stage | Seconds | Approx. |
| --- | ---: | ---: |
| Variant summary | `518.632` | `8m 39s` |
| Consequence bucketization | `999.975` | `16m 40s` |
| Schema validation | `0.031` | `<1s` |
| Variant diff | `4.081` | `4s` |
| Consequence diff | `490.070` | `8m 10s` |
| Total compare-only | `2012.789` | `33m 33s` |

Observed result:

- variant tier: `equal=true`
- variant joined equal rows: `4096123`
- variant `left_only=0`, `right_only=0`, `unequal=0`
- consequence tier: `equal=false`
- consequence joined equal rows: `36914657`
- consequence `left_only=392`, `right_only=392`, `unequal=0`

What this means:

- at the variant level VEP and `vepyr` are fully equal
- at the consequence level there are still a small number of rows that exist only on one side
- those rows are not random noise: the sampled mismatch artifact collapses cleanly into `392` transcript-level `left_only/right_only` pairs
- there are still no cases where the same normalized consequence key exists on both sides with different payload values inside one joined key

Why this run still takes a long time:

- the dominant cost is not the final diff itself, but preprocessing the annotated VCFs into bounded-memory bucket artifacts
- `consequence_bucketization` alone is still about `16m 40s`
- even in `fast` mode, only a tiny minority of consequence buckets are proven equal by precheck, so the final diff still has to run exact compare for most buckets
- the inputs are very large annotated VCFs, so parsing `CSQ`, normalizing transcript consequences, spilling to disk, and compacting buckets dominate wall-clock time
- the run deliberately used a conservative `1024 MB` memory budget, which reduces RAM pressure but also reduces concurrency

Per-chromosome view from the same summary:

- all chromosomes are `variant.equal=true`
- only chromosome `1` is fully equal on the consequence tier
- the largest consequence mismatch counts are on chromosomes `2` (`88/88`), `9` (`54/54`), `6` (`44/44`), `3` (`39/39`), and `14` (`37/37`)
- the per-chromosome timings below are attributed timings from the summary, not independent wall-clock slices, so they are best used for relative weight and hotspot analysis

| Chr | Variant | Consequence | Left/Right only | Attributed total | Variant summary | Consequence bucketization | Consequence diff |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | `yes` | `yes` | `0 / 0` | `161.248s` | `40.951s` | `80.516s` | `39.459s` |
| 2 | `yes` | `no` | `88 / 88` | `166.191s` | `41.951s` | `83.156s` | `40.754s` |
| 3 | `yes` | `no` | `39 / 39` | `158.916s` | `36.532s` | `81.940s` | `40.157s` |
| 4 | `yes` | `no` | `10 / 10` | `135.640s` | `38.908s` | `64.712s` | `31.714s` |
| 5 | `yes` | `no` | `8 / 8` | `110.990s` | `33.478s` | `51.842s` | `25.407s` |
| 6 | `yes` | `no` | `44 / 44` | `129.702s` | `34.435s` | `63.752s` | `31.244s` |
| 7 | `yes` | `no` | `9 / 9` | `114.174s` | `29.694s` | `56.538s` | `27.708s` |
| 8 | `yes` | `no` | `26 / 26` | `107.026s` | `28.519s` | `52.536s` | `25.747s` |
| 9 | `yes` | `no` | `54 / 54` | `83.531s` | `22.298s` | `40.976s` | `20.082s` |
| 10 | `yes` | `no` | `8 / 8` | `105.041s` | `27.028s` | `52.212s` | `25.588s` |
| 11 | `yes` | `no` | `12 / 12` | `100.479s` | `26.187s` | `49.719s` | `24.367s` |
| 12 | `yes` | `no` | `8 / 8` | `95.336s` | `25.046s` | `47.040s` | `23.053s` |
| 13 | `yes` | `no` | `7 / 7` | `61.341s` | `20.438s` | `27.342s` | `13.400s` |
| 14 | `yes` | `no` | `37 / 37` | `68.388s` | `16.865s` | `34.488s` | `16.902s` |
| 15 | `yes` | `no` | `6 / 6` | `68.004s` | `15.850s` | `34.917s` | `17.112s` |
| 16 | `yes` | `no` | `13 / 13` | `62.537s` | `15.619s` | `31.404s` | `15.391s` |
| 17 | `yes` | `no` | `3 / 3` | `65.283s` | `13.722s` | `34.530s` | `16.923s` |
| 18 | `yes` | `no` | `1 / 1` | `56.589s` | `15.116s` | `27.753s` | `13.601s` |
| 19 | `yes` | `no` | `5 / 5` | `55.707s` | `11.484s` | `29.618s` | `14.515s` |
| 20 | `yes` | `no` | `11 / 11` | `41.981s` | `11.003s` | `20.731s` | `10.160s` |
| 21 | `yes` | `no` | `2 / 2` | `29.526s` | `7.067s` | `15.035s` | `7.368s` |
| 22 | `yes` | `no` | `1 / 1` | `35.130s` | `6.440s` | `19.220s` | `9.419s` |

What the `392 / 392 / 0` mismatches actually are:

- they are not missing transcripts versus extra transcripts in the old sense
- each sampled mismatch row on the VEP side has a matching transcript-level partner on the `vepyr` side with the same `chrom`, `pos`, `ref`, `alt`, and `Feature`
- the mismatch is therefore usually a payload drift for the same transcript consequence, not a different transcript set
- the top low-level field signatures are:
  - only `HGVSp`: `121` pairs
  - only `Consequence`: `72`
  - only `HGNC_ID`: `63`
  - only `HGVSc`: `54`
  - `Consequence + IMPACT`: `40`
- the top semantic classes are:
  - `hgvs_payload_drift`: `114`
  - `missing_hgvs`: `84`
  - `consequence_reclassification`: `79`
  - `missing_hgnc_id`: `63`
  - `consequence_and_impact_reclassification`: `55`

In practice that means the `392` mismatch pairs are mostly:

- rows where both sides point at the same transcript consequence, but `HGVSp` uses different protein coordinates or wording
- rows where VEP has `HGVSc` or `HGVSp` filled and `vepyr` leaves it empty
- rows where VEP and `vepyr` classify the same transcript consequence differently, for example `inframe_insertion` versus `inframe_insertion&start_retained_variant`
- rows where the same transcript consequence gets a different `IMPACT`, for example `HIGH` versus `LOW`
- rows where `HGNC_ID` is missing on one side, most often present in `vepyr` and empty in VEP

Representative raw `CSQ` examples confirm that these are already present in the source annotated VCFs:

- `missing_hgvs`:
  - VEP contains `ENST00000944442.1:c.-111_-110ins...`
  - `vepyr` emits the same row with an empty `HGVSc`
- `consequence_and_impact_reclassification`:
  - VEP: `frameshift_variant&splice_region_variant|HIGH`
  - `vepyr`: `splice_region_variant|LOW`
- `consequence_reclassification`:
  - VEP: `inframe_insertion`
  - `vepyr`: `inframe_insertion&start_retained_variant`
- `missing_hgnc_id`:
  - VEP: empty `HGNC_ID`
  - `vepyr`: populated `HGNC:...`
- `hgvs_payload_drift`:
  - both sides agree on `Consequence`, `IMPACT`, and `HGVSc`
  - but `HGVSp` refers to different duplicated amino-acid coordinates

Follow-up optimization note:

- after the older retained full run, the consequence compare path was tightened further so bucketed compare no longer writes unused per-bucket TSV files, can use cheap bucket metadata to skip the expensive eager precheck for buckets that are already provably different, and now uses in-process worker threads instead of spawning separate compare subprocesses
- the fresh full rerun above confirms those changes on the full pipeline: `consequence_diff` is now `490.070s` instead of the older retained `761.191s`, which is about `35.6%` faster
- on the same machine and with the same full-golden `vepyr.annotated.vcf`, a one-side `materialize_consequence_buckets(...)` benchmark at `bucket_count=256` and `chunk_variants=8192` first dropped from `708.518s` to `654.883s` after switching temporary bucket part writes to lighter parquet settings (`lz4`, no statistics), and then to `447.560s` after buffering several reduced chunks before each flush
- that latest bucketization benchmark is about `36.8%` faster than the original one-side baseline, and it also cuts temporary file churn dramatically because bucket parts are flushed around every `65536` input variants instead of every `8192`
- the fresh full rerun also reflects that in end-to-end compare wall-clock:
  - older retained compare-only baseline: `49m 00s`
  - fresh rerun on current code: `33m 33s`

Latest retained mismatch diagnosis:

- all `784` sampled consequence mismatch rows collapse cleanly into `392` transcript-level `left_only/right_only` pairs
- there are `0` unpaired rows in that sampled mismatch artifact
- the top differing fields are still `HGVSp`, `Consequence`, `HGNC_ID`, `IMPACT`, and `HGVSc`
- the top semantic drift classes are:
  - `hgvs_payload_drift`: `114`
  - `missing_hgvs`: `84`
  - `consequence_reclassification`: `79`
  - `missing_hgnc_id`: `63`
  - `consequence_and_impact_reclassification`: `55`
- raw `CSQ` extraction on representative cases confirms that these differences already exist in the source annotated VCF payloads, not just in downstream normalization
- the current debug tooling can now go from mismatch TSV -> semantic category summary -> raw left/right `CSQ` examples without rescanning the whole mismatch space manually

Most useful artifacts for debugging the golden mismatch:

- [summary.json](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast-per-chrom-20260330/summary.json)
- [summary.md](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast-per-chrom-20260330/summary.md)
- [consequence_mismatches.tsv](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast-per-chrom-20260330/consequence_mismatches.tsv)
- [consequence_diff.parquet](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast-per-chrom-20260330/consequence_diff.parquet)
- [runtime/compare.progress.log](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/runs/full-golden-compare-fast-per-chrom-20260330/runtime/compare.progress.log)
- `/tmp/vepyr-diffly-golden-mismatch-analysis.json`
- `/tmp/vepyr-diffly-golden-raw-csq.json`

Logs and structured diagnostics added for this:

- `runtime/compare.progress.log`: live stage and per-bucket compare progress for the full run
- `consequence_mismatch_analysis.json`: grouped transcript-level mismatch summary, field counts, semantic categories, and representative examples
- `raw-csq-examples.json`: raw left/right `CSQ` payloads for representative cases in each semantic category

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
python -m vepyr_diffly.cli run --preset ensembl_everything --input-vcf /path/to/input.vcf --output-dir runs/demo --sample-first-n 1000 --chromosomes 1,2
python -m vepyr_diffly.cli compare-existing --preset ensembl_everything --left-vcf /path/to/vep.annotated.vcf --right-vcf /path/to/vepyr.annotated.vcf --output-dir runs/compare-only --chromosomes 1,2 --compare-mode fast --memory-budget-mb 1024
python -m vepyr_diffly.cli annotate-vepyr --input-vcf /path/to/prepared_input.vcf --output-vcf /path/to/vepyr.annotated.vcf
python -m vepyr_diffly.cli benchmark-compare --left-vcf /path/to/vep.annotated.vcf --right-vcf /path/to/vepyr.annotated.vcf --output-json /tmp/benchmark.json
python -m vepyr_diffly.cli inspect-run --run-dir runs/demo
python -m vepyr_diffly.cli analyze-consequence-mismatches --run-dir runs/full-golden-compare-fast --output-json /tmp/mismatch-analysis.json
python -m vepyr_diffly.cli extract-mismatch-csq-examples --run-dir runs/full-golden-fresh --analysis-json /tmp/mismatch-analysis.json --output-json /tmp/raw-csq-examples.json
```

`analyze-consequence-mismatches` reads `consequence_mismatches.tsv`, groups `left_only/right_only` rows back into transcript-level pairs, counts the fields that drift most often, and emits representative examples. It is the fastest way to turn a retained run into a concrete CSQ drift report.

`extract-mismatch-csq-examples` takes that analysis and enriches a few representative cases per semantic category with the original raw `CSQ=` entries from both annotated VCFs. This is the shortest path from "which categories drift?" to "what exactly did VEP and vepyr emit for the same transcript consequence?".

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

`summary.json` and `summary.md` now include both:

- overall results for the processed subset as a whole
- per-chromosome variant/consequence equality and counts
- attributed per-chromosome stage timings for normalization and diff

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
