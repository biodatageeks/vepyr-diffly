# Reproducibility package — master's thesis on VEPYR

Thesis: *Równoległa oraz rozproszona implementacja rozwiązania do adnotacji
wykrytych zmian genomicznych* (Parallel and Distributed Implementation of a
System for Annotating Detected Genomic Variants), Łukasz Jezapkowicz,
Warsaw University of Technology, 2026.

This directory holds the small raw artifacts behind the numbers reported in the
thesis (chapter 5 and appendix A): per-run timing and memory records, summaries,
command logs, concordance reports and the measurement scripts. Large outputs
(annotated VCFs, Parquet caches, source data) are not included; they are listed
in `LARGE_ARTIFACTS.tsv`. `MANIFEST.sha256` lists a SHA-256 checksum for every
file in this directory.

**Status (2026-10-23):** the October 2026 completion run (`results/completion-2026-10`,
steps S0–S11) and VEP full-genome repeats 3–5 are still running; their results
will be added to this directory when they finish.

## Environment

- Hardware: MacBook Air, Apple M5 (4 performance + 6 efficiency cores), 32 GiB RAM.
  One additional test (independent processes) was originally run on an Apple M1
  machine and repeated on the M5 in October 2026 (see `results/completion-2026-10`).
- OS: macOS 27.0 (annotation series, September 2026), 27.0.1 (cache-build series, October 2026).
- Ensembl VEP 116.0 in Docker (`ensemblorg/ensembl-vep:release_116.0`) under Colima
  (Virtualization.framework, virtiofs, 10 vCPU); VEP cache 116, merged profile, GRCh38.
- Input: GIAB HG002 v4.2.1, chromosomes 1–22, normalised with `bcftools norm -m -both`
  (4,048,342 VCF records → 4,096,123 alleles).
- VEPYR used for annotation and concordance: vepyr `479d8f46eee2` (0.7.0 plus
  post-release fixes), datafusion-bio-functions `06fb397d6b47` (v0.22.1),
  datafusion-bio-formats `00487b9edf53` (v1.13.0); recorded in the `build`
  section of every E2E report.
- Cache-build experiments: code on the branches `thesis/cache-build-parallel`
  (biodatageeks/vepyr, biodatageeks/datafusion-bio-functions) and
  `thesis/cache-build-distributed` (biodatageeks/vepyr). See the READMEs in
  `scripts/thesis/` on those branches for provenance caveats (the measured
  native binary was built from an uncommitted state on 2026-09-30; the
  container image version is inferred).

## Map: thesis → artifacts

Times are means of five runs unless stated otherwise; "repeat 1" of the
annotation series is the single run of 2026-09-24, repeats 2–5 were run on
2026-09-25 (VEPYR) and 2026-09-26/28 (VEP).

| Thesis item (LaTeX label) | Artifacts | Script |
|---|---|---|
| VEPYR worker scaling chr1/chr22 (`tab:vepyr-scaling`) | `results/annotation-2026-09/5repeats_m5/` (`summary_5repeats.json`, `raw_runs.jsonl`); repeat 1: `chapter5-m5-20260924T061203Z/benchmarks/vepyr_chr*/` | `scripts/run_5repeats_benchmark.py`, `scripts/measure_once.py` |
| VEP `--fork` scaling chr1/chr22 (section `sec:speedup-investigation`) | `results/annotation-2026-09/vep_5repeats_m5/`; repeat 1: `chapter5-m5-.../benchmarks/vep/`; second chr1 no-fork value: `vep_chr1_warmup/` | `scripts/run_vep_5repeats_suite.py` |
| Plugin overhead O_x and VEP vs VEPYR per plugin (`tab:plugin-performance`, `tab:plugin-vep-vepyr`, `fig:plugin-overhead`) | `5repeats_m5/summary_5repeats.json` (`chr*_plugins`), `vep_5repeats_m5/summary_vep_5repeats.json` (`chr*_plugins`) | `scripts/fig_plugin_overhead.py` (figure) |
| All-plugins scaling (`tab:scaling-synthesis`) | `results/annotation-2026-10/plugin_scaling_m5_5repeats_20261013/runs.jsonl` | `repeat_plugin_scaling_m5.py` (thesis repository) |
| Full-input annotation (`tab:full-genome-performance`) | VEP: `results/vep-full-genome-2026-10/` (repeats 2–5, `summary.json`) and `chapter5-m5-.../benchmarks/vep/full_forknone.*` (repeat 1); VEPYR: `results/completion-2026-10/m5_completion_20261024/S5/` | `scripts/run_vep_full_repeats.sh`, `scripts/m5_completion_batch.py` |
| E2E concordance chr1/chr22, baseline and all plugins (`tab:variant-results`, `tab:unit-accounting`) | `results/concordance-2026-09/e2e-reports/fast_chr{1,22}_merged{,_plugins}_116_*` (JSON report, empty mismatch ledger, summaries); log: `chapter5-m5-.../correctness/correctness.log` | `run_comparison.py` in biodatageeks/vepyr `e2e-testing/scripts` |
| Single-plugin field comparison (`tab:plugin-field-results`) | `chapter5-m5-.../correctness/chr{1,22}_individual_plugins/benchmark_plugins_summary.{json,md}` | `vepyr-diffly benchmark-plugins` (this repository) |
| Cache build: dynamic queue, native sequential, containers (`tab:cache-lpt`) | `results/cache-build-2026-10/parallel_full_cache/queue_20_10_repeats/`, `full_core_5x/`, `full_cache_20/` (`distributed_5_summary.json`, `baseline_*_summary.json`) | branches `thesis/cache-build-*` |
| Full build including translation | `results/completion-2026-10/m5_completion_20261024/S11/` (repeat 1 from 2026-10-19 and repeats 2–5) | `full_e2e_all.py` on `thesis/cache-build-parallel` |
| Determinism `workers=1` vs `workers=8`, container start-up overhead, independent-process test, per-plugin CSQ field counts, cache sizes, full-genome E2E concordance, VEP peak RSS on chr1/chr22 | `results/completion-2026-10/m5_completion_20261024/S4, S3, S6, S7, S8, S10, S9` (`result.json`, raw JSONL, logs) | `scripts/m5_completion_batch.py` |

## Notes and deviations (also stated in the thesis)

- No warm-up runs were excluded; configurations ran in a fixed order.
- VEP time is the wall time of `docker run` measured with `/usr/bin/time -p`
  on the host (includes container start); VEPYR time is the wall time of the
  annotation call measured with `time.perf_counter` inside an already running
  Python process.
- VEP peak memory is the `VmHWM` of the VEP (perl) process inside the container,
  sampled with `docker exec`.
- Annotated VCFs from timing runs were deleted after each run; record counts and
  body checksums were kept.
- Paths in scripts and JSON files refer to the author's machine.
