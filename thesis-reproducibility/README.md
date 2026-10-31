# Reproducibility package — master's thesis on VEPYR

Thesis: *Równoległa oraz rozproszona implementacja rozwiązania do adnotacji
wykrytych zmian genomicznych* (Parallel and Distributed Implementation of a
System for Annotating Detected Genomic Variants), Łukasz Jezapkowicz,
Warsaw University of Technology, 2026.

This directory holds the small raw artifacts behind the numbers reported in the
thesis (chapter 5 and appendix A): per-run timing and memory records, summaries,
command logs, concordance reports and the measurement scripts. Large outputs
(annotated VCFs, Parquet caches, source data) are not included. `LARGE_ARTIFACTS.tsv`
lists only the large files that were skipped while copying the result directories;
caches and source data were never copied. `MANIFEST.sha256` lists a SHA-256 checksum
for every file in this directory.

**Status (2026-10-25):** complete. The October 2026 completion run
(`results/completion-2026-10`, steps S0–S11), VEP full-genome repeats 2–6, the
five-plugin references built on the M5 and the re-measured VEPYR plugin series have
finished; all numbers in the thesis refer to artifacts in this directory.

## Environment

- Hardware: MacBook Air, Apple M5 (4 performance + 6 efficiency cores), 32 GiB RAM.
  One additional test (independent processes) was originally run on an Apple M1
  machine and repeated on the M5 in October 2026 (see `results/completion-2026-10`).
- OS: macOS 27.0 (September 2026 runs), 27.0.1 (all October 2026 runs: cache builds, VEP
  full-genome repeats 2–6, completion run, references and plugin series of 2026-10-25).
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

Times are means of five runs. In the chromosome 1/22 annotation series (VEP and
VEPYR worker scaling, VEP plugins), "repeat 1" is the single run of 2026-09-24;
repeats 2–5 were run on 2026-09-25 (VEPYR) and 2026-09-26/28 (VEP). The VEPYR plugin
series on chromosomes 1/22 used in the thesis was measured again on 2026-10-25, all
five repeats with `measure_once.py` (`results/plugin-series-2026-10-25/`); the
September VEPYR plugin series stays in `results/annotation-2026-09/5repeats_m5/`.
VEP full-genome time and peak RSS are means over repeats 2–6 (2026-10-23/25,
`summary_repeats_2_6.json`); repeat 1 (2026-09-25 00:59–08:53 UTC, 28,469.05 s, no
RSS sampling) is kept but not averaged. VEPYR on the full input was run five times
on 2026-10-24 (step S5).

| Thesis item (LaTeX label) | Artifacts | Script |
|---|---|---|
| VEPYR worker scaling chr1/chr22 (`tab:vepyr-scaling`) | `results/annotation-2026-09/5repeats_m5/` (`summary_5repeats.json`, `raw_runs.jsonl`); repeat 1: `chapter5-m5-20260924T061203Z/benchmarks/vepyr_chr*/` | `scripts/run_5repeats_benchmark.py`, `scripts/measure_once.py` |
| VEP `--fork` scaling chr1/chr22 (section `sec:speedup-investigation`) | `results/annotation-2026-09/vep_5repeats_m5/`; repeat 1: `chapter5-m5-.../benchmarks/vep/`; second chr1 no-fork value: `vep_chr1_warmup/` | `scripts/run_vep_5repeats_suite.py` |
| Plugin overhead O_x and VEP vs VEPYR per plugin (`tab:plugin-performance`, `tab:plugin-vep-vepyr`, `fig:plugin-overhead`) | VEPYR: `results/plugin-series-2026-10-25/` (`raw_runs.jsonl`, `summary.json`, `environment.txt`, `run.log`); VEP: `results/annotation-2026-09/vep_5repeats_m5/summary_vep_5repeats.json` (`chr*_plugins`) | `scripts/rerun_plugin_series_5x.py` (calls `scripts/measure_once.py`), `scripts/after_repeat6_plugin_series.sh`, `scripts/fig_plugin_overhead.py` (figure, reads `raw_runs.jsonl`, writes the PDF used in the thesis) |
| Interleaved VEP/VEPYR series (D5, D6; `tab:interleaved`, `sec:interleaved`, interleaved rows of `tab:ratio-ci` and `tab:repro-dispersion`): baseline chr1 and chr22, 2026-10-31, five repeats per chromosome, tool order alternating per repeat; wall time, CPU time (VEP: container cgroup `cpu.stat`; VEPYR: `getrusage` user + sys) and output body SHA-256 of every run. All 20 bodies equal the VEP reference (chr1 `e9703e5e...`, chr22 `19efe3c7...`) | `results/interleaved-2026-10-31/` (`chr{1,22}/runs.jsonl`, `summary.json`, `chain.log`, `power.log`, `vep.time.txt` of the last VEP run) | `scripts/interleaved_vep_vepyr.py`, `scripts/interleaved_analysis.py` (summary, ratio CI as in `scripts/ratio_ci.py`) |
| Plugin cache of known provenance (C1, `sec:plugin-cache-rebuild`): the five plugin caches rebuilt for chr1 and chr22 on 2026-10-31 with vepyr 0.7.0 (479d8f4) from the vepyr-plugins v0.2.0 (7736bbe) source manifests, source MD5 verified (AlphaMissense: re-compressed BGZF file, digest differs from the upstream one by design, `verify_source='warn'`); sorted-row digests of all 10 plugin-chromosome files equal to the 2026-08-27 copy used in the thesis; annotation with all five plugins gives the same body on the new and on the old cache (chr22 `d0bafc91...`, chr1 `ce615dc8...`); this body differs from the VEP reference only in the order of the CSQ plugin fields (chr22: 50,861 records equal after mapping the 124 CSQ fields by name) | `results/plugin-cache-rebuild-2026-10-31/` (`c1_summary.json`, `c1-old_summary.json`, `c1_manifest_*.json`, `chr22_csq_field_order_check.json`, build logs) | `scripts/plugin_cache_rebuild.py` (modes `c1`, `c1-old`), `scripts/compare_cache_content.py`, `scripts/csq_field_order_check.py` |
| All-plugins scaling (`tab:scaling-synthesis`) | `results/annotation-2026-10/plugin_scaling_m5_5repeats_20261013/runs.jsonl` (use rows with `returncode` 0; 25 failed chr22 runs were repeated) | `scripts/repeat_plugin_scaling_m5.py`, which calls the snapshot of the VEPYR E2E scripts in `scripts/full_correctness_m5_20261012/` (`run_comparison.py`; annotation time from `comparison/annotate.py`, `time.time`, 0.1 s). On the author's machine the snapshot sat in `full_correctness_m5_20261012/scripts/` next to the runner script |
| Full-input annotation (`tab:full-genome-performance`) | VEP: `results/vep-full-genome-2026-10/` (repeats 2–6 in `summary.json`, means in `summary_repeats_2_6.json`; repeat 1, not averaged: `chapter5-m5-.../benchmarks/vep/full_forknone.*`); VEPYR: `results/completion-2026-10/m5_completion_20261024/S5/` | `scripts/run_vep_full_repeats.sh`, `scripts/m5_references_repeat6.sh` (repeat 6 and the 2–6 summary), `scripts/m5_completion_batch.py` |
| E2E concordance chr1/chr22, baseline and all plugins (`tab:variant-results`, `tab:unit-accounting`) | baseline: `results/completion-2026-10/m5_completion_20261024/S10/reports/fast_chr{1,22}_merged_116_*` (2026-10-24, against the full-genome VEP output of repeat 2); all plugins: `results/references-2026-10-25/e2e_reports/` (2026-10-25, against five-plugin references built on the M5 that day). The September comparisons against downloaded references (`results/concordance-2026-09/`, VEP headers 2026-07-27 and 2026-08-26) gave the same counts | `run_comparison.py` in biodatageeks/vepyr `e2e-testing/scripts` |
| Five-plugin VEP references chr1/chr22 (appendix A) | `results/references-2026-10-25/` (`chr*_reference_*`: body SHA-256, record count, VEP header, plugin code and source identities, build log; `environment.txt`, `run.log`) | `generate_vep_plugin_references.sh` in biodatageeks/vepyr `e2e-testing/scripts`, driven by `scripts/m5_references_repeat6.sh` |
| Single-plugin field comparison (`tab:plugin-field-results`) | `results/references-2026-10-25/chr{1,22}_individual_plugins/benchmark_plugins_summary.{json,md}` (2026-10-25; the September run in `chapter5-m5-.../correctness/` gave the same rows) | `vepyr-diffly benchmark-plugins` (this repository) |
| Cache build: dynamic queue, native sequential, containers (`tab:cache-lpt`) | `results/cache-build-2026-10/parallel_full_cache/queue_20_10_repeats/`, `full_core_5x/`, `full_cache_20/` (`distributed_5_summary.json`, `baseline_*_summary.json`) | branches `thesis/cache-build-*` |
| Full build including translation | `results/completion-2026-10/m5_completion_20261024/S11/` (repeat 1 from 2026-10-19 and repeats 2–5) | `full_e2e_all.py` on `thesis/cache-build-parallel` |
| Queue analysis: makespan vs lower bound, cost-estimate correlation per entity, Gantt chart (`fig:queue-gantt`, ch. 5) | `results/queue-analysis-2026-10-29/queue_analysis.json`, from the queue summary above | `scripts/queue_analysis.py [OUT_PDF]` |
| Annotation scaling and time-RSS figure (`fig:annotation-scaling`, ch. 5) | dispersion JSON above and `completion-2026-10/.../S5/result.json` | `scripts/fig_annotation_scaling.py OUT_PDF` |
| Negative control of the comparators (`tab:negative-control`): 9 injected difference classes + 3 equivalence classes, 20 variants each, E2E comparator and vepyr-diffly before (`4cc02d4`) and after (`1d657ce`, branch `fix/csq-allele-common-prefix`) the CSQ-allele fix | `results/negative-control-2026-10-29/` (`diffly_*/summary.{json,tsv}`, per-class `result.json`, E2E ledgers, diffly summaries) | `scripts/negative_control.py`, `scripts/negative_control_summary.py` |
| Annotation scaling diagnosis (`tab:scaling-diagnosis`): wall/CPU time, busy cores, RSS, block reads, VEP_PROFILE stage sums and output checksums for workers 1-10 (chr1 3x, chr1-22 1x; 2026-10-29). chr1-22 body SHA-256 equals the VEP repeat-2 body (5fee6224...) | `results/profiling-2026-10-29/` (`runs.jsonl`, `stderr_chr1_w*_profile.txt`) | `scripts/profile_annotation.py` |
| Scheduling-policy ablation by simulation (E3, `tab:queue-policies`): measured task times of the five p=10 queue runs replayed under 9 dispatch orders / static splits; simulated submitted order vs measured build time 0.996-0.9995 | `results/queue-policies-2026-10-30/policies_p10.json` | `scripts/simulate_queue_policies.py` |
| Edge-case representations (B8, `sec:limitations`): 17 chr22 records (multi-allelic SNV, `*`, `<NON_REF>`, MNV, complex delins, REF mismatch, ALT `.`, duplicate, left/right-aligned indels in a homopolymer, 59-bp deletion, symbolic DEL/DUP/INS); VEPYR workers=1 and workers=2 (no VEP run yet) | `results/edge-variants-2026-10-29/` (`edge_variants.vcf`, `cases.json`, `summary.json`, outputs, stderr) | `scripts/edge_variants.py` |
| Unusual contigs (B19, `sec:limitations`): VEPYR 0.7.0 on chr22, chrX, chrM, chrUn_KI270742v1, chrFAKE and on the same mitochondrial SNV spelled MT; records in/out, CSQ entries, consequences (no VEP run) | `results/edge-contigs-2026-10-29/` (`summary.json`, inputs, annotated VCFs, stderr) | `scripts/edge_contigs.py` |
| Dispersion of all timing configurations: median, SD (n-1), min-max, CV (appendix A, `tab:repro-dispersion`) | `results/dispersion-2026-10-29/` (`dispersion.json`, `dispersion.tsv`), computed from the raw files listed above | `scripts/dispersion_stats.py` |
| Determinism `workers=1` vs `workers=8`, container start-up overhead, independent-process test, per-plugin CSQ field counts, cache sizes, full-genome E2E concordance, VEP peak RSS on chr1/chr22 | `results/completion-2026-10/m5_completion_20261024/S4, S3, S6, S7, S8, S10, S9` (`result.json`, raw JSONL, logs) | `scripts/m5_completion_batch.py` |

## Notes and deviations (also stated in the thesis)

- No warm-up runs were excluded; configurations ran in a fixed order.
- VEP time is the wall time of `docker run` measured with `/usr/bin/time -p`
  on the host (includes container start); VEPYR time is the wall time of the
  annotation call measured with `time.perf_counter` inside an already running
  Python process.
- VEP peak memory is the `VmHWM` of the VEP (perl) process inside the container,
  sampled with `docker exec` (every 30 s in the full-genome series, every 5 s in
  step S9). It was not recorded in full-genome repeat 1, which is therefore left out
  of the full-genome means, and on chromosomes 1/22 it was measured for the baseline
  scenario only (S9, 2026-10-25).
- Memory units: the thesis reports peak RSS in GiB. VEPYR records store
  `ru_maxrss` in bytes and a field named `max_rss_gb` that is bytes / 1024^3
  (i.e. GiB); VEP records store VmHWM in kB, with `peak_rss_gb` = kB × 1024 / 10^9
  (decimal GB) and, in S9, also `peak_rss_gib`.
- The September VEPYR plugin series on chr1/chr22 (`5repeats_m5`, `chr*_plugins`) mixed
  two tools: repeat 1 was measured with `vepyr-diffly benchmark-plugins` (time including
  Python interpreter start-up, peak RSS sampled over the process tree in MiB and stored
  as MiB/1000 under `max_rss_gb`), repeats 2–5 with `measure_once.py`. The thesis
  therefore uses the 2026-10-25 series, in which all five repeats use `measure_once.py`
  with the same settings as repeats 2–5 (same vepyr build installed 2026-09-23,
  `workers=1`, same per-plugin cache directories, same order). Its core cache is the one
  rebuilt by step S2 with the same code, and its inputs are the chr1/chr22 slices
  recreated by step S0 with the September command (the September files had been
  deleted). `scripts/plugin_rss_gib.py` reproduces the RSS means of the
  September series. RSS from the worker-scaling series is not used in the thesis (its
  repeat-1 entries in `summary_5repeats.json` do not all match the `*.metrics.json`
  records, which store KiB).
- Step S9 was first attempted without Docker (failed run kept in
  `S9/phase1_no_docker/`), and S7 likewise (`S7/result_phase1_no_docker.json`);
  both were re-run by `scripts/m5_phase2.sh`.
- Annotated VCFs from timing runs were deleted after each run; record counts and
  body checksums were kept.
- The five-plugin references use the VEP 116 plugin code except CADD, which is the
  corrected `CADD.pm` from commit `7a1f6450fd12` of the `mwiewior/VEP_plugins` fork (the release-116 plugin reads
  the HGVS-shifted allele and misses valid annotations); the generator pins the SHA-256 of
  every plugin file. Plugin sources were served from the local `plugin_input/` directory
  with `rclone serve http`. The parity gate script (`verify_parity_gate.py`) covers the
  core CSQ fields only and refuses the `merged_plugins` profile (`gate.log`); concordance
  of that profile rests on the E2E reports and their empty mismatch ledgers.
- All VEPYR runs with plugins (plugin series, all-plugins scaling, E2E `merged_plugins`,
  single-plugin field comparisons, S6, S8) used one plugin cache, built on 2026-08-27 and
  downloaded from the team's data store with the other inputs. Its manifests record neither
  the code version that built it nor the identities of the source files (they predate the
  `url`/`md5` fields of the `vepyr-plugins` manifests), so its agreement with the source files
  used by VEP rests on the plugin field comparisons.
- Full-genome parity gate (S10): the PASS reported in the thesis is recorded in
  `S10/gate_original_dir.log` ("PARITY GATE PASSED … 22 contigs", all counters 0), run on
  the original reports directory. `S10/result.json` shows `gate_exit: 1` and `S10/gate.log`
  two FAILs because the gate was first run on the copied reports, whose mismatch-ledger
  paths no longer matched; all 22 reports have zero differences and empty ledgers.
- Repeat 1 of the chr1/chr22 VEPYR worker-scaling series: use the unrounded times in
  `chapter5-m5-.../benchmarks/vepyr_chr*/workers*.metrics.json`. `summary_5repeats.json`
  stores them rounded to two decimals, which changes chr22 S_2 from 1.4246 (thesis: 1,42)
  to 1.4255.
- Paths in scripts and JSON files refer to the author's machine.
