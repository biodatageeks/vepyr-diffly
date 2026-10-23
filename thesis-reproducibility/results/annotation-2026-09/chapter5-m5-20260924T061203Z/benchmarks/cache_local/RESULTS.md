# Cache build benchmark — 1,22

- created: 2026-09-25T02:41:31+0200
- knob swept: `partitions`
- platform: macOS-27.0-arm64-arm-64bit-Mach-O
- git: main@812c4b14ab
- raw cache: /Users/lukaszjezapkowicz/workspace/data_vepyr/runs/chapter5-m5-20260924T061203Z/benchmarks/cache_local/116_GRCh38_merged

| partitions | median s | min s | max s | speedup | RSS GiB | rows |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 824.34 | 824.34 | 824.34 | 1.0x | 21.86 | 188875139 |
| 2 | 914.98 | 914.98 | 914.98 | 0.901x | 22.26 | 188875139 |
| 4 | 961.95 | 961.95 | 961.95 | 0.857x | 22.04 | 188875139 |
| 8 | 866.53 | 866.53 | 866.53 | 0.951x | 23.08 | 188875139 |
