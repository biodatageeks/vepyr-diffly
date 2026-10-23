# Cache build benchmark — 20,21,22

- created: 2026-09-25T02:54:06+0200
- knob swept: `build_concurrency`
- platform: macOS-27.0-arm64-arm-64bit-Mach-O
- git: unknown@unknown
- raw cache: /Users/lukaszjezapkowicz/workspace/data_vepyr/runs/chapter5-m5-20260924T061203Z/benchmarks/cache_lpt/116_GRCh38_merged

| build_concurrency | median s | min s | max s | speedup | RSS GiB | rows |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 220.39 | 220.39 | 220.39 | 1.0x | 21.08 | 67663043 |
| 2 | 168.44 | 168.44 | 168.44 | 1.308x | 20.89 | 67663043 |
| 4 | 136.68 | 136.68 | 136.68 | 1.612x | 20.64 | 67663043 |
| 8 | 169.13 | 169.13 | 169.13 | 1.303x | 20.51 | 67663043 |
