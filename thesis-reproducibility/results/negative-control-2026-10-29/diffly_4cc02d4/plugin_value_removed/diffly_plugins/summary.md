# Run Summary

- Preset: `ensembl_everything`
- Input VCF: `/Users/lukaszjezapkowicz/workspace/data_vepyr/runs/thesis_extra_20261029/B1/base_chr22_all_plugins.vcf`
- Left annotated VCF: `base_chr22_all_plugins.vcf`
- Right annotated VCF: `out/plugin_value_removed/modified.vcf`
- Sample first N: `full`
- Requested chromosomes: `all`
- Effective chromosomes: `all`
- Progress log: `/Users/lukaszjezapkowicz/workspace/data_vepyr/runs/thesis_extra_20261029/B1/out/plugin_value_removed/diffly_plugins/runtime/compare.progress.log`
- Compare mode: `fast`
- Memory budget MB: `auto`
- Fingerprint only: `False`

## Resource Plan

```json
{
  "memory_budget_mb": 2048,
  "bucket_count": 8,
  "variant_chunk_rows": 50000,
  "consequence_chunk_rows": 16384,
  "compare_workers": 4,
  "parallelize_sides": false,
  "use_bucketed_variant": true,
  "use_bucketed_consequence": true
}
```

## Timings

```text
variant_summary_seconds=4.485s
consequence_bucketization_seconds=5.323s
schema_validation_seconds=0.004s
variant_diff_seconds=0.081s
consequence_diff_seconds=1.841s
```

## Variant Tier

```text
Bucketized diffly variant comparison
bucket_count=8
buckets_compared=8
compare_mode=fast
fingerprint_only=False
precheck_equal_buckets=8
exact_diff_run_buckets=0
left_rows=50861
right_rows=50861
joined_equal_rows=50861
left_only_rows=0
right_only_rows=0
unequal_rows=0
```

## Consequence Tier

```text
Bucketized diffly consequence comparison
workers=4
bucket_count=8
buckets_compared=8
compare_mode=fast
fingerprint_only=False
precheck_equal_buckets=0
exact_diff_run_buckets=8
mismatch_tsv_row_limit=2000
left_rows=59043
right_rows=59043
joined_equal_rows=59007
left_only_rows=36
right_only_rows=36
unequal_rows=0
```

## Per Chromosome

```json
{
  "per_chromosome": {
    "22": {
      "variant": {
        "left_rows": 50861,
        "right_rows": 50861,
        "joined_equal_rows": 50861,
        "left_only_rows": 0,
        "right_only_rows": 0,
        "unequal_rows": 0,
        "equal": true
      },
      "consequence": {
        "left_rows": 59043,
        "right_rows": 59043,
        "joined_equal_rows": 59007,
        "left_only_rows": 36,
        "right_only_rows": 36,
        "unequal_rows": 0,
        "equal": false
      },
      "timings": {
        "variant_summary_seconds": 4.485,
        "consequence_bucketization_seconds": 5.323,
        "consequence_summary_seconds": 0.0,
        "variant_diff_seconds": 0.081,
        "consequence_diff_seconds": 1.841
      }
    }
  },
  "requested": [],
  "effective": []
}
```
