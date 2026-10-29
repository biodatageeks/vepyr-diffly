# Run Summary

- Preset: `ensembl_everything`
- Input VCF: `/Users/lukaszjezapkowicz/workspace/data_vepyr/runs/thesis_extra_20261029/B1/base_chr22_all_plugins.vcf`
- Left annotated VCF: `base_chr22_all_plugins.vcf`
- Right annotated VCF: `out_fixed/csq_entry_added/modified.vcf`
- Sample first N: `full`
- Requested chromosomes: `all`
- Effective chromosomes: `all`
- Progress log: `/Users/lukaszjezapkowicz/workspace/data_vepyr/runs/thesis_extra_20261029/B1/out_fixed/csq_entry_added/diffly_full/runtime/compare.progress.log`
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
variant_summary_seconds=4.455s
consequence_bucketization_seconds=14.025s
schema_validation_seconds=0.005s
variant_diff_seconds=0.130s
consequence_diff_seconds=73.166s
```

## Variant Tier

```text
Bucketized diffly variant comparison
bucket_count=8
buckets_compared=8
compare_mode=fast
fingerprint_only=False
precheck_equal_buckets=0
exact_diff_run_buckets=8
left_rows=50861
right_rows=50861
joined_equal_rows=50841
left_only_rows=0
right_only_rows=0
unequal_rows=20
```

## Consequence Tier

```text
Bucketized diffly consequence comparison
workers=4
bucket_count=8
buckets_compared=8
compare_mode=fast
fingerprint_only=False
precheck_equal_buckets=1
exact_diff_run_buckets=7
mismatch_tsv_row_limit=2000
left_rows=1242828
right_rows=1242848
joined_equal_rows=1242828
left_only_rows=0
right_only_rows=20
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
        "joined_equal_rows": 50841,
        "left_only_rows": 0,
        "right_only_rows": 0,
        "unequal_rows": 20,
        "equal": false
      },
      "consequence": {
        "left_rows": 1242828,
        "right_rows": 1242848,
        "joined_equal_rows": 1242828,
        "left_only_rows": 0,
        "right_only_rows": 20,
        "unequal_rows": 0,
        "equal": false
      },
      "timings": {
        "variant_summary_seconds": 4.455,
        "consequence_bucketization_seconds": 14.025,
        "consequence_summary_seconds": 0.0,
        "variant_diff_seconds": 0.13,
        "consequence_diff_seconds": 73.166
      }
    }
  },
  "requested": [],
  "effective": []
}
```
