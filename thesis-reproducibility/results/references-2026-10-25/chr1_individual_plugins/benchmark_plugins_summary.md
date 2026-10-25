# vepyr plugin benchmark

- input: `/Users/lukaszjezapkowicz/workspace/data_vepyr/runs/m5_references_20261025/input/HG002_normalized_chr1.vcf`
- VEP reference: `/Users/lukaszjezapkowicz/workspace/data_vepyr/runs/m5_references_20261025/input/HG002_chr1_5plugins_vep116_caddfix.vcf`
- variants: 323430 (SNV: 271740, indel/other: 51690)

## Summary

| config | wall (s) | cpu (s) | peak RSS (MB) | throughput (var/s) | cache size | output size (Δ base) | plugin fields equal | mismatches (L/R/U) |
|---|---|---|---|---|---|---|---|---|
| base | 18.3 | 29.6 | 3696 | 17678.6 | 0.0 MB | 2302.4 MB (—) | n/a | n/a |
| clinvar | 19.2 | 30.5 | 3569 | 16816.4 | 77.3 MB | 2350.7 MB (+48.2 MB) | True | 0/0/0 |
| spliceai | 27.7 | 38.1 | 3731 | 11668.2 | 24298.3 MB | 2435.6 MB (+133.2 MB) | True | 0/0/0 |
| cadd | 30.5 | 40.4 | 3701 | 10587.9 | 70841.7 MB | 2385.4 MB (+83.0 MB) | True | 0/0/0 |
| alphamissense | 19.2 | 29.9 | 3593 | 16820.8 | 544.7 MB | 2313.4 MB (+11.0 MB) | True | 0/0/0 |
| dbnsfp | 20.1 | 30.7 | 3618 | 16069.5 | 3843.6 MB | 2413.0 MB (+110.6 MB) | True | 0/0/0 |

## Per-config detail

### base

**Consequence types (top 10):**

| type | count |
|---|---|
| intron_variant | 4446360 |
| non_coding_transcript_variant | 970287 |
| NMD_transcript_variant | 518556 |
| upstream_gene_variant | 495927 |
| downstream_gene_variant | 492922 |
| 3_prime_UTR_variant | 81574 |
| intergenic_variant | 56088 |
| non_coding_transcript_exon_variant | 39087 |
| synonymous_variant | 20254 |
| 5_prime_UTR_variant | 19230 |

### clinvar

**Consequence types (top 10):**

| type | count |
|---|---|
| intron_variant | 4446360 |
| non_coding_transcript_variant | 970287 |
| NMD_transcript_variant | 518556 |
| upstream_gene_variant | 495927 |
| downstream_gene_variant | 492922 |
| 3_prime_UTR_variant | 81574 |
| intergenic_variant | 56088 |
| non_coding_transcript_exon_variant | 39087 |
| synonymous_variant | 20254 |
| 5_prime_UTR_variant | 19230 |

**ClinVar_CLNSIG** (categorical, n=3464, top values):

| value | count |
|---|---|
| Benign | 3255 |
| Benign/Likely_benign | 72 |
| Likely_benign | 55 |
| Uncertain_significance | 27 |
| Conflicting_classifications_of_pathogenicity | 19 |
| association | 9 |
| not_provided | 8 |
| risk_factor | 6 |
| drug_response | 4 |
| Benign&association | 4 |

### spliceai

**Consequence types (top 10):**

| type | count |
|---|---|
| intron_variant | 4446360 |
| non_coding_transcript_variant | 970287 |
| NMD_transcript_variant | 518556 |
| upstream_gene_variant | 495927 |
| downstream_gene_variant | 492922 |
| 3_prime_UTR_variant | 81574 |
| intergenic_variant | 56088 |
| non_coding_transcript_exon_variant | 39087 |
| synonymous_variant | 20254 |
| 5_prime_UTR_variant | 19230 |

**SpliceAI_pred_DS_AG** (numeric, n=101440, min=0.0, max=0.99, mean=0.001):

| bucket | count |
|---|---|
| [0.0, 0.198] | 101380 |
| [0.198, 0.396] | 49 |
| [0.396, 0.594] | 7 |
| [0.594, 0.792] | 1 |
| [0.792, 0.99] | 3 |

### cadd

**Consequence types (top 10):**

| type | count |
|---|---|
| intron_variant | 4446360 |
| non_coding_transcript_variant | 970287 |
| NMD_transcript_variant | 518556 |
| upstream_gene_variant | 495927 |
| downstream_gene_variant | 492922 |
| 3_prime_UTR_variant | 81574 |
| intergenic_variant | 56088 |
| non_coding_transcript_exon_variant | 39087 |
| synonymous_variant | 20254 |
| 5_prime_UTR_variant | 19230 |

**CADD_PHRED** (numeric, n=320746, min=0.001, max=40.0, mean=2.7001):

| bucket | count |
|---|---|
| [0.001, 8.0008] | 297397 |
| [8.0008, 16.0006] | 20498 |
| [16.0006, 24.0004] | 2762 |
| [24.0004, 32.0002] | 82 |
| [32.0002, 40.0] | 7 |

### alphamissense

**Consequence types (top 10):**

| type | count |
|---|---|
| intron_variant | 4446360 |
| non_coding_transcript_variant | 970287 |
| NMD_transcript_variant | 518556 |
| upstream_gene_variant | 495927 |
| downstream_gene_variant | 492922 |
| 3_prime_UTR_variant | 81574 |
| intergenic_variant | 56088 |
| non_coding_transcript_exon_variant | 39087 |
| synonymous_variant | 20254 |
| 5_prime_UTR_variant | 19230 |

**am_class** (categorical, n=626, top values):

| value | count |
|---|---|
| likely_benign | 589 |
| likely_pathogenic | 19 |
| ambiguous | 18 |

### dbnsfp

**Consequence types (top 10):**

| type | count |
|---|---|
| intron_variant | 4446360 |
| non_coding_transcript_variant | 970287 |
| NMD_transcript_variant | 518556 |
| upstream_gene_variant | 495927 |
| downstream_gene_variant | 492922 |
| 3_prime_UTR_variant | 81574 |
| intergenic_variant | 56088 |
| non_coding_transcript_exon_variant | 39087 |
| synonymous_variant | 20254 |
| 5_prime_UTR_variant | 19230 |

**REVEL_score** (numeric, n=1574, min=0.001, max=0.842, mean=0.1123):

| bucket | count |
|---|---|
| [0.001, 0.1692] | 1254 |
| [0.1692, 0.3374] | 264 |
| [0.3374, 0.5056] | 34 |
| [0.5056, 0.6738] | 16 |
| [0.6738, 0.842] | 6 |

