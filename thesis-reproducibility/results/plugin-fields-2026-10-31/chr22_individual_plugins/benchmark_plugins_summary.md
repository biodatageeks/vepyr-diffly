# vepyr plugin benchmark

- input: `/Users/lukaszjezapkowicz/workspace/data_vepyr/runs/plugin_fields_fixed_20261031/input/HG002_normalized_chr22.vcf`
- VEP reference: `/Users/lukaszjezapkowicz/workspace/data_vepyr/runs/plugin_fields_fixed_20261031/input/HG002_chr22_5plugins_vep116_caddfix.vcf`
- variants: 50861 (SNV: 42845, indel/other: 8016)

## Summary

| config | wall (s) | cpu (s) | peak RSS (MB) | throughput (var/s) | cache size | output size (Δ base) | plugin fields equal | mismatches (L/R/U) |
|---|---|---|---|---|---|---|---|---|
| base | 4.5 | 4.5 | 1120 | 11305.0 | 0.0 MB | 494.7 MB (—) | n/a | n/a |
| clinvar | 3.5 | 4.4 | 1179 | 14699.7 | 77.3 MB | 507.4 MB (+12.7 MB) | True | 0/0/0 |
| spliceai | 5.0 | 5.9 | 1272 | 10268.7 | 24298.3 MB | 523.8 MB (+29.1 MB) | True | 0/0/0 |
| cadd | 5.1 | 6.1 | 1315 | 9990.4 | 70841.7 MB | 512.9 MB (+18.2 MB) | True | 0/0/0 |
| alphamissense | 3.3 | 4.5 | 1102 | 15264.4 | 544.7 MB | 497.1 MB (+2.4 MB) | True | 0/0/0 |
| dbnsfp | 3.7 | 4.6 | 1190 | 13768.5 | 3843.6 MB | 520.0 MB (+25.3 MB) | True | 0/0/0 |

## Per-config detail

### base

**Consequence types (top 10):**

| type | count |
|---|---|
| intron_variant | 918255 |
| non_coding_transcript_variant | 218655 |
| upstream_gene_variant | 135827 |
| downstream_gene_variant | 130749 |
| NMD_transcript_variant | 121254 |
| 3_prime_UTR_variant | 20222 |
| non_coding_transcript_exon_variant | 11707 |
| intergenic_variant | 6491 |
| 5_prime_UTR_variant | 5322 |
| synonymous_variant | 4892 |

### clinvar

**Consequence types (top 10):**

| type | count |
|---|---|
| intron_variant | 918255 |
| non_coding_transcript_variant | 218655 |
| upstream_gene_variant | 135827 |
| downstream_gene_variant | 130749 |
| NMD_transcript_variant | 121254 |
| 3_prime_UTR_variant | 20222 |
| non_coding_transcript_exon_variant | 11707 |
| intergenic_variant | 6491 |
| 5_prime_UTR_variant | 5322 |
| synonymous_variant | 4892 |

**ClinVar_CLNSIG** (categorical, n=977, top values):

| value | count |
|---|---|
| Benign | 934 |
| Uncertain_significance | 10 |
| drug_response | 10 |
| Likely_benign | 9 |
| Benign/Likely_benign | 5 |
| risk_factor | 2 |
| Likely_benign&drug_response&other | 2 |
| Conflicting_classifications_of_pathogenicity | 2 |
| Benign&other | 1 |
| other | 1 |

### spliceai

**Consequence types (top 10):**

| type | count |
|---|---|
| intron_variant | 918255 |
| non_coding_transcript_variant | 218655 |
| upstream_gene_variant | 135827 |
| downstream_gene_variant | 130749 |
| NMD_transcript_variant | 121254 |
| 3_prime_UTR_variant | 20222 |
| non_coding_transcript_exon_variant | 11707 |
| intergenic_variant | 6491 |
| 5_prime_UTR_variant | 5322 |
| synonymous_variant | 4892 |

**SpliceAI_pred_DS_AG** (numeric, n=18200, min=0.0, max=0.96, mean=0.0014):

| bucket | count |
|---|---|
| [0.0, 0.192] | 18180 |
| [0.192, 0.384] | 14 |
| [0.384, 0.576] | 3 |
| [0.576, 0.768] | 2 |
| [0.768, 0.96] | 1 |

### cadd

**Consequence types (top 10):**

| type | count |
|---|---|
| intron_variant | 918255 |
| non_coding_transcript_variant | 218655 |
| upstream_gene_variant | 135827 |
| downstream_gene_variant | 130749 |
| NMD_transcript_variant | 121254 |
| 3_prime_UTR_variant | 20222 |
| non_coding_transcript_exon_variant | 11707 |
| intergenic_variant | 6491 |
| 5_prime_UTR_variant | 5322 |
| synonymous_variant | 4892 |

**CADD_PHRED** (numeric, n=50461, min=0.001, max=36.0, mean=2.3415):

| bucket | count |
|---|---|
| [0.001, 7.2008] | 47212 |
| [7.2008, 14.4006] | 2721 |
| [14.4006, 21.6004] | 490 |
| [21.6004, 28.8002] | 32 |
| [28.8002, 36.0] | 6 |

### alphamissense

**Consequence types (top 10):**

| type | count |
|---|---|
| intron_variant | 918255 |
| non_coding_transcript_variant | 218655 |
| upstream_gene_variant | 135827 |
| downstream_gene_variant | 130749 |
| NMD_transcript_variant | 121254 |
| 3_prime_UTR_variant | 20222 |
| non_coding_transcript_exon_variant | 11707 |
| intergenic_variant | 6491 |
| 5_prime_UTR_variant | 5322 |
| synonymous_variant | 4892 |

**am_class** (categorical, n=139, top values):

| value | count |
|---|---|
| likely_benign | 133 |
| ambiguous | 4 |
| likely_pathogenic | 2 |

### dbnsfp

**Consequence types (top 10):**

| type | count |
|---|---|
| intron_variant | 918255 |
| non_coding_transcript_variant | 218655 |
| upstream_gene_variant | 135827 |
| downstream_gene_variant | 130749 |
| NMD_transcript_variant | 121254 |
| 3_prime_UTR_variant | 20222 |
| non_coding_transcript_exon_variant | 11707 |
| intergenic_variant | 6491 |
| 5_prime_UTR_variant | 5322 |
| synonymous_variant | 4892 |

**REVEL_score** (numeric, n=355, min=0.007, max=0.765, mean=0.1391):

| bucket | count |
|---|---|
| [0.007, 0.1586] | 246 |
| [0.1586, 0.3102] | 63 |
| [0.3102, 0.4618] | 35 |
| [0.4618, 0.6134] | 8 |
| [0.6134, 0.765] | 3 |

