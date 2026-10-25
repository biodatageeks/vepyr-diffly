# Fast Annotation Report: chr1 (parquet, release 116, profile merged_plugins)

**Date:** 2026-09-24 07:27
**Variants:** 323,430 (HG002 GRCh38, bcftools norm -m -both)
**Backend:** parquet only
**Total annotation time:** 18s (0.3 min)
**Aggregate rate:** 17,869 variants/s
**Build:** branch `master` @ `479d8f46eee2`, bio-functions rev `06fb397d6b47`
**Cargo.lock SHA-256:** `7060c5c616a65ea52a585d6678c1b82248a33f66ff3fd03f74b4bc4576560ade`
**datafusion-bio-format-core:** `00487b9edf53cbf36918d7d22a03a79fd0acf863`; effective source `git+https://github.com/biodatageeks/datafusion-bio-formats.git?tag=v1.13.0#00487b9edf53cbf36918d7d22a03a79fd0acf863`
**datafusion-bio-format-ensembl-cache:** `00487b9edf53cbf36918d7d22a03a79fd0acf863`; effective source `git+https://github.com/biodatageeks/datafusion-bio-formats.git?tag=v1.13.0#00487b9edf53cbf36918d7d22a03a79fd0acf863`
**datafusion-bio-format-vcf:** `00487b9edf53cbf36918d7d22a03a79fd0acf863`; effective source `git+https://github.com/biodatageeks/datafusion-bio-formats.git?tag=v1.13.0#00487b9edf53cbf36918d7d22a03a79fd0acf863`
**datafusion-bio-function-vep:** `06fb397d6b474d594066230c17904f373e257f76`; effective source `git+https://github.com/biodatageeks/datafusion-bio-functions.git?tag=v0.22.1#06fb397d6b474d594066230c17904f373e257f76`
**VEP reference:** VEP `116.0`, API `116`, cache `116`, Ensembl `116.c0cf13d`, variation `116.2fb834b`
**Validated cache:** release `116`, source `merged` (contig-local Parquet metadata)
**Native target:** VEP `116.0`, semantics `116`, cache `116`

## Headline

- **124 / 124 CSQ fields at 100% match** (0 value mismatches)
- **0 fields** with mismatches, **0 total** across CSQ entries
- **49 fields FIXED** to 0 vs previous benchmark (was 2,864 mismatches)
- All mismatches traced to root cause classes with upstream issues filed

## Root Cause Classification & Issue Tracker

| # | Root Cause | Mismatches | Fields Affected | Upstream Issue | Status |
|---|-----------|-----------|-----------------|---------------|--------|
| 1 | `stop_retained_variant` false positive on inframe ops | 0 | Consequence | [#90](https://github.com/biodatageeks/datafusion-bio-functions/issues/90), [#117](https://github.com/biodatageeks/datafusion-bio-functions/issues/117), [#113](https://github.com/biodatageeks/datafusion-bio-functions/pull/113) | FIXED |
| 2 | `stop_gained` extra on frameshift | 0 | Consequence | [#114](https://github.com/biodatageeks/datafusion-bio-functions/issues/114) | FIXED |
| 3 | `stop_lost` false positive | 0 | Consequence | [#115](https://github.com/biodatageeks/datafusion-bio-functions/issues/115) | FIXED |
| 4 | `stop_lost` missing on frameshift past stop codon | 0 | Consequence | [#90](https://github.com/biodatageeks/datafusion-bio-functions/issues/90) | FIXED |
| 5 | Inframe/frameshift disagree at CDS boundary | 0 | Consequence | [#117](https://github.com/biodatageeks/datafusion-bio-functions/issues/117) | FIXED |
| 6 | Incomplete terminal codon: IMPACT/HGVSp residual (Xaa vs Ter, missing p.Ter=) | ~0 IMPACT + 0 HGVSp | IMPACT, HGVSp | [#130](https://github.com/biodatageeks/datafusion-bio-functions/issues/130) | FIXED |
| 7 | `stop_gained` missing on frameshift/inframe_deletion | 0 | Consequence | [#116](https://github.com/biodatageeks/datafusion-bio-functions/issues/116) | FIXED |
| 8 | `incomplete_terminal_codon` companion terms | 0 | Consequence | [#101](https://github.com/biodatageeks/datafusion-bio-functions/issues/101) | FIXED |
| 9 | HGVSc/HGVS_OFFSET on non-coding + UTR indels | ~0 + 0 | HGVSc, HGVS_OFFSET | [#112](https://github.com/biodatageeks/datafusion-bio-functions/issues/112) | FIXED |
| 10 | HGNC_ID false-positive propagation | ~0 | HGNC_ID | [#108](https://github.com/biodatageeks/datafusion-bio-functions/issues/108) | FIXED |
| 11 | CDS/protein fields missing at CDS boundary | ~0 | CDS_position, Protein_position, Amino_acids, Codons, DOMAINS | [#118](https://github.com/biodatageeks/datafusion-bio-functions/issues/118) | FIXED |
| 12 | miRNA dedup (stem repeated in VEP) | 0 | miRNA | [#100](https://github.com/biodatageeks/datafusion-bio-functions/issues/100) | FIXED |
| 13 | `mature_miRNA_variant` emitted outside a mature product | 0 | Consequence | [#90](https://github.com/biodatageeks/datafusion-bio-functions/issues/90) | FIXED |
| 14 | `frameshift_variant` missing at CDS boundary | 0 | Consequence | [#117](https://github.com/biodatageeks/datafusion-bio-functions/issues/117) | FIXED |
| 15 | `protein_altering_variant` not emitted for complex inframe changes | 0 | Consequence | [#124](https://github.com/biodatageeks/datafusion-bio-functions/issues/124) | FIXED |
| 16 | `start_retained_variant` missing alongside `start_lost` | 0 | Consequence | [#125](https://github.com/biodatageeks/datafusion-bio-functions/issues/125) | FIXED |

## Per-Chromosome Performance

| Chrom | Variants | Time (s) | Rate (v/s) |
|-------|----------|----------|------------|
| chr1 | 323,430 | 18.1 | 17,869 |
| **TOTAL** | **323,430** | **18.1** | **17,869** |

## Variant Coverage

| Metric | Value |
|--------|-------|
| Variants compared | 323,430 |
| CSQ entry count match | 323,430 |
| CSQ entry count mismatch | 0 |
| Only in vepyr | 0 |
| Only in VEP | 0 |
| CSQ entries only in vepyr | 0 |
| CSQ entries only in VEP | 0 |
| Uncapped mismatch-ledger rows | 0 |

## Field Equality Shapes

| Shape | CSQ field comparisons |
|-------|----------------------:|
| `both_empty` | 443,066,763 |
| `both_nonempty_equal` | 262,352,745 |
| `vepyr_empty_only` | 0 |
| `vep_empty_only` | 0 |
| `both_nonempty_unequal` | 0 |

## Field-Level Mismatches: NEW vs OLD Benchmark

| Field | NEW (this run) | OLD (benchmark) | Delta | Status |
|-------|---------------|-----------------|-------|--------|
| HGVSc | 0 | 317 | -317 | FIXED |
| Feature | 0 | 260 | -260 | FIXED |
| ENSP | 0 | 242 | -242 | FIXED |
| HGVSp | 0 | 200 | -200 | FIXED |
| INTRON | 0 | 178 | -178 | FIXED |
| Consequence | 0 | 175 | -175 | FIXED |
| UNIPARC | 0 | 173 | -173 | FIXED |
| TREMBL | 0 | 117 | -117 | FIXED |
| cDNA_position | 0 | 91 | -91 | FIXED |
| APPRIS | 0 | 88 | -88 | FIXED |
| CCDS | 0 | 84 | -84 | FIXED |
| TSL | 0 | 83 | -83 | FIXED |
| IMPACT | 0 | 83 | -83 | FIXED |
| HGNC_ID | 0 | 78 | -78 | FIXED |
| EXON | 0 | 65 | -65 | FIXED |
| CDS_position | 0 | 64 | -64 | FIXED |
| Protein_position | 0 | 64 | -64 | FIXED |
| DOMAINS | 0 | 61 | -61 | FIXED |
| UNIPROT_ISOFORM | 0 | 51 | -51 | FIXED |
| FLAGS | 0 | 45 | -45 | FIXED |
| Codons | 0 | 44 | -44 | FIXED |
| BIOTYPE | 0 | 44 | -44 | FIXED |
| Amino_acids | 0 | 44 | -44 | FIXED |
| SWISSPROT | 0 | 43 | -43 | FIXED |
| CANONICAL | 0 | 30 | -30 | FIXED |
| MANE | 0 | 22 | -22 | FIXED |
| MANE_SELECT | 0 | 22 | -22 | FIXED |
| Gene | 0 | 15 | -15 | FIXED |
| SYMBOL | 0 | 15 | -15 | FIXED |
| STRAND | 0 | 14 | -14 | FIXED |
| HGVS_OFFSET | 0 | 11 | -11 | FIXED |
| PolyPhen | 0 | 8 | -8 | FIXED |
| SIFT | 0 | 8 | -8 | FIXED |
| DISTANCE | 0 | 7 | -7 | FIXED |
| SYMBOL_SOURCE | 0 | 4 | -4 | FIXED |
| gnomADg_SAS_AF | 0 | 1 | -1 | FIXED |
| gnomADg_ASJ_AF | 0 | 1 | -1 | FIXED |
| gnomADg_EAS_AF | 0 | 1 | -1 | FIXED |
| gnomADg_MID_AF | 0 | 1 | -1 | FIXED |
| gnomADg_AMI_AF | 0 | 1 | -1 | FIXED |
| gnomADg_AMR_AF | 0 | 1 | -1 | FIXED |
| gnomADg_FIN_AF | 0 | 1 | -1 | FIXED |
| gnomADg_AFR_AF | 0 | 1 | -1 | FIXED |
| MAX_AF_POPS | 0 | 1 | -1 | FIXED |
| gnomADg_REMAINING_AF | 0 | 1 | -1 | FIXED |
| gnomADg_NFE_AF | 0 | 1 | -1 | FIXED |
| gnomADg_AF | 0 | 1 | -1 | FIXED |
| MAX_AF | 0 | 1 | -1 | FIXED |
| miRNA | 0 | 1 | -1 | FIXED |

**Total mismatches: 0** (was 2,864, delta -2,864)

### Fields FIXED (previously had mismatches, now 0): 49 fields

**Consequence** (175), **IMPACT** (83), **SYMBOL** (15), **Gene** (15), **Feature** (260), **BIOTYPE** (44), **EXON** (65), **INTRON** (178), **HGVSc** (317), **HGVSp** (200), **cDNA_position** (91), **CDS_position** (64), **Protein_position** (64), **Amino_acids** (44), **Codons** (44), **DISTANCE** (7), **STRAND** (14), **FLAGS** (45), **SYMBOL_SOURCE** (4), **HGNC_ID** (78), **CANONICAL** (30), **MANE** (22), **MANE_SELECT** (22), **TSL** (83), **APPRIS** (88), **CCDS** (84), **ENSP** (242), **SWISSPROT** (43), **TREMBL** (117), **UNIPARC** (173), **UNIPROT_ISOFORM** (51), **SIFT** (8), **PolyPhen** (8), **DOMAINS** (61), **miRNA** (1), **HGVS_OFFSET** (11), **gnomADg_AF** (1), **gnomADg_AFR_AF** (1), **gnomADg_AMI_AF** (1), **gnomADg_AMR_AF** (1), **gnomADg_ASJ_AF** (1), **gnomADg_EAS_AF** (1), **gnomADg_FIN_AF** (1), **gnomADg_MID_AF** (1), **gnomADg_NFE_AF** (1), **gnomADg_REMAINING_AF** (1), **gnomADg_SAS_AF** (1), **MAX_AF** (1), **MAX_AF_POPS** (1)

## Remaining Mismatch Details

## Per-Chromosome Mismatch Breakdown

| Chrom | Variants | CSQ Match | Consequence | HGVSc | HGVSp | IMPACT | HGNC_ID | Other |
|-------|----------|-----------|-------------|-------|-------|--------|---------|-------|
| chr1 | 323,430 | 323,430 | 0 | 0 | 0 | 0 | 0 | 0 |
