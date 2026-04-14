# GitHub Issue Check: `datafusion-bio-formats#137`

Original issue link: https://github.com/biodatageeks/datafusion-bio-formats/issues/137

Checked on: `2026-04-14`

## Verdict

Issue intent is satisfied in the current local workspace.

Full plugin sources are too large to use as the default proof path, so validation intentionally uses:

- chromosome-scoped builds,
- preview-sized builds,
- fresh clean cache roots,
- schema inspection,
- automated tests,
- round-trip verification,
- runtime annotation smoke verification.

This is enough because it exercises the same codepaths that matter for the issue:

- source loading,
- plugin-specific parsing,
- parquet generation,
- Fjall generation,
- runtime lookup,
- annotation field population.

## Scope Used For Validation

Because full plugin inputs are very large, this check intentionally uses:

- `--chromosomes 1`
- `--preview-rows 1000`

That keeps the validation local and reproducible while still exercising:

- source loading,
- plugin-specific parsing,
- parquet writing,
- Fjall writing,
- current cache layout,
- runtime-facing schema.

## What The Issue Required In Practice

The issue-driven work in this repo and the sibling local checkouts was about introducing and validating a real plugin-cache pipeline for the current 5 plugins:

- `clinvar`
- `spliceai`
- `cadd`
- `alphamissense`
- `dbnsfp`

The concrete expectations that were checked here are:

1. plugin build support exists for all 5 plugins,
2. cache output layout matches the version-root layout instead of the old wrapper-specific layout,
3. generated parquet schema is plugin-specific and consistent with the intended runtime contract,
4. heavy indexed sources can be sliced chromosome-wise using `tabix`,
5. local-source builds work without download-time coupling,
6. automated tests covering the build flow pass,
7. dedicated round-trip verification passes for all 5 plugins,
8. dedicated annotation smoke verification passes for all 5 plugins on a fresh clean cache root.

## Repositories / Files That Implement The Issue

The implementation is spread across the local sibling repos:

- `vepyr`
  - `vepyr/src/vepyr/__init__.py`
  - `vepyr/src/lib.rs`
  - `vepyr/src/plugin_convert.rs`
- `datafusion-bio-functions`
  - `datafusion-bio-functions/datafusion/bio-function-vep/src/plugin_cache_builder.rs`
  - `datafusion-bio-functions/datafusion/bio-function-vep/src/plugin_lookup.rs`
  - `datafusion-bio-functions/datafusion/bio-function-vep/src/plugin.rs`
- `datafusion-bio-formats`
  - `datafusion-bio-formats/datafusion/bio-format-vep-plugin/src/lib.rs`
- orchestration in this repo
  - `scripts/build_chr_cache.py`
  - `scripts/create_plugin_indexes.py`
  - `scripts/plugin_round_trip_test.py`
  - `scripts/plugin_annotation_smoke_test.py`

## Step-By-Step Evidence

### 1. The current cache layout matches the intended version-root contract

Current plugin/core cache layout:

```text
.cache/vepyr_cache/115_GRCh38_vep/
  alphamissense/
  alphamissense.fjall/
  cadd/
  cadd.fjall/
  clinvar/
  clinvar.fjall/
  dbnsfp/
  dbnsfp.fjall/
  spliceai/
  spliceai.fjall/
  exon/
  motif/
  regulatory/
  transcript/
  translation_core/
  translation_sift/
  translation_sift.fjall/
  variation/
  variation.fjall/
```

This matters because the old wrapper-added `parquet/<version>/...` layout was explicitly superseded. The cache is now under:

- `<cache_root>/<version>/<plugin>`
- `<cache_root>/<version>/<plugin>.fjall`

which is the intended issue-aligned structure.

### 1a. A fresh clean reduced-scope cache root was also built successfully

To avoid relying only on the reused local cache root, a fresh clean cache was built under:

```text
/tmp/vepyr-issue-check-cache-y
```

Command used:

```bash
cd /Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly
source .venv/bin/activate
python scripts/build_chr_cache.py \
  --cache-dir /tmp/vepyr-issue-check-cache-y \
  --chromosomes Y \
  --plugins clinvar,spliceai,cadd,alphamissense,dbnsfp \
  --preview-rows 1000 \
  --clean-plugin-output \
  --skip-install
```

Key core output:

```text
variation: chrY 2.6M rows (total: 2.6M)
transcript: chrY 68 rows (total: 68)
exon: chrY 375 rows (total: 375)
translation: chrY core=33 sift=33
regulatory: 0 rows
motif: 0 rows
building core fjall caches
```

This is important because it proves the issue check is not based only on stale cache artifacts from the main local cache root.

### 2. All 5 plugin caches were generated locally on the reduced validation path

Command used:

```bash
cd /Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly
source .venv/bin/activate
python scripts/build_chr_cache.py \
  --only-plugins \
  --plugins clinvar,spliceai,cadd,alphamissense,dbnsfp \
  --chromosomes 1 \
  --preview-rows 1000 \
  --clean-plugin-output \
  --skip-install
```

Key output:

```text
preview fallback for clinvar: tabix index not found for .../plugins/clinvar.vcf.gz
clinvar (preview=1000): .gz=39.7KB, raw=373.9KB, parquet=2.0MB, preview_prep=0.0s, convert=0.2s, total=0.2s

tabix source slice for spliceai: using preview from spliceai_scores.masked.snv.ensembl_mane.grch38.110.vcf.gz
spliceai (preview=1000): .gz=9.1KB, raw=66.2KB, parquet=2.0MB, preview_prep=0.0s, convert=0.2s, total=0.2s

tabix source slice for cadd_snv: using preview from whole_genome_SNVs.tsv.gz
tabix source slice for cadd_indel: using preview from gnomad.genomes.r4.0.indel.tsv.gz
cadd (snv+indel) (preview=1000): .gz=19.9KB, raw=68.6KB, parquet=2.0MB, preview_prep=0.1s, convert=0.2s, total=0.3s

preview fallback for alphamissense: tabix index not found for .../plugins/AlphaMissense_hg38.tsv.gz
alphamissense (preview=1000): .gz=9.0KB, raw=69.0KB, parquet=2.0MB, preview_prep=0.0s, convert=0.2s, total=0.2s

tabix source slice for dbnsfp: using preview from dbNSFP5.3.1a_grch38.gz
dbnsfp (preview=1000): .gz=250.1KB, raw=1.7MB, parquet=2.0MB, preview_prep=0.2s, convert=0.2s, total=0.4s
```

Why this is important:

- all 5 plugins completed successfully,
- the heaviest indexed sources already use `tabix` on the chromosome-scoped path,
- the build script now exposes meaningful timing and size information per plugin.

The same was also re-verified on the fresh clean `chrY` cache root:

```text
clinvar: chrY 675 rows (total: 675)
spliceai: chrY 1000 rows (total: 1000)
cadd: chrY 2000 rows (total: 2000)
alphamissense: chrY 1000 rows (total: 1000)
dbnsfp: chrY 1000 rows (total: 1000)

clinvar (preview=1000): ... total=6.7s
spliceai (preview=1000): ... total=0.4s
cadd (preview=1000): ... total=0.5s
alphamissense (preview=1000): ... total=59.5s
dbnsfp (preview=1000): ... total=0.6s
```

This clean-cache run exercised:

- fresh core parquet generation,
- fresh core Fjall generation,
- fresh plugin parquet generation,
- fresh plugin Fjall generation,
- all five plugins in one brand-new cache root.

### 3. The generated parquet schemas match the intended plugin-specific outputs

Fresh schema inspection of the generated plugin parquet files:

```text
[clinvar] chr1.parquet
chrom,pos,ref,alt,clnsig,clnrevstat,clndn,clnvc,clnvi,af_esp,af_exac,af_tgp

[spliceai] chr1.parquet
chrom,pos,ref,alt,symbol,ds_ag,ds_al,ds_dg,ds_dl,dp_ag,dp_al,dp_dg,dp_dl

[cadd] chr1.parquet
chrom,pos,ref,alt,raw_score,phred_score

[alphamissense] chr1.parquet
chrom,pos,ref,alt,genome,uniprot_id,transcript_id,protein_variant,am_pathogenicity,am_class

[dbnsfp] chr1.parquet
chrom,pos,ref,alt,sift4g_score,sift4g_pred,polyphen2_hdiv_score,polyphen2_hvar_score,lrt_score,lrt_pred,mutationtaster_score,mutationtaster_pred,fathmm_score,fathmm_pred,provean_score,provean_pred,vest4_score,metasvm_score,metasvm_pred,metalr_score,metalr_pred,revel_score,gerp_rs,phylop100way,phylop30way,phastcons100way,phastcons30way,siphy_29way,cadd_raw,cadd_phred
```

Why this matters:

- each plugin writes its own plugin-specific schema,
- `cadd` is exposed as one logical plugin output,
- `dbnsfp` writes the widened issue-aligned field set rather than a placeholder/minimal schema,
- the generated cache is not just “files on disk”; it has the correct runtime-facing columns.

### 4. Indexed heavy sources are actually being used through `tabix`

The reduced build above confirms that the optimized path is active for:

- `spliceai`
- `cadd_snv`
- `cadd_indel`
- `dbnsfp`

Observed output:

```text
tabix source slice for spliceai: using preview from ...
tabix source slice for cadd_snv: using preview from ...
tabix source slice for cadd_indel: using preview from ...
tabix source slice for dbnsfp: using preview from ...
```

This matters because the heavy source files are:

- `dbNSFP5.3.1a_grch38.gz`
- `whole_genome_SNVs.tsv.gz`
- `gnomad.genomes.r4.0.indel.tsv.gz`
- `spliceai_scores.masked.snv.ensembl_mane.grch38.110.vcf.gz`

Without this optimization, chromosome-scoped validation would still degrade into much larger linear scans.

### 5. Automated build-path tests are green

Command used:

```bash
cd /Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly
source .venv/bin/activate
pytest -q tests/test_build_chr_cache_script.py vepyr/tests/test_build_cache_plugins.py
```

Result:

```text
25 passed in 0.07s
```

This is important because it covers the orchestration and plugin-build layer, including the newer cache layout and current script contract.

### 6. The helper script for round-trip validation was found to be stale and was corrected

During verification, `scripts/plugin_round_trip_test.py` still called an old private helper signature:

```text
TypeError: _prepare_preview_source() missing 1 required positional argument: 'temp_paths'
```

That was a real local regression in the helper script, not in plugin-cache generation itself.

It was corrected by:

- passing the new `chromosomes` argument,
- making schema validation tolerant of `string` vs `string_view`, which is expected in current Arrow output.

This matters because the issue check should distinguish:

- actual product/build-path failures,
- from stale validation helper failures.

### 6a. Round-trip validation now passes for all 5 plugins

After correcting the helper, the following command passed:

```bash
cd /Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly
source .venv/bin/activate
python -u scripts/plugin_round_trip_test.py \
  --cache-dir /tmp/vepyr-plugin-roundtrip-cache-y \
  --plugins clinvar,spliceai,cadd,alphamissense,dbnsfp \
  --preview-rows 1000 \
  --skip-install
```

Observed result:

```text
[clinvar] round-trip: PASSED
[spliceai] round-trip: PASSED
[cadd] round-trip: PASSED
[alphamissense] round-trip: PASSED
[dbnsfp] round-trip: PASSED
round-trip check PASSED
```

This is stronger evidence than schema inspection alone because it confirms:

- conversion succeeded,
- parquet files were readable,
- expected row counts matched actual parquet row counts,
- schema validation passed for every plugin.

### 7. A separate annotation smoke attempt exposed a pre-existing corrupted core cache artifact

When validating annotation against the existing local cache, the failure was:

```text
RuntimeError: Stream: Parquet error: Parquet error: Invalid Parquet file. Corrupt footer
```

Follow-up inspection identified corrupted existing core parquet files:

```text
.cache/vepyr_cache/115_GRCh38_vep/transcript/chr1.parquet
ArrowInvalid('Parquet file size is 4 bytes, smaller than the minimum file footer (8 bytes)')

.cache/vepyr_cache/115_GRCh38_vep/variation/chr1.parquet
ArrowInvalid('Parquet magic bytes not found in footer. Either the file is corrupted or this is not a parquet file.')
```

This does **not** invalidate the plugin-cache issue itself, because:

- the plugin parquet/fjall generation succeeded,
- the corruption is in the local core cache artifacts from an earlier interrupted run,
- it is not evidence that the plugin issue implementation is missing.

### 7a. Fresh-cache annotation smoke initially exposed a runtime bug outside plugin-cache generation

After building the clean `chrY` cache root, annotation smoke was re-run against that fresh cache:

```bash
cd /Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly
source .venv/bin/activate
python -u scripts/plugin_annotation_smoke_test.py \
  --cache-root /tmp/vepyr-issue-check-cache-y \
  --plugins clinvar,spliceai,cadd,alphamissense,dbnsfp \
  --skip-install
```

The fresh clean cache root reached the real runtime path, but direct runtime debugging then exposed a concrete failure:

```text
Annotation stream error: Error during planning:
table 'datafusion.public.__vep_partitioned_transcript_chrY' not found
```

This was diagnosed as a runtime issue in ephemeral DataFusion table naming for mixed-case chromosome identifiers like `chrY`, not a plugin-cache generation issue. The clean cache root contained the required parquet files, including:

```text
/tmp/vepyr-issue-check-cache-y/115_GRCh38_vep/transcript/chrY.parquet
```

The runtime fix was applied in:

- `datafusion-bio-functions/datafusion/bio-function-vep/src/partitioned_cache.rs`

by sanitizing per-chrom ephemeral table names to lowercase-safe identifiers before registration.

### 7b. After the runtime fix, direct annotation and the full annotation smoke helper both pass

The runtime fix was rebuilt with:

```bash
cd /Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/vepyr
source ../.venv/bin/activate
maturin develop --release
```

Minimal direct runtime verification on the fresh clean cache root then succeeded for:

- no plugins,
- `dbnsfp` only,
- all 5 plugins.

Representative result:

```text
no_plugins annotate_s 0.199 collect_s 0.751 rows 1
dbnsfp_only annotate_s 0.012 collect_s 0.728 rows 1
all_plugins annotate_s 0.013 collect_s 0.793 rows 1
```

Finally, the full operator helper also passed:

```bash
cd /Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly
source .venv/bin/activate
python -u scripts/plugin_annotation_smoke_test.py \
  --cache-root /tmp/vepyr-issue-check-cache-y \
  --plugins clinvar,spliceai,cadd,alphamissense,dbnsfp \
  --skip-install
```

Observed result:

```text
[clinvar] annotate: variant=chrY:2786976 T>C populated=clnrevstat,clndn,clnvc,clnvi
[spliceai] annotate: variant=chrY:2786856 G>A populated=ds_dg,ds_dl
[cadd] annotate: variant=chrY:10001 C>A populated=raw_score
[alphamissense] annotate: variant=chrY:284191 G>A populated=uniprot_id,transcript_id,protein_variant,am_class
[dbnsfp] annotate: variant=chrY:2786989 C>A populated=polyphen2_hdiv_score
annotation smoke test PASSED
```

This is the strongest current evidence because:

- the cache was freshly built,
- the runtime bug was found and fixed,
- the full operator helper completed successfully,
- all 5 plugins populated real runtime fields from the generated cache.

## Why The Issue Should Be Considered Satisfied

The issue was about having a real plugin-cache implementation path, not just stub plumbing.

That is now demonstrably true because:

1. all 5 target plugins build locally from real sources,
2. the build writes real parquet plus Fjall outputs,
3. the cache layout matches the intended version-root structure,
4. heavy indexed plugin sources are actually sliced with `tabix`,
5. parquet schema reflects the intended plugin field sets,
6. the automated build-path tests are green,
7. the dedicated plugin round-trip validation now passes for all 5 plugins,
8. the dedicated plugin annotation smoke validation now passes for all 5 plugins on a fresh clean cache root.

## What Was Not Claimed

This check does **not** claim the following:

- a full all-chromosome/full-row plugin-cache build was completed in this session,
- the existing local core cache is clean,
- every helper script was already up to date before validation.

Those are separate operational concerns.

## Practical Conclusion

For the scope of issue `datafusion-bio-formats#137`, the local implementation is in place and working:

- plugin cache generation exists,
- it supports all 5 plugins,
- current reduced builds produce correct-looking cache artifacts,
- the cache layout and schema are aligned with the intended design.

The remaining work is operational hardening:

- rebuild/refresh corrupted core cache artifacts before relying on annotation smoke from that older cache root,
- keep helper validation scripts in sync with the evolving build API,
- run larger chromosome/full builds only when time and disk budget allow.
