#!/usr/bin/env bash
# Single-plugin field comparison repeated with the corrected vepyr-diffly (2026-10-31).
#
# The comparison in tab:plugin-field-results (results/references-2026-10-25) ran vepyr-diffly
# 4cc02d4, whose CSQ allele rule dropped the entries of indels whose VEP allele trims more than
# the first shared base (362 variants on chr22, 2545 on chr1; negative control, B1). The VEP
# five-plugin references of 2026-10-25 are still on disk (output/116/plugins, body SHA-256 as in
# references-2026-10-25), so this repeats the 2026-10-25 command unchanged except for the
# vepyr-diffly code: branch fix/csq-allele-common-prefix (1d657ce), same VEPYR (.venv-native,
# 0.7.0 @ 479d8f4), same core and plugin caches, same inputs, compare-mode fast.
#
# Usage: rerun_benchmark_plugins_fixed.sh OUT_DIR DIFFLY_FIX_CHECKOUT
set -uo pipefail
D=/Users/lukaszjezapkowicz/workspace/data_vepyr
OUT=$1
FIX=$2
PYN=$D/.venv-native/bin/python3
FASTA=$D/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa
mkdir -p "$OUT/input"
echo "vepyr-diffly: $(git -C "$FIX" rev-parse HEAD) $(git -C "$FIX" status --short | wc -l | tr -d ' ') local changes" > "$OUT/environment.txt"
"$PYN" -c 'import vepyr; print("vepyr", vepyr.__version__)' >> "$OUT/environment.txt"
for c in 22 1; do
  ref=$D/output/116/plugins/HG002_chr${c}_5plugins_vep116_caddfix.vcf.gz
  echo "chr$c reference body sha256 $(bgzip -dc "$ref" | grep -v '^#' | shasum -a 256 | cut -d' ' -f1)" >> "$OUT/environment.txt"
  bcftools view -Ov -o "$OUT/input/HG002_chr${c}_5plugins_vep116_caddfix.vcf" "$ref"
  bgzip -dc "$D/runs/m5_completion_20261024/inputs/HG002_normalized_chr${c}.vcf.gz" > "$OUT/input/HG002_normalized_chr${c}.vcf"
  echo "$(date -u +%FT%TZ) benchmark-plugins chr$c: start" >> "$OUT/run.log"
  PYTHONPATH=$FIX/src "$PYN" -m vepyr_diffly.cli benchmark-plugins --preset ensembl_everything \
    --input-vcf "$OUT/input/HG002_normalized_chr${c}.vcf" \
    --reference-vcf "$OUT/input/HG002_chr${c}_5plugins_vep116_caddfix.vcf" \
    --output-dir "$OUT/chr${c}_individual_plugins" --cache-dir "$D/cache/116_GRCh38_merged" \
    --plugin-cache-root "$D/plugin_cache" --reference-fasta "$FASTA" --vepyr-python "$PYN" --compare-mode fast \
    > "$OUT/chr${c}_benchmark_plugins.log" 2>&1
  echo "$(date -u +%FT%TZ) benchmark-plugins chr$c: exit $?" >> "$OUT/run.log"
  rm -f "$OUT"/input/*.vcf
done
