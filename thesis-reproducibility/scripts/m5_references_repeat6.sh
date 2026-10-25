#!/usr/bin/env bash
# Added 2026-10-25 after the thesis audit ("every result from this M5, mean of five runs").
#
# 1. Regenerates the VEP 116 five-plugin references for chr1 and chr22 on this machine with
#    generate_vep_plugin_references.sh (pinned plugin code incl. the corrected CADD.pm, plugin
#    sources served from plugin_input/ over HTTP by rclone, as the script expects). The earlier
#    references (VEP headers 2026-07-27 / 2026-08-26) were downloaded, not built here.
# 2. With Colima stopped: E2E comparison of VEPYR 0.7.0 (479d8f4) against them
#    (run_comparison.py --profile merged_plugins, chr1 and chr22, then the parity gate) and the
#    per-plugin field comparison (vepyr-diffly benchmark-plugins, compare-mode fast). The
#    baseline chr1/chr22 comparison against an M5 reference already exists (S10).
# 3. Starts Colima again and runs a sixth full-genome VEP repeat with the same script and
#    directory as repeats 2-5, so that wall time and peak RSS come from the same five runs (2-6).
#
# Usage: nohup caffeinate -is tools/m5_references_repeat6.sh > runs/m5_references_20261025/run.out 2>&1 &
set -uo pipefail
D=/Users/lukaszjezapkowicz/workspace/data_vepyr
OUT=$D/runs/m5_references_20261025
SCRIPTS=$D/tools/vepyr/e2e-testing/scripts
REPORTS=$D/tools/vepyr/e2e-testing/reports
PYN=$D/.venv-native/bin/python3
TARGET=$D/output/116/plugins
FULL=$D/runs/vep_full_repeats_m5_20261023
FASTA=$D/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa
SOURCE_URL=http://127.0.0.1:18080
mkdir -p "$OUT/input"
LOG=$OUT/run.log
log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$LOG"; }

{
  echo "started_utc: $(date -u +%FT%TZ)"
  sw_vers
  sysctl -n machdep.cpu.brand_string hw.memsize
  grep -E '^(cpu|memory|disk|vmType|mountType):' ~/.colima/default/colima.yaml
  echo "vepyr: $(git -C "$D/tools/vepyr" rev-parse HEAD) $(git -C "$D/tools/vepyr" status --short | wc -l | tr -d ' ') local changes"
  echo "vepyr-diffly: $(git -C "$D/tools/vepyr-diffly" rev-parse HEAD)"
  "$PYN" -c 'import vepyr; print("vepyr python package", vepyr.__version__)'
  shasum -a 256 "$D"/runs/m5_completion_20261024/S7/plugins/*.pm
  df -h "$D" | tail -1
} > "$OUT/environment.txt" 2>&1

# --- 1. five-plugin references --------------------------------------------------------------
log "starting Colima"
colima start >> "$LOG" 2>&1
docker info > /dev/null 2>&1 || { log "docker not available after colima start, aborting"; exit 1; }
echo "image_id: $(docker image inspect --format '{{.Id}}' ensemblorg/ensembl-vep:release_116.0)" >> "$OUT/environment.txt"

rclone serve http "$D/plugin_input" --addr 127.0.0.1:18080 --read-only > "$OUT/rclone.log" 2>&1 &
RCLONE=$!
trap 'kill $RCLONE 2> /dev/null' EXIT
for _ in $(seq 30); do
  curl -sf --head "$SOURCE_URL/clinvar/clinvar.vcf.gz" > /dev/null && break
  sleep 1
done

mkdir -p "$TARGET/plugin_code"
cp -p "$D"/runs/m5_completion_20261024/S7/plugins/*.pm "$TARGET/plugin_code/"

cd "$OUT" || exit 1
# htslib keeps a downloaded remote .tbi in the working directory; fetch each one once here so the
# two parallel builds reuse it instead of downloading the same file concurrently.
for p in clinvar/clinvar.vcf.gz spliceai/spliceai_scores.masked.snv.ensembl_mane.grch38.110.vcf.gz \
         alphamissense/AlphaMissense_hg38.bgz.tsv.gz dbnsfp/dbNSFP5.3.1a_grch38.gz \
         cadd/whole_genome_SNVs.tsv.gz cadd/gnomad.genomes.r4.0.indel.tsv.gz; do
  tabix -l "$SOURCE_URL/$p" > /dev/null || { log "cannot read index of $p over HTTP, aborting"; exit 1; }
done
log "references chr1 chr22: start"
/usr/bin/time -p -o "$OUT/references.time.txt" env DATA_VEPYR_DIR="$D" VEP_PLUGIN_SOURCE_URL="$SOURCE_URL" \
  VEP_NORMALIZED_VCF="$D/input/HG002_normalized.vcf.gz" VEP_REFERENCE_JOBS=2 \
  "$SCRIPTS/generate_vep_plugin_references.sh" 1 22 > "$OUT/references.out" 2>&1
RC=$?
log "references chr1 chr22: exit $RC"
kill $RCLONE 2> /dev/null
[ "$RC" -eq 0 ] || { log "reference generation failed, stopping (see references.out)"; colima stop >> "$LOG" 2>&1; exit 1; }

for c in 1 22; do
  f=$TARGET/HG002_chr${c}_5plugins_vep116_caddfix.vcf.gz
  cp -p "$f.plugins" "$OUT/chr${c}_reference.plugins"
  cp -p "$TARGET/logs/chr${c}.log" "$OUT/chr${c}_reference_build.log"
  bcftools view -h "$f" | grep -E '^##(VEP|VEP-command-line)=' > "$OUT/chr${c}_reference_header.txt"
  bgzip -dc "$f" | grep -v '^#' | tee >(wc -l | tr -d ' ' > "$OUT/chr${c}_reference_records.txt") \
    | shasum -a 256 | cut -d' ' -f1 > "$OUT/chr${c}_reference_body.sha256"
  shasum -a 256 "$f" > "$OUT/chr${c}_reference_file.sha256"
done
log "reference checksums recorded; stopping Colima"
colima stop >> "$LOG" 2>&1

# --- 2. comparisons against the new references (Colima stopped) -----------------------------
mkdir -p "$OUT/reports_backup_before" "$OUT/e2e_reports"
cp -p "$REPORTS"/fast_chr1_merged_plugins_116_* "$REPORTS"/fast_chr22_merged_plugins_116_* "$OUT/reports_backup_before/"
for c in 1 22; do
  log "E2E merged_plugins chr$c: start"
  DATA_VEPYR_DIR="$D" "$PYN" "$SCRIPTS/run_comparison.py" --release 116 --profile merged_plugins --chroms "$c" \
    --workers 4 --bgzf --force >> "$OUT/run_comparison.log" 2>&1
  log "E2E merged_plugins chr$c: exit $?"
done
DATA_VEPYR_DIR="$D" "$PYN" "$SCRIPTS/verify_parity_gate.py" --release 116 --profile merged_plugins --chroms 1 22 \
  > "$OUT/gate.log" 2>&1
log "parity gate merged_plugins chr1 chr22: exit $?"
for c in 1 22; do
  cp -p "$REPORTS/fast_chr${c}_merged_plugins_116_report.json" "$REPORTS/fast_chr${c}_merged_plugins_116_mismatches.jsonl" \
    "$OUT/e2e_reports/"
done

for c in 1 22; do
  bcftools view -Ov -o "$OUT/input/HG002_chr${c}_5plugins_vep116_caddfix.vcf" "$TARGET/HG002_chr${c}_5plugins_vep116_caddfix.vcf.gz"
  bgzip -dc "$D/runs/m5_completion_20261024/inputs/HG002_normalized_chr${c}.vcf.gz" > "$OUT/input/HG002_normalized_chr${c}.vcf"
  log "benchmark-plugins chr$c: start"
  PYTHONPATH=$D/tools/vepyr-diffly/src "$PYN" -m vepyr_diffly.cli benchmark-plugins --preset ensembl_everything \
    --input-vcf "$OUT/input/HG002_normalized_chr${c}.vcf" \
    --reference-vcf "$OUT/input/HG002_chr${c}_5plugins_vep116_caddfix.vcf" \
    --output-dir "$OUT/chr${c}_individual_plugins" --cache-dir "$D/cache/116_GRCh38_merged" \
    --plugin-cache-root "$D/plugin_cache" --reference-fasta "$FASTA" --vepyr-python "$PYN" --compare-mode fast \
    > "$OUT/chr${c}_benchmark_plugins.log" 2>&1
  log "benchmark-plugins chr$c: exit $?"
done
rm -f "$OUT"/input/*.vcf

# --- 3. sixth full-genome VEP repeat --------------------------------------------------------
log "VEP full-genome repeat 6: starting Colima"
colima start >> "$LOG" 2>&1
docker info > /dev/null 2>&1 || { log "docker not available after colima start, aborting"; exit 1; }
cp -p "$FULL/environment.txt" "$FULL/environment_repeats_4_5.txt"
cp -p "$FULL/summary.json" "$FULL/summary_repeats_1_5.json"
FIRST_REPEAT=6 LAST_REPEAT=6 OUT=$FULL "$D/tools/run_vep_full_repeats.sh" >> "$LOG" 2>&1

# Wall time and peak RSS over the same five runs (repeats 2-6, all with RSS sampling).
python3 -I - "$FULL" << 'EOF' | tee -a "$LOG"
import json, statistics, sys
from pathlib import Path

full = Path(sys.argv[1])
runs = [r for r in json.loads((full / "summary.json").read_text())["runs"] if 2 <= r["repeat"] <= 6]
ok = [r for r in runs if r["exit_code"] == 0 and r["wall_s"] and r.get("peak_rss_gb")]
wall = [r["wall_s"] for r in ok]
rss_gib = [r["peak_rss_gb"] * 1e9 / 1024**3 for r in ok]
summary = {
    "repeats": [r["repeat"] for r in ok],
    "n": len(ok),
    "mean_wall_s": round(statistics.mean(wall), 2),
    "min_wall_s": min(wall),
    "max_wall_s": max(wall),
    "rsd_wall_pct": round(100 * statistics.stdev(wall) / statistics.mean(wall), 2),
    "mean_peak_rss_gib": round(statistics.mean(rss_gib), 3),
    "min_peak_rss_gib": round(min(rss_gib), 3),
    "max_peak_rss_gib": round(max(rss_gib), 3),
    "body_sha256_identical": len({r["output_body"] for r in ok}) == 1,
    "output_records": sorted({r["output_records"] for r in ok}),
}
(full / "summary_repeats_2_6.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
EOF
colima stop >> "$LOG" 2>&1
log "done"
