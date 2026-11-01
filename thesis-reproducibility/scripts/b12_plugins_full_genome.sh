#!/usr/bin/env bash
# Plugin scenarios on chromosomes 1-22 (thesis item B12, 2026-10-31).
#
# 1. VEP 116 five-plugin references for the autosomes not built yet (chr2-chr21; chr1 and chr22
#    exist since 2026-10-25 and are skipped by the generator's validation), built exactly as on
#    2026-10-25: generate_vep_plugin_references.sh, pinned plugin code incl. the corrected CADD.pm,
#    plugin sources served from plugin_input/ over HTTP by rclone, one VEP process per chromosome
#    (no --fork), VEP_REFERENCE_JOBS chromosomes at a time.
# 2. E2E comparison of VEPYR 0.7.0 (479d8f4) with all five plugins against them for chr1-22
#    (run_comparison.py --profile merged_plugins, plugin cache = the copy used in all measurements),
#    as on 2026-10-25 for chr1 and chr22. The profile has one reference file per chromosome, so
#    run_comparison.py is called once per chromosome.
# The first run (2026-11-01) built the references and stopped at the E2E step, which then passed all
# 22 chromosomes in one call; the E2E step was rerun on its own with B12_SKIP_REFERENCES=1.
# Usage: [B12_SKIP_REFERENCES=1] b12_plugins_full_genome.sh OUT_DIR [JOBS, default 8]
set -uo pipefail
D=/Users/lukaszjezapkowicz/workspace/data_vepyr
OUT=$1; JOBS=${2:-8}
SCRIPTS=$D/tools/vepyr/e2e-testing/scripts
REPORTS=$D/tools/vepyr/e2e-testing/reports
PYN=$D/.venv-native/bin/python3
mkdir -p "$OUT/e2e_reports"
LOG=$OUT/run.log
log() { echo "$(date -u +%FT%TZ) $* $(pmset -g batt | sed -n 1p)" | tee -a "$LOG"; }
{
  echo "started_utc: $(date -u +%FT%TZ)"
  sw_vers; grep -E '^(cpu|memory|disk|vmType|mountType):' ~/.colima/default/colima.yaml
  echo "vepyr: $(git -C "$D/tools/vepyr" rev-parse HEAD) $(git -C "$D/tools/vepyr" status --short | wc -l | tr -d ' ') local changes"
  "$PYN" -c 'import vepyr; print("vepyr python package", vepyr.__version__)'
  docker image inspect --format '{{.Id}}' ensemblorg/ensembl-vep:release_116.0
  df -h "$D" | tail -1
} >> "$OUT/environment.txt" 2>&1

if [ "${B12_SKIP_REFERENCES:-0}" != 1 ]; then
rclone serve http "$D/plugin_input" --addr 127.0.0.1:18080 --read-only > "$OUT/rclone.log" 2>&1 &
RCLONE=$!
sleep 3
log "references chr1-22 (existing ones skipped): start, jobs=$JOBS"
/usr/bin/time -p -o "$OUT/references.time.txt" env VEP_PLUGIN_SOURCE_URL=http://127.0.0.1:18080 VEP_REFERENCE_JOBS=$JOBS \
  "$SCRIPTS/generate_vep_plugin_references.sh" $(seq 1 22) > "$OUT/references.out" 2>&1
log "references: exit $?"
kill $RCLONE

for c in $(seq 1 22); do
  f=$D/output/116/plugins/HG002_chr${c}_5plugins_vep116_caddfix.vcf.gz
  echo "chr$c $(bgzip -dc "$f" | grep -v '^#' | tee >(wc -l | tr -d ' ' > "$OUT/.n") | shasum -a 256 | cut -d' ' -f1) $(cat "$OUT/.n")" >> "$OUT/reference_bodies.txt"
done
rm -f "$OUT/.n"
fi

mkdir -p "$OUT/reports_backup_before"
cp -p "$REPORTS"/fast_chr*_merged_plugins_116_* "$OUT/reports_backup_before/" 2>/dev/null
log "E2E merged_plugins chr1-22: start"
for c in $(seq 1 22); do
  echo "=== chr$c $(date -u +%FT%TZ)" >> "$OUT/run_comparison.log"
  DATA_VEPYR_DIR="$D" "$PYN" "$SCRIPTS/run_comparison.py" --release 116 --profile merged_plugins --chroms $c \
    --workers 4 --bgzf --force >> "$OUT/run_comparison.log" 2>&1
  log "E2E chr$c: exit $?"
done
DATA_VEPYR_DIR="$D" "$PYN" "$SCRIPTS/verify_parity_gate.py" --release 116 --profile merged_plugins --chroms $(seq 1 22) \
  > "$OUT/gate.log" 2>&1
log "parity gate merged_plugins chr1-22: exit $?"
for c in $(seq 1 22); do
  cp -p "$REPORTS/fast_chr${c}_merged_plugins_116_report.json" "$REPORTS/fast_chr${c}_merged_plugins_116_mismatches.jsonl" "$OUT/e2e_reports/" 2>/dev/null
done
log "done"
