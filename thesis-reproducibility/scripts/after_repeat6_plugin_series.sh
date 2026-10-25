#!/usr/bin/env bash
# Added 2026-10-25: waits for m5_references_repeat6.sh to finish (VEP full-genome repeat 6,
# Colima stopped afterwards), then re-measures the VEPYR plugin series on chr1/chr22 (5x, all
# repeats with measure_once.py) so that no two measurements overlap.
#
# Usage: nohup caffeinate -is tools/after_repeat6_plugin_series.sh > runs/plugin_series_m5_20261025.out 2>&1 &
set -uo pipefail
D=/Users/lukaszjezapkowicz/workspace/data_vepyr
OUT=$D/runs/plugin_series_m5_20261025
WAIT_LOG=$D/runs/m5_references_20261025/run.log
log() { echo "$(date -u +%FT%TZ) $*"; }

until grep -qE '^[0-9TZ:-]+ done$' "$WAIT_LOG"; do sleep 60; done
log "repeat-6 script finished"
if colima status > /dev/null 2>&1; then
  log "Colima still running, stopping it"
  colima stop
fi
sleep 120
pgrep -fl 'measure_once\.py|vepyr_diffly|run_vep_full_repeats|m5_references_repeat6' && { log "other measurement processes found, aborting"; exit 1; }

mkdir -p "$OUT"
{
  echo "started_utc: $(date -u +%FT%TZ)"
  sw_vers
  sysctl -n machdep.cpu.brand_string hw.memsize
  echo "colima: $(colima status 2>&1 | tail -1)"
  "$D/.venv-native/bin/python3" -c 'import vepyr; print("vepyr python package", vepyr.__version__)'
  ls -lT "$D"/.venv-native/lib/python3.14/site-packages/vepyr-0.7.0.dist-info/RECORD
  echo "vepyr: $(git -C "$D/tools/vepyr" rev-parse HEAD)"
  for c in 1 22; do
    f=$D/runs/m5_completion_20261024/inputs/HG002_normalized_chr$c.vcf.gz
    echo "chr$c input body sha256: $(bcftools view -H "$f" | shasum -a 256 | cut -d' ' -f1)"
  done
  for p in clinvar spliceai cadd alphamissense dbnsfp; do
    echo "plugin_cache/plugin/$p -> $(readlink "$D/runs/chapter5-m5-20260924T061203Z/correctness/chr1_individual_plugins/$p/plugin_cache/plugin/$p")"
  done
  shasum -a 256 "$D/tools/measure_once.py" "$D/tools/rerun_plugin_series_5x.py"
} > "$OUT/environment.txt" 2>&1

log "plugin series: start"
python3 -I "$D/tools/rerun_plugin_series_5x.py" "$OUT" > "$OUT/run.log" 2>&1
log "plugin series: exit $?"
