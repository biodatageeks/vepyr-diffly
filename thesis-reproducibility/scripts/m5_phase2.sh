#!/usr/bin/env bash
# Phase 2 of the M5 completion run (added 2026-10-24 after the VEP series stopped
# at repeat 4 for lack of disk space caused by host swap from the 28 GiB Colima VM).
#
# Waits until phase 1 (m5_completion_batch.py) has finished — it runs the
# non-Docker steps with Colima stopped, so VEPYR timings are not disturbed by
# swapping — then starts Colima with its saved configuration (10 vCPU, 28 GiB)
# and runs, sequentially: the Docker steps S3 and S7, VEP full-genome repeats
# 4 and 5 (same script and output directory as repeats 2–3), and S9.
#
# Usage: nohup caffeinate -is tools/m5_phase2.sh > runs/m5_completion_20261024/phase2.out 2>&1 &
set -uo pipefail
D=/Users/lukaszjezapkowicz/workspace/data_vepyr
OUT=$D/runs/m5_completion_20261024
LOG=$OUT/phase2.log
log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$LOG"; }

log "phase 2 waiting for phase 1 to finish"
until python3 -I -c "import json,sys; sys.exit(0 if json.load(open('$OUT/status.json')).get('phase')=='finished' else 1)"; do
  sleep 120
done
log "phase 1 finished"

# S10 (full-genome E2E) failed in phase 1 on a contig-range argument; rerun it before
# Colima starts, so it does not compete with the VM for memory.
log "S10 rerun (Colima stopped)"
python3 -I "$D/tools/m5_completion_batch.py" --only S10 >> "$OUT/phase2_batch.out" 2>&1

# S7 was marked done in phase 1 although every VEP call failed without Docker: reset it.
python3 -I - "$OUT" <<'PY'
import json, shutil, sys
from pathlib import Path
out = Path(sys.argv[1]); st = json.loads((out / "status.json").read_text())
for s in st["steps"]:
    if s["id"] == "S7":
        s.update(state="pending", detail="reset in phase 2: phase-1 run had no Docker")
(out / "status.json").write_text(json.dumps(st, indent=2, ensure_ascii=False))
r = out / "S7/result.json"
if r.exists():
    shutil.move(r, out / "S7/result_phase1_no_docker.json")
PY

log "starting Colima"
colima start >> "$LOG" 2>&1
docker info > /dev/null 2>&1 || { log "docker not available after colima start, aborting"; exit 1; }

log "S3 and S7 (Docker steps of phase 1)"
python3 -I "$D/tools/m5_completion_batch.py" --only S3 S7 >> "$OUT/phase2_batch.out" 2>&1

log "VEP full-genome repeats 4-5"
cp -p "$D/runs/vep_full_repeats_m5_20261023/environment.txt" "$D/runs/vep_full_repeats_m5_20261023/environment_repeats_2_3.txt"
FIRST_REPEAT=4 LAST_REPEAT=5 OUT=$D/runs/vep_full_repeats_m5_20261023 "$D/tools/run_vep_full_repeats.sh" >> "$LOG" 2>&1

log "S9 (VEP peak RSS on chr22 and chr1)"
python3 -I "$D/tools/m5_completion_batch.py" --only S9 >> "$OUT/phase2_batch.out" 2>&1

log "phase 2 done"
