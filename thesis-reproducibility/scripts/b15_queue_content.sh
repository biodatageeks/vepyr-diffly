#!/usr/bin/env bash
# Content of the queue-built cache against the single-process build (thesis item B15, 2026-11-01).
#
# 1. One full dynamic-queue build (E5: full_queue_e5.py on vepyr branch thesis/cache-build-full-queue,
#    10 processes, 20 groups, six raw entities incl. translation, per-task outputs merged into one
#    cache root), run with the .venv-lpt interpreter (the queue's vepyr build, a wheel built from
#    tools/vepyr-lpt on branch thesis/cache-build-parallel).
# 2. compare_cache_content.py: per (entity, contig) Parquet file of the per-task outputs, row count
#    and order-independent row digests against the core cache built by VEPYR 0.7.0 in one process
#    (cache/116_GRCh38_merged, S2, used in all measurements).
# 3. Annotation of the HG002 chr1 slice (baseline, workers=1) with the merged queue cache, using
#    VEPYR 0.7.0 (.venv-native); the body SHA-256 must equal the native chr1 baseline body.
# The build time is recorded by full_queue_e5.py but is not a measurement of the thesis.
#
# The first attempt (2026-11-01 15:58 UTC) ran with the Colima VM up (28 of 32 GiB) and little free
# disk; the host ran out of memory and swap and rebooted (watchdog panic). Run with the VM stopped.
# A guard samples swap and memory every 10 s into host.log and stops the build (whole process tree)
# when used swap exceeds SWAP_LIMIT_MB (default 8192).
# Usage: [SWAP_LIMIT_MB=8192] b15_queue_content.sh OUT_DIR VEPYR_E5_CHECKOUT
set -uo pipefail
D=/Users/lukaszjezapkowicz/workspace/data_vepyr
OUT=$1; E5=$2; LIMIT=${SWAP_LIMIT_MB:-8192}
PKG=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$OUT"
log() { echo "$(date -u +%FT%TZ) $* $(pmset -g batt | sed -n 1p)" | tee -a "$OUT/run.log"; }
swap_mb() { sysctl -n vm.swapusage | awk '{print $6}' | tr -d M | cut -d. -f1; }
tree() { local p; for p in $(pgrep -P "$1"); do tree "$p"; done; echo "$1"; }
guard() {  # $1 = pid to watch
  while kill -0 "$1" 2>/dev/null; do
    s=$(swap_mb)
    echo "$(date -u +%FT%TZ) swap_used_mb=$s $(memory_pressure -Q 2>/dev/null | tail -1 | tr -d '\n')" >> "$OUT/host.log"
    if [ "$s" -gt "$LIMIT" ]; then
      log "GUARD: swap ${s} MB > ${LIMIT} MB, stopping pid $1 and its children"
      kill -TERM $(tree "$1") 2>/dev/null; sleep 5; kill -KILL $(tree "$1") 2>/dev/null
      return
    fi
    sleep 10
  done
}
{
  echo "started_utc: $(date -u +%FT%TZ)"
  echo "vepyr-e5: $(git -C "$E5" rev-parse HEAD) $(git -C "$E5" branch --show-current) $(git -C "$E5" status --short | wc -l | tr -d ' ') local changes"
  echo "queue vepyr wheel: $(cat "$D"/.venv-lpt/lib/python3.14/site-packages/vepyr-*.dist-info/direct_url.json)"
  echo "vepyr-lpt checkout: $(git -C "$D/tools/vepyr-lpt" rev-parse HEAD) $(git -C "$D/tools/vepyr-lpt" branch --show-current)"
  "$D/.venv-native/bin/python3" -c 'import vepyr; print("native vepyr", vepyr.__version__)'
  echo "colima: $(colima status 2>&1 | tail -1)"
  sysctl hw.memsize vm.swapusage
  df -h "$D" | tail -1
} > "$OUT/environment.txt" 2>&1

log "E5 queue build: start (swap limit ${LIMIT} MB)"
"$D/.venv-lpt/bin/python" "$E5/scripts/thesis/cache-build-parallel/full_queue_e5.py" "$OUT/e5" > "$OUT/e5.out" 2>&1 &
E5PID=$!
guard $E5PID &
wait $E5PID; rc=$?
log "E5 queue build: exit $rc"
[ $rc -eq 0 ] || { log "stop: build failed"; exit 1; }

log "content comparison: start"
"$D/.venv-native/bin/python3" "$PKG/scripts/compare_cache_content.py" "$OUT/e5/current" "$D/cache/116_GRCh38_merged" \
  "$OUT/content.json" > "$OUT/content.out" 2>&1
log "content comparison: exit $?"

log "chr1 annotation on the merged queue cache: start"
"$D/.venv-native/bin/python3" - "$D/runs/m5_completion_20261024/inputs/HG002_normalized_chr1.vcf.gz" "$OUT/chr1.vcf" \
  "$OUT/e5/assembled/116_GRCh38_merged" "$D/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa" > "$OUT/chr1_annotation.json" 2> "$OUT/chr1_annotation.err" <<'EOF'
import hashlib, json, sys, time, vepyr
vcf, out, cache, fasta = sys.argv[1:5]
t0 = time.perf_counter()
vepyr.annotate(vcf=vcf, cache_dir=cache, everything=True, hgvs=True, reference_fasta=fasta, workers=1,
               output_vcf=out, show_progress=False)
dt = time.perf_counter() - t0
h = hashlib.sha256(); n = 0
with open(out, "rb") as f:
    for line in f:
        if not line.startswith(b"#"):
            h.update(line); n += 1
print(json.dumps({"annotation_seconds": dt, "body_sha256": h.hexdigest(), "records": n, "cache": cache}))
EOF
log "chr1 annotation: exit $? $(cat "$OUT/chr1_annotation.json")"
rm -f "$OUT/chr1.vcf"
log "done"
