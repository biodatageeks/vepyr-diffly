#!/usr/bin/env bash
# B15, second part (2026-11-01): content check of the five queue entities of an E5 run whose
# translation tasks were stopped.
#
# In the run with 2 processes (b15_queue_content.sh) all 100 tasks of the five entities of the
# measured queue (variation, transcript, exon, regulatory, motif; 20 groups) finished; the guard then
# stopped the build when the first two translation tasks (chr1 and chr2 groups, submitted last) ran
# at the same time and swap grew by 13 GB within 30 s. This script:
# 1. checks that the 100 five-entity tasks have status done and moves the stopped translation tasks
#    out of the run directory (kept in killed_translation_tasks/);
# 2. compare_cache_content.py on the five entities against the single-process VEPYR 0.7.0 cache;
# 3. assembles the five entities with the E5 assemble() (hard links, merged chrom_manifest.json) and
#    adds translation_core/translation_sift of the reference cache (hard links), because the queue
#    did not build them;
# 4. annotates the HG002 chr1 slice (baseline, workers=1) with this cache; the body SHA-256 must equal
#    the native chr1 baseline body.
# Usage: b15_finish.sh OUT_DIR VEPYR_E5_CHECKOUT   (OUT_DIR of b15_queue_content.sh)
set -uo pipefail
D=/Users/lukaszjezapkowicz/workspace/data_vepyr
OUT=$1; E5=$2
PKG=$(cd "$(dirname "$0")/.." && pwd)
REF=$D/cache/116_GRCh38_merged
RUN=$OUT/e5/current
log() { echo "$(date -u +%FT%TZ) $* $(pmset -g batt | sed -n 1p)" | tee -a "$OUT/run.log"; }
caffeinate -ims -w $$ &

"$D/.venv-native/bin/python3" - "$RUN" "$OUT/e5/killed_translation_tasks" > "$OUT/tasks.json" <<'EOF' || { log "stop: task check failed"; exit 1; }
import collections, json, shutil, sys
from pathlib import Path
run, killed = Path(sys.argv[1]), Path(sys.argv[2])
rows = [json.loads(p.read_text()) for p in sorted(run.glob("task_*/state.json"))]
done = collections.Counter(r["entity"] for r in rows if r["status"] == "done")
stopped = [r for r in rows if r["status"] != "done"]
assert all(done[e] == 20 for e in ("variation", "transcript", "exon", "regulatory", "motif")), done
assert all(r["entity"] == "translation" for r in stopped), stopped
killed.mkdir(exist_ok=True)
for r in stopped:
    shutil.move(str(run / f"task_{r['task']:03d}"), str(killed / f"task_{r['task']:03d}"))
print(json.dumps({"done": done, "stopped": [{k: r[k] for k in ("task", "group", "entity", "contigs")} for r in stopped],
                  "task_seconds_sum": sum(r.get("seconds", 0) for r in rows if r["status"] == "done")}, indent=1))
EOF
log "tasks: $(tr -d '\n ' < "$OUT/tasks.json" | cut -c1-200)"

log "content comparison (five entities): start"
"$D/.venv-native/bin/python3" "$PKG/scripts/compare_cache_content.py" "$RUN" "$REF" "$OUT/content.json" > "$OUT/content.out" 2>&1
log "content comparison: exit $?"

A=$OUT/e5/assembled/116_GRCh38_merged
rm -rf "$OUT/e5/assembled"; mkdir -p "$A"
"$D/.venv-lpt/bin/python" - "$E5/scripts/thesis/cache-build-parallel" "$RUN" "$A" > "$OUT/assemble.json" <<'EOF'
import json, sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
import full_queue_e5 as q
print(json.dumps(q.assemble(Path(sys.argv[2]), Path(sys.argv[3])), indent=1))
EOF
log "assemble five entities: exit $?"
for e in translation_core translation_sift; do
  mkdir -p "$A/$e" && for f in "$REF/$e"/*; do ln "$f" "$A/$e/"; done
done
log "translation entities linked from the reference cache: $(ls "$A/translation_core" | wc -l | tr -d ' ') + $(ls "$A/translation_sift" | wc -l | tr -d ' ') files"

log "chr1 annotation on the queue cache: start"
"$D/.venv-native/bin/python3" - "$D/runs/m5_completion_20261024/inputs/HG002_normalized_chr1.vcf.gz" "$OUT/chr1.vcf" \
  "$A" "$D/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa" > "$OUT/chr1_annotation.json" 2> "$OUT/chr1_annotation.err" <<'EOF'
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
