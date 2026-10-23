#!/usr/bin/env bash
# Repeats 2..5 of the full-genome Ensembl VEP 116 baseline (no --fork) on HG002,
# with the same Docker command as repeat 1 (2026-09-24, runs/chapter5-m5-20260924T061203Z,
# benchmarks/vep/full_forknone.*). Runs sequentially; each run takes ~8 h.
#
# Per repeat it records: exact command, wall time (/usr/bin/time -p, outside Docker),
# exit code, VEP stdout/stderr, peak RSS of the VEP process (VmHWM sampled every 30 s
# via `docker exec`, so the measured command itself is unchanged), record count and a
# SHA-256 of the VCF body (header lines differ between runs). The ~27 GB output VCF is
# deleted after hashing to keep disk usage at one output at a time.
#
# Usage: nohup caffeinate -is tools/run_vep_full_repeats.sh > /dev/null 2>&1 &
set -uo pipefail

D=/Users/lukaszjezapkowicz/workspace/data_vepyr
IMAGE=ensemblorg/ensembl-vep:release_116.0
FIRST_REPEAT=${FIRST_REPEAT:-2}
LAST_REPEAT=${LAST_REPEAT:-5}
REPEAT1_TIME=$D/runs/chapter5-m5-20260924T061203Z/benchmarks/vep/full_forknone.time.txt
OUT=${OUT:-$D/runs/vep_full_repeats_m5_$(date +%Y%m%dT%H%M%S)}
mkdir -p "$OUT"
LOG=$OUT/series.log
log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$LOG"; }

{
  echo "started_utc: $(date -u +%FT%TZ)"
  sw_vers
  sysctl -n machdep.cpu.brand_string hw.memsize
  grep -E '^(cpu|memory|disk|vmType|mountType):' ~/.colima/default/colima.yaml
  echo "image_id: $(docker image inspect --format '{{.Id}}' "$IMAGE")"
  df -h "$D" | tail -1
} > "$OUT/environment.txt" 2>&1

# Reads VmHWM/VmRSS (kB) of the VEP perl process inside the container.
SAMPLER='for d in /proc/[0-9]*; do c=$(cat $d/comm 2>/dev/null); case "$c" in perl|vep) awk -v c="$c" "/^VmHWM|^VmRSS/{printf \"%s%s \", \$1, \$2} END{print c}" $d/status;; esac; done'

for i in $(seq "$FIRST_REPEAT" "$LAST_REPEAT"); do
  RUN=$OUT/repeat_$i
  NAME=vep_full_repeat_$i
  mkdir -p "$RUN"
  docker rm -f "$NAME" > /dev/null 2>&1
  FREE_GB=$(df -g "$D" | awk 'NR==2 {print $4}')
  if [ "$FREE_GB" -lt 35 ]; then
    log "repeat $i: only ${FREE_GB} GB free (need ~30 GB for the output VCF), stopping"
    break
  fi

  CMD=(docker run --rm --name "$NAME" --user 501:20 --env HOME=/tmp
       -v "$D/homo_sapiens_merged/116_GRCh38:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro"
       -v "$D/input:/input:ro"
       -v "$RUN:/output"
       "$IMAGE" vep --dir /opt/vep/.vep --cache --merged --offline --assembly GRCh38
       --input_file /input/HG002_normalized.vcf.gz
       --output_file /output/HG002_vep_forknone.vcf
       --vcf --force_overwrite --no_stats --everything --hgvs
       --fasta /input/Homo_sapiens.GRCh38.dna.primary_assembly.fa)
  printf '%q ' /usr/bin/time -p -o "$RUN/time.txt" "${CMD[@]}" > "$RUN/command.txt"

  log "repeat $i: start"
  date -u +%FT%TZ > "$RUN/started_utc.txt"
  /usr/bin/time -p -o "$RUN/time.txt" "${CMD[@]}" > "$RUN/stdout.txt" 2> "$RUN/stderr.txt" &
  RUNPID=$!

  (
    sleep 20
    while kill -0 "$RUNPID" 2> /dev/null; do
      s=$(docker exec "$NAME" sh -c "$SAMPLER" 2> /dev/null)
      [ -n "$s" ] && echo "$(date -u +%FT%TZ) $s" >> "$RUN/rss_samples.txt"
      sleep 30
    done
  ) &
  SAMPLERPID=$!

  wait "$RUNPID"; RC=$?
  wait "$SAMPLERPID" 2> /dev/null
  echo "$RC" > "$RUN/exit_code.txt"
  date -u +%FT%TZ > "$RUN/finished_utc.txt"
  log "repeat $i: exit $RC, $(grep '^real' "$RUN/time.txt" 2> /dev/null)"

  VCF=$RUN/HG002_vep_forknone.vcf
  if [ -f "$VCF" ]; then
    stat -f %z "$VCF" > "$RUN/output_bytes.txt"
    grep -v '^#' "$VCF" | tee >(wc -l | tr -d ' ' > "$RUN/output_records.txt") \
      | shasum -a 256 | cut -d' ' -f1 > "$RUN/output_body.sha256"
    rm -f "$VCF"
    log "repeat $i: body sha256 $(cat "$RUN/output_body.sha256"), output removed"
  fi
done

python3 -I - "$OUT" "$REPEAT1_TIME" << 'EOF' | tee -a "$LOG"
import json, re, statistics, sys
from pathlib import Path

out, repeat1 = Path(sys.argv[1]), Path(sys.argv[2])

def real(path):
    m = re.search(r"^real\s+([0-9.]+)", path.read_text(), re.M)
    return float(m.group(1)) if m else None

runs = [{"repeat": 1, "wall_s": real(repeat1), "exit_code": 0, "source": str(repeat1)}]
for d in sorted(out.glob("repeat_*"), key=lambda p: int(p.name.split("_")[1])):
    r = {"repeat": int(d.name.split("_")[1]), "source": str(d)}
    r["wall_s"] = real(d / "time.txt") if (d / "time.txt").exists() else None
    r["exit_code"] = int((d / "exit_code.txt").read_text()) if (d / "exit_code.txt").exists() else None
    hwm = [int(x) for x in re.findall(r"VmHWM:(\d+)", (d / "rss_samples.txt").read_text())] \
        if (d / "rss_samples.txt").exists() else []
    r["peak_rss_gb"] = round(max(hwm) * 1024 / 1e9, 3) if hwm else None
    for name in ("output_body.sha256", "output_records.txt", "output_bytes.txt"):
        f = d / name
        r[name.split(".")[0]] = f.read_text().strip() if f.exists() else None
    runs.append(r)

ok = [r["wall_s"] for r in runs if r["exit_code"] == 0 and r["wall_s"]]
summary = {
    "runs": runs,
    "n_successful": len(ok),
    "mean_wall_s": round(statistics.mean(ok), 2) if ok else None,
    "min_wall_s": min(ok) if ok else None,
    "max_wall_s": max(ok) if ok else None,
    "body_sha256_identical_repeats_2_5": len({r.get("output_body") for r in runs[1:]}) == 1,
    "peak_rss_gb_max": max((r["peak_rss_gb"] for r in runs[1:] if r.get("peak_rss_gb")), default=None),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps({k: v for k, v in summary.items() if k != "runs"}, indent=2))
EOF
log "done"
