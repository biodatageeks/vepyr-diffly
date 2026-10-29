#!/usr/bin/env bash
# Stability of the VEP 116 five-plugin reference (thesis item B17, 2026-10-30).
#
# 1. Rebuilds the chr22 five-plugin reference a second time with the same builder
#    (build_vep_plugin_reference.sh), image, pinned plugin code (corrected CADD.pm) and plugin
#    sources served by rclone, as on 2026-10-25, into a separate work directory (the reference
#    in output/116/plugins is not touched). The SHA-256 of the record body is compared with the
#    reference (results/references-2026-10-25/chr22_reference_body.sha256).
# 2. Annotates the same input with the same VEP command plus --fork 10 and compares again.
# Both runs are timed (one run each: not a timing series).
#
# Usage: vep_plugin_stability.sh OUT_DIR   (Colima running; OUT_DIR under data_vepyr, which is
#        the only directory mounted into the container)
set -euo pipefail
D=/Users/lukaszjezapkowicz/workspace/data_vepyr
OUT=$(cd "$1" 2>/dev/null && pwd || { mkdir -p "$1" && cd "$1" && pwd; })
case "$OUT" in "$D"/*) ;; *) echo "OUT_DIR must be under $D" >&2; exit 1 ;; esac
PKG=$(cd "$(dirname "$0")/.." && pwd)
BUILDER=$D/tools/vepyr/e2e-testing/scripts/build_vep_plugin_reference.sh
PLUGINS=$D/output/116/plugins/plugin_code
REF=$D/output/116/plugins/HG002_chr22_5plugins_vep116_caddfix.vcf.gz
IMAGE=ensemblorg/ensembl-vep:release_116.0
WORK=$OUT/work
SOURCE_URL=http://127.0.0.1:18080
mkdir -p "$WORK/output"

body() { bgzip -dc "$1" | grep -v '^#' | shasum -a 256 | cut -d' ' -f1; }

{
  echo "started_utc: $(date -u +%FT%TZ)"
  sw_vers -productVersion
  grep -E '^(cpu|memory|vmType|mountType):' ~/.colima/default/colima.yaml
  echo "image_id: $(docker image inspect --format '{{.Id}} {{.Architecture}}' $IMAGE)"
  echo "builder: $(git -C "$D/tools/vepyr" log -1 --format=%h -- e2e-testing/scripts/build_vep_plugin_reference.sh)"
  shasum -a 256 "$PLUGINS"/*.pm
} > "$OUT/environment.txt"

rclone serve http "$D/plugin_input" --addr 127.0.0.1:18080 --read-only > "$OUT/rclone.log" 2>&1 &
RCLONE=$!
trap 'kill $RCLONE 2> /dev/null' EXIT
for _ in $(seq 30); do curl -sf --head "$SOURCE_URL/clinvar/clinvar.vcf.gz" > /dev/null && break; sleep 1; done

# --- 1. second build of the reference -------------------------------------------------------
RUN2=$WORK/output/HG002_chr22_5plugins_run2.vcf
/usr/bin/time -p -o "$OUT/run2.time.txt" env DATA_VEPYR_DIR="$D" VEP_PLUGIN_DIR="$PLUGINS" \
  VEP_PLUGIN_SOURCE_URL="$SOURCE_URL" VEP_NORMALIZED_VCF="$D/input/HG002_normalized.vcf.gz" \
  VEP_OUTPUT_VCF="$RUN2" VEP_KEEP_PLAIN=0 "$BUILDER" 22 "$WORK" > "$OUT/run2.log" 2>&1
kill $RCLONE 2> /dev/null

# --- 2. same command with --fork 10 (slices and input from the work directory) ---------------
DBNSFP_COLS=$(sed -n '/^DBNSFP_COLS="/,/"$/p' "$BUILDER" | tr -d '\\\n' | sed 's/^DBNSFP_COLS="//; s/"$//')
S=/data/${WORK#"$D"/}/slices
FORK=$WORK/output/HG002_chr22_5plugins_fork10.vcf
/usr/bin/time -p -o "$OUT/fork10.time.txt" docker run --rm --user "$(id -u):$(id -g)" \
  -v "$D":/data -v "$PLUGINS":/plugins:ro $IMAGE \
  vep --cache --cache_version 116 --dir_cache /data --offline --merged \
      --everything --no_stats --force_overwrite --vcf \
      --fasta /data/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
      --input_file "/data/${WORK#"$D"/}/input/HG002_norm_chr22.vcf.gz" \
      --output_file "/data/${FORK#"$D"/}" \
      --dir_plugins /plugins --fork 10 \
      --custom "$S/clinvar_chr22.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNVC,CLNVI" \
      --plugin "SpliceAI,snv=$S/spliceai_chr22.vcf.gz,indel=$S/spliceai_chr22.vcf.gz" \
      --plugin "CADD,snv=$S/cadd_snv_chr22.tsv.gz,indels=$S/cadd_indel_chr22.tsv.gz" \
      --plugin "AlphaMissense,file=$S/alphamissense_chr22.tsv.gz" \
      --plugin "dbNSFP,$S/dbNSFP5.3.1a_grch38_chr22.gz,$DBNSFP_COLS" > "$OUT/fork10.log" 2>&1
bgzip -f "$FORK"

# --- 3. checksums ----------------------------------------------------------------------------
python3 -I - "$OUT" "$(cat "$PKG/results/references-2026-10-25/chr22_reference_body.sha256")" \
  "$(body "$REF")" "$(body "$RUN2.gz")" "$(body "$FORK.gz")" \
  "$(bgzip -dc "$RUN2.gz" | grep -vc '^#')" "$(bgzip -dc "$FORK.gz" | grep -vc '^#')" << 'EOF'
import json, sys
from pathlib import Path
out, recorded, ref, run2, fork, n2, nf = sys.argv[1:]
t = lambda n: float(next(l.split()[1] for l in (Path(out) / n).read_text().splitlines() if l.startswith("real")))
res = {"reference_body_sha256_recorded": recorded, "reference_body_sha256_now": ref,
       "run2_body_sha256": run2, "fork10_body_sha256": fork, "run2_records": int(n2), "fork10_records": int(nf),
       "run2_identical": run2 == recorded, "fork10_identical": fork == recorded,
       "run2_seconds_incl_slicing": t("run2.time.txt"), "fork10_seconds": t("fork10.time.txt")}
(Path(out) / "summary.json").write_text(json.dumps(res, indent=1) + "\n")
print(json.dumps(res, indent=1))
EOF
