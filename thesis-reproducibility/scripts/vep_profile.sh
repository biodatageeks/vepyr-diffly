#!/usr/bin/env bash
# Stage profile of VEP 116 on chr22, baseline scenario (thesis item D12, 2026-10-31).
#
# build: image vep-nytprof:116 = ensemblorg/ensembl-vep:release_116.0 + Ubuntu's libdevel-nytprof-perl.
# run:   the command of the interleaved series (one process, no --fork, merged cache, --everything --hgvs
#        --fasta) under perl -d:NYTProf with subroutine profiling only (stmts=0), then a per-subroutine
#        dump (calls, inclusive and exclusive seconds) read with Devel::NYTProf::Data, and the output body
#        SHA-256, which must equal the unprofiled output (19efe3c7...).
# Usage: vep_profile.sh build OUT_DIR | vep_profile.sh run OUT_DIR
set -euo pipefail
D=/Users/lukaszjezapkowicz/workspace/data_vepyr
MODE=$1; OUT=$2
mkdir -p "$OUT"
if [ "$MODE" = build ]; then
  printf 'FROM ensemblorg/ensembl-vep:release_116.0\nUSER root\nRUN apt-get update -qq && apt-get install -y -qq libdevel-nytprof-perl > /dev/null && rm -rf /var/lib/apt/lists/*\nUSER vep\n' > "$OUT/Dockerfile"
  docker build -t vep-nytprof:116 "$OUT" > "$OUT/build.log" 2>&1
  docker image inspect --format '{{.Id}}' vep-nytprof:116 > "$OUT/image_id.txt"
  docker run --rm --entrypoint perl vep-nytprof:116 -MDevel::NYTProf::Data -e 'print "$Devel::NYTProf::VERSION\n"' >> "$OUT/image_id.txt"
  exit 0
fi
INP=$D/runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz
DUMP='use Devel::NYTProf::Data; my $p = Devel::NYTProf::Data->new({filename => $ARGV[0], quiet => 1});
my $s = $p->subname_subinfo_map; print "sub\tcalls\tincl_s\texcl_s\n";
for my $n (sort keys %$s) { my $i = $s->{$n}; next unless $i->calls; printf "%s\t%d\t%.6f\t%.6f\n", $n, $i->calls, $i->incl_time, $i->excl_time; }'
date -u +%FT%TZ > "$OUT/started_utc.txt"
/usr/bin/time -p -o "$OUT/time.txt" docker run --rm --user 501:20 --env HOME=/tmp \
  --env NYTPROF=file=/output/nytprof.out:stmts=0:addpid=0 \
  -v "$D/homo_sapiens_merged/116_GRCh38:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro" \
  -v "$D/input:/ref:ro" -v "$(dirname "$INP"):/input:ro" -v "$OUT:/output" --entrypoint perl vep-nytprof:116 \
  -d:NYTProf /opt/vep/src/ensembl-vep/vep \
  --dir /opt/vep/.vep --cache --merged --offline --assembly GRCh38 \
  --input_file "/input/$(basename "$INP")" --output_file /output/vep.vcf --vcf --force_overwrite \
  --no_stats --everything --hgvs --fasta /ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa > "$OUT/run.log" 2>&1
grep -v '^#' "$OUT/vep.vcf" | shasum -a 256 | cut -d' ' -f1 > "$OUT/body.sha256"
rm -f "$OUT/vep.vcf"
docker run --rm --user 501:20 -v "$OUT:/output" --entrypoint perl vep-nytprof:116 -e "$DUMP" /output/nytprof.out > "$OUT/subs.tsv"
