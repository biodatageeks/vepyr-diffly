#!/usr/bin/env python3
"""Where VEP 116 spends its time (thesis item D12, 2026-10-30): Devel::NYTProf profile of
one VEP process (no --fork) on the HG002 chr22 baseline, with the options of the timing
series (--cache --merged --offline --everything --hgvs --fasta --vcf --no_stats).

Image vep-nytprof:116 = ensemblorg/ensembl-vep:release_116.0 + the distribution package
libdevel-nytprof-perl (same system perl). Only the subroutine profiler runs (stmts=0), which
keeps the profile small and the slowdown lower. The profile is reduced inside the container with
Devel::NYTProf::Data to one row per subroutine (calls, exclusive and inclusive seconds);
the raw profile is deleted. Reported: profiled wall time (NYTProf slows Perl down, so only
shares are used), exclusive time per package, the 40 most expensive subroutines, and the
inclusive time of the stage entry points listed in STAGES (those present in the profile).

Usage: vep_nytprof.py OUT_DIR   (Colima running)
"""
from __future__ import annotations

import gzip
import json
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
BASE = "ensemblorg/ensembl-vep:release_116.0"
TAG = "vep-nytprof:116"
CHR22 = D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz"
FASTA = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
DUMP = r'''
use strict; use warnings; use Devel::NYTProf::Data;
my $p = Devel::NYTProf::Data->new({ filename => $ARGV[0], quiet => 1 });
my $subs = $p->subname_subinfo_map;
for my $name (sort keys %$subs) {
    my $si = $subs->{$name};
    next unless $si->calls;
    printf "%s\t%d\t%.6f\t%.6f\n", $name, $si->calls, $si->excl_time, $si->incl_time;
}
'''
# Entry points of VEP stages (inclusive time); names from ensembl-vep / ensembl-variation 116.
STAGES = [
    "Bio::EnsEMBL::VEP::Runner::run",
    "Bio::EnsEMBL::VEP::Runner::_buffer_to_output",
    "Bio::EnsEMBL::VEP::InputBuffer::next",
    "Bio::EnsEMBL::VEP::Parser::VCF::create_VariationFeatures",
    "Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript::annotate_InputBuffer",
    "Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript::annotate_InputBuffer",
    "Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation::annotate_InputBuffer",
    "Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat::annotate_InputBuffer",
    "Bio::EnsEMBL::VEP::AnnotationSource::Cache::get_features_by_regions_uncached",
    "Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript::deserialize_from_file",
    "Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript::read_from_file",
    "Bio::EnsEMBL::VEP::OutputFactory::get_all_lines_by_InputBuffer",
    "Bio::EnsEMBL::VEP::OutputFactory::VCF::get_all_lines_by_InputBuffer",
    "Bio::EnsEMBL::VEP::OutputFactory::get_all_output_hashes_by_InputBuffer",
    "Bio::EnsEMBL::Variation::TranscriptVariationAllele::hgvs_transcript",
    "Bio::EnsEMBL::Variation::TranscriptVariationAllele::hgvs_protein",
    "Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele::get_all_OverlapConsequences",
    "Bio::EnsEMBL::Variation::TranscriptVariationAllele::get_all_OverlapConsequences",
    "Bio::EnsEMBL::VEP::Haplo::Runner::run",
    "Storable::retrieve",
    "Storable::nstore",
    "Bio::EnsEMBL::Slice::subseq",
    "Bio::DB::HTS::Faidx::get_sequence_no_length",
]


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    img = out / "image"
    img.mkdir(exist_ok=True)
    (img / "Dockerfile").write_text(f"FROM {BASE}\nUSER root\nRUN apt-get update -qq && apt-get install -y -qq "
                                    "libdevel-nytprof-perl > /dev/null && rm -rf /var/lib/apt/lists/*\nUSER vep\n")
    subprocess.run(["docker", "build", "-q", "-t", TAG, str(img)], check=True, capture_output=True)
    (out / "dump.pl").write_text(DUMP)
    vep_args = ["--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
                "--input_file", f"/input/{CHR22.name}", "--output_file", "/output/vep.vcf", "--vcf",
                "--force_overwrite", "--no_stats", "--everything", "--hgvs", "--fasta", f"/ref/{FASTA}"]
    script = ('perl -d:NYTProf "$(command -v vep)" "$@" && '
              'perl /output/dump.pl /output/nytprof.out | gzip -c > /output/subs.tsv.gz && rm -f /output/nytprof.out')
    cmd = ["/usr/bin/time", "-p", "-o", str(out / "time.txt"), "docker", "run", "--rm", "--user", "501:20",
           "--env", "HOME=/tmp", "--env", "NYTPROF=file=/output/nytprof.out:addpid=0:stmts=0",
           "-v", f"{D}/homo_sapiens_merged/116_GRCh38:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
           "-v", f"{D}/input:/ref:ro", "-v", f"{CHR22.parent}:/input:ro", "-v", f"{out}:/output",
           TAG, "sh", "-c", script, "vep", *vep_args]
    (out / "command.txt").write_text(" ".join(cmd) + "\n")
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    (out / "run.log").write_text(p.stdout + p.stderr)
    if p.returncode:
        print(p.stderr[-2000:])
        return 1
    (out / "vep.vcf").unlink(missing_ok=True)

    rows = []
    with gzip.open(out / "subs.tsv.gz", "rt") as f:
        for line in f:
            name, calls, excl, incl = line.rstrip("\n").split("\t")
            rows.append((name, int(calls), float(excl), float(incl)))
    total = sum(r[2] for r in rows)
    pkg = defaultdict(float)
    for name, _, excl, _ in rows:
        pkg[name.rsplit("::", 1)[0] if "::" in name else "main"] += excl
    by_name = {r[0]: r for r in rows}
    res = {"image_id": subprocess.run(["docker", "image", "inspect", "--format", "{{.Id}}", TAG],
                                      capture_output=True, text=True).stdout.strip(),
           "profiled_wall_s": time.time() - t0, "sum_exclusive_s": total, "subroutines": len(rows),
           "packages_by_exclusive": [{"package": k, "seconds": v, "share": v / total}
                                     for k, v in sorted(pkg.items(), key=lambda x: -x[1])[:30]],
           "top_subs_by_exclusive": [{"sub": n, "calls": c, "excl_s": e, "incl_s": i, "excl_share": e / total}
                                     for n, c, e, i in sorted(rows, key=lambda r: -r[2])[:40]],
           "stages_inclusive": {s: {"calls": by_name[s][1], "incl_s": by_name[s][3], "incl_share": by_name[s][3] / total}
                                for s in STAGES if s in by_name}}
    (out / "summary.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: res[k] for k in ("profiled_wall_s", "sum_exclusive_s", "subroutines")}))
    for s, v in res["stages_inclusive"].items():
        print(f"{v['incl_share']:6.1%}  {v['incl_s']:9.1f}s  {s}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
