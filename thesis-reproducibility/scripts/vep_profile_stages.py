#!/usr/bin/env python3
"""Stages of VEP from its NYTProf subroutine profile (vep_profile.sh), thesis item D12.

Inclusive time of the top-level subroutines of the VEP 116 run loop
(Runner::_buffer_to_output per input buffer), as shares of Runner::run:

  transcripts       AnnotationType::Transcript::annotate_InputBuffer  (transcripts from the cache,
                    TranscriptVariation objects for every overlapping transcript)
  known_variants    AnnotationSource::Cache::VariationTabix::annotate_InputBuffer  (tabix lookups)
  csq_entries       OutputFactory::get_all_VariationFeatureOverlapAllele_output_hashes  (the CSQ
                    entries: consequence predicates, HGVS and the other fields are computed lazily here)
    of which        BaseVariationFeatureOverlapAllele::get_all_OverlapConsequences  (SO predicates),
                    TranscriptVariationAllele::hgvs_transcript + hgvs_protein  (HGVS); these two are
                    inclusive times and may overlap, so no remainder is derived from them
  vcf_formatting    OutputFactory::VCF::output_hash_to_vcf_info_chunk
  rest              Runner::run minus the four stages above (input parsing, buffers, other sources)

Exclusive time per sub is not used: most of it lands in shared helpers (interval trees, coordinate
mappers) called from several stages. NYTProf adds a per-call overhead, so the profiled run is about
three times slower than an unprofiled one and the shares are approximate.

Usage: vep_profile_stages.py SUBS_TSV OUT_JSON
"""
import csv, json, sys
from pathlib import Path

TOP = "Bio::EnsEMBL::VEP::Runner::run"
STAGES = {
    "transcripts": ["Bio::EnsEMBL::VEP::AnnotationType::Transcript::annotate_InputBuffer"],
    "known_variants": ["Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix::annotate_InputBuffer"],
    "csq_entries": ["Bio::EnsEMBL::VEP::OutputFactory::get_all_VariationFeatureOverlapAllele_output_hashes"],
    "vcf_formatting": ["Bio::EnsEMBL::VEP::OutputFactory::VCF::output_hash_to_vcf_info_chunk"],
}
WITHIN_CSQ = {
    "consequence_predicates": ["Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele::get_all_OverlapConsequences"],
    "hgvs": ["Bio::EnsEMBL::Variation::TranscriptVariationAllele::hgvs_transcript",
             "Bio::EnsEMBL::Variation::TranscriptVariationAllele::hgvs_protein"],
}


def main():
    subs = {r["sub"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
    incl = lambda names: sum(float(subs[n]["incl_s"]) for n in names)
    total = float(subs[TOP]["incl_s"])
    res = {"total_s": total, "stages": {}, "within_csq_entries": {},
           "csq_entries_calls": int(subs["Bio::EnsEMBL::VEP::OutputFactory::VariationFeatureOverlapAllele_to_output_hash"]["calls"])}
    for k, names in STAGES.items():
        res["stages"][k] = {"seconds": incl(names), "share": incl(names) / total, "subs": names}
    res["stages"]["rest"] = {"seconds": total - sum(v["seconds"] for v in res["stages"].values())}
    res["stages"]["rest"]["share"] = res["stages"]["rest"]["seconds"] / total
    for k, names in WITHIN_CSQ.items():
        res["within_csq_entries"][k] = {"seconds": incl(names), "share_of_total": incl(names) / total, "subs": names}
    Path(sys.argv[2]).write_text(json.dumps(res, indent=1) + "\n")
    print(f"Runner::run {total:.1f} s; CSQ entries built: {res['csq_entries_calls']}")
    for k, v in res["stages"].items():
        print(f"  {k:16s} {v['seconds']:8.1f} s  {100 * v['share']:5.1f}%")
    for k, v in res["within_csq_entries"].items():
        print(f"    of csq_entries: {k:22s} {v['seconds']:7.1f} s  {100 * v['share_of_total']:5.1f}% of total")


if __name__ == "__main__":
    main()
