from pathlib import Path

from vepyr_diffly.mismatch_analysis import (
    analyze_consequence_mismatches,
    extract_csq_examples_from_analysis,
)


def test_analyze_consequence_mismatches_groups_pairs_and_unpaired(tmp_path: Path) -> None:
    mismatches_tsv = tmp_path / "consequence_mismatches.tsv"
    mismatches_tsv.write_text(
        "\t".join(
            [
                "chrom",
                "pos",
                "ref",
                "alt",
                "Feature",
                "Consequence",
                "IMPACT",
                "HGNC_ID",
                "HGVSc",
                "HGVSp",
                "diff_kind",
                "duplicate_count_left",
                "duplicate_count_right",
            ]
        )
        + "\n"
        + "\n".join(
            [
                "chr1\t10\tA\tT\tENST1\tmissense_variant\tMODERATE\tHGNC:1\tc.1A>T\tp.Lys1Asn\tleft_only\t1\t",
                "chr1\t10\tA\tT\tENST1\tmissense_variant\tMODERATE\tHGNC:1\tc.1A>T\t.\tright_only\t\t1",
                "chr1\t20\tG\tGA\tENST2\tframeshift_variant\tHIGH\tHGNC:2\tc.2dup\t.\tleft_only\t1\t",
                "chr1\t20\tG\tGA\tENST2\tsplice_region_variant\tLOW\tHGNC:2\tc.2dup\t.\tright_only\t\t1",
                "chr2\t30\tC\tG\tENST3\tdownstream_gene_variant\tMODIFIER\t.\t.\t.\tleft_only\t1\t",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    payload = analyze_consequence_mismatches(mismatches_tsv)

    assert payload["total_rows"] == 5
    assert payload["paired_keys"] == 2
    assert payload["unpaired_rows"] == 1
    assert payload["left_only_rows"] == 3
    assert payload["right_only_rows"] == 2
    assert {"field": "HGVSp", "paired_keys": 1} in payload["top_differing_fields"]
    assert {"field": "Consequence", "paired_keys": 1} in payload["top_differing_fields"]
    assert {"field": "IMPACT", "paired_keys": 1} in payload["top_differing_fields"]
    assert {"category": "missing_hgvs", "paired_keys": 1} in payload["top_semantic_categories"]
    assert {
        "category": "consequence_and_impact_reclassification",
        "paired_keys": 1,
    } in payload["top_semantic_categories"]
    assert payload["top_field_signatures"][0]["paired_keys"] == 1
    assert len(payload["paired_examples"]) == 2
    assert "missing_hgvs" in payload["category_examples"]
    assert "consequence_and_impact_reclassification" in payload["category_examples"]
    assert len(payload["unpaired_examples"]) == 1


def test_extract_csq_examples_from_analysis(tmp_path: Path) -> None:
    analysis = {
        "category_examples": {
            "missing_hgvs": [
                {
                    "key": {
                        "chrom": "chr1",
                        "pos": "10",
                        "ref": "A",
                        "alt": "T",
                        "Feature": "ENST1",
                    },
                    "diff_fields": ["HGVSc"],
                    "categories": ["missing_hgvs"],
                }
            ]
        }
    }
    left_vcf = tmp_path / "left.vcf"
    right_vcf = tmp_path / "right.vcf"
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    left_vcf.write_text(
        header
        + "chr1\t10\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|GENE|ENSG1|Transcript|ENST1|protein_coding|||ENST1:c.1A>T|ENSP1:p.Lys1Asn\n",
        encoding="utf-8",
    )
    right_vcf.write_text(
        header
        + "chr1\t10\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|GENE|ENSG1|Transcript|ENST1|protein_coding||||ENSP1:p.Lys1Asn\n",
        encoding="utf-8",
    )

    payload = extract_csq_examples_from_analysis(
        analysis,
        left_vcf=left_vcf,
        right_vcf=right_vcf,
        per_category_limit=3,
    )

    example = payload["categories"]["missing_hgvs"][0]
    assert example["left_raw_csq"] == [
        "T|missense_variant|MODERATE|GENE|ENSG1|Transcript|ENST1|protein_coding|||ENST1:c.1A>T|ENSP1:p.Lys1Asn"
    ]
    assert example["right_raw_csq"] == [
        "T|missense_variant|MODERATE|GENE|ENSG1|Transcript|ENST1|protein_coding||||ENSP1:p.Lys1Asn"
    ]
