from pathlib import Path

from vepyr_diffly.normalize import normalize_annotated_vcf


def test_normalize_annotated_vcf_builds_variant_and_consequence_tables() -> None:
    fixture = Path(__file__).parent / "fixtures" / "annotated_left.vcf"
    result = normalize_annotated_vcf(fixture)

    assert result.variant.height == 2
    assert result.consequence.height == 2
    assert result.csq_fields[:4] == ["Allele", "Consequence", "IMPACT", "SYMBOL"]
    assert result.variant.columns == [
        "chrom",
        "pos",
        "ref",
        "alt",
        "record_count",
        "consequence_count",
        "ids",
        "filters",
    ]
    assert result.variant.schema["record_count"].is_integer()
    assert result.variant.schema["consequence_count"].is_integer()
    assert result.consequence.schema["duplicate_count"].is_integer()
