from pathlib import Path

from vepyr_diffly.normalize import (
    materialize_consequence_buckets,
    materialize_variant_summary,
    normalize_alt_for_csq_allele,
    normalize_annotated_vcf,
)


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


def test_materialize_consequence_buckets_writes_bucket_parts(tmp_path: Path) -> None:
    fixture = Path(__file__).parent / "fixtures" / "annotated_left.vcf"
    result = normalize_annotated_vcf(fixture)
    bucket_root = tmp_path / "buckets"

    buckets = materialize_consequence_buckets(
        vcf_path=fixture,
        csq_fields=result.csq_fields,
        bucket_root=bucket_root,
        side_label="VEP",
        bucket_count=8,
    )

    assert buckets
    written_parts = sorted(bucket_root.glob("bucket-*/*.parquet"))
    assert written_parts
    assert sorted(bucket_root.glob("bucket-*/bucket.meta.json"))


def test_materialize_variant_summary_writes_bucketed_variant_artifacts(tmp_path: Path) -> None:
    fixture = Path(__file__).parent / "fixtures" / "annotated_left.vcf"
    variant_path = tmp_path / "left.variant.parquet"
    bucket_root = tmp_path / "variant-buckets"

    materialize_variant_summary(
        vcf_path=fixture,
        variant_path=variant_path,
        bucket_root=bucket_root,
        side_label="VEP",
        bucket_count=8,
        chunk_variants=1,
    )

    assert variant_path.exists()
    assert sorted(bucket_root.glob("bucket-*/bucket.parquet"))
    assert sorted(bucket_root.glob("bucket-*/bucket.meta.json"))


def test_normalize_alt_for_csq_allele_handles_insertions_and_deletions() -> None:
    assert normalize_alt_for_csq_allele("G", "GGTTT") == "GTTT"
    assert normalize_alt_for_csq_allele("G", "GTTTT") == "TTTT"
    assert normalize_alt_for_csq_allele("CAAAACAAAAACA", "CAAAACA") == "AAAACA"
    assert normalize_alt_for_csq_allele("CAAAACAAAAACA", "C") == "-"


def test_normalize_annotated_vcf_matches_csq_to_the_correct_multiallelic_alt() -> None:
    fixture = Path(__file__).parent / "fixtures" / "annotated_multiallelic.vcf"
    result = normalize_annotated_vcf(fixture)

    variant_counts = {
        (row["pos"], row["alt"]): row["consequence_count"] for row in result.variant.to_dicts()
    }
    assert variant_counts[(100, "GGTTT")] == 1
    assert variant_counts[(100, "GTTTT")] == 1
    assert variant_counts[(200, "CAAAACA")] == 1
    assert variant_counts[(200, "C")] == 1

    consequence_rows = {
        (row["pos"], row["alt"], row["Allele"]) for row in result.consequence.to_dicts()
    }
    assert (100, "GGTTT", "GTTT") in consequence_rows
    assert (100, "GTTTT", "TTTT") in consequence_rows
    assert (200, "CAAAACA", "AAAACA") in consequence_rows
    assert (200, "C", "-") in consequence_rows
    assert (100, "GGTTT", "TTTT") not in consequence_rows
    assert (200, "C", "AAAACA") not in consequence_rows
