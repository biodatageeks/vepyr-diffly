from pathlib import Path
import polars as pl

from vepyr_diffly.chromosomes import parse_chromosome_selection
from vepyr_diffly.normalize import (
    _extract_csq_entries_from_info_bytes,
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


def test_materialize_consequence_buckets_respects_selected_csq_field_indexes(tmp_path: Path) -> None:
    fixture = Path(__file__).parent / "fixtures" / "annotated_left.vcf"
    result = normalize_annotated_vcf(fixture)
    bucket_root = tmp_path / "buckets"
    selected_fields = ["IMPACT", "SYMBOL"]
    field_indexes = {field: index for index, field in enumerate(result.csq_fields)}

    materialize_consequence_buckets(
        vcf_path=fixture,
        csq_fields=selected_fields,
        field_indexes=field_indexes,
        bucket_root=bucket_root,
        side_label="VEP",
        bucket_count=8,
    )

    bucket_paths = sorted(bucket_root.glob("bucket-*/bucket.parquet"))
    frame = pl.read_parquet([str(path) for path in bucket_paths]).sort(["pos", "alt"])
    assert frame.select("IMPACT").to_series().to_list() == ["MODERATE", "LOW"]
    assert frame.select("SYMBOL").to_series().to_list() == ["GENE1", "GENE2"]


def test_materialize_consequence_buckets_can_drop_empty_plugin_rows(tmp_path: Path) -> None:
    fixture = tmp_path / "plugin_only.vcf"
    fixture.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from test. Format: Allele|Consequence|am_pathogenicity|am_class">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "1\t100\t.\tG\tA\t.\t.\tCSQ=A|downstream_gene_variant|.|.,A|missense_variant|0.3296|likely_benign",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    bucket_root = tmp_path / "buckets"

    materialize_consequence_buckets(
        vcf_path=fixture,
        csq_fields=["am_pathogenicity", "am_class"],
        field_indexes={"am_pathogenicity": 2, "am_class": 3},
        drop_empty_csq_rows=True,
        bucket_root=bucket_root,
        side_label="VEP",
        bucket_count=8,
    )

    frame = pl.read_parquet([str(path) for path in sorted(bucket_root.glob("bucket-*/bucket.parquet"))])
    assert frame.select("am_pathogenicity").to_series().to_list() == ["0.3296"]
    assert frame.select("am_class").to_series().to_list() == ["likely_benign"]


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


def test_extract_csq_entries_from_info_bytes_keeps_unescaped_semicolons_inside_csq() -> None:
    info = (
        b"AC=1;CSQ=T|foo|HP:0001252;_HP:0001776|single_nucleotide_variant|,"
        b"T|bar|HP:1;_HP:2|snv|;SOMETHING=1"
    )

    entries = _extract_csq_entries_from_info_bytes(info)

    assert entries == [
        "T|foo|HP:0001252;_HP:0001776|single_nucleotide_variant|",
        "T|bar|HP:1;_HP:2|snv|",
    ]


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


def test_normalize_annotated_vcf_can_filter_standard_chromosome_aliases() -> None:
    fixture = Path(__file__).parent / "fixtures" / "annotated_multi_chrom_left.vcf"
    _, aliases = parse_chromosome_selection("1")

    result = normalize_annotated_vcf(fixture, chromosome_aliases=aliases)

    assert result.variant.height == 1
    assert result.variant.item(0, "chrom") == "1"
    assert result.consequence.height == 1
