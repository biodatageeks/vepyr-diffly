from pathlib import Path
import json

from vepyr_diffly.chromosomes import parse_chromosome_selection
from vepyr_diffly.compare import (
    _metadata_proves_mismatch,
    compare_bucketed_consequence_tier,
    compare_bucketed_variant_tier,
    compare_tier,
)
from vepyr_diffly.normalize import (
    VARIANT_KEY,
    materialize_consequence_buckets,
    materialize_variant_summary,
    normalize_annotated_vcf,
)


def test_compare_tier_detects_variant_differences(tmp_path: Path) -> None:
    fixture_dir = Path(__file__).parent / "fixtures"
    left = normalize_annotated_vcf(fixture_dir / "annotated_left.vcf")
    right = normalize_annotated_vcf(fixture_dir / "annotated_right.vcf")

    artifact = compare_tier(
        name="variant",
        left=left.variant,
        right=right.variant,
        primary_key=VARIANT_KEY,
        left_name="VEP",
        right_name="vepyr",
        diff_frame_path=tmp_path / "variant.parquet",
        mismatches_tsv_path=tmp_path / "variant.tsv",
    )

    assert artifact.tier.equal is False
    assert artifact.tier.right_only_rows == 1
    assert artifact.tier.unequal_rows == 0
    assert artifact.tier.diff_frame_path.exists()
    assert artifact.tier.mismatches_tsv_path.exists()


def test_compare_tier_decodes_percent_escaped_string_values(tmp_path: Path) -> None:
    import polars as pl

    left = pl.DataFrame(
        {
            "chrom": ["1"],
            "pos": [100],
            "ref": ["A"],
            "alt": ["T"],
            "ClinVar_CLNDN": ["HP:0001252%3B_HP:0001776"],
        }
    )
    right = pl.DataFrame(
        {
            "chrom": ["1"],
            "pos": [100],
            "ref": ["A"],
            "alt": ["T"],
            "ClinVar_CLNDN": ["HP:0001252;_HP:0001776"],
        }
    )

    artifact = compare_tier(
        name="consequence",
        left=left,
        right=right,
        primary_key=["chrom", "pos", "ref", "alt"],
        left_name="VEP",
        right_name="vepyr",
        diff_frame_path=tmp_path / "consequence.parquet",
        mismatches_tsv_path=tmp_path / "consequence.tsv",
    )

    assert artifact.tier.equal is True
    assert artifact.tier.joined_equal_rows == 1
    assert artifact.tier.unequal_rows == 0


def test_compare_bucketed_consequence_tier_detects_consequence_differences(
    tmp_path: Path,
) -> None:
    fixture_dir = Path(__file__).parent / "fixtures"
    left_bucket_dir = tmp_path / "left-buckets"
    right_bucket_dir = tmp_path / "right-buckets"

    left = normalize_annotated_vcf(fixture_dir / "annotated_left.vcf")
    right = normalize_annotated_vcf(fixture_dir / "annotated_right.vcf")
    materialize_consequence_buckets(
        vcf_path=fixture_dir / "annotated_left.vcf",
        csq_fields=left.csq_fields,
        bucket_root=left_bucket_dir,
        side_label="VEP",
        bucket_count=8,
    )
    materialize_consequence_buckets(
        vcf_path=fixture_dir / "annotated_right.vcf",
        csq_fields=right.csq_fields,
        bucket_root=right_bucket_dir,
        side_label="vepyr",
        bucket_count=8,
    )

    artifact = compare_bucketed_consequence_tier(
        left_bucket_dir=left_bucket_dir,
        right_bucket_dir=right_bucket_dir,
        csq_fields=left.csq_fields,
        left_name="VEP",
        right_name="vepyr",
        diff_frame_path=tmp_path / "consequence.parquet",
        mismatches_tsv_path=tmp_path / "consequence.tsv",
        bucket_count=8,
        max_workers=2,
    )

    assert artifact.tier.equal is False
    assert artifact.tier.right_only_rows == 2
    assert artifact.tier.diff_frame_path.exists()
    assert artifact.tier.mismatches_tsv_path.exists()


def test_compare_bucketed_variant_tier_detects_variant_differences(tmp_path: Path) -> None:
    fixture_dir = Path(__file__).parent / "fixtures"
    left_bucket_dir = tmp_path / "left-variant-buckets"
    right_bucket_dir = tmp_path / "right-variant-buckets"
    left_variant_path = tmp_path / "left.variant.parquet"
    right_variant_path = tmp_path / "right.variant.parquet"

    materialize_variant_summary(
        vcf_path=fixture_dir / "annotated_left.vcf",
        variant_path=left_variant_path,
        bucket_root=left_bucket_dir,
        side_label="VEP",
        bucket_count=8,
        chunk_variants=1,
    )
    materialize_variant_summary(
        vcf_path=fixture_dir / "annotated_right.vcf",
        variant_path=right_variant_path,
        bucket_root=right_bucket_dir,
        side_label="vepyr",
        bucket_count=8,
        chunk_variants=1,
    )

    artifact = compare_bucketed_variant_tier(
        left_bucket_dir=left_bucket_dir,
        right_bucket_dir=right_bucket_dir,
        left_name="VEP",
        right_name="vepyr",
        diff_frame_path=tmp_path / "variant.parquet",
        mismatches_tsv_path=tmp_path / "variant.tsv",
        bucket_count=8,
        compare_mode="fast",
    )

    assert artifact.tier.equal is False
    assert artifact.tier.right_only_rows == 1
    assert artifact.tier.diff_frame_path.exists()
    assert artifact.tier.mismatches_tsv_path.exists()


def test_compare_bucketed_consequence_tier_fast_mode_skips_exact_diff_for_equal_buckets(
    tmp_path: Path,
) -> None:
    fixture_dir = Path(__file__).parent / "fixtures"
    left_bucket_dir = tmp_path / "left-buckets"
    right_bucket_dir = tmp_path / "right-buckets"

    normalized = normalize_annotated_vcf(fixture_dir / "annotated_left.vcf")
    materialize_consequence_buckets(
        vcf_path=fixture_dir / "annotated_left.vcf",
        csq_fields=normalized.csq_fields,
        bucket_root=left_bucket_dir,
        side_label="VEP",
        bucket_count=8,
        chunk_variants=2,
    )
    materialize_consequence_buckets(
        vcf_path=fixture_dir / "annotated_left.vcf",
        csq_fields=normalized.csq_fields,
        bucket_root=right_bucket_dir,
        side_label="vepyr",
        bucket_count=8,
        chunk_variants=2,
    )

    artifact = compare_bucketed_consequence_tier(
        left_bucket_dir=left_bucket_dir,
        right_bucket_dir=right_bucket_dir,
        csq_fields=normalized.csq_fields,
        left_name="VEP",
        right_name="vepyr",
        diff_frame_path=tmp_path / "consequence.parquet",
        mismatches_tsv_path=tmp_path / "consequence.tsv",
        bucket_count=8,
        max_workers=2,
        compare_mode="fast",
    )

    assert artifact.tier.equal is True
    assert (
        artifact.tier.details["precheck_equal_buckets"] == artifact.tier.details["buckets_compared"]
    )
    assert artifact.tier.details["exact_diff_run_buckets"] == 0
    assert artifact.tier.diff_frame_path.exists()
    assert artifact.tier.mismatches_tsv_path.exists()


def test_compare_bucketed_consequence_tier_can_ignore_duplicate_count(tmp_path: Path) -> None:
    left_bucket_dir = tmp_path / "left-buckets"
    right_bucket_dir = tmp_path / "right-buckets"
    (left_bucket_dir / "bucket-0000").mkdir(parents=True)
    (right_bucket_dir / "bucket-0000").mkdir(parents=True)

    left_frame = {
        "chrom": ["1"],
        "pos": [100],
        "ref": ["A"],
        "alt": ["T"],
        "ClinVar": ["123"],
        "ClinVar_CLNSIG": ["Benign"],
        "duplicate_count": [5],
    }
    right_frame = {
        "chrom": ["1"],
        "pos": [100],
        "ref": ["A"],
        "alt": ["T"],
        "ClinVar": ["123"],
        "ClinVar_CLNSIG": ["Benign"],
        "duplicate_count": [1],
    }

    import polars as pl

    pl.DataFrame(left_frame).write_parquet(left_bucket_dir / "bucket-0000" / "bucket.parquet")
    pl.DataFrame(right_frame).write_parquet(right_bucket_dir / "bucket-0000" / "bucket.parquet")

    artifact = compare_bucketed_consequence_tier(
        left_bucket_dir=left_bucket_dir,
        right_bucket_dir=right_bucket_dir,
        csq_fields=["ClinVar", "ClinVar_CLNSIG"],
        left_name="VEP",
        right_name="vepyr",
        diff_frame_path=tmp_path / "consequence.parquet",
        mismatches_tsv_path=tmp_path / "consequence.tsv",
        bucket_count=1,
        compare_mode="fast",
        compare_duplicate_count=False,
    )

    assert artifact.tier.equal is True
    assert artifact.tier.unequal_rows == 0
    assert artifact.tier.left_only_rows == 0
    assert artifact.tier.right_only_rows == 0


def test_bucket_metadata_can_shortcut_obvious_mismatch(tmp_path: Path) -> None:
    fixture_dir = Path(__file__).parent / "fixtures"
    left_bucket_dir = tmp_path / "left-buckets"
    right_bucket_dir = tmp_path / "right-buckets"

    left = normalize_annotated_vcf(fixture_dir / "annotated_left.vcf")
    materialize_consequence_buckets(
        vcf_path=fixture_dir / "annotated_left.vcf",
        csq_fields=left.csq_fields,
        bucket_root=left_bucket_dir,
        side_label="VEP",
        bucket_count=8,
    )
    materialize_consequence_buckets(
        vcf_path=fixture_dir / "annotated_left.vcf",
        csq_fields=left.csq_fields,
        bucket_root=right_bucket_dir,
        side_label="vepyr",
        bucket_count=8,
    )

    bucket_path = next(left_bucket_dir.glob("bucket-*/bucket.parquet"))
    peer_path = right_bucket_dir / bucket_path.parent.name / "bucket.parquet"
    peer_meta = peer_path.with_suffix(".meta.json")
    payload = json.loads(peer_meta.read_text(encoding="utf-8"))
    payload["row_count"] += 1
    peer_meta.write_text(json.dumps(payload) + "\n", encoding="utf-8")

    assert _metadata_proves_mismatch([bucket_path], [peer_path]) is True


def test_compare_bucketed_variant_tier_tracks_per_chromosome_results(tmp_path: Path) -> None:
    fixture_dir = Path(__file__).parent / "fixtures"
    left_bucket_dir = tmp_path / "left-variant-buckets"
    right_bucket_dir = tmp_path / "right-variant-buckets"
    left_variant_path = tmp_path / "left.variant.parquet"
    right_variant_path = tmp_path / "right.variant.parquet"
    _, chr1_aliases = parse_chromosome_selection("1")

    materialize_variant_summary(
        vcf_path=fixture_dir / "annotated_multi_chrom_left.vcf",
        variant_path=left_variant_path,
        bucket_root=left_bucket_dir,
        side_label="VEP",
        bucket_count=8,
        chunk_variants=1,
        chromosome_aliases=chr1_aliases,
    )
    materialize_variant_summary(
        vcf_path=fixture_dir / "annotated_multi_chrom_right.vcf",
        variant_path=right_variant_path,
        bucket_root=right_bucket_dir,
        side_label="vepyr",
        bucket_count=8,
        chunk_variants=1,
        chromosome_aliases=chr1_aliases,
    )

    artifact = compare_bucketed_variant_tier(
        left_bucket_dir=left_bucket_dir,
        right_bucket_dir=right_bucket_dir,
        left_name="VEP",
        right_name="vepyr",
        diff_frame_path=tmp_path / "variant.parquet",
        mismatches_tsv_path=tmp_path / "variant.tsv",
        bucket_count=8,
        compare_mode="fast",
    )

    assert artifact.tier.equal is True
    assert artifact.tier.per_chromosome["1"]["joined_equal_rows"] == 1
