from pathlib import Path

from vepyr_diffly.compare import compare_bucketed_consequence_tier, compare_tier
from vepyr_diffly.normalize import (
    VARIANT_KEY,
    materialize_consequence_buckets,
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
