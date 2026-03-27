from pathlib import Path

from vepyr_diffly.compare import compare_tier
from vepyr_diffly.normalize import VARIANT_KEY, normalize_annotated_vcf


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

