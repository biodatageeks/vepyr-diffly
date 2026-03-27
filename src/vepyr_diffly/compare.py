from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl
from diffly import compare_frames

from .models import TierResult


@dataclass(frozen=True)
class ComparisonArtifacts:
    diff_frame: pl.DataFrame
    mismatches: pl.DataFrame
    tier: TierResult


def _write_tsv(frame: pl.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(frame.write_csv(separator="\t"), encoding="utf-8")


def compare_tier(
    *,
    name: str,
    left: pl.DataFrame,
    right: pl.DataFrame,
    primary_key: list[str],
    left_name: str,
    right_name: str,
    diff_frame_path: Path,
    mismatches_tsv_path: Path,
) -> ComparisonArtifacts:
    comparison = compare_frames(left, right, primary_key=primary_key)
    left_only = comparison.left_only()
    right_only = comparison.right_only()
    unequal = comparison.joined_unequal()
    joined_equal = comparison.joined_equal()
    summary = comparison.summary(
        top_k_column_changes=10,
        show_sample_primary_key_per_change=True,
        left_name=left_name,
        right_name=right_name,
    ).format(pretty=False)

    diff_frame = pl.concat(
        [
            left_only.with_columns(pl.lit("left_only").alias("diff_kind")),
            right_only.with_columns(pl.lit("right_only").alias("diff_kind")),
            unequal.with_columns(pl.lit("unequal").alias("diff_kind")),
        ],
        how="diagonal_relaxed",
    )
    diff_frame.write_parquet(diff_frame_path)
    _write_tsv(diff_frame, mismatches_tsv_path)

    tier = TierResult(
        name=name,
        summary=summary,
        equal=comparison.equal(),
        left_only_rows=left_only.height,
        right_only_rows=right_only.height,
        unequal_rows=unequal.height,
        joined_equal_rows=joined_equal.height,
        diff_frame_path=diff_frame_path,
        mismatches_tsv_path=mismatches_tsv_path,
    )
    return ComparisonArtifacts(diff_frame=diff_frame, mismatches=diff_frame, tier=tier)

