from __future__ import annotations

from dataclasses import asdict, dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
from decimal import Decimal, InvalidOperation
import json
import os
from pathlib import Path
import re
import shutil
import threading
from typing import Callable
from urllib.parse import unquote

import polars as pl
from diffly import compare_frames

from .models import TierResult
from .normalize import CONSEQUENCE_BUCKET_COUNT, VARIANT_KEY
from .progress import ProgressReporter


MISMATCH_TSV_ROW_LIMIT = 2000
_NUMERIC_STRING_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


def _normalize_compare_string(value: str | None) -> str | None:
    if value is None:
        return None
    decoded = unquote(value)
    if not _NUMERIC_STRING_RE.match(decoded):
        return decoded
    try:
        decimal = Decimal(decoded)
    except InvalidOperation:
        return decoded
    if decimal.is_zero():
        return "0"
    return format(decimal.normalize(), "f")


@dataclass(frozen=True)
class ComparisonArtifacts:
    diff_frame_path: Path
    mismatches_tsv_path: Path
    tier: TierResult


@dataclass(frozen=True)
class BucketCompareResult:
    bucket_id: int
    left_rows: int
    right_rows: int
    joined_equal_rows: int
    left_only_rows: int
    right_only_rows: int
    unequal_rows: int
    diff_frame_path: Path
    mismatches_tsv_path: Path
    exact_diff_run: bool = True
    precheck_equal: bool = False
    per_chromosome: dict[str, dict[str, int]] | None = None


def _merge_per_chromosome(
    aggregate: dict[str, dict[str, int]],
    incoming: dict[str, dict[str, int]] | None,
) -> None:
    if not incoming:
        return
    for chrom, counts in incoming.items():
        bucket = aggregate.setdefault(
            chrom,
            {
                "left_rows": 0,
                "right_rows": 0,
                "joined_equal_rows": 0,
                "left_only_rows": 0,
                "right_only_rows": 0,
                "unequal_rows": 0,
            },
        )
        for key, value in counts.items():
            bucket[key] = bucket.get(key, 0) + int(value)


def _group_counts_by_chrom(frame: pl.DataFrame) -> dict[str, int]:
    if frame.is_empty():
        return {}
    grouped = frame.group_by("chrom").agg(pl.len().alias("rows"))
    return {str(row["chrom"]): int(row["rows"]) for row in grouped.to_dicts()}


def _build_per_chromosome_counts(
    *,
    left_frame: pl.DataFrame,
    right_frame: pl.DataFrame,
    key: list[str],
) -> dict[str, dict[str, int]]:
    value_columns = [column for column in left_frame.columns if column not in key]
    chroms = set(left_frame.get_column("chrom").to_list()) | set(right_frame.get_column("chrom").to_list())
    if not chroms:
        return {}

    left_only = left_frame.join(right_frame.select(key), on=key, how="anti")
    right_only = right_frame.join(left_frame.select(key), on=key, how="anti")
    joined = left_frame.join(right_frame, on=key, how="inner", suffix="_right")
    if value_columns:
        equal_expr = pl.all_horizontal(
            *(pl.col(column) == pl.col(f"{column}_right") for column in value_columns)
        )
        joined_equal = joined.filter(equal_expr)
        joined_unequal = joined.filter(~equal_expr)
    else:
        joined_equal = joined
        joined_unequal = joined.head(0)

    left_counts = _group_counts_by_chrom(left_frame)
    right_counts = _group_counts_by_chrom(right_frame)
    joined_equal_counts = _group_counts_by_chrom(joined_equal)
    left_only_counts = _group_counts_by_chrom(left_only)
    right_only_counts = _group_counts_by_chrom(right_only)
    unequal_counts = _group_counts_by_chrom(joined_unequal)
    output: dict[str, dict[str, int]] = {}
    for chrom in sorted(str(value) for value in chroms):
        output[chrom] = {
            "left_rows": left_counts.get(chrom, 0),
            "right_rows": right_counts.get(chrom, 0),
            "joined_equal_rows": joined_equal_counts.get(chrom, 0),
            "left_only_rows": left_only_counts.get(chrom, 0),
            "right_only_rows": right_only_counts.get(chrom, 0),
            "unequal_rows": unequal_counts.get(chrom, 0),
        }
    return output


def _normalize_compare_eager_frame(frame: pl.DataFrame, key: list[str]) -> pl.DataFrame:
    preserved_columns = {"chrom", "pos", "ref", "alt"}
    string_columns = [
        column
        for column, dtype in frame.schema.items()
        if column not in preserved_columns and dtype == pl.String
    ]
    if not string_columns:
        return frame
    return frame.with_columns(
        pl.col(column)
        .cast(pl.String)
        .map_elements(_normalize_compare_string, return_dtype=pl.String)
        .alias(column)
        for column in string_columns
    )


def _sink_lazy_frame(frame: pl.LazyFrame, path: Path, *, as_csv: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if as_csv and hasattr(frame, "sink_csv"):
        frame.sink_csv(str(path), separator="\t")
        return
    if not as_csv and hasattr(frame, "sink_parquet"):
        frame.sink_parquet(str(path))
        return
    eager = frame.collect()
    if as_csv:
        path.write_text(eager.write_csv(separator="\t"), encoding="utf-8")
    else:
        eager.write_parquet(path)


def _sink_sampled_tsv(frame: pl.LazyFrame, path: Path, *, row_limit: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    sampled = frame.limit(row_limit).collect()
    path.write_text(sampled.write_csv(separator="\t"), encoding="utf-8")


def _write_empty_diff_artifacts(
    *,
    diff_frame_path: Path,
    mismatches_tsv_path: Path,
    schema: dict[str, pl.DataType],
) -> None:
    empty = pl.DataFrame(schema={**schema, "diff_kind": pl.String})
    diff_frame_path.parent.mkdir(parents=True, exist_ok=True)
    mismatches_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    empty.write_parquet(diff_frame_path)
    mismatches_tsv_path.write_text(empty.write_csv(separator="\t"), encoding="utf-8")


def compare_tier(
    *,
    name: str,
    left: pl.DataFrame | pl.LazyFrame,
    right: pl.DataFrame | pl.LazyFrame,
    primary_key: list[str],
    left_name: str,
    right_name: str,
    diff_frame_path: Path,
    mismatches_tsv_path: Path,
    reporter: ProgressReporter | None = None,
) -> ComparisonArtifacts:
    left_eager = left.collect() if isinstance(left, pl.LazyFrame) else left
    right_eager = right.collect() if isinstance(right, pl.LazyFrame) else right
    left_eager = _normalize_compare_eager_frame(left_eager, primary_key)
    right_eager = _normalize_compare_eager_frame(right_eager, primary_key)
    if reporter is not None:
        reporter.stage(f"{name}: initializing diffly comparison")
    comparison = compare_frames(left_eager.lazy(), right_eager.lazy(), primary_key=primary_key)

    if reporter is not None:
        reporter.stage(f"{name}: counting left rows")
    left_rows = comparison.num_rows_left()
    if reporter is not None:
        reporter.stage(f"{name}: counting right rows")
    right_rows = comparison.num_rows_right()
    if reporter is not None:
        reporter.stage(f"{name}: counting joined equal rows")
    joined_equal_rows = comparison.num_rows_joined_equal()
    if reporter is not None:
        reporter.stage(f"{name}: counting left-only rows")
    left_only_rows = comparison.num_rows_left_only()
    if reporter is not None:
        reporter.stage(f"{name}: counting right-only rows")
    right_only_rows = comparison.num_rows_right_only()
    if reporter is not None:
        reporter.stage(f"{name}: counting unequal rows")
    unequal_rows = comparison.num_rows_joined_unequal()
    equal = left_only_rows == 0 and right_only_rows == 0 and unequal_rows == 0

    if reporter is not None:
        reporter.log(
            f"{name}: counts left={left_rows} right={right_rows} "
            f"joined_equal={joined_equal_rows} left_only={left_only_rows} "
            f"right_only={right_only_rows} unequal={unequal_rows}"
        )
        reporter.stage(f"{name}: generating human-readable summary")
    summary = comparison.summary(
        top_k_column_changes=10,
        show_sample_primary_key_per_change=True,
        left_name=left_name,
        right_name=right_name,
    ).format(pretty=False)

    diff_frame = pl.concat(
        [
            comparison.left_only(lazy=True).with_columns(pl.lit("left_only").alias("diff_kind")),
            comparison.right_only(lazy=True).with_columns(pl.lit("right_only").alias("diff_kind")),
            comparison.joined_unequal(lazy=True).with_columns(pl.lit("unequal").alias("diff_kind")),
        ],
        how="diagonal_relaxed",
    )
    if reporter is not None:
        reporter.stage(
            f"{name}: writing diff parquet → {diff_frame_path}",
            tracked_paths=[diff_frame_path],
        )
    _sink_lazy_frame(diff_frame, diff_frame_path, as_csv=False)
    if reporter is not None:
        reporter.stage(
            f"{name}: writing mismatches tsv → {mismatches_tsv_path}",
            tracked_paths=[mismatches_tsv_path],
        )
    _sink_sampled_tsv(diff_frame, mismatches_tsv_path, row_limit=MISMATCH_TSV_ROW_LIMIT)

    tier = TierResult(
        name=name,
        summary=summary,
        equal=equal,
        left_only_rows=left_only_rows,
        right_only_rows=right_only_rows,
        unequal_rows=unequal_rows,
        joined_equal_rows=joined_equal_rows,
        diff_frame_path=diff_frame_path,
        mismatches_tsv_path=mismatches_tsv_path,
        per_chromosome=_build_per_chromosome_counts(
            left_frame=left_eager,
            right_frame=right_eager,
            key=primary_key,
        ),
    )
    return ComparisonArtifacts(
        diff_frame_path=diff_frame_path,
        mismatches_tsv_path=mismatches_tsv_path,
        tier=tier,
    )


def _variant_schema() -> dict[str, pl.DataType]:
    return {
        "chrom": pl.String,
        "pos": pl.Int64,
        "ref": pl.String,
        "alt": pl.String,
        "record_count": pl.Int64,
        "consequence_count": pl.Int64,
        "ids": pl.String,
        "filters": pl.String,
    }


def _scan_single_bucket(paths: list[Path], schema: dict[str, pl.DataType]) -> pl.LazyFrame:
    if not paths:
        return pl.DataFrame(schema=schema).lazy()
    return pl.scan_parquet([str(path) for path in paths])


def _load_bucket_metadata(path: Path) -> dict[str, int] | None:
    meta_path = path.with_suffix(".meta.json")
    if not meta_path.exists():
        return None
    payload = json.loads(meta_path.read_text(encoding="utf-8"))
    return {str(key): int(value) for key, value in payload.items()}


def _metadata_proves_mismatch(left_paths: list[Path], right_paths: list[Path]) -> bool:
    if len(left_paths) != 1 or len(right_paths) != 1:
        return False
    if left_paths[0].name != "bucket.parquet" or right_paths[0].name != "bucket.parquet":
        return False
    left_meta = _load_bucket_metadata(left_paths[0])
    right_meta = _load_bucket_metadata(right_paths[0])
    if left_meta is None or right_meta is None:
        return False
    for payload in (left_meta, right_meta):
        payload.pop("file_size_bytes", None)
    return left_meta != right_meta


def _precheck_bucket_frames(
    *,
    bucket_id: int,
    left_paths: list[Path],
    right_paths: list[Path],
    key: list[str],
    schema: dict[str, pl.DataType],
    exact_counts_on_mismatch: bool = False,
) -> BucketCompareResult:
    if _metadata_proves_mismatch(left_paths, right_paths):
        return BucketCompareResult(
            bucket_id=bucket_id,
            left_rows=0,
            right_rows=0,
            joined_equal_rows=0,
            left_only_rows=0,
            right_only_rows=0,
            unequal_rows=0,
            diff_frame_path=Path(),
            mismatches_tsv_path=Path(),
            exact_diff_run=False,
            precheck_equal=False,
        )
    left_eager = _normalize_compare_eager_frame(
        _scan_single_bucket(left_paths, schema).sort(key).collect(),
        key,
    )
    right_eager = _normalize_compare_eager_frame(
        _scan_single_bucket(right_paths, schema).sort(key).collect(),
        key,
    )
    left_rows = left_eager.height
    right_rows = right_eager.height
    precheck_equal = left_eager.equals(right_eager)
    left_only_rows = 0
    right_only_rows = 0
    unequal_rows = 0
    joined_equal_rows = left_rows if precheck_equal else 0
    if not precheck_equal and exact_counts_on_mismatch:
        comparison = compare_frames(left_eager.lazy(), right_eager.lazy(), primary_key=key)
        left_only_rows = comparison.num_rows_left_only()
        right_only_rows = comparison.num_rows_right_only()
        unequal_rows = comparison.num_rows_joined_unequal()
    return BucketCompareResult(
        bucket_id=bucket_id,
        left_rows=left_rows,
        right_rows=right_rows,
        joined_equal_rows=joined_equal_rows,
        left_only_rows=left_only_rows,
        right_only_rows=right_only_rows,
        unequal_rows=unequal_rows,
        diff_frame_path=Path(),
        mismatches_tsv_path=Path(),
        exact_diff_run=False,
        precheck_equal=precheck_equal,
        per_chromosome=_build_per_chromosome_counts(
            left_frame=left_eager,
            right_frame=right_eager,
            key=key,
        ),
    )


def _compare_direct_bucket(
    *,
    bucket_id: int,
    left_paths: list[Path],
    right_paths: list[Path],
    key: list[str],
    schema: dict[str, pl.DataType],
    diff_frame_path: Path,
    mismatches_tsv_path: Path | None = None,
) -> BucketCompareResult:
    left_eager = _normalize_compare_eager_frame(_scan_single_bucket(left_paths, schema).collect(), key)
    right_eager = _normalize_compare_eager_frame(
        _scan_single_bucket(right_paths, schema).collect(),
        key,
    )
    comparison = compare_frames(
        left_eager.lazy(),
        right_eager.lazy(),
        primary_key=key,
    )
    left_rows = comparison.num_rows_left()
    right_rows = comparison.num_rows_right()
    joined_equal_rows = comparison.num_rows_joined_equal()
    left_only_rows = comparison.num_rows_left_only()
    right_only_rows = comparison.num_rows_right_only()
    unequal_rows = comparison.num_rows_joined_unequal()
    diff_frame = pl.concat(
        [
            comparison.left_only(lazy=True).with_columns(pl.lit("left_only").alias("diff_kind")),
            comparison.right_only(lazy=True).with_columns(pl.lit("right_only").alias("diff_kind")),
            comparison.joined_unequal(lazy=True).with_columns(pl.lit("unequal").alias("diff_kind")),
        ],
        how="diagonal_relaxed",
    )
    _sink_lazy_frame(diff_frame, diff_frame_path, as_csv=False)
    if mismatches_tsv_path is not None:
        _sink_sampled_tsv(diff_frame, mismatches_tsv_path, row_limit=MISMATCH_TSV_ROW_LIMIT)
    return BucketCompareResult(
        bucket_id=bucket_id,
        left_rows=left_rows,
        right_rows=right_rows,
        joined_equal_rows=joined_equal_rows,
        left_only_rows=left_only_rows,
        right_only_rows=right_only_rows,
        unequal_rows=unequal_rows,
        diff_frame_path=diff_frame_path,
        mismatches_tsv_path=mismatches_tsv_path or Path(),
        exact_diff_run=True,
        precheck_equal=(left_only_rows == 0 and right_only_rows == 0 and unequal_rows == 0),
        per_chromosome=_build_per_chromosome_counts(
            left_frame=left_eager,
            right_frame=right_eager,
            key=key,
        ),
    )


def _consequence_schema(
    csq_fields: list[str], *, include_duplicate_count: bool = True
) -> dict[str, pl.DataType]:
    schema: dict[str, pl.DataType] = {
        "chrom": pl.String,
        "pos": pl.Int64,
        "ref": pl.String,
        "alt": pl.String,
    }
    for field in csq_fields:
        schema[field] = pl.String
    if include_duplicate_count:
        schema["duplicate_count"] = pl.Int64
    return schema


def _aggregate_bucket_frame(
    paths: list[Path],
    csq_fields: list[str],
    *,
    include_duplicate_count: bool = True,
) -> pl.LazyFrame:
    if not paths:
        return pl.DataFrame(
            schema=_consequence_schema(csq_fields, include_duplicate_count=include_duplicate_count)
        ).lazy()
    if len(paths) == 1 and paths[0].name == "bucket.parquet":
        frame = pl.scan_parquet(str(paths[0]))
        if include_duplicate_count:
            return frame
        return frame.select("chrom", "pos", "ref", "alt", *csq_fields)
    group_key = ["chrom", "pos", "ref", "alt", *csq_fields]
    grouped = pl.scan_parquet([str(path) for path in paths]).group_by(group_key)
    if include_duplicate_count:
        return grouped.agg(pl.col("duplicate_count").sum().cast(pl.Int64).alias("duplicate_count"))
    return grouped.agg()


def _count_lazy_rows(frame: pl.LazyFrame) -> int:
    return frame.select(pl.len()).collect().item()


def _collect_sorted_bucket_frame(
    paths: list[Path],
    csq_fields: list[str],
    *,
    include_duplicate_count: bool = True,
) -> pl.DataFrame:
    key = ["chrom", "pos", "ref", "alt", *csq_fields]
    return _aggregate_bucket_frame(
        paths,
        csq_fields,
        include_duplicate_count=include_duplicate_count,
    ).sort(key).collect()


def _count_bucket_left_only(
    left_frame: pl.LazyFrame,
    right_frame: pl.LazyFrame,
    key: list[str],
) -> int:
    return left_frame.join(right_frame, on=key, how="anti").select(pl.len()).collect().item()


def _count_bucket_unequal(
    left_frame: pl.LazyFrame,
    right_frame: pl.LazyFrame,
    key: list[str],
) -> int:
    return (
        left_frame.join(right_frame, on=key, how="inner", suffix="_right")
        .filter(pl.col("duplicate_count") != pl.col("duplicate_count_right"))
        .select(pl.len())
        .collect()
        .item()
    )


def precheck_single_bucket(
    *,
    bucket_id: int,
    left_paths: list[Path],
    right_paths: list[Path],
    csq_fields: list[str],
    exact_counts_on_mismatch: bool = False,
    compare_duplicate_count: bool = True,
) -> BucketCompareResult:
    key = ["chrom", "pos", "ref", "alt", *csq_fields]
    if _metadata_proves_mismatch(left_paths, right_paths):
        return BucketCompareResult(
            bucket_id=bucket_id,
            left_rows=0,
            right_rows=0,
            joined_equal_rows=0,
            left_only_rows=0,
            right_only_rows=0,
            unequal_rows=0,
            diff_frame_path=Path(),
            mismatches_tsv_path=Path(),
            exact_diff_run=False,
            precheck_equal=False,
        )
    left_eager = _normalize_compare_eager_frame(
        _collect_sorted_bucket_frame(
        left_paths,
        csq_fields,
        include_duplicate_count=compare_duplicate_count,
    ),
        key,
    )
    right_eager = _normalize_compare_eager_frame(
        _collect_sorted_bucket_frame(
        right_paths,
        csq_fields,
        include_duplicate_count=compare_duplicate_count,
    ),
        key,
    )
    left_rows = left_eager.height
    right_rows = right_eager.height
    precheck_equal = left_eager.equals(right_eager)
    left_only_rows = 0
    right_only_rows = 0
    unequal_rows = 0
    joined_equal_rows = left_rows if precheck_equal else 0
    if not precheck_equal and exact_counts_on_mismatch:
        left_frame = left_eager.lazy()
        right_frame = right_eager.lazy()
        left_only_rows = _count_bucket_left_only(left_frame, right_frame, key)
        right_only_rows = _count_bucket_left_only(right_frame, left_frame, key)
        if compare_duplicate_count and left_only_rows == 0 and right_only_rows == 0:
            unequal_rows = _count_bucket_unequal(left_frame, right_frame, key)
    return BucketCompareResult(
        bucket_id=bucket_id,
        left_rows=left_rows,
        right_rows=right_rows,
        joined_equal_rows=joined_equal_rows,
        left_only_rows=left_only_rows,
        right_only_rows=right_only_rows,
        unequal_rows=unequal_rows,
        diff_frame_path=Path(),
        mismatches_tsv_path=Path(),
        exact_diff_run=False,
        precheck_equal=precheck_equal,
        per_chromosome=_build_per_chromosome_counts(
            left_frame=left_eager,
            right_frame=right_eager,
            key=key,
        ),
    )


def compare_single_bucket(
    *,
    bucket_id: int,
    left_paths: list[Path],
    right_paths: list[Path],
    csq_fields: list[str],
    diff_frame_path: Path,
    mismatches_tsv_path: Path | None = None,
    compare_duplicate_count: bool = True,
) -> BucketCompareResult:
    key = ["chrom", "pos", "ref", "alt", *csq_fields]
    left_eager = _normalize_compare_eager_frame(
        _aggregate_bucket_frame(
        left_paths,
        csq_fields,
        include_duplicate_count=compare_duplicate_count,
    ).collect(),
        key,
    )
    right_eager = _normalize_compare_eager_frame(
        _aggregate_bucket_frame(
        right_paths,
        csq_fields,
        include_duplicate_count=compare_duplicate_count,
    ).collect(),
        key,
    )
    comparison = compare_frames(
        left_eager.lazy(),
        right_eager.lazy(),
        primary_key=key,
    )
    left_rows = comparison.num_rows_left()
    right_rows = comparison.num_rows_right()
    joined_equal_rows = comparison.num_rows_joined_equal()
    left_only_rows = comparison.num_rows_left_only()
    right_only_rows = comparison.num_rows_right_only()
    unequal_rows = comparison.num_rows_joined_unequal()

    diff_frame = pl.concat(
        [
            comparison.left_only(lazy=True).with_columns(pl.lit("left_only").alias("diff_kind")),
            comparison.right_only(lazy=True).with_columns(pl.lit("right_only").alias("diff_kind")),
            comparison.joined_unequal(lazy=True).with_columns(pl.lit("unequal").alias("diff_kind")),
        ],
        how="diagonal_relaxed",
    )
    _sink_lazy_frame(diff_frame, diff_frame_path, as_csv=False)
    if mismatches_tsv_path is not None:
        _sink_lazy_frame(diff_frame, mismatches_tsv_path, as_csv=True)
    return BucketCompareResult(
        bucket_id=bucket_id,
        left_rows=left_rows,
        right_rows=right_rows,
        joined_equal_rows=joined_equal_rows,
        left_only_rows=left_only_rows,
        right_only_rows=right_only_rows,
        unequal_rows=unequal_rows,
        diff_frame_path=diff_frame_path,
        mismatches_tsv_path=mismatches_tsv_path or Path(),
        exact_diff_run=True,
        precheck_equal=(left_only_rows == 0 and right_only_rows == 0 and unequal_rows == 0),
        per_chromosome=_build_per_chromosome_counts(
            left_frame=left_eager,
            right_frame=right_eager,
            key=key,
        ),
    )


def _discover_bucket_jobs(
    *,
    left_bucket_dir: Path,
    right_bucket_dir: Path,
    bucket_count: int,
) -> list[tuple[int, list[Path], list[Path]]]:
    jobs: list[tuple[int, list[Path], list[Path]]] = []
    for bucket_id in range(bucket_count):
        left_dir = left_bucket_dir / f"bucket-{bucket_id:04d}"
        right_dir = right_bucket_dir / f"bucket-{bucket_id:04d}"
        left_compacted = left_dir / "bucket.parquet"
        right_compacted = right_dir / "bucket.parquet"
        left_paths = (
            [left_compacted] if left_compacted.exists() else sorted(left_dir.glob("*.parquet"))
        )
        right_paths = (
            [right_compacted] if right_compacted.exists() else sorted(right_dir.glob("*.parquet"))
        )
        if left_paths or right_paths:
            jobs.append((bucket_id, left_paths, right_paths))
    return jobs


def compare_bucketed_variant_tier(
    *,
    left_bucket_dir: Path,
    right_bucket_dir: Path,
    left_name: str,
    right_name: str,
    diff_frame_path: Path,
    mismatches_tsv_path: Path,
    reporter: ProgressReporter | None = None,
    bucket_count: int = CONSEQUENCE_BUCKET_COUNT,
    compare_mode: str = "fast",
    fingerprint_only: bool = False,
) -> ComparisonArtifacts:
    temp_root = diff_frame_path.parent / ".variant_bucket_compare"
    if temp_root.exists():
        shutil.rmtree(temp_root)
    temp_diff_dir = temp_root / "diff"
    temp_diff_dir.mkdir(parents=True, exist_ok=True)

    jobs = _discover_bucket_jobs(
        left_bucket_dir=left_bucket_dir,
        right_bucket_dir=right_bucket_dir,
        bucket_count=bucket_count,
    )
    if reporter is not None:
        reporter.stage(
            f"variant: comparing {len(jobs)} buckets",
            tracked_paths=[left_bucket_dir, right_bucket_dir],
        )

    if not jobs:
        _write_empty_diff_artifacts(
            diff_frame_path=diff_frame_path,
            mismatches_tsv_path=mismatches_tsv_path,
            schema=_variant_schema(),
        )
        tier = TierResult(
            name="variant",
            summary="No variant rows were produced by either side.\n",
            equal=True,
            left_only_rows=0,
            right_only_rows=0,
            unequal_rows=0,
            joined_equal_rows=0,
            diff_frame_path=diff_frame_path,
            mismatches_tsv_path=mismatches_tsv_path,
        )
        return ComparisonArtifacts(
            diff_frame_path=diff_frame_path,
            mismatches_tsv_path=mismatches_tsv_path,
            tier=tier,
        )

    aggregate = {
        "left_rows": 0,
        "right_rows": 0,
        "joined_equal_rows": 0,
        "left_only_rows": 0,
        "right_only_rows": 0,
        "unequal_rows": 0,
    }
    per_chromosome: dict[str, dict[str, int]] = {}
    diff_paths: list[Path] = []
    precheck_equal_buckets = 0
    exact_diff_run_buckets = 0
    key = VARIANT_KEY
    schema = _variant_schema()

    for index, (bucket_id, left_paths, right_paths) in enumerate(jobs, start=1):
        if compare_mode == "fast":
            precheck = _precheck_bucket_frames(
                bucket_id=bucket_id,
                left_paths=left_paths,
                right_paths=right_paths,
                key=key,
                schema=schema,
                exact_counts_on_mismatch=fingerprint_only,
            )
            if precheck.precheck_equal or fingerprint_only:
                result = precheck
            else:
                result = _compare_direct_bucket(
                    bucket_id=bucket_id,
                    left_paths=left_paths,
                    right_paths=right_paths,
                    key=key,
                    schema=schema,
                    diff_frame_path=temp_diff_dir / f"bucket-{bucket_id:04d}.parquet",
                )
        else:
            result = _compare_direct_bucket(
                bucket_id=bucket_id,
                left_paths=left_paths,
                right_paths=right_paths,
                key=key,
                schema=schema,
                diff_frame_path=temp_diff_dir / f"bucket-{bucket_id:04d}.parquet",
            )
        aggregate["left_rows"] += result.left_rows
        aggregate["right_rows"] += result.right_rows
        aggregate["joined_equal_rows"] += result.joined_equal_rows
        aggregate["left_only_rows"] += result.left_only_rows
        aggregate["right_only_rows"] += result.right_only_rows
        aggregate["unequal_rows"] += result.unequal_rows
        _merge_per_chromosome(per_chromosome, result.per_chromosome)
        if result.precheck_equal:
            precheck_equal_buckets += 1
        if result.exact_diff_run and result.diff_frame_path.exists():
            diff_paths.append(result.diff_frame_path)
            exact_diff_run_buckets += 1
        if reporter is not None:
            reporter.log(
                f"variant: completed bucket {bucket_id:04d} ({index}/{len(jobs)}) "
                f"left={result.left_rows} right={result.right_rows} "
                f"unequal={result.unequal_rows} left_only={result.left_only_rows} "
                f"right_only={result.right_only_rows} exact_diff_run={'yes' if result.exact_diff_run else 'no'}"
            )

    if diff_paths and not fingerprint_only:
        _sink_lazy_frame(
            pl.scan_parquet([str(path) for path in diff_paths]), diff_frame_path, as_csv=False
        )
        _sink_sampled_tsv(
            pl.scan_parquet([str(path) for path in diff_paths]),
            mismatches_tsv_path,
            row_limit=MISMATCH_TSV_ROW_LIMIT,
        )
    else:
        _write_empty_diff_artifacts(
            diff_frame_path=diff_frame_path,
            mismatches_tsv_path=mismatches_tsv_path,
            schema=schema,
        )

    equal = (
        aggregate["left_only_rows"] == 0
        and aggregate["right_only_rows"] == 0
        and aggregate["unequal_rows"] == 0
    )
    summary = "\n".join(
        [
            "Bucketized diffly variant comparison",
            f"bucket_count={bucket_count}",
            f"buckets_compared={len(jobs)}",
            f"compare_mode={compare_mode}",
            f"fingerprint_only={fingerprint_only}",
            f"precheck_equal_buckets={precheck_equal_buckets}",
            f"exact_diff_run_buckets={exact_diff_run_buckets}",
            f"left_rows={aggregate['left_rows']}",
            f"right_rows={aggregate['right_rows']}",
            f"joined_equal_rows={aggregate['joined_equal_rows']}",
            f"left_only_rows={aggregate['left_only_rows']}",
            f"right_only_rows={aggregate['right_only_rows']}",
            f"unequal_rows={aggregate['unequal_rows']}",
            "",
        ]
    )
    tier = TierResult(
        name="variant",
        summary=summary,
        equal=equal,
        left_only_rows=aggregate["left_only_rows"],
        right_only_rows=aggregate["right_only_rows"],
        unequal_rows=aggregate["unequal_rows"],
        joined_equal_rows=aggregate["joined_equal_rows"],
        diff_frame_path=diff_frame_path,
        mismatches_tsv_path=mismatches_tsv_path,
        details={
            "compare_mode": compare_mode,
            "fingerprint_only": fingerprint_only,
            "precheck_equal_buckets": precheck_equal_buckets,
            "exact_diff_run_buckets": exact_diff_run_buckets,
            "bucket_count": bucket_count,
            "buckets_compared": len(jobs),
            "workers": 1,
        },
        per_chromosome=per_chromosome,
    )
    shutil.rmtree(temp_root, ignore_errors=True)
    return ComparisonArtifacts(
        diff_frame_path=diff_frame_path,
        mismatches_tsv_path=mismatches_tsv_path,
        tier=tier,
    )


def _shard_bucket_ids(bucket_ids: list[int], worker_count: int) -> list[list[int]]:
    shards: list[list[int]] = [[] for _ in range(worker_count)]
    for index, bucket_id in enumerate(bucket_ids):
        shards[index % worker_count].append(bucket_id)
    return [shard for shard in shards if shard]


def compare_bucketed_consequence_tier(
    *,
    left_bucket_dir: Path,
    right_bucket_dir: Path,
    csq_fields: list[str],
    left_name: str,
    right_name: str,
    diff_frame_path: Path,
    mismatches_tsv_path: Path,
    reporter: ProgressReporter | None = None,
    bucket_count: int = CONSEQUENCE_BUCKET_COUNT,
    max_workers: int | None = None,
    compare_mode: str = "fast",
    fingerprint_only: bool = False,
    compare_duplicate_count: bool = True,
) -> ComparisonArtifacts:
    temp_root = diff_frame_path.parent / ".consequence_bucket_compare"
    if temp_root.exists():
        shutil.rmtree(temp_root)
    temp_diff_dir = temp_root / "diff"
    temp_diff_dir.mkdir(parents=True, exist_ok=True)

    jobs = _discover_bucket_jobs(
        left_bucket_dir=left_bucket_dir,
        right_bucket_dir=right_bucket_dir,
        bucket_count=bucket_count,
    )

    if reporter is not None:
        reporter.stage(
            f"consequence: comparing {len(jobs)} buckets with diffly",
            tracked_paths=[left_bucket_dir, right_bucket_dir],
        )

    if not jobs:
        _write_empty_diff_artifacts(
            diff_frame_path=diff_frame_path,
            mismatches_tsv_path=mismatches_tsv_path,
            schema=_consequence_schema(csq_fields, include_duplicate_count=compare_duplicate_count),
        )
        tier = TierResult(
            name="consequence",
            summary="No consequence rows were produced by either side.\n",
            equal=True,
            left_only_rows=0,
            right_only_rows=0,
            unequal_rows=0,
            joined_equal_rows=0,
            diff_frame_path=diff_frame_path,
            mismatches_tsv_path=mismatches_tsv_path,
        )
        return ComparisonArtifacts(
            diff_frame_path=diff_frame_path,
            mismatches_tsv_path=mismatches_tsv_path,
            tier=tier,
        )

    cpu_count = max(1, os.cpu_count() or 1)
    worker_count = max_workers or min(len(jobs), max(2, cpu_count - 1))
    shards = _shard_bucket_ids([bucket_id for bucket_id, _, _ in jobs], worker_count)
    if reporter is not None:
        reporter.log(f"consequence: launching {len(shards)} compare worker threads")

    aggregate = {
        "left_rows": 0,
        "right_rows": 0,
        "joined_equal_rows": 0,
        "left_only_rows": 0,
        "right_only_rows": 0,
        "unequal_rows": 0,
    }
    per_chromosome: dict[str, dict[str, int]] = {}
    diff_paths: list[Path] = []
    completed = 0
    exact_diff_run_buckets = 0
    precheck_equal_buckets = 0
    progress_lock = threading.Lock()

    def _log_bucket_progress(
        completed_in_shard: int,
        total_in_shard: int,
        result: BucketCompareResult,
    ) -> None:
        nonlocal completed
        if reporter is None:
            return
        with progress_lock:
            completed += 1
            reporter.log(
                f"worker: completed bucket {result.bucket_id:04d} "
                f"({completed_in_shard}/{total_in_shard} in shard)"
            )
            reporter.log(
                f"consequence: completed bucket {result.bucket_id:04d} "
                f"({completed}/{len(jobs)}) left={result.left_rows} right={result.right_rows} "
                f"unequal={result.unequal_rows} left_only={result.left_only_rows} "
                f"right_only={result.right_only_rows} exact_diff_run={'yes' if result.exact_diff_run else 'no'}"
            )

    with ThreadPoolExecutor(max_workers=len(shards)) as executor:
        future_map = {
            executor.submit(
                compare_bucket_shard,
                left_bucket_dir=left_bucket_dir,
                right_bucket_dir=right_bucket_dir,
                bucket_ids=shard,
                csq_fields=csq_fields,
                temp_diff_dir=temp_diff_dir,
                temp_tsv_dir=None,
                compare_mode=compare_mode,
                fingerprint_only=fingerprint_only,
                compare_duplicate_count=compare_duplicate_count,
                on_bucket_complete=_log_bucket_progress,
            ): shard
            for shard in shards
        }
        for future in as_completed(future_map):
            shard_results = future.result()
            for result in shard_results:
                aggregate["left_rows"] += result.left_rows
                aggregate["right_rows"] += result.right_rows
                aggregate["joined_equal_rows"] += result.joined_equal_rows
                aggregate["left_only_rows"] += result.left_only_rows
                aggregate["right_only_rows"] += result.right_only_rows
                aggregate["unequal_rows"] += result.unequal_rows
                _merge_per_chromosome(per_chromosome, result.per_chromosome)
                if result.exact_diff_run and result.diff_frame_path.exists():
                    diff_paths.append(result.diff_frame_path)
                    exact_diff_run_buckets += 1
                if result.precheck_equal:
                    precheck_equal_buckets += 1

    if reporter is not None:
        reporter.stage(
            f"consequence: merging bucket diff artifacts → {diff_frame_path}",
            tracked_paths=[temp_diff_dir],
        )
    if diff_paths and not fingerprint_only:
        _sink_lazy_frame(
            pl.scan_parquet([str(path) for path in diff_paths]), diff_frame_path, as_csv=False
        )
        if reporter is not None:
            reporter.stage(
                f"consequence: writing merged mismatches tsv → {mismatches_tsv_path}",
                tracked_paths=[diff_frame_path],
            )
        _sink_sampled_tsv(
            pl.scan_parquet([str(path) for path in diff_paths]),
            mismatches_tsv_path,
            row_limit=MISMATCH_TSV_ROW_LIMIT,
        )
    else:
        _write_empty_diff_artifacts(
            diff_frame_path=diff_frame_path,
            mismatches_tsv_path=mismatches_tsv_path,
            schema=_consequence_schema(csq_fields, include_duplicate_count=compare_duplicate_count),
        )

    equal = (
        aggregate["left_only_rows"] == 0
        and aggregate["right_only_rows"] == 0
        and aggregate["unequal_rows"] == 0
    )
    summary = "\n".join(
        [
            "Bucketized diffly consequence comparison",
            f"workers={len(shards)}",
            f"bucket_count={bucket_count}",
            f"buckets_compared={len(jobs)}",
            f"compare_mode={compare_mode}",
            f"fingerprint_only={fingerprint_only}",
            f"precheck_equal_buckets={precheck_equal_buckets}",
            f"exact_diff_run_buckets={exact_diff_run_buckets}",
            f"mismatch_tsv_row_limit={MISMATCH_TSV_ROW_LIMIT}",
            f"left_rows={aggregate['left_rows']}",
            f"right_rows={aggregate['right_rows']}",
            f"joined_equal_rows={aggregate['joined_equal_rows']}",
            f"left_only_rows={aggregate['left_only_rows']}",
            f"right_only_rows={aggregate['right_only_rows']}",
            f"unequal_rows={aggregate['unequal_rows']}",
            "",
        ]
    )
    tier = TierResult(
        name="consequence",
        summary=summary,
        equal=equal,
        left_only_rows=aggregate["left_only_rows"],
        right_only_rows=aggregate["right_only_rows"],
        unequal_rows=aggregate["unequal_rows"],
        joined_equal_rows=aggregate["joined_equal_rows"],
        diff_frame_path=diff_frame_path,
        mismatches_tsv_path=mismatches_tsv_path,
        details={
            "compare_mode": compare_mode,
            "fingerprint_only": fingerprint_only,
            "precheck_equal_buckets": precheck_equal_buckets,
            "exact_diff_run_buckets": exact_diff_run_buckets,
            "bucket_count": bucket_count,
            "buckets_compared": len(jobs),
            "workers": len(shards),
        },
        per_chromosome=per_chromosome,
    )
    shutil.rmtree(temp_root, ignore_errors=True)
    return ComparisonArtifacts(
        diff_frame_path=diff_frame_path,
        mismatches_tsv_path=mismatches_tsv_path,
        tier=tier,
    )


def compare_bucket_shard(
    *,
    left_bucket_dir: Path,
    right_bucket_dir: Path,
    bucket_ids: list[int],
    csq_fields: list[str],
    temp_diff_dir: Path,
    temp_tsv_dir: Path | None,
    compare_mode: str = "fast",
    fingerprint_only: bool = False,
    compare_duplicate_count: bool = True,
    on_bucket_complete: Callable[[int, int, BucketCompareResult], None] | None = None,
) -> list[BucketCompareResult]:
    results: list[BucketCompareResult] = []
    total_buckets = len(bucket_ids)
    for index, bucket_id in enumerate(bucket_ids, start=1):
        left_dir = left_bucket_dir / f"bucket-{bucket_id:04d}"
        right_dir = right_bucket_dir / f"bucket-{bucket_id:04d}"
        left_compacted = left_dir / "bucket.parquet"
        right_compacted = right_dir / "bucket.parquet"
        left_paths = (
            [left_compacted] if left_compacted.exists() else sorted(left_dir.glob("*.parquet"))
        )
        right_paths = (
            [right_compacted] if right_compacted.exists() else sorted(right_dir.glob("*.parquet"))
        )
        if not left_paths and not right_paths:
            continue
        if compare_mode == "fast":
            precheck = precheck_single_bucket(
                bucket_id=bucket_id,
                left_paths=left_paths,
                right_paths=right_paths,
                csq_fields=csq_fields,
                exact_counts_on_mismatch=fingerprint_only,
                compare_duplicate_count=compare_duplicate_count,
            )
            if precheck.precheck_equal or fingerprint_only:
                result = precheck
            else:
                result = compare_single_bucket(
                    bucket_id=bucket_id,
                    left_paths=left_paths,
                    right_paths=right_paths,
                    csq_fields=csq_fields,
                    diff_frame_path=temp_diff_dir / f"bucket-{bucket_id:04d}.parquet",
                    mismatches_tsv_path=None
                    if temp_tsv_dir is None
                    else temp_tsv_dir / f"bucket-{bucket_id:04d}.tsv",
                    compare_duplicate_count=compare_duplicate_count,
                )
        else:
            result = compare_single_bucket(
                bucket_id=bucket_id,
                left_paths=left_paths,
                right_paths=right_paths,
                csq_fields=csq_fields,
                diff_frame_path=temp_diff_dir / f"bucket-{bucket_id:04d}.parquet",
                mismatches_tsv_path=None
                if temp_tsv_dir is None
                else temp_tsv_dir / f"bucket-{bucket_id:04d}.tsv",
                compare_duplicate_count=compare_duplicate_count,
            )
        results.append(result)
        if on_bucket_complete is not None:
            on_bucket_complete(index, total_buckets, result)
    return results


def dump_bucket_shard_summary(results: list[BucketCompareResult], summary_path: Path) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps({"results": [asdict(result) for result in results]}, default=str, indent=2)
        + "\n",
        encoding="utf-8",
    )
