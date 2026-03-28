from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Callable

import polars as pl
from diffly import compare_frames

from .models import TierResult
from .normalize import CONSEQUENCE_BUCKET_COUNT
from .progress import ProgressReporter


MISMATCH_TSV_ROW_LIMIT = 2000


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
    if reporter is not None:
        reporter.stage(f"{name}: initializing diffly comparison")
    comparison = compare_frames(left, right, primary_key=primary_key)

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
    )
    return ComparisonArtifacts(
        diff_frame_path=diff_frame_path,
        mismatches_tsv_path=mismatches_tsv_path,
        tier=tier,
    )


def _consequence_schema(csq_fields: list[str]) -> dict[str, pl.DataType]:
    schema: dict[str, pl.DataType] = {
        "chrom": pl.String,
        "pos": pl.Int64,
        "ref": pl.String,
        "alt": pl.String,
    }
    for field in csq_fields:
        schema[field] = pl.String
    schema["duplicate_count"] = pl.Int64
    return schema


def _aggregate_bucket_frame(paths: list[Path], csq_fields: list[str]) -> pl.LazyFrame:
    if not paths:
        return pl.DataFrame(schema=_consequence_schema(csq_fields)).lazy()
    group_key = ["chrom", "pos", "ref", "alt", *csq_fields]
    return (
        pl.scan_parquet([str(path) for path in paths])
        .group_by(group_key)
        .agg(pl.col("duplicate_count").sum().cast(pl.Int64).alias("duplicate_count"))
    )


def compare_single_bucket(
    *,
    bucket_id: int,
    left_paths: list[Path],
    right_paths: list[Path],
    csq_fields: list[str],
    diff_frame_path: Path,
    mismatches_tsv_path: Path,
) -> BucketCompareResult:
    key = ["chrom", "pos", "ref", "alt", *csq_fields]
    comparison = compare_frames(
        _aggregate_bucket_frame(left_paths, csq_fields),
        _aggregate_bucket_frame(right_paths, csq_fields),
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
        mismatches_tsv_path=mismatches_tsv_path,
    )


def _discover_bucket_jobs(
    *,
    left_bucket_dir: Path,
    right_bucket_dir: Path,
    bucket_count: int,
) -> list[tuple[int, list[Path], list[Path]]]:
    jobs: list[tuple[int, list[Path], list[Path]]] = []
    for bucket_id in range(bucket_count):
        left_paths = sorted((left_bucket_dir / f"bucket-{bucket_id:04d}").glob("*.parquet"))
        right_paths = sorted((right_bucket_dir / f"bucket-{bucket_id:04d}").glob("*.parquet"))
        if left_paths or right_paths:
            jobs.append((bucket_id, left_paths, right_paths))
    return jobs


def _shard_bucket_ids(bucket_ids: list[int], worker_count: int) -> list[list[int]]:
    shards: list[list[int]] = [[] for _ in range(worker_count)]
    for index, bucket_id in enumerate(bucket_ids):
        shards[index % worker_count].append(bucket_id)
    return [shard for shard in shards if shard]


def _current_python() -> str:
    return sys.executable


def _load_shard_results(summary_path: Path) -> list[BucketCompareResult]:
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    return [
        BucketCompareResult(
            bucket_id=item["bucket_id"],
            left_rows=item["left_rows"],
            right_rows=item["right_rows"],
            joined_equal_rows=item["joined_equal_rows"],
            left_only_rows=item["left_only_rows"],
            right_only_rows=item["right_only_rows"],
            unequal_rows=item["unequal_rows"],
            diff_frame_path=Path(item["diff_frame_path"]),
            mismatches_tsv_path=Path(item["mismatches_tsv_path"]),
        )
        for item in payload["results"]
    ]


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
) -> ComparisonArtifacts:
    temp_root = diff_frame_path.parent / ".consequence_bucket_compare"
    if temp_root.exists():
        shutil.rmtree(temp_root)
    temp_diff_dir = temp_root / "diff"
    temp_tsv_dir = temp_root / "tsv"
    temp_summary_dir = temp_root / "summary"
    temp_diff_dir.mkdir(parents=True, exist_ok=True)
    temp_tsv_dir.mkdir(parents=True, exist_ok=True)
    temp_summary_dir.mkdir(parents=True, exist_ok=True)

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
            schema=_consequence_schema(csq_fields),
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
        reporter.log(
            f"consequence: launching {len(shards)} compare worker processes"
        )

    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1]) + (
        os.pathsep + env["PYTHONPATH"] if "PYTHONPATH" in env else ""
    )
    processes: list[tuple[subprocess.Popen[str], Path]] = []
    for shard_index, shard in enumerate(shards):
        summary_path = temp_summary_dir / f"shard-{shard_index:03d}.json"
        command = [
            _current_python(),
            "-m",
            "vepyr_diffly.worker",
            "compare-bucket-shard",
            "--left-bucket-dir",
            str(left_bucket_dir),
            "--right-bucket-dir",
            str(right_bucket_dir),
            "--bucket-ids",
            ",".join(str(bucket_id) for bucket_id in shard),
            "--csq-fields-json",
            json.dumps(csq_fields),
            "--temp-diff-dir",
            str(temp_diff_dir),
            "--temp-tsv-dir",
            str(temp_tsv_dir),
            "--summary-path",
            str(summary_path),
            "--progress-log",
            "" if reporter is None else str(reporter.log_path),
        ]
        processes.append(
            (
                subprocess.Popen(
                    command,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=env,
                ),
                summary_path,
            )
        )

    aggregate = {
        "left_rows": 0,
        "right_rows": 0,
        "joined_equal_rows": 0,
        "left_only_rows": 0,
        "right_only_rows": 0,
        "unequal_rows": 0,
    }
    diff_paths: list[Path] = []
    completed = 0

    for process, summary_path in processes:
        _, stderr = process.communicate()
        if process.returncode != 0:
            raise RuntimeError(
                f"bucket compare worker failed with exit code {process.returncode}: {stderr}"
            )
        shard_results = _load_shard_results(summary_path)
        for result in shard_results:
            aggregate["left_rows"] += result.left_rows
            aggregate["right_rows"] += result.right_rows
            aggregate["joined_equal_rows"] += result.joined_equal_rows
            aggregate["left_only_rows"] += result.left_only_rows
            aggregate["right_only_rows"] += result.right_only_rows
            aggregate["unequal_rows"] += result.unequal_rows
            diff_paths.append(result.diff_frame_path)
            completed += 1
            if reporter is not None:
                reporter.log(
                    f"consequence: completed bucket {result.bucket_id:04d} "
                    f"({completed}/{len(jobs)}) left={result.left_rows} right={result.right_rows} "
                    f"unequal={result.unequal_rows} left_only={result.left_only_rows} "
                    f"right_only={result.right_only_rows}"
                )

    if reporter is not None:
        reporter.stage(
            f"consequence: merging bucket diff artifacts → {diff_frame_path}",
            tracked_paths=[temp_diff_dir],
        )
    if diff_paths:
        _sink_lazy_frame(pl.scan_parquet([str(path) for path in diff_paths]), diff_frame_path, as_csv=False)
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
            schema=_consequence_schema(csq_fields),
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
            f"buckets_compared={len(jobs)}",
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
    temp_tsv_dir: Path,
    on_bucket_complete: Callable[[int, int, BucketCompareResult], None] | None = None,
) -> list[BucketCompareResult]:
    results: list[BucketCompareResult] = []
    total_buckets = len(bucket_ids)
    for index, bucket_id in enumerate(bucket_ids, start=1):
        left_paths = sorted((left_bucket_dir / f"bucket-{bucket_id:04d}").glob("*.parquet"))
        right_paths = sorted((right_bucket_dir / f"bucket-{bucket_id:04d}").glob("*.parquet"))
        if not left_paths and not right_paths:
            continue
        result = (
            compare_single_bucket(
                bucket_id=bucket_id,
                left_paths=left_paths,
                right_paths=right_paths,
                csq_fields=csq_fields,
                diff_frame_path=temp_diff_dir / f"bucket-{bucket_id:04d}.parquet",
                mismatches_tsv_path=temp_tsv_dir / f"bucket-{bucket_id:04d}.tsv",
            )
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
