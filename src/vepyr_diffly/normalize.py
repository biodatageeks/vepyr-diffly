from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import shutil
from typing import Iterator

import polars as pl

from .progress import ProgressReporter
from .vcf_io import _scan_annotated_vcf_text, parse_csq_header, scan_annotated_vcf


@dataclass(frozen=True)
class NormalizedTables:
    variant: pl.DataFrame
    consequence: pl.DataFrame
    csq_fields: list[str]


VARIANT_KEY = ["chrom", "pos", "ref", "alt"]
SMALL_VCF_TEXT_FALLBACK_BYTES = 4 * 1024 * 1024
STREAMING_CONSEQUENCE_THRESHOLD_BYTES = 8 * 1024 * 1024
CONSEQUENCE_BUCKET_COUNT = 256
CONSEQUENCE_BUCKET_CHUNK_VARIANTS = 200_000
CONSEQUENCE_BUFFERED_VARIANTS_TARGET = 65_536
TEMP_PARQUET_COMPRESSION = "lz4"
TEMP_PARQUET_STATISTICS = False

BASE_VCF_SCHEMA = {
    "chrom": pl.String,
    "pos": pl.Int64,
    "id": pl.String,
    "ref": pl.String,
    "alt": pl.String,
    "qual": pl.String,
    "filter": pl.String,
    "csq": pl.String,
}

VARIANT_CHUNK_SCHEMA = {
    "chrom": pl.String,
    "pos": pl.Int64,
    "ref": pl.String,
    "alt": pl.String,
    "id": pl.String,
    "filter": pl.String,
    "csq_entries": pl.List(pl.String),
}

CONSEQUENCE_CHUNK_SCHEMA = {
    "chrom": pl.String,
    "pos": pl.Int64,
    "ref": pl.String,
    "alt": pl.String,
    "csq_entries": pl.List(pl.String),
}


def _invalid_string_expr(column: str) -> pl.Expr:
    value = pl.col(column).cast(pl.String)
    return pl.col(column).is_null() | value.is_in(["", "None"])


def _clean_string_expr(column: str) -> pl.Expr:
    return (
        pl.when(_invalid_string_expr(column))
        .then(pl.lit("."))
        .otherwise(pl.col(column).cast(pl.String))
        .alias(column)
    )


def _invalid_csq_expr() -> pl.Expr:
    value = pl.col("csq").cast(pl.String)
    return pl.col("csq").is_null() | value.is_in(["", ".", "None"])


def _base_variant_rows(source: pl.LazyFrame) -> pl.LazyFrame:
    return (
        source.with_columns(
            pl.col("chrom").cast(pl.String).alias("chrom"),
            pl.col("pos").cast(pl.Int64).alias("pos"),
            pl.col("ref").cast(pl.String).alias("ref"),
            pl.col("alt").cast(pl.String).str.split(",").alias("alt"),
            _clean_string_expr("id"),
            _clean_string_expr("qual"),
            _clean_string_expr("filter"),
            pl.col("csq").cast(pl.String).alias("csq"),
        )
        .explode("alt")
        .with_columns(_clean_string_expr("alt"))
        .filter(~pl.col("alt").is_in([".", ""]))
    )


def _normalized_allele_expr() -> pl.Expr:
    ref = pl.col("ref").cast(pl.String)
    alt = pl.col("alt").cast(pl.String)
    trimmed = (
        pl.when(
            (ref.str.len_chars() > 0)
            & (alt.str.len_chars() > 0)
            & (ref.str.slice(0, 1) == alt.str.slice(0, 1))
        )
        .then(alt.str.slice(1))
        .otherwise(alt)
    )
    return pl.when(trimmed == "").then(pl.lit("-")).otherwise(trimmed).alias("normalized_allele")


def normalize_alt_for_csq_allele(ref: str, alt: str) -> str:
    ref_value = ref or ""
    alt_value = alt or ""
    if ref_value and alt_value and ref_value[0] == alt_value[0]:
        alt_trimmed = alt_value[1:]
    else:
        alt_trimmed = alt_value
    if alt_trimmed == "":
        return "-"
    return alt_trimmed


def build_variant_summary_lazy(source: pl.LazyFrame) -> pl.LazyFrame:
    consequence_counts = (
        _base_variant_rows(source)
        .with_columns(
            _normalized_allele_expr(),
            pl.when(_invalid_csq_expr())
            .then(pl.lit([], dtype=pl.List(pl.String)))
            .otherwise(pl.col("csq").str.split(","))
            .alias("csq_entry"),
        )
        .explode("csq_entry")
        .filter(pl.col("csq_entry").is_not_null())
        .with_columns(
            pl.col("csq_entry").cast(pl.String).str.split("|").alias("_csq_parts"),
            pl.col("csq_entry")
            .cast(pl.String)
            .str.split("|")
            .list.get(0, null_on_oob=True)
            .cast(pl.String)
            .alias("csq_allele"),
        )
        .filter(pl.col("csq_allele") == pl.col("normalized_allele"))
        .group_by(VARIANT_KEY)
        .agg(pl.len().cast(pl.Int64).alias("consequence_count"))
    )
    variant_base = (
        _base_variant_rows(source)
        .group_by(VARIANT_KEY)
        .agg(
            pl.len().cast(pl.Int64).alias("record_count"),
            pl.col("id").sort().str.join("|").alias("ids"),
            pl.col("filter").sort().str.join("|").alias("filters"),
        )
    )
    return (
        variant_base.join(consequence_counts, on=VARIANT_KEY, how="left")
        .with_columns(
            pl.col("consequence_count").fill_null(0).cast(pl.Int64).alias("consequence_count")
        )
        .select(*VARIANT_KEY, "record_count", "consequence_count", "ids", "filters")
    )


def _csq_field_exprs(csq_fields: list[str]) -> list[pl.Expr]:
    return [
        (
            pl.when(
                pl.col("_csq_parts").list.get(index, null_on_oob=True).is_null()
                | (
                    pl.col("_csq_parts")
                    .list.get(index, null_on_oob=True)
                    .cast(pl.String)
                    .is_in(["", "None"])
                )
            )
            .then(pl.lit("."))
            .otherwise(pl.col("_csq_parts").list.get(index, null_on_oob=True).cast(pl.String))
            .alias(field_name)
        )
        for index, field_name in enumerate(csq_fields)
    ]


def build_consequence_summary_lazy(source: pl.LazyFrame, csq_fields: list[str]) -> pl.LazyFrame:
    consequence_rows = (
        _base_variant_rows(source)
        .with_columns(
            pl.when(_invalid_csq_expr())
            .then(pl.lit([], dtype=pl.List(pl.String)))
            .otherwise(pl.col("csq").str.split(","))
            .alias("csq_entry")
        )
        .explode("csq_entry")
        .filter(pl.col("csq_entry").is_not_null())
        .with_columns(pl.col("csq_entry").cast(pl.String).str.split("|").alias("_csq_parts"))
    )
    group_key = VARIANT_KEY + csq_fields
    return (
        consequence_rows.select(*VARIANT_KEY, *_csq_field_exprs(csq_fields))
        .with_columns(_normalized_allele_expr())
        .filter(pl.col("Allele") == pl.col("normalized_allele"))
        .drop("normalized_allele")
        .group_by(group_key)
        .agg(pl.len().cast(pl.Int64).alias("duplicate_count"))
    )


def _materialize_lazyframe(frame: pl.LazyFrame, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(frame, "sink_parquet"):
        frame.sink_parquet(str(target))
        return
    frame.collect().write_parquet(target)


def _count_rows(path: Path) -> int:
    return pl.scan_parquet(path).select(pl.len()).collect().item()


def _extract_csq_from_info(info: str) -> str:
    for item in info.split(";"):
        if item.startswith("CSQ="):
            return item[4:]
    return "."


def _extract_csq_from_info_bytes(info: bytes) -> str:
    marker = b"CSQ="
    start = info.find(marker)
    if start == -1:
        return "."
    start += len(marker)
    end = info.find(b";", start)
    if end == -1:
        end = len(info)
    return info[start:end].decode("utf-8")


def _extract_csq_entries_from_info_bytes(info: bytes) -> list[str]:
    marker = b"CSQ="
    start = info.find(marker)
    if start == -1:
        return []
    start += len(marker)
    end = info.find(b";", start)
    if end == -1:
        end = len(info)
    payload = info[start:end]
    if not payload:
        return []
    return [item.decode("utf-8") for item in payload.split(b",") if item]


def _count_vcf_records(vcf_path: Path) -> int:
    total = 0
    with vcf_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                total += 1
    return total


def recommend_consequence_bucket_count(*vcf_paths: Path) -> int:
    max_size = max(path.stat().st_size for path in vcf_paths)
    cpu_count = max(1, os.cpu_count() or 1)
    if max_size < 64 * 1024 * 1024:
        return 8
    if max_size < 512 * 1024 * 1024:
        return min(64, max(16, cpu_count * 2))
    if max_size < 2 * 1024 * 1024 * 1024:
        return min(128, max(32, cpu_count * 4))
    if max_size < 8 * 1024 * 1024 * 1024:
        return min(256, max(64, cpu_count * 8))
    return min(512, max(CONSEQUENCE_BUCKET_COUNT, cpu_count * 16))


def recommend_consequence_chunk_variants(*vcf_paths: Path) -> int:
    max_size = max(path.stat().st_size for path in vcf_paths)
    if max_size < 256 * 1024 * 1024:
        return CONSEQUENCE_BUCKET_CHUNK_VARIANTS
    if max_size < 2 * 1024 * 1024 * 1024:
        return 300_000
    if max_size < 8 * 1024 * 1024 * 1024:
        return 400_000
    return 500_000


def _empty_variant_schema(
    *, include_bucket: bool = False, compacted: bool = False
) -> dict[str, pl.DataType]:
    schema: dict[str, pl.DataType] = {
        "chrom": pl.String,
        "pos": pl.Int64,
        "ref": pl.String,
        "alt": pl.String,
    }
    if include_bucket:
        schema["bucket"] = pl.Int64
    schema["record_count"] = pl.Int64
    schema["consequence_count"] = pl.Int64
    if compacted:
        schema["ids"] = pl.String
        schema["filters"] = pl.String
    else:
        schema["id_values"] = pl.List(pl.String)
        schema["filter_values"] = pl.List(pl.String)
    return schema


def _empty_variant_frame(*, include_bucket: bool = False, compacted: bool = False) -> pl.DataFrame:
    return pl.DataFrame(
        schema=_empty_variant_schema(include_bucket=include_bucket, compacted=compacted)
    )


def _empty_consequence_schema(
    csq_fields: list[str], *, include_bucket: bool = False
) -> dict[str, pl.DataType]:
    schema: dict[str, pl.DataType] = {
        "chrom": pl.String,
        "pos": pl.Int64,
        "ref": pl.String,
        "alt": pl.String,
    }
    if include_bucket:
        schema["bucket"] = pl.Int64
    for field in csq_fields:
        schema[field] = pl.String
    schema["duplicate_count"] = pl.Int64
    return schema


def _empty_consequence_frame(
    csq_fields: list[str], *, include_bucket: bool = False
) -> pl.DataFrame:
    return pl.DataFrame(schema=_empty_consequence_schema(csq_fields, include_bucket=include_bucket))


def _iter_vcf_record_chunks(
    vcf_path: Path,
    *,
    chunk_variants: int,
) -> Iterator[tuple[int, pl.DataFrame]]:
    chroms: list[str] = []
    positions: list[int] = []
    refs: list[str] = []
    alts: list[str] = []
    csq_entries: list[list[str]] = []
    processed = 0
    with vcf_path.open("rb", buffering=8 * 1024 * 1024) as handle:
        for raw_line in handle:
            if raw_line.startswith(b"#"):
                continue
            fields = raw_line.rstrip(b"\n").split(b"\t", 8)
            if len(fields) < 8:
                continue
            chroms.append(fields[0].decode("utf-8"))
            positions.append(int(fields[1]))
            refs.append(fields[3].decode("utf-8"))
            alts.append(fields[4].decode("utf-8"))
            csq_entries.append(_extract_csq_entries_from_info_bytes(fields[7]))
            processed += 1
            if len(chroms) >= chunk_variants:
                yield (
                    processed,
                    pl.DataFrame(
                        {
                            "chrom": chroms,
                            "pos": positions,
                            "ref": refs,
                            "alt": alts,
                            "csq_entries": csq_entries,
                        },
                        schema=CONSEQUENCE_CHUNK_SCHEMA,
                    ),
                )
                chroms = []
                positions = []
                refs = []
                alts = []
                csq_entries = []
    if chroms:
        yield (
            processed,
            pl.DataFrame(
                {
                    "chrom": chroms,
                    "pos": positions,
                    "ref": refs,
                    "alt": alts,
                    "csq_entries": csq_entries,
                },
                schema=CONSEQUENCE_CHUNK_SCHEMA,
            ),
        )


def _iter_variant_record_chunks(
    vcf_path: Path,
    *,
    chunk_variants: int,
) -> Iterator[tuple[int, pl.DataFrame]]:
    chroms: list[str] = []
    positions: list[int] = []
    refs: list[str] = []
    alts: list[str] = []
    ids: list[str] = []
    filters: list[str] = []
    csq_entries: list[list[str]] = []
    processed = 0
    with vcf_path.open("rb", buffering=8 * 1024 * 1024) as handle:
        for raw_line in handle:
            if raw_line.startswith(b"#"):
                continue
            fields = raw_line.rstrip(b"\n").split(b"\t", 8)
            if len(fields) < 8:
                continue
            chroms.append(fields[0].decode("utf-8"))
            positions.append(int(fields[1]))
            refs.append(fields[3].decode("utf-8"))
            alts.append(fields[4].decode("utf-8"))
            ids.append(fields[2].decode("utf-8"))
            filters.append(fields[6].decode("utf-8"))
            csq_entries.append(_extract_csq_entries_from_info_bytes(fields[7]))
            processed += 1
            if len(chroms) >= chunk_variants:
                yield (
                    processed,
                    pl.DataFrame(
                        {
                            "chrom": chroms,
                            "pos": positions,
                            "ref": refs,
                            "alt": alts,
                            "id": ids,
                            "filter": filters,
                            "csq_entries": csq_entries,
                        },
                        schema=VARIANT_CHUNK_SCHEMA,
                    ),
                )
                chroms = []
                positions = []
                refs = []
                alts = []
                ids = []
                filters = []
                csq_entries = []
    if chroms:
        yield (
            processed,
            pl.DataFrame(
                {
                    "chrom": chroms,
                    "pos": positions,
                    "ref": refs,
                    "alt": alts,
                    "id": ids,
                    "filter": filters,
                    "csq_entries": csq_entries,
                },
                schema=VARIANT_CHUNK_SCHEMA,
            ),
        )


def _bucket_expr(group_key: list[str], bucket_count: int) -> pl.Expr:
    return (pl.struct(group_key).hash(seed=0) % pl.lit(bucket_count)).cast(pl.Int64).alias("bucket")


def _aggregate_consequence_chunk(
    chunk: pl.DataFrame,
    csq_fields: list[str],
    *,
    bucket_count: int,
) -> pl.DataFrame:
    if chunk.is_empty():
        return _empty_consequence_frame(csq_fields, include_bucket=True)

    group_key = VARIANT_KEY + csq_fields
    return (
        chunk.lazy()
        .with_columns(
            pl.col("alt").cast(pl.String).str.split(",").alias("alt"),
        )
        .explode("alt")
        .with_columns(_clean_string_expr("alt"))
        .filter(~pl.col("alt").is_in([".", ""]))
        .explode("csq_entries")
        .rename({"csq_entries": "csq_entry"})
        .filter(pl.col("csq_entry").is_not_null())
        .with_columns(
            pl.col("csq_entry").cast(pl.String).str.split("|").alias("_csq_parts"),
            _normalized_allele_expr(),
            pl.col("csq_entry")
            .cast(pl.String)
            .str.split("|")
            .list.get(0, null_on_oob=True)
            .cast(pl.String)
            .alias("_csq_allele"),
        )
        .filter(pl.col("_csq_allele") == pl.col("normalized_allele"))
        .select(*VARIANT_KEY, *_csq_field_exprs(csq_fields))
        .with_columns(_bucket_expr(group_key, bucket_count))
        .group_by(["bucket", *group_key])
        .agg(pl.len().cast(pl.Int64).alias("duplicate_count"))
        .collect()
    )


def _aggregate_variant_chunk(
    chunk: pl.DataFrame,
    *,
    bucket_count: int,
) -> pl.DataFrame:
    if chunk.is_empty():
        return _empty_variant_frame(include_bucket=True)

    consequence_counts = (
        chunk.lazy()
        .with_columns(
            pl.col("alt").cast(pl.String).str.split(",").alias("alt"),
        )
        .explode("alt")
        .with_columns(
            _clean_string_expr("alt"),
            _clean_string_expr("id"),
            _clean_string_expr("filter"),
        )
        .filter(~pl.col("alt").is_in([".", ""]))
        .explode("csq_entries")
        .rename({"csq_entries": "csq_entry"})
        .filter(pl.col("csq_entry").is_not_null())
        .with_columns(
            _normalized_allele_expr(),
            pl.col("csq_entry")
            .cast(pl.String)
            .str.split("|")
            .list.get(0, null_on_oob=True)
            .cast(pl.String)
            .alias("_csq_allele"),
        )
        .filter(pl.col("_csq_allele") == pl.col("normalized_allele"))
        .group_by(VARIANT_KEY)
        .agg(pl.len().cast(pl.Int64).alias("consequence_count"))
    )
    variant_base = (
        chunk.lazy()
        .with_columns(
            pl.col("alt").cast(pl.String).str.split(",").alias("alt"),
        )
        .explode("alt")
        .with_columns(
            _clean_string_expr("alt"),
            _clean_string_expr("id"),
            _clean_string_expr("filter"),
        )
        .filter(~pl.col("alt").is_in([".", ""]))
        .group_by(VARIANT_KEY)
        .agg(
            pl.len().cast(pl.Int64).alias("record_count"),
            pl.col("id").sort().alias("id_values"),
            pl.col("filter").sort().alias("filter_values"),
        )
    )
    return (
        variant_base.join(consequence_counts, on=VARIANT_KEY, how="left")
        .with_columns(
            pl.col("consequence_count").fill_null(0).cast(pl.Int64).alias("consequence_count"),
            _bucket_expr(VARIANT_KEY, bucket_count),
        )
        .select(
            "bucket",
            *VARIANT_KEY,
            "record_count",
            "consequence_count",
            "id_values",
            "filter_values",
        )
        .collect()
    )


def _write_bucket_chunk_parts(
    chunk_frame: pl.DataFrame,
    *,
    bucket_root: Path,
    part_index: int,
) -> int:
    if chunk_frame.is_empty():
        return 0
    partitioned = chunk_frame.partition_by(
        "bucket",
        as_dict=True,
        include_key=True,
        maintain_order=False,
    )
    written = 0
    for bucket_key, bucket_frame in partitioned.items():
        if isinstance(bucket_key, tuple):
            bucket_value = int(bucket_key[0])
        else:
            bucket_value = int(bucket_key)
        bucket_dir = bucket_root / f"bucket-{bucket_value:04d}"
        bucket_dir.mkdir(parents=True, exist_ok=True)
        bucket_frame.drop("bucket").write_parquet(
            bucket_dir / f"part-{part_index:05d}.parquet",
            compression=TEMP_PARQUET_COMPRESSION,
            statistics=TEMP_PARQUET_STATISTICS,
        )
        written += 1
    return written


def _merge_buffered_consequence_frames(
    frames: list[pl.DataFrame],
    *,
    csq_fields: list[str],
) -> pl.DataFrame:
    if not frames:
        return _empty_consequence_frame(csq_fields, include_bucket=True)
    if len(frames) == 1:
        return frames[0]
    group_key = ["bucket", *VARIANT_KEY, *csq_fields]
    return (
        pl.concat(frames, how="vertical_relaxed")
        .group_by(group_key)
        .agg(pl.col("duplicate_count").sum().cast(pl.Int64).alias("duplicate_count"))
    )


def _flush_consequence_buffer(
    *,
    frames: list[pl.DataFrame],
    csq_fields: list[str],
    bucket_root: Path,
    part_index: int,
) -> tuple[int, set[int]]:
    merged = _merge_buffered_consequence_frames(frames, csq_fields=csq_fields)
    _write_bucket_chunk_parts(
        merged,
        bucket_root=bucket_root,
        part_index=part_index,
    )
    if merged.is_empty():
        return 0, set()
    non_empty_buckets = {int(value) for value in merged.get_column("bucket").unique().to_list()}
    return 1, non_empty_buckets


def _compact_bucket_parts(
    *,
    bucket_dir: Path,
    csq_fields: list[str],
) -> None:
    part_paths = sorted(bucket_dir.glob("part-*.parquet"))
    if not part_paths:
        return
    group_key = VARIANT_KEY + csq_fields
    compacted_path = bucket_dir / "bucket.parquet"
    (
        pl.scan_parquet([str(path) for path in part_paths])
        .group_by(group_key)
        .agg(pl.col("duplicate_count").sum().cast(pl.Int64).alias("duplicate_count"))
        .sink_parquet(str(compacted_path))
    )
    for part_path in part_paths:
        part_path.unlink()
    _write_bucket_metadata(
        compacted_path=compacted_path,
        count_column="duplicate_count",
    )


def _compact_variant_bucket_parts(
    *,
    bucket_dir: Path,
) -> None:
    part_paths = sorted(bucket_dir.glob("part-*.parquet"))
    if not part_paths:
        return
    compacted_path = bucket_dir / "bucket.parquet"
    compacted = (
        pl.read_parquet([str(path) for path in part_paths])
        .group_by(VARIANT_KEY)
        .agg(
            pl.col("record_count").sum().cast(pl.Int64).alias("record_count"),
            pl.col("consequence_count").sum().cast(pl.Int64).alias("consequence_count"),
            pl.col("id_values").flatten().sort().alias("id_values"),
            pl.col("filter_values").flatten().sort().alias("filter_values"),
        )
        .with_columns(
            pl.col("id_values").list.join("|").alias("ids"),
            pl.col("filter_values").list.join("|").alias("filters"),
        )
        .select(*VARIANT_KEY, "record_count", "consequence_count", "ids", "filters")
    )
    compacted.write_parquet(compacted_path)
    for part_path in part_paths:
        part_path.unlink()
    _write_bucket_metadata(
        compacted_path=compacted_path,
        count_column="record_count",
        extra_sum_columns=["consequence_count"],
    )


def _write_bucket_metadata(
    *,
    compacted_path: Path,
    count_column: str,
    extra_sum_columns: list[str] | None = None,
) -> None:
    extra_sum_columns = extra_sum_columns or []
    frame = pl.scan_parquet(str(compacted_path))
    select_exprs: list[pl.Expr] = [pl.len().cast(pl.Int64).alias("row_count")]
    select_exprs.append(pl.col(count_column).sum().cast(pl.Int64).alias(f"{count_column}_sum"))
    for column in extra_sum_columns:
        select_exprs.append(pl.col(column).sum().cast(pl.Int64).alias(f"{column}_sum"))
    stats = frame.select(*select_exprs).collect().row(0, named=True)
    payload = {
        "row_count": int(stats["row_count"]),
        f"{count_column}_sum": int(stats[f"{count_column}_sum"]),
        "file_size_bytes": int(compacted_path.stat().st_size),
    }
    for column in extra_sum_columns:
        payload[f"{column}_sum"] = int(stats[f"{column}_sum"])
    (compacted_path.with_suffix(".meta.json")).write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )


def _merge_compacted_buckets(
    *,
    bucket_root: Path,
    target_path: Path,
    empty_schema: dict[str, pl.DataType],
) -> None:
    bucket_paths = sorted(bucket_root.glob("bucket-*/bucket.parquet"))
    target_path.parent.mkdir(parents=True, exist_ok=True)
    if not bucket_paths:
        pl.DataFrame(schema=empty_schema).write_parquet(target_path)
        return
    if hasattr(pl.LazyFrame, "sink_parquet"):
        pl.scan_parquet([str(path) for path in bucket_paths]).sink_parquet(str(target_path))
        return
    pl.read_parquet([str(path) for path in bucket_paths]).write_parquet(target_path)


def materialize_consequence_buckets(
    *,
    vcf_path: Path,
    csq_fields: list[str],
    bucket_root: Path,
    reporter: ProgressReporter | None = None,
    side_label: str,
    bucket_count: int = CONSEQUENCE_BUCKET_COUNT,
    chunk_variants: int = CONSEQUENCE_BUCKET_CHUNK_VARIANTS,
    total_variants: int | None = None,
) -> list[int]:
    if bucket_root.exists():
        shutil.rmtree(bucket_root)
    bucket_root.mkdir(parents=True, exist_ok=True)
    if total_variants is None:
        total_variants = _count_vcf_records(vcf_path)
    if reporter is not None:
        reporter.stage(
            f"{side_label}: bucketizing consequence rows → {bucket_root}",
            tracked_paths=[bucket_root],
        )
        reporter.log(f"{side_label}: consequence total input variants={total_variants}")

    non_empty_buckets: set[int] = set()
    part_index = 0
    buffered_frames: list[pl.DataFrame] = []
    buffered_variants = 0
    flush_target_variants = max(chunk_variants, CONSEQUENCE_BUFFERED_VARIANTS_TARGET)
    for processed, chunk in _iter_vcf_record_chunks(
        vcf_path,
        chunk_variants=chunk_variants,
    ):
        chunk_frame = _aggregate_consequence_chunk(
            chunk,
            csq_fields,
            bucket_count=bucket_count,
        )
        if not chunk_frame.is_empty():
            buffered_frames.append(chunk_frame)
            buffered_variants += chunk.height
        if buffered_variants >= flush_target_variants and buffered_frames:
            written_parts, flushed_buckets = _flush_consequence_buffer(
                frames=buffered_frames,
                csq_fields=csq_fields,
                bucket_root=bucket_root,
                part_index=part_index,
            )
            non_empty_buckets.update(flushed_buckets)
            part_index += written_parts
            buffered_frames = []
            buffered_variants = 0
        if reporter is not None:
            percentage = (processed / total_variants * 100.0) if total_variants else 0.0
            reporter.log(
                f"{side_label}: consequence progress {processed}/{total_variants} "
                f"variants ({percentage:.2f}%), chunk_parts={part_index}, "
                f"buckets={len(non_empty_buckets)}, chunk_variants={chunk_variants}, "
                f"buffered_variants={buffered_variants}"
            )

    if buffered_frames:
        written_parts, flushed_buckets = _flush_consequence_buffer(
            frames=buffered_frames,
            csq_fields=csq_fields,
            bucket_root=bucket_root,
            part_index=part_index,
        )
        non_empty_buckets.update(flushed_buckets)
        part_index += written_parts

    if non_empty_buckets:
        if reporter is not None:
            reporter.stage(
                f"{side_label}: compacting consequence bucket parts",
                tracked_paths=[bucket_root],
            )
        for bucket_id in sorted(non_empty_buckets):
            _compact_bucket_parts(
                bucket_dir=bucket_root / f"bucket-{bucket_id:04d}",
                csq_fields=csq_fields,
            )

    if reporter is not None:
        reporter.log(
            f"{side_label}: bucketized consequence rows into {len(non_empty_buckets)} buckets"
        )
    return sorted(non_empty_buckets)


def materialize_variant_buckets(
    *,
    vcf_path: Path,
    bucket_root: Path,
    variant_path: Path,
    reporter: ProgressReporter | None = None,
    side_label: str,
    bucket_count: int = CONSEQUENCE_BUCKET_COUNT,
    chunk_variants: int = 10_000,
    total_variants: int | None = None,
) -> list[int]:
    if bucket_root.exists():
        shutil.rmtree(bucket_root)
    bucket_root.mkdir(parents=True, exist_ok=True)
    if total_variants is None:
        total_variants = _count_vcf_records(vcf_path)
    if reporter is not None:
        reporter.stage(
            f"{side_label}: bucketizing variant rows → {bucket_root}",
            tracked_paths=[bucket_root],
        )
        reporter.log(f"{side_label}: variant total input variants={total_variants}")

    non_empty_buckets: set[int] = set()
    part_index = 0
    for processed, chunk in _iter_variant_record_chunks(
        vcf_path,
        chunk_variants=chunk_variants,
    ):
        chunk_frame = _aggregate_variant_chunk(
            chunk,
            bucket_count=bucket_count,
        )
        _write_bucket_chunk_parts(
            chunk_frame,
            bucket_root=bucket_root,
            part_index=part_index,
        )
        if not chunk_frame.is_empty():
            non_empty_buckets.update(
                int(value) for value in chunk_frame.get_column("bucket").unique().to_list()
            )
        part_index += 1
        if reporter is not None:
            percentage = (processed / total_variants * 100.0) if total_variants else 0.0
            reporter.log(
                f"{side_label}: variant progress {processed}/{total_variants} "
                f"variants ({percentage:.2f}%), chunk_parts={part_index}, "
                f"buckets={len(non_empty_buckets)}, chunk_variants={chunk_variants}"
            )

    if non_empty_buckets:
        if reporter is not None:
            reporter.stage(
                f"{side_label}: compacting variant bucket parts",
                tracked_paths=[bucket_root],
            )
        for bucket_id in sorted(non_empty_buckets):
            _compact_variant_bucket_parts(bucket_dir=bucket_root / f"bucket-{bucket_id:04d}")

    _merge_compacted_buckets(
        bucket_root=bucket_root,
        target_path=variant_path,
        empty_schema=_empty_variant_schema(compacted=True),
    )
    if reporter is not None:
        reporter.log(f"{side_label}: variant summary rows={_count_rows(variant_path)}")
        reporter.log(f"{side_label}: bucketized variant rows into {len(non_empty_buckets)} buckets")
    return sorted(non_empty_buckets)


def _eager_normalize_fallback(vcf_path: Path, csq_fields: list[str]) -> NormalizedTables:
    source = _scan_annotated_vcf_text(vcf_path)
    variant = build_variant_summary_lazy(source.lazy()).collect()
    consequence = build_consequence_summary_lazy(source.lazy(), csq_fields).collect()
    return NormalizedTables(variant=variant, consequence=consequence, csq_fields=csq_fields)


def materialize_variant_summary(
    *,
    vcf_path: Path,
    variant_path: Path,
    bucket_root: Path | None = None,
    reporter: ProgressReporter | None = None,
    side_label: str,
    bucket_count: int = CONSEQUENCE_BUCKET_COUNT,
    chunk_variants: int = 10_000,
    total_variants: int | None = None,
) -> list[str]:
    csq_fields = parse_csq_header(vcf_path)
    if reporter is not None:
        reporter.log(f"{side_label}: parsed CSQ header with {len(csq_fields)} fields")
    if vcf_path.stat().st_size <= SMALL_VCF_TEXT_FALLBACK_BYTES:
        if reporter is not None:
            reporter.stage(
                f"{side_label}: using eager text normalization for small file {vcf_path.name}"
            )
        eager = _eager_normalize_fallback(vcf_path, csq_fields)
        eager.variant.write_parquet(variant_path)
        if bucket_root is not None:
            if bucket_root.exists():
                shutil.rmtree(bucket_root)
            bucket_root.mkdir(parents=True, exist_ok=True)
            partitioned = eager.variant.with_columns(
                _bucket_expr(VARIANT_KEY, bucket_count)
            ).partition_by("bucket", as_dict=True, include_key=True, maintain_order=False)
            for bucket_key, bucket_frame in partitioned.items():
                bucket_value = int(bucket_key[0] if isinstance(bucket_key, tuple) else bucket_key)
                bucket_dir = bucket_root / f"bucket-{bucket_value:04d}"
                bucket_dir.mkdir(parents=True, exist_ok=True)
                compacted_path = bucket_dir / "bucket.parquet"
                bucket_frame.drop("bucket").write_parquet(compacted_path)
                _write_bucket_metadata(
                    compacted_path=compacted_path,
                    count_column="record_count",
                    extra_sum_columns=["consequence_count"],
                )
        if reporter is not None:
            reporter.log(f"{side_label}: variant summary rows={eager.variant.height}")
        return csq_fields

    if bucket_root is not None:
        materialize_variant_buckets(
            vcf_path=vcf_path,
            bucket_root=bucket_root,
            variant_path=variant_path,
            reporter=reporter,
            side_label=side_label,
            bucket_count=bucket_count,
            chunk_variants=chunk_variants,
            total_variants=total_variants,
        )
        return csq_fields

    if reporter is not None:
        reporter.stage(
            f"{side_label}: materializing variant summary → {variant_path}",
            tracked_paths=[variant_path],
        )
    try:
        _materialize_lazyframe(
            build_variant_summary_lazy(scan_annotated_vcf(vcf_path)), variant_path
        )
    except Exception as exc:
        if reporter is not None:
            reporter.log(
                f"{side_label}: lazy variant normalization failed, falling back to eager text parser ({exc})"
            )
        eager = _eager_normalize_fallback(vcf_path, csq_fields)
        eager.variant.write_parquet(variant_path)
    if reporter is not None:
        reporter.log(f"{side_label}: variant summary rows={_count_rows(variant_path)}")
    return csq_fields


def materialize_consequence_summary(
    *,
    vcf_path: Path,
    consequence_path: Path,
    csq_fields: list[str],
    reporter: ProgressReporter | None = None,
    side_label: str,
) -> None:
    if vcf_path.stat().st_size <= SMALL_VCF_TEXT_FALLBACK_BYTES:
        eager = _eager_normalize_fallback(vcf_path, csq_fields)
        consequence_path.parent.mkdir(parents=True, exist_ok=True)
        eager.consequence.write_parquet(consequence_path)
        if reporter is not None:
            reporter.log(f"{side_label}: consequence summary rows={eager.consequence.height}")
        return

    if reporter is not None:
        reporter.stage(
            f"{side_label}: materializing consequence summary → {consequence_path}",
            tracked_paths=[consequence_path],
        )
    try:
        _materialize_lazyframe(
            build_consequence_summary_lazy(scan_annotated_vcf(vcf_path), csq_fields),
            consequence_path,
        )
    except Exception as exc:
        if reporter is not None:
            reporter.log(
                f"{side_label}: lazy consequence normalization failed, falling back to eager text parser ({exc})"
            )
        eager = _eager_normalize_fallback(vcf_path, csq_fields)
        eager.consequence.write_parquet(consequence_path)
    if reporter is not None:
        reporter.log(f"{side_label}: consequence summary rows={_count_rows(consequence_path)}")


def materialize_normalized_tables(
    *,
    vcf_path: Path,
    variant_path: Path,
    consequence_path: Path,
    reporter: ProgressReporter | None = None,
    side_label: str,
) -> list[str]:
    csq_fields = materialize_variant_summary(
        vcf_path=vcf_path,
        variant_path=variant_path,
        reporter=reporter,
        side_label=side_label,
    )
    if vcf_path.stat().st_size >= STREAMING_CONSEQUENCE_THRESHOLD_BYTES:
        bucket_root = consequence_path.parent / f"{consequence_path.stem}.buckets"
        materialize_consequence_buckets(
            vcf_path=vcf_path,
            csq_fields=csq_fields,
            bucket_root=bucket_root,
            reporter=reporter,
            side_label=side_label,
        )
        return csq_fields

    materialize_consequence_summary(
        vcf_path=vcf_path,
        consequence_path=consequence_path,
        csq_fields=csq_fields,
        reporter=reporter,
        side_label=side_label,
    )
    return csq_fields


def normalize_annotated_vcf(vcf_path: Path) -> NormalizedTables:
    csq_fields = parse_csq_header(vcf_path)
    if vcf_path.stat().st_size <= SMALL_VCF_TEXT_FALLBACK_BYTES:
        return _eager_normalize_fallback(vcf_path, csq_fields)
    try:
        variant = build_variant_summary_lazy(scan_annotated_vcf(vcf_path)).collect()
        consequence = build_consequence_summary_lazy(
            scan_annotated_vcf(vcf_path), csq_fields
        ).collect()
        return NormalizedTables(
            variant=variant,
            consequence=consequence,
            csq_fields=csq_fields,
        )
    except Exception:
        return _eager_normalize_fallback(vcf_path, csq_fields)
