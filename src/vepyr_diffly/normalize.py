from __future__ import annotations

from dataclasses import dataclass
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
    return (
        pl.struct(["ref", "alt"])
        .map_elements(
            lambda value: normalize_alt_for_csq_allele(
                value["ref"],
                value["alt"],
            ),
            return_dtype=pl.String,
        )
        .alias("normalized_allele")
    )


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
            .alias("csq_entry")
        )
        .explode("csq_entry")
        .filter(pl.col("csq_entry").is_not_null())
        .with_columns(
            pl.col("csq_entry").cast(pl.String).str.split("|").alias("_csq_parts"),
            pl.col("csq_entry").cast(pl.String).str.split("|").list.get(0, null_on_oob=True).cast(pl.String).alias("csq_allele"),
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
            .otherwise(
                pl.col("_csq_parts").list.get(index, null_on_oob=True).cast(pl.String)
            )
            .alias(field_name)
        )
        for index, field_name in enumerate(csq_fields)
    ]


def build_consequence_summary_lazy(
    source: pl.LazyFrame, csq_fields: list[str]
) -> pl.LazyFrame:
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
    rows: list[dict[str, object]] = []
    processed = 0
    with vcf_path.open("rb", buffering=8 * 1024 * 1024) as handle:
        for raw_line in handle:
            if raw_line.startswith(b"#"):
                continue
            fields = raw_line.rstrip(b"\n").split(b"\t", 8)
            if len(fields) < 8:
                continue
            rows.append(
                {
                    "chrom": fields[0].decode("utf-8"),
                    "pos": int(fields[1]),
                    "id": fields[2].decode("utf-8"),
                    "ref": fields[3].decode("utf-8"),
                    "alt": fields[4].decode("utf-8"),
                    "qual": fields[5].decode("utf-8"),
                    "filter": fields[6].decode("utf-8"),
                    "csq": _extract_csq_from_info_bytes(fields[7]),
                }
            )
            processed += 1
            if len(rows) >= chunk_variants:
                yield processed, pl.DataFrame(rows, schema=BASE_VCF_SCHEMA)
                rows = []
    if rows:
        yield processed, pl.DataFrame(rows, schema=BASE_VCF_SCHEMA)


def _bucket_expr(group_key: list[str], bucket_count: int) -> pl.Expr:
    return (pl.struct(group_key).hash(seed=0) % pl.lit(bucket_count)).cast(pl.Int64).alias(
        "bucket"
    )


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
            pl.when(_invalid_csq_expr())
            .then(pl.lit([], dtype=pl.List(pl.String)))
            .otherwise(pl.col("csq").cast(pl.String).str.split(","))
            .alias("csq_entry"),
        )
        .explode("alt")
        .with_columns(_clean_string_expr("alt"))
        .filter(~pl.col("alt").is_in([".", ""]))
        .explode("csq_entry")
        .filter(pl.col("csq_entry").is_not_null())
        .with_columns(pl.col("csq_entry").cast(pl.String).str.split("|").alias("_csq_parts"))
        .select(*VARIANT_KEY, *_csq_field_exprs(csq_fields))
        .with_columns(_normalized_allele_expr())
        .filter(pl.col("Allele") == pl.col("normalized_allele"))
        .drop("normalized_allele")
        .with_columns(_bucket_expr(group_key, bucket_count))
        .group_by(["bucket", *group_key])
        .agg(pl.len().cast(pl.Int64).alias("duplicate_count"))
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
        bucket_frame.drop("bucket").write_parquet(bucket_dir / f"part-{part_index:05d}.parquet")
        written += 1
    return written


def materialize_consequence_buckets(
    *,
    vcf_path: Path,
    csq_fields: list[str],
    bucket_root: Path,
    reporter: ProgressReporter | None = None,
    side_label: str,
    bucket_count: int = CONSEQUENCE_BUCKET_COUNT,
) -> list[int]:
    if bucket_root.exists():
        shutil.rmtree(bucket_root)
    bucket_root.mkdir(parents=True, exist_ok=True)
    total_variants = _count_vcf_records(vcf_path)
    if reporter is not None:
        reporter.stage(
            f"{side_label}: bucketizing consequence rows → {bucket_root}",
            tracked_paths=[bucket_root],
        )
        reporter.log(f"{side_label}: consequence total input variants={total_variants}")

    non_empty_buckets: set[int] = set()
    part_index = 0
    for processed, chunk in _iter_vcf_record_chunks(
        vcf_path,
        chunk_variants=CONSEQUENCE_BUCKET_CHUNK_VARIANTS,
    ):
        chunk_frame = _aggregate_consequence_chunk(
            chunk,
            csq_fields,
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
                f"{side_label}: consequence progress {processed}/{total_variants} "
                f"variants ({percentage:.2f}%), chunk_parts={part_index}, "
                f"buckets={len(non_empty_buckets)}"
            )

    if reporter is not None:
        reporter.log(
            f"{side_label}: bucketized consequence rows into {len(non_empty_buckets)} buckets"
        )
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
    reporter: ProgressReporter | None = None,
    side_label: str,
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
        if reporter is not None:
            reporter.log(f"{side_label}: variant summary rows={eager.variant.height}")
        return csq_fields

    if reporter is not None:
        reporter.stage(
            f"{side_label}: materializing variant summary → {variant_path}",
            tracked_paths=[variant_path],
        )
    try:
        _materialize_lazyframe(build_variant_summary_lazy(scan_annotated_vcf(vcf_path)), variant_path)
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
