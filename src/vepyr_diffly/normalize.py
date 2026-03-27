from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from .vcf_io import _scan_annotated_vcf_text, parse_csq_header, scan_annotated_vcf


@dataclass(frozen=True)
class NormalizedTables:
    variant: pl.DataFrame
    consequence: pl.DataFrame
    csq_fields: list[str]


VARIANT_KEY = ["chrom", "pos", "ref", "alt"]


def _clean_string(value: str | None) -> str:
    if value is None or value == "":
        return "."
    return value


def _split_csv(value: str | None) -> list[str]:
    if value is None or value in {"", "."}:
        return []
    return [item for item in value.split(",") if item]


def _build_csq_rows(
    row: dict[str, object], csq_fields: list[str], field_count: int
) -> list[dict[str, object]]:
    csq_value = row.get("csq")
    csq_entries = [] if csq_value in {None, "", "."} else str(csq_value).split(",")
    rows: list[dict[str, object]] = []
    for entry in csq_entries:
        values = entry.split("|")
        if len(values) < field_count:
            values.extend([""] * (field_count - len(values)))
        payload: dict[str, object] = {
            "chrom": row["chrom"],
            "pos": row["pos"],
            "ref": row["ref"],
            "alt": row["alt"],
        }
        for field_name, value in zip(csq_fields, values, strict=False):
            payload[field_name] = _clean_string(value)
        rows.append(payload)
    return rows


def _iter_variant_rows(frame: pl.DataFrame) -> Iterable[dict[str, object]]:
    for row in frame.iter_rows(named=True):
        for alt in _split_csv(str(row["alt"])):
            yield {
                "chrom": str(row["chrom"]),
                "pos": int(row["pos"]),
                "ref": str(row["ref"]),
                "alt": alt,
                "id": _clean_string(str(row["id"])),
                "qual": _clean_string(str(row["qual"])),
                "filter": _clean_string(str(row["filter"])),
                "csq": row["csq"],
            }


def normalize_annotated_vcf(vcf_path: Path) -> NormalizedTables:
    csq_fields = parse_csq_header(vcf_path)
    try:
        source = scan_annotated_vcf(vcf_path).collect()
    except Exception:
        source = _scan_annotated_vcf_text(vcf_path)

    variant_rows = list(_iter_variant_rows(source))
    variant_frame = pl.DataFrame(variant_rows or [], schema={
        "chrom": pl.String,
        "pos": pl.Int64,
        "ref": pl.String,
        "alt": pl.String,
        "id": pl.String,
        "qual": pl.String,
        "filter": pl.String,
        "csq": pl.String,
    })
    if variant_frame.is_empty():
        consequence_frame = pl.DataFrame()
    else:
        consequence_rows: list[dict[str, object]] = []
        for row in variant_rows:
            consequence_rows.extend(_build_csq_rows(row, csq_fields, len(csq_fields)))
        consequence_frame = (
            pl.DataFrame(consequence_rows)
            if consequence_rows
            else pl.DataFrame(schema={**{k: variant_frame.schema[k] for k in VARIANT_KEY}, **{name: pl.String for name in csq_fields}})
        )

    variant_summary = (
        variant_frame.with_columns(
            pl.col("csq")
            .map_elements(lambda value: 0 if value in {None, "", "."} else len(str(value).split(",")), return_dtype=pl.Int64)
            .alias("csq_entry_count")
        )
        .group_by(VARIANT_KEY)
        .agg(
            pl.len().alias("record_count"),
            pl.col("csq_entry_count").sum().alias("consequence_count"),
            pl.col("id").sort().str.join("|").alias("ids"),
            pl.col("filter").sort().str.join("|").alias("filters"),
        )
        .sort(VARIANT_KEY)
        .with_columns(
            pl.col("record_count").cast(pl.Int64),
            pl.col("consequence_count").cast(pl.Int64),
        )
    )

    if consequence_frame.is_empty():
        consequence_summary = pl.DataFrame(
            schema={**{key: variant_summary.schema.get(key, pl.String) for key in VARIANT_KEY}}
        )
    else:
        group_key = VARIANT_KEY + csq_fields
        consequence_summary = (
            consequence_frame.group_by(group_key)
            .agg(pl.len().alias("duplicate_count"))
            .sort(group_key)
            .with_columns(pl.col("duplicate_count").cast(pl.Int64))
        )

    return NormalizedTables(
        variant=variant_summary,
        consequence=consequence_summary,
        csq_fields=csq_fields,
    )
