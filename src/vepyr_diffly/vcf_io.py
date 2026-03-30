from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import polars as pl

CSQ_DESCRIPTION_RE = re.compile(r'Format: ([^"]+)')


def _is_list_dtype(dtype: Any) -> bool:
    return getattr(dtype, "base_type", lambda: dtype)() == pl.List


def parse_csq_header(vcf_path: Path) -> list[str]:
    with vcf_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("##INFO=<ID=CSQ"):
                continue
            match = CSQ_DESCRIPTION_RE.search(line)
            if match is None:
                raise ValueError(f"CSQ header found in {vcf_path} but format could not be parsed")
            return [field.strip() for field in match.group(1).split("|")]
    raise ValueError(f"VCF {vcf_path} does not contain a CSQ header")


def _scan_vcf_polars_bio(vcf_path: Path) -> pl.LazyFrame:
    import polars_bio as pb

    scan_vcf = getattr(pb, "scan_vcf", None)
    if scan_vcf is None:
        raise RuntimeError("polars-bio does not expose scan_vcf")
    for kwargs in (
        {"info_fields": ["CSQ"]},
        {"include_info": True},
        {},
    ):
        try:
            frame = scan_vcf(str(vcf_path), **kwargs)
        except TypeError:
            continue
        if isinstance(frame, pl.DataFrame):
            return frame.lazy()
        return frame
    frame = scan_vcf(str(vcf_path))
    if isinstance(frame, pl.DataFrame):
        return frame.lazy()
    return frame


def _normalize_scan_columns(frame: pl.LazyFrame) -> pl.LazyFrame:
    schema = frame.collect_schema()
    rename_map: dict[str, str] = {}
    for source, target in (
        ("#CHROM", "chrom"),
        ("CHROM", "chrom"),
        ("chromosome", "chrom"),
        ("POS", "pos"),
        ("position", "pos"),
        ("start", "start"),
        ("end", "end"),
        ("REF", "ref"),
        ("ALT", "alt"),
        ("ID", "id"),
        ("FILTER", "filter"),
        ("QUAL", "qual"),
    ):
        if source in schema and target not in schema:
            rename_map[source] = target
    if rename_map:
        frame = frame.rename(rename_map)
    return frame


def _extract_csq_expr(schema: dict[str, Any]) -> pl.Expr:
    if "CSQ" in schema:
        dtype = schema["CSQ"]
        if _is_list_dtype(dtype):
            return pl.col("CSQ").list.join(",").alias("csq")
        return pl.col("CSQ").cast(pl.String)
    if "INFO" in schema:
        return (
            pl.col("INFO")
            .cast(pl.String)
            .str.extract(r"(?:^|;)CSQ=([^;]+)", group_index=1)
            .alias("csq")
        )
    info_csq_candidates = [name for name in schema if name.lower() in {"info.csq", "info_csq"}]
    if info_csq_candidates:
        candidate = info_csq_candidates[0]
        dtype = schema[candidate]
        if _is_list_dtype(dtype):
            return pl.col(candidate).list.join(",").alias("csq")
        return pl.col(candidate).cast(pl.String).alias("csq")
    raise ValueError("could not find CSQ or INFO column in VCF scan")


def scan_annotated_vcf(vcf_path: Path) -> pl.LazyFrame:
    frame = _normalize_scan_columns(_scan_vcf_polars_bio(vcf_path))
    schema = dict(frame.collect_schema())
    csq_expr = _extract_csq_expr(schema)
    if "pos" in schema:
        pos_expr = pl.col("pos").cast(pl.Int64)
    elif "start" in schema:
        # polars-bio exposes VCF coordinates as genomic interval starts.
        pos_expr = (pl.col("start").cast(pl.Int64) + 1).alias("pos")
    else:
        raise ValueError("could not find POS/start column in VCF scan")
    return frame.with_columns(
        pl.col("chrom").cast(pl.String).alias("chrom"),
        pos_expr.alias("pos"),
        pl.col("ref").cast(pl.String).alias("ref"),
        pl.col("alt").cast(pl.String).alias("alt"),
        (pl.col("id").cast(pl.String) if "id" in schema else pl.lit(".")).alias("id"),
        (pl.col("qual").cast(pl.String) if "qual" in schema else pl.lit(".")).alias("qual"),
        (pl.col("filter").cast(pl.String) if "filter" in schema else pl.lit(".")).alias("filter"),
        csq_expr.alias("csq"),
    ).select("chrom", "pos", "id", "ref", "alt", "qual", "filter", "csq")


def _scan_annotated_vcf_text(vcf_path: Path) -> pl.DataFrame:
    rows: list[dict[str, object]] = []
    with vcf_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = fields[7]
            csq = "."
            for item in info.split(";"):
                if item.startswith("CSQ="):
                    csq = item[4:]
                    break
            rows.append(
                {
                    "chrom": fields[0],
                    "pos": int(fields[1]),
                    "id": fields[2],
                    "ref": fields[3],
                    "alt": fields[4],
                    "qual": fields[5],
                    "filter": fields[6],
                    "csq": csq,
                }
            )
    return pl.DataFrame(
        rows,
        schema={
            "chrom": pl.String,
            "pos": pl.Int64,
            "id": pl.String,
            "ref": pl.String,
            "alt": pl.String,
            "qual": pl.String,
            "filter": pl.String,
            "csq": pl.String,
        },
    )
