from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


DEFAULT_KEY_FIELDS = ["chrom", "pos", "ref", "alt", "Feature"]
IGNORED_DIFF_FIELDS = {"diff_kind", "side", "duplicate_count_left", "duplicate_count_right"}
MISSING_VALUE_MARKERS = {"", "."}


def _diff_kind(row: dict[str, str]) -> str:
    return row.get("diff_kind") or row.get("side") or ""


def _is_missing(value: str) -> bool:
    return value in MISSING_VALUE_MARKERS


def _classify_pair(
    left: dict[str, str],
    right: dict[str, str],
    diff_fields: list[str],
) -> list[str]:
    categories: list[str] = []

    if any(field in diff_fields for field in ("HGVSc", "HGVSp")):
        hgvs_missing_only = True
        for field in ("HGVSc", "HGVSp"):
            if field not in diff_fields:
                continue
            if not (_is_missing(left.get(field, "")) ^ _is_missing(right.get(field, ""))):
                hgvs_missing_only = False
                break
        categories.append("missing_hgvs" if hgvs_missing_only else "hgvs_payload_drift")

    if "HGNC_ID" in diff_fields:
        categories.append(
            "missing_hgnc_id"
            if _is_missing(left.get("HGNC_ID", "")) ^ _is_missing(right.get("HGNC_ID", ""))
            else "hgnc_id_drift"
        )

    if "Consequence" in diff_fields and "IMPACT" in diff_fields:
        categories.append("consequence_and_impact_reclassification")
    elif "Consequence" in diff_fields:
        categories.append("consequence_reclassification")
    elif "IMPACT" in diff_fields:
        categories.append("impact_reclassification")

    if any(
        field in diff_fields
        for field in ("CDS_position", "Protein_position", "Amino_acids", "Codons", "DOMAINS")
    ):
        categories.append("protein_annotation_drift")

    if any(field in diff_fields for field in ("SIFT", "PolyPhen")):
        categories.append("predictor_drift")

    if not categories:
        categories.append("other_payload_drift")
    return categories


def _build_pair_example(
    key_fields: list[str],
    group_key: tuple[str, ...],
    diff_fields: list[str],
    categories: list[str],
    left: dict[str, str],
    right: dict[str, str],
    fieldnames: list[str],
) -> dict[str, Any]:
    return {
        "key": dict(zip(key_fields, group_key, strict=False)),
        "diff_fields": diff_fields,
        "categories": categories,
        "left": {
            field: left.get(field, "")
            for field in ["Consequence", "IMPACT", "HGNC_ID", "HGVSc", "HGVSp"]
            if field in fieldnames
        },
        "right": {
            field: right.get(field, "")
            for field in ["Consequence", "IMPACT", "HGNC_ID", "HGVSc", "HGVSp"]
            if field in fieldnames
        },
    }


def analyze_consequence_mismatches(
    mismatches_tsv: Path,
    *,
    key_fields: list[str] | None = None,
    top_n: int = 10,
    example_limit: int = 8,
) -> dict[str, Any]:
    key_fields = key_fields or DEFAULT_KEY_FIELDS
    with mismatches_tsv.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    grouped: dict[tuple[str, ...], list[dict[str, str]]] = {}
    for row in rows:
        key = tuple(row.get(field, "") for field in key_fields)
        grouped.setdefault(key, []).append(row)

    field_counts: dict[str, int] = {}
    signature_counts: dict[str, int] = {}
    category_counts: dict[str, int] = {}
    paired_examples: list[dict[str, Any]] = []
    category_examples: dict[str, list[dict[str, Any]]] = {}
    unpaired_examples: list[dict[str, Any]] = []
    paired_keys = 0
    left_only_rows = 0
    right_only_rows = 0
    unpaired_rows = 0

    for group_key, group_rows in grouped.items():
        left_rows = [row for row in group_rows if _diff_kind(row) == "left_only"]
        right_rows = [row for row in group_rows if _diff_kind(row) == "right_only"]
        left_only_rows += len(left_rows)
        right_only_rows += len(right_rows)

        if len(left_rows) == 1 and len(right_rows) == 1 and len(group_rows) == 2:
            paired_keys += 1
            left = left_rows[0]
            right = right_rows[0]
            diff_fields = [
                field
                for field in fieldnames
                if field not in IGNORED_DIFF_FIELDS and left.get(field, "") != right.get(field, "")
            ]
            for field in diff_fields:
                field_counts[field] = field_counts.get(field, 0) + 1
            signature = ",".join(diff_fields) if diff_fields else "<no_field_diff>"
            signature_counts[signature] = signature_counts.get(signature, 0) + 1
            categories = _classify_pair(left, right, diff_fields)
            for category in categories:
                category_counts[category] = category_counts.get(category, 0) + 1
            if len(paired_examples) < example_limit:
                paired_examples.append(
                    _build_pair_example(
                        key_fields, group_key, diff_fields, categories, left, right, fieldnames
                    )
                )
            for category in categories:
                category_bucket = category_examples.setdefault(category, [])
                if len(category_bucket) < example_limit:
                    category_bucket.append(
                        _build_pair_example(
                            key_fields,
                            group_key,
                            diff_fields,
                            categories,
                            left,
                            right,
                            fieldnames,
                        )
                    )
            continue

        unpaired_rows += len(group_rows)
        if len(unpaired_examples) < example_limit:
            unpaired_examples.append(
                {
                    "key": dict(zip(key_fields, group_key, strict=False)),
                    "row_count": len(group_rows),
                    "diff_kinds": sorted(_diff_kind(row) for row in group_rows),
                }
            )

    top_differing_fields = [
        {"field": field, "paired_keys": count}
        for field, count in sorted(field_counts.items(), key=lambda item: (-item[1], item[0]))[:top_n]
    ]
    top_signatures = [
        {"fields": signature.split(",") if signature != "<no_field_diff>" else [], "paired_keys": count}
        for signature, count in sorted(
            signature_counts.items(),
            key=lambda item: (-item[1], item[0]),
        )[:top_n]
    ]
    top_categories = [
        {"category": category, "paired_keys": count}
        for category, count in sorted(
            category_counts.items(),
            key=lambda item: (-item[1], item[0]),
        )[:top_n]
    ]

    return {
        "mismatches_tsv": str(mismatches_tsv),
        "key_fields": key_fields,
        "total_rows": len(rows),
        "grouped_keys": len(grouped),
        "paired_keys": paired_keys,
        "left_only_rows": left_only_rows,
        "right_only_rows": right_only_rows,
        "unpaired_rows": unpaired_rows,
        "top_differing_fields": top_differing_fields,
        "top_field_signatures": top_signatures,
        "top_semantic_categories": top_categories,
        "paired_examples": paired_examples,
        "category_examples": category_examples,
        "unpaired_examples": unpaired_examples,
    }


def write_mismatch_analysis(output_path: Path, payload: dict[str, Any]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def extract_csq_examples_from_analysis(
    analysis: dict[str, Any],
    *,
    left_vcf: Path,
    right_vcf: Path,
    per_category_limit: int = 3,
) -> dict[str, Any]:
    wanted: dict[str, dict[tuple[str, str, str, str], set[str]]] = {}
    for category, examples in analysis.get("category_examples", {}).items():
        selected = examples[:per_category_limit]
        if not selected:
            continue
        category_cases: dict[tuple[str, str, str, str], set[str]] = {}
        for example in selected:
            key = example["key"]
            locus = (
                str(key["chrom"]),
                str(key["pos"]),
                str(key["ref"]),
                str(key["alt"]),
            )
            category_cases.setdefault(locus, set()).add(str(key["Feature"]))
        wanted[category] = category_cases

    left_entries = _collect_csq_entries(left_vcf, wanted)
    right_entries = _collect_csq_entries(right_vcf, wanted)
    categories: dict[str, list[dict[str, Any]]] = {}
    for category, examples in analysis.get("category_examples", {}).items():
        selected = examples[:per_category_limit]
        if not selected:
            continue
        enriched: list[dict[str, Any]] = []
        for example in selected:
            key = example["key"]
            locus = (
                str(key["chrom"]),
                str(key["pos"]),
                str(key["ref"]),
                str(key["alt"]),
            )
            feature = str(key["Feature"])
            enriched.append(
                {
                    **example,
                    "left_raw_csq": left_entries.get((category, locus, feature), []),
                    "right_raw_csq": right_entries.get((category, locus, feature), []),
                }
            )
        categories[category] = enriched

    return {
        "left_vcf": str(left_vcf),
        "right_vcf": str(right_vcf),
        "per_category_limit": per_category_limit,
        "categories": categories,
    }


def write_csq_examples(output_path: Path, payload: dict[str, Any]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _collect_csq_entries(
    vcf_path: Path,
    wanted: dict[str, dict[tuple[str, str, str, str], set[str]]],
) -> dict[tuple[str, tuple[str, str, str, str], str], list[str]]:
    matches: dict[tuple[str, tuple[str, str, str, str], str], list[str]] = {}
    if not wanted:
        return matches
    by_locus: dict[tuple[str, str, str, str], list[tuple[str, set[str]]]] = {}
    for category, loci in wanted.items():
        for locus, features in loci.items():
            by_locus.setdefault(locus, []).append((category, features))

    with vcf_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 8:
                continue
            locus = (columns[0], columns[1], columns[3], columns[4])
            needed = by_locus.get(locus)
            if not needed:
                continue
            csq_payload = _extract_csq_payload(columns[7])
            if not csq_payload:
                continue
            for entry in csq_payload.split(","):
                for category, features in needed:
                    for feature in features:
                        if feature in entry:
                            matches.setdefault((category, locus, feature), []).append(entry)
    return matches


def _extract_csq_payload(info: str) -> str:
    for item in info.split(";"):
        if item.startswith("CSQ="):
            return item[4:]
    return ""
