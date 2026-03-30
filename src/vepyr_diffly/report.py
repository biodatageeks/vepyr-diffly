from __future__ import annotations

import json
from typing import Any

try:
    from rich.console import Console
    from rich.table import Table
except ModuleNotFoundError:  # pragma: no cover

    class Console:  # type: ignore[override]
        def print(self, value: object = "") -> None:
            print(value)

    class Table:  # pragma: no cover
        def __init__(self, title: str | None = None) -> None:
            self.title = title
            self.columns: list[str] = []
            self.rows: list[list[str]] = []

        def add_column(self, name: str) -> None:
            self.columns.append(name)

        def add_row(self, *values: str) -> None:
            self.rows.append(list(values))

        def __str__(self) -> str:
            lines = []
            if self.title:
                lines.append(self.title)
            if self.columns:
                lines.append("\t".join(self.columns))
            lines.extend("\t".join(row) for row in self.rows)
            return "\n".join(lines)


from .models import RunArtifacts, RuntimeConfig, TierResult


def print_run_summary(
    *,
    console: Console,
    config: RuntimeConfig,
    variant: TierResult,
    consequence: TierResult,
    left_vcf: str | Any,
    right_vcf: str | Any,
    progress_log_path: str | Any,
    resource_plan: dict[str, Any] | None = None,
    timings: dict[str, float] | None = None,
    chromosome_summary: dict[str, Any] | None = None,
) -> None:
    console.print(f"Preset: [bold]{config.preset.name}[/bold]")
    console.print(f"Input: {config.input_vcf}")
    console.print(f"Left annotated VCF: {left_vcf}")
    console.print(f"Right annotated VCF: {right_vcf}")
    console.print(
        f"Sample first N: {config.sample_first_n if config.sample_first_n is not None else 'full'}"
    )
    console.print(
        "Chromosomes: "
        + (
            ",".join(config.selected_chromosomes)
            if config.selected_chromosomes
            else "all"
        )
    )
    console.print(f"Progress log: {progress_log_path}")
    if resource_plan:
        console.print(
            "Resource plan: "
            + ", ".join(f"{name}={value}" for name, value in resource_plan.items())
        )
    if timings:
        console.print(
            "Timings (s): " + ", ".join(f"{name}={value:.3f}" for name, value in timings.items())
        )
    console.print("")

    table = Table(title="Comparison Summary")
    table.add_column("Tier")
    table.add_column("Equal")
    table.add_column("Left only")
    table.add_column("Right only")
    table.add_column("Unequal")
    table.add_column("Joined equal")
    for tier in (variant, consequence):
        table.add_row(
            tier.name,
            "yes" if tier.equal else "no",
            str(tier.left_only_rows),
            str(tier.right_only_rows),
            str(tier.unequal_rows),
            str(tier.joined_equal_rows),
        )
    console.print(table)
    if chromosome_summary and chromosome_summary.get("per_chromosome"):
        console.print("")
        chrom_table = Table(title="Per-Chromosome Summary")
        chrom_table.add_column("Chrom")
        chrom_table.add_column("Variant equal")
        chrom_table.add_column("Variant left/right/uneq")
        chrom_table.add_column("Consequence equal")
        chrom_table.add_column("Consequence left/right/uneq")
        chrom_table.add_column("Attributed timings (s)")
        for chrom, payload in chromosome_summary["per_chromosome"].items():
            variant_payload = payload["variant"]
            consequence_payload = payload["consequence"]
            timing_payload = payload.get("timings", {})
            chrom_table.add_row(
                chrom,
                "yes" if variant_payload["equal"] else "no",
                f"{variant_payload['left_only_rows']}/{variant_payload['right_only_rows']}/{variant_payload['unequal_rows']}",
                "yes" if consequence_payload["equal"] else "no",
                f"{consequence_payload['left_only_rows']}/{consequence_payload['right_only_rows']}/{consequence_payload['unequal_rows']}",
                ", ".join(
                    f"{name}={value:.3f}" for name, value in timing_payload.items() if value > 0
                )
                or "n/a",
            )
        console.print(chrom_table)


def write_run_summary(
    *,
    config: RuntimeConfig,
    artifacts: RunArtifacts,
    variant: TierResult,
    consequence: TierResult,
    left_vcf: str | Any,
    right_vcf: str | Any,
    resource_plan: dict[str, Any] | None = None,
    timings: dict[str, float] | None = None,
    chromosome_summary: dict[str, Any] | None = None,
) -> None:
    payload: dict[str, Any] = {
        "preset": config.preset.name,
        "input_vcf": str(config.input_vcf),
        "annotated_left_vcf": str(left_vcf),
        "annotated_right_vcf": str(right_vcf),
        "sample_first_n": config.sample_first_n,
        "requested_chromosomes": config.chromosome_filter_raw,
        "effective_chromosomes": config.selected_chromosomes,
        "progress_log_path": str(artifacts.progress_log_path),
        "compare_mode": config.compare_mode,
        "bucket_count": config.compare_bucket_count,
        "compare_workers": config.compare_workers,
        "memory_budget_mb": config.memory_budget_mb,
        "resource_plan": resource_plan or {},
        "fingerprint_only": config.fingerprint_only,
        "timings": timings or {},
        "tiers": {
            "variant": {
                "equal": variant.equal,
                "left_only_rows": variant.left_only_rows,
                "right_only_rows": variant.right_only_rows,
                "unequal_rows": variant.unequal_rows,
                "joined_equal_rows": variant.joined_equal_rows,
                "diff_frame_path": str(variant.diff_frame_path),
                "mismatches_tsv_path": str(variant.mismatches_tsv_path),
                "details": variant.details,
                "per_chromosome": variant.per_chromosome,
            },
            "consequence": {
                "equal": consequence.equal,
                "left_only_rows": consequence.left_only_rows,
                "right_only_rows": consequence.right_only_rows,
                "unequal_rows": consequence.unequal_rows,
                "joined_equal_rows": consequence.joined_equal_rows,
                "diff_frame_path": str(consequence.diff_frame_path),
                "mismatches_tsv_path": str(consequence.mismatches_tsv_path),
                "details": consequence.details,
                "per_chromosome": consequence.per_chromosome,
            },
        },
        "chromosomes": chromosome_summary or {},
    }
    artifacts.summary_json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# Run Summary",
        "",
        f"- Preset: `{config.preset.name}`",
        f"- Input VCF: `{config.input_vcf}`",
        f"- Left annotated VCF: `{left_vcf}`",
        f"- Right annotated VCF: `{right_vcf}`",
        f"- Sample first N: `{config.sample_first_n if config.sample_first_n is not None else 'full'}`",
        f"- Requested chromosomes: `{config.chromosome_filter_raw or 'all'}`",
        f"- Effective chromosomes: `{','.join(config.selected_chromosomes) if config.selected_chromosomes else 'all'}`",
        f"- Progress log: `{artifacts.progress_log_path}`",
        f"- Compare mode: `{config.compare_mode}`",
        f"- Memory budget MB: `{config.memory_budget_mb if config.memory_budget_mb is not None else 'auto'}`",
        f"- Fingerprint only: `{config.fingerprint_only}`",
        "",
    ]
    if resource_plan:
        lines.extend(
            [
                "## Resource Plan",
                "",
                "```json",
                json.dumps(resource_plan, indent=2),
                "```",
                "",
            ]
        )
    if timings:
        lines.extend(
            [
                "## Timings",
                "",
                "```text",
                *(f"{name}={value:.3f}s" for name, value in timings.items()),
                "```",
                "",
            ]
        )
    lines.extend(
        [
            "## Variant Tier",
            "",
            "```text",
            variant.summary.rstrip(),
            "```",
            "",
            "## Consequence Tier",
            "",
            "```text",
            consequence.summary.rstrip(),
            "```",
            "",
        ]
    )
    if chromosome_summary and chromosome_summary.get("per_chromosome"):
        lines.extend(
            [
                "## Per Chromosome",
                "",
                "```json",
                json.dumps(chromosome_summary, indent=2),
                "```",
                "",
            ]
        )
    artifacts.summary_md_path.write_text("\n".join(lines), encoding="utf-8")
