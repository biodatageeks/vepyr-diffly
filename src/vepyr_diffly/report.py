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
) -> None:
    console.print(f"Preset: [bold]{config.preset.name}[/bold]")
    console.print(f"Input: {config.input_vcf}")
    console.print(f"Sample first N: {config.sample_first_n if config.sample_first_n is not None else 'full'}")
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


def write_run_summary(
    *,
    config: RuntimeConfig,
    artifacts: RunArtifacts,
    variant: TierResult,
    consequence: TierResult,
) -> None:
    payload: dict[str, Any] = {
        "preset": config.preset.name,
        "input_vcf": str(config.input_vcf),
        "sample_first_n": config.sample_first_n,
        "tiers": {
            "variant": {
                "equal": variant.equal,
                "left_only_rows": variant.left_only_rows,
                "right_only_rows": variant.right_only_rows,
                "unequal_rows": variant.unequal_rows,
                "joined_equal_rows": variant.joined_equal_rows,
                "diff_frame_path": str(variant.diff_frame_path),
                "mismatches_tsv_path": str(variant.mismatches_tsv_path),
            },
            "consequence": {
                "equal": consequence.equal,
                "left_only_rows": consequence.left_only_rows,
                "right_only_rows": consequence.right_only_rows,
                "unequal_rows": consequence.unequal_rows,
                "joined_equal_rows": consequence.joined_equal_rows,
                "diff_frame_path": str(consequence.diff_frame_path),
                "mismatches_tsv_path": str(consequence.mismatches_tsv_path),
            },
        },
    }
    artifacts.summary_json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    lines = [
        f"# Run Summary",
        "",
        f"- Preset: `{config.preset.name}`",
        f"- Input VCF: `{config.input_vcf}`",
        f"- Sample first N: `{config.sample_first_n if config.sample_first_n is not None else 'full'}`",
        "",
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
    artifacts.summary_md_path.write_text("\n".join(lines), encoding="utf-8")
