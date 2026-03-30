from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path


def test_compare_existing_command_reuses_annotated_vcfs(tmp_path: Path) -> None:
    fixture_dir = Path(__file__).resolve().parents[1] / "fixtures"
    run_dir = tmp_path / "compare-existing"

    completed = subprocess.run(
        [
            str(Path(__file__).resolve().parents[2] / ".venv" / "bin" / "python"),
            "-m",
            "vepyr_diffly.cli",
            "compare-existing",
            "--preset",
            "ensembl_everything",
            "--left-vcf",
            str(fixture_dir / "annotated_left.vcf"),
            "--right-vcf",
            str(fixture_dir / "annotated_right.vcf"),
            "--output-dir",
            str(run_dir),
            "--memory-budget-mb",
            "256",
        ],
        check=False,
        capture_output=True,
        text=True,
        env={
            **os.environ,
            "PYTHONPATH": str(Path(__file__).resolve().parents[2] / "src"),
        },
    )

    assert completed.returncode == 0, completed.stderr
    summary_path = run_dir / "summary.json"
    assert summary_path.exists()
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["annotated_left_vcf"].endswith("annotated_left.vcf")
    assert payload["annotated_right_vcf"].endswith("annotated_right.vcf")
    assert payload["resource_plan"]["memory_budget_mb"] == 256
    assert (run_dir / "runtime" / "compare.progress.log").exists()
