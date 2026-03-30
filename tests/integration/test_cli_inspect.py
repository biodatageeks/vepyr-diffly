import json
from pathlib import Path

from vepyr_diffly.cli import main


def test_inspect_run_prints_existing_summary(tmp_path: Path, capsys: object) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "summary.json").write_text(
        json.dumps({"preset": "ensembl_everything"}), encoding="utf-8"
    )

    exit_code = main(["inspect-run", "--run-dir", str(run_dir)])

    assert exit_code == 0
    captured = capsys.readouterr()
    assert "ensembl_everything" in captured.out
