from __future__ import annotations

import gzip
import importlib.util
import subprocess
import sys
from pathlib import Path


def _load_module():
    module_path = (
        Path(__file__).resolve().parents[1] / "scripts" / "create_plugin_indexes.py"
    )
    spec = importlib.util.spec_from_file_location("create_plugin_indexes", module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    previous = sys.modules.get("create_plugin_indexes")
    sys.modules["create_plugin_indexes"] = module
    try:
        spec.loader.exec_module(module)
    finally:
        if previous is None:
            sys.modules.pop("create_plugin_indexes", None)
        else:
            sys.modules["create_plugin_indexes"] = previous
    return module


def _write_bgzf_stub(path: Path) -> None:
    # Minimal empty BGZF block.
    path.write_bytes(
        bytes.fromhex("1f8b08040000000000ff0600424302001b00030000000000000000")
    )


def test_is_bgzf_distinguishes_bgzf_from_plain_gzip(tmp_path: Path) -> None:
    module = _load_module()
    bgzf_path = tmp_path / "bgzf.tsv.gz"
    gzip_path = tmp_path / "plain.tsv.gz"

    _write_bgzf_stub(bgzf_path)
    with gzip.open(gzip_path, "wt", encoding="utf-8") as handle:
        handle.write("#Chrom\tPos\n1\t10\n")

    assert module.is_bgzf(bgzf_path) is True
    assert module.is_bgzf(gzip_path) is False


def test_resolve_inputs_expands_cadd_to_two_files(tmp_path: Path) -> None:
    module = _load_module()
    (tmp_path / "whole_genome_SNVs.tsv.gz").write_text("fixture", encoding="utf-8")
    (tmp_path / "gnomad.genomes.r4.0.indel.tsv.gz").write_text("fixture", encoding="utf-8")

    resolved = module.resolve_inputs(tmp_path, ["cadd"])

    assert [item[1].label for item in resolved] == ["cadd_snv", "cadd_indel"]
    assert [item[2].name for item in resolved] == [
        "whole_genome_SNVs.tsv.gz",
        "gnomad.genomes.r4.0.indel.tsv.gz",
    ]


def test_bgzf_sibling_path_uses_bgz_extension() -> None:
    module = _load_module()
    assert module.bgzf_sibling_path(Path("clinvar.vcf.gz")) == Path("clinvar.vcf.bgz")
    assert module.bgzf_sibling_path(Path("AlphaMissense_hg38.tsv.gz")) == Path(
        "AlphaMissense_hg38.tsv.bgz"
    )


def test_tabix_command_builds_expected_args() -> None:
    module = _load_module()
    module.ensure_tool = lambda name: f"/usr/bin/{name}"

    assert module.tabix_command(Path("clinvar.vcf.gz"), "vcf", force=False) == [
        "/usr/bin/tabix",
        "-p",
        "vcf",
        "clinvar.vcf.gz",
    ]
    assert module.tabix_command(Path("dbNSFP5.3.1a_grch38.gz"), "tsv", force=True) == [
        "/usr/bin/tabix",
        "-f",
        "-s",
        "1",
        "-b",
        "2",
        "-e",
        "2",
        "-c",
        "#",
        "dbNSFP5.3.1a_grch38.gz",
    ]


def test_run_with_progress_invokes_subprocess(monkeypatch, capsys) -> None:
    module = _load_module()

    class FakeProc:
        def __init__(self) -> None:
            self._polls = 0
            self.stderr = type("Stderr", (), {"read": lambda self: ""})()

        def poll(self):
            self._polls += 1
            return 0 if self._polls >= 2 else None

    monkeypatch.setattr(module.subprocess, "Popen", lambda *args, **kwargs: FakeProc())
    monkeypatch.setattr(module.time, "sleep", lambda _: None)

    elapsed = module.run_with_progress(
        ["/usr/bin/tabix", "example.gz"],
        label="dbnsfp",
        action="tabix example.gz",
        progress_interval_seconds=0.0,
    )

    assert elapsed >= 0.0
    output = capsys.readouterr().out
    assert "[dbnsfp] tabix example.gz: still running" in output
    assert "[dbnsfp] tabix example.gz: done in" in output


def test_run_with_progress_raises_on_failure(monkeypatch) -> None:
    module = _load_module()

    class FakeProc:
        def __init__(self) -> None:
            self.stderr = type("Stderr", (), {"read": lambda self: "boom"})()

        def poll(self):
            return 1

    monkeypatch.setattr(module.subprocess, "Popen", lambda *args, **kwargs: FakeProc())

    try:
        module.run_with_progress(
            ["/usr/bin/tabix", "example.gz"],
            label="dbnsfp",
            action="tabix example.gz",
        )
    except subprocess.CalledProcessError as exc:
        assert exc.returncode == 1
        assert exc.stderr == "boom"
    else:
        raise AssertionError("expected CalledProcessError")
