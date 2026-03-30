from __future__ import annotations

from pathlib import Path
import types

from vepyr_diffly import vepyr_runner


def test_import_vepyr_prefers_current_environment(monkeypatch) -> None:
    current_module = types.SimpleNamespace(annotate=lambda **_: None)

    def fake_import_module(name: str):
        assert name == "vepyr"
        return current_module

    monkeypatch.setattr(vepyr_runner.importlib, "import_module", fake_import_module)
    monkeypatch.setattr(vepyr_runner, "_repo_vepyr_src", lambda: Path("/nonexistent"))

    imported = vepyr_runner._import_vepyr()

    assert imported is current_module


def test_import_vepyr_falls_back_to_repo_checkout(monkeypatch, tmp_path: Path) -> None:
    current_module_error = ModuleNotFoundError("vepyr._core")
    repo_module = types.SimpleNamespace(annotate=lambda **_: None)
    calls = {"count": 0}

    def fake_import_module(name: str):
        assert name == "vepyr"
        calls["count"] += 1
        if calls["count"] == 1:
            raise current_module_error
        return repo_module

    repo_src = tmp_path / "vepyr" / "src"
    repo_src.mkdir(parents=True)
    monkeypatch.setattr(vepyr_runner.importlib, "import_module", fake_import_module)
    monkeypatch.setattr(vepyr_runner, "_repo_vepyr_src", lambda: repo_src)
    monkeypatch.delitem(vepyr_runner.sys.modules, "vepyr", raising=False)

    imported = vepyr_runner._import_vepyr()

    assert imported is repo_module
    assert calls["count"] == 2
