from __future__ import annotations

import tomllib
from pathlib import Path

from .models import Preset


def load_presets(presets_path: Path | None = None) -> dict[str, Preset]:
    path = presets_path or Path(__file__).resolve().parents[2] / "presets" / "presets.toml"
    payload = tomllib.loads(path.read_text(encoding="utf-8"))
    raw_presets = payload.get("presets", {})
    presets: dict[str, Preset] = {}
    for name, raw in raw_presets.items():
        presets[name] = Preset(
            name=name,
            enabled=bool(raw["enabled"]),
            description=str(raw["description"]),
            species=str(raw["species"]),
            assembly=str(raw["assembly"]),
            cache_flavor=str(raw["cache_flavor"]),
            normalization_policy=str(raw["normalization_policy"]),
            supports_full_run=bool(raw["supports_full_run"]),
            supports_sampling=bool(raw["supports_sampling"]),
            vep_args=[str(item) for item in raw.get("vep_args", [])],
            vepyr_args=[str(item) for item in raw.get("vepyr_args", [])],
        )
    return presets


def get_preset(name: str, presets_path: Path | None = None) -> Preset:
    presets = load_presets(presets_path)
    try:
        return presets[name]
    except KeyError as exc:
        available = ", ".join(sorted(presets))
        raise KeyError(f"unknown preset '{name}'. Available presets: {available}") from exc
