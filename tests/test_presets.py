from vepyr_diffly.presets import get_preset, load_presets


def test_load_presets_contains_active_v1_preset() -> None:
    presets = load_presets()
    assert "ensembl_everything" in presets
    assert presets["ensembl_everything"].enabled is True


def test_get_preset_returns_expected_cache_flavor() -> None:
    preset = get_preset("ensembl_everything")
    assert preset.cache_flavor == "ensembl"
