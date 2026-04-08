from __future__ import annotations

import importlib.util
from pathlib import Path
import pytest


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "build_chr_cache.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("build_chr_cache", SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_resolve_requested_chromosomes_accepts_csv_and_positional() -> None:
    module = _load_module()

    result = module.resolve_requested_chromosomes(["chr1", "Y"], "chrY,2")

    assert result == ["1", "Y", "2"]


def test_resolve_local_cache_source_matches_runtime_layout() -> None:
    module = _load_module()

    result = module.resolve_local_cache_source(
        vep_cache_dir=Path("/tmp/.vep"),
        species="homo_sapiens",
        release=115,
        assembly="GRCh38",
        cache_flavor="vep",
    )

    assert result == Path("/tmp/.vep/homo_sapiens/115_GRCh38")


def test_parse_plugins_deduplicates_and_normalizes_case() -> None:
    module = _load_module()

    result = module.parse_plugins("ClinVar, spliceai,clinvar , CADD")

    assert result == ["clinvar", "spliceai", "cadd"]


def test_resolve_plugin_source_returns_resolved_file(tmp_path: Path) -> None:
    module = _load_module()
    source = tmp_path / "clinvar.vcf.gz"
    source.write_text("fixture", encoding="utf-8")

    result = module.resolve_plugin_source("clinvar", source)

    assert result == source.resolve()


def test_resolve_plugin_source_requires_configured_path() -> None:
    module = _load_module()

    with pytest.raises(SystemExit, match="missing local source for plugin 'dbnsfp'"):
        module.resolve_plugin_source("dbnsfp", None)


def test_resolve_plugin_source_requires_existing_file(tmp_path: Path) -> None:
    module = _load_module()

    with pytest.raises(SystemExit, match="source file not found"):
        module.resolve_plugin_source("spliceai", tmp_path / "missing.vcf.gz")


def test_resolve_plugin_sources_requires_both_cadd_inputs(tmp_path: Path) -> None:
    module = _load_module()
    snv = tmp_path / "whole_genome_SNVs.tsv.gz"
    snv.write_text("fixture", encoding="utf-8")

    args = module.argparse.Namespace(
        clinvar_source=None,
        spliceai_source=None,
        cadd_snv_source=snv,
        cadd_indel_source=None,
        alphamissense_source=None,
        dbnsfp_source=None,
    )

    with pytest.raises(SystemExit, match="missing local source for plugin 'cadd_indel'"):
        module.resolve_plugin_sources(args, ["cadd"])


def test_resolve_plugin_sources_returns_split_cadd_mapping(tmp_path: Path) -> None:
    module = _load_module()
    snv = tmp_path / "whole_genome_SNVs.tsv.gz"
    indel = tmp_path / "gnomad.genomes.r4.0.indel.tsv.gz"
    snv.write_text("fixture", encoding="utf-8")
    indel.write_text("fixture", encoding="utf-8")

    args = module.argparse.Namespace(
        clinvar_source=None,
        spliceai_source=None,
        cadd_snv_source=snv,
        cadd_indel_source=indel,
        alphamissense_source=None,
        dbnsfp_source=None,
    )

    result = module.resolve_plugin_sources(args, ["cadd"])

    assert result == {
        "cadd": {
            "snv": snv.resolve(),
            "indel": indel.resolve(),
        }
    }
