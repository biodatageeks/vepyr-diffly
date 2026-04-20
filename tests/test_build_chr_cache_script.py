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


def test_resolve_requested_chromosomes_accepts_csv_only() -> None:
    module = _load_module()

    result = module.resolve_requested_chromosomes("chrY,2")

    assert result == ["Y", "2"]


def test_resolve_requested_chromosomes_returns_empty_for_missing_value() -> None:
    module = _load_module()

    result = module.resolve_requested_chromosomes(None)

    assert result == []


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


def test_resolve_partitioned_cache_dir_matches_runtime_layout() -> None:
    module = _load_module()

    result = module.resolve_partitioned_cache_dir(
        cache_dir=Path("/tmp/vepyr_cache"),
        release=115,
        assembly="GRCh38",
        cache_flavor="vep",
    )

    assert result == Path("/tmp/vepyr_cache/115_GRCh38_vep")


def test_parse_args_rejects_positional_chromosomes() -> None:
    module = _load_module()

    with pytest.raises(SystemExit):
        module.parse_args(["1"])


def test_resolve_legacy_partitioned_cache_dir_matches_previous_layout() -> None:
    module = _load_module()

    result = module.resolve_legacy_partitioned_cache_dir(
        cache_dir=Path("/tmp/vepyr_cache"),
        release=115,
        assembly="GRCh38",
        cache_flavor="vep",
    )

    assert result == Path("/tmp/vepyr_cache/parquet/115_GRCh38_vep")


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


def test_remove_old_layout_cache_removes_recognized_legacy_artifacts(tmp_path: Path) -> None:
    module = _load_module()
    legacy_partitioned = tmp_path / "parquet" / "115_GRCh38_vep"
    root_plugin_dir = tmp_path / "clinvar"
    root_plugin_fjall = tmp_path / "clinvar.fjall"
    variation_fjall = tmp_path / "variation.fjall"
    sift_fjall = tmp_path / "translation_sift.fjall"
    current_cache_dir = tmp_path / "115_GRCh38_vep"

    for path in (
        legacy_partitioned,
        root_plugin_dir,
        root_plugin_fjall,
        variation_fjall,
        sift_fjall,
        current_cache_dir / "clinvar",
        current_cache_dir / "clinvar.fjall",
    ):
        path.mkdir(parents=True, exist_ok=True)

    removed = module.remove_old_layout_cache(
        cache_dir=tmp_path,
        release=115,
        assembly="GRCh38",
        cache_flavor="vep",
        plugins=["clinvar"],
    )

    assert set(removed) == {
        legacy_partitioned,
        root_plugin_dir,
        root_plugin_fjall,
        variation_fjall,
        sift_fjall,
    }
    assert not legacy_partitioned.exists()
    assert not root_plugin_dir.exists()
    assert not root_plugin_fjall.exists()
    assert not variation_fjall.exists()
    assert not sift_fjall.exists()
    assert (current_cache_dir / "clinvar").exists()
    assert (current_cache_dir / "clinvar.fjall").exists()


def test_remove_current_plugin_outputs_removes_only_requested_plugin_artifacts(
    tmp_path: Path,
) -> None:
    module = _load_module()
    partitioned_cache_dir = tmp_path / "115_GRCh38_vep"
    clinvar_dir = partitioned_cache_dir / "clinvar"
    clinvar_fjall = partitioned_cache_dir / "clinvar.fjall"
    spliceai_dir = partitioned_cache_dir / "spliceai"
    variation_dir = partitioned_cache_dir / "variation"

    for path in (clinvar_dir, clinvar_fjall, spliceai_dir, variation_dir):
        path.mkdir(parents=True, exist_ok=True)

    removed = module.remove_current_plugin_outputs(
        partitioned_cache_dir=partitioned_cache_dir,
        plugins=["clinvar"],
    )

    assert set(removed) == {clinvar_dir, clinvar_fjall}
    assert not clinvar_dir.exists()
    assert not clinvar_fjall.exists()
    assert spliceai_dir.exists()
    assert variation_dir.exists()


def test_create_core_preview_cache_truncates_selected_chromosome(tmp_path: Path) -> None:
    module = _load_module()
    source_root = tmp_path / "source"
    chrom1 = source_root / "1"
    chrom2 = source_root / "2"
    chrom1.mkdir(parents=True)
    chrom2.mkdir(parents=True)
    (source_root / "info.txt").write_text("info", encoding="utf-8")
    (source_root / "chr_synonyms.txt").write_text("syn", encoding="utf-8")

    import gzip

    with gzip.open(chrom1 / "1-1000000.gz", "wt", encoding="utf-8") as handle:
        handle.write("a\nb\nc\n")
    with gzip.open(chrom1 / "1000001-2000000.gz", "wt", encoding="utf-8") as handle:
        handle.write("d\ne\nf\n")
    with gzip.open(chrom1 / "1-1000000_reg.gz", "wt", encoding="utf-8") as handle:
        handle.write("r1\nr2\nr3\n")
    with gzip.open(chrom1 / "1000001-2000000_reg.gz", "wt", encoding="utf-8") as handle:
        handle.write("r4\nr5\nr6\n")
    with gzip.open(chrom1 / "all_vars.gz", "wt", encoding="utf-8") as handle:
        handle.write("v1\nv2\nv3\n")
    with gzip.open(chrom2 / "1-1000000.gz", "wt", encoding="utf-8") as handle:
        handle.write("x1\nx2\n")

    preview_root = Path(module.create_core_preview_cache(source_root, 2, ["1"]))
    try:
        assert (preview_root / "info.txt").exists()
        assert (preview_root / "chr_synonyms.txt").exists()
        assert (preview_root / "1" / "1-1000000.gz").exists()
        assert not (preview_root / "1" / "1000001-2000000.gz").exists()
        assert not (preview_root / "2").exists()

        with gzip.open(preview_root / "1" / "1-1000000.gz", "rt", encoding="utf-8") as handle:
            assert handle.read().splitlines() == ["a", "b", "c"]
        with gzip.open(preview_root / "1" / "1-1000000_reg.gz", "rt", encoding="utf-8") as handle:
            assert handle.read().splitlines() == ["r1", "r2", "r3"]
        with gzip.open(preview_root / "1" / "all_vars.gz", "rt", encoding="utf-8") as handle:
            assert handle.read().splitlines() == ["v1", "v2", "v3"]
    finally:
        import shutil

        shutil.rmtree(preview_root, ignore_errors=True)


def test_create_preview_source_filters_requested_chromosomes_for_tsv(tmp_path: Path) -> None:
    module = _load_module()
    source = tmp_path / "dbnsfp.tsv.gz"

    import gzip

    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write("chr\tpos(1-based)\tref\talt\n")
        handle.write("1\t10\tA\tG\n")
        handle.write("Y\t20\tC\tT\n")
        handle.write("Y\t21\tG\tA\n")

    preview, _ = module._create_preview_source(source, "dbnsfp", 2, ["Y"])
    try:
        with gzip.open(preview, "rt", encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        assert lines == [
            "chr\tpos(1-based)\tref\talt",
            "chrY\t20\tC\tT",
            "chrY\t21\tG\tA",
        ]
    finally:
        preview.unlink(missing_ok=True)


def test_create_preview_source_skips_tsv_preamble_before_header(tmp_path: Path) -> None:
    module = _load_module()
    source = tmp_path / "cadd.tsv.gz"

    import gzip

    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write("##CADD metadata\n")
        handle.write("#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n")
        handle.write("1\t10\tA\tG\t0.1\t1.0\n")
        handle.write("1\t11\tC\tT\t0.2\t2.0\n")

    preview, _ = module._create_preview_source(source, "cadd_snv", 1, None)
    try:
        with gzip.open(preview, "rt", encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        assert lines == [
            "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED",
            "chr1\t10\tA\tG\t0.1\t1.0",
        ]
    finally:
        preview.unlink(missing_ok=True)


def test_create_preview_source_filters_requested_chromosomes_for_vcf(tmp_path: Path) -> None:
    module = _load_module()
    source = tmp_path / "clinvar.vcf.gz"

    import gzip

    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        handle.write("chr1\t10\t.\tA\tG\t.\t.\t.\n")
        handle.write("chrY\t20\t.\tC\tT\t.\t.\t.\n")

    preview, _ = module._create_preview_source(source, "clinvar", 1, ["Y"])
    try:
        with gzip.open(preview, "rt", encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        assert lines == [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chrY\t20\t.\tC\tT\t.\t.\t.",
        ]
    finally:
        preview.unlink(missing_ok=True)


def test_prepare_preview_source_prefers_tabix_for_full_chromosome_scoped_build(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _load_module()
    source = tmp_path / "dbNSFP5.3.1a_grch38.gz"
    source.write_text("fixture", encoding="utf-8")
    temp_paths: list[Path] = []
    tabix_preview = tmp_path / "dbnsfp.tabix-preview.gz"
    tabix_preview.write_text("preview", encoding="utf-8")

    def fake_extract(original, plugin_name, preview_rows, chromosomes):
        assert original == source
        assert plugin_name == "dbnsfp"
        assert preview_rows is None
        assert chromosomes == ["10"]
        return tabix_preview, 123

    def fail_create(*args, **kwargs):
        raise AssertionError("fallback path should not be used when tabix succeeds")

    monkeypatch.setattr(module, "_extract_tabix_preview_source", fake_extract)
    monkeypatch.setattr(module, "_create_preview_source", fail_create)

    preview, stats = module._prepare_preview_source(
        "dbnsfp",
        source,
        None,
        ["10"],
        temp_paths,
    )

    assert preview == tabix_preview
    assert temp_paths == [tabix_preview]
    assert stats == {
        "compressed_bytes": tabix_preview.stat().st_size,
        "uncompressed_bytes": 123,
    }


def test_build_plugin_caches_uses_builder_limit_for_global_preview_when_sorted_opt_in_enabled(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _load_module()
    source = tmp_path / "dbNSFP5.3.1a_grch38.gz"
    source.write_text("fixture", encoding="utf-8")
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    calls: list[tuple] = []

    class FakeVepyr:
        def build_plugin(
            self,
            plugin,
            source_path,
            cache_dir_value,
            *,
            build_fjall=True,
            partitions=8,
            chromosomes=None,
            assume_sorted_input=False,
            preview_rows=None,
        ):
            calls.append(
                (
                    plugin,
                    source_path,
                    cache_dir_value,
                    build_fjall,
                    partitions,
                    chromosomes,
                    assume_sorted_input,
                    preview_rows,
                )
            )

    def fail_prepare(*args, **kwargs):
        raise AssertionError("global preview should use builder row limit instead of temp preview")

    monkeypatch.setattr(module, "_prepare_preview_source", fail_prepare)
    monkeypatch.setattr(module, "_parquet_size", lambda _path: 0)

    module.build_plugin_caches(
        vepyr_module=FakeVepyr(),
        cache_dir=cache_dir,
        plugins=["dbnsfp"],
        plugin_sources={"dbnsfp": source},
        partitions=8,
        chromosomes=None,
        preview_rows=1000,
        assume_sorted_input=True,
    )

    assert calls == [
        (
            "dbnsfp",
            str(source),
            str(cache_dir),
            True,
            8,
            None,
            True,
            1000,
        )
    ]


def test_build_plugin_caches_falls_back_to_temp_preview_without_sorted_opt_in(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _load_module()
    source = tmp_path / "dbNSFP5.3.1a_grch38.gz"
    source.write_text("fixture", encoding="utf-8")
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    temp_preview = tmp_path / "dbnsfp-preview.gz"
    temp_preview.write_text("preview", encoding="utf-8")

    calls: list[tuple] = []

    class FakeVepyr:
        def build_plugin(
            self,
            plugin,
            source_path,
            cache_dir_value,
            *,
            build_fjall=True,
            partitions=8,
            chromosomes=None,
            assume_sorted_input=False,
            preview_rows=None,
        ):
            calls.append(
                (
                    plugin,
                    source_path,
                    cache_dir_value,
                    build_fjall,
                    partitions,
                    chromosomes,
                    assume_sorted_input,
                    preview_rows,
                )
            )

    monkeypatch.setattr(
        module,
        "_prepare_preview_source",
        lambda *args, **kwargs: (temp_preview, {"compressed_bytes": 1, "uncompressed_bytes": 2}),
    )
    monkeypatch.setattr(module, "_parquet_size", lambda _path: 0)

    module.build_plugin_caches(
        vepyr_module=FakeVepyr(),
        cache_dir=cache_dir,
        plugins=["dbnsfp"],
        plugin_sources={"dbnsfp": source},
        partitions=8,
        chromosomes=None,
        preview_rows=1000,
        assume_sorted_input=False,
    )

    assert calls == [
        (
            "dbnsfp",
            str(temp_preview),
            str(cache_dir),
            True,
            8,
            None,
            False,
            None,
        )
    ]


def test_build_plugin_caches_keeps_temp_preview_path_for_cadd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _load_module()
    snv = tmp_path / "whole_genome_SNVs.tsv.gz"
    indel = tmp_path / "gnomad.genomes.r4.0.indel.tsv.gz"
    snv.write_text("fixture", encoding="utf-8")
    indel.write_text("fixture", encoding="utf-8")
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    snv_preview = tmp_path / "cadd-snv-preview.gz"
    indel_preview = tmp_path / "cadd-indel-preview.gz"
    snv_preview.write_text("preview", encoding="utf-8")
    indel_preview.write_text("preview", encoding="utf-8")

    calls: list[tuple] = []

    class FakeVepyr:
        def build_plugin(
            self,
            plugin,
            source_path,
            cache_dir_value,
            *,
            build_fjall=True,
            partitions=8,
            chromosomes=None,
            assume_sorted_input=False,
            preview_rows=None,
        ):
            calls.append(
                (
                    plugin,
                    source_path,
                    cache_dir_value,
                    build_fjall,
                    partitions,
                    chromosomes,
                    assume_sorted_input,
                    preview_rows,
                )
            )

    def fake_prepare(plugin_name, source_path, preview_rows, chromosomes, temp_paths):
        if plugin_name == "cadd_snv":
            temp_paths.append(snv_preview)
            return snv_preview, {"compressed_bytes": 1, "uncompressed_bytes": 2}
        if plugin_name == "cadd_indel":
            temp_paths.append(indel_preview)
            return indel_preview, {"compressed_bytes": 3, "uncompressed_bytes": 4}
        raise AssertionError(f"unexpected plugin preview request: {plugin_name}")

    monkeypatch.setattr(module, "_prepare_preview_source", fake_prepare)
    monkeypatch.setattr(module, "_parquet_size", lambda _path: 0)

    module.build_plugin_caches(
        vepyr_module=FakeVepyr(),
        cache_dir=cache_dir,
        plugins=["cadd"],
        plugin_sources={"cadd": {"snv": snv, "indel": indel}},
        partitions=8,
        chromosomes=None,
        preview_rows=1000,
        assume_sorted_input=True,
    )

    assert calls == [
        (
            "cadd",
            {"snv": str(snv_preview), "indel": str(indel_preview)},
            str(cache_dir),
            True,
            8,
            None,
            True,
            None,
        )
    ]


def test_build_plugin_caches_can_skip_plugin_fjall(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    module = _load_module()
    source = tmp_path / "AlphaMissense_hg38.tsv.gz"
    source.write_text("fixture", encoding="utf-8")
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    calls: list[tuple] = []

    class FakeVepyr:
        def build_plugin(
            self,
            plugin,
            source_path,
            cache_dir_value,
            *,
            build_fjall=True,
            partitions=8,
            chromosomes=None,
            assume_sorted_input=False,
            preview_rows=None,
        ):
            calls.append(
                (
                    plugin,
                    source_path,
                    cache_dir_value,
                    build_fjall,
                    partitions,
                    chromosomes,
                    assume_sorted_input,
                    preview_rows,
                )
            )

    monkeypatch.setattr(module, "_prepare_preview_source", lambda *args, **kwargs: (source, None))
    monkeypatch.setattr(module, "_parquet_size", lambda _path: 0)

    module.build_plugin_caches(
        vepyr_module=FakeVepyr(),
        cache_dir=cache_dir,
        plugins=["alphamissense"],
        plugin_sources={"alphamissense": source},
        partitions=8,
        chromosomes=None,
        preview_rows=None,
        assume_sorted_input=False,
        build_fjall=False,
    )

    assert calls == [
        (
            "alphamissense",
            str(source),
            str(cache_dir),
            False,
            8,
            None,
            False,
            None,
        )
    ]
