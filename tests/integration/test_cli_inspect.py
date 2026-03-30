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


def test_analyze_consequence_mismatches_prints_summary_and_writes_json(
    tmp_path: Path, capsys: object
) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "consequence_mismatches.tsv").write_text(
        "\t".join(
            [
                "chrom",
                "pos",
                "ref",
                "alt",
                "Feature",
                "Consequence",
                "IMPACT",
                "HGNC_ID",
                "HGVSc",
                "HGVSp",
                "diff_kind",
                "duplicate_count_left",
                "duplicate_count_right",
            ]
        )
        + "\n"
        + "\n".join(
            [
                "chr1\t10\tA\tT\tENST1\tmissense_variant\tMODERATE\tHGNC:1\tc.1A>T\tp.Lys1Asn\tleft_only\t1\t",
                "chr1\t10\tA\tT\tENST1\tmissense_variant\tMODERATE\tHGNC:1\tc.1A>T\t.\tright_only\t\t1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_json = tmp_path / "analysis.json"

    exit_code = main(
        [
            "analyze-consequence-mismatches",
            "--run-dir",
            str(run_dir),
            "--output-json",
            str(output_json),
        ]
    )

    assert exit_code == 0
    captured = capsys.readouterr()
    assert '"paired_keys": 1' in captured.out
    assert output_json.exists()
    payload = json.loads(output_json.read_text(encoding="utf-8"))
    assert payload["paired_keys"] == 1


def test_extract_mismatch_csq_examples_prints_and_writes_json(
    tmp_path: Path, capsys: object
) -> None:
    run_dir = tmp_path / "run"
    runtime_dir = run_dir / "runtime"
    runtime_dir.mkdir(parents=True)
    analysis = {
        "category_examples": {
            "missing_hgvs": [
                {
                    "key": {
                        "chrom": "chr1",
                        "pos": "10",
                        "ref": "A",
                        "alt": "T",
                        "Feature": "ENST1",
                    },
                    "diff_fields": ["HGVSc"],
                    "categories": ["missing_hgvs"],
                }
            ]
        }
    }
    analysis_path = run_dir / "consequence_mismatch_analysis.json"
    analysis_path.write_text(json.dumps(analysis), encoding="utf-8")
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    (runtime_dir / "vep.annotated.vcf").write_text(
        header
        + "chr1\t10\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|GENE|ENSG1|Transcript|ENST1|protein_coding|||ENST1:c.1A>T|ENSP1:p.Lys1Asn\n",
        encoding="utf-8",
    )
    (runtime_dir / "vepyr.annotated.vcf").write_text(
        header
        + "chr1\t10\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|GENE|ENSG1|Transcript|ENST1|protein_coding||||ENSP1:p.Lys1Asn\n",
        encoding="utf-8",
    )
    output_json = tmp_path / "raw-csq.json"

    exit_code = main(
        [
            "extract-mismatch-csq-examples",
            "--run-dir",
            str(run_dir),
            "--output-json",
            str(output_json),
        ]
    )

    assert exit_code == 0
    captured = capsys.readouterr()
    assert '"missing_hgvs"' in captured.out
    payload = json.loads(output_json.read_text(encoding="utf-8"))
    example = payload["categories"]["missing_hgvs"][0]
    assert "ENST1:c.1A>T" in example["left_raw_csq"][0]
    assert "ENST1:c.1A>T" not in example["right_raw_csq"][0]
