#!/usr/bin/env python3
"""VEP-side checks run after the queue-scaling series (thesis items B19, D15, B9; 2026-10).

Subcommands (all use the VEP 116 image and options of the timing series: --cache --merged
--offline --everything --hgvs --fasta, --vcf, --no_stats):

  edge OUT_DIR       VEP on the two inputs of edge_contigs.py (unusual contigs, chrM vs MT);
                     records in/out, CSQ entries, consequences, VEP warnings (B19).
  forkprocs OUT_DIR  VEP on chr22 with --fork 1 and --fork 4; `docker top` sampled every
                     0.2 s, maximum number of concurrent vep processes (D15). No timing claim.
  b8 OUT_DIR         VEP on the edge-case input of edge_variants.py; per-record CSQ entries,
                     consequences and warnings next to the VEPYR summary (B8).
  sexmt OUT_DIR      ClinVar variants on X (random 5000 + all in PAR1/PAR2), Y (all) and MT
                     (all), contigs written as chrX, chrY and MT; VEP and VEPYR (workers=1)
                     annotate the same input and the E2E comparator of the vepyr repository
                     compares them (B9).

Environment: PYN (.venv-native python with vepyr 0.7.0), E2E_DIR (e2e-testing/scripts).
"""
from __future__ import annotations

import gzip
import json
import os
import random
import re
import subprocess
import sys
import time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
PYN = os.environ.get("PYN", str(D / ".venv-native/bin/python3"))
E2E_DIR = Path(os.environ.get("E2E_DIR", str(D / "tools/vepyr/e2e-testing/scripts")))
VEP_CACHE = D / "homo_sapiens_merged/116_GRCh38"
VEPYR_CACHE = D / "cache/116_GRCh38_merged"
FASTA = D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
IMAGE = "ensemblorg/ensembl-vep:release_116.0"
CLINVAR = D / "plugin_input/clinvar/clinvar.vcf.gz"
CHR22 = D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz"
EDGE = Path(__file__).resolve().parent.parent / "results/edge-contigs-2026-10-29"
EDGEV = Path(__file__).resolve().parent.parent / "results/edge-variants-2026-10-29"
PAR = {"X": [(10001, 2781479), (155701383, 156030895)]}
SEED = 20261030


def vep_cmd(inp: Path, out_dir: Path, out_name: str, name: str | None = None, fork: int | None = None) -> list[str]:
    cmd = ["docker", "run", "--rm", "--user", "501:20", "--env", "HOME=/tmp"]
    if name:
        cmd += ["--name", name]
    cmd += ["-v", f"{VEP_CACHE}:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
            "-v", f"{inp}:/input/{inp.name}:ro", "-v", f"{FASTA.parent}:/refdir:ro", "-v", f"{out_dir}:/output",
            IMAGE, "vep", "--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
            "--input_file", f"/input/{inp.name}", "--output_file", f"/output/{out_name}",
            "--vcf", "--force_overwrite", "--no_stats", "--everything", "--hgvs",
            "--fasta", f"/refdir/{FASTA.name}"]
    if fork:
        cmd += ["--fork", str(fork)]
    return cmd


def run_vepyr(inp: Path, out: Path) -> subprocess.CompletedProcess:
    child = ("import sys, vepyr; vepyr.annotate(vcf=sys.argv[1], cache_dir=sys.argv[3], everything=True, hgvs=True, "
             "reference_fasta=sys.argv[4], workers=1, output_vcf=sys.argv[2], show_progress=False)")
    return subprocess.run([PYN, "-c", child, str(inp), str(out), str(VEPYR_CACHE), str(FASTA)],
                          capture_output=True, text=True)


def records(path: Path) -> dict:
    op = gzip.open if path.suffix == ".gz" else open
    out = {}
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            csq = next((x[4:] for x in c[7].split(";") if x.startswith("CSQ=")), "")
            ents = [e.split("|") for e in csq.split(",")] if csq else []
            out[f"{c[0]}:{c[1]}:{c[3]}:{c[4]}"] = {"csq_entries": len(ents),
                                                   "consequences": sorted({t for e in ents for t in e[1].split("&")})}
    return out


def cmd_edge(out: Path) -> dict:
    res = {}
    for name in ("edge", "mt"):
        inp = EDGE / f"{name}.vcf"
        p = subprocess.run(vep_cmd(inp, out, f"{name}.vep.vcf"), capture_output=True, text=True)
        (out / f"{name}.vep.stderr.txt").write_text(p.stdout + p.stderr)
        warn = out / f"{name}.vep.vcf_warnings.txt"
        got = records(out / f"{name}.vep.vcf") if p.returncode == 0 else {}
        res[name] = {"returncode": p.returncode, "records_in": len(records(inp)), "records_out": len(got),
                     "records": got, "warnings": warn.read_text().splitlines() if warn.exists() else []}
    return res


def cmd_b8(out: Path) -> dict:
    inp = EDGEV / "edge_variants.vcf"
    p = subprocess.run(vep_cmd(inp, out, "edge_variants.vep.vcf"), capture_output=True, text=True)
    (out / "edge_variants.vep.stderr.txt").write_text(p.stdout + p.stderr)
    warn = out / "edge_variants.vep.vcf_warnings.txt"
    vep = records(out / "edge_variants.vep.vcf") if p.returncode == 0 else {}
    vepyr = records(EDGEV / "vepyr_w1.vcf")
    cases = json.loads((EDGEV / "cases.json").read_text())
    rows = []
    for c in cases:
        k = f"chr22:{c['pos']}:{c['ref']}:{c['alt']}"
        rows.append({"case": c["case"], "key": k, "vep": vep.get(k, "missing from output"),
                     "vepyr": vepyr.get(k, "missing from output")})
    return {"returncode": p.returncode, "records_vep": len(vep), "records_vepyr": len(vepyr), "cases": rows,
            "warnings": warn.read_text().splitlines() if warn.exists() else []}


def cmd_forkprocs(out: Path) -> dict:
    res = {}
    for fork in (1, 4):
        name = f"vep_forkprocs_{fork}"
        proc = subprocess.Popen(vep_cmd(CHR22, out, f"chr22_fork{fork}.vcf", name=name, fork=fork),
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        samples = []
        while proc.poll() is None:
            top = subprocess.run(["docker", "top", name, "-o", "pid,ppid,args"], capture_output=True, text=True)
            if top.returncode == 0:
                rows = [ln.split(None, 2) for ln in top.stdout.splitlines()[1:]]
                vep = [r for r in rows if len(r) == 3 and re.search(r"(^|/|\s)vep(\s|$)", r[2])]
                pids = {r[0] for r in vep}
                children = [r for r in vep if r[1] in pids]
                samples.append({"t": time.time(), "vep_processes": len(vep), "vep_children": len(children)})
            time.sleep(0.2)
        log = proc.stdout.read()
        (out / f"fork{fork}.log").write_text(log)
        (out / f"fork{fork}_samples.jsonl").write_text("".join(json.dumps(s) + "\n" for s in samples))
        (out / f"chr22_fork{fork}.vcf").unlink(missing_ok=True)
        res[f"fork{fork}"] = {"returncode": proc.returncode, "samples": len(samples),
                              "max_vep_processes": max((s["vep_processes"] for s in samples), default=None),
                              "max_vep_children": max((s["vep_children"] for s in samples), default=None)}
    return res


def build_sexmt_input(out: Path) -> tuple[Path, dict, dict]:
    rng = random.Random(SEED)
    rows = {"X": [], "Y": [], "MT": []}
    header = []
    with gzip.open(CLINVAR, "rt") as f:
        for line in f:
            if line.startswith("##"):
                if line.startswith(("##fileformat", "##reference")):
                    header.append(line)
                continue
            if line.startswith("#"):
                continue
            c = line.split("\t", 6)
            if c[0] in rows and c[4] != "." and re.fullmatch(r"[ACGT]+", c[3]) and re.fullmatch(r"[ACGT,]+", c[4]):
                rows[c[0]].append((int(c[1]), c[3], c[4]))
    par = {r for r in rows["X"] if any(a <= r[0] <= b for a, b in PAR["X"])}
    rest = sorted({r for r in rows["X"]} - par)
    pick_x = sorted(par | set(rng.sample(rest, 5000)))
    sel = {"chrX": pick_x, "chrY": sorted(set(rows["Y"])), "MT": sorted(set(rows["MT"]))}
    inp = out / "clinvar_sex_mt.vcf"
    with inp.open("w") as f:
        f.write("##fileformat=VCFv4.2\n" + "".join(f"##contig=<ID={k}>\n" for k in sel))
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, rs in sel.items():
            for pos, ref, alt in rs:
                f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")
    counts = {k: len(v) for k, v in sel.items()}
    counts["chrX_PAR"] = len(par)
    return inp, sel, counts


def cmd_sexmt(out: Path) -> dict:
    inp, sel, counts = build_sexmt_input(out)
    t0 = time.time()
    pv = subprocess.run(vep_cmd(inp, out, "vep.vcf"), capture_output=True, text=True)
    t1 = time.time()
    py = run_vepyr(inp, out / "vepyr.vcf")
    t2 = time.time()
    (out / "vep.stderr.txt").write_text(pv.stdout + pv.stderr)
    (out / "vepyr.stderr.txt").write_text(py.stderr)
    res = {"input_records": counts, "vep_returncode": pv.returncode, "vepyr_returncode": py.returncode,
           "seconds_vep": t1 - t0, "seconds_vepyr": t2 - t1}
    if pv.returncode == 0 and py.returncode == 0:
        rv, ry = records(out / "vep.vcf"), records(out / "vepyr.vcf")
        res["records_vep"], res["records_vepyr"] = len(rv), len(ry)
        res["per_contig"] = {ch: {"vep": sum(k.startswith(ch + ":") for k in rv),
                                  "vepyr": sum(k.startswith(ch + ":") for k in ry)} for ch in sel}
        sys.path.insert(0, str(E2E_DIR))
        from comparison import compare  # noqa: E402
        import contextlib, io
        with contextlib.redirect_stdout(io.StringIO()):
            res["e2e"] = compare.compare_vcfs(str(out / "vepyr.vcf"), str(out / "vep.vcf"), "sex-mt",
                                              mismatch_ledger_path=str(out / "e2e_ledger.jsonl"))
    return res


def main() -> int:
    sub, out = sys.argv[1], Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    res = {"edge": cmd_edge, "b8": cmd_b8, "forkprocs": cmd_forkprocs, "sexmt": cmd_sexmt}[sub](out)
    (out / "summary.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps(res, indent=1, default=str)[:4000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
