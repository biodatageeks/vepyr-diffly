#!/usr/bin/env python3
"""ClinVar and CADD values from echtvar and vcfanno against the VEPYR plugins (thesis item C9).

echtvar and vcfanno attach variant-level values only (no matching by transcript), so the
comparison covers ClinVar (CLNSIG, CLNREVSTAT, CLNDN, CLNVC, CLNVI) and CADD (raw score,
PHRED) on the HG002 chr22 or chr1 baseline:
  sources  the same files as the VEPYR plugin caches: ClinVar 20260706 (contig renamed to
           chrN) and CADD v1.7 SNV + gnomAD r4.0 indel TSVs, cut by tabix and written as one
           VCF (INFO CADD_RAW, CADD_PHRED; vcfanno matches TSV files by position only, VCF
           files by REF/ALT);
  vcfanno  v0.3.9, native osx_aarch64 binary, ops "self", -p 1 and -p 10;
  echtvar  v0.2.4, aarch64 musl binary in a debian container in the Colima VM; encode timed
           separately; floats stored with multiplier 1e6 (its integer encoding);
  VEPYR    0.7.0, plugins cadd and clinvar from the plugin cache of the plugin series,
           workers=1, no --everything (fewest other fields).
Per variant and field: equal (numbers within 1e-6), different, only in VEPYR, only in the
tool, missing in both. Wall time and peak RSS (/usr/bin/time -l, natively; cgroup
memory.peak in the container) of each run; one run each.

Usage: c9_echtvar_vcfanno.py OUT_DIR CHROM   (CHROM = chr22 or chr1; Colima running)
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
PYN = D / ".venv-native/bin/python3"
PI = D / "plugin_input"
INPUT = {"chr22": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz",
         "chr1": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr1.vcf.gz"}
CACHE = D / "cache/116_GRCh38_merged"
PLUGIN_ROOT = D / "cache/plugin_cache_116"
VCFANNO_URL = "https://github.com/brentp/vcfanno/releases/download/v0.3.9/vcfanno_osx_aarch64"
ECHTVAR_URL = "https://github.com/brentp/echtvar/releases/download/v0.2.4/echtvar-v0.2.4-aarch64-unknown-linux-musl"
CTR = "debian:bookworm-slim"
CLN = ["CLNSIG", "CLNREVSTAT", "CLNDN", "CLNVC", "CLNVI"]
CADD = ["CADD_RAW", "CADD_PHRED"]
VEPYR_FIELD = {**{f: f"ClinVar_{f}" for f in CLN}, "CADD_RAW": "CADD_RAW", "CADD_PHRED": "CADD_PHRED"}
VEPYR = r'''
import json, resource, sys, time, vepyr
t0 = time.perf_counter()
vepyr.annotate(vcf=sys.argv[1], cache_dir=sys.argv[2], plugin_cache_root=sys.argv[3], plugins=["cadd", "clinvar"],
               workers=1, output_vcf=sys.argv[4], show_progress=False)
print(json.dumps({"seconds": time.perf_counter() - t0, "max_rss_gib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3}))
'''


def sh(cmd: str) -> None:
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def fetch(url: str, path: Path) -> str:
    if not path.exists():
        path.write_bytes(urllib.request.urlopen(url).read())
        path.chmod(0o755)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources(src: Path, chrom: str) -> tuple[Path, Path]:
    n = chrom.removeprefix("chr")
    cln, cadd = src / f"clinvar_{chrom}.vcf.gz", src / f"cadd_{chrom}.vcf.gz"
    if not cln.exists():
        (src / "rename.txt").write_text(f"{n}\t{chrom}\n")
        sh(f"bcftools view -r {n} {PI}/clinvar/clinvar.vcf.gz | bcftools annotate --rename-chrs {src}/rename.txt "
           f"-Oz -o {cln} && tabix -p vcf {cln}")
    if not cadd.exists():
        hdr = ("##fileformat=VCFv4.2\n" f"##contig=<ID={chrom}>\n"
               '##INFO=<ID=CADD_RAW,Number=1,Type=Float,Description="CADD v1.7 RawScore">\n'
               '##INFO=<ID=CADD_PHRED,Number=1,Type=Float,Description="CADD v1.7 PHRED">\n'
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        (src / "cadd_header.txt").write_text(hdr)
        awk = f"awk -F'\\t' -v OFS='\\t' '{{print \"{chrom}\",$2,\".\",$3,$4,\".\",\".\",\"CADD_RAW=\"$5\";CADD_PHRED=\"$6}}'"
        sh(f"tabix {PI}/cadd/whole_genome_SNVs.tsv.gz {n} | {awk} > {src}/snv.txt")
        sh(f"tabix {PI}/cadd/gnomad.genomes.r4.0.indel.tsv.gz {n} | {awk} > {src}/indel.txt")
        sh(f"(cat {src}/cadd_header.txt; sort -m -k2,2n -s {src}/snv.txt {src}/indel.txt) | bgzip -c > {cadd} "
           f"&& tabix -p vcf {cadd} && rm {src}/snv.txt {src}/indel.txt")
    return cln, cadd


def timed(cmd: list[str], stdout=None) -> tuple[float, float | None, subprocess.CompletedProcess]:
    t = Path(os.environ.get("TMPDIR", "/tmp")) / f"c9_time_{os.getpid()}.txt"
    t0 = time.perf_counter()
    p = subprocess.run(["/usr/bin/time", "-l", "-o", str(t), *cmd], stdout=stdout, stderr=subprocess.PIPE, text=True)
    dt = time.perf_counter() - t0
    m = re.search(r"(\d+)\s+maximum resident set size", t.read_text()) if t.exists() else None
    return dt, int(m.group(1)) / 1024**3 if m else None, p


def info_values(path: Path, prefix: str, fields: list[str]) -> dict:
    out = {}
    with (gzip.open if path.suffix == ".gz" else open)(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 8)
            kv = dict(x.split("=", 1) if "=" in x else (x, "") for x in c[7].split(";"))
            out[f"{c[0]}:{c[1]}:{c[3]}:{c[4]}"] = {fl: kv.get(prefix + fl, "") for fl in fields}
    return out


def vepyr_values(path: Path, fields: list[str]) -> tuple[dict, int]:
    out, inconsistent, idx = {}, 0, None
    with open(path) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                names = line.split("Format: ")[1].rstrip('">\n').split("|")
                idx = {fl: names.index(VEPYR_FIELD[fl]) for fl in fields}
            if line[0] == "#":
                continue
            c = line.split("\t", 8)
            csq = next((x[4:] for x in c[7].split(";") if x.startswith("CSQ=")), "")
            ents = [e.split("|") for e in csq.split(",")] if csq else []
            row = {}
            for fl, i in idx.items():
                vals = {e[i] for e in ents if len(e) > i and e[i]}
                inconsistent += len(vals) > 1
                row[fl] = sorted(vals)[0] if vals else ""
            out[f"{c[0]}:{c[1]}:{c[3]}:{c[4]}"] = row
    return out, inconsistent


def same(a: str, b: str) -> bool:
    if a == b:
        return True
    try:
        return abs(float(a) - float(b)) <= 1e-6 * max(1.0, abs(float(a)))
    except ValueError:  # list separators: VEP CSQ uses '&', VCF INFO ',' (or %2C escaped)
        norm = lambda s: s.replace("%2C", ",").replace("&", ",")  # noqa: E731
        return norm(a) == norm(b)


def compare(ref: dict, tool: dict, fields: list[str]) -> dict:
    res = {}
    for fl in fields:
        c = dict.fromkeys(("equal", "different", "only_vepyr", "only_tool", "none"), 0)
        ex = []
        for k, r in ref.items():
            a, b = r.get(fl, ""), tool.get(k, {}).get(fl, "")
            a, b = ("" if a == "." else a), ("" if b == "." else b)
            key = "none" if not a and not b else "only_vepyr" if not b else "only_tool" if not a else \
                "equal" if same(a, b) else "different"
            c[key] += 1
            if key in ("different", "only_vepyr", "only_tool") and len(ex) < 5:
                ex.append({"variant": k, "vepyr": a, "tool": b})
        res[fl] = {**c, "examples": ex}
    return res


def main() -> int:
    out, chrom = Path(sys.argv[1]).resolve(), sys.argv[2]
    src, bins = out / "src", out / "bin"
    src.mkdir(parents=True, exist_ok=True); bins.mkdir(exist_ok=True)
    res = {"chrom": chrom, "vcfanno_sha256": fetch(VCFANNO_URL, bins / "vcfanno"),
           "echtvar_sha256": fetch(ECHTVAR_URL, bins / "echtvar")}
    inp = INPUT[chrom]
    t0 = time.time()
    cln, cadd = sources(src, chrom)
    res["source_prep_s"] = time.time() - t0
    fields = CLN + CADD

    v = out / f"{chrom}_vepyr.vcf"
    p = subprocess.run([str(PYN), "-c", VEPYR, str(inp), str(CACHE), str(PLUGIN_ROOT), str(v)],
                       capture_output=True, text=True)
    if p.returncode:
        raise SystemExit(p.stderr[-2000:])
    res["vepyr"] = json.loads(p.stdout.strip().splitlines()[-1])
    ref, res["vepyr"]["csq_entries_with_differing_values"] = vepyr_values(v, fields)

    (out / "vcfanno.toml").write_text(
        f'[[annotation]]\nfile="{cln}"\nfields={json.dumps(CLN)}\nnames={json.dumps(["va_" + f for f in CLN])}\n'
        f'ops={json.dumps(["self"] * len(CLN))}\n\n[[annotation]]\nfile="{cadd}"\nfields={json.dumps(CADD)}\n'
        f'names={json.dumps(["va_" + f for f in CADD])}\nops={json.dumps(["self"] * len(CADD))}\n')
    for procs in (1, 10):
        o = out / f"{chrom}_vcfanno_p{procs}.vcf"
        with open(o, "w") as fo:
            dt, rss, p = timed([str(bins / "vcfanno"), "-p", str(procs), str(out / "vcfanno.toml"), str(inp)], stdout=fo)
        (out / f"vcfanno_p{procs}.log").write_text(p.stderr)
        res[f"vcfanno_p{procs}"] = {"seconds": dt, "max_rss_gib": rss, "returncode": p.returncode,
                                    "compare": compare(ref, info_values(o, "va_", fields), fields)}
        o.unlink()

    (out / "echt_clinvar.json").write_text(json.dumps(
        [{"field": f, "alias": f"echt_{f}", "missing_string": "."} for f in CLN]))
    (out / "echt_cadd.json").write_text(json.dumps(
        [{"field": f, "alias": f"echt_{f}", "multiplier": 1000000, "zigzag": True, "missing_value": -2147483648}
         for f in CADD]))
    run = ["docker", "run", "--rm", "--user", "501:20", "-v", f"{out}:/w", "-v", f"{src}:/src:ro",
           "-v", f"{inp.parent}:/in:ro", CTR, "sh", "-c"]
    wrap = '; rc=$?; echo "memory_peak $(cat /sys/fs/cgroup/memory.peak)"; exit $rc'
    for name, cmd in (("encode_clinvar", f"/w/bin/echtvar encode /w/clinvar.echtvar.zip /w/echt_clinvar.json /src/{cln.name}"),
                      ("encode_cadd", f"/w/bin/echtvar encode /w/cadd.echtvar.zip /w/echt_cadd.json /src/{cadd.name}"),
                      ("anno", f"/w/bin/echtvar anno -e /w/clinvar.echtvar.zip -e /w/cadd.echtvar.zip /in/{inp.name} "
                               f"/w/{chrom}_echtvar.vcf.gz")):
        t1 = time.perf_counter()
        p = subprocess.run(run + [cmd + wrap], capture_output=True, text=True)
        m = re.search(r"memory_peak (\d+)", p.stdout)
        res[f"echtvar_{name}"] = {"seconds": time.perf_counter() - t1, "returncode": p.returncode,
                                  "cgroup_memory_peak_gib": int(m.group(1)) / 1024**3 if m else None,
                                  "stderr_tail": p.stderr[-800:]}
    e = out / f"{chrom}_echtvar.vcf.gz"
    if e.exists():
        res["echtvar_anno"]["compare"] = compare(ref, info_values(e, "echt_", fields), fields)
        e.unlink()
    v.unlink()
    (out / f"{chrom}_summary.json").write_text(json.dumps(res, indent=1) + "\n")
    for k in ("vcfanno_p1", "vcfanno_p10", "echtvar_anno"):
        if "compare" in res.get(k, {}):
            print(k, {f: {x: y for x, y in r.items() if x != "examples"} for f, r in res[k]["compare"].items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
