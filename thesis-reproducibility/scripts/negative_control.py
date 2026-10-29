#!/usr/bin/env python3
"""Negative control of the two comparators used in the thesis (2026-10-29).

Takes one annotated VCF (VEPYR 0.7.0, chr22, all five plugins), writes copies with a
known number of injected differences of one class each, and compares every copy with
the original using
  * the E2E comparator of the vepyr repository (e2e-testing/scripts/comparison,
    compare_vcfs; tab:variant-results), and
  * vepyr-diffly compare-existing (full and --compare-only-plugins;
    tab:plugin-field-results).
The original plays the VEP side, the modified copy the VEPYR side. A comparator passes
a class when it reports exactly the injected differences and nothing else; the
equivalence classes (number formatting, order inside '&' lists, header lines) must
not be reported as strict differences.

Usage: negative_control.py BASE_VCF OUT_DIR
Environment: E2E_DIR (vepyr e2e-testing/scripts), DIFFLY_SRC (vepyr-diffly/src),
DIFFLY_PYTHON (interpreter with polars), K (instances per class, default 20).
"""
from __future__ import annotations

import json
import os
import random
import re
import subprocess
import sys
import time
from pathlib import Path

K = int(os.environ.get("K", "20"))
SEED = 20261029
E2E_DIR = Path(os.environ.get("E2E_DIR", "/Users/lukaszjezapkowicz/workspace/data_vepyr/tools/vepyr/e2e-testing/scripts"))
DIFFLY_SRC = Path(os.environ.get("DIFFLY_SRC", "/Users/lukaszjezapkowicz/workspace/data_vepyr/tools/vepyr-diffly/src"))
DIFFLY_PY = os.environ.get("DIFFLY_PYTHON", "/Users/lukaszjezapkowicz/workspace/data_vepyr/.venv-native/bin/python3")
PLUGINS = ["clinvar", "spliceai", "cadd", "alphamissense", "dbnsfp"]
CSQ_RE = re.compile(r"(^|;)CSQ=([^;\t]*)")


def read_vcf(path: Path):
    header, rows = [], []
    with path.open() as f:
        for line in f:
            (header if line.startswith("#") else rows).append(line.rstrip("\n"))
    fmt = next(h for h in header if h.startswith("##INFO=<ID=CSQ"))
    fields = fmt.split("Format: ")[1].rstrip('">').split("|")
    return header, rows, fields


def csq_of(row: str) -> str:
    return CSQ_RE.search(row.split("\t")[7]).group(2)


def with_csq(row: str, csq: str) -> str:
    cols = row.split("\t")
    cols[7] = CSQ_RE.sub(lambda m: f"{m.group(1)}CSQ={csq}", cols[7], count=1)
    return "\t".join(cols)


def key(row: str) -> str:
    c = row.split("\t")
    return f"{c[0]}:{c[1]}:{c[3]}:{c[4]}"


def bump_last_digit(s: str) -> str | None:
    m = re.search(r"(\d)(?!.*\d)", s)
    if not m:
        return None
    d = str((int(m.group(1)) + 1) % 10)
    return s[: m.start()] + d + s[m.end():]


class Mutator:
    """Each class picks fresh rows (never reused) and returns the expected effect."""

    def __init__(self, header, rows, fields):
        self.header, self.rows, self.f = header, rows, fields
        self.ix = {n: i for i, n in enumerate(fields)}
        self.rng = random.Random(SEED)
        self.used: set[int] = set()

    def pick(self, pred, k=K):
        cand = [i for i, r in enumerate(self.rows) if i not in self.used and pred(r)]
        sel = sorted(self.rng.sample(cand, k))
        self.used.update(sel)
        return sel

    def entries(self, row):
        return [e.split("|") for e in csq_of(row).split(",")]

    def join(self, ents):
        return ",".join("|".join(e) for e in ents)

    def apply(self, cls):
        rows = list(self.rows)
        header = list(self.header)
        exp = {"variants": K}
        F = self.ix
        if cls == "variant_removed":
            sel = self.pick(lambda r: True)
            for i in reversed(sel):
                del rows[i]
        elif cls == "variant_added":
            sel = self.pick(lambda r: len(r.split("\t")[3]) == 1 and len(r.split("\t")[4]) == 1)
            for i in reversed(sel):
                c = rows[i].split("\t")
                alt = next(b for b in "ACGT" if b not in (c[3], c[4]))
                c[4] = alt
                rows.insert(i + 1, "\t".join(c))
        elif cls == "csq_entry_removed":
            sel = self.pick(lambda r: len(self.entries(r)) >= 2)
            for i in sel:
                rows[i] = with_csq(rows[i], self.join(self.entries(rows[i])[:-1]))
        elif cls == "csq_entry_added":
            sel = self.pick(lambda r: self.entries(r)[0][F["Feature"]] != "")
            for i in sel:
                e = self.entries(rows[i])
                extra = list(e[0])
                extra[F["Feature"]] = "ENST99999999999"
                rows[i] = with_csq(rows[i], self.join(e + [extra]))
        elif cls == "consequence_changed":
            sel = self.pick(lambda r: True)
            for i in sel:
                e = self.entries(rows[i])
                e[0][F["Consequence"]] = "synonymous_variant" if e[0][F["Consequence"]] != "synonymous_variant" else "intron_variant"
                rows[i] = with_csq(rows[i], self.join(e))
            exp["field_mismatches"] = {"Consequence": K}
        elif cls == "hgvs_changed":
            sel = self.pick(lambda r: any(bump_last_digit(x[F["HGVSc"]]) for x in self.entries(r)))
            for i in sel:
                e = self.entries(rows[i])
                j = next(j for j, x in enumerate(e) if bump_last_digit(x[F["HGVSc"]]))
                e[j][F["HGVSc"]] = bump_last_digit(e[j][F["HGVSc"]])
                rows[i] = with_csq(rows[i], self.join(e))
            exp["field_mismatches"] = {"HGVSc": K}
        elif cls == "plugin_value_changed":
            fld = "CADD_PHRED"
            sel = self.pick(lambda r: self.entries(r)[0][F[fld]] not in ("",))
            n = 0
            for i in sel:
                e = self.entries(rows[i])
                for x in e:
                    if x[F[fld]]:
                        x[F[fld]] = f"{float(x[F[fld]]) + 1.0:g}"
                        n += 1
                rows[i] = with_csq(rows[i], self.join(e))
            exp["field_mismatches"] = {fld: n}
            exp["plugin"] = "cadd"
        elif cls == "plugin_value_removed":
            fld = next(n for n in self.f if n.startswith("ClinVar") and "CLNSIG" in n)
            sel = self.pick(lambda r: self.entries(r)[0][F[fld]] != "")
            n = 0
            for i in sel:
                e = self.entries(rows[i])
                for x in e:
                    if x[F[fld]]:
                        x[F[fld]] = ""
                        n += 1
                rows[i] = with_csq(rows[i], self.join(e))
            exp["field_mismatches"] = {fld: n}
            exp["plugin"] = "clinvar"
        elif cls == "csq_order_swapped":
            def two_distinct(r):
                e = self.entries(r)
                sig = lambda x: (x[F["Allele"]], x[F["Feature"]], x[F["Consequence"]])  # noqa: E731
                return len(e) >= 2 and sig(e[0]) != sig(e[1])
            sel = self.pick(two_distinct)
            for i in sel:
                e = self.entries(rows[i])
                e[0], e[1] = e[1], e[0]
                rows[i] = with_csq(rows[i], self.join(e))
        elif cls == "equiv_number_format":
            fld = "CADD_PHRED"
            sel = self.pick(lambda r: "." in self.entries(r)[0][F[fld]])
            n = 0
            for i in sel:
                e = self.entries(rows[i])
                for x in e:
                    if "." in x[F[fld]]:
                        x[F[fld]] += "0"
                        n += 1
                rows[i] = with_csq(rows[i], self.join(e))
            exp["format_only"] = {fld: n}
            exp["plugin"] = "cadd"
        elif cls == "equiv_amp_order":
            fld = "Existing_variation"
            sel = self.pick(lambda r: any(len(set(x[F[fld]].split("&"))) > 1 for x in self.entries(r)))
            n = 0
            for i in sel:
                e = self.entries(rows[i])
                for x in e:
                    parts = x[F[fld]].split("&")
                    if len(set(parts)) > 1:
                        x[F[fld]] = "&".join(reversed(parts))
                        n += 1
                rows[i] = with_csq(rows[i], self.join(e))
            exp["order_only"] = {fld: n}
        elif cls == "equiv_header":
            header.insert(1, "##negative_control=header line added; not a data difference")
            header = [h.replace("##VEP=", "##VEP_EDITED=") if h.startswith("##VEP=") else h for h in header]
            exp = {"variants": 0}
        elif cls != "identity":
            raise ValueError(cls)
        else:
            exp = {"variants": 0}
        return header, rows, exp


def write_vcf(path: Path, header, rows):
    with path.open("w") as f:
        f.write("\n".join(header) + "\n")
        f.write("\n".join(rows) + "\n")


def run_e2e(ref: Path, test: Path, ledger: Path) -> dict:
    sys.path.insert(0, str(E2E_DIR))
    from comparison import compare  # noqa: E402
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        res = compare.compare_vcfs(str(test), str(ref), "negative-control", mismatch_ledger_path=str(ledger))
    return res


def run_diffly(ref: Path, test: Path, out: Path, only_plugins: bool) -> dict:
    cmd = [DIFFLY_PY, "-m", "vepyr_diffly.cli", "compare-existing", "--preset", "ensembl_everything",
           "--left-vcf", str(ref), "--right-vcf", str(test), "--output-dir", str(out),
           "--plugins", ",".join(PLUGINS)]
    if only_plugins:
        cmd.append("--compare-only-plugins")
    env = {**os.environ, "PYTHONPATH": str(DIFFLY_SRC)}
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    (out / "stdout.txt").write_text(p.stdout + p.stderr) if out.exists() else None
    if p.returncode != 0:
        return {"returncode": p.returncode, "tail": (p.stdout + p.stderr)[-2000:]}
    return json.loads((out / "summary.json").read_text())


CLASSES = ["identity", "variant_removed", "variant_added", "csq_entry_removed", "csq_entry_added",
           "consequence_changed", "hgvs_changed", "plugin_value_changed", "plugin_value_removed",
           "csq_order_swapped", "equiv_number_format", "equiv_amp_order", "equiv_header"]


def main() -> int:
    base, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    header, rows, fields = read_vcf(base)
    mut = Mutator(header, rows, fields)
    only = os.environ.get("CLASSES")
    results = []
    for cls in CLASSES:
        if only and cls not in only.split(","):
            continue
        d = out / cls
        d.mkdir(exist_ok=True)
        h, r, exp = mut.apply(cls)
        test = d / "modified.vcf"
        write_vcf(test, h, r)
        t0 = time.time()
        e2e = run_e2e(base, test, d / "e2e_ledger.jsonl")
        t1 = time.time()
        dif_full = run_diffly(base, test, d / "diffly_full", only_plugins=False)
        t2 = time.time()
        dif_plug = run_diffly(base, test, d / "diffly_plugins", only_plugins=True)
        t3 = time.time()
        item = {"class": cls, "expected": exp, "e2e": e2e, "diffly_full": dif_full, "diffly_plugins": dif_plug,
                "seconds": {"e2e": t1 - t0, "diffly_full": t2 - t1, "diffly_plugins": t3 - t2}}
        (d / "result.json").write_text(json.dumps(item, indent=1, default=str) + "\n")
        test.unlink()
        results.append({"class": cls, "expected": exp})
        print(cls, json.dumps(exp), flush=True)
    (out / "classes.json").write_text(json.dumps(results, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
