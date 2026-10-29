#!/usr/bin/env python3
"""VEPYR in the Colima VM, in a container like VEP (thesis item D7, 2026-10-30).

In the timing series VEPYR ran natively on macOS and VEP in Docker inside the Colima VM
(vz, virtiofs, 10 CPU, 28 GiB). Here VEPYR runs in the same VM with the same kind of mounts
(cache, input and FASTA read-only over virtiofs, output to a virtiofs directory).

  build  Linux aarch64 wheel of vepyr 479d8f4 (git archive of the studied commit, Cargo.lock
         --locked) built with maturin 1.15.0 --release, as the native wheel, in
         rust:1.91-bookworm (rust-toolchain.toml: 1.91); image vepyr-vm:479d8f4 =
         python:3.14-slim-bookworm + the wheel + the native venv's pyarrow, polars, tqdm.
  run    chr22 and chr1 baseline (everything, hgvs, FASTA) with workers 1 and 10, REPEATS
         rounds, order shuffled in each round: wall time of `docker run` (/usr/bin/time -p,
         as for VEP), annotation time inside the container (perf_counter around
         vepyr.annotate, as measure_once.py), CPU (getrusage), peak RSS (ru_maxrss) and the
         body SHA-256, which must equal the native output (chr1: e9703e5e...).

Usage: vepyr_in_vm.py build OUT_DIR
       vepyr_in_vm.py run OUT_DIR [REPEATS]
"""
from __future__ import annotations

import json
import random
import shlex
import subprocess
import sys
import time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
COMMIT = "479d8f4"
TAG = f"vepyr-vm:{COMMIT}"
RUST = "rust:1.91-bookworm"
PYIMG = "python:3.14-slim-bookworm"
DEPS = "pyarrow==25.0.1 polars==1.44.2 tqdm==4.70.1"
INPUT = {"chr22": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz",
         "chr1": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr1.vcf.gz"}
NATIVE_BODY = {"chr1": "e9703e5e71f9e85438575d0382a6ea797b668b09db799a9e7825e602c1f81ffa"}
FASTA = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
SEED = 20261030
CHILD = r'''
import hashlib, json, resource, sys, time, vepyr
vcf, out, workers = sys.argv[1], sys.argv[2], int(sys.argv[3])
t0 = time.perf_counter()
vepyr.annotate(vcf=vcf, cache_dir="/cache", everything=True, hgvs=True, reference_fasta="/ref/%s",
               workers=workers, output_vcf=out, show_progress=False)
dt = time.perf_counter() - t0
u = resource.getrusage(resource.RUSAGE_SELF)
h = hashlib.sha256(); n = 0
with open(out, "rb") as f:
    for line in f:
        if not line.startswith(b"#"):
            h.update(line); n += 1
print(json.dumps({"annotation_seconds": dt, "user_s": u.ru_utime, "sys_s": u.ru_stime,
                  "max_rss_gib": u.ru_maxrss * 1024 / 1024**3, "body_sha256": h.hexdigest(), "records": n,
                  "vepyr": vepyr.__version__}))
''' % FASTA


def sh(cmd: list[str], log: Path) -> None:
    with open(log, "a") as f:
        f.write("$ " + shlex.join(cmd) + "\n"); f.flush()
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)


def build(out: Path) -> None:
    src, log = out / "src", out / "build.log"
    if not src.exists():
        src.mkdir(parents=True)
        subprocess.run(f"git -C {D}/tools/vepyr archive {COMMIT} | tar -x -C {shlex.quote(str(src))}", shell=True,
                       check=True)
    t0 = time.time()
    sh(["docker", "run", "--rm", "-v", f"{src}:/io", "-v", "vepyr-vm-cargo:/usr/local/cargo/registry", "-w", "/io",
        RUST, "bash", "-c",
        "apt-get update -qq && apt-get install -y -qq python3-venv cmake clang pkg-config protobuf-compiler"
        " > /dev/null && python3 -m venv /tmp/m && /tmp/m/bin/pip install -q maturin==1.15.0"
        " && /tmp/m/bin/maturin build --release --locked --compatibility linux -i python3 -o /io/dist"], log)
    wheel = next((src / "dist").glob("vepyr-*.whl"))
    img = out / "image"
    img.mkdir(exist_ok=True)
    (img / wheel.name).write_bytes(wheel.read_bytes())
    (img / "Dockerfile").write_text(f"FROM {PYIMG}\nCOPY {wheel.name} /tmp/\n"
                                    f"RUN pip install --no-cache-dir {DEPS} /tmp/{wheel.name}\n")
    sh(["docker", "build", "-t", TAG, str(img)], log)
    info = {"commit": COMMIT, "wheel": wheel.name, "build_seconds": time.time() - t0,
            "wheel_sha256": subprocess.run(["shasum", "-a", "256", str(wheel)], capture_output=True,
                                           text=True).stdout.split()[0],
            "image_id": subprocess.run(["docker", "image", "inspect", "--format", "{{.Id}}", TAG],
                                       capture_output=True, text=True).stdout.strip(),
            "base_images": {i: subprocess.run(["docker", "image", "inspect", "--format", "{{.Id}}", i],
                                              capture_output=True, text=True).stdout.strip() for i in (RUST, PYIMG)}}
    (out / "build.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info, indent=1))


def run_one(out: Path, chrom: str, workers: int, rep: int) -> dict:
    work = out / "work"
    work.mkdir(exist_ok=True)
    inp, t = INPUT[chrom], out / "time_tmp.txt"
    cmd = ["/usr/bin/time", "-p", "-o", str(t), "docker", "run", "--rm", "--user", "501:20", "--env", "HOME=/tmp",
           "-v", f"{D}/cache/116_GRCh38_merged:/cache:ro", "-v", f"{D}/input:/ref:ro",
           "-v", f"{inp.parent}:/input:ro", "-v", f"{work}:/output", TAG,
           "python", "-c", CHILD, f"/input/{inp.name}", "/output/vepyr.vcf", str(workers)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    (work / "vepyr.vcf").unlink(missing_ok=True)
    r = {"chrom": chrom, "workers": workers, "repeat": rep, "returncode": p.returncode,
         "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
         "docker_run_wall_s": float(next(l.split()[1] for l in t.read_text().splitlines() if l.startswith("real")))}
    if p.returncode == 0:
        r.update(json.loads(p.stdout.strip().splitlines()[-1]))
        if chrom in NATIVE_BODY:
            r["body_equal_native"] = r["body_sha256"] == NATIVE_BODY[chrom]
    else:
        r["stderr_tail"] = p.stderr[-1500:]
    return r


def run(out: Path, reps: int) -> int:
    (out / "environment.json").write_text(json.dumps({
        "image_id": subprocess.run(["docker", "image", "inspect", "--format", "{{.Id}}", TAG], capture_output=True,
                                   text=True).stdout.strip(),
        "colima": [l for l in (Path.home() / ".colima/default/colima.yaml").read_text().splitlines()
                   if l.split(":")[0] in ("cpu", "memory", "vmType", "mountType")]}, indent=1) + "\n")
    rng = random.Random(SEED)
    for chrom in ("chr22", "chr1"):  # warm-up, excluded
        run_one(out, chrom, 1, 0)
    for rep in range(1, reps + 1):
        cells = [(c, w) for c in ("chr22", "chr1") for w in (1, 10)]
        rng.shuffle(cells)
        for chrom, w in cells:
            r = run_one(out, chrom, w, rep)
            with (out / "runs.jsonl").open("a") as f:
                f.write(json.dumps(r) + "\n")
            print(json.dumps({k: v for k, v in r.items() if k != "stderr_tail"}), flush=True)
            if r["returncode"]:
                print(r["stderr_tail"])
                return 1
            time.sleep(5)
    return 0


def main() -> int:
    mode, out = sys.argv[1], Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if mode == "build":
        build(out)
        return 0
    if mode == "run":
        return run(out, int(sys.argv[3]) if len(sys.argv) > 3 else 5)
    raise SystemExit(__doc__)


if __name__ == "__main__":
    raise SystemExit(main())
