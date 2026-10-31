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
  stage  copy the chr1 and chr22 inputs, their core cache fragments (+ chrom_manifest.json) and the FASTA
         into the docker volume VOLUME (ext4 inside the VM), so that no file is read over virtiofs.
  run-local / grid-local  as run / grid, with cache, FASTA, input and output in that volume
         (separates the cost of virtiofs from the cost of the VM itself).
  grid   item D13: chr1 baseline with workers p in {1, 2, 4, 8, 10} on c in {1, 2, 4, 10} VM
         CPUs (docker --cpuset-cpus 0..c-1), REPEATS rounds in shuffled order. The engine's
         threads are not bounded by `workers` alone (spawn_blocking and std::thread besides the
         tokio runtime), so the core count is limited from outside; workers=p on one CPU is
         the single-core time of the partitioned (tabix-indexed) path.

Usage: vepyr_in_vm.py build OUT_DIR
       vepyr_in_vm.py stage OUT_DIR
       vepyr_in_vm.py run-local|grid-local OUT_DIR [REPEATS]
       vepyr_in_vm.py run OUT_DIR [REPEATS]
       vepyr_in_vm.py grid OUT_DIR [REPEATS]
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
VOLUME = "vepyr-local-479d8f4"
LOCAL = False  # set by the *-local modes
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


def run_one(out: Path, chrom: str, workers: int, rep: int, cpus: int | None = None) -> dict:
    work = out / "work"
    work.mkdir(exist_ok=True)
    inp, t = INPUT[chrom], out / "time_tmp.txt"
    limit = ["--cpuset-cpus", f"0-{cpus - 1}"] if cpus else []
    if LOCAL:
        mounts = [a for sub, dst, ro in (("cache", "/cache", 1), ("ref", "/ref", 1), ("input", "/input", 1), ("out", "/output", 0))
                  for a in ("--mount", f"type=volume,src={VOLUME},dst={dst},volume-subpath={sub}" + (",readonly" if ro else ""))]
    else:
        mounts = ["-v", f"{D}/cache/116_GRCh38_merged:/cache:ro", "-v", f"{D}/input:/ref:ro",
                  "-v", f"{inp.parent}:/input:ro", "-v", f"{work}:/output"]
    cmd = ["/usr/bin/time", "-p", "-o", str(t), "docker", "run", "--rm", "--user", "501:20", "--env", "HOME=/tmp",
           *limit, *mounts, TAG,
           "python", "-c", CHILD, f"/input/{inp.name}", "/output/vepyr.vcf", str(workers)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    (work / "vepyr.vcf").unlink(missing_ok=True)
    r = {"chrom": chrom, "workers": workers, "cpus": cpus, "repeat": rep, "returncode": p.returncode,
         "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
         "docker_run_wall_s": float(next(l.split()[1] for l in t.read_text().splitlines() if l.startswith("real")))}
    if p.returncode == 0:
        r.update(json.loads(p.stdout.strip().splitlines()[-1]))
        if chrom in NATIVE_BODY:
            r["body_equal_native"] = r["body_sha256"] == NATIVE_BODY[chrom]
    else:
        r["stderr_tail"] = p.stderr[-1500:]
    return r


def stage(out: Path) -> None:
    """Copy inputs, chr1/chr22 cache fragments and the FASTA into the volume; record what was copied."""
    cache = D / "cache/116_GRCh38_merged"
    ents = sorted(p.name for p in cache.iterdir() if p.is_dir())
    files = [f"{e}/{n}" for e in ents for n in ("chrom_manifest.json", "chr1.parquet", "chr22.parquet") if (cache / e / n).exists()]
    script = ("set -e; rm -rf /v/cache /v/ref /v/input /v/out; mkdir -p /v/cache /v/ref /v/input /v/out; "
              + "; ".join(f"mkdir -p /v/cache/{e}" for e in ents) + "; "
              + "; ".join(f"cp /c/{f} /v/cache/{f}" for f in files) + "; "
              + f"cp /r/{FASTA} /r/{FASTA}.fai /v/ref/; "
              + "; ".join(f"cp /i/{p.name} /i/{p.name}.tbi /v/input/" for p in INPUT.values())
              + "; chown -R 501:20 /v/out; du -sb /v/cache /v/ref /v/input")
    t0 = time.time()
    p = subprocess.run(["docker", "run", "--rm", "-v", f"{VOLUME}:/v", "-v", f"{cache}:/c:ro", "-v", f"{D}/input:/r:ro",
                        "-v", f"{INPUT['chr1'].parent}:/i:ro", TAG, "sh", "-c", script], capture_output=True, text=True, check=True)
    info = {"volume": VOLUME, "files": files, "seconds": time.time() - t0, "du": p.stdout.split("\n")}
    (out / "stage.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info))


def run(out: Path, reps: int, grid: bool = False) -> int:
    (out / "environment.json").write_text(json.dumps({
        "mode": f"local volume {VOLUME}" if LOCAL else "virtiofs mounts",
        "image_id": subprocess.run(["docker", "image", "inspect", "--format", "{{.Id}}", TAG], capture_output=True,
                                   text=True).stdout.strip(),
        "colima": [l for l in (Path.home() / ".colima/default/colima.yaml").read_text().splitlines()
                   if l.split(":")[0] in ("cpu", "memory", "vmType", "mountType")]}, indent=1) + "\n")
    rng = random.Random(SEED)
    for chrom in ("chr22", "chr1"):  # warm-up, excluded
        run_one(out, chrom, 1, 0)
    for rep in range(1, reps + 1):
        if grid:
            cells = [("chr1", w, c) for w in (1, 2, 4, 8, 10) for c in (1, 2, 4, 10)]
        else:
            cells = [(c, w, None) for c in ("chr22", "chr1") for w in (1, 10)]
        rng.shuffle(cells)
        for chrom, w, cpus in cells:
            r = run_one(out, chrom, w, rep, cpus)
            with (out / ("grid_runs.jsonl" if grid else "runs.jsonl")).open("a") as f:
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
    if mode == "stage":
        stage(out)
        return 0
    if mode in ("run-local", "grid-local"):
        global LOCAL
        LOCAL = True
        mode = mode[:-6]
    if mode in ("run", "grid"):
        return run(out, int(sys.argv[3]) if len(sys.argv) > 3 else (5 if mode == "run" else 3), mode == "grid")
    raise SystemExit(__doc__)


if __name__ == "__main__":
    raise SystemExit(main())
