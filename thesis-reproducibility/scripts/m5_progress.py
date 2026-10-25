#!/usr/bin/env python3
"""Progress report for the thesis measurement plan (VEP repeats + M5 completion batch).

Prints a compact status and rewrites the block between the STAN markers in
PLAN_DOKONCZENIA.md, so the plan file always shows the current state.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
VEP = D / "runs/vep_full_repeats_m5_20261023"
OUT = D / "runs/m5_completion_20261024"
PLAN = Path("/Users/lukaszjezapkowicz/Desktop/magisterka/praca/PLAN_DOKONCZENIA.md")
VEP_MEAN_S = 28469.05  # repeat 1 (2026-09-24); replaced by the mean of finished repeats when available


def ts(s: str) -> datetime:
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


def fmt_dur(sec: float) -> str:
    sec = max(0, int(sec))
    h, m = divmod(sec // 60, 60)
    return f"{h} h {m:02d} min" if h else f"{m} min"


def fmt_local(dt: datetime) -> str:
    return dt.astimezone().strftime("%a %d.%m %H:%M")  # system time zone (CEST/CET switch on 25.10)


def vep_status(now: datetime) -> tuple[list[str], datetime | None]:
    lines = []
    walls = []
    current = None
    for i in range(2, 6):
        d = VEP / f"repeat_{i}"
        if not d.exists():
            continue
        t = d / "time.txt"
        m = re.search(r"^real\s+([0-9.]+)", t.read_text(), re.M) if t.exists() else None
        if m:
            rc = (d / "exit_code.txt").read_text().strip() if (d / "exit_code.txt").exists() else "?"
            sha = (d / "output_body.sha256").read_text().strip()[:12] if (d / "output_body.sha256").exists() else "—"
            peak = max((int(x) for x in re.findall(r"VmHWM:(\d+)", (d / "rss_samples.txt").read_text())), default=0) \
                if (d / "rss_samples.txt").exists() else 0
            walls.append(float(m.group(1)))
            lines.append(f"  - powtórzenie {i}: zakończone, {float(m.group(1)):.0f} s ({fmt_dur(float(m.group(1)))}), "
                         f"kod {rc}, szczyt VmHWM {peak * 1024 / 1e9:.2f} GB, sha256 treści {sha}…")
        elif (d / "started_utc.txt").exists():
            current = (i, ts((d / "started_utc.txt").read_text().strip()))
    mean = sum(walls) / len(walls) if walls else VEP_MEAN_S
    finished = (VEP / "repeat_5/exit_code.txt").exists() and (VEP / "summary.json").exists()
    eta_end = None
    if current and not finished:
        i, start = current
        el = (now - start).total_seconds()
        vcf = VEP / f"repeat_{i}/HG002_vep_forknone.vcf"
        size = vcf.stat().st_size / 1024**3 if vcf.exists() else 0  # repeat 1 output: 27.40 GiB
        remaining_repeats = 5 - i
        eta_end = start + timedelta(seconds=mean * (remaining_repeats + 1) + 300 * (remaining_repeats + 1))
        lines.append(f"  - powtórzenie {i}: W TOKU od {fmt_local(start)}, upłynęło {fmt_dur(el)}, "
                     f"wynik ~{size:.1f}/27,4 GiB ({100 * size / 27.4:.0f}%), "
                     f"szac. koniec tego powtórzenia {fmt_local(start + timedelta(seconds=mean))}")
    head = "zakończone" if finished else (f"w toku, szac. koniec wszystkich: {fmt_local(eta_end)}" if eta_end else "w toku")
    return [f"- **Powtórzenia VEP (pełny genom, 2–5)**: {head}"] + lines, eta_end


def batch_status(now: datetime, vep_end: datetime | None) -> list[str]:
    st_file = OUT / "status.json"
    if not st_file.exists():
        return ["- **Kolejka pomiarów M5**: nie uruchomiona"]
    st = json.loads(st_file.read_text())
    alive = False
    if st.get("pid"):
        alive = subprocess.run(["ps", "-p", str(st["pid"])], capture_output=True).returncode == 0
    phase = st.get("phase")
    lines = [f"- **Kolejka pomiarów M5**: faza `{phase}`, proces {'działa' if alive else 'NIE DZIAŁA'}"]
    t = now if phase != "waiting_for_vep_reruns" else (vep_end or now)
    remaining = 0.0
    for s in st["steps"]:
        state = s["state"]
        extra = ""
        if state == "running":
            el = (now - ts(s["started"])).total_seconds() / 60
            left = max(s["eta_min"] - el, 1)
            remaining += left
            extra = f", trwa {fmt_dur(el * 60)} (szac. {s['eta_min']} min)"
        elif state == "pending":
            remaining += s["eta_min"]
            extra = f", szac. {s['eta_min']} min"
        elif state == "done":
            extra = f", {s.get('minutes', '?')} min"
        elif state in ("failed", "skipped"):
            extra = f" — {s.get('detail', '')[:160]}"
        mark = {"done": "[x]", "running": "[>]", "failed": "[!]", "skipped": "[-]"}.get(state, "[ ]")
        lines.append(f"  - {mark} {s['id']} {s['name']}{extra}")
    if phase != "finished":
        lines.append(f"  - szacowany koniec kolejki: {fmt_local(t + timedelta(minutes=remaining))} "
                     f"(pozostało ok. {fmt_dur(remaining * 60)} pracy kolejki)")
    return lines


def phase2_status(now: datetime) -> list[str]:
    log = OUT / "phase2.log"
    if not log.exists():
        return ["- **Faza 2** (S10, potem Colima + S3, S7, powtórzenia VEP 4–5, S9): nie uruchomiona"]
    lines = log.read_text().strip().splitlines()
    last = lines[-1] if lines else ""
    alive = subprocess.run(["pgrep", "-f", "m5_phase2.sh"], capture_output=True).returncode == 0
    out = [f"- **Faza 2** (S10, potem Colima + S3, S7, powtórzenia VEP 4–5, S9): proces {'działa' if alive else 'nie działa'}; "
           f"ostatni wpis: `{last[:140]}`"]
    if "waiting for phase 1" in last:
        out.append("  - czeka na koniec fazy 1; szacunkowo potrwa ok. 23 h od startu (S10 ~2 h, S3+S7 ~25 min, VEP 2 × ~7,9 h, S9 ~4 h)")
    return out


def main() -> None:
    now = datetime.now(timezone.utc)
    free = shutil.disk_usage(D).free / 1e9
    vl, vep_end = vep_status(now)
    bl = batch_status(now, vep_end)
    block = [f"_Stan na {fmt_local(now)} (aktualizowane automatycznie przez tools/m5_progress.py)_", "", *vl, *bl, *phase2_status(now),
             f"- Wolne miejsce na dysku: {free:.0f} GB"]
    text = "\n".join(block)
    print(text)
    if PLAN.exists():
        p = PLAN.read_text()
        new = re.sub(r"(<!-- STAN:START -->\n).*?(\n<!-- STAN:END -->)", lambda m: m.group(1) + text + m.group(2), p, flags=re.S)
        PLAN.write_text(new)


if __name__ == "__main__":
    main()
