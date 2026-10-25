"""Figure fig:plugin-overhead: O_x per plugin on chr1/chr22 from the 2026-10-25 plugin series.

Usage: python3 fig_plugin_overhead.py <raw_runs.jsonl> <output.png>
"""
import json, statistics as st, sys
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
src, out = sys.argv[1], sys.argv[2]
runs = [json.loads(l) for l in open(src) if l.strip()]
order = [("clinvar","ClinVar"),("alphamissense","AlphaMissense"),("dbnsfp","dbNSFP"),("spliceai","SpliceAI"),("cadd","CADD")]
def mean_t(ch, cfg):
    return st.mean(r["annotation_seconds"] for r in runs if r["chrom"] == ch and r["config"] == cfg)
def ox(ch):
    tb = mean_t(ch, "baseline")
    return [round(100*(mean_t(ch, k)-tb)/tb, 1) for k,_ in order]
v1, v22 = ox("chr1"), ox("chr22")
print(v1, v22)
fig, ax = plt.subplots(figsize=(8.9, 3.9), dpi=200)
x = range(len(order)); w = 0.36
b1 = ax.bar([i-w/2 for i in x], v1, w, color="#b35c0b", label="chr1")
b2 = ax.bar([i+w/2 for i in x], v22, w, color="#356b8a", label="chr22")
for bars in (b1, b2):
    for b in bars:
        ax.annotate(f"{b.get_height():.1f}".replace(".", ","), (b.get_x()+b.get_width()/2, b.get_height()), xytext=(0,3), textcoords="offset points", ha="center", va="bottom", fontsize=8)
ax.set_xticks(list(x)); ax.set_xticklabels([n for _,n in order], rotation=20, ha="right")
ax.set_ylabel(r"Narzut $O_x$ [%]")
ax.set_ylim(0, max(75, 10*(int(max(v1+v22))//10+2)))
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",")))
ax.axhline(0, color="#555", lw=0.8)
ax.yaxis.grid(True, color="#e6e6e6"); ax.set_axisbelow(True)
ax.legend(loc="upper left", ncol=2, frameon=False)
fig.tight_layout(); fig.savefig(out)
