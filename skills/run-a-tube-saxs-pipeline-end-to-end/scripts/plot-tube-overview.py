#!/usr/bin/env python
"""Overview figure for a batch of tube-pipeline runs: <root>/_summary/overview.png

Panels: (a) subtracted I(q) of every sample, (b) Kratky, (c) I(0) per sample on a log axis
(and, when the sample names form a numbered series, I0 vs series order), (d) Rg - auto-Guinier
vs IFT real space - with a marker for runs whose IFT passed the trust gates.

Only reads numbers the pipeline already produced (tables/*.json, tables/*.csv) and the .dat
profiles; it fits nothing.
"""
import argparse, csv, glob, json, os, re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument("root", help="directory holding one sub-directory per run")
args = ap.parse_args()
root = os.path.abspath(os.path.expanduser(args.root))
out = os.path.join(root, "_summary")
os.makedirs(out, exist_ok=True)


def load_dat(path):
    d = np.loadtxt(path, comments="#", ndmin=2)
    return d[:, 0], d[:, 1], d[:, 2]


runs = []
for d in sorted(glob.glob(os.path.join(root, "*"))):
    sj, dat = os.path.join(d, "summary.json"), os.path.join(
        d, "profiles", "03_subtracted", "subtracted.dat")
    if not (os.path.isdir(d) and os.path.exists(sj) and os.path.exists(dat)):
        continue
    s = json.load(open(sj))
    runs.append(dict(name=s.get("sample", os.path.basename(d)), dir=d, s=s, dat=dat))

fig = plt.figure(figsize=(15, 10))
ax_i, ax_k, ax_c, ax_r = (fig.add_subplot(2, 2, i) for i in (1, 2, 3, 4))
for r in runs:
    q, I, E = load_dat(r["dat"])
    m = (I > 0) & np.isfinite(I)
    ax_i.errorbar(q[m], I[m], yerr=E[m], fmt=".-", ms=2, lw=0.7, label=r["name"])
    ax_k.plot(q[m], (q ** 2 * I)[m], ".-", ms=2, lw=0.7, label=r["name"])
ax_i.set_xscale("log"); ax_i.set_yscale("log")
ax_i.set_xlabel("q (1/A)"); ax_i.set_ylabel("I"); ax_i.set_title("subtracted (log-log)")
ax_i.legend(fontsize=7, ncol=2)
ax_k.set_xlabel("q"); ax_k.set_ylabel("q^2 I"); ax_k.set_title("Kratky")
ax_k.legend(fontsize=7, ncol=2)

names = [r["name"] for r in runs]
i0 = [(r["s"].get("guinier_auto") or {}).get("rg") and
      float(json.load(open(os.path.join(r["dir"], "tables", "guinier_results.json")))["auto"]["i0"])
      for r in runs]
ax_c.bar(range(len(runs)), i0)
ax_c.set_yscale("log"); ax_c.set_xticks(range(len(runs)))
ax_c.set_xticklabels(names, rotation=60, ha="right", fontsize=8)
ax_c.set_ylabel("I(0) from auto-Guinier"); ax_c.set_title("I(0) per sample (log)")

xa = [float((r["s"].get("guinier_auto") or {}).get("rg", np.nan)) for r in runs]
xj = [float((r["s"].get("ift") or {}).get("runs") and next(
    (x.get("rg_realspace", np.nan) for x in r["s"]["ift"]["runs"]
     if x.get("tag") == r["s"]["ift"].get("chosen") and "failed" not in x), np.nan)) for r in runs]
tr = [bool((r["s"].get("ift") or {}).get("trusted")) for r in runs]
ax_r.plot(range(len(runs)), xa, "o-", label="auto-Guinier Rg")
ax_r.plot(range(len(runs)), xj, "s--", label="IFT real-space Rg", alpha=0.8)
for k, t in enumerate(tr):
    if t:
        ax_r.plot(k, xa[k], "o", ms=11, mfc="none", mec="green", mew=1.5)
ax_r.set_xticks(range(len(runs)))
ax_r.set_xticklabels(names, rotation=60, ha="right", fontsize=8)
ax_r.set_ylabel("Rg (A)"); ax_r.set_title("Rg (green ring = IFT passed trust gates)")
ax_r.legend(fontsize=8)

plt.tight_layout()
png = os.path.join(out, "overview.png")
plt.savefig(png, dpi=110)
print(png)
