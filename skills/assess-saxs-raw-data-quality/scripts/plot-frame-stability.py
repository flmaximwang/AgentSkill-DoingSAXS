#!/usr/bin/env python
"""把 raw-qc-frames.py 的数字画成"逐帧稳定性"图。

    <RAW python> plot-frame-stability.py <root> [--qc <root>/_rawqc]

读 <qc>/<sample>/qc.json，每个样品一个小格：每帧的 I(低 q) 除以该 run 的中位数（蓝=样品，红=control），
绿色右轴是透射。蓝线平在 1 = 帧可复现；爬坡 = 漂移（辐照损伤、沉降、毛细管积垢）。
标题里带该样品的结论与三个关键数（上翘、缩放因子、对比度）——图可以直接拿去讨论。
"""
from __future__ import annotations

import argparse, glob, json, os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LAB = {"好": "GOOD", "可用": "OK", "勉强可用": "WARN", "不可用": "BAD"}   # 默认字体没有中文


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("root", help="the root given to raw-qc-frames.py (for --qc default)")
    ap.add_argument("--qc", default=None, help="default <root>/_rawqc")
    ap.add_argument("--out", default=None, help="default <qc>/raw_qc.png")
    ap.add_argument("--ncol", type=int, default=3)
    args = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(args.root))
    qc = os.path.abspath(os.path.expanduser(args.qc)) if args.qc else os.path.join(root, "_rawqc")

    files = sorted(glob.glob(os.path.join(qc, "*", "qc.json")))
    if not files:
        raise SystemExit("no qc.json under %s - run raw-qc-frames.py first" % qc)
    n, ncol = len(files), max(1, args.ncol)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 2.4 * nrow), squeeze=False)
    for ax, f in zip(axes.ravel(), files):
        d = json.load(open(f))
        info, frames = d["info"], d["frames"]
        for kind, col in (("sample", "tab:blue"), ("control", "tab:red")):
            v = np.array([r["i_lo"] for r in frames if r["kind"] == kind], float)
            if v.size:
                ax.plot(np.arange(v.size), v / np.median(v), "o-", ms=2.5, lw=0.8, color=col, label=kind)
        tx = np.array([r["tx"] for r in frames if r["kind"] == "sample" and r["tx"] == r["tx"]], float)
        if tx.size:
            ax2 = ax.twinx()
            ax2.plot(np.arange(tx.size), tx / np.median(tx), "-", lw=0.8, color="green", alpha=0.6)
            ax2.set_ylim(0.8, 1.2)
            ax2.tick_params(labelsize=6, colors="green")
        ax.axhline(1, color="k", lw=0.5, ls=":")
        ax.axhspan(0.85, 1.15, color="k", alpha=0.05)
        ax.set_ylim(0.4, 1.6)
        ax.set_title("%s | %s | upturn=%.2f scale=%.3f contr=%.1f%%"
                     % (info["sample"], LAB.get(info["verdict"], info["verdict"]),
                        info["lowq_upturn"], info["control_scale"], info["contrast_lowq"] * 100),
                     fontsize=8)
        ax.tick_params(labelsize=6)
        ax.set_xlabel("frame", fontsize=7)
        ax.set_ylabel("I(lo q) / run median", fontsize=7)
        ax.legend(fontsize=6, loc="lower left")
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    plt.tight_layout()
    png = os.path.abspath(os.path.expanduser(args.out)) if args.out else os.path.join(qc, "raw_qc.png")
    plt.savefig(png, dpi=105)
    print(png)


if __name__ == "__main__":
    main()
