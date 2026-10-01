#!/usr/bin/env python3
"""画 model-vs-fit 的对比图：I(q) 叠加 + 标准化残差 + 判据 χ² 条形图。

读 <产物目录>/comparison.json 与其中各条曲线的 .fit 文件（以及 data_used.dat），
写 <产物目录>/model-vs-fit.png。只读产物、不重算拟合（重算在 model-vs-fit.py 里）。
图内文字用英文 —— matplotlib 默认字体没有中文字形，中文会画成方框（控制台消息仍可用中文）。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_sasdata(path: Path) -> np.ndarray:
    rows = []
    for line in Path(path).read_text(errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith(("#", "!", ";", ">")):
            continue
        parts = s.split()
        if len(parts) < 3:
            continue
        try:
            vals = [float(x) for x in parts[:4]]
        except ValueError:
            continue
        rows.append(vals)
    width = min(len(r) for r in rows)
    return np.array([r[:width] for r in rows], dtype=float)


def align_values(qf, yf, qt):
    """与 model-vs-fit.py 同一套取值规则：目标点能落在源网格上就直接取，否则插值。"""
    qf, yf, qt = np.asarray(qf, float), np.asarray(yf, float), np.asarray(qt, float)
    if len(qf) >= 2:
        j = np.clip(np.searchsorted(qf, qt), 1, len(qf) - 1)
        left, right = qf[j - 1], qf[j]
        k = np.where(np.abs(qt - left) <= np.abs(right - qt), j - 1, j)
        if np.all(np.abs(qf[k] - qt) <= 1e-6 * np.maximum(np.abs(qt), 1e-6)):
            return yf[k]
    out = np.interp(qt, qf, yf)
    pos = yf > 0
    if pos.sum() >= 2:
        qp, yp = qf[pos], yf[pos]
        m = (qt >= qp[0]) & (qt <= qp[-1])
        if m.any():
            out[m] = np.exp(np.interp(qt[m], qp, np.log(yp)))
    return out


def main() -> int:
    outdir = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    payload = json.loads((outdir / "comparison.json").read_text())
    data = load_sasdata(Path(payload["data_used"]))
    q, I, err = data[:, 0], data[:, 1], data[:, 2]

    curves = [c for c in payload["curves"]
              if c.get("fit_file") and Path(c["fit_file"]).exists()]
    cmap = plt.get_cmap("tab10")
    fig = plt.figure(figsize=(11.5, 9.5))
    gs = fig.add_gridspec(3, 1, height_ratios=[2.2, 1.4, 1.3], hspace=0.35)
    ax1, ax2, ax3 = (fig.add_subplot(gs[i]) for i in range(3))

    m = err > 0
    ax1.errorbar(q[m], I[m], yerr=err[m], fmt="o", ms=2.0, lw=0.5, alpha=0.45, color="black",
                 label=f"data: {Path(payload['data_original']).name}", zorder=1)
    for i, c in enumerate(curves):
        arr = load_sasdata(Path(c["fit_file"]))
        if arr.shape[1] < 4:            # 3 列文件（如 GNOM .out）没有拟合列，画不了
            continue
        col = cmap(i % 10)
        v = c.get("primary_chi2")
        stem = Path(c["fit_file"]).stem          # 用文件名做图例：ASCII，且不会缺字形
        lab = (f"{stem}  (chi2={v:.2f}, {c.get('primary_source')})" if v is not None else stem)
        ax1.plot(arr[:, 0], arr[:, 3], "-", color=col, lw=1.4, label=lab, zorder=3)
        ax2.plot(q[m], (I[m] - align_values(arr[:, 0], arr[:, 3], q[m])) / err[m],
                 "-", color=col, lw=1.0, alpha=0.85)

    for ax in (ax1, ax2):
        ax.set_xlim(q.min(), q.max())
        if payload.get("qmin"):
            ax.axvline(payload["qmin"], color="grey", ls=":", lw=1)
        if payload.get("qmax"):
            ax.axvline(payload["qmax"], color="grey", ls=":", lw=1)
    ax1.set_yscale("log")
    ax1.set_ylabel("I(q)  [arbitrary units]")
    ax1.set_title("CRYSOL model curves vs the fits already present in processed/", fontsize=11)
    ax1.legend(fontsize=8, loc="best")
    ax1.grid(alpha=0.25)

    ax2.axhline(0, color="black", lw=0.8)
    for y in (-1, 1):
        ax2.axhline(y, color="grey", ls="--", lw=0.6)
    for y in (-3, 3):
        ax2.axhline(y, color="red", ls=":", lw=0.6)
    ax2.set_ylabel("(I_data - I_fit) / sigma")
    ax2.set_xlabel("q  [1/A]")
    ax2.grid(alpha=0.25)

    labels, vals, colors = [], [], []
    for i, c in enumerate(curves):
        v = c.get("primary_chi2")
        if v is None:
            continue
        labels.append(f"{Path(c['fit_file']).stem[:34]} ({c.get('primary_source')})")
        vals.append(v)
        colors.append(cmap(i % 10))
    ypos = np.arange(len(vals))
    ax3.barh(ypos, vals, color=colors, alpha=0.85)
    for y, v in zip(ypos, vals):
        ax3.text(v, y, f" {v:.2f}", va="center", fontsize=8)
    ax3.set_yticks(ypos, labels, fontsize=8)
    ax3.invert_yaxis()
    band = payload.get("chi2_band_for_n_data") or [None, None]
    if band[0]:
        ax3.axvspan(band[0], band[1], color="green", alpha=0.15,
                    label=f"chi2_n 95% band (n={payload.get('chi2_band_for_n')})")
        ax3.legend(fontsize=8)
    ax3.axvline(1.0, color="black", ls=":", lw=1)
    ax3.set_xlabel("decisive reduced chi2 (log scale; see README for which column)")
    if vals:                                    # 对数轴：一条 χ²=155 的曲线不会把 1–5 的差异压平
        ax3.set_xscale("log")
        ax3.set_xlim(min(vals) * 0.6, max(vals) * 1.8)
    ax3.grid(alpha=0.25, axis="x")

    fig.suptitle("Model -> CRYSOL curve vs existing fits (like-for-like comparison)",
                 fontsize=12, y=0.995)
    out = outdir / "model-vs-fit.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"{out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
