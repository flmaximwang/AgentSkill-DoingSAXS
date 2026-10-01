#!/usr/bin/env python
"""Write a human-readable README.md for one (or many) embed/ result directories.

    /Applications/BioXTASRAW/bin/python write-embed-readme.py <embed-dir> [<embed-dir> ...]

Reads <embed-dir>/embed_results.json (written by embed-model.py) and writes
<embed-dir>/README.md: one-line conclusion, the numbers with their criteria, the file map,
what was NOT done, and how to reproduce.  Nothing is recomputed here - if a number is not
in the JSON it is not in the README.
"""
from __future__ import annotations

import json
import os
import sys


def f(x, nd=2, dash="-"):
    return dash if x is None else ("%.*f" % (nd, x))


def build(d):
    s = d.get("sample", "?")
    inp = d.get("input", {}) or {}
    bead = d.get("bead_model") or {}
    bf = d.get("beads_fit") or {}
    df = d.get("density_fit") or {}
    L = []
    A = L.append

    # ---- one-line conclusion
    concl = []
    if bf.get("nsd") is not None:
        concl.append("珠模型分支：NSD %s（%s 口径）" % (f(bf["nsd"], 3), bf.get("selection")))
    if df.get("best_correlation") is not None:
        concl.append("电子云分支：correlation %s @ %s Å" % (f(df["best_correlation"], 4), df.get("best_resolution")))
        if df.get("pose_spread_max_A") and df["pose_spread_max_A"] > 5:
            concl.append("但至少有两个分数接近的不同摆法（姿态散布 %s Å）" % f(df["pose_spread_max_A"], 1))
    if df.get("status"):
        concl.append("电子云分支：%s" % df["status"])
    if not concl:
        concl.append("两条分支都没跑出结果（见下）")
    A("# %s —— 高分辨模型嵌入珠模型 / 电子云" % s)
    A("")
    A("> 由 `embed-model.py` 生成的分析 + 本 README ｜ 模型：`%s` ｜ 结果目录：`%s`"
      % (d.get("model", "?"), d.get("out_dir", ".")))
    A("")
    A("## 一句话结论")
    A("")
    A("；".join(concl) + "。")
    if d.get("warning"):
        A("")
        A("⚠️ **告警**：%s" % d["warning"])
    A("")

    # ---- numbers
    A("## 1. 关键结果")
    A("")
    A("| 项目 | 数值 | 怎么看 |")
    A("|---|---|---|")
    A("| 用的 IFT run | %s（trusted=%s） | 上游闸门选中的那条；代建珠模型只用它的 q 窗口 |"
      % (inp.get("tag", "-"), inp.get("trusted")))
    A("| P(r) 的 Dmax / Rg | %s Å / %s Å | 珠模型的 Rg/Dmax 应当与它对得上 |"
      % (f(inp.get("dmax"), 1), f(inp.get("rg"), 1)))
    A("| 电子云 | `%s` | DENSS 的最终图（不是 current/support） |"
      % os.path.basename(inp.get("density_map") or "-"))
    g = bead.get("gnom")
    if g:
        A("| GNOM | χ² %s，总估计 %s（%s），Rg(实空间) %s Å | χ² 数百 = P(r) 本身拟合差，珠模型跟着不可信 |"
          % (f(g.get("chi_sq")), f(g.get("total_estimate")), g.get("quality", "-"), f(g.get("rg_out"))))
    dm = bead.get("dammif")
    if dm:
        runs = dm.get("runs") or []
        if runs:
            A("| DAMMIF ×%d（%s） | χ² %s–%s，Rg %s Å，Dmax %s–%s Å | 各模型 Rg 应该彼此接近、并与 P(r) 的 Rg 对得上 |"
              % (dm.get("n_models"), dm.get("mode"),
                 f(min(r["chi2"] for r in runs)), f(max(r["chi2"] for r in runs)),
                 f(runs[0]["rg"]), f(min(r["dmax"] for r in runs), 1), f(max(r["dmax"] for r in runs), 1)))
        da = bead.get("damaver") or {}
        if da.get("mean_nsd") is not None:
            A("| DAMAVER | 平均 NSD %s ± %s | 这是**这批模型彼此**像不像（≠ 本目录的嵌入 NSD） |"
              % (f(da["mean_nsd"]), f(da.get("std_nsd"))))
        A("| 共识珠模型 | `%s`（%s 个 bead） | damfilt = 滤过平均（最可能模型） |"
          % (os.path.basename(bead.get("bead_model", "-")), bead.get("bead_model_beads", "-")))
    if bf.get("nsd") is not None:
        table = bf.get("nsd_by_selection") or {}
        A("| **CIFSUP NSD** | **%s**（%s；method=NSD） | 官方刻度只有一条：→0 理想，>1 表示系统性不同；**必须连口径一起报** |"
          % (f(bf["nsd"], 3), bf.get("selection")))
        if table:
            A("| NSD 的其他口径 | %s | 同一个模型换 selection 能差一个量级，别跨口径比较 |"
              % "，".join("%s %s" % (k, f(v, 3)) for k, v in table.items()))
    if df.get("best_correlation") is not None:
        A("| **fitmap** | **correlation %s**，cam %s，overlap %s（%s Å，%s metric，search %s seed %s） | 相关值随 R 变，别单点引用；overlap 与标度有关 |"
          % (f(df["best_correlation"], 4), f(df.get("best_correlation_about_mean"), 4),
             f(df.get("overlap")), f(df.get("best_resolution"), 1), df.get("metric"),
             df.get("search"), df.get("seed")))
        fits = df.get("fits") or []
        if fits:
            A("| R 阶梯 | %s | 相关值平不平不重要，**姿态**稳不稳才重要 |"
              % "；".join("%s Å → %s（%s 个唯一解）" % (f(x["resolution"], 0), f(x.get("correlation"), 4), x.get("n_unique_fits"))
                          for x in fits if x.get("correlation") is not None))
        if df.get("pose_rmsd_vs_best_A"):
            A("| 姿态散布 | %s（最大 %s Å） | >5 Å 且有分数接近的解 → 不能只报一个摆法 |"
              % ("，".join("%s: %s Å" % (k, f(v, 1)) for k, v in df["pose_rmsd_vs_best_A"].items()),
                 f(df.get("pose_spread_max_A"), 1)))
        if df.get("ambiguity"):
            A("")
            A("> **歧义提示**：%s" % df["ambiguity"])
    A("")

    # ---- files
    A("## 2. 这个目录里有什么")
    A("")
    A("- `embed_results.json` —— 本 README 里所有数字的机器可读版（含 `warning` / `ambiguity`）")
    if bf.get("aligned_model"):
        A("- `%s`（同内容 `.cif` 在旁）—— CIFSUP 叠好的模型（**在珠模型坐标系里**）；"
          "NSD 写在 cif 的 `score` 行 / PDB 的 REMARK 265" % os.path.basename(bf["aligned_model"]))
    if os.path.isdir(os.path.join(d.get("out_dir", "."), "bead")):
        A("- `bead/` —— 代建珠模型的全部中间件：GNOM `.out`、DAMMIF ×N 模型与日志、DAMAVER 共识模型与 `-distances.txt`（两两 NSD）")
    if df.get("fits"):
        A("- `fitmap_r<RR>.csv` / `.log` —— 每个 resolution 下 fitmap 的**全部唯一解**：旋转/平移矩阵 + correlation / cam / overlap / average_map_value / steps / shift / angle")
        A("- `*_in_density_r<RR>.pdb` —— 该 resolution 下最好的那个摆法（**在电子云坐标系里**）")
    figs = d.get("figures") or {}
    if not isinstance(figs, dict):
        figs = {}

    def _path(x):
        # figures.<panel> is {png, pse, render_log} (new) or a path string (old runs)
        return os.path.basename(x.get("png", "") if isinstance(x, dict) else (x or ""))

    beads_png = _path(figs.get("beads") or figs.get("bead")) or "embed_beads.png"
    dens_png = _path(figs.get("density")) or "embed_density.png"
    A("- `%s`（珠模型坐标系）/ `%s`（电子云坐标系）—— 两张图**故意分开**：两个靶子各自居中在各自原点，相对取向没有数据约束，混画会得到\"模型戳出珠球\"的假象"
      % (beads_png, dens_png))
    pses = [(k, os.path.basename(v["pse"])) for k, v in figs.items()
            if isinstance(v, dict) and v.get("pse")]
    if pses:
        A("- **`%s` —— 上面两张图的 PyMOL 会话（每张图一个）：样品 = 黄色 cartoon，包络 = 半透明白色**；"
          "直接用 `pymol <文件>.pse` 打开即可改视角/换 representation，成品图不必重渲"
          % "` / `".join(n for _, n in pses))
    for k, v in figs.items():
        if isinstance(v, dict) and v.get("status"):
            A("- ⚠️ 面板 %s 没出图：`%s`" % (k, v["status"]))
    rl = [l for v in figs.values() if isinstance(v, dict) for l in (v.get("render_log") or [])]
    if rl:
        A("- 渲染参数（density 的等值面水平与包络体积、beads 的球半径）：%s"
          % "；".join("`%s`" % l.strip() for l in rl))
    A("")
    A("## 3. 判据（本目录的数字该怎么读）")
    A("")
    A("| 数字 | 好 | 要小心 | 不可用 |")
    A("|---|---|---|---|")
    A("| CIFSUP NSD（method=NSD） | → 0（理想叠合） | 0.7–1.0：能摆进去，但珠模型/系综会影响第 2 位 | **> 1：两个物体系统性地不同**（官方刻度） |")
    A("| fitmap correlation | 越接近 1 越像 | 0.85–0.92（本机同类样品实测档）——**相关值本身不判对错** | — |")
    A("| 姿态散布（各 R 的最好摆法之间的 CA-RMSD） | < 5 Å：摆法基本唯一 | > 5 Å：存在多个分数相当的摆法 | — |")
    A("| 唯一解个数 | 少 | 多（同一图里 9–23 个）说明包络对摆放的约束弱 | — |")
    A("")
    A("**本 skill 的分数不是拟合优度**：NSD / correlation 只说\"形状能不能摆进去、摆得贴不贴\"。")
    A("模型配不配这条曲线要另做 CRYSOL/PDB2SAS 直接拟合（见 `fit-a-high-resolution-model-to-data`）。")
    A("")
    A("## 4. 这份结果**没有**做的事")
    A("")
    A("- 没有判模型配不配数据（无 χ² 拟合优度；CRYSOL 是另一条路）")
    A("- 没有做任何精修/柔性拟合（只有刚体叠合）")
    A("- 没有测电子云的分辨率（单模型无 FSC；用 R 阶梯 + 姿态散布代替）")
    A("- 没有在珠模型与电子云之间建立相对取向（两者坐标系独立）")
    A("")
    A("## 5. 怎么复核 / 重跑")
    A("")
    A("```bash")
    A("# 只重做嵌入（珠模型已在 embed/bead/ 里时不会重建）")
    A("/Applications/BioXTASRAW/bin/python embed-model.py %s \\" % (inp.get("sample_dir") or "<sample-dir>"))
    A("    --model %s --out-dir %s" % (d.get("model", "<model.pdb>"), d.get("out_dir", "<out>")))
    A("")
    A("# 单独复核 CIFSUP 的分数")
    A("grep -m1 '^score' %s" % os.path.basename(bf.get("aligned_model_cif") or "xxx_in_beads.cif"))
    A("")
    A("# 单独看 fitmap 的全部解（第 25 列 = correlation，按它降序）")
    A("sort -k25 -g -r %s | head" % (os.path.basename((df.get("fits") or [{}])[0].get("csv") or "fitmap_r15.csv")))
    A("```")
    A("")
    return "\n".join(L)


def main():
    for d in sys.argv[1:]:
        d = os.path.abspath(os.path.expanduser(d))
        j = os.path.join(d, "embed_results.json")
        if not os.path.exists(j):
            jj = sorted(__import__("glob").glob(os.path.join(d, "*", "embed_results.json")))
            for x in jj:
                _one(x)
            continue
        _one(j)
    return 0


def _one(j):
    d = json.load(open(j))
    d.setdefault("out_dir", os.path.dirname(j))
    txt = build(d)
    out = os.path.join(os.path.dirname(j), "README.md")
    open(out, "w").write(txt)
    print("wrote", out)


if __name__ == "__main__":
    sys.exit(main())
