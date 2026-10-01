#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把一次 SEC-SAXS 处理的结果目录写成一份"给人看"的 README.md。

设计原则（照着用户的原话来的："我看不懂你现在的结果"）：
  · 只用**产物里现成的东西**（run_meta.json / tables/*.csv / series/ranges.json / 目录清单），
    不重新拟合、不臆造数字——README 里的每个数都能在某个文件里找到原文。
  · 结论先行：开头就写"这次算出来什么、可信到什么程度、下一步看哪三个文件"。
  · 每个目录逐条说"是什么 / 怎么打开 / 是不是最终产物"。
  · 把"判读红线"和"本次告警"显式写出来（rg=-1 是失败哨兵、χ² 多大算差、NSD 多大算一致、
    没 ATSAS 所以没珠模、没监视器所以没归一化……）。

用法（独立跑，可用在已有结果目录上）：
  python write-results-readme.py --out <结果目录> [--language zh]
管线内会被 run-raw-sec-pipeline.py 在最后自动调用。
"""
from __future__ import annotations

import csv
import glob
import json
import os

SENTINEL = -1.0          # RAW 的失败哨兵值：rg / i0 / mw 给 -1.0 表示该区间没收敛


# ----------------------------------------------------------------- 小工具
def _fmt(x, nd=2, dash="—"):
    """数字格式化：失败哨兵、NaN、None 一律显示成破折号或"未收敛"。"""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return dash
    if v != v or v is None:                       # NaN
        return dash
    if abs(v + 1.0) < 1e-9:                       # RAW 哨兵
        return "未收敛"
    if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e5):
        return f"{v:.2e}"
    return f"{v:.{nd}f}"


def _read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="", errors="ignore") as fh:
        return list(csv.DictReader(fh))


def _read_json(path, default=None):
    if not os.path.exists(path):
        return default
    try:
        with open(path, errors="ignore") as fh:
            return json.load(fh)
    except Exception:
        return default


def _count(path, pattern="*"):
    return len(glob.glob(os.path.join(path, pattern)))


def _parse_damaver(models_dir):
    """DAMAVER 自己的口径写在 *-distances.txt 里（"Mean value of nsd"），比 API 返回值更可信。"""
    hits = glob.glob(os.path.join(models_dir, "*-distances.txt"))
    if not hits:
        return None
    txt = open(sorted(hits)[0], errors="ignore").read()
    out = {}
    for line in txt.splitlines():
        low = line.lower()
        if "mean value of nsd" in low:
            out["mean"] = line.split(":")[-1].strip()
        elif "standard deviation of nsd" in low:
            out["std"] = line.split(":")[-1].strip()
    return out or None


def _tree_lines(out, depth=1):
    lines = []
    for root, dirs, files in os.walk(out):
        dirs[:] = sorted(d for d in dirs if not d.startswith("."))
        rel = os.path.relpath(root, out)
        if rel == ".":
            continue
        level = rel.count(os.sep) + 1
        if level > depth:
            continue
        n = len(files)
        ex = sorted(os.path.splitext(f)[1] for f in files if not f.startswith("."))
        kinds = ", ".join(sorted({e for e in ex if e})) or "—"
        lines.append(f"{rel}/  ({n} 个文件{'；'+kinds if kinds else ''})")
    return lines


# ----------------------------------------------------------------- 数据汇总
def collect(out):
    meta = _read_json(os.path.join(out, "run_meta.json"), {}) or {}
    args = meta.get("args", {}) or {}
    ranges = _read_json(os.path.join(out, "series", "ranges.json"), {}) or {}
    d = {
        "out": out,
        "meta": meta,
        "args": args,
        "ranges": ranges,
        "guinier": _read_csv(os.path.join(out, "tables", "guinier_multi_range.csv")),
        "ift": _read_csv(os.path.join(out, "tables", "ift_summary.csv")),
        "mw": _read_csv(os.path.join(out, "tables", "mw.csv")),
        "frames_params": _read_csv(os.path.join(out, "tables", "frame_params.csv")),
        "n_integrated": _count(os.path.join(out, "profiles", "01_integrated")),
        "n_subtracted": _count(os.path.join(out, "profiles", "03_subtracted")),
        "models": sorted(os.path.basename(p) for p in glob.glob(os.path.join(out, "models", "*"))),
        "video": sorted(os.path.basename(p) for p in glob.glob(os.path.join(out, "video", "*"))),
        "ifts": sorted(os.path.basename(p) for p in glob.glob(os.path.join(out, "ifts", "*"))),
        "reports": sorted(os.path.basename(p) for p in glob.glob(os.path.join(out, "reports", "*"))),
        "norm": sorted(os.path.basename(p) for p in glob.glob(os.path.join(out, "norm", "*"))),
        "tables": sorted(os.path.basename(p) for p in glob.glob(os.path.join(out, "tables", "*"))),
    }
    d["prefix"] = meta.get("prefix") or os.path.basename(os.path.normpath(out))
    d["denss_rows"] = meta.get("denss_rows", [])
    d["dammif_rows"] = meta.get("dammif_rows", [])
    d["damaver_row"] = meta.get("damaver_row", [])
    d["damaver_txt"] = _parse_damaver(os.path.join(out, "models"))
    return d


def _guinier_table(rows):
    if not rows:
        return ["（没有 tables/guinier_multi_range.csv）"]
    out = ["| 区间 | Rg (Å) | Rg 误差 | I(0) | qRg 范围 | r² | 判读 |",
           "|---|---|---|---|---|---|---|"]
    for r in rows:
        rg = float(r.get("rg", "nan") or "nan")
        rsq = float(r.get("r_sqr", "nan") or "nan")
        if rg < 0 or rsq < 0:
            verdict = "**不收敛**（哨兵值，别用）"
        elif rsq >= 0.99:
            verdict = "好"
        elif rsq >= 0.95:
            verdict = "可用"
        else:
            verdict = "差（低 q 被污染或区间选得不合适）"
        out.append(
            f"| {r.get('range_label','?')} | {_fmt(rg)} | {_fmt(r.get('rg_err'))} | {_fmt(r.get('i0'))} | "
            f"{_fmt(r.get('qRg_min'),2)}–{_fmt(r.get('qRg_max'),2)} | {_fmt(rsq,5)} | {verdict} |")
    return out


def _ift_table(rows):
    if not rows:
        return ["（没有 tables/ift_summary.csv）"]
    out = ["| 方法 | Dmax (Å) | Rg (Å) | χ² | log α | 说明 |", "|---|---|---|---|---|---|"]
    note = {"BIFT": "贝叶斯 IFT（RAW 原生），Dmax 自己扫出来的",
            "GNOM": "ATSAS GNOM（IFTM），形状重建的输入",
            "DIFT": "RAW 原生 DIFT，Dmax 是显式给的那个（DENSS 的输入）"}
    for r in rows:
        m = r.get("method", "?")
        out.append(f"| {m} | {_fmt(r.get('dmax'))} | {_fmt(r.get('rg'))} | {_fmt(r.get('chi_sq'),2)} | "
                   f"{_fmt(r.get('log_alpha'),2)} | {note.get(m,'')} |")
    return out


def _mw_table(rows):
    if not rows:
        return ["（没有 tables/mw.csv）"]
    out = ["| 方法 | 分子量 (kDa) | 备注 |", "|---|---|---|"]
    extra = {"Vc": "按体积不变假设（Vc，密堆效率 0.3）",
             "Vp": "按 Porod 体积（Vp）",
             "Bayesian": "ATSAS DATMW 的贝叶斯估计（需 ATSAS）",
             "Datclass": "ATSAS DATCLASS 的分类 + 分子量（需 ATSAS）"}
    for r in rows:
        out.append(f"| {r.get('method','?')} | {_fmt(r.get('mw'),1)} | {extra.get(r.get('method'),'')} |")
    return out


def _shape_block(d):
    out = []
    # --- DENSS（电子云）---
    out.append("**电子云重建（RAW 原生 DENSS，输出 `models/*_denss.mrc`）**")
    if d["denss_rows"]:
        out += ["", "| χ² | Rg (Å) | support 体积 (Å³) | 盒子边长 (Å) |", "|---|---|---|---|"]
        for r in d["denss_rows"]:
            out.append(f"| {_fmt(r[1],2)} | {_fmt(r[2],1)} | {_fmt(r[3],0)} | {_fmt(r[4],1)} |")
    elif any("_denss" in n and n.endswith(".mrc") for n in d["models"]):
        out.append("")
        out.append("（在 models/ 里，但本次 run_meta.json 没记数——看 `models/*_denss.log` 尾部的 "
                   "`Chi^2`/`Rg`/`Support` 行）")
    else:
        out.append("")
        out.append("**没有出电子云**（`models/` 里无 `*_denss.mrc`）——要么没跑 `shape` 步，要么 DENSS 失败；"
                   "失败原因看 `models/*_denss.log`（`.mrc` 只写出一部分时也会留 `_current.mrc`）。")
    # --- DAMMIF（珠模）---
    out += ["", "**珠模（ATSAS DAMMIF，dummy-atom，输出 `models/*_dammif_*.cif`）**"]
    if d["dammif_rows"]:
        out += ["", "| 模型 | χ² | Rg (Å) | Dmax (Å) | MW (kDa) |", "|---|---|---|---|---|"]
        for i, r in enumerate(d["dammif_rows"], 1):
            out.append(f"| #{i} | {_fmt(r[0],2)} | {_fmt(r[1],1)} | {_fmt(r[2],1)} | {_fmt(r[3],0)} |")
    elif not any("_dammif_" in n and n.endswith((".cif", ".pdb")) for n in d["models"]):
        out.append("")
        out.append("**没有出珠模**（`models/` 里无 `*_dammif_*.cif`）：DAMMIF 是 ATSAS 可执行文件的外壳，"
                   "要么没装 ATSAS / 没给 `--atsas-dir`，要么这一步失败了（看日志里的 `DAMMIF #n 失败`）。")
    elif not any("_dammif_" in n and n.endswith(".png") for n in d["models"]):
        out.append("")
        out.append(f"（本次 run_meta.json 没记每个模型的 χ²/Rg/Dmax —— 见 `models/` 下 "
                   f"`*_dammif_*.log` 与 `*_dammif_*.fit`）")
    rep = None
    for n in d["models"]:
        if n.endswith("-summary.txt"):
            try:
                for ln in open(os.path.join(d["out"], "models", n), errors="ignore"):
                    if "Most representative" in ln:
                        rep = ln.split(":")[-1].strip()
            except Exception:
                pass
    if d["damaver_txt"]:
        m, sd = d["damaver_txt"].get("mean", "?"), d["damaver_txt"].get("std", "?")
        out += ["", f"**DAMAVER 一致性**：平均 NSD = **{m}** ± {sd}"
                    "（这是 DAMAVER 自己在 `*-distances.txt` 里给的口径：模型两两差异，越小越一致；"
                    "同一类形状一般 ≲2，明显更大或分成多个 cluster = 形状还没定）"]
        if rep:
            out += [f"**代表模型**：`models/{rep}`（在 `*-global-summary.txt` 里被标为 most representative）"]
    elif d["damaver_row"]:
        out += ["", f"**DAMAVER**：平均 NSD = {_fmt(d['damaver_row'][0],3)}"
                    f" ± {_fmt(d['damaver_row'][1] if len(d['damaver_row'])>1 else None,3)}"]
    return out


def _warnings(d):
    w = []
    st = d["meta"].get("settings", {}) or {}
    if not st.get("EnableNormalization", False):
        w.append("**本次没有做逐帧归一化**（`EnableNormalization=False`）：该系列没有线站逐帧 txt / 监视器，"
                 "所以束流起伏会留在数据里——低 q 与绝对刻度都别太当真。")
    if not d["models"]:
        w.append("**`models/` 是空的**：形状重建没跑（`--steps` 里没有 `shape`，或整步跳过）。")
    for r in d["guinier"]:
        if float(r.get("rg", 0) or 0) < 0:
            w.append(f"Guinier 区间 `{r.get('range_label')}` 未收敛（RAW 用 `rg=-1` 当失败哨兵，不是真的负值）"
                     "——通常是 q_min 太低、低 q 上翘把拟合带跑。")
    rg = None
    for r in d["guinier"]:
        if r.get("range_label") == "auto":
            rg = float(r.get("rg", "nan") or "nan")
    for r in d["ift"]:
        try:
            dd = float(r.get("dmax"))
        except (TypeError, ValueError):
            continue
        if dd != dd or rg is None or rg != rg or rg <= 0:
            continue
        if dd / rg <= 4:
            continue
        others = [f"{x.get('method')} {float(x.get('dmax')):.0f} Å" for x in d["ift"]
                  if x is not r and x.get("dmax") not in (None, "", "nan")
                  and float(x.get("dmax") or 0) > 0 and float(x.get("dmax")) / rg <= 4]
        w.append(f"**{r.get('method')} 的 Dmax（{dd:.0f} Å）相对 Guinier Rg（{rg:.1f} Å）大了 {dd/rg:.1f} 倍**"
                 "（球状蛋白一般 2.5–3 倍）：这一支的 P(r)/重建被低 q 拖长，只能当示意。"
                 + (f" 同一份数据里 **{'、'.join(others)}** 才在正常量级（形状重建用后者）。" if others else
                    " 先回去看低 q 是否该裁（`--trim-qmin`）或显式给 Dmax（`--ift-dmax`）。"))
    if d["dammif_rows"]:
        chis = [float(r[0]) for r in d["dammif_rows"] if r]
        if chis and max(chis) > 5:
            w.append(f"**珠模 χ² 偏大（最高 {max(chis):.1f}）**：这不是重建参数问题，是**曲线本身还不够干净**"
                     "（弱峰 / 未归一化 / 低 q 污染）。χ² 单位是拟合残差，接近 1–3 才算像样。")
    return w


# ----------------------------------------------------------------- 主函数
def write_results_readme(out_dir, prefix=None, language="zh"):
    d = collect(out_dir)
    prefix = prefix or d["prefix"]
    meta, args, ranges = d["meta"], d["args"], d["ranges"]
    st = meta.get("settings", {}) or {}
    rng_b = ranges.get("buffer") or []
    rng_s = ranges.get("sample") or []
    fmt_rng = lambda rs: "、".join(f"第 {a}–{b} 帧" for a, b in rs) if rs else "（未记录）"

    rg_auto = None
    for r in d["guinier"]:
        if r.get("range_label") == "auto":
            rg_auto = float(r.get("rg", "nan") or "nan")
    mw_main = None
    for r in d["mw"]:
        if r.get("method") == "Vc":
            mw_main = float(r.get("mw", "nan") or "nan")

    lines = []
    A = lines.append
    A(f"# {prefix} —— SEC-SAXS 处理结果")
    A("")
    A("> 这份 README 由 `run-a-sec-saxs-pipeline-end-to-end` 的 `write-results-readme.py` 从**产物本身**"
      "读出来（不重新拟合、不臆造数字）。生成时间取自 `run_meta.json`。")
    A("")
    A("## 0. 三十秒结论")
    A("")
    concl = []
    if rg_auto is not None and rg_auto > 0:
        concl.append(f"Rg ≈ **{rg_auto:.1f} Å**")
    if mw_main is not None and mw_main > 0:
        concl.append(f"分子量（Vc）≈ **{mw_main:.0f} kDa**")
    if d["dammif_rows"]:
        best = min(float(r[0]) for r in d["dammif_rows"] if r)
        concl.append(f"珠模 χ² ≈ **{best:.2f}**")
    if d["damaver_txt"]:
        concl.append(f"DAMAVER 平均 NSD **{d['damaver_txt'].get('mean','?')}**")
    elif d["damaver_row"]:
        concl.append(f"DAMAVER 平均 NSD **{_fmt(d['damaver_row'][0],3)}**")
    A("· " + "；".join(concl) if concl else "· （结果目录里还没有可汇总的关键数字）")
    A("")
    A(f"· 输入：`{meta.get('input','?')}`，共 **{meta.get('n_frames', d['n_integrated'])} 帧**"
      f"；扣减用 buffer {fmt_rng(rng_b)}，样品区 {fmt_rng(rng_s)}。")
    A(f"· 处理：全程 RAW（`{meta.get('cfg','?')}` 的设置档 + 命令行参数），"
      f"运行时间戳 `{meta.get('timestamp','?')}`。")
    A("")
    A("**只看三个文件就够**：")
    A("")
    A("| 想回答 | 打开 |")
    A("|---|---|")
    A(f"| 整体对不对？ | `reports/{d['reports'][0] if d['reports'] else prefix+'_raw_report.pdf'}`（RAW 报告：Guinier/IFT/序列图都在里面） |")
    A("| 选帧选得对不对？ | `series/series_plot.png`（洗脱曲线 + 选区）+ `series/ranges.json` |")
    A("| 形状重建长什么样？ | `models/*_dammif_01.png`（拟合图）；平均模型 `models/*_damaver-global-damaver.cif` |")
    A("")
    A("## 1. 目录导航")
    A("")
    A("| 路径 | 里面是什么 | 怎么打开 / 注意 |")
    A("|---|---|---|")
    A(f"| `profiles/01_integrated/` | 每一帧的 I(q)（径向积分后的曲线，含逐帧归一化） "
      f"| 文本 .dat，3 列 q/I/err；{d['n_integrated']} 条 |")
    A(f"| `profiles/02_buffer/` | buffer（空池/基线）平均曲线 | `buffer_avg.dat` |")
    A(f"| `profiles/03_subtracted/` | 每一帧扣减 buffer 之后的 I(q) | 逐帧 .dat，{d['n_subtracted']} 条；系列图的数据源 |")
    A(f"| `profiles/04_sample/` | 样品区平均后的曲线（下游分析用的就是它） | `sample_avg*.dat`，多区间拟合基于它 |")
    A(f"| `profiles/06_guinier/` | 每个 qRg 区间的 Guinier 拟合（原始点+拟合线） | 文本 .dat，画图看直线段 |")
    A(f"| `tables/` | 汇总表（见第 2 节） | CSV，Excel 直接开 |")
    A(f"| `series/` | 序列对象与洗脱曲线 | `.hdf5`（RAW 打开）、`series_plot.png`、`ranges.json` |")
    A(f"| `ifts/` | P(r) 反变换结果 | `*_bift.ift`/`*_gnom.out`/`*_denss_ift.ift`（RAW 打开） |")
    A(f"| `models/` | 3D 重建：电子云(.mrc) + 珠模(.cif) + DAMAVER 平均 | `.mrc` 用 ChimeraX/PyMOL 开；`.cif` 同上；日志同目录 |")
    A(f"| `reports/` | RAW 生成的 PDF 报告 | 双击打开 |")
    A(f"| `video/` | 裁剪区（束斑附近）的归一化逐帧视频 | .mp4；`*.centroid.csv` 是质心/强度轨迹 |")
    A(f"| `norm/` | 逐帧归一化因子表（若做了归一化） | `normalization_factors.csv` |")
    A(f"| `run_meta.json` | 本次运行的参数与关键结果（机器可读） | 复现时照它抄参数 |")
    A("")
    A("## 2. 关键数字（都从产物里现读）")
    A("")
    A("### 2.1 Guinier 多区间（判断 Rg 稳不稳）")
    A("")
    lines += _guinier_table(d["guinier"])
    A("")
    A("> 多区间的意义：**每个区间都给一个 Rg**，如果它们彼此接近，说明 Rg 可信；如果随区间乱跳，"
      "说明低 q 被污染（上翘/聚集）或者浓度序列有问题。")
    A("")
    A("### 2.2 IFT / P(r)")
    A("")
    lines += _ift_table(d["ift"])
    A("")
    A("### 2.3 分子量")
    A("")
    lines += _mw_table(d["mw"])
    A("")
    A("### 2.4 形状重建")
    A("")
    lines += _shape_block(d)
    A("")
    A("## 3. 本次用的参数（复现用）")
    A("")
    A("```")
    A("run-raw-sec-pipeline.py \\")
    A(f"  --series-dir {meta.get('input','<帧目录>')} \\")
    A(f"  --out-dir {out_dir} \\")
    A(f"  --cfg {meta.get('cfg','<设置档>')} \\")
    if rng_b:
        A(f'  --buffer-range "{";".join(f"{a},{b}" for a,b in rng_b)}" \\')
    if rng_s:
        A(f"  --sample-range {rng_s[0][0]},{rng_s[0][1]} \\")
    if args.get("trim_qmin"):
        A(f"  --trim-qmin {args['trim_qmin']} \\")
    if args.get("ift_dmax"):
        A(f"  --ift-dmax {args['ift_dmax']} \\")
    if args.get("no_header_normalization"):
        A("  --no-header-normalization \\")
    if args.get("atsas_dir"):
        A(f"  --atsas-dir {args['atsas_dir']} \\")
    A(f"  --steps {','.join(meta.get('steps', []))} \\")
    A(f"  --model-engine {meta.get('model_engine','auto')} --n-models {meta.get('n_models',4)}")
    A("```")
    A("")
    A(f"* RAW 端设置：`ImageHdrFormat={st.get('ImageHdrFormat','None')}`，"
      f"`EnableNormalization={st.get('EnableNormalization',False)}`，"
      f"`NormalizationList={st.get('NormalizationList','—')}`；ATSAS：`{st.get('ATSASDir') or '未指定'}`。")
    A("")
    A("## 4. 怎么判读（红线）")
    A("")
    A("| 数字 | 好 | 差 | 说明 |")
    A("|---|---|---|---|")
    A("| Guinier r² | ≥0.99 很好；≥0.95 可用 | <0.9 或为负 | 负 r² / Rg=-1 是 RAW 的**失败哨兵**，不是物理值 |")
    A("| 多区间 Rg 一致性 | 各区间差 <10% | 随区间单调乱变 | 乱变=低 q 污染（上翘/聚集） |")
    A("| BIFT/GNOM 的 Dmax/Rg | 2.5–3 倍（球状） | >4 倍 | 大了说明 P(r) 被低 q 拖长，重建只能当示意 |")
    A("| 珠模 χ² | 1–3 | >5 | χ² 大是**数据**问题，先修曲线再谈参数 |")
    A("| DAMAVER 平均 NSD | ≲2（同一类形状） | 明显更大 / 分成多个 cluster | 模型分簇 = 形状还没定 |")
    A("| DENSS support 体积 | 与 MW 自洽 | 明显过大 | 盒子边长 ≈ oversampling×Dmax |")
    A("")
    A("## 5. 本次的告警 / 未完成项")
    A("")
    ws = _warnings(d)
    if ws:
        for x in ws:
            A(f"- {x}")
    else:
        A("- 无（本流程能自动判的项都过了）")
    A("")
    A("## 6. 这些缩写")
    A("")
    A("- **Rg**：回转半径（Å）。**I(0)**：零角散射强度。**qRg**：无量纲化 q，Guinier 区一般要求 qRg≲1.3。")
    A("- **IFT/P(r)**：把 I(q) 反变换成实空间距离分布，Dmax 是它的支撑上限。")
    A("- **Vc/Vp/Bayesian/Datclass**：四种独立的分子量估计（前两个只需曲线，后两个用 ATSAS）。")
    A("- **DENSS**：电子云（密度图）重建；**DAMMIF**：dummy-atom 珠模；**DAMAVER**：把多个珠模对齐平均、给一致性 NSD。")
    A("- **χ²**：重建模型与实验曲线的拟合残差（≈1 表示在误差范围内）。")
    A("")
    A("---")
    A("")
    A("### 附：目录树（一层）")
    A("")
    A("```")
    for t in _tree_lines(out_dir, depth=2):
        A(t)
    A("```")
    A("")

    path = os.path.join(out_dir, "README.md")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    return path


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="把 SEC-SAXS 结果目录写成 README.md")
    ap.add_argument("--out", required=True, help="结果目录（含 run_meta.json）")
    ap.add_argument("--prefix", default=None, help="覆盖前缀（默认取 run_meta.json / 目录名）")
    a = ap.parse_args()
    print("→", write_results_readme(a.out, a.prefix))
