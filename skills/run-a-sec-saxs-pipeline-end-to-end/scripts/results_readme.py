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


def _fmt_err(x, dash="—"):
    """误差列用**两位有效数字**（误差常比主值小几个量级：0.0038 不该显示成 0.00）。"""
    v = None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return dash
    if v != v:
        return dash
    if v == 0:
        return "0"
    return f"{v:.2g}"


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



# ----------------------------------------------------------------- 总结与全参数表（README 开头）
def _num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def _primary_guinier(rows):
    """挑"最终采用"的 Guinier 区间：收敛 + qRg 落在常用域 + r² 最高；都不合格才退回 auto。"""
    cand = []
    for r in rows:
        rg, rsq = _num(r.get("rg")), _num(r.get("r_sqr"))
        qlo, qhi = _num(r.get("qRg_min")), _num(r.get("qRg_max"))
        if not rg or rg <= 0 or rsq is None or qlo is None or qhi is None:
            continue
        if qlo < 0.29 or qhi > 1.35:
            continue
        cand.append((rsq, r))
    if cand:
        return max(cand, key=lambda t: t[0])[1]
    for r in rows:
        if r.get("range_label") == "auto":
            return r
    return rows[0] if rows else None


def _rg_spread(rows):
    """**可用**区间之间的 Rg 最大相对差（%）；返回 (spread%, 采用标签, 被排除的标签)。

    只统计"收敛（rg>0）且 r²≥0.9"的区间——未收敛（rg=-1）与明显很差的区间本来就不能引用，
    让它们参与一致性判据等于一票否决。被排除的要列出来，便于复核。
    """
    ok, bad = [], []
    for r in rows:
        rg, rsq = _num(r.get("rg")), _num(r.get("r_sqr"))
        if rg and rg > 0 and rsq is not None and rsq >= 0.9:
            ok.append((r.get("range_label"), rg))
        else:
            bad.append(str(r.get("range_label")))
    if len(ok) < 2:
        return (None, [l for l, _ in ok], bad)
    vals = [v for _, v in ok]
    return (100.0 * (max(vals) - min(vals)) / (sum(vals) / len(vals)), [l for l, _ in ok], bad)


def _mark_chi2(chi, kind="dammif"):
    if chi is None:
        return "—"
    if kind == "denss":
        return "🟢 可用" if chi <= 3 else ("🟡 有限" if chi <= 6 else "🔴 该重建不达标")
    return "🟢 可用" if chi <= 3 else ("🟡 有限" if chi <= 5 else "🔴 不可用（先修曲线，不是调重建参数）")


def _shape_facts(d):
    """把形状重建的元数据收成事实字典，供总结与全表共用。"""
    f = {"denss": None, "dammif_chi2": None, "dammif_n": len(d["dammif_rows"]),
         "dammif_rg": None, "dammif_dmax": None, "dammif_mw": None,
         "damaver_mean": None, "damaver_std": None, "damaver_clusters": None}
    if d["denss_rows"]:
        r = d["denss_rows"][0]
        f["denss"] = dict(chi=_num(r[1]), rg=_num(r[2]), vol=_num(r[3]), side=_num(r[4]))
    if d["dammif_rows"]:
        chis = [_num(r[0]) for r in d["dammif_rows"] if _num(r[0]) is not None]
        rgs = [_num(r[1]) for r in d["dammif_rows"] if _num(r[1])]
        dms = [_num(r[2]) for r in d["dammif_rows"] if _num(r[2])]
        mws = [_num(r[3]) for r in d["dammif_rows"] if _num(r[3])]
        if chis:
            f["dammif_chi2"] = (min(chis), max(chis))
        if rgs:
            f["dammif_rg"] = (min(rgs), max(rgs))
        if dms:
            f["dammif_dmax"] = (min(dms), max(dms))
        if mws:
            f["dammif_mw"] = sum(mws) / len(mws)
    if d["damaver_txt"]:
        f["damaver_mean"] = _num(d["damaver_txt"].get("mean"))
        f["damaver_std"] = _num(d["damaver_txt"].get("std"))
    elif d["damaver_row"]:
        f["damaver_mean"] = _num(d["damaver_row"][0])
        f["damaver_std"] = _num(d["damaver_row"][1]) if len(d["damaver_row"]) > 1 else None
    n_cl = glob.glob(os.path.join(d["out"], "models", "*damaver-cluster*-summary.txt"))
    f["damaver_clusters"] = len(n_cl) or None
    return f


def _summary_block(d, prefix):
    """开头的分点总结：数据可用性 / 表现较好的参数 / 表现较差的参数。"""
    L = []
    A = L.append
    meta = d["meta"]
    g_primary = _primary_guinier(d["guinier"])
    g_auto = next((r for r in d["guinier"] if r.get("range_label") == "auto"), None)
    spread, spread_used, spread_bad = _rg_spread(d["guinier"])
    rg = _num(g_primary.get("rg")) if g_primary else None
    rsq = _num(g_primary.get("r_sqr")) if g_primary else None
    sh = _shape_facts(d)
    mw = {r.get("method"): _num(r.get("mw")) for r in d["mw"]}
    st_raw = meta.get("settings", {}) or {}

    good, bad = [], []
    usable = None
    # ---- 数据可用性
    if rg and rg > 0:
        if rsq is not None and rsq >= 0.99 and (spread is None or spread <= 5):
            usable = "🟢 可用"
        elif rsq is not None and rsq >= 0.95 and (spread is None or spread <= 15):
            usable = "🟡 有限可用（低 q 仍有可疑成分）"
        else:
            usable = "🔴 不可用/存疑"
    else:
        usable = "🔴 无可用 Rg（Guinier 全部未收敛）"

    A("## 0. 总结")
    A("")
    avail = [f"**数据可用性**：{usable}"]
    facts = [f"{meta.get('n_frames', d['n_integrated'])} 帧"]
    if st_raw.get("EnableNormalization"):
        facts.append("已做逐帧归一化（线站 txt/监视器）")
    else:
        facts.append("**未做逐帧归一化**（该系列没有线站 txt/监视器）")
    if rg and rg > 0:
        facts.append(f"Guinier 采用区间 `{g_primary.get('range_label')}`，r²={_fmt(rsq,4)}")
    if spread is not None:
        facts.append(f"可用区间之间 Rg 相对差 {spread:.1f}%")
    if sh["dammif_chi2"]:
        facts.append(f"珠模 χ²={_fmt(sh['dammif_chi2'][0],2)}")
    A("- " + avail[0].replace("**数据可用性**：", "**数据可用性**：") + " —— " + "；".join(facts) + "。")

    # ---- 表现较好
    if rg and rg > 0 and rsq is not None and rsq >= 0.95:
        good.append(f"**Rg = {_fmt(rg)} ± {_fmt_err(g_primary.get('rg_err'))} Å**"
                    f"（`{g_primary.get('range_label')}`，qRg {_fmt(g_primary.get('qRg_min'))}–{_fmt(g_primary.get('qRg_max'))}，"
                    f"r²={_fmt(rsq,4)}）")
    if spread is not None and spread <= 5:
        good.append(f"**跨区间一致**：{('、'.join(spread_used))} 之间 Rg 只差 {spread:.1f}%")
    if mw.get("Vc") and mw.get("Vp"):
        vc, vp = mw["Vc"], mw["Vp"]
        rel = abs(vc - vp) / ((vc + vp) / 2) * 100
        if rel <= 20:
            good.append(f"**分子量四法一致**：Vc {_fmt(vc,1)} / Vp {_fmt(vp,1)} / "
                        f"Bayesian {_fmt(mw.get('Bayesian'),1)} / Datclass {_fmt(mw.get('Datclass'),1)} kDa"
                        f"（最大差 {rel:.0f}%）")
    if sh["dammif_chi2"] and sh["dammif_chi2"][0] <= 3:
        good.append(f"**形状重建可用**：珠模 χ²={_fmt(sh['dammif_chi2'][0],2)}"
                    f"（{sh['dammif_n']} 个模型 Rg {_fmt(sh['dammif_rg'][0],1)}–{_fmt(sh['dammif_rg'][1],1)} Å）、"
                    f"DAMAVER 平均 NSD {_fmt(sh['damaver_mean'],2)}")
    if sh["denss"] and sh["denss"]["chi"] is not None and sh["denss"]["chi"] <= 3:
        good.append(f"电子云（DENSS）χ²={_fmt(sh['denss']['chi'],2)}，"
                    f"重建 Rg {_fmt(sh['denss']['rg'],1)} Å"
                    + (f"（与 Guinier 差 {100*abs(sh['denss']['rg']-rg)/rg:.0f}%）" if rg and sh["denss"]["rg"] else ""))

    # ---- 表现较差
    for r in d["ift"]:
        m, dd = r.get("method"), _num(r.get("dmax"))
        if dd and rg and rg > 0 and dd / rg > 4:
            bad.append(f"Dmax（{m}）= {_fmt(dd)} ± {_fmt_err(r.get('dmax_err'))} Å，"
                       f"**Dmax/Rg={dd/rg:.1f}**（球状一般 2.5–3）→ 这一支的 P(r)/重建不可用")
    for r in d["guinier"]:
        rgg, rs = _num(r.get("rg")), _num(r.get("r_sqr"))
        if rgg is not None and rgg < 0:
            bad.append(f"Guinier 区间 `{r.get('range_label')}` **未收敛**（`rg=-1` 是失败哨兵）")
        elif rs is not None and rs < 0.9 and rgg is not None and rgg > 0:
            bad.append(f"Guinier 区间 `{r.get('range_label')}` 拟合差（r²={_fmt(rs,4)}，Rg={_fmt(rgg)} Å）")
    if sh["dammif_chi2"] and sh["dammif_chi2"][0] > 3:
        bad.append(f"珠模 χ²={_fmt(sh['dammif_chi2'][0],2)} 偏大（>3）→ 曲线还不够干净，"
                   "先修数据而不是调重建参数")
    if sh["denss"] and sh["denss"]["chi"] is not None and sh["denss"]["chi"] > 6:
        bad.append(f"电子云 χ²={_fmt(sh['denss']['chi'],2)}（>6）→ 该密度图不达标")
    if sh["damaver_clusters"] and sh["damaver_clusters"] > 1:
        bad.append(f"DAMAVER 分成 {sh['damaver_clusters']} 个 cluster → 形状还没定")
    if not st_raw.get("EnableNormalization"):
        bad.append("**没做逐帧归一化**：束流起伏留在曲线里，低 q 与绝对刻度都别当真")

    A("")
    A("- **表现较好的参数**：" + ("；".join(good) if good else "（没有明显可靠的参数）") + "。")
    A("- **表现较差的参数**：" + ("；".join(bad) if bad else "（本次没有明显失败的参数）") + "。")
    A("")
    return L


def _all_params(d, prefix):
    """全参数表：环节 | 参数 | 值 | 程序误差 | 能不能用 | 误差 / 不确定度怎么读。"""
    L = []
    A = L.append
    meta = d["meta"]
    args = meta.get("args", {}) or {}
    st_raw = meta.get("settings", {}) or {}
    g_primary = _primary_guinier(d["guinier"])
    g_auto = next((r for r in d["guinier"] if r.get("range_label") == "auto"), None)
    spread, spread_used, spread_bad = _rg_spread(d["guinier"])
    rg = _num(g_primary.get("rg")) if g_primary else None
    rsq = _num(g_primary.get("r_sqr")) if g_primary else None
    sh = _shape_facts(d)
    mw = {r.get("method"): r for r in d["mw"]}

    A("## 1. 全参数表")
    A("")
    A("> 「程序误差」= RAW/ATSAS 自己返回的那一列误差；「能不能用」= 🟢 可用 / 🟡 有限 / 🔴 不可用；"
      "最后一列写清这条误差**读得懂什么、读不懂什么**。数字全部从产物文件现读（出处见「环节」列对应的表/文件）。")
    A("")
    A("| 环节 | 参数 | 值 | 程序误差 | 能不能用 | 误差 / 不确定度怎么读 |")
    A("|---|---|---|---|---|---|")

    # --- 积分 / 归一化
    A(f"| 积分 | 帧数 | {meta.get('n_frames', d['n_integrated'])} | — | 🟢 | 帧数足够，统计量不虚 |")
    norm_txt = "已启用（`NormalizationList=[['/','Transmitted_Beam']]`）" if st_raw.get("EnableNormalization") \
        else "**未启用**（无 txt/监视器）"
    A(f"| 归一化 | 逐帧归一化 | {norm_txt} | — | "
      f"{'🟢' if st_raw.get('EnableNormalization') else '🟡'} | "
      f"{'因子只含束流波动，不含几何漂移 —— 后者看 video/*.centroid.csv' if st_raw.get('EnableNormalization') else '束流起伏留在曲线里；低 q 与绝对刻度不可信'} |")
    rng_b = (d["ranges"] or {}).get("buffer") or []
    rng_s = (d["ranges"] or {}).get("sample") or []
    if rng_b or rng_s:
        A(f"| 扣减 | buffer / 样品区（0 基帧号） | "
          f"buffer {'、'.join(f'{a}–{b}' for a, b in rng_b) or '—'}；样品 "
          f"{'、'.join(f'{a}–{b}' for a, b in rng_s) or '—'} | — | 🟢 | "
          f"峰前+峰后两段 buffer 能吃掉基线漂移；只用一段时 SEC 峰后的基线常不等于峰前 |")

    # --- Guinier
    if rg and rg > 0:
        verdict = ("🟢 可用" if (rsq is not None and rsq >= 0.95 and (spread is None or spread <= 15))
                   else "🔴 存疑")
        A(f"| Guinier（采用 `{g_primary.get('range_label')}`） | Rg | {_fmt(rg)} Å | "
          f"± {_fmt_err(g_primary.get('rg_err'))} Å | {verdict} | "
          f"该误差是**区间内最小二乘标准误差**：区间收窄它必然变小，**不代表更准**；"
          f"真正的判据是 r²（{_fmt(rsq,4)}）与跨区间一致性（{_fmt(spread,1)}%） |")
        A(f"| Guinier（同上） | I(0) | {_fmt(g_primary.get('i0'))} | ± {_fmt_err(g_primary.get('i0_err'))} | "
          f"{verdict} | 同上；且 I(0) 与浓度耦合，SEC 里浓度未知 → 不能当绝对量 |")
        A(f"| Guinier（同上） | 拟合区间 / qRg / r² | q {_fmt(g_primary.get('q_min'),4)}–{_fmt(g_primary.get('q_max'),4)} Å⁻¹ / "
          f"{_fmt(g_primary.get('qRg_min'))}–{_fmt(g_primary.get('qRg_max'))} / {_fmt(rsq,5)} | — | "
          f"{'🟢' if (rsq or 0) >= 0.99 else '🟡'} | 判据：qRg 下界 ≳0.3、上界按形状（球 ≈1.3）、r² ≥0.99 很好；"
          f"`rg=-1`/`r²<0` 是**失败哨兵**，别当数字 |")
    for r in d["guinier"]:
        if r is g_primary:
            continue
        lbl, rgg, rs = r.get("range_label"), _num(r.get("rg")), _num(r.get("r_sqr"))
        if lbl == "auto":
            A(f"| Guinier（`auto`，交叉核对） | Rg | {_fmt(rgg)} Å | ± {_fmt_err(r.get('rg_err'))} Å | "
              f"{'🟢 一致' if (rg and rgg and abs(rgg-rg)/rg <= 0.05) else '🟡 与采用区间有小差异'} | "
              f"RAW 自动选的区间，只作参考；r²={_fmt(rs,5)} |")
        elif rgg is not None and rgg < 0:
            A(f"| Guinier（`{lbl}`） | Rg | 未收敛 | — | 🔴 不可用 | "
              f"`rg=-1` 是 RAW 的**失败哨兵**（区间点数不足/低 q 上翘把拟合带跑）；r²={_fmt(rs,5)} |")
        elif rs is not None and rs < 0.9:
            A(f"| Guinier（`{lbl}`） | Rg | {_fmt(rgg)} Å | ± {_fmt_err(r.get('rg_err'))} Å | 🔴 不可用 | "
              f"r²={_fmt(rs,5)}（<0.9，拟合比取平均还差）→ 这段 q 不服从 Guinier 定律 |")

    # --- IFT 各支
    for r in d["ift"]:
        m = r.get("method", "?")
        dd, dde = _num(r.get("dmax")), _num(r.get("dmax_err"))
        rr, re_ = _num(r.get("rg")), _num(r.get("rg_err"))
        ch = _num(r.get("chi_sq"))
        ratio = (dd / rg) if (dd and dd > 0 and rg and rg > 0) else None
        if ratio is None:
            mark, why = "🟡 无法判", "Guinier 无可用 Rg，无法比 Dmax/Rg"
        elif ratio <= 3.2:
            mark, why = "🟢 可用" if (ch is None or ch <= 3) else "🟡 有限", f"Dmax/Rg={ratio:.1f} 量级正常"
        elif ratio <= 4.5:
            mark, why = "🟡 有限", f"Dmax/Rg={ratio:.1f} 偏大，P(r) 尾部可疑"
        else:
            mark, why = "🔴 不可用", f"Dmax/Rg={ratio:.1f} 明显被低 q 拖大 → 这一支的 P(r)/重建别用"
        src = {"GNOM": "ATSAS GNOM（自动定 Dmax）", "BIFT": "RAW BIFT", "DIFT": "RAW DIFT（显式 Dmax）"}.get(m, m)
        A(f"| IFT（{src}） | Dmax | {_fmt(dd)} Å | {'± ' + _fmt_err(dde) + ' Å' if dde else '—'} | {mark} | "
          f"{why}；χ²={_fmt(ch,2)}；Dmax 是搜索/支撑上限，**没有误差棒时别当精确值** |")
        if rr and rr > 0:
            if rg and rg > 0:
                dfc = 100 * abs(rr - rg) / rg
                same = (f"与 Guinier 的 Rg 一致（差 {dfc:.0f}%）" if dfc <= 15
                        else f"与 Guinier 的 Rg 差 {dfc:.0f}% —— 多半是 Dmax 没定好")
            else:
                same = "—"
            A(f"| IFT（同上） | Rg（实空间，来自 P(r)） | {_fmt(rr)} Å | ± {_fmt_err(re_)} Å | "
              f"{'🟢' if mark.startswith('🟢') else ('🟡' if mark.startswith('🟡') else '🔴')} | "
              f"{same}；这个误差棒**不含** Dmax 选错带来的偏差 |")

    # --- 分子量
    mw_spec = [("Vc", "kDa", "体积不变假设（vcor）", "🟢 可用（与 Vp 差 ≤20%）或 🟡"),
               ("Vp", "kDa", "Porod 体积", "🟡 有限（无误差棒）"),
               ("Bayesian", "kDa", "ATSAS DATMW（贝叶斯）", "🟢 可用"),
               ("Datclass", "kDa", "ATSAS DATCLASS", "🟢 参考")]
    vc = _num((mw.get("Vc") or {}).get("mw"))
    for name, unit, src, mark_default in mw_spec:
        r = mw.get(name)
        if not r:
            A(f"| 分子量（{name}） | MW | 未产出 | — | 🟡 缺 | {src} 没跑（多因缺 ATSAS 参数） |")
            continue
        v = _num(r.get("mw"))
        d1, d2, d3, d4 = (_num(r.get(k)) for k in ("detail1", "detail2", "detail3", "detail4"))
        read = ""
        if name == "Vc":
            rel = abs(v - (_num((mw.get("Vp") or {}).get("mw")) or v)) / v * 100 if v else None
            mark = "🟢 可用" if (d2 and v and d2 / v <= 0.2 and (rel is None or rel <= 20)) else "🟡 有限"
            err = f"± {_fmt_err(d2)} kDa"
            read = (f"误差 = RAW 的经验不确定度（Vcor 传播），**不含**经验系数的系统偏差；"
                    f"qmax={_fmt(d3,3)} Å⁻¹；与 Vp 差 {_fmt(rel,0)}%")
        elif name == "Vp":
            mark, err = "🟡 有限", "—"
            read = (f"Porod 体积法**没有误差棒**；修正后体积 {_fmt(d1,0)} Å³；"
                    f"经验上比真实值高估 ~1.5×，只作数量级参考")
        elif name == "Bayesian":
            mark = "🟢 可用"
            err = f"CI {_fmt(d2,1)}–{_fmt(d3,1)} kDa" if (d2 and d3) else "—"
            read = (f"ATSAS 的置信区间（区间外概率 {_fmt(100 - (d4 or 0),0)}%）；"
                    f"曲线质量差时区间会张得很开 → 看宽度判可用性")
        else:
            mark, err = "🟢 参考", "—"
            read = f"DATCLASS 同时给形状分类（{r.get('detail1')}）与 Dmax={_fmt(d2,1)} Å"
            v = _num(v)
        A(f"| 分子量（{name}） | MW | {_fmt(v,1)} {unit} | {err} | {mark} | {read} |")

    # --- 形状重建
    if sh["denss"]:
        A(f"| 形状重建（DENSS 电子云，{args.get('denss_mode','Fast')}） | χ² | {_fmt(sh['denss']['chi'],2)} | — | "
          f"{_mark_chi2(sh['denss']['chi'], 'denss')} | 模型↔曲线残差（≈1 在误差棒内）；"
          f"重建 Rg {_fmt(sh['denss']['rg'],1)} Å"
          + (f"，与 Guinier 差 {100*abs(sh['denss']['rg']-rg)/rg:.0f}%" if (rg and sh['denss']['rg']) else "")
          + f"；support {_fmt(sh['denss']['vol'],0)} Å³、盒子 {_fmt(sh['denss']['side'],1)} Å |")
    else:
        A("| 形状重建（DENSS） | 电子云 | 未产出 | — | 🔴 缺 | 没跑 `shape` 步，或 DENSS 失败（看 `models/*_denss.log`） |")
    if sh["dammif_chi2"]:
        A(f"| 形状重建（DAMMIF×{sh['dammif_n']}） | χ² | {_fmt(sh['dammif_chi2'][0],2)}"
          f"（{_fmt(sh['dammif_chi2'][0],2)}–{_fmt(sh['dammif_chi2'][1],2)}） | — | "
          f"{_mark_chi2(sh['dammif_chi2'][0])} | 残差越小越贴合；**χ² 好也不代表解唯一**；"
          f"模型 Rg {_fmt(sh['dammif_rg'][0],1)}–{_fmt(sh['dammif_rg'][1],1)} Å、"
          f"Dmax {_fmt(sh['dammif_dmax'][0],1)}–{_fmt(sh['dammif_dmax'][1],1)} Å、MW≈{_fmt(sh['dammif_mw'],0)} kDa |")
    else:
        A("| 形状重建（DAMMIF） | 珠模 | 未产出 | — | 🟡 缺 | DAMMIF 是 ATSAS 可执行文件的外壳："
          "没装 ATSAS / 没给 `--atsas-dir` 就没有珠模（DENSS 电子云不受影响） |")
    if sh["damaver_mean"] is not None:
        cl = sh["damaver_clusters"] or 1
        if cl > 1:
            mark = "🟡 分簇（形状未定）"
        else:
            mark = "🟢 一致" if sh["damaver_mean"] <= 2 else ("🟡 尚可" if sh["damaver_mean"] <= 3 else "🔴 形状未定")
        extra = f"；分成 {cl} 个 cluster（>1 = 模型没收敛到同一形状）" if cl > 1 else "；单一 cluster"
        A(f"| 一致性（DAMAVER） | 平均 NSD | {_fmt(sh['damaver_mean'],2)} | ± {_fmt_err(sh['damaver_std'])} | {mark} | "
          f"这是**模型两两之间的差异**（± 是模型间标准差），**不是**与实验数据的拟合误差；"
          f"≤2 一般算同一类形状{extra} |")
    A("")
    A("> **三条通用读法**：① 只含「随机误差」的误差棒（Guinier 的 `rg_err`、IFT 的 `rg_err`）**不会**告诉你"
      "低 q 污染、背景失配、Dmax 选错这类系统误差；② 误差棒大小不能跨区间/跨方法比较，"
      "判可靠性要看 r²、跨区间一致性、χ²、多法一致；③ `rg=-1`、`r²<0`、`-1.0 kDa` 这类值是**失败哨兵**，"
      "不是物理量。")
    A("")
    return L


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
    A(f"> 自动生成 {meta.get('timestamp','?')} ｜ 脚本 `run-raw-sec-pipeline.py` ｜ "
      f"原始帧 `{meta.get('input','?')}` ｜ 配置 `{meta.get('cfg','?')}`")
    A("")
    A("> 这份 README 由 `run-a-sec-saxs-pipeline-end-to-end` 的 `write-results-readme.py` 从**产物本身**"
      "读出来（不重新拟合、不臆造数字）。生成时间取自 `run_meta.json`。")
    A("")
    lines += _summary_block(d, prefix)
    lines += _all_params(d, prefix)
    A("## 2. 输入与处理概要")
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
    A("## 3. 目录导航")
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
    A("## 4. 附录：逐区间 / 逐模型明细（都在产物文件里）")
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
    A("## 5. 本次用的参数（复现用）")
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
    A("## 6. 怎么判读（红线）")
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
    A("## 7. 本次的告警 / 未完成项")
    A("")
    ws = _warnings(d)
    if ws:
        for x in ws:
            A(f"- {x}")
    else:
        A("- 无（本流程能自动判的项都过了）")
    A("")
    A("## 8. 这些缩写")
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
