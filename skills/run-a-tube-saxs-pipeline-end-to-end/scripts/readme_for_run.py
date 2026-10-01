#!/usr/bin/env python
"""Write the human-readable README.md of ONE tube-pipeline result folder.

Reads <out>/summary.json (plus the files actually present) and writes <out>/README.md that
answers, in this order: 这批数据能不能用 → 关键数字 → 这个文件夹里每个文件是什么 →
每个数字的判据 → 这次没做什么/为什么 → 怎么自己复核。

Called automatically at the end of run-raw-tube-pipeline.py; can also be run alone on an
existing result folder to (re)generate the README without re-processing anything:

    python readme_for_run.py <out-dir> [<out-dir> ...]
"""
from __future__ import annotations

import argparse
import fnmatch
import numpy as np
import glob
import json
import os
import sys
import time

# (路径或通配, 中文说明, 什么时候看)  —— 顺序就是建议的阅读顺序
FILE_DOC = [
    ("qc.png", "四联诊断图：扣减曲线 / Kratky / 多区间 Rg 与 χ² / P(r)", "先看这张"),
    ("README.md", "本文件：结果说明与判据", "正在看"),
    ("reports/raw_report.pdf", "BioXTAS RAW 自己出的报告（Guinier 拟合、IFT、MW 的原始输出）", "要核对数字时"),
    ("profiles/03_subtracted/subtracted.dat", "**最终扣减曲线**：q、I(q)、误差三列，纯文本", "要拿去画图/喂别的软件时用这个"),
    ("profiles/02_sample/sample_avg.dat", "样品帧的平均曲线（未扣背景）", "想检查扣减前后的差别时"),
    ("profiles/01_control/control_avg.dat", "对照（buffer）帧的平均曲线", "同上"),
    ("profiles/01_control/control_avg_scaled.dat", "对照曲线按高 q 窗缩放后（真正参与扣减的那份）", "怀疑背景不匹配时"),
    ("profiles/04_guinier/recommended_range.dat", "被推荐的那段 Guinier 范围内的曲线", "复核 Rg 时"),
    ("profiles/04_guinier/guinier_*.dat", "每个候选 q 区间各一份曲线（文件名即 q 范围）", "想看区间怎么影响 Rg 时"),
    ("tables/guinier_multi_range.csv", "多区间 Guinier 全表：Rg/I0/误差/qRg/χ²_red/曲率/左右半段差 + 每条闸门", "判 Rg 信不信时（核心表）"),
    ("tables/guinier_results.json", "上表的机器可读版（含 auto 与 recommended）", "脚本用"),
    ("tables/ift_summary.csv", "各 IFT 引擎（BIFT 两次起跑 + GNOM）的 Dmax/Rg_real/chisq 与三条可信闸门", "判 P(r) 信不信时（核心表）"),
    ("tables/mw.csv", "分子量：Porod 体积法（Vp）与浓度无关的 Vc 法", "看分子量时"),
    ("tables/shape_results.json", "形状重建全表：电子云（DENSS）+ 每个珠模（DAMMIF）+ DAMAVER 一致性", "看 3D 时"),
    ("ifts/pr.dat", "**本次采用的 P(r)**（r、P(r)、误差；文件头写明来自哪个引擎）", "画 P(r) / 写论文时用这个"),
    ("ifts/pr_gnom.dat", "GNOM 的 P(r)（珠模就是用它算的）", "对比 GNOM 与 BIFT 的 P(r) 时"),
    ("ifts/bift.ift", "BIFT 的 IFT 解（RAW 原生，可被其它程序读）", "要复核 BIFT 解时"),
    ("ifts/gnom.out", "GNOM 的 IFT 解（ATSAS 格式，DAMMIF 的输入就是它）", "要重跑珠模时"),
    ("ifts/ift_fit.dat", "IFT 拟合曲线：实测 I(q) vs 拟合 I(q)", "看 IFT 拟合得好不好时"),
    ("ifts/*_untrusted_*", "**不可信**的 IFT 解（闸门没过，仅留档，不要引用）", "排查为什么没有 P(r) 时"),
    ("models/denss.mrc", "**电子云密度图**（DENSS，可用 ChimeraX/PyMOL 打开）", "要形状/内部空腔时"),
    ("models/denss_support.mrc", "DENSS 的支撑掩膜（模型边界）", "看模型大小是否合理时"),
    ("models/denss.log", "DENSS 自己的日志（每次迭代的 χ²/Rg）", "怀疑没收敛时"),
    ("models/dammif_*", "**珠模**（ATSAS DAMMIF，dummy-atom；`.cif` 直接拖进 ChimeraX/PyMOL）", "要形状、要投稿图时"),
    ("models/dammif_*.fit", "珠模对数据的拟合曲线与 χ²（RAW/PyMOL 都可看）", "判珠模拟合好不好时"),
    ("models/dammif_*.log", "DAMMIF 自己的日志（每次起跑的 χ²、珠子数）", "珠模 χ² 大时看这里"),
    ("models/damaver-*", "**多个珠模的一致性平均与代表模型**（`-global-damaver.cif` 就是代表模型）", "要给一个模型时用这个"),
    ("models/damaver-distances.txt", "模型两两差异（DAMAVER 自己的口径，里面的 Mean value of nsd 就是平均 NSD）", "判形状定没定时（核心）"),
    ("norm/frame_qc.csv", "逐帧质检：每帧低 q 电平、透射、与本 run 中位数的偏离、是否离群", "怀疑某一帧坏时"),
    ("<样品>_workspace.hdf5", "RAW 的 workspace 文件（File → Open Workspace 可打开，所有曲线都在里面）", "想在 GUI 里交互式复核时"),
    ("summary.json", "本 README 里所有数字的机器可读版", "脚本/批处理用"),
    ("frames/", "逐帧的 1D 曲线（只在用了 --save-frames 时才有）", "要逐帧看时"),
]

RULES_TEXT = """\
| 数字 | 好 | 勉强 | 差 / 不可信 |
|---|---|---|---|
| 高 q 窗对照缩放因子 | 1.00±0.02 | ±0.05 | ±0.05 以上：背景与样品不匹配，扣减后会有常数残留 |
| 低 q 对比度（扣减后 ÷ 背景电平） | >5% | 2–5% | <1%：信号埋在背景里，低 q 什么都不能说 |
| auto-Guinier R² | >0.99 | 0.95–0.99 | <0.95：这段不是 Guinier 区 |
| 多区间闸门 | 有区间全过 | 无区间全过但 Rg 漂移不大 | Rg 在区间间漂 >±20%：低 q 不服从单一 Guinier 定律 |
| IFT 三条闸门 | 全过（trusted） | 只差 Rg 那条（见 rg_tol） | 有 Dmax 跑到搜索域顶/超 4.5·Rg：**假解**，不许报 P(r) |
| Dmax / Rg | 2.5–4 | 4–4.5 | >4.5：P(r) 里多半是大颗粒尾巴 |
| Vc 分子量 | 与预期相符 | — | 明显偏大 >2×：多半是聚集体（Vp 法高估是常态，Vp≈1.5×Vc 属正常） |
| DENSS χ² | 越低越好 | — | 只在 IFT 可信时才有意义；模型 Rg 应与 P(r) 的 Rg 接近（差 <10%） |
| 珠模 χ²（DAMMIF） | 1–3 | 3–5 | >5：**是数据问题不是参数问题**——多半是低 q 污染（束挡光晕/背景残留），先修曲线再跑 |
| 珠模 Rg / Dmax | 与 P(r) 的差 <10% | 10–20% | 差很多：模型和数据说的不是同一个东西，别混着报 |
| DAMAVER 平均 NSD | ≲2（同一类形状） | — | 明显更大、或模型分成多个 cluster：**形状还没定**，加模型数再跑 |
"""


def _fmt_sub_contrast(v):
    """pipeline 的 contrast_lowq = I_扣减(q_min) / I_样品(q_min)，是个分数，报成百分数。"""
    return "—" if v is None else _fmt(v * 100, "%")


def _fmt(x, unit="", nd=2):
    if x is None:
        return "—"
    try:
        return (("%." + str(nd) + "f") % float(x)) + unit
    except (TypeError, ValueError):
        return str(x)


def _one_liner(sm, present):
    """一句话结论：先给"能不能用"，再说不能用的部分与原因。只依据 summary.json 判状态。"""
    ift = sm.get("ift") or {}
    sub = sm.get("subtraction") or {}
    auto = sm.get("guinier_auto") or {}
    rec = sm.get("guinier_recommended")
    shape = sm.get("shape")
    denss = shape.get("denss") if isinstance(shape, dict) else None
    shape_msg = denss if isinstance(denss, str) else None        # 失败/跳过时这里是一句话
    denss = denss if isinstance(denss, dict) else None
    _dm = shape.get("dammif") if isinstance(shape, dict) else None
    dammif_msg = _dm if isinstance(_dm, str) else None
    dammif_ok = ([r for r in _dm if isinstance(r, dict) and "failed" not in r]
                 if isinstance(_dm, list) else [])
    _dv = shape.get("damaver") if isinstance(shape, dict) else None
    damaver = _dv if isinstance(_dv, dict) else None
    runs = ift.get("runs") or []
    ran, trusted = bool(runs), bool(ift.get("trusted"))
    crow = next((r for r in runs if r.get("tag") == ift.get("chosen") and "failed" not in r), None)

    c = sub.get("contrast_lowq")
    weak = c is not None and float(c) < 0.02
    rg_auto = auto.get("rg")
    rg_real = crow.get("rg_realspace") if crow else None

    bits = []
    # ① 先给"能不能用"——对比度太低时不能因为闸门过了就说可用
    if weak:
        bits.append("**这批数据基本不能用**：低 q 对比度只有 %s（信号埋在背景里），"
                    "下面的 Rg/P(r) 都可能是在拟合背景" % _fmt_sub_contrast(c))
    elif rec and all((rec.get("gates") or {}).values()):
        bits.append("低 q 与整体都可用")
    elif auto:
        bits.append("整体可用，但**低 q 不服从单一 Guinier 定律**（多分散或残留背景），"
                    "报 Rg 必须带上 q 区间")
    # ② 形状是不是像大颗粒（并指出与 P(r) 的矛盾）
    if rg_auto and float(rg_auto) > 40:
        extra = ""
        if rg_real:
            extra = ("；同一份数据的 P(r) 给 Rg(实空间) = %s Å，两者差得远 —— "
                     "低 q 还有额外成分，别只报 auto-Guinier 那个数" % _fmt(rg_real, "", 1))
        bits.append("auto-Guinier 的 Rg = %s Å 偏大（大颗粒/聚集或多分散的迹象）%s"
                    % (_fmt(rg_auto, "", 0), extra))
    # ③ P(r)
    if not ran:
        bits.append("本次运行**没有跑 IFT/P(r)**（`--steps` 里没包含）")
    elif not trusted:
        bits.append("**P(r) 不可信**（IFT 闸门没过）→ 不报 P(r)、不建 3D，只能用高 q 的形状信息")
    else:
        bits.append("P(r) 算出来了：Dmax = %s Å，Rg(实空间) = %s Å%s"
                    % (_fmt(crow.get("dmax") if crow else None, "", 1), _fmt(rg_real, "", 1),
                       "（但对比度太低，别当结论用）" if weak else ""))
    # ④ 形状重建：电子云（DENSS）与珠模（DAMMIF/DAMAVER）
    shapes = []
    if denss:
        shapes.append("电子云(DENSS)已出")
    elif shape_msg:
        shapes.append("电子云(DENSS)**失败**（%s）" % str(shape_msg)[:60])
    if dammif_ok:
        chis = min(float(r["chisq"]) for r in dammif_ok)
        shapes.append("珠模 %d 个（最好 χ² = %s%s）"
                      % (len(dammif_ok), _fmt(chis, "", 1),
                         "，先修曲线" if chis > 5 else ""))
    elif dammif_msg:
        shapes.append("没出珠模（%s）" % str(dammif_msg)[:40])
    if damaver:
        shapes.append("平均 NSD %s%s"
                      % (_fmt(damaver.get("mean_nsd"), "", 3),
                         ("（%d 个 cluster）" % len(damaver.get("clusters") or [])
                          if (damaver.get("clusters") or []) else "")))
    if shapes:
        bits.append("；".join(shapes))
    elif ran and trusted:
        bits.append("形状重建没跑（`--model-engine none`）")
    return "；".join(bits) + "。"


def write_readme(out, extra_note=None):
    out = os.path.abspath(os.path.expanduser(out))
    sj = os.path.join(out, "summary.json")
    if not os.path.exists(sj):
        return None
    sm = json.load(open(sj))
    present = sorted(os.path.relpath(p, out) for p in glob.glob(os.path.join(out, "**", "*"), recursive=True)
                     if os.path.isfile(p))
    rel = set(present)
    sample = sm.get("sample") or os.path.basename(out)
    sub = sm.get("subtraction") or {}
    auto = sm.get("guinier_auto") or {}
    rec = sm.get("guinier_recommended") or {}
    ift = sm.get("ift") or {}
    runs = ift.get("runs") or []
    chosen = ift.get("chosen")
    crow = next((r for r in runs if r.get("tag") == chosen and "failed" not in r), None)
    mw = sm.get("mw") or {}
    sh = sm.get("shape")
    denss = sh.get("denss") if isinstance(sh, dict) else None
    shape_msg = denss if isinstance(denss, str) else None      # 失败/跳过时这里是一句话
    denss = denss if isinstance(denss, dict) else None
    dammif = sh.get("dammif") if isinstance(sh, dict) else None
    dammif_msg = dammif if isinstance(dammif, str) else None   # "skipped: ..." 之类
    dammif = [r for r in dammif if isinstance(r, dict)] if isinstance(dammif, list) else []
    dammif_ok = [r for r in dammif if "failed" not in r]
    damaver = sh.get("damaver") if isinstance(sh, dict) else None
    damaver_msg = damaver if isinstance(damaver, str) else None
    damaver = damaver if isinstance(damaver, dict) else None

    L = []
    A = L.append
    A("# %s —— 管式/静态 SAXS 处理结果" % sample)
    A("")
    A("> 自动生成 %s ｜ 脚本 `run-raw-tube-pipeline.py` ｜ 原始帧 `%s` ｜ 配置 `%s`"
      % (time.strftime("%Y-%m-%d %H:%M"), sm.get("sample_dir", "?"), sm.get("cfg", "?")))
    A("")
    A("## 一句话结论")
    A("")
    A(_one_liner(sm, rel))
    A("")
    A("## 1. 关键结果")
    A("")
    A("| 项目 | 数值 | 怎么看 |")
    A("|---|---|---|")
    A("| 用哪些帧当背景 | %s | 文件名前缀 = 样品名的是样品，其余都是背景；每次运行前先确认 |"
      % (", ".join(sorted((sm.get("control_runs") or sub.get("control_runs") or {}).keys()))
         or "见 norm/frame_qc.csv"))
    A("| 高 q 缩放因子 | %s（在 q 0.30–0.44 Å⁻¹ 定标） | 偏离 1 超过 5%% 说明背景与样品不匹配，"
      "扣减后会有常数残留 |" % _fmt(sub.get("control_scale_factor"), "", 4))
    A("| 低 q 对比度 | %s（扣减曲线 q_min 处 ÷ 样品同点） | <2%% 时低 q 不可信 |"
      % _fmt_sub_contrast(sub.get("contrast_lowq")))
    A("| auto-Guinier | Rg = %s Å（q %s–%s，R² = %s） | RAW 自动选的区间，仅作参考，"
      "不要直接引用 |"
      % (_fmt(auto.get("rg"), "", 1), _fmt(auto.get("qmin"), "", 4),
         _fmt(auto.get("qmax"), "", 4), _fmt(auto.get("r_sqr"), "", 3)))
    if rec:
        A("| 推荐 Guinier 区间 | q %s–%s，Rg = %s ± %s Å，qRg ≤ %s，R² = %s，%s | "
          "这是本流程推荐引用的 Rg，报告里必须写 q 区间 |"
          % (_fmt(rec.get("qmin"), "", 4), _fmt(rec.get("qmax"), "", 4),
             _fmt(rec.get("rg"), "", 2), _fmt(rec.get("rg_err"), "", 2),
             _fmt(rec.get("qrg_max")), _fmt(rec.get("r_sqr"), "", 4),
             "全过闸门" if all((rec.get("gates") or {}).values()) else
             "**未全过闸门**：" + ",".join(k for k, v in (rec.get("gates") or {}).items() if not v)))
    else:
        A("| 推荐 Guinier 区间 | **没有区间通过全部闸门** | 低 q 不服从单一 Guinier 定律，"
          "Rg 只能给范围，见上面那张核心表 |")
    if crow:
        ok = bool(ift.get("trusted"))
        de = crow.get("dmax_err")
        de_s = ("" if de is None or not np.isfinite(float(de))
                else " ± %s" % _fmt(de, "", 1))          # GNOM 不给 Dmax 误差 -> 别印 ±nan
        A("| P(r) / IFT | %s%s | %s |"
          % ("Dmax = %s%s Å，Rg(实空间) = %s Å，chisq = %s（引擎/起跑点：%s）"
             % (_fmt(crow.get("dmax"), "", 1), de_s,
                _fmt(crow.get("rg_realspace"), "", 1), _fmt(crow.get("chisq"), "", 2), chosen),
             "" if ok else " ——**不可信，不要引用**（闸门没过，见第 4 节）",
             "三闸门（Rg 对得上 / Dmax ≤ 4.5·Rg / Dmax 未越出搜索网格）**全过才算结果**"))
    elif ift.get("runs"):
        A("| P(r) / IFT | **没解出来**（BIFT 两次起跑都失败） | 见 `tables/ift_summary.csv` |")
    else:
        A("| P(r) / IFT | **本次没跑**（`--steps` 里没包含 ift） | 重跑时不加 `--steps` |")
    vp = mw.get("Vp_porod")
    if isinstance(vp, dict) and vp.get("mw") is not None:
        aux = (vp.get("aux") or [])
        A("| 分子量（Porod 体积法） | %s kDa%s | 会高估 ~1.5×，只作数量级参考 |"
          % (_fmt(vp.get("mw"), "", 1),
             ("（Porod 体积 %s Å³）" % _fmt(aux[0], "", 0)) if aux else ""))
    vc = mw.get("Vc")
    if isinstance(vc, dict) and vc.get("mw") is not None:
        A("| 分子量（Vc 法，浓度无关） | %s kDa | **这个更可信**；比单体预期大 2 倍以上多为聚集体 |"
          % _fmt(vc.get("mw"), "", 1))
    if denss:
        A("| 电子云（DENSS） | χ² = %s，模型 Rg = %s Å，支撑体积 = %s Å³ | "
          "模型 Rg 应与 P(r) 的 Rg 接近；χ² 只在 IFT 可信时才有意义 |"
          % (_fmt(denss.get("chi2"), "", 2), _fmt(denss.get("rg_model"), "", 1),
             _fmt(denss.get("support_volume"), "", 0)))
    elif shape_msg:
        A("| 电子云（DENSS） | **没出**：%s | 见第 4 节 |" % str(shape_msg)[:120])
    else:
        A("| 电子云（DENSS） | **本次没跑** | 想出就重跑时去掉 `--model-engine none` |")
    if dammif_ok:
        best_d = min(dammif_ok, key=lambda r: r["chisq"])
        rg_p = crow.get("rg_realspace") if crow else None
        drift = ("；与 P(r) 的 Rg 差 %s%%" % _fmt(abs(best_d["rg"] - rg_p) / rg_p * 100, "", 0)
                 if rg_p else "")
        A("| 珠模（DAMMIF ×%d） | 最好的 χ² = %s；Rg = %s Å；Dmax = %s Å；MW = %s kDa%s | "
          "珠模 χ² **1–3 才算拟合好**；>5 说明曲线（尤其低 q）还不干净，先修数据 |"
          % (len(dammif_ok), _fmt(best_d["chisq"], "", 1), _fmt(best_d["rg"], "", 1),
             _fmt(best_d["dmax"], "", 1), _fmt(best_d["mw"], "", 0), drift))
        bad = [r for r in dammif if "failed" in r]
        if bad:
            A("| 珠模失败次数 | %d 次（%s） | 见 `tables/shape_results.json` 与日志 |"
              % (len(bad), str(bad[0].get("failed"))[:80]))
    elif dammif_msg:
        A("| 珠模（DAMMIF） | **没出**：%s | 装了 ATSAS 就能出；见第 4 节 |" % str(dammif_msg)[:120])
    if damaver:
        cl = damaver.get("clusters") or []
        A("| 珠模一致性（DAMAVER，%s 个模型） | 平均 NSD = %s ± %s；%d 个 cluster%s | "
          "NSD ≲2 说明这些模型是同一类形状；**分成多个 cluster = 形状还没定** |"
          % (_fmt(damaver.get("n_models"), "", 0), _fmt(damaver.get("mean_nsd"), "", 3),
             _fmt(damaver.get("stdev_nsd"), "", 3), len(cl),
             ("；代表模型 `%s`" % damaver.get("rep_model")) if damaver.get("rep_model") else ""))
    elif damaver_msg:
        A("| 珠模一致性（DAMAVER） | %s | — |" % str(damaver_msg)[:120])
    A("")
    A("## 2. 这个文件夹里有什么（按建议阅读顺序）")
    A("")
    for pat, desc, when in FILE_DOC:
        if pat.startswith("<"):
            hit = [p for p in rel if p.endswith("workspace.hdf5")]
            if not hit:
                continue
            shown = hit[0]
        elif "*" in pat:
            hits = sorted(p for p in rel if fnmatch.fnmatch(p, pat))
            if not hits:
                continue
            shown = hits[0] + ("（共 %d 个）" % len(hits) if len(hits) > 1 else "")
        else:
            if pat not in rel:
                continue
            shown = pat
        A("- `%s` —— %s（%s）" % (shown, desc, when))
    A("")
    missing = []
    has_models = bool([p for p in rel if p.startswith("models/")])
    if not ift.get("runs"):
        missing.append("- **本次运行没有跑 IFT/P(r) 节点**（`--steps` 里没包含 ift/mw/shape）："
                       "重跑时不加 `--steps` 即可")
    elif not ift.get("trusted"):
        missing.append("- **没有可信的 P(r)**：IFT 两条起跑都没过闸门（见下表），"
                       "按设计不报 P(r)、不拿它建 3D；只留 `ifts/*_untrusted_*` 供目视")
    if shape_msg:
        missing.append("- **电子云（DENSS）这次失败了**：%s —— 常见于 IFT 解本身不稳或数据在低 q/高 q "
                       "有坏点；可先 `--model-engine none` 把前面的结果定下来，再单独试 "
                       "`--denss-mode Slow`" % str(shape_msg)[:160])
    if dammif_msg:
        missing.append("- **没有珠模（DAMMIF）**：%s。珠模需要 ATSAS（GNOM + DAMMIF + DAMAVER）"
                       "并且只吃 GNOM 的 IFT（`--ift-engine both` 时自动跑 GNOM）"
                       % str(dammif_msg)[:160])
    elif dammif_ok and min(r["chisq"] for r in dammif_ok) > 5:
        missing.append("- **珠模 χ² 偏大（最好 %s）**：这不是 DAMMIF 参数问题，是**曲线本身还不够"
                       "干净**——低 q 的束挡光晕/背景残留会把 χ² 顶到几百（本机实测：同一份数据把"
                       "起算 q 从 0.0064 提到 0.0567 后 χ² 从 556 掉到 1.4）。回去修扣减与 q 窗口，"
                       "别在 DAMMIF 里调参" % _fmt(min(r["chisq"] for r in dammif_ok), "", 1))
    if has_models and not isinstance(sh, dict):
        missing.append("- `models/` 里有文件，但**本次运行没有 3D 记录**（summary.json 里没有）"
                       "：那是更早一次运行留下的，重跑才与本文一致")
    elif not has_models and ift.get("trusted"):
        missing.append("- **没有 `models/`（3D 模型）**：本次 `--model-engine none`，"
                       "想建模型就重跑时去掉这个参数（DENSS 几分钟，DAMMIF 需要 ATSAS）")
    if not [p for p in rel if p.startswith("frames/")]:
        missing.append("- 没有 `frames/`：没加 `--save-frames`（逐帧曲线默认不落盘）")
    A("## 3. 每个数字的判据")
    A("")
    A(RULES_TEXT)
    A("## 4. 这次没做的 / 不能信的")
    A("")
    if missing:
        L.extend(missing)
    else:
        A("- （没有跳过任何节点）")
    if ift.get("runs"):
        A("")
        A("IFT 两次起跑的全部结果：")
        A("")
        A("| 起跑点 | q 范围 | Dmax (Å) | Rg_real (Å) | chisq | Rg 闸门 | Dmax/Rg 闸门 | 网格闸门 | 可信 |")
        A("|---|---|---|---|---|---|---|---|---|")
        for r in ift["runs"]:
            if "failed" in r:
                A("| %s | — | — | — | — | — | — | — | 失败：%s |" % (r.get("tag"), r.get("failed")))
                continue
            g = r.get("gates") or {}
            A("| %s | %s–%s | %s | %s | %s | %s | %s | %s | %s |"
              % (r.get("tag"), _fmt(r.get("qmin"), "", 4), _fmt(r.get("qmax"), "", 4),
                 _fmt(r.get("dmax"), "", 1), _fmt(r.get("rg_realspace"), "", 1),
                 _fmt(r.get("chisq"), "", 2),
                 "过" if g.get("rg_vs_guinier") else "**不过**",
                 "过" if g.get("dmax_over_rg") else "**不过**",
                 "过" if g.get("dmax_within_grid") else "**不过**",
                 "是" if r.get("trusted") else "**否**"))
        A("")
        A("（BIFT 的 Dmax 只是搜索网格，优化步可以跑出去；跑到几百 Å 还配一个很小的 chisq 就是假解。）")
    if extra_note:
        A("")
        A(extra_note)
    A("")
    A("## 5. 想自己复核")
    A("")
    A("- **在 RAW 里交互式看**：`File → Open Workspace` 打开 `%s_workspace.hdf5`，曲线、Guinier、IFT 都在里面。"
      % sample)
    A("- **核对本文件的数字**：`tables/*.csv`、`summary.json`；RAW 自己的原始输出在 `reports/raw_report.pdf`。")
    A("- **重跑**：`run-raw-tube-pipeline.py --sample-dir <原始目录> --cfg <当天.cfg> --out-dir <这里>`"
      "（加 `--denss-mode Slow` 得到更细的 3D，加 `--qmin 0.010` 可跳过光束挡边缘的点）。")
    A("- **低 q 到底能不能信 / 上翘是真还是假**：走 `AgentSkill-DoingSAXS` 里的 "
      "`assess-saxs-raw-data-quality`（空白−空白对照、背景形状失配、2D 差分三道检验）。")
    A("")
    path = os.path.join(out, "README.md")
    with open(path, "w") as fh:
        fh.write("\n".join(L))
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("out_dirs", nargs="+", help="result folder(s) containing summary.json")
    args = ap.parse_args()
    for d in args.out_dirs:
        p = write_readme(d)
        print(p or "no summary.json in %s" % d)


if __name__ == "__main__":
    main()
