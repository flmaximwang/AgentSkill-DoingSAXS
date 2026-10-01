#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""结果目录的**唯一**约定与**唯一**实现（SEC 与管式共用）。

两个入口用它：
  · write-readme.py            → 把结果目录写成给人看的 `README.md`
  · verify-results-folder.py   → 验收：结构齐不齐 + README 里的数字与产物对不对得上

只读产物（`run_meta.json` / `summary.json` / `tables/*.csv` / `models/` / 目录清单），
**不重新拟合、不臆造数字**：README 里出现的每个数字都能在某个文件里找到原文。

两种管线的产物字段名不一样（`rg` vs `Rg`、`q_min` vs `qmin`、`mw` vs `MW_kDa`…），
本模块用**别名表**把它们归一到同一套 canonical 字段，再套同一个章节模板 —— 这就是
「两边的 README 长得一样、关键结果能逐条对齐」的机制所在。
"""
from __future__ import annotations

import csv
import glob
import json
import os
import time

VERSION = "1.0.0"

MODE_LABEL = {"sec": "SEC-SAXS", "tube": "管式/静态 SAXS"}
PRODUCER = {"sec": "run-raw-sec-pipeline.py", "tube": "run-raw-tube-pipeline.py"}
SKILL_NAME = "write-saxs-results-readme"

# =============================================================== 结果目录契约
# (级别, 路径/通配, 说明)：must = 缺了就没法交付；should = 缺了要写明原因；info = 有就说
STRUCTURE = {
    "sec": [
        ("must", "run_meta.json", "本次运行的参数与关键结果（机器可读）"),
        ("must", "profiles/01_integrated", "逐帧积分曲线（含逐帧归一化）"),
        ("must", "profiles/03_subtracted", "逐帧扣减曲线（系列图的数据源）"),
        ("must", "profiles/04_sample", "样品区平均曲线（下游分析用的就是它）"),
        ("must", "tables/guinier_multi_range.csv", "多区间 Guinier 全表"),
        ("must", "tables/ift_summary.csv", "各 IFT 引擎的 Dmax/Rg/χ²"),
        ("must", "tables/mw.csv", "分子量（Vc / Vp / 贝叶斯 / DATCLASS）"),
        ("should", "series", "序列对象、洗脱曲线与区间记录（ranges.json）"),
        ("should", "profiles/02_buffer", "buffer 平均曲线"),
        ("should", "profiles/06_guinier", "每个 q 区间的 Guinier 拟合副本"),
        ("should", "reports", "RAW 自己出的 PDF 报告"),
        ("should", "ifts", "IFT 解（P(r) 与 RAW/ATSAS 原文件）"),
        ("info", "tables/frame_params.csv", "逐帧 Rg/I0/Vc-MW/SEC"),
        ("info", "tables/frames_integrated.csv", "逐帧透射与低 q 强度"),
        ("info", "models", "3D 重建（电子云 .mrc / 珠模 .cif / DAMAVER）"),
        ("info", "video", "归一化裁剪视频与束斑质心轨迹"),
        ("info", "norm", "逐帧归一化因子表"),
    ],
    "tube": [
        ("must", "summary.json", "本次运行的参数、关键结果与逐节点状态（机器可读）"),
        ("must", "tables/guinier_multi_range.csv", "多区间 Guinier 全表 + 每条闸门"),
        ("must", "tables/ift_summary.csv", "各引擎/起跑点的 Dmax/Rg/χ² + 三条闸门"),
        ("must", "tables/mw.csv", "分子量（Vp / Vc / 贝叶斯 / DATCLASS）"),
        ("must", "profiles/01_control", "对照（buffer）平均曲线"),
        ("must", "profiles/02_sample", "样品平均曲线"),
        ("must", "profiles/03_subtracted", "扣减曲线（最终交付的那条）"),
        ("must", "profiles/04_guinier", "每个候选区间一份截断后的曲线"),
        ("must", "ifts", "IFT 解与 P(r)"),
        ("should", "tables/shape_results.json", "形状重建全表（电子云 / 每个珠模 / DAMAVER）"),
        ("should", "qc.png", "四联诊断图（扣减曲线 / Kratky / 多区间 Rg 与 χ² / P(r)）"),
        ("should", "reports", "RAW 自己出的 PDF 报告"),
        ("should", "models", "3D 重建（电子云 .mrc / 珠模 .cif / DAMAVER）"),
        ("info", "norm/frame_qc.csv", "逐帧质检（低 q 电平 / 透射 / 离群）"),
        ("info", "frames", "逐帧 1D 曲线（仅在 --save-frames 时）"),
    ],
}

# 「这个文件夹里有什么」：按**建议阅读顺序**排列（* = 两种模式通用）
FILE_DOC = [
    ("*", "README.md", "本文件：结论、关键结果、判据与复核方式", "正在看"),
    ("tube", "qc.png", "四联诊断图：扣减曲线 / Kratky / 多区间 Rg 与 χ² / P(r)", "先看这张"),
    ("*", "reports/*", "RAW 自己出的报告（Guinier、IFT、MW 的原始输出）", "要核对数字时"),
    ("*", "profiles/03_subtracted/subtracted.dat", "**最终扣减曲线**：q、I(q)、误差三列，纯文本",
     "画图 / 喂别的软件时用这个"),
    ("sec", "profiles/03_subtracted/*_sub.dat", "每一帧扣减 buffer 之后的 I(q)", "想自己重选区间时"),
    ("tube", "profiles/02_sample/sample_avg.dat", "样品帧的平均曲线（未扣背景）", "想检查扣减前后差别时"),
    ("tube", "profiles/01_control/control_avg.dat", "对照（buffer）帧的平均曲线", "同上"),
    ("tube", "profiles/01_control/control_avg_scaled.dat", "对照曲线按高 q 窗缩放后（真正参与扣减的那份）",
     "怀疑背景不匹配时"),
    ("sec", "profiles/04_sample/sample_avg.dat", "样品区平均后的曲线（下游分析基于它）", "复核 Rg / I(0) 时"),
    ("sec", "profiles/02_buffer/buffer_avg.dat", "buffer 平均曲线", "怀疑背景不匹配时"),
    ("sec", "profiles/01_integrated/*", "逐帧积分曲线", "要逐帧看时"),
    ("*", "tables/guinier_multi_range.csv",
     "多区间 Guinier 全表：Rg/I0/误差/qRg/χ²_red/曲率 + 每条闸门", "判 Rg 信不信时（核心表）"),
    ("tube", "tables/guinier_results.json", "上表的机器可读版（含 auto 与 recommended）", "脚本用"),
    ("*", "profiles/04_guinier/*", "每个候选 q 区间各一份曲线（文件名即 q 范围）",
     "想看区间怎么影响 Rg 时"),
    ("sec", "profiles/06_guinier/*", "每个区间的 Guinier 拟合（原始点 + 拟合线）", "想画拟合直线时"),
    ("*", "tables/ift_summary.csv", "各 IFT 引擎（BIFT/GNOM…）的 Dmax/Rg/chisq 与可信闸门",
     "判 P(r) 信不信时（核心表）"),
    ("*", "ifts/pr.dat", "**本次采用的 P(r)**（r、P(r)、误差；文件头写明引擎）", "画 P(r) / 写论文时用这个"),
    ("*", "ifts/pr_gnom.dat", "GNOM 的 P(r)（珠模就是用它算的）", "对比 GNOM 与 BIFT 时"),
    ("*", "ifts/*.ift", "BIFT 的 IFT 解（RAW 原生）", "要复核 BIFT 解时"),
    ("*", "ifts/*.out", "GNOM 的 IFT 解（ATSAS 格式，DAMMIF 的输入）", "要重跑珠模时"),
    ("*", "ifts/ift_fit*.dat", "IFT 拟合曲线：实测 I(q) vs 拟合 I(q)", "看 IFT 拟合得好不好时"),
    ("tube", "ifts/*untrusted*", "**不可信**的 IFT 解（闸门没过，仅留档，不要引用）", "排查为什么没有 P(r) 时"),
    ("*", "tables/mw.csv", "分子量：Porod 体积法（Vp）与浓度无关的 Vc 法", "看分子量时"),
    ("*", "models/*_denss.mrc", "**电子云密度图**（DENSS，可用 ChimeraX/PyMOL 打开）", "要形状 / 内部空腔时"),
    ("*", "models/*support*.mrc", "DENSS 的支撑掩膜（模型边界）", "看模型大小是否合理时"),
    ("*", "models/*_denss.log", "DENSS 自己的日志（每次迭代的 χ²/Rg）", "怀疑没收敛时"),
    ("*", "models/*dammif*.cif", "**珠模**（ATSAS DAMMIF，dummy-atom；直接拖进 ChimeraX/PyMOL）",
     "要形状、要投稿图时"),
    ("*", "models/*dammif*.fit", "珠模对数据的拟合曲线与 χ²", "判珠模拟合好不好时"),
    ("*", "models/*dammif*.log", "DAMMIF 自己的日志（每次起跑的 χ²、珠子数）", "珠模 χ² 大时看这里"),
    ("*", "models/*damaver-global-damaver.cif", "**多个珠模的一致性平均代表模型**", "要给一个模型时用这个"),
    ("*", "models/*-distances.txt", "模型两两差异（里面的 Mean value of nsd 就是平均 NSD）", "判形状定没定时（核心）"),
    ("tube", "tables/shape_results.json", "形状重建全表（电子云 + 每个珠模 + DAMAVER 一致性）", "看 3D 时"),
    ("tube", "norm/frame_qc.csv", "逐帧质检：低 q 电平、透射、与本 run 中位数的偏离、是否离群",
     "怀疑某一帧坏时"),
    ("sec", "norm/*.csv", "逐帧归一化因子（Transmitted_Beam 等）", "怀疑某一帧坏时"),
    ("*", "*workspace.hdf5", "RAW 的 workspace 文件（File → Open Workspace 可交互复核所有曲线）",
     "想在 GUI 里复核时"),
    ("*", "*.hdf5", "序列对象（RAW 打开）", "想在 GUI 里复核时"),
    ("sec", "series/series_plot.png", "洗脱曲线 + 选区阴影", "查选帧选得对不对时"),
    ("sec", "series/ranges.json", "脚本使用的 buffer / 样品区间（0 基帧号）", "复核选帧时"),
    ("sec", "video/*.mp4", "裁剪区（束斑附近）的归一化逐帧视频", "看束位/亮度漂移时"),
    ("sec", "tables/frame_params.csv", "逐帧 Rg / I0 / Vc-MW", "看某个帧异不异常时"),
    ("*", "summary.json", "本 README 里所有数字的机器可读版（管式）", "脚本 / 批处理用"),
    ("*", "run_meta.json", "本 README 里所有数字的机器可读版（SEC）", "脚本 / 批处理用"),
    ("tube", "frames/*", "逐帧的 1D 曲线（只在用了 --save-frames 时才有）", "要逐帧看时"),
]

# 「每个数字的判据」——两种模式共用一张表
RULES = [
    ("高 q 窗对照缩放因子", "1.00±0.02", "±0.05", "±0.05 以上：背景与样品不匹配，扣减后会有常数残留"),
    ("低 q 对比度（扣减后 ÷ 背景电平）", ">5%", "2–5%", "<1%：信号埋在背景里，低 q 什么都不能说"),
    ("Guinier R²", ">0.99", "0.95–0.99", "<0.95：这段不是 Guinier 区（<0.9 或为负 = 比取平均还差）"),
    ("多区间 Rg 一致性", "区间间差 <10%", "10–20%", "随区间乱变 / 漂 >±20%：低 q 不服从单一 Guinier 定律"),
    ("IFT 三条闸门", "全过（trusted）", "只差 rg 那条（见 rg_tol）", "有 Dmax 跑到搜索域顶 / 超 4.5·Rg：**假解**，不许报 P(r)"),
    ("Dmax / Rg", "2.5–4", "4–4.5", ">4.5：P(r) 里多半是大颗粒尾巴"),
    ("Vc 分子量", "与预期相符", "—", "明显偏大 >2×：多半是聚集体（Vp 法高估是常态，Vp≈1.5×Vc 属正常）"),
    ("电子云 DENSS χ²", "越低越好", "—", "只在 IFT 可信时才有意义；模型 Rg 应与 P(r) 的 Rg 接近（差 <10%）"),
    ("DENSS support 体积", "与 MW 自洽", "—", "明显过大（盒子边长 ≫ 3·Dmax）说明 IFT 被低 q 拖长"),
    ("珠模 χ²（DAMMIF）", "1–3", "3–5", ">5：**是数据问题不是参数问题**——先修曲线（低 q 污染 / 背景残留）再跑"),
    ("珠模 Rg / Dmax", "与 P(r) 的差 <10%", "10–20%", "差很多：模型和数据说的不是同一个东西，别混着报"),
    ("DAMAVER 平均 NSD", "≲2（同一类形状）", "2–3", "明显更大、或模型分成多个 cluster：**形状还没定**，加模型数再跑"),
]

ABBREV = [
    "**Rg**：回转半径（Å）。**I(0)**：零角散射强度。**qRg**：无量纲化 q，Guinier 区一般要求 qRg ≲ 1.3。",
    "**IFT / P(r)**：把 I(q) 反变换成实空间距离分布，Dmax 是它的支撑上限。",
    "**Vc / Vp / Bayesian / Datclass**：四种独立的分子量估计（前两个只需曲线，后两个用 ATSAS）。",
    "**DENSS**：电子云（密度图）重建；**DAMMIF**：dummy-atom 珠模；**DAMAVER**：把多个珠模对齐平均并给一致性 NSD。",
    "**χ²**：重建模型与实验曲线的拟合残差（≈1 表示在误差棒范围内）；**NSD**：模型两两之间的归一化空间差异。",
    "**哨兵值**：`rg = -1`、`r² < 0`、`-1.0 kDa` 是 RAW 的**失败标记**，不是物理量。",
]

# 章节顺序 = 对齐契约：两种模式都必须有这 8 节，顺序一致
SECTION_ORDER = [
    "## 0. 结论速览",
    "## 1. 最终拟合参数（明细表）",
    "## 2. 关键结果",
    "## 3. 这个文件夹里有什么（按建议阅读顺序）",
    "## 4. 每个数字的判据",
    "## 5. 这次没做的 / 不能信的",
    "## 6. 本次用的参数（复现用）",
    "## 7. 想自己复核",
]


# =============================================================== 格式化小工具
def _num(x):
    """转 float；None / NaN / 非数一律 None（调用方再决定怎么显示）。"""
    if x is None or x == "":
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def _fmt(x, unit="", nd=2):
    if x is None:
        return "—"
    v = _num(x)
    if v is None:
        return x if isinstance(x, str) else "—"
    if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e5):
        s = "%.2e" % v
    else:
        s = ("%." + str(nd) + "f") % v
    return s + unit


def _fmt_err(x, nd=2, dash="—"):
    """误差列用两位有效数字（误差常比主值小几个量级：0.0038 不该显示成 0.00）。"""
    v = _num(x)
    if v is None:
        return dash
    if v == 0:
        return "0"
    return "%.*g" % (nd, v)


def _finite(x):
    v = _num(x)
    return v is not None and v == v and abs(v) != float("inf")


def _fmt_mw(x):
    """分子量显示：inf / -1（RAW 的失败哨兵）都要说明白，不能当数字印出来。"""
    v = _num(x)
    if v is None:
        return "—"
    if abs(v) == float("inf"):
        return "inf（发散）"
    if v < 0:
        return "未收敛（哨兵 -1）"
    return _fmt(v, "", 1)


def _pct(x, nd=1):
    v = _num(x)
    return "—" if v is None else ("%." + str(nd) + "f%%") % (v * 100)


def _color(s):
    return (s.replace("✅", "🟢").replace("⚠️", "🟡").replace("❌", "🔴"))


def _gate(v):
    """summary.json 里的布尔可能被 default=str 写成 "False"（numpy 标量的历史遗留）。
    字符串 "False" 在 Python 里是**真**值 —— 不显式归一化就会把"没过闸门"读成"全过"。"""
    if isinstance(v, str):
        return v.strip().lower() not in ("false", "0", "no", "none", "")
    return bool(v)


def _cell(s):
    """markdown 表格单元：竖线会破表，换行会断路。"""
    return str(s).replace("|", "\\|").replace("\n", " ").strip()


# =============================================================== 字段别名
# 两种管线的表头不一样，这里统一（canonical 名 → 可能的列名，按优先级）
ALIAS = {
    "label": ("label", "range_label", "tag"),
    "rg": ("rg", "Rg"),
    "rg_err": ("rg_err", "Rg_err"),
    "i0": ("i0", "I0"),
    "i0_err": ("i0_err", "I0_err"),
    "qmin": ("qmin", "q_min"),
    "qmax": ("qmax", "q_max"),
    "qrg_min": ("qrg_min", "qRg_min"),
    "qrg_max": ("qrg_max", "qRg_max"),
    "r2": ("r2", "r_sqr", "R2"),
    "chi2_red": ("chi2_red",),
    "curvature": ("curvature",),
    "rg_split_rel": ("rg_split_rel",),
    "idx_min": ("idx_min",),
    "idx_max": ("idx_max",),
    "dmax": ("dmax", "Dmax"),
    "dmax_err": ("dmax_err", "Dmax_err"),
    "rg_real": ("rg_real", "rg_realspace", "Rg_realspace", "rg"),   # SEC 表里 `rg` = 实空间 Rg
    "rg_real_err": ("rg_real_err", "rg_err", "Rg_err"),
    "chisq": ("chisq", "chi_sq"),
    "method": ("method", "engine"),
    "mw": ("mw", "MW_kDa"),
    "pass_rg": ("pass_rg",),
    "pass_dmax_over_rg": ("pass_dmax_over_rg",),
    "pass_dmax_within_grid": ("pass_dmax_within_grid",),
    "trusted": ("trusted",),
}


def field(row, key, default=None):
    for name in ALIAS.get(key, (key,)):
        if name in row and row[name] not in (None, ""):
            return row[name]
    return default


def _read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="", errors="ignore") as fh:
        return list(csv.DictReader(fh))


def _read_mw_csv(path):
    """`tables/mw.csv` 的读取：**不用 DictReader**。

    管式旧版把 `aux` 按 Python 列表直接写进一行（内含逗号 → 列被撑开、
    表头对不上数据），DictReader 会把它切碎。这里按位置读：第 1 列 = 方法、第 2 列 = MW，
    其余按 CSV 行还原成 aux 文本。SEC 的 6 列格式（method,mw,detail1..4）同样吃。"""
    if not os.path.exists(path):
        return []
    out = []
    with open(path, newline="", errors="ignore") as fh:
        rd = csv.reader(fh)
        try:
            header = next(rd)
        except StopIteration:
            return []
        ncol = len(header)
        for fields in rd:
            if not fields:
                continue
            row = {"method": fields[0], "mw": fields[1] if len(fields) > 1 else None}
            if ncol >= 6 and len(fields) >= 6 and not fields[2].lstrip().startswith("["):
                for i in range(2, 6):
                    row["detail%d" % (i - 1)] = fields[i]
            else:
                rest = ",".join(fields[2:]).rstrip(",").strip()
                if rest:
                    row["aux"] = rest
                if fields[-1] and not fields[-1].lstrip().startswith(("[", "'", '"')) \
                        and fields[-1].strip():
                    row["note"] = fields[-1]
            out.append(row)
    return out


def _read_json(path, default=None):
    if not os.path.exists(path):
        return default
    try:
        with open(path, errors="ignore") as fh:
            return json.load(fh)
    except Exception:
        return default


def _rel_files(out):
    return sorted(os.path.relpath(p, out)
                  for p in glob.glob(os.path.join(out, "**", "*"), recursive=True)
                  if os.path.isfile(p) and not os.path.basename(p).startswith("."))


# =============================================================== 统一选值策略
# SEC 的 run_meta.json 把表存成**位置参数列表**（不是 dict），列序与同目录 CSV 一致。
_SEC_G_KEYS = ("label", "rg", "i0", "rg_err", "i0_err", "qmin", "qmax", "qrg_min", "qrg_max", "r2")
_SEC_I_KEYS = ("method", "dmax", "dmax_err", "rg_real", "rg_err", "chisq", "log_alpha", "evidence")


def _as_dict(row, keys):
    if isinstance(row, dict):
        return row
    try:
        return dict(zip(keys, list(row)))
    except TypeError:
        return {}


def guinier_rows(raw_rows):
    """归一成 canonical 行：label / rg / rg_err / i0 / i0_err / qmin / qmax / qRg / r2。"""
    rows = []
    for r in raw_rows or []:
        r = _as_dict(r, _SEC_G_KEYS)
        rg = _num(field(r, "rg"))
        label = str(field(r, "label") or "").strip()
        if label in ("", "?"):
            lo, hi = _num(field(r, "idx_min")), _num(field(r, "idx_max"))
            label = ("idx %s-%s" % (_fmt(lo, "", 0), _fmt(hi, "", 0))
                     if lo is not None else "?")
        rows.append(dict(
            label=label, rg=rg,
            rg_err=_num(field(r, "rg_err")), i0=_num(field(r, "i0")),
            i0_err=_num(field(r, "i0_err")), qmin=_num(field(r, "qmin")),
            qmax=_num(field(r, "qmax")), qrg_min=_num(field(r, "qrg_min")),
            qrg_max=_num(field(r, "qrg_max")), r2=_num(field(r, "r2")),
            chi2_red=_num(field(r, "chi2_red")),
            split_rel=_num(field(r, "rg_split_rel")),
            gates={k: _gate(field(r, "pass_" + k)) for k in
                   ("n", "rg", "qRg", "chi2", "stable") if field(r, "pass_" + k) is not None},
            sentinel=bool(rg is not None and rg < 0)))
    return rows


def declared_verdict(declared):
    """产物自己声明的推荐区间怎么读：'ok'（闸门全过）/ 'rejected'（有闸门没过）/ None（没声明）。"""
    if not isinstance(declared, dict) or not declared:
        return None
    rg = _num(declared.get("rg") if declared.get("rg") is not None
              else field(declared, "rg"))
    if rg is None or rg <= 0:
        return "rejected"
    gates = declared.get("gates") or {}
    if gates and not all(_gate(v) for v in gates.values()):
        return "rejected"
    return "ok"


def pick_adopted_guinier(rows, declared=None):
    """「采用区间」的**唯一**口径（两种模式同一条）：

    ① 产物自己声明了推荐区间**且闸门全过**（管式的 `guinier_recommended`）→ 用它；
    ② 产物声明了推荐区间但**有闸门没过**（它的 `reason` 里写着 NOT recommended）→ **不采用**，
       退回通用规则，并在 README 里写明这句话；
    ③ 通用规则：收敛（rg>0）+ qRg 落在 [0.28, 1.35] + r² 最高；
    ④ 都没有 → None（**没有可引用的 Rg**，不许报数）。
    """
    if declared_verdict(declared) == "ok":
        d = dict(declared)
        d["_source"] = "producer"
        d["_label"] = str(declared.get("label") or field(declared, "label") or "recommended")
        d["_rg"] = _num(d.get("rg") if d.get("rg") is not None else field(d, "rg"))
        d["_r2"] = _num(d.get("r_sqr") if d.get("r_sqr") is not None else field(d, "r2"))
        d["_qmin"], d["_qmax"] = _num(d.get("qmin")), _num(d.get("qmax"))
        d["_rg_err"] = _num(d.get("rg_err"))
        d["_i0"], d["_i0_err"] = _num(d.get("i0")), _num(d.get("i0_err"))
        d["_qrg_min"] = _num(field(d, "qrg_min"))
        d["_qrg_max"] = _num(field(d, "qrg_max"))
        d["_split_rel"] = _num(field(d, "rg_split_rel"))
        if d["_rg"] is not None and d["_rg"] > 0:
            return d
    cand = [r for r in rows if r["rg"] and r["rg"] > 0 and r["r2"] is not None
            and (r["qrg_min"] is None or r["qrg_min"] >= 0.28)
            and (r["qrg_max"] is None or r["qrg_max"] <= 1.35)]
    if not cand:
        return None
    r = max(cand, key=lambda t: t["r2"])
    d = dict(r)
    d["_source"] = "auto-rule"
    d["_label"] = r["label"]
    d["_rg"], d["_r2"] = r["rg"], r["r2"]
    d["_qmin"], d["_qmax"] = r["qmin"], r["qmax"]
    d["_rg_err"], d["_i0"], d["_i0_err"] = r["rg_err"], r["i0"], r["i0_err"]
    d["_qrg_min"] = r["qrg_min"]
    d["_qrg_max"] = r["qrg_max"]
    d["_split_rel"] = r.get("split_rel")
    return d


def guinier_spread(rows):
    """可用区间之间的 Rg 最大相对差（%）。

    只统计「收敛（rg>0）且 r² ≥ 0.9」的区间——未收敛（哨兵 -1）与明显很差的区间本来就不能引用，
    让它们参与一致性判据等于一票否决。被排除的**列出来**，便于复核。
    """
    ok, bad = [], []
    for r in rows:
        if r["rg"] and r["rg"] > 0 and r["r2"] is not None and r["r2"] >= 0.9:
            ok.append((r["label"], r["rg"]))
        else:
            bad.append(r["label"])
    if len(ok) < 2:
        return None, [l for l, _ in ok], bad
    vals = [v for _, v in ok]
    return (100.0 * (max(vals) - min(vals)) / (sum(vals) / len(vals)),
            [l for l, _ in ok], bad)


def ift_rows(raw_rows, rg_guinier=None):
    """归一成 canonical IFT 行，并**统一**算三条闸门（管式产物自带的闸门优先）。

    幂等：已经归一过的行（带 `_canonical`）只补算缺失的闸门，不会把产物给的闸门覆盖掉。
    """
    rows = []
    for r in raw_rows or []:
        r = _as_dict(r, _SEC_I_KEYS)                 # SEC 的表是位置列表，先变 dict
        if r.get("_canonical"):                      # 再进来一次（例如补算闸门）
            d = dict(r)
            if rg_guinier:
                g = d.setdefault("gates", {})
                if g.get("rg_vs_guinier") is None and d.get("rg_real"):
                    g["rg_vs_guinier"] = abs(d["rg_real"] - rg_guinier) / rg_guinier <= 0.15
                if g.get("dmax_over_rg") is None and d.get("dmax"):
                    g["dmax_over_rg"] = d["dmax"] / rg_guinier <= 4.5
            rows.append(d)
            continue
        if "failed" in r:
            rows.append(dict(failed=str(r.get("failed")), tag=str(r.get("tag") or "?"),
                             engine=str(field(r, "method") or "?"), _canonical=True))
            continue
        dmax = _num(field(r, "dmax"))
        rr = _num(field(r, "rg_real"))
        ch = _num(field(r, "chisq"))
        gates = {}
        if field(r, "pass_rg") is not None:
            gates = {"rg_vs_guinier": _gate(field(r, "pass_rg")),
                     "dmax_over_rg": _gate(field(r, "pass_dmax_over_rg")),
                     "dmax_within_grid": _gate(field(r, "pass_dmax_within_grid"))}
        else:                                     # SEC：自己算（网格闸门无法判，标 None）
            if rg_guinier and rr:
                gates["rg_vs_guinier"] = abs(rr - rg_guinier) / rg_guinier <= 0.15
            if dmax and rg_guinier:
                gates["dmax_over_rg"] = dmax / rg_guinier <= 4.5
            gates["dmax_within_grid"] = None
        rows.append(dict(
            tag=str(r.get("tag") or field(r, "label") or field(r, "method") or "?"),
            engine=str(field(r, "method") or "?"),
            qmin=_num(field(r, "qmin")), qmax=_num(field(r, "qmax")),
            dmax=dmax, dmax_err=_num(field(r, "dmax_err")),
            rg_real=rr, rg_err=_num(field(r, "rg_real_err")), chisq=ch,
            gates=gates, trusted=(_gate(field(r, "trusted")) if field(r, "trusted") is not None
                                  else None),
            total_est=_num(r.get("total_est")), alpha=_num(r.get("alpha")),
            quality=r.get("quality"), _canonical=True))
    return rows


def pick_adopted_ift(rows, declared_trusted=None, chosen=None):
    """「采用那一支 P(r)」的**唯一**口径：先看产物声明的 chosen/trusted，再看闸门，再看 χ²。"""
    ok = [r for r in rows if "failed" not in r]
    if not ok:
        return None
    if chosen:
        hit = next((r for r in ok if r["tag"] == chosen or r["engine"] == chosen), None)
        if hit and (hit.get("trusted") is not False):
            return hit
    trusted = [r for r in ok if r.get("trusted") is True]
    if trusted:
        return min(trusted, key=lambda r: (r["chisq"] if r["chisq"] is not None else 9e9))
    passing = [r for r in ok if all(v is True for k, v in r["gates"].items()
                                    if k != "dmax_within_grid" and v is not None)]
    if passing:
        return min(passing, key=lambda r: (r["chisq"] if r["chisq"] is not None else 9e9))
    return None


# =============================================================== 收集：SEC
def collect_sec(out):
    meta = _read_json(os.path.join(out, "run_meta.json"), {}) or {}
    args = meta.get("args") or {}
    st = meta.get("settings") or {}
    ranges = _read_json(os.path.join(out, "series", "ranges.json"), {}) or {}
    f = {
        "mode": "sec",
        "sample": meta.get("prefix") or os.path.basename(os.path.normpath(out)),
        "input": meta.get("input") or args.get("series_dir"),
        "cfg": meta.get("cfg") or args.get("cfg"),
        "generated_at": meta.get("timestamp") or time.strftime("%Y-%m-%d %H:%M:%S"),
        "n_frames": meta.get("n_frames"),
        "steps": meta.get("steps") or [],
        "args": args, "raw_settings": st,
        "command": args.get("_command"),
        "normalization": {
            "on": bool(st.get("EnableNormalization")),
            "how": "BL19U2 逐帧 header txt 的 Transmitted_Beam",
            "detail": "ImageHdrFormat=%s，NormalizationList=%s" % (
                st.get("ImageHdrFormat", "None"), st.get("NormalizationList", "—")),
        },
        "subtraction": dict(
            kind="series（连续洗脱）",
            control_runs=None,
            buffer=ranges.get("buffer") or [],
            sample=ranges.get("sample") or [],
            control_scale_factor=None, scale_window=None, contrast_lowq=None,
            qmin=_num((meta.get("dmax") is not None and None) or None), qmax=None,
            n_points=None, n_frame_outliers=None),
        "guinier_rows": guinier_rows(meta.get("guinier_rows")
                                     or _read_csv(os.path.join(out, "tables",
                                                               "guinier_multi_range.csv"))),
        "ift_raw": ift_rows(meta.get("ift_rows")
                            or _read_csv(os.path.join(out, "tables", "ift_summary.csv"))),
    }
    # 扣减/分析窗信息在表里（SEC 的 run_meta 不带 subtraction 块）
    gi = _read_csv(os.path.join(out, "tables", "frames_integrated.csv"))
    if gi:
        f["subtraction"]["qmin"] = _num(gi[0].get("q_min"))
        f["subtraction"]["qmax"] = _num(gi[0].get("q_max"))
        f["subtraction"]["n_points"] = len(gi)
    # 分子量
    f["mw"] = _mw_from_rows(_read_mw_csv(os.path.join(out, "tables", "mw.csv")), kind="sec")
    f["shape"] = _shape_from_sec(out, meta)
    f["denss_mode"] = args.get("denss_mode")
    f["review_target"] = None
    return f


def _mw_from_rows(rows, kind):
    """SEC 把 RAW 的返回值写成 `detail1..4` 列；管式写成 `aux` 列表。两者的顺序**相同**：

    `mw_vc` → (mw, vcor, mw_err, qmax)；`mw_vp` → (mw, pvol_cor, pvol, qmax)；
    `mw_bayes` → (mw, mw_prob, ci_lower, ci_upper, ci_prob)；`mw_datclass` → (mw, 形状分类, …)。
    """
    mw = {}
    for r in rows:
        m = str(field(r, "method") or "").strip()
        v = _num(field(r, "mw"))
        aux = r.get("aux")
        if isinstance(aux, str):
            aux = _parse_aux(aux)
        aux = aux or []

        def slot(i):
            d = _num(r.get("detail%d" % (i + 1)))
            return d if d is not None else _num(aux[i] if len(aux) > i else None)

        def stext(i):
            d = r.get("detail%d" % (i + 1))
            if d not in (None, ""):
                return d
            return aux[i] if len(aux) > i else None

        if m == "Vc":
            mw["Vc"] = dict(mw=v, err=slot(1), vcor=slot(0), qmax=slot(2))
        elif m in ("Vp", "Vp_porod"):
            mw["Vp"] = dict(mw=v, vol=slot(0), pvol=slot(1), qmax=slot(2))
        elif m in ("Bayesian", "datmw_bayes"):
            mw["Bayesian"] = dict(mw=v, prob=slot(0), ci=(slot(1), slot(2)), ci_prob=slot(3))
        elif m in ("Datclass", "datclass"):
            mw["Datclass"] = dict(mw=v, cls=stext(0), second=slot(1))
        elif v is None and m:
            mw[m] = dict(mw=None, note=str(r.get("note") or ""))
    return mw


def _parse_aux(s):
    """把 mw.csv 的 aux 文本解析成列表。

    合法 JSON 直接用；`inf`/`nan`/`'compact'` 这类不是合法 JSON 的走回退：
    逐个 token 试 float（`inf` → float('inf')），试不动就当字符串留着。
    """
    s = str(s).strip()
    if not s:
        return []
    try:
        return json.loads(s.replace("'", '"'))
    except Exception:
        pass
    out = []
    for tok in s.strip("[]").split(","):
        t = tok.strip().strip("'\"")
        if not t:
            continue
        try:
            out.append(float(t))
        except ValueError:
            out.append(t)
    return out


def _shape_from_sec(out, meta):
    denss, dammif, damaver = None, [], None
    dr = meta.get("denss_rows") or []
    if dr and isinstance(dr[0], list) and len(dr[0]) >= 5:
        # SEC 的 denss_rows = [密度数组字符串, chi², Rg_model, support 体积, 盒子边长]
        denss = dict(chi2=_num(dr[0][1]), rg_model=_num(dr[0][2]),
                     support_volume=_num(dr[0][3]), side=_num(dr[0][4]))
    elif isinstance(meta.get("shape"), dict):
        pass
    dm = meta.get("dammif_rows") or []
    for i, r in enumerate(dm, 1):
        if isinstance(r, list) and len(r) >= 4:
            dammif.append(dict(model="dammif_%02d" % i, chisq=_num(r[0]), rg=_num(r[1]),
                               dmax=_num(r[2]), mw=_num(r[3])))
        elif isinstance(r, dict):
            dammif.append(r)
    dv = meta.get("damaver_row")
    if isinstance(dv, list) and dv:
        ncl = glob.glob(os.path.join(out, "models", "*damaver-cluster*-summary.txt"))
        damaver = dict(mean_nsd=_num(dv[0]), stdev_nsd=_num(dv[1]) if len(dv) > 1 else None,
                       clusters=ncl or None, rep_model=None)
    return dict(denss=denss, dammif=dammif, damaver=damaver)


# =============================================================== 收集：管式
def _tube_ift_rows(out, ift):
    """管式的 IFT 行：**优先读 `tables/ift_summary.csv`**（它带 engine 列与三条闸门列），
    summary.json 只用来补那些失败、没进表的起跑点。"""
    rows = _read_csv(os.path.join(out, "tables", "ift_summary.csv"))
    runs = [r for r in (ift.get("runs") or []) if isinstance(r, dict)]
    if not rows:
        rows = [dict(r) for r in runs]
        for r in rows:
            if not r.get("engine"):
                tag = str(r.get("tag") or "").lower()
                r["engine"] = "GNOM" if "gnom" in tag else ("BIFT" if "start" in tag else "?")
    else:
        have = {str(r.get("tag")) for r in rows}
        for r in runs:
            if "failed" in r and str(r.get("tag")) not in have:
                rows.append(r)
    return ift_rows(rows)


def collect_tube(out):
    sm = _read_json(os.path.join(out, "summary.json"), {}) or {}
    sub = sm.get("subtraction") or {}
    rec = sm.get("guinier_recommended")
    if isinstance(rec, dict):
        rec = dict(rec)
        rec.setdefault("label", "recommended")
    ift = sm.get("ift") or {}
    ift = ift if isinstance(ift, dict) else {}
    sh = sm.get("shape")
    sh = sh if isinstance(sh, dict) else {}
    denss = sh.get("denss")
    denss = denss if isinstance(denss, dict) else None
    dm = sh.get("dammif")
    dammif = [r for r in dm if isinstance(r, dict) and "failed" not in r] if isinstance(dm, list) else []
    dammif_bad = [r for r in dm if isinstance(r, dict) and "failed" in r] if isinstance(dm, list) else []
    dv = sh.get("damaver")
    dv = dv if isinstance(dv, dict) else None
    n_models = sm.get("n_models")
    f = {
        "mode": "tube",
        "sample": sm.get("sample") or os.path.basename(os.path.normpath(out)),
        "input": sm.get("sample_dir"),
        "cfg": sm.get("cfg"),
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "n_frames": sm.get("n_frames"),
        "n_frame_outliers": sm.get("n_frame_outliers"),
        "steps": [s for s in ("integrate", "average", "scale", "subtract", "guinier", "ift",
                              "mw", "shape", "report", "workspace")],
        "args": {}, "raw_settings": {},
        "command": sm.get("command"),
        "atsas_dir": sm.get("atsas_dir"), "has_atsas": sm.get("has_atsas"),
        "normalization": {
            "on": True,   # 管式：RAW 逐帧读 BL19U2 txt（cfg 里两个开关由脚本强制打开）
            "how": "BL19U2 逐帧 header txt 的 Transmitted_Beam",
            "detail": "见 norm/frame_qc.csv 的逐帧 TB/SR",
        },
        "subtraction": dict(
            kind="管式（目录内样品 + 夹着它的对照）",
            control_runs=sorted((sm.get("control_runs") or {}).keys()) or None,
            buffer=[], sample=[],
            control_scale_factor=_num(sub.get("control_scale_factor")),
            scale_window=sub.get("scale_window"),
            contrast_lowq=_num(sub.get("contrast_lowq")),
            qmin=_num(sub.get("qmin")), qmax=_num(sub.get("qmax")),
            n_points=sub.get("n_points"),
            sample_over_control_midq=_num(sub.get("sample_over_control_midq")),
            negative_kratky_frac=_num(sub.get("negative_kratky_frac")),
            n_frame_outliers=sm.get("n_frame_outliers")),
        "guinier_rows": guinier_rows(
            _read_csv(os.path.join(out, "tables", "guinier_multi_range.csv"))),
        "ift_raw": _tube_ift_rows(out, ift),
        "ift_chosen": ift.get("chosen"), "ift_declared_trusted": ift.get("trusted"),
        "gnom": ift.get("gnom") if isinstance(ift.get("gnom"), dict) else None,
        "denss_mode": (denss or {}).get("mode") if denss else None,
        "mw": _mw_from_rows(_read_mw_csv(os.path.join(out, "tables", "mw.csv")), kind="tube"),
        "mw_note": (sm.get("mw") or {}).get("note") if isinstance(sm.get("mw"), dict) else None,
        "shape": dict(denss=denss, dammif=dammif, dammif_bad=dammif_bad,
                      damaver=dv, shape_msg=(sh.get("denss") if isinstance(sh.get("denss"), str)
                                             else None),
                      dammif_msg=(dm if isinstance(dm, str) else None),
                      damaver_msg=(dv.get("_msg") if isinstance(dv, dict) else
                                   (dv if isinstance(dv, str) else None)),
                      engine=sh.get("engine"), n_models=n_models),
        "auto_declared": sm.get("guinier_auto"),
        "rec_declared": rec,
    }
    return f


def collect(out, mode=None):
    mode = mode or detect_mode(out)
    if mode == "sec":
        f = collect_sec(out)
    elif mode == "tube":
        f = collect_tube(out)
    else:
        return None
    f["out"] = os.path.abspath(out)
    f["files"] = _rel_files(out)
    f["mode_label"] = MODE_LABEL[mode]
    f["producer"] = PRODUCER[mode]
    return _finish(f)


def detect_mode(out):
    if os.path.exists(os.path.join(out, "run_meta.json")):
        return "sec"
    if os.path.exists(os.path.join(out, "summary.json")):
        return "tube"
    return None


def _finish(f):
    """把原始行归一成 canonical 事实：采用区间 / 一致性 / 采用的那支 P(r)。"""
    rows = f["guinier_rows"]
    adopted = pick_adopted_guinier(rows, f.get("rec_declared"))
    spread, used, excluded = guinier_spread(rows)
    rg_for_gates = adopted["_rg"] if adopted else None
    if f["mode"] == "tube" and f.get("auto_declared") and not rg_for_gates:
        rg_for_gates = _num((f["auto_declared"] or {}).get("rg")) or None
    f["ift"] = ift_rows(f["ift_raw"], rg_guinier=rg_for_gates)
    f["ift_adopted"] = pick_adopted_ift(f["ift"], f.get("ift_declared_trusted"),
                                        f.get("ift_chosen"))
    if f.get("ift_declared_trusted") is not None:      # 管式：产物自己给了 trusted
        f["ift_trusted"] = any(r.get("trusted") is True for r in f["ift"]
                               if "failed" not in r)
    else:                                              # SEC：用本模块统一算的三条闸门
        ad_ift = f["ift_adopted"]
        f["ift_trusted"] = bool(ad_ift and ad_ift.get("dmax") and all(
            v is True for k, v in (ad_ift.get("gates") or {}).items()
            if k != "dmax_within_grid" and v is not None))
    rec = f.get("rec_declared")
    f["rec_verdict"] = declared_verdict(rec)
    f["rec_reason"] = (rec.get("reason") if isinstance(rec, dict) else None)
    f["guinier_adopted"] = adopted
    f["guinier_spread"] = spread
    f["guinier_used"] = used
    f["guinier_excluded"] = excluded
    used_rgs = [r["rg"] for r in rows
                if r["label"] in used and r["rg"] and r["rg"] > 0]
    f["guinier_used_rg"] = (min(used_rgs), max(used_rgs)) if used_rgs else None
    auto = None
    if f["mode"] == "tube" and isinstance(f.get("auto_declared"), dict):
        auto = f["auto_declared"]
    else:
        auto = next((r for r in rows if r["label"] == "auto"), None)
    f["guinier_auto"] = auto
    return f


# =============================================================== 判语与告警
def flags(f):
    """整体不合格的硬伤（不否决每一行，但引用数字前必须先看）。"""
    out = []
    sub, mw, sh = f["subtraction"], f["mw"], f["shape"]
    cq = sub.get("contrast_lowq")
    if cq is not None and cq < 0.02:
        out.append("低 q 对比度只有 %s（<2%%）" % _pct(cq))
    rg = f["guinier_adopted"]["_rg"] if f["guinier_adopted"] else None
    if rg and rg > 40:
        out.append("Rg ≈ %s Å 像多聚/聚集" % _fmt(rg, "", 0))
    gnom = f.get("gnom") or {}
    if gnom and "suspicious" in str(gnom.get("quality") or "").lower():
        out.append("GNOM 判语 %s（TE %s / α %s）" % (gnom.get("quality"),
                  _fmt(gnom.get("total_est"), "", 3), _fmt(gnom.get("alpha"), "", 0)))
    ad = f["ift_adopted"]
    if ad and ad.get("chisq") is not None and ad["chisq"] > 10:
        out.append("IFT χ² = %s（≫1）" % _fmt(ad["chisq"], "", 1))
    if f["guinier_adopted"] is None and f["guinier_rows"]:
        out.append("**没有 q 区间同时通过四条 Guinier 闸门**（见 tables/guinier_multi_range.csv）")
    if f.get("rec_verdict") == "rejected":
        out.append("**产物自己把推荐区间标成 NOT recommended**（%s）"
                   % str(f.get("rec_reason") or "闸门未全过")[:160])
    chis = [r["chisq"] for r in sh.get("dammif") or [] if r.get("chisq") is not None]
    if chis and min(chis) > 20:
        out.append("珠模 χ² = %s（≫5）" % _fmt(min(chis), "", 1))
    dv = sh.get("damaver")
    if isinstance(dv, dict) and len(dv.get("clusters") or []) > 1:
        out.append("DAMAVER 分成 %d 个 cluster（形状还没定）" % len(dv["clusters"]))
    if chis and 5 < min(chis) <= 20:
        out.append("珠模 χ² = %s 偏大 → 是数据问题不是参数问题，先修曲线" % _fmt(min(chis), "", 1))
    return out


def missing_items(f):
    """本次没做的 / 不能信的（每条都写清原因与下一步）。"""
    out, rel = [], set(f["files"])
    sub, sh, mw = f["subtraction"], f["shape"], f["mw"]
    if not f["normalization"]["on"]:
        out.append("**没有做逐帧归一化**：该系列没有线站逐帧 txt / 监视器，束流起伏留在曲线里 —— "
                   "低 q 与绝对刻度都别太当真（本次按 `--no-header-normalization` 跑）。")
    if not f["ift"]:
        out.append("**本次运行没有跑 IFT/P(r) 节点**：重跑时把 ift 加到 `--steps` 即可。")
    elif not f["ift_trusted"]:
        out.append("**没有可信的 P(r)**：各支 IFT 都没过闸门（见下表），按设计不报 P(r)、"
                   "不拿它建 3D；只留 IFT 解供目视。")
    msg = sh.get("shape_msg")
    if msg:
        out.append("**电子云（DENSS）这次失败了**：`%s` —— 常见于 IFT 解本身不稳或曲线在低 q/高 q "
                   "有坏点；先把前面的结果定下来，再单独试 `--denss-mode Slow`。" % str(msg)[:160])
    elif sh.get("denss") is None and sh.get("engine") not in (None, "none"):
        out.append("**没有出电子云**：`models/` 里没有 `*denss*.mrc`（要么没跑 shape 步，"
                   "要么 DENSS 失败，看 `models/*denss.log`）。")
    dmsg = sh.get("dammif_msg")
    if dmsg:
        out.append("**没有珠模（DAMMIF）**：`%s`。珠模需要 ATSAS（GNOM + DAMMIF + DAMAVER），"
                   "而且 DAMMIF 只吃 GNOM 的 IFT。" % str(dmsg)[:160])
    elif not sh.get("dammif") and sh.get("engine") == "none":
        out.append("**本次 `--model-engine none`**：没有做任何 3D 重建。")
    chis = [r["chisq"] for r in sh.get("dammif") or [] if r.get("chisq") is not None]
    if chis and min(chis) > 5:
        out.append("**珠模 χ² 偏大（最好 %s）**：这不是 DAMMIF 参数问题，是**曲线本身还不够干净**"
                   "——低 q 的束挡光晕/背景残留会把 χ² 顶到几百（本机实测：把起算 q 从 0.0064 提到 "
                   "0.0567 后 χ² 从 556 掉到 1.4）。回去修扣减与 q 窗口，别在 DAMMIF 里调参。"
                   % _fmt(min(chis), "", 1))
    bad = sh.get("dammif_bad") or []
    if bad:
        out.append("**%d 个珠模起跑失败**：`%s`（见 `tables/shape_results.json` 与 models/ 下的日志）。"
                   % (len(bad), str(bad[0].get("failed"))[:200]))
    if mw.get("Vp") and mw.get("Vc") and mw["Vp"].get("mw") and mw["Vc"].get("mw"):
        if mw["Vp"]["mw"] > 2 * mw["Vc"]["mw"]:
            out.append("**Vp 分子量是 Vc 的 %.1f 倍**：Porod 体积法高估 ~1.5× 是常态，但超过 2× "
                       "要怀疑聚集体或低 q 污染。" % (mw["Vp"]["mw"] / mw["Vc"]["mw"]))
    if f["mode"] == "sec":
        if not [p for p in rel if p.startswith("models/")]:
            out.append("**没有 `models/`（3D 模型）**：本次 `--model-engine none`，"
                       "想建模型就重跑时去掉这个参数（DENSS 几分钟，DAMMIF 需要 ATSAS）。")
    else:
        if not [p for p in rel if p.startswith("frames/")]:
            out.append("没有 `frames/`：没加 `--save-frames`（逐帧曲线默认不落盘）。")
    return out


# =============================================================== 章节
def _s0(f, keys):
    """0. 结论速览：结论先行 + 好/差两组（与 §1 同一张表，保证口径一致）。"""
    rows = param_rows(f)
    good = [r for r in rows if r[4].startswith("🟢")]
    bad = [r for r in rows if not r[4].startswith("🟢")]
    L = ["## 0. 结论速览", ""]
    fl = flags(f)
    if fl:
        L += ["> 🟡 **需要警惕**（不否定整份结果，但引用数字前先看这几条）："
              + "；".join(x.replace("**", "") for x in fl) + "。", ""]
    avail, facts = _availability(f)
    L.append("- **数据可用性**：%s" % avail)
    L.append("- **表现较好的参数**：")
    for r in (good or [("", "", "（本次没有环节处在可用档）", "", "", "")]):
        if good:
            L.append("  - 🟢 %s %s —— %s" % (r[0], r[1], _why(r)))
    if not good:
        L.append("  - （本次没有环节处在可用档）")
    L.append("- **表现较差的参数**：")
    if bad:
        for r in bad:
            L.append("  - %s %s %s —— %s" % (r[4].split(" ")[0], r[0], r[1], _why(r)))
    else:
        L.append("  - （本次所有环节都在可用档）")
    L.append("- **只想先看结论，打开这几个就够**：")
    for line in _first_files(f):
        L.append("  - %s" % line)
    if facts:
        L.append("")
        L.append("- 本次数据的事实：%s" % "；".join(facts) + "。")
    L.append("")
    keys.extend(_key_pairs(rows))
    return L


def _why(r, limit=72):
    """§0 里的短因由：只取「误差怎么读」里第一句（完整说明在 §1 的表里）。"""
    txt = r[5].split("｜", 1)[-1]
    for c in ("🟢 ", "🟡 ", "🔴 "):
        txt = txt.replace(c, "")
    txt = txt.strip()
    for sep in ("；", "。", "（"):
        if sep in txt and len(txt) > limit:
            txt = txt.split(sep, 1)[0]
    return txt if len(txt) <= limit else txt[:limit].rstrip() + "…"


def _availability(f):
    """数据可用性：三档结论 + 支撑事实（两种模式**同一口径**）。

    🟢 = 采用区间 R² ≥ 0.99、跨区间 Rg 不漂（≤5%）、且 IFT 有可信的一支；
    🟡 = 有可引用的 Rg，但上面任一条打折（区间间漂移 / R² 略低 / 产物标了 NOT recommended）；
    🔴 = 没有可引用的 Rg（全部未收敛或没有区间过闸门），或 R² < 0.95。
    """
    ad = f["guinier_adopted"]
    rsq = ad["_r2"] if ad else None
    spread = f["guinier_spread"]
    if ad is None:
        return "🔴 不可用/存疑（没有可引用的 Rg：Guinier 全部未收敛或没有区间过闸门）", []
    if rsq is not None and rsq < 0.95:
        tag = "🔴 不可用/存疑（采用区间 R² < 0.95）"
    elif rsq is not None and rsq >= 0.99 and (spread is None or spread <= 5) \
            and (f["ift_trusted"] or not f["ift"]):
        tag = "🟢 可用"
    elif rsq is not None and rsq >= 0.99 and (spread is None or spread <= 15):
        tag = ("🟡 有限可用（采用区间 R² = %s，但跨区间 Rg 漂 %.1f%%）"
               % (_fmt(rsq, "", 4), spread) if spread is not None
               else "🟡 有限可用（低 q 仍有可疑成分）")
    elif f["rec_verdict"] == "rejected":
        tag = ("🟡 有限可用（**Rg 不要直接引用**：产物把它标成 NOT recommended，区间之间 "
               "Rg 漂 %.1f%%）" % (spread or 0.0))
    elif rsq is not None and rsq >= 0.95:
        tag = "🟡 有限可用（低 q 仍有可疑成分）"
    else:
        tag = "🔴 不可用/存疑"
    facts = []
    if f.get("n_frames"):
        facts.append("**%s 帧**" % f["n_frames"])
    facts.append("归一化：%s" % ("已做逐帧归一化（%s）" % f["normalization"]["how"]
                             if f["normalization"]["on"] else "**未做**"))
    facts.append("Guinier 采用区间 `%s`（q %s–%s Å⁻¹），r² = %s"
                 % (ad["_label"], _fmt(ad.get("_qmin"), "", 4), _fmt(ad.get("_qmax"), "", 4),
                    _fmt(rsq, "", 4)))
    if spread is not None:
        facts.append("可用区间之间 Rg 相对差 %.1f%%" % spread)
    if f["ift_adopted"] and f["ift_adopted"].get("dmax") is not None:
        facts.append("IFT 采用 `%s`，Dmax = %s Å" % (f["ift_adopted"]["tag"],
                    _fmt(f["ift_adopted"]["dmax"], "", 1)))
    chis = [r["chisq"] for r in f["shape"].get("dammif") or [] if r.get("chisq") is not None]
    if chis:
        facts.append("珠模 χ² = %s" % _fmt(min(chis), "", 2))
    return tag, facts


def _first_files(f):
    """「先看这三个」：按模式给出最省事的三个入口（只列真实存在的）。"""
    rel = set(f["files"])
    outs = []
    if f["mode"] == "tube" and "qc.png" in rel:
        outs.append("**整体对不对？** → `qc.png`（四联诊断图：曲线 / Kratky / Rg / P(r)）")
    rep = sorted(p for p in rel if p.startswith("reports/") and p.endswith(".pdf"))
    if rep:
        outs.append("**RAW 的原始输出？** → `%s`（Guinier / IFT / 序列图都在里面）" % rep[0])
    sub = [p for p in rel if p.startswith("profiles/03_subtracted/") and p.endswith(".dat")]
    if sub:
        outs.append("**最终扣减曲线？** → `%s`（要拿去画图 / 喂别的软件用这个）" % sorted(sub)[0])
    if f["mode"] == "sec" and "series/series_plot.png" in rel:
        outs.append("**选帧选得对不对？** → `series/series_plot.png` + `series/ranges.json`")
    md = [p for p in rel if p.startswith("models/") and p.endswith((".mrc", ".cif"))]
    if md:
        outs.append("**3D 长什么样？** → `%s`" % sorted(md)[0])
    return outs[:3] or ["（本目录里没有可推荐的入口文件，见第 5 节）"]


def param_rows(f):
    """§1 的表：环节 / 参数 / 值 / 程序误差 / 能不能用 / 误差·不确定度怎么读。

    两种模式**同一套行、同一顺序**；有数据就出，没有就整行不出现。
    """
    rows, sub, mw, sh = [], f["subtraction"], f["mw"], f["shape"]
    ad = f["guinier_adopted"]
    cq = sub.get("contrast_lowq")
    weak = cq is not None and cq < 0.02
    # ---- Guinier
    if ad:
        gates = ad.get("gates") or {}
        passed = ad.get("_source") == "producer"
        split = _num(ad.get("_split_rel", ad.get("rg_split_rel")))
        ok = "🟢 可用（闸门全过）" if passed else "🟢 可用（按 qRg 与 r² 挑出的采用区间）"
        rows.append(("Guinier 拟合", "Rg", "%s Å" % _fmt(ad["_rg"], "", 2),
                     "± %s Å" % _fmt_err(ad.get("_rg_err") or ad.get("rg_err")), ok,
                     "程序误差是**区间内最小二乘标准误差**，区间收窄它必然变小，**不代表更准**；"
                     "真不确定度看**区间间漂移**（左右半段相对差 %s）与 R²。报 Rg 必须一起写 "
                     "q 区间 %s–%s Å⁻¹" % (_pct(split) if split is not None else "—",
                                          _fmt(ad.get("_qmin", ad.get("qmin")), "", 4),
                                          _fmt(ad.get("_qmax", ad.get("qmax")), "", 4))))
        _ie = ad.get("_i0_err", ad.get("i0_err"))
        rows.append(("Guinier 拟合", "I(0)",
                     _fmt(ad.get("_i0", ad.get("i0")), "", 1),
                     ("± %s" % _fmt_err(_ie)) if _num(_ie) is not None else "—",
                     "🟡 只作相对量" if weak else "🟢 可用于序列内比较",
                     "统计误差同上；**系统校验**看稀释序列：I(0) 应随浓度成比例，并与 Vc-MW 自洽。"
                     "浓度未知时它没有绝对意义"))
        rows.append(("Guinier 拟合", "R² / χ²_red / qRg",
                     "%s / %s / %s–%s" % (_fmt(ad.get("_r2"), "", 4),
                                          _fmt(_num(ad.get("chi2_red")), "", 3),
                                          _fmt(ad.get("_qrg_min"), "", 2),
                                          _fmt(ad.get("_qrg_max"), "", 2)),
                     "—", "🟢 区间形状自洽" if passed else "🟢 采用区间（R² 最高）",
                     "R²>0.99 且 χ²_red≈1 才说明这段真的是 Guinier 区；χ²_red≫1 = 误差被低估或"
                     "形状不对（SEC 的表里没有 χ²_red 一列时显示 —）"))
    elif f["guinier_auto"]:
        a = f["guinier_auto"]
        rows.append(("Guinier 拟合", "Rg（auto，仅参考）",
                     "%s Å" % _fmt(a.get("rg"), "", 2), "—", "🔴 没有区间通过闸门",
                     "auto 区间只是为了给 IFT 一个出发点；**本样品没有可引用的 Rg**。"
                     "见 `tables/guinier_multi_range.csv` 里每条闸门为什么没过"))
    # ---- IFT
    adf = f["ift_adopted"]
    if adf:
        de = adf.get("dmax_err")
        de_s = ("± %s Å" % _fmt_err(de)) if _num(de) is not None else "程序不给"
        rg_ref = ad["_rg"] if ad else None
        ratio = (adf["dmax"] / rg_ref) if (adf.get("dmax") and rg_ref) else None
        rows.append(("IFT / P(r)", "Dmax",
                     "%s Å（引擎 %s）" % (_fmt(adf.get("dmax"), "", 1), adf.get("engine")),
                     de_s, "🟢 可信（闸门全过）" if adf.get("trusted") is not False
                     else "🟡 条件可信",
                     "**真不确定度 = 换起算 q / Dmax 网格后它怎么变**（见第 5 节的全表）；"
                     "程序误差只反映单次拟合内部一致性。本次 Dmax/Rg = %s"
                     % (_fmt(ratio, "", 2) if ratio else "—")))
        if _num(adf.get("rg_real")) is not None:
            drift = (abs(adf["rg_real"] - rg_ref) / rg_ref * 100) if rg_ref else None
            rows.append(("IFT / P(r)", "Rg（实空间）",
                         "%s Å" % _fmt(adf.get("rg_real"), "", 2),
                         "± %s Å" % _fmt_err(adf.get("rg_err")),
                         "🟢 与 Guinier 对得上" if (drift is not None and drift <= 15)
                         else "🟡 与 Guinier 差得多",
                         ("和 Guinier 的 Rg 差 %.0f%%（<15%% 才算对得上）；差得多说明 q 窗口或 "
                          "Dmax 不对，而不是「更精确的 Rg」" % drift) if drift is not None
                         else "本样品没有可引用的 Guinier Rg，只能靠 Rg < Dmax/2 之类的自洽性打量"))
        if _num(adf.get("chisq")) is not None:
            ch = adf["chisq"]
            rows.append(("IFT / P(r)", "χ²（拟合优度）", _fmt(ch, "", 3), "—",
                         "🟢 合理" if 0.3 <= ch <= 3 else "🟡 偏离 1 较远",
                         "≈1 = 误差标定到位；≫1 → 误差被低估或形状不对；**<1 常见于误差被高估**，"
                         "不代表拟合更好"))
        if adf.get("total_est") is not None or adf.get("alpha") is not None:
            q = str(adf.get("quality") or "").lower()
            rows.append(("GNOM", "Total estimate / α",
                         "%s / %s" % (_fmt(adf.get("total_est"), "", 3), _fmt(adf.get("alpha"), "", 2)),
                         "—",
                         "🔴 **GNOM 判语：%s**" % adf.get("quality") if "suspicious" in q else
                         ("🟢 %s" % adf.get("quality") if "reasonable" in q
                          else "🟡 %s" % (adf.get("quality") or "无判语")),
                         "GNOM 自己的解质量指标（TE<1 一般算合理）+ 自动定的惩罚权重 α"
                         "（α 到几百以上说明解被「拽」得很凶）。它们只说明「这个解自洽」，"
                         "**不证明解唯一**"))
    # ---- 分子量
    vc, vp = mw.get("Vc") or {}, mw.get("Vp") or {}
    if vc.get("mw") is not None:
        rel = (abs(vp["mw"] - vc["mw"]) / vc["mw"] * 100) \
            if (_finite(vp.get("mw")) and _finite(vc.get("mw")) and vc["mw"]) else None
        rows.append(("分子量", "Vc（浓度无关）", "%s kDa" % _fmt_mw(vc["mw"]),
                     "± %s kDa" % _fmt_err(vc.get("err")) if vc.get("err") else "± ~10%",
                     "🟡 偏大需查（像聚集/多聚）"
                     if (weak or (ad and ad["_rg"] and ad["_rg"] > 40)) else "🟢 本次最可信的一条",
                     "误差 = RAW 的经验不确定度（Vcor 传播），**不含**经验系数的系统偏差；"
                     "**浓度无关**所以可直接用。与 Vp 差 %s；比单体预期大 2 倍以上多为聚集体"
                     % (_fmt(rel, "", 0) + "%" if rel is not None else "—")))
    if vp.get("mw") is not None:
        vol = vp.get("vol")
        rows.append(("分子量", "Vp（Porod 体积）", "%s kDa" % _fmt_mw(vp["mw"]), "—",
                     "🟡 只作数量级",
                     "**系统性高估 ~1.5×**是常态（已知偏差，不是误差）；修正后的 Porod 体积 "
                     "%s Å³ 可用来算 P/V 比值判断紧密程度" % _fmt(vol, "", 0)))
    bay = mw.get("Bayesian") or {}
    if bay.get("mw") is not None:
        ci = bay.get("ci") or (None, None)
        err = ("CI %s–%s kDa" % (_fmt(ci[0], "", 1), _fmt(ci[1], "", 1))
               if ci[0] is not None and ci[1] is not None else "—")
        cp = _num(bay.get("ci_prob"))
        if cp is not None and abs(cp) != float("inf"):
            err += "（区间外概率 %s）" % _fmt(100 - cp, "", 1)
        if not _finite(bay.get("mw")):
            err = "发散（程序没给有限值）"
        rows.append(("分子量", "DATMW（贝叶斯）", "%s kDa" % _fmt_mw(bay["mw"]), err,
                     "🟢 与 Vc 同量级才安心",
                     "ATSAS 用 P(r) + Vc 做的贝叶斯推断；曲线质量差时区间会张得很开 —— "
                     "**看区间宽度判可用性**，与 Vc 差得多说明 P(r) 不可信或体系不单一"))
    dc = mw.get("Datclass") or {}
    if dc.get("mw") is not None:
        sec = _num(dc.get("second"))
        sec_txt = ("%s Å" % _fmt(sec, "", 1)) if (_finite(sec) and sec > 0) \
            else "—（本次没给有效值）"
        rows.append(("分子量", "DATCLASS", "%s kDa" % _fmt_mw(dc["mw"]), "—",
                     "🟢 参考" if _finite(dc["mw"]) and _num(dc["mw"]) > 0 else "🔴 本次没给出结果",
                     "ATSAS 按曲线形状分类（本次判为 `%s`，第二返回值 %s）；"
                     "高 q 截得太短时会以 `insufficient data` 跳过"
                     % (dc.get("cls") or "?", sec_txt)))
    # ---- 珠模
    dm = sh.get("dammif") or []
    if dm:
        best = min(dm, key=lambda r: (r.get("chisq") if r.get("chisq") is not None else 9e9))
        ch = _num(best.get("chisq"))
        rows.append(("珠模 DAMMIF（%d 个）" % len(dm), "χ²", _fmt(ch, "", 2), "—",
                     "🟢 拟合好" if (ch is not None and ch <= 3) else
                     ("🟡 偏大 → 数据问题" if ch is not None and ch <= 5 else "🔴 不可用"),
                     "1–3 才算拟合好；>5 先修曲线（低 q 污染 / 背景残留），"
                     "**不是 DAMMIF 参数问题**。χ² 好也不代表解唯一"))
        rgp = (f["ift_adopted"] or {}).get("rg_real") or (ad["_rg"] if ad else None)
        drift = (abs(best["rg"] - rgp) / rgp * 100) if (best.get("rg") and rgp) else None
        rows.append(("珠模 DAMMIF", "Rg / Dmax / MW",
                     "%s / %s Å / %s kDa" % (_fmt(best.get("rg"), "", 1),
                                             _fmt(best.get("dmax"), "", 1),
                                             _fmt(best.get("mw"), "", 0)),
                     "—（系综内看 SD）",
                     "🟢 与 P(r) 一致" if (drift is not None and drift <= 10)
                     else "🟡 与 P(r) 不一致",
                     "模型 Rg 应与 P(r) 的 Rg 接近（本次差 %s）；Dmax 一般比 P(r) 的略大（珠子外壳）"
                     % (_fmt(drift, "%", 0) if drift is not None else "—")))
    dv = sh.get("damaver")
    if isinstance(dv, dict) and dv.get("mean_nsd") is not None:
        ncl = len(dv.get("clusters") or []) or 1
        nsd = _num(dv["mean_nsd"])
        rows.append(("珠模一致性 DAMAVER", "平均 NSD", _fmt(nsd, "", 3),
                     "± %s" % _fmt_err(dv.get("stdev_nsd")),
                     ("🟡 分簇（形状未定）" if ncl > 1 else
                      ("🟢 同一类形状" if nsd is not None and nsd <= 2 else
                       ("🟡 尚可" if nsd is not None and nsd <= 3 else "🔴 形状未定"))),
                     "这是**模型两两之间的差异**，**不是**与实验数据的拟合误差；≤2 一般算同一类形状。"
                     "本次 %d 个 cluster%s" % (ncl, "（>1 = 模型没收敛到同一形状）" if ncl > 1 else "")))
    # ---- 电子云
    de = sh.get("denss")
    if isinstance(de, dict) and de.get("chi2") is not None:
        rgp = (f["ift_adopted"] or {}).get("rg_real") or (ad["_rg"] if ad else None)
        d2 = (abs(de["rg_model"] - rgp) / rgp * 100) if (de.get("rg_model") and rgp) else None
        ch = _num(de["chi2"])
        rows.append(("电子云 DENSS", "χ² / 模型 Rg",
                     "%s / %s Å" % (_fmt(ch, "", 2), _fmt(de.get("rg_model"), "", 1)), "—",
                     "🔴 无意义（IFT 不可信）" if not f["ift_trusted"] else
                     ("🟢 与 P(r) 一致" if (d2 is not None and d2 <= 10) else "🟡 与 P(r) 不一致"),
                     "模型 Rg 应与 P(r) 的 Rg 差 <10%%（本次差 %s）；χ² 只在 IFT 可信时才有意义。"
                     "support %s Å³、盒子 %s Å"
                     % (_fmt(d2, "%", 0) if d2 is not None else "—",
                        _fmt(de.get("support_volume"), "", 0), _fmt(de.get("side"), "", 1))))
    return rows


def _key_pairs(rows):
    """验收用：README 里必须出现的「关键数字」原文（label, 显示串）。"""
    out = []
    for r in rows:
        if r[2] not in ("", "—", None):
            out.append(("%s %s" % (r[0], r[1]), str(r[2])))
    return out


def _s1(f, keys):
    rows = param_rows(f)
    L = ["## 1. 最终拟合参数（明细表）", "",
         "> **怎么读**：先看「能不能用」这一列（🟢 好 / 🟡 需要警惕 / 🔴 坏）—— 它综合了闸门"
         "（值本身自不自洽）、跨引擎/跨起跑点一致性、以及与别的独立判据（稀释序列、MW）对不对得上。",
         "> 「程序误差」一列只是**拟合随机误差**，它**永远不会**告诉你背景扣错了、低 q 被束挡污染了、"
         "或样品是多分散的 —— 那类**系统误差**写在最后一列。", "",
         "| 环节 | 参数 | 值 | 程序误差 | 能不能用 | 误差 / 不确定度怎么读 |",
         "|---|---|---|---|---|---|"]
    for r in rows:
        L.append("| " + " | ".join(_cell(x) for x in r) + " |")
    L.append("")
    return L


def _s2(f, keys):
    """2. 关键结果：一行一个关键结论，每行带「怎么看」（两种模式同一顺序）。"""
    sub, mw, sh = f["subtraction"], f["mw"], f["shape"]
    ad, adf = f["guinier_adopted"], f["ift_adopted"]
    L = ["## 2. 关键结果", "", "| 项目 | 数值 | 怎么看 |", "|---|---|---|"]
    rows = []

    def add(item, value, how):
        rows.append((item, str(value), how))

    # 数据 / 归一化
    if f["mode"] == "tube":
        ctrl = "、".join(sub.get("control_runs") or []) or "见 norm/frame_qc.csv"
        add("用哪些 run 当对照", "%s（样品 + 对照共 %s 帧）" % (ctrl, _fmt(f.get("n_frames"), "", 0)),
            "文件名前缀 = 样品名的是样品，其余都是对照；每次运行前先确认")
        add("高 q 缩放因子", "%s（在 q %s–%s Å⁻¹ 定标）"
            % (_fmt(sub.get("control_scale_factor"), "", 4),
               _fmt((sub.get("scale_window") or [None, None])[0], "", 2),
               _fmt((sub.get("scale_window") or [None, None])[1], "", 2)),
            "偏离 1 超过 5% 说明对照与样品不匹配，扣减后会有常数残留")
    else:
        add("扣减用的帧区间（0 基）", "buffer %s；样品 %s"
            % ("、".join("%d–%d" % (a, b) for a, b in (sub.get("buffer") or [])) or "—",
               "、".join("%d–%d" % (a, b) for a, b in (sub.get("sample") or [])) or "—"),
            "峰前 + 峰后两段 buffer 能吃掉基线漂移；只用一段时峰后的基线常不等于峰前")
    add("归一化", "已做逐帧归一化（%s）" % f["normalization"]["how"] if f["normalization"]["on"]
        else "**未做**（该系列没有线站 txt / 监视器）",
        "归一化只含束流起伏，**不含几何漂移**（后者看 video/*.centroid.csv）"
        if f["normalization"]["on"] else "束流起伏留在曲线里；低 q 与绝对刻度都不可信")
    if sub.get("contrast_lowq") is not None:
        add("低 q 对比度", "%s（扣减曲线 q_min 处 ÷ 样品同点）" % _pct(sub["contrast_lowq"]),
            "<2% 时低 q 不可信；>5% 才算干净")
    if sub.get("qmin") is not None:
        add("分析窗", "q %s–%s Å⁻¹（%s 点）"
            % (_fmt(sub["qmin"], "", 4), _fmt(sub["qmax"], "", 4), _fmt(sub.get("n_points"), "", 0)),
            "管式按 SNR≥2 定上界；SEC 用 I/σ 掉到 ~2 的位置截高 q。低 q 落在这里的"
            "束挡光晕要单独判（--qmin）")
    # Guinier
    a = f["guinier_auto"]
    if a:
        add("auto-Guinier（仅参考）",
            "Rg = %s Å（q %s–%s，R² = %s）"
            % (_fmt(a.get("rg"), "", 1), _fmt(a.get("qmin"), "", 4), _fmt(a.get("qmax"), "", 4),
               _fmt(a.get("r_sqr") if a.get("r_sqr") is not None else a.get("r2"), "", 3)),
            "RAW 自动选的区间，**不要直接引用**")
    if ad:
        gates = ad.get("gates") or {}
        bad = [k for k, v in gates.items() if v is False]
        add("采用的 Guinier 区间",
            "q %s–%s Å⁻¹：Rg = %s ± %s Å，qRg ≤ %s，R² = %s"
            % (_fmt(ad.get("_qmin", ad.get("qmin")), "", 4),
               _fmt(ad.get("_qmax", ad.get("qmax")), "", 4),
               _fmt(ad["_rg"], "", 2), _fmt(ad.get("_rg_err", ad.get("rg_err")), "", 2),
               _fmt(ad.get("_qrg_max", ad.get("qrg_max")), "", 2), _fmt(ad["_r2"], "", 4)),
            ("这是本流程推荐引用的 Rg；报告里**必须一起写 q 区间**。"
             + ("未过闸门：%s" % "、".join(bad) if bad else "闸门全过"))
            if ad.get("_source") == "producer" else
            ("**产物没有给出可采用的推荐区间**"
             + ("（它把候选区间标成 NOT recommended：%s）" % str(f.get("rec_reason"))[:150]
                if f.get("rec_verdict") == "rejected" else "")
             + "，本流程按「收敛 + qRg 落在 [0.28, 1.35] + R² 最高」挑的；"
               "报告里**必须一起写 q 区间**，并说明 Rg 随区间会漂多少（见下一行）"))
    else:
        add("采用的 Guinier 区间", "**没有区间通过全部闸门**",
            "低 q 不服从单一 Guinier 定律，Rg 只能给范围，见 `tables/guinier_multi_range.csv`")
    if f["guinier_spread"] is not None or f["guinier_excluded"]:
        used = f["guinier_used"]
        lo_hi = f.get("guinier_used_rg")
        used_txt = ("、".join(used) if len(used) <= 4 else
                    "%d 个区间（Rg %s–%s Å）" % (len(used), _fmt(lo_hi[0], "", 1),
                                                _fmt(lo_hi[1], "", 1)) if lo_hi else
                    "%d 个区间" % len(used))
        add("跨区间一致性", "可用区间之间 Rg 相对差 %s（参与：%s）%s"
            % (_fmt(f["guinier_spread"], "%", 1), used_txt or "—",
               ("；**排除**：" + "、".join(f["guinier_excluded"])) if f["guinier_excluded"] else ""),
            "只统计收敛且 R² ≥ 0.9 的区间 —— 未收敛（哨兵 -1）与明显很差的区间本来就不能引用，"
            "让它们参与一致性判据等于一票否决")
    # IFT
    if adf:
        add("P(r) / IFT 采用的那一支",
            "引擎/起跑点 `%s`：Dmax = %s Å，Rg(实空间) = %s Å，chisq = %s"
            % (adf.get("tag"), _fmt(adf.get("dmax"), "", 1), _fmt(adf.get("rg_real"), "", 1),
               _fmt(adf.get("chisq"), "", 2)),
            "三闸门（Rg 对得上 / Dmax ≤ 4.5·Rg / Dmax 未越出搜索网格）"
            + ("全过才算结果" if adf.get("trusted") is not False
               else "**未全过，不要引用**（见第 5 节全表）"))
    elif f["ift"]:
        add("P(r) / IFT", "**没有可信的一支**（见第 5 节全表）",
            "按设计不报 P(r)、不拿它建 3D")
    else:
        add("P(r) / IFT", "**本次没跑**（`--steps` 里没包含 ift）", "重跑时加上 ift 节点")
    # MW
    for name, key, how in (("Vc 法·浓度无关", "Vc", "**这个更可信**；比单体预期大 2 倍以上多为聚集体"),
                           ("Vp·Porod 体积", "Vp", "系统性高估 ~1.5×，只作数量级参考"),
                           ("DATMW 贝叶斯", "Bayesian", "看置信区间宽度判可用性"),
                           ("DATCLASS", "Datclass", "同时给曲线形状分类；高 q 截太短时会被跳过")):
        d = mw.get(key) or {}
        if d.get("mw") is not None:
            add("分子量（%s）" % name, "%s kDa" % _fmt_mw(d["mw"]), how)
    # 3D
    de = sh.get("denss")
    if isinstance(de, dict):
        add("电子云（DENSS）", "χ² = %s，模型 Rg = %s Å，support 体积 = %s Å³"
            % (_fmt(de.get("chi2"), "", 2), _fmt(de.get("rg_model"), "", 1),
               _fmt(de.get("support_volume"), "", 0)),
            "模型 Rg 应与 P(r) 的 Rg 接近；χ² 只在 IFT 可信时才有意义")
    elif sh.get("shape_msg"):
        add("电子云（DENSS）", "**没出**：%s" % str(sh["shape_msg"])[:120], "见第 5 节")
    elif sh.get("engine") == "none":
        add("电子云（DENSS）", "**本次没跑**", "重跑时去掉 `--model-engine none`")
    dm = sh.get("dammif") or []
    if dm:
        best = min(dm, key=lambda r: (r.get("chisq") if r.get("chisq") is not None else 9e9))
        rgp = (adf or {}).get("rg_real")
        drift = ("；与 P(r) 的 Rg 差 %s%%"
                 % _fmt(abs(best["rg"] - rgp) / rgp * 100, "", 0)) if (rgp and best.get("rg")) else ""
        add("珠模（DAMMIF ×%d）" % len(dm),
            "最好的 χ² = %s；Rg = %s Å；Dmax = %s Å；MW = %s kDa%s"
            % (_fmt(best.get("chisq"), "", 1), _fmt(best.get("rg"), "", 1),
               _fmt(best.get("dmax"), "", 1), _fmt(best.get("mw"), "", 0), drift),
            "珠模 χ² **1–3 才算拟合好**；>5 说明曲线（尤其低 q）还不干净，先修数据")
    elif sh.get("dammif_msg"):
        add("珠模（DAMMIF）", "**没出**：%s" % str(sh["dammif_msg"])[:120], "见第 5 节")
    dv = sh.get("damaver")
    if isinstance(dv, dict) and dv.get("mean_nsd") is not None:
        cl = dv.get("clusters") or []
        add("珠模一致性（DAMAVER，%s 个模型）" % _fmt(dv.get("n_models") or len(dm), "", 0),
            "平均 NSD = %s ± %s；%d 个 cluster%s"
            % (_fmt(dv.get("mean_nsd"), "", 3), _fmt(dv.get("stdev_nsd"), "", 3), len(cl),
               "；代表模型 `%s`" % dv["rep_model"] if dv.get("rep_model") else ""),
            "NSD ≲2 说明这些模型是同一类形状；**分成多个 cluster = 形状还没定**")

    for item, value, how in rows:
        L.append("| %s | %s | %s |" % (_cell(item), _cell(value), _cell(how)))
    L.append("")
    keys.extend([("%s（§2）" % item, value) for item, value, _ in rows
                 if value not in ("", "—") and "没出" not in value and "没跑" not in value])
    return L


def _s3(f, keys):
    """3. 文件清单：无序列表（memo 式条目一律用列表，不写成段落）。"""
    rel = f["files"]
    used, L = set(), ["## 3. 这个文件夹里有什么（按建议阅读顺序）", ""]
    for mode, pat, desc, when in FILE_DOC:
        if mode not in ("*", f["mode"]):
            continue
        if pat.startswith("*"):
            hits = sorted(p for p in rel if p.endswith(pat.lstrip("*")))
        elif "*" in pat:
            import fnmatch
            hits = sorted(p for p in rel if fnmatch.fnmatch(p, pat))
        else:
            hits = [pat] if pat in rel else []
        if not hits:
            continue
        shown = hits[0]
        more = "，共 %d 个同名文件" % len(hits) if len(hits) > 1 else ""
        used.update(hits)
        L.append("- `%s` —— %s%s（%s）" % (shown, desc, more, when))
    orphan = [p for p in rel if p not in used and not p.endswith("README.md")]
    if orphan:
        L += ["", "> 还有 %d 个文件没在上面的清单里（多为各节点的中间产物）：%s%s。"
              % (len(orphan), "、".join("`%s`" % p for p in orphan[:6]),
                 " 等" if len(orphan) > 6 else "")]
    L.append("")
    return L


def _s4(f, keys):
    L = ["## 4. 每个数字的判据", "",
         "| 数字 | 好 | 勉强 | 差 / 不可信 |", "|---|---|---|---|"]
    for r in RULES:
        L.append("| " + " | ".join(_cell(x) for x in r) + " |")
    L.append("")
    return L


def _s5(f, keys):
    L = ["## 5. 这次没做的 / 不能信的", ""]
    mi = missing_items(f)
    if mi:
        L += ["- %s" % x for x in mi]
    else:
        L.append("- （没有跳过任何节点）")
    if f["ift"]:
        L += ["", "IFT 各起跑点/引擎的全部结果：", "",
              "| 起跑点 | 引擎 | q 范围 | Dmax (Å) | Rg_real (Å) | chisq | Rg 闸门 | Dmax/Rg 闸门 | 网格闸门 | 可信 |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for r in f["ift"]:
            if "failed" in r:
                L.append("| %s | %s | — | — | — | — | — | — | — | 失败：%s |"
                         % (_cell(r.get("tag")), _cell(r.get("engine")),
                            _cell(str(r.get("failed"))[:80])))
                continue
            g = r.get("gates") or {}

            def mark(k):
                v = g.get(k)
                return "—" if v is None else ("过" if v else "**不过**")
            qr = ("%s–%s" % (_fmt(r["qmin"], "", 4), _fmt(r["qmax"], "", 4))
                  if (r["qmin"] is not None and r["qmax"] is not None) else "—")
            L.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |"
                     % (_cell(r["tag"]), _cell(r["engine"]), qr,
                        _fmt(r["dmax"], "", 1),
                        _fmt(r["rg_real"], "", 1), _fmt(r["chisq"], "", 2),
                        mark("rg_vs_guinier"), mark("dmax_over_rg"), mark("dmax_within_grid"),
                        ("是" if r.get("trusted") is True else
                         ("**否**" if r.get("trusted") is False else "—"))))
        L += ["", "（BIFT 的 Dmax 只是搜索网格，优化步可以跑出去；跑到几百 Å 还配一个很小的 chisq "
              "就是**假解**。GNOM 的 Dmax 是它自己定的，没有误差棒。）"]
    L.append("")
    return L


def _s6(f, keys):
    L = ["## 6. 本次用的参数（复现用）", ""]
    cmd = f.get("command")
    if f["mode"] == "tube":
        L += ["```", (cmd or "%s --sample-dir %s --cfg %s --out-dir %s"
                      % (PRODUCER["tube"], f.get("input") or "<帧目录>",
                         f.get("cfg") or "<设置档>", f["out"])), "```", ""]
        L.append("> ATSAS：`%s`；逐帧归一化由 RAW 读 BL19U2 header txt 完成（不是自写算法）。"
                 % (f.get("atsas_dir") or "未指定（本次无 GNOM/DAMMIF/DAMAVER）"))
    else:
        a = f.get("args") or {}
        lines = ["%s \\" % PRODUCER["sec"],
                 "  --series-dir %s \\" % (f.get("input") or "<帧目录>"),
                 "  --out-dir %s \\" % f["out"],
                 "  --cfg %s \\" % (f.get("cfg") or "<设置档>")]
        for key, flag in (("buffer_range", "--buffer-range"), ("sample_range", "--sample-range"),
                          ("baseline", "--baseline"), ("trim_qmin", "--trim-qmin"),
                          ("trim_qmax", "--trim-qmax"), ("ift_dmax", "--ift-dmax"),
                          ("atsas_dir", "--atsas-dir"), ("guinier_ranges", "--guinier-ranges")):
            v = a.get(key)
            if v not in (None, "", "none"):
                lines.append('  %s "%s" \\' % (flag, v) if "range" in key else
                             "  %s %s \\" % (flag, v))
        if a.get("no_header_normalization"):
            lines.append("  --no-header-normalization \\")
        lines.append("  --steps %s \\" % ",".join(f.get("steps") or []))
        lines.append("  --model-engine %s --n-models %s"
                     % (a.get("model_engine", "auto"), a.get("n_models", 4)))
        L += ["```", "\n".join(lines), "```", ""]
        st = f.get("raw_settings") or {}
        L.append("> RAW 端设置：`ImageHdrFormat=%s`，`EnableNormalization=%s`，"
                 "`NormalizationList=%s`；ATSAS：`%s`。"
                 % (st.get("ImageHdrFormat", "None"), st.get("EnableNormalization", False),
                    st.get("NormalizationList", "—"), st.get("ATSSDir") or st.get("ATSASDir")
                    or "未指定"))
    L.append("")
    return L


def _s7(f, keys):
    rel = set(f["files"])
    ws = sorted(p for p in rel if p.endswith(".hdf5"))
    L = ["## 7. 想自己复核", ""]
    if ws:
        L.append("- **在 RAW 里交互式看**：`File → Open Workspace` 打开 `%s`，"
                 "曲线、Guinier、IFT 都在里面。" % ws[0])
    tables = sorted(p for p in rel if p.startswith("tables/") and p.endswith((".csv", ".json")))
    if tables:
        L.append("- **核对本文件的数字**：`%s`，机器可读版是 `%s`。"
                 % ("`、`".join(tables[:4]), "summary.json" if f["mode"] == "tube"
                    else "run_meta.json"))
    L.append("- **重跑**：见第 6 节的命令。本 README 由 `%s/scripts/write-readme.py` 生成"
             "（只读产物、不重算），所以写完要和 `%s` 的时间戳对一下 —— 避免出现「README 比数据旧」。"
             % (SKILL_NAME, "summary.json" if f["mode"] == "tube" else "run_meta.json"))
    L.append("- **低 q 到底能不能信 / 上翘是真还是假**：走 `assess-saxs-raw-data-quality`"
             "（空白−空白对照、背景形状失配、2D 差分三道检验）。")
    L.append("- **判据细节**：Rg 取点 → `assess-guinier-fit-quality`；P(r)/Dmax → "
             "`compute-and-validate-p-of-r`；MW 方法 → `choose-a-molecular-weight-method`；"
             "重建评估 → `evaluate-a-shape-reconstruction`。")
    L += ["", "### 这些缩写", ""]
    L += ["- %s" % x for x in ABBREV]
    L += ["", "### 目录树（两层）", "", "```"]
    L += _tree_lines(f["out"], depth=2)
    L += ["```", ""]
    return L


def _tree_lines(out, depth=2):
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
        ex = sorted({os.path.splitext(f)[1] for f in files if not f.startswith(".")} - {""})
        kinds = ", ".join(ex) or "—"
        lines.append("%s/  (%d 个文件；%s)" % (rel, n, kinds))
    return lines or ["(空)"]


def build(f):
    """返回 (markdown 行, 关键数字列表)。关键数字 = README 里必须出现的原文。"""
    keys = []
    L = ["# %s —— %s 处理结果" % (f["sample"], f["mode_label"]), "",
         "> 自动生成 %s ｜ 产生它的脚本 `%s` ｜ 原始帧 `%s` ｜ 配置 `%s`"
         % (f["generated_at"], f["producer"], f.get("input") or "?", f.get("cfg") or "?"), "",
         "> 这份 README 由 `%s` 的 `write-readme.py` 从**产物本身**读出来（不重新拟合、不臆造数字）；"
         "每个数字都能在 `tables/` / `%s` / `models/` 里找到原文。"
         "SEC 与管式两条流水线共用同一套章节与表格，验收命令："
         "`verify-results-folder.py <目录>`。" % (SKILL_NAME,
                                             "summary.json" if f["mode"] == "tube"
                                             else "run_meta.json"), ""]
    L += _s0(f, keys)
    L += _s1(f, keys)
    L += _s2(f, keys)
    L += _s3(f, keys)
    L += _s4(f, keys)
    L += _s5(f, keys)
    L += _s6(f, keys)
    L += _s7(f, keys)
    return L, keys


def render(out, mode=None):
    f = collect(out, mode)
    if f is None:
        return None, None, None
    lines, keys = build(f)
    return "\n".join(lines) + "\n", keys, f


def write_readme(out, mode=None, dry=False):
    """写 <out>/README.md；返回 (path|None, keys, facts)。"""
    text, keys, f = render(out, mode)
    if text is None:
        return None, None, None
    path = os.path.join(os.path.abspath(out), "README.md")
    if not dry:
        with open(path, "w") as fh:
            fh.write(text)
    return path, keys, f


# =============================================================== 验收
def structure_check(out, mode=None):
    mode = mode or detect_mode(out)
    if mode is None:
        return None, [dict(level="must", path="run_meta.json 或 summary.json", status="🔴",
                           note="两个标记文件都没有：这不像本流程产出的结果目录")]
    rows = []
    for level, pat, note in STRUCTURE[mode]:
        p = os.path.join(out, pat)
        if os.path.exists(p):
            if os.path.isdir(p):
                n = len([x for x in glob.glob(os.path.join(p, "*")) if not os.path.basename(x)
                         .startswith(".")])
                st = "🟢" if n else ("🔴" if level == "must" else "🟡")
                extra = "%d 个文件" % n if n else "空目录（本次没落盘）"
            else:
                st, extra = "🟢", os.path.getsize(p) and "存在" or "存在（0 字节）"
        else:
            st, extra = ("🔴" if level == "must" else "🟡"), {
                "must": "缺失（交付不完整）", "should": "缺失（第 5 节要写明原因）",
                "info": "本次没有（可接受）"}[level]
        rows.append(dict(level=level, path=pat, status=st, note=("%s；%s" % (extra, note))))
    return mode, rows


def _ndec(txt):
    """CSV 文本里的小数位数（用来判断"两个数在这个精度下是否相等"）。"""
    t = str(txt).strip()
    if not t or any(c in t.lower() for c in ("e", "n")) or "." not in t:
        return None
    return len(t.split(".")[1])


def _same(a, b, b_text=None):
    """同一物理量的两个读数是否一致。CSV 里是四舍五入后的文本 —— 按它的精度比，
    否则 48.41485（json）与 48.41（csv）会被误判成"自相矛盾"。"""
    a, b = _num(a), _num(b)
    if a is None or b is None:
        return True
    nd = _ndec(b_text) if b_text is not None else None
    if nd is not None:
        return round(a, nd) == round(b, nd)
    return abs(a - b) <= 1e-9 * max(1.0, abs(a), abs(b))


def cross_source_checks(out, mode, f):
    """两条独立路径读同一个数：不一致 = 产物自相矛盾。"""
    res = []
    if mode == "sec":
        meta = _read_json(os.path.join(out, "run_meta.json"), {}) or {}
        g = _read_csv(os.path.join(out, "tables", "guinier_multi_range.csv"))
        grow = meta.get("guinier_rows") or []
        if g and grow:
            n = min(len(g), len(grow))
            bad = []
            for i in range(n):
                for j, k in ((1, "rg"), (2, "i0"), (5, "qmin"), (6, "qmax"), (9, "r2")):
                    b = field(g[i], k)
                    if not _same(grow[i][j], b, b):
                        bad.append("%s 列 %s：run_meta=%s vs csv=%s" % (grow[i][0], k,
                                                                       grow[i][j], b))
            res.append(("run_meta.json ↔ tables/guinier_multi_range.csv",
                        "🟢" if not bad else "🔴", "%d 行逐列一致" % n if not bad
                        else "；".join(bad[:3])))
        i_g = _read_csv(os.path.join(out, "tables", "ift_summary.csv"))
        irow = meta.get("ift_rows") or []
        if i_g and irow:
            bad = []
            for i, r in enumerate(irow):
                if i >= len(i_g):
                    break
                for j, k in ((1, "dmax"), (3, "rg_real"), (5, "chisq")):
                    b = field(i_g[i], k)
                    if not _same(r[j], b, b):
                        bad.append("%s 的 %s：run_meta=%s vs csv=%s" % (r[0], k, r[j], b))
            res.append(("run_meta.json ↔ tables/ift_summary.csv",
                        "🟢" if not bad else "🔴", "%d 行一致" % min(len(i_g), len(irow))
                        if not bad else "；".join(bad[:3])))
    else:
        sm = _read_json(os.path.join(out, "summary.json"), {}) or {}
        g = _read_csv(os.path.join(out, "tables", "guinier_multi_range.csv"))
        rec = sm.get("guinier_recommended")
        if isinstance(rec, dict) and g:
            cands = [r for r in g
                     if _same(field(r, "qmin"), rec.get("qmin"), field(r, "qmin"))
                     and (rec.get("qmax") is None
                          or _same(field(r, "qmax"), rec.get("qmax"), field(r, "qmax")))]
            hit = min(cands, key=lambda r: abs((_num(field(r, "rg")) or 0)
                                               - (_num(rec.get("rg")) or 0))) if cands else None
            if hit is None:
                res.append(("summary.json ↔ tables/guinier_multi_range.csv", "🟡",
                            "推荐区间的 qmin 在表里没找到同一条（可能被重算过）"))
            else:
                ok = _same(rec.get("rg"), field(hit, "rg"), field(hit, "rg"))
                res.append(("summary.json ↔ tables/guinier_multi_range.csv",
                            "🟢" if ok else "🔴",
                            "推荐区间 Rg 一致（按 csv 的精度比）" if ok
                            else "推荐区间 Rg 不一致：json=%s vs csv=%s"
                                 % (rec.get("rg"), field(hit, "rg"))))
        ift = (sm.get("ift") or {}).get("runs") or []
        ic = _read_csv(os.path.join(out, "tables", "ift_summary.csv"))
        if ift and ic:
            bad = []
            for r in ift:
                hit = next((x for x in ic if x.get("tag") == r.get("tag")), None)
                if hit is None:
                    bad.append("`%s` 在 csv 里没有" % r.get("tag"))
                    continue
                for k, j in (("dmax", "dmax"), ("rg_realspace", "rg_real"), ("chisq", "chisq")):
                    b = field(hit, j)
                    if not _same(r.get(k), b, b):
                        bad.append("`%s` 的 %s：json=%s vs csv=%s" % (r.get("tag"), k,
                                                                     r.get(k), b))
            res.append(("summary.json ↔ tables/ift_summary.csv",
                        "🟢" if not bad else "🔴",
                        "%d 支一致" % len(ift) if not bad else "；".join(bad[:3])))
        mw_json = sm.get("mw") or {}
        mw_csv = {str(field(r, "method") or ""): _num(field(r, "mw"))
                  for r in _read_mw_csv(os.path.join(out, "tables", "mw.csv"))}
        if isinstance(mw_json, dict) and mw_csv:
            bad = []
            for k, v in mw_json.items():
                if not isinstance(v, dict) or v.get("mw") is None:
                    continue
                got = mw_csv.get(k)
                if got is None:
                    continue
                if not _same(v["mw"], got, got):
                    bad.append("%s：json=%s vs csv=%s" % (k, v["mw"], got))
            res.append(("summary.json ↔ tables/mw.csv", "🟢" if not bad else "🔴",
                        "%d 法一致" % len(mw_json) if not bad else "；".join(bad[:3])))
    return res


def staleness_check(out):
    """README 是不是比数据旧：**管线最后一步才写 README**，所以任何关键表比它新 ⇒ 这张表被改过。

    关键表（`tables/*.csv`、`summary.json` / `run_meta.json`）比 README 新 = 🔴（先重写 README 再说）；
    其它产物（图、视频、日志）比 README 新 = 🟡（不影响数字，但也提示重跑一次）。
    """
    rp = os.path.join(out, "README.md")
    if not os.path.exists(rp):
        return ("README 与数据的时间戳", "🔴", "没有 README.md")
    t_readme = os.path.getmtime(rp)
    key, other = [], []
    for pat in ("tables/*.csv", "tables/*.json"):
        key += glob.glob(os.path.join(out, pat))
    for nm in ("summary.json", "run_meta.json"):
        if os.path.exists(os.path.join(out, nm)):
            key.append(os.path.join(out, nm))
    for pat in ("models/*", "qc.png", "reports/*", "ifts/*", "profiles/04_sample/*.dat",
                "profiles/06_guinier/*", "series/*"):
        other += glob.glob(os.path.join(out, pat))
    newer_key = [(os.path.relpath(p, out), os.path.getmtime(p) - t_readme)
                 for p in key if os.path.getmtime(p) > t_readme + 60]
    newer_other = [(os.path.relpath(p, out), os.path.getmtime(p) - t_readme)
                   for p in other if os.path.getmtime(p) > t_readme + 60]
    if newer_key:
        newer_key.sort(key=lambda t: -t[1])
        return ("README 与数据的时间戳", "🔴",
                "关键表比 README 新 %.0f 分钟：`%s` —— 先跑 write-readme.py 重写（README 里的数字"
                "可能是旧版本）" % (newer_key[0][1] / 60, newer_key[0][0]))
    if newer_other:
        newer_other.sort(key=lambda t: -t[1])
        return ("README 与数据的时间戳", "🟡",
                "有 %d 个非关键产物比 README 新（如 `%s`）；数字不受影响，但可重跑一次让它同步"
                % (len(newer_other), newer_other[0][0]))
    return ("README 与数据的时间戳", "🟢", "README 比所有产物新（= 管线跑完的结果）")


def readme_check(out, text, keys, f):
    """README 本身：章节齐不齐、关键数字在不在。"""
    res = []
    heads = [l.strip() for l in text.splitlines() if l.startswith("## ")]
    missing = [h for h in SECTION_ORDER if h not in heads]
    ordered = [h for h in heads if h in SECTION_ORDER]
    res.append(("README 章节（%d 节）" % len(SECTION_ORDER), "🟢" if not missing else "🔴",
                "顺序与规范一致" if not missing and ordered == SECTION_ORDER
                else ("缺：" + "、".join(missing) if missing
                      else "顺序不对：%s" % "、".join(ordered))))
    absent = [("%s = %s" % (lab, val)) for lab, val in keys if str(val) not in text]
    res.append(("README ↔ 产物（关键数字 %d 条）" % len(keys),
                "🟢" if not absent else "🔴",
                "全部能在 README 里找到原文" if not absent
                else "有 %d 条在 README 里找不到：%s" % (len(absent), "；".join(absent[:4]))))
    doc = [p for p in f["files"] if "`%s`" % p in text]
    res.append(("README ↔ 目录清单", "🟢" if doc else "🟡",
                "清单里提到 %d 个真实文件" % len(doc)))
    res.append(staleness_check(out))
    return res
