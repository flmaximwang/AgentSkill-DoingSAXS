#!/usr/bin/env python
"""Assess the quality of the RAW frames of a batch SAXS dataset - before any fit.

For every sample directory under <root> it integrates the frames exactly the way the
production pipeline does (same RAW settings, forced transmission normalization) and looks
at the frame-level invariants that decide whether any downstream fit can be trusted:

  信号     低 q 对比度（样品 vs 自己的 control）、单帧信噪比与可用 q 上限、误差是否诚实
  稳定性   20 帧一个 run 内的电平漂移、帧间起伏、离群帧（低 q 与 q=0.2 两套判据）
  背景     control 电平、run 内漂移（毛细管污染）、run 间一致
  形状提示 I(q->0)/I(0.2)、扣减后低 q 上翘（真散射还是假信号见 attribute-lowq-upturn.py）

输出（默认 <root>/_rawqc）：
    raw_quality.md / raw_quality.csv      每样品一行 + 结论 + 理由（阈值写在里面）
    <sample>/frame_qc.csv, <sample>/qc.json
    raw_qc.png                            另跑 plot-frame-stability.py

    <RAW python> raw-qc-frames.py <root> --cfg <day>.cfg [--qmin 0.010]

判据全部写在 RULES 里，结论是可以被逐条反驳的，不是印象。
"""
from __future__ import annotations

import argparse, glob, json, os, sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qc_common as qc                                    # noqa: E402

RULES = dict(
    contrast_fatal=0.005,      # 低 q 对比度低于 0.5 % 直接不可用
    contrast_warn=0.02,        # 低于 2 % 报警（信号埋在背景里）
    ctl_drift_warn=0.03,       # control 在 run 内漂移 >3 %
    ctl_spread_warn=0.05,      # control 帧间起伏 >5 %
    sam_spread_warn=0.05,      # 样品帧间起伏 >5 %
    sam_drift_warn=0.10,       # 样品 run 内漂移 >10 %（辐照损伤/沉降）
    outlier_warn=0.20,         # 离群帧占比 >20 %
    jump_warn=0.15,            # 相邻帧跳变 >15 %（q=0.2）
    qmax_warn=0.15,            # 单帧 SNR>=2 的 q 上限低于 0.15
    err_ratio_warn=4.0,        # 单点帧间起伏是 RAW 报的 sigma 的 >4 倍
    upturn_warn=1.20,          # 扣减后 I(q1)/I(2 q1) > 1.2
    factor_warn=0.05,          # control 缩放因子偏离 1 超过 5 %
    run_ratio_warn=0.05,       # 两次 control run 之间低 q 电平差 >5 %
)


def frame_row(p, kind):
    q, I, E = p.getQ(), p.getI(), p.getErr()
    h = p.getParameter("counters") or {}
    lo = q < 0.02                       # low-q level bin (the 1/Rg region)
    i2 = int(np.argmin(abs(q - 0.2)))
    k1 = 1
    k2 = int(np.argmin(abs(q - 2 * q[k1])))
    try:
        tx = float(h.get("Transmitted_Beam", "nan"))
    except ValueError:
        tx = float("nan")
    return dict(file=os.path.basename(p.getParameter("filename") or "?"), kind=kind,
                run=str(h.get("Run_Number", "")), frame=str(h.get("Frame_Number", "")),
                tx=tx, i_lo=float(np.mean(I[lo])), e_lo=float(np.mean(E[lo])),
                i_02=float(I[i2]), e_02=float(E[i2]),
                lowq_upturn_raw=float(I[k1] / I[k2]),   # raw frames: ~q^-4 parasitic tail
                i_single_lowq=float(I[k1]), e_single_lowq=float(E[k1]))


def stats(v):
    v = np.asarray(v, float)
    if v.size == 0:
        return dict(n=0, med=float("nan"), rel_spread=float("nan"))
    med = float(np.median(v))
    return dict(n=int(v.size), med=med, min=float(v.min()), max=float(v.max()),
                rel_spread=float(np.std(v) / med) if med else float("nan"))


def drift_pct(v):
    """total drift across a run from a straight-line fit, as a fraction of the mean."""
    v = np.asarray(v, float)
    if v.size < 3:
        return float("nan")
    k = np.polyfit(np.arange(v.size), v, 1)[0]
    return float(k * (v.size - 1) / np.mean(v))


def assess(sample_dir, s, qmin):
    key = os.path.basename(sample_dir.rstrip("/"))
    files = sorted(glob.glob(os.path.join(sample_dir, "*.tif")))
    if not files:
        return None, "no .tif"
    sk = qc.guess_sample_key(files, key)
    if sk is None:
        return None, "no sample run (dir name does not match any run prefix)"
    sam_files, ctl_runs, _ = qc.split_runs(files, sk)
    ctl_files = [f for v in ctl_runs.values() for f in v]
    if not ctl_files:
        return None, "no control run next to the sample"

    sam = qc.average(sam_files, s)
    ctl = qc.average(ctl_files, s)
    rows = ([frame_row(p, "sample") for p in _load(sam_files, s)] +
            [frame_row(p, "control") for p in _load(ctl_files, s)])
    S = [r for r in rows if r["kind"] == "sample"]
    C = [r for r in rows if r["kind"] == "control"]

    by_run = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by_run[r["run"]][r["kind"]].append(r["i_lo"])
    run_spread, run_drift, ctl_spread, ctl_drift = [], [], [], []
    runmed = {"sample": [], "control": []}
    for run, kinds in sorted(by_run.items(), key=lambda kv: str(kv[0])):
        for kind, sp, dr in (("sample", run_spread, run_drift),
                             ("control", ctl_spread, ctl_drift)):
            v = kinds.get(kind) or []
            if len(v) >= 3:
                sp.append(stats(v)["rel_spread"])
                dr.append(drift_pct(v))
                runmed[kind].append(float(np.median(v)))

    def _run_ratio(v):
        return (max(v) / min(v) - 1) if len(v) > 1 and min(v) > 0 else float("nan")

    # frame-to-frame outliers at q=0.2 (where the particle scatters): the low-q bin is
    # dominated by the parasitic background shared with the control and hides sample-side
    # excursions (beam / injection events) the control run cannot cancel.
    n_out_mid, max_jump, worst = 0, 0.0, []
    for run in sorted({r["run"] for r in S}, key=str):
        v = [r for r in S if r["run"] == run]
        med = float(np.median([r["i_02"] for r in v]))
        prev = None
        for r in v:
            r["dev_midq"] = r["i_02"] / med - 1.0
            if abs(r["dev_midq"]) > 0.10:
                n_out_mid += 1
                worst.append("%s %+.0f%%" % (r["file"], 100 * r["dev_midq"]))
            if prev is not None:
                max_jump = max(max_jump, abs(r["i_02"] / prev - 1.0))
            prev = r["i_02"]
    s_med = stats([r["i_lo"] for r in S])["med"]
    n_out_low = sum(1 for r in S if abs(r["i_lo"] / s_med - 1) > 0.15) if s_med else 0

    def _err_ratio(Ikey, Ekey):        # single q point vs RAW's sigma for that point
        I = np.array([r[Ikey] for r in S], float)
        E = np.array([r[Ekey] for r in S], float)
        m = E > 0
        return float(np.std(I[m]) / np.mean(E[m])) if m.any() else float("nan")

    dif, factor = qc.subtract(sam, ctl)
    qd, Id, Ed = dif.getQ(), dif.getI(), dif.getErr()
    i0 = int(np.argmin(abs(qd - qmin))) if qmin > 0 else 0
    qd, Id, Ed = qd[i0:], Id[i0:], Ed[i0:]
    lo = qd < 0.02
    contrast_lo = float(np.mean(Id[lo]) / np.mean(ctl.getI()[i0:][lo]))
    mid = (qd > 0.25) & (qd < 0.35)
    contrast_mid = float(np.mean(Id[mid]) / np.mean(ctl.getI()[i0:][mid])) if mid.any() else float("nan")
    k1, k2 = 1, int(np.argmin(abs(qd - 2 * qd[1])))
    upturn = float(Id[k1] / Id[k2])
    snr = Id / np.where(Ed > 0, Ed, np.inf)
    snr_one = snr / np.sqrt(max(len(S), 1))
    qmax_one = float(qd[np.argmax(snr_one < 2)]) if np.any(snr_one < 2) else float(qd[-1])
    qmax_avg = float(qd[np.argmax(snr < 2)]) if np.any(snr < 2) else float(qd[-1])

    info = dict(sample=key, sample_run=sk, n_sample=len(S), n_control=len(C),
                control_runs={k: len(v) for k, v in ctl_runs.items()},
                i_lo_sample=stats([r["i_lo"] for r in S]),
                i_lo_control=stats([r["i_lo"] for r in C]),
                run_spread_sample=[round(x, 4) for x in run_spread],
                run_drift_sample=[round(x, 4) for x in run_drift],
                run_spread_control=[round(x, 4) for x in ctl_spread],
                run_drift_control=[round(x, 4) for x in ctl_drift],
                runmed_lo_sample=[round(x, 2) for x in runmed["sample"]],
                runmed_lo_control=[round(x, 2) for x in runmed["control"]],
                ctl_run_ratio=_r(_run_ratio(runmed["control"])),
                sam_run_ratio=_r(_run_ratio(runmed["sample"])),
                n_outliers_lowq=n_out_low, n_outliers_midq=n_out_mid,
                outlier_frac_midq=round(n_out_mid / max(len(S), 1), 3),
                max_jump_midq=round(max_jump, 3), worst_frames_midq=worst[:6],
                err_ratio_lowq=_r(_err_ratio("i_single_lowq", "e_single_lowq")),
                err_ratio_midq=_r(_err_ratio("i_02", "e_02")),
                tx=stats([r["tx"] for r in S if r["tx"] == r["tx"]]),
                control_scale=round(factor, 4), qmin_used=qmin,
                contrast_lowq=round(contrast_lo, 5), contrast_midq=_r(contrast_mid),
                lowq_upturn=round(upturn, 3), snr_lowq_average=round(float(np.mean(snr[lo])), 1),
                qmax_SNR2_single=round(qmax_one, 4), qmax_SNR2_average=round(qmax_avg, 4))

    R = RULES
    fatal, warn = [], []
    if contrast_lo < R["contrast_fatal"]:
        fatal.append("低 q 对比度 %.2f%% < %.1f%%：信号埋在背景里"
                     % (contrast_lo * 100, R["contrast_fatal"] * 100))
    elif contrast_lo < R["contrast_warn"]:
        warn.append("低 q 对比度只有 %.2f%%" % (contrast_lo * 100))
    for name, val, lim, msg in (
            ("background", max(ctl_drift + [0]), R["ctl_drift_warn"],
             "背景在 run 内漂移 %.1f%%（毛细管污染/束流变化）"),
            ("ctl spread", max(ctl_spread + [0]), R["ctl_spread_warn"], "背景帧间起伏 %.1f%%"),
            ("sample spread", max(run_spread + [0]), R["sam_spread_warn"], "样品帧间起伏 %.1f%%"),
            ("sample drift", max(run_drift + [0]), R["sam_drift_warn"],
             "样品在 run 内漂移 %.1f%%（辐照损伤/沉降）")):
        if val > lim:
            warn.append(msg % (val * 100))
    if n_out_low / max(len(S), 1) > R["outlier_warn"]:
        warn.append(">15%% 偏离中位数的帧占 %.0f%%（低 q 判据）"
                    % (100 * n_out_low / max(len(S), 1)))
    if n_out_mid / max(len(S), 1) > R["outlier_warn"]:
        warn.append("q=0.2 处偏离 run 中位数 >10%% 的帧占 %.0f%%（名单见 qc.json）"
                    % (100 * n_out_mid / max(len(S), 1)))
    if max_jump > R["jump_warn"]:
        warn.append("相邻帧最大跳变 %.0f%%（q=0.2）" % (100 * max_jump))
    if upturn > R["upturn_warn"]:
        warn.append("扣减后低 q 仍在上翘（I(q1)/I(2q1)=%.2f）：先跑 attribute-lowq-upturn.py，"
                    "不要直接当相互作用/聚集的证据" % upturn)
    if abs(factor - 1.0) > R["factor_warn"]:
        warn.append("control 缩放因子 %.3f 偏离 1 达 %.1f%%（背景通量/位置不匹配）"
                    % (factor, abs(factor - 1) * 100))
    for nm, val, msg in (("ctl run", info["ctl_run_ratio"],
                          "两次 control run 之间低 q 电平差 %.1f%%（单个缩放因子消不掉）"),
                         ("sam run", info["sam_run_ratio"],
                          "样品在两次 run 之间低 q 电平差 %.1f%%")):
        if val is not None and val > R["run_ratio_warn"]:
            warn.append(msg % (val * 100))
    if info["err_ratio_midq"] and info["err_ratio_midq"] > R["err_ratio_warn"]:
        warn.append("单帧在 q=0.2 的帧间起伏是 RAW 报的 sigma 的 %.1f 倍（通量归一化残差）"
                    % info["err_ratio_midq"])
    if qmax_one < R["qmax_warn"]:
        warn.append("单帧 SNR>=2 只到 q=%.3f（该浓度下 q 范围不够）" % qmax_one)
    info["warnings"], info["fatal"] = warn, fatal
    info["verdict"] = ("不可用" if fatal else "勉强可用" if len(warn) >= 2 else
                       "可用" if len(warn) == 1 else "好")
    return info, rows


def _r(x):
    x = None if x is None else x
    return None if x is None or x != x else round(float(x), 4)


def _load(files, s):
    profs, _ = qc.raw.load_and_integrate_images(sorted(files), settings=s)
    return profs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("root", help="directory holding one sub-directory per sample")
    ap.add_argument("--cfg", required=True, help="RAW .cfg for this dataset")
    ap.add_argument("--hdr-format", default=qc.HDR_FORMAT,
                    help="RAW image-header reader (per-frame .txt must sit next to the tif)")
    ap.add_argument("--qmin", type=float, default=0.0, metavar="Q",
                    help="[1/A] drop q below this before judging the subtracted curve; "
                         "0.010 skips the beam-stop-edge points (see the reference)")
    ap.add_argument("--out", default=None, help="default <root>/_rawqc")
    args = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(args.root))
    out = os.path.abspath(os.path.expanduser(args.out)) if args.out else os.path.join(root, "_rawqc")
    os.makedirs(out, exist_ok=True)
    s = qc.load_settings(args.cfg, args.hdr_format)

    infos = []
    for d in sorted(glob.glob(os.path.join(root, "*"))):
        if not os.path.isdir(d) or os.path.basename(d).startswith("_"):
            continue
        if not glob.glob(os.path.join(d, "*.tif")):
            continue
        print("== %s" % os.path.basename(d), flush=True)
        info, rows = assess(d, s, args.qmin)
        if info is None:
            print("   skipped: %s" % rows, flush=True)
            continue
        infos.append(info)
        sd = os.path.join(out, info["sample"])
        os.makedirs(sd, exist_ok=True)
        with open(os.path.join(sd, "qc.json"), "w") as fh:
            json.dump(dict(info=info, frames=rows), fh, indent=1, default=str)
        with open(os.path.join(sd, "frame_qc.csv"), "w") as fh:
            fh.write("file,kind,run,frame,tx,i_lo,e_lo,i_02,e_02,lowq_upturn_raw\n")
            for r in rows:
                fh.write("%s,%s,%s,%s,%.6f,%.6g,%.6g,%.6g,%.6g,%.3f\n"
                         % (r["file"], r["kind"], r["run"], r["frame"], r["tx"], r["i_lo"],
                            r["e_lo"], r["i_02"], r["e_02"], r["lowq_upturn_raw"]))
        print("   %s | contrast=%.2f%% scale=%.3f upturn=%.2f qmax1=%.3f snr_lo=%s | %s"
              % (info["verdict"], info["contrast_lowq"] * 100, info["control_scale"],
                 info["lowq_upturn"], info["qmax_SNR2_single"], info["snr_lowq_average"],
                 "; ".join(info["fatal"] + info["warnings"]) or "no issues"), flush=True)

    with open(os.path.join(out, "raw_quality.csv"), "w") as fh:
        cols = ["sample", "verdict", "n_sample", "n_control", "contrast_lowq", "contrast_midq",
                "control_scale", "lowq_upturn", "snr_lowq_average", "qmax_SNR2_single",
                "qmax_SNR2_average", "err_ratio_lowq", "err_ratio_midq", "outlier_frac_midq",
                "run_spread_sample", "run_drift_sample", "run_spread_control",
                "run_drift_control", "issues"]
        fh.write(",".join(cols) + "\n")
        for i in infos:
            fh.write(",".join(str(x) for x in [
                i["sample"], i["verdict"], i["n_sample"], i["n_control"], i["contrast_lowq"],
                i["contrast_midq"], i["control_scale"], i["lowq_upturn"],
                i["snr_lowq_average"], i["qmax_SNR2_single"], i["qmax_SNR2_average"],
                i["err_ratio_lowq"], i["err_ratio_midq"], i["outlier_frac_midq"],
                "|".join(str(x) for x in i["run_spread_sample"]),
                "|".join(str(x) for x in i["run_drift_sample"]),
                "|".join(str(x) for x in i["run_spread_control"]),
                "|".join(str(x) for x in i["run_drift_control"]),
                "; ".join(i["fatal"] + i["warnings"])]) + "\n")

    with open(os.path.join(out, "raw_quality.md"), "w") as fh:
        fh.write("# 原始数据质量评估（逐样品）\n\n判据：%s\n\n" % json.dumps(RULES, ensure_ascii=False))
        fh.write("| 样品 | 结论 | 对比度(低q) | 缩放因子 | 扣减后上翘 | SNR均值(lowq) | "
                 "单帧qmax | 平均qmax | 帧间/sigma(低q) | 帧间/sigma(0.2) | 问题 |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|---|---|\n")
        for i in infos:
            fh.write("| %s | %s | %.2f%% | %.3f | %.3f | %s | %.4f | %.4f | %s | %s | %s |\n"
                     % (i["sample"], i["verdict"], i["contrast_lowq"] * 100, i["control_scale"],
                        i["lowq_upturn"], i["snr_lowq_average"], i["qmax_SNR2_single"],
                        i["qmax_SNR2_average"], i["err_ratio_lowq"], i["err_ratio_midq"],
                        "; ".join(i["fatal"] + i["warnings"]) or "—"))
    print("wrote %s" % os.path.join(out, "raw_quality.md"))


if __name__ == "__main__":
    main()
