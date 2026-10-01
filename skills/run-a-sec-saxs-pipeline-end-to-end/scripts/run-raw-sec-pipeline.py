#!/usr/bin/env python3
"""端到端跑一条 SEC-SAXS 系列：**全部处理都在 RAW 里做**，本脚本只设参数、收产物。

设计原则（用户明确要求）：
  * 不自己写积分/拟合/P(r)/重建算法——只用 bioxtasraw.RAWAPI 的函数（它们就是 RAW GUI 面板背后的实现）；
  * 逐帧归一化在 RAW 内部完成（settings: ImageHdrFormat='BL19U2, SSRF' + EnableNormalization +
    NormalizationList=[['/','Transmitted_Beam']]，每帧除以该帧 header txt 里的 Transmitted_Beam），
    因此**不需要写归一化后的 tif**；
  * 每个节点都落 `.dat`：逐帧积分 → buffer 平均 → 逐帧扣减 → 样品区平均 → （可选）基线校正 →
    每个 Guinier 区间的 profile 副本 → IFT(.ift/.out)；
  * 拟合好坏要能被复核：逐帧 Rg/I0/MW、多区间 Guinier（各自 Rg/I0/误差/qRg 边界/r²）、IFT 的
    χ²/α/evidence、平台一致性，全部写成表 + 图 + RAW 自带 PDF 报告。

典型用法（先跑 emit-bl19u2-header-txt.py 生成 work 目录）：
  run-raw-sec-pipeline.py --work-dir PROC/bsa/work --out-dir PROC/bsa --cfg data/20261001.cfg

**多峰（≥2 个洗脱峰）**：默认 --multi-peak auto —— 认到 ≥2 个峰就**逐峰**建子目录，各自用
"本峰紧邻的前后两段 buffer + 本峰窗口"重跑扣减与下游；认峰结果总是画成图
（series/sec_peaks.png、series/sec_peaks_zoom.png）给人复核。RAW 自己的 find_sample_range
只认最大的那个峰（SASCalc.findSampleRange 里 argmax），所以 ≥2 峰的数据用原流程会静默丢掉别的峰。
先看有几个峰：
  run-raw-sec-pipeline.py ... --steps integrate,peaks        # 只出峰检测产物（秒级判读）
改参数只看峰（不重跑积分）：find-sec-peaks.py --products <产物目录> [--sweep]
人工定下来的区间再跑全流程：--peak-ranges "..." [--peak-buffers "..."]

注意：必须在**不是 RAW 源码目录**的地方运行（源码树里的 bioxtasraw/ 会遮蔽 site-packages 里
编译好的那个，缺 sascalc_exts 扩展会直接 ImportError）。
"""
import argparse
import copy
import csv
import glob
import json
import os
import sys
import time
import traceback

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import bioxtasraw.RAWAPI as raw
except ModuleNotFoundError as exc:  # 最常见：从 RAW 源码目录里跑
    raise SystemExit(f"无法 import bioxtasraw.RAWAPI（{exc}）。\n"
                     "· 用装了 RAW 的 python 跑（如 /Applications/BioXTASRAW/bin/python）；\n"
                     "· 并且不要站在 RAW 源码目录里跑（源码树会遮蔽 site-packages 里的编译扩展）。")

try:
    from bioxtasraw.SASExceptions import NoATSASError
except Exception:  # pragma: no cover
    class NoATSASError(Exception):
        pass

# ------------------------------------------------------------------- 结果 README
# README.md 由**另一个技能**统一生成（`write-saxs-results-readme`），本管线不自己写模板：
# 两条流水线（SEC / 管式）套同一份实现 → 章节、表格列、判据、采用口径逐条对齐。
README_SKILL = "write-saxs-results-readme"


def find_readme_writer(explicit=None):
    """找共享生成器 `write-readme.py`：① `--readme-script` 显式给；
    ② 同类目下的同名技能目录（安装后 <profile>/skills/<category>/ 下互为兄弟）；
    ③ 各 profile 里搜一遍。找不到返回 None —— README 是交付物不是计算步骤，
    不能因为它没装就让整条管线失败（只在日志里明确告警）。"""
    if explicit:
        return explicit if os.path.isfile(explicit) else None
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(os.path.dirname(here), README_SKILL, "scripts", "write-readme.py"),
                 os.path.join(os.path.dirname(os.path.dirname(here)), README_SKILL,
                              "scripts", "write-readme.py")):
        if os.path.isfile(cand):
            return cand
    hits = sorted(glob.glob(os.path.expanduser(
        "~/.hermes/**/%s/scripts/write-readme.py" % README_SKILL), recursive=True))
    return hits[0] if hits else None


def write_results_readme(out_dir, mode, script=None):
    """调共享实现写 `<out_dir>/README.md`（同解释器内 import；那份实现只用标准库）。

    独立补写/验收（不重算）：`python <writer> <目录>`、`python verify-results-folder.py <目录>`。
    """
    writer = find_readme_writer(script)
    if writer is None:
        raise RuntimeError(
            "找不到 %s：hermes skills install flmaximwang/AgentSkill-DoingSAXS/skills/%s "
            "--category saxs -y" % (README_SKILL, README_SKILL))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "saxs_readme_common", os.path.join(os.path.dirname(writer), "readme_common.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    path, keys, _f = mod.write_readme(out_dir, mode)
    return path, len(keys or [])


# 多峰识别与逐峰区间划分（RAW 自己只认最大的那个峰，见 sec_peaks 模块头）
try:
    import sec_peaks as sp
except Exception:                                     # pragma: no cover
    sp = None

FMT = argparse.ArgumentDefaultsHelpFormatter
STEPS = ["integrate", "peaks", "series", "guinier", "ift", "mw", "shape", "report"]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ----------------------------------------------------------------- RAW 设置
def make_settings(cfg, atsas_dir=None, header_normalization=True):
    st = raw.load_settings(cfg)
    if header_normalization:
        st.set('ImageHdrFormat', 'BL19U2, SSRF')  # 认 <帧名>.txt 作为 header
        st.set('EnableNormalization', True)       # 打开图像归一化
        st.set('NormalizationList', [['/', 'Transmitted_Beam']])  # 除以 header 里的 Transmitted_Beam
    else:
        # 该系列没有逐帧 txt（线站没给监视器/日志）→ 关掉，避免 RAW 去找不存在的 <帧名>.txt
        st.set('ImageHdrFormat', 'None')
        st.set('EnableNormalization', False)
    if atsas_dir:
        st.set('ATSASDir', atsas_dir)
    return st


def qindex(profile, q):
    """q → 该 profile 当前 q 向量里的整数下标。"""
    return int(np.argmin(np.abs(profile.getQ() - q)))


def save_dat(profile, name, folder):
    os.makedirs(folder, exist_ok=True)
    return raw.save_profile(profile, name, folder)


# ----------------------------------------------------------------- 各步骤
def step_integrate(files, st, out, prefix, limit):
    log(f"积分 {len(files)} 帧（RAW 径向积分，含逐帧 header 归一化）")
    profiles, imgs = raw.load_and_integrate_images(files, st)
    fr_dir = os.path.join(out, "profiles", "01_integrated")
    os.makedirs(fr_dir, exist_ok=True)
    rows = []
    for p in profiles:
        stem = os.path.splitext(p.getParameter('filename'))[0]
        save_dat(p, stem + ".dat", fr_dir)
        c = p.getParameter('counters') or {}
        tb = c.get('Transmitted_Beam')
        rows.append((stem, tb if tb is not None else "",
                     float(p.getI()[0]), float(p.getQ()[0]), float(p.getQ()[-1]),
                     float(np.trapezoid(p.getI(), p.getQ()))))
    with open(os.path.join(out, "tables", "frames_integrated.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["stem", "Transmitted_Beam_as_read", "I_q0", "q_min", "q_max", "int_I_dq"])
        w.writerows(rows)
    log(f"  → {len(profiles)} 条曲线 → {fr_dir}/*.dat + tables/frames_integrated.csv")
    return profiles


def clip_ranges(rngs, n_frames, what):
    """把手工给的区间收进 [0, n_frames-1]：RAW 的区间是 0 基闭区间，end 写大了会 IndexError
    （实测：445 帧的系列给 '400,445' → SECM.averageFrames 里 list index out of range）。"""
    for r in rngs:
        if r[1] > n_frames - 1:
            log(f"  ！{what} 区 {r} 超出帧范围（共 {n_frames} 帧，最大下标 {n_frames-1}）→ 收到 {n_frames-1}")
            r[1] = n_frames - 1
        if r[0] < 0:
            log(f"  ！{what} 区 {r} 起点 < 0 → 收到 0")
            r[0] = 0
    return rngs


def step_ranges_and_subtraction(profiles, st, out, prefix, args, buffers=None, sample=None,
                                baseline_ranges=None, series=None, det=None, peak=None):
    """buffer 区 → 扣减 →（可选）基线 → sample 区。返回样品平均曲线与 series。

    ``buffers`` / ``sample`` 给了就用给定的（0 基闭区间；buffers 可以是多段 [[s,e],...]）——
    多峰模式就是这么逐峰复用本函数的；都不给 = 旧的单峰行为（args 里手工区间优先，否则 RAW 自动）。
    """
    series = series if series is not None else raw.profiles_to_series(profiles, st)

    # ---- buffer 区
    if buffers:
        b_rng = [[int(a), int(b)] for a, b in buffers]
        clip_ranges(b_rng, len(profiles), "buffer")
        b_ok = True
    elif args.buffer_range:
        b_rng = [[int(x) for x in part.split(",")] for part in args.buffer_range.split(";")]
        clip_ranges(b_rng, len(profiles), "buffer")
        b_ok = True
    else:
        ok, s, e = raw.find_buffer_range(series)
        if s is None or e is None:
            raise SystemExit(f"RAW 没找到 buffer 区（success={ok}）：系列里可能没有可识别的洗脱峰"
                             "（例如只取了峰前的帧）。用 --buffer-range start,end 手工指定（0 基帧号），"
                             "或让它按 --peaks/--peak-ranges 逐峰自己挑 buffer。")
        b_ok, b_rng = bool(ok), [[int(s), int(e)]]
    log(f"buffer 区 {b_rng}（{'给定' if (buffers or args.buffer_range) else 'RAW 自动找'}, success={b_ok}）")
    (sub_profiles, rg, rger, i0, i0er, vcmw, vcmwer, vpmw) = raw.set_buffer_range(series, b_rng)

    buf_dir = os.path.join(out, "profiles", "02_buffer")
    buf_sasms = [s for r in b_rng for s in series.getSASMList(r[0], r[1], 'unsub')]
    buf_avg = raw.average(buf_sasms)
    save_dat(buf_avg, "buffer_avg.dat", buf_dir)

    sub_dir = os.path.join(out, "profiles", "03_subtracted")
    for p in sub_profiles:
        stem = os.path.splitext(p.getParameter('filename'))[0]
        save_dat(p, stem + "_sub.dat", sub_dir)
    log(f"  → buffer 平均 + {len(sub_profiles)} 条扣减曲线 已存 .dat")

    # ---- 基线校正（可选）
    pt = 'sub'
    if args.baseline != 'none':
        try:
            if args.baseline == 'integral':
                start_found, end_found, start_range, end_range = raw.find_baseline_range(series)
                bl_type = 'Integral'
            else:
                bl_spec = baseline_ranges or args.baseline_ranges or "0,20;1800,1990"
                start_range, end_range = bl_spec.split(";")
                start_range = [int(x) for x in start_range.split(",")]
                end_range = [int(x) for x in end_range.split(",")]
                bl_type = 'Linear'
            res = raw.set_baseline_correction(series, start_range, end_range, bl_type)
            bl_profiles, bl_corr, bl_fit = res[0], res[-2], res[-1]
            pt = 'baseline'
            bl_dir = os.path.join(out, "profiles", "05_baseline")
            for p in bl_profiles:
                stem = os.path.splitext(p.getParameter('filename'))[0]
                save_dat(p, stem + "_bl.dat", bl_dir)
            save_dat(raw.average(series.getSASMList(0, len(profiles) - 1, 'baseline')),
                     "baseline_allframes_avg.dat", bl_dir)
            log(f"  → {bl_type} 基线校正完成（区间 {start_range}–{end_range}）→ 后续 profile_type='baseline'")
        except Exception as exc:
            log(f"  ！基线校正失败，退回 'sub'：{type(exc).__name__}: {exc}")

    # ---- sample 区 + 样品平均
    if sample:
        s_rng = [[int(sample[0]), int(sample[1])]]
        clip_ranges(s_rng, len(profiles), "sample")
        s_ok = True
        s_valid = None
    elif args.sample_range:
        s0, s1 = (int(x) for x in args.sample_range.split(","))
        s_ok, s_rng = True, [[s0, s1]]
        clip_ranges(s_rng, len(profiles), "sample")
        s_valid = None
    else:
        ok, s, e = raw.find_sample_range(series, profile_type=pt)
        if s is None or e is None:
            raise SystemExit(f"RAW 没找到样品区（success={ok}）。用 --sample-range start,end 手工指定（0 基帧号）。")
        s_ok, s_rng = bool(ok), [[int(s), int(e)]]
        s_valid = None
    log(f"样品区 {s_rng[0]}（{'给定' if (sample or args.sample_range) else 'RAW 自动找'}, success={s_ok}, profile_type={pt}）")
    sample_profile = raw.set_sample_range(series, s_rng, profile_type=pt)
    sample_dir = os.path.join(out, "profiles", "04_sample")
    save_dat(sample_profile, "sample_avg.dat", sample_dir)

    # ---- 可选：给下游（IFT/MW）用的低 q 裁剪。低 q 被寄生散射/聚集污染时，
    # BIFT/DENSS 会给出离谱的 Dmax（本机实测未裁时 Dmax 417 Å / Rg 149 Å）。裁剪用 RAW 自己的 setQrange，
    # 不是自写拟合：只是把"哪些点参与"变成显式参数。
    if args.trim_qmin or args.trim_qmax:
        old_range = sample_profile.getQrange()
        i0i = qindex(sample_profile, args.trim_qmin) if args.trim_qmin else 0      # 相对当前 getQ() 的下标
        if args.trim_qmax:                                                          # 上界：模型表达不了高 q 细结构时用
            i1i = qindex(sample_profile, args.trim_qmax) or len(sample_profile.getQ())
        else:
            i1i = old_range[1] - old_range[0]
        sample_profile = copy.deepcopy(sample_profile)
        # setQrange 用绝对下标 q[start:end]
        sample_profile.setQrange((old_range[0] + i0i, old_range[0] + i1i))
        save_dat(sample_profile, "sample_avg_qmin%s_qmax%s.dat" % (
            "%.4f" % args.trim_qmin if args.trim_qmin else "unset",
            "%.4f" % args.trim_qmax if args.trim_qmax else "unset"), sample_dir)
        log(f"下游（Guinier 表/IFT/MW/重建）改用裁剪后曲线：q {sample_profile.getQ()[0]:.4f}–"
            f"{sample_profile.getQ()[-1]:.4f} 1/A（{len(sample_profile.getQ())} 点）")

    # ---- 逐帧参数（注意：SECM 的 getRg/getI0/getVcMW/getVpMW 返回 (值, 误差) 两个数组）
    frames = np.asarray(series.getFrames())
    rg_a, rger_a = (np.asarray(x, dtype=float) for x in series.getRg())
    i0_a, i0er_a = (np.asarray(x, dtype=float) for x in series.getI0())
    vc_a, vcer_a = (np.asarray(x, dtype=float) for x in series.getVcMW())
    vp_a, vper_a = (np.asarray(x, dtype=float) for x in series.getVpMW())
    int_a = np.asarray(series.getIntI('unsub'), dtype=float)
    with open(os.path.join(out, "tables", "frame_params.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "rg", "rg_err", "i0", "i0_err", "vc_mw", "vc_mw_err", "vp_mw", "int_unsub"])
        for i in range(len(frames)):
            w.writerow([float(frames[i]), float(rg_a[i]), float(rger_a[i]), float(i0_a[i]), float(i0er_a[i]),
                        float(vc_a[i]), float(vcer_a[i]), float(vp_a[i]), float(int_a[i])])
    try:
        plot_series(out, frames, int_a, rg_a, i0_a, vc_a, vp_a, b_rng[0], s_rng[0], det=det)
    except Exception as exc:
        log(f"  ！series 图失败：{type(exc).__name__}: {exc}")

    ranges = dict(buffer=b_rng, sample=s_rng, profile_type=pt,
                  sample_range_success=s_ok, buffer_range_success=b_ok)   # type: dict
    if det is not None:
        ranges['detected_peaks'] = [dict(index=p['index'], apex=p['apex'], window=list(p['window']),
                                         buffers=[list(x) for x in p.get('buffer_ranges', [])])
                                    for p in det.get('peaks', [])]
    if peak is not None:
        ranges['this_is_peak'] = peak
    with open(os.path.join(out, "series", "ranges.json"), "w") as fh:
        json.dump(ranges, fh, indent=2)
    return series, sample_profile, ranges


def plot_series(out, frames, int_a, rg_a, i0_a, vc_a, vp_a, b_rng, s_rng, det=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    axes[0].plot(frames, int_a, lw=0.8)
    axes[0].set_ylabel("integrated I (unsubtracted)")
    axes[1].plot(frames, rg_a, lw=0.8)
    axes[1].set_ylabel("Rg (A)")
    axes[2].plot(frames, i0_a, lw=0.8)
    axes[2].set_ylabel("I(0)")
    axes[3].plot(frames, vc_a, lw=0.8, label="Vc MW")
    axes[3].plot(frames, vp_a, lw=0.8, label="Vp MW")
    axes[3].set_ylabel("MW (kDa)")
    axes[3].set_xlabel("frame")
    axes[3].legend(fontsize=8)
    for ax in axes:
        for rng, c, lab in ((b_rng, 'tab:blue', 'buffer'), (s_rng, 'tab:green', 'sample')):
            ax.axvspan(rng[0], rng[1], color=c, alpha=0.15)
        if det is not None:                    # 多峰时：把**所有**检测到的峰窗都标出来
            for p in det.get('peaks', []):
                cc = sp.PEAK_COLORS[p['index'] % len(sp.PEAK_COLORS)] if sp else 'tab:gray'
                ax.axvspan(p['window'][0], p['window'][1], color=cc, alpha=0.10)
                ax.axvline(p['apex'], color=cc, lw=0.9, ls='--', alpha=0.9)
                ax.annotate(f"P{p['index']+1}", (p['apex'], ax.get_ylim()[1]),
                            fontsize=8, color=cc, ha='left', va='top')
        ax.grid(alpha=0.3)
    fig.suptitle("SEC-SAXS series: chromatogram, Rg, I(0), MW (shaded = buffer / sample ranges)")
    fig.tight_layout()
    os.makedirs(os.path.join(out, "series"), exist_ok=True)
    fig.savefig(os.path.join(out, "series", "series_plot.png"), dpi=130)
    plt.close(fig)


def step_guinier(sample_profile, st, out, args):
    """多区间 Guinier：每个区间一份 profile 副本（qrange 截到该区间）+ 一张汇总表/图。

    判据（交给用户复核）：qRg_min ≳0.3（避开低 q 寄生散射）、qRg_max 按形状 0.65–1.3、r² 与 Rg 误差大小、
    以及不同区间 Rg 是否稳定。
    """
    auto = raw.auto_guinier(sample_profile, settings=st)
    rg_a = float(auto[0])
    log(f"auto_guinier: Rg={auto[0]:.2f} A (err {auto[2]:.2f}), q={auto[4]:.4f}–{auto[5]:.4f}, "
        f"qRg={auto[6]:.2f}–{auto[7]:.2f}, r²={auto[10]:.5f}")

    q = sample_profile.getQ()
    if args.guinier_ranges:
        rngs = []
        for part in args.guinier_ranges.split(","):
            lo, hi = part.split(":")
            rngs.append((f"q{float(lo):.4f}-{float(hi):.4f}", float(lo), float(hi)))
    else:  # 以 auto Rg 为锚，铺一条跨判据边界的阶梯，方便你判断"哪段才合适"
        rngs = [(f"qRg{lo}-{hi}", lo / rg_a, hi / rg_a)
                for lo, hi in ((0.3, 0.6), (0.4, 0.8), (0.5, 1.0), (0.6, 1.3))]
        rngs = [(lab, max(q[0], lo), min(q[-1], hi)) for lab, lo, hi in rngs]

    rows = [("auto", rg_a, float(auto[1]), float(auto[2]), float(auto[3]), float(auto[4]),
             float(auto[5]), float(auto[6]), float(auto[7]), float(auto[10]))]
    fit_dir = os.path.join(out, "profiles", "06_guinier")
    report_profiles = []
    for lab, qlo, qhi in rngs:
        i0i, i1i = qindex(sample_profile, qlo), qindex(sample_profile, qhi)
        p2 = copy.deepcopy(sample_profile)
        rg, i00, rge, i0e, qmin, qmax, qrgmin, qrgmax, r2 = raw.guinier_fit(p2, i0i, i1i, settings=st)
        p2.setQrange((i0i, i1i + 1))          # 让 RAW 报告只画这个区间
        save_dat(p2, f"guinier_{lab}.dat", fit_dir)
        report_profiles.append(p2)
        rows.append((lab, float(rg), float(i00), float(rge), float(i0e), float(qmin), float(qmax),
                     float(qrgmin), float(qrgmax), float(r2)))
        log(f"  {lab:>16}: Rg={rg:7.2f}±{rge:5.2f} A  I0={i00:9.3g}  qRg={qrgmin:.2f}–{qrgmax:.2f}  r²={r2:.5f}")

    tbl = os.path.join(out, "tables", "guinier_multi_range.csv")
    os.makedirs(os.path.dirname(tbl), exist_ok=True)
    with open(tbl, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["range_label", "rg", "i0", "rg_err", "i0_err", "q_min", "q_max",
                    "qRg_min", "qRg_max", "r_sqr"])
        w.writerows(rows)
    log(f"  → {tbl}")

    # 多区间同框图（用 RAW 返回的 Rg/I0 画模型线，不自己拟合）
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.5, 5))
        qq, ii, ee = sample_profile.getQ(), sample_profile.getI(), sample_profile.getErr()
        ax.errorbar(qq, ii, yerr=ee, fmt='.', ms=3, lw=0.7, color='0.5', label='sample avg')
        qf = np.linspace(qq[0], qq[-1], 400)
        for (lab, rg, i00, *_rest) in [(r[0], r[1], r[2]) for r in rows]:
            ax.plot(qf, i00 * np.exp(-(qf ** 2) * rg ** 2 / 3.0), lw=1.2, label=f"{lab}: Rg={rg:.1f}")
        ax.set_yscale('log')
        ax.set_xlabel("q (1/A)")
        ax.set_ylabel("I(q)")
        ax.set_title("Guinier model lines from RAW's fits, overlaid on the sample profile")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "tables", "guinier_multi_range.png"), dpi=130)
        plt.close(fig)
    except Exception as exc:
        log(f"  ！Guinier 汇总图失败：{type(exc).__name__}: {exc}")
    return rows, report_profiles


def write_ift_summary(out, rows):
    with open(os.path.join(out, "tables", "ift_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["method", "dmax", "dmax_err", "rg", "rg_err", "chi_sq", "log_alpha", "evidence"])
        w.writerows(rows)


def step_ift(sample_profile, st, out, prefix, atsas_dir, ift_dmax=None):
    """IFT：BIFT（RAW 原生）+ （可选）显式 Dmax 的 DIFT（给 DENSS 用）+ GNOM（需 ATSAS）。"""
    ifts, rows = [], []
    try:
        b = raw.bift(sample_profile, settings=st, single_proc=True)
        ift = b[0]
        raw.save_ift(ift, f"{prefix}_bift.ift", os.path.join(out, "ifts"))
        ifts.append(ift)
        # bift 返回顺序：0 ift / 1 dmax / 2 rg / 3 i0 / 4 dmax_err / 5 rg_err / 6 i0_err /
        #                7 chi_sq / 8 log_alpha / 9 ld_err / 10 evidence / 11 evidence_err
        rows.append(("BIFT", float(b[1]), float(b[4]), float(b[2]), float(b[5]),
                     float(b[7]), float(b[8]), float(b[10])))
        log(f"BIFT : Dmax={b[1]:.1f} A  Rg={b[2]:.2f} A  chi²={b[7]:.1f}  log_alpha={b[8]:.2f}  evidence={b[10]:.1f}")
    except Exception as exc:
        log(f"  ！BIFT 失败：{type(exc).__name__}: {exc}")
        ift = None

    # 显式 Dmax 的 DIFT：DENSS 的推荐输入。DENSS 按 dmax 建盒子；Dmax 被低 q 拖大 → 盒子里密度太稀，
    # 收缩包络塌成空 support → DENSS 报 `IndexError: index -1 is out of bounds ... labeled_support == feature`
    # （DENSS.py:1931 在 num_features==0 时 sums 长度为 0）。修法是回去修 IFT 输入，不是调 DENSS 参数。
    denss_ift_obj = None
    if ift_dmax:
        try:
            d = raw.denss_ift(sample_profile, dmax=float(ift_dmax))
            denss_ift_obj = d[0]
            raw.save_ift(denss_ift_obj, f"{prefix}_denss_ift.ift", os.path.join(out, "ifts"))
            ifts.append(denss_ift_obj)
            rows.append(("DIFT", float(d[1]), float('nan'), float(d[2]), float(d[4]),
                         float(d[6]), float('nan'), float('nan')))
            log(f"DIFT : Dmax={d[1]:.1f} A（显式）Rg={d[2]:.2f} A chi²={d[6]:.2f} alpha={d[7]:.2f} → 给 DENSS 用")
        except Exception as exc:
            log(f"  ！DIFT 失败：{type(exc).__name__}: {exc}")

    dmax = None
    if atsas_dir:
        try:
            # 截断后的曲线（尤其裁掉高 q）常让 auto_dmax 返回 -1 → GNOM 直接报 'rmax got -1'。
            # 有显式 --ift-dmax 就用它，否则才让 RAW 自己定。
            dmax = float(ift_dmax) if ift_dmax else raw.auto_dmax(sample_profile)
            if not dmax or dmax <= 0:
                raise RuntimeError(f"auto_dmax 返回 {dmax}（曲线被截断后常见）→ 请用 --ift-dmax 显式指定 Dmax")
            g = raw.gnom(sample_profile, dmax)
            gnom_ift = g[0]
            raw.save_ift(gnom_ift, f"{prefix}_gnom.out", os.path.join(out, "ifts"))
            # gnom 返回顺序：0 ift / 1 dmax / 2 rg / 3 i0 / 4 rg_err / 5 i0_err / 6 chi_sq / ...
            rows.append(("GNOM", float(dmax), float('nan'), float(g[2]), float(g[4]),
                         float(g[6]), float('nan'), float('nan')))
            log(f"GNOM : Dmax={dmax:.1f} A  Rg={g[2]:.2f}±{g[4]:.2f} A  chi²={g[6]:.2f}  → {prefix}_gnom.out")
            ifts = ifts + [gnom_ift]
            ift = gnom_ift
            write_ift_summary(out, rows)      # 别在 early return 前漏掉这张表
            return (ift, ifts, rows, dmax, denss_ift_obj)
        except NoATSASError as exc:
            log(f"  ！GNOM 需要 ATSAS（未装/未指定 --atsas-dir）：{exc}")
        except Exception as exc:
            log(f"  ！GNOM 失败：{type(exc).__name__}: {exc}")
    write_ift_summary(out, rows)
    return (ift, ifts, rows, dmax, denss_ift_obj)


def step_mw(sample_profile, st, out, ift, atsas_dir):
    rows = []
    for name, fn, keys in (("Vc", raw.mw_vc, ("mw_vc", "vcor", "mw_err", "qmax")),
                           ("Vp", raw.mw_vp, ("mw_vp", "pvol_cor", "pvol", "qmax"))):
        try:
            res = fn(sample_profile, settings=st)
            rows.append((name, *[float(x) if isinstance(x, (int, float, np.floating)) else str(x) for x in res]))
            log(f"MW {name}: {res[0]:.1f} kDa (details {res[1:]})")
        except Exception as exc:
            log(f"  ！MW {name} 失败：{type(exc).__name__}: {exc}")
    # 注意：mw_bayes / mw_datclass 没有 settings 参数（签名：profile, rg, i0, first, atsas_dir, ...）
    for name, fn in (("Bayesian", raw.mw_bayes), ("Datclass", raw.mw_datclass)):
        try:
            res = fn(sample_profile, atsas_dir=atsas_dir)
            rows.append((name, *[float(x) if isinstance(x, (int, float, np.floating)) else str(x)
                                 for x in res[:5]]))   # mw + 最多 4 个 detail（Bayes 的 CI 就在里面）
            log(f"MW {name}: {res[0]:.1f} kDa")
        except Exception as exc:
            log(f"  ！MW {name} 跳过（需 ATSAS）：{type(exc).__name__}: {exc}")
    with open(os.path.join(out, "tables", "mw.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["method", "mw", "detail1", "detail2", "detail3", "detail4"])
        w.writerows(rows)
    return rows


def step_shape(ift, ifts, out, prefix, args, atsas_dir, model_ift=None):
    """形状重建：**电子云（RAW 原生 DENSS）总跑**；**珠模（DAMMIF）仅在有 ATSAS 时跑**。

    两者是同一份 IFT 的两种重建：DENSS 给电子密度图（.mrc），DAMMIF 给 dummy-atom 珠模（.pdb），
    互不替代——所以默认两个都要，缺 ATSAS 时退化成"只出电子云"。
    """
    mdir = os.path.join(out, "models")
    os.makedirs(mdir, exist_ok=True)
    engine = args.model_engine
    if engine == 'none':
        log("形状重建：按参数跳过")
        return
    if engine == 'auto':
        do_denss, do_dammif = True, bool(atsas_dir)
    else:
        do_denss = engine in ('denss', 'both')
        do_dammif = engine in ('dammif', 'both')
    rows = {"denss": [], "dammif": [], "damaver": []}
    if do_denss:
        denss_in = model_ift if model_ift is not None else ift
        if denss_in is None:
            log("  ！DENSS 需要一个 IFT，当前没有 → 跳过")
        try:
            res = raw.denss(denss_in, f"{prefix}_denss", mdir, mode=args.denss_mode)
            log(f"  DENSS（电子云）: chi²={res[1]:.2f} Rg={res[2]:.1f} A "
                f"support_vol={res[3]:.0f} side={res[4]:.1f} A  mode={args.denss_mode}")
            rows["denss"].append([str(res[0]), float(res[1]), float(res[2]), float(res[3]), float(res[4])])
        except Exception as exc:
            log(f"  ！DENSS 失败：{type(exc).__name__}: {exc}")
    if do_dammif:
        if not atsas_dir:
            log("  ！珠模 DAMMIF 需要 ATSAS（未指定 --atsas-dir）→ 本次只出电子云（DENSS）")
            return
        if ift is None:
            log("  ！DAMMIF 需要 GNOM 的 IFTM，当前没有 → 跳过珠模")
            return
        files = []
        for i in range(args.n_models):
            try:
                res = raw.dammif(ift, f"{prefix}_dammif_{i+1:02d}", mdir, mode=args.dammif_mode,
                                 symmetry=args.symmetry, atsas_dir=atsas_dir,
                                 model_format=args.model_format)
                log(f"  DAMMIF #{i+1}: chi²={res[0]:.2f} Rg={res[1]:.1f} Dmax={res[2]:.1f} MW={res[3]:.0f}")
                rows["dammif"].append([float(res[0]), float(res[1]), float(res[2]), float(res[3])])
                # dammif 写出的模型名是 <prefix>-1.<model_format>；damaver 只要文件名、不要路径
                files.append(f"{prefix}_dammif_{i+1:02d}-1.{args.model_format}")
            except Exception as exc:
                log(f"  ！DAMMIF #{i+1} 失败：{type(exc).__name__}: {exc}")
                break
        if len(files) >= 2:
            try:
                a = raw.damaver(files, f"{prefix}_damaver", mdir,
                                model_format=args.model_format, atsas_dir=atsas_dir)
                # 返回是 (mean_nsd, stdev_nsd, include_list, ...)：第 2 个是标准差，别当一致性指标
                log(f"  DAMAVER: 平均 NSD={a[0]:.3f} ± {a[1]:.3f}"
                    f"（代表模型见 models/{prefix}_damaver-global-summary.txt）")
                rows["damaver"] = [float(a[0]), float(a[1])]
            except Exception as exc:
                log(f"  ！DAMAVER 失败：{type(exc).__name__}: {exc}")
    return rows


def step_report(profiles, sample_profile, report_profiles, ifts, series, out, prefix):
    try:
        os.makedirs(os.path.join(out, "reports"), exist_ok=True)
        proflist = [sample_profile] + report_profiles
        raw.save_report(f"{prefix}_raw_report.pdf", os.path.join(out, "reports"), proflist, ifts, [series])
        log(f"  → RAW 报告 reports/{prefix}_raw_report.pdf（Guinier/IFT/系列图都在里面）")
    except Exception as exc:
        log(f"  ！RAW 报告失败：{type(exc).__name__}: {exc}\n{traceback.format_exc()[-500:]}")


# ------------------------------------------------------------- 多峰识别与逐峰分析
def peak_kwargs(args):
    """args → sec_peaks.detect_peaks 的参数（默认值的唯一来源是 sec_peaks）。"""
    return dict(q_range=tuple(float(x) for x in args.peak_q_range.split(",")),
                smooth_window=args.peak_smooth_window,
                baseline_window=args.peak_baseline_window,
                min_height=args.peak_min_height, min_prominence=args.peak_min_prominence,
                min_snr=args.peak_min_snr, min_width=args.peak_min_width,
                max_peaks=args.peak_max_n, window_mode=args.peak_window,
                buffer_mode=args.peak_buffer, buffer_min_len=args.peak_buffer_min,
                buffer_max_len=args.peak_buffer_max)


def chromatogram_and_mask(profiles, series, st, args):
    """认峰用的色谱图 + "无束流帧"掩码。返回 (y, mask, 说明, sub_profiles|None)。

    **两条都是低 q 窗口 [qlo,qhi] 的积分强度**（``SASM.getIofQRange``），不是总强度：
    总强度被基线漂移主导（实测 BSA 2000 帧只有 ~10% 起伏），任何阈值都会找出十几个假峰。
    * ``sub``（默认）：先用（命令行给的或 RAW 自动找的）全局 buffer 扣减一次，再取扣减后曲线
      —— 组分分得开，代价是依赖这次 buffer 选得对（选错就看图能看出来）；
    * ``unsub``：直接用未扣减曲线 —— 不依赖 buffer，代价是基线漂移更大。
    """
    mask = sp.no_beam_mask(np.array([p.getTotalI() for p in profiles]))
    q_range = tuple(float(x) for x in args.peak_q_range.split(","))
    bufs = sp.parse_ranges(args.buffer_range) if args.buffer_range else None
    if bufs is None:
        ok, s, e = raw.find_buffer_range(series)
        if s is not None and e is not None:
            bufs = [[int(s), int(e)]]
            log(f"认峰用的全局 buffer：RAW 自动找 {bufs[0]}（success={ok}）")
        else:
            log(f"  ！RAW 没找到 buffer 区（success={ok}）→ 认峰改用**未扣减**曲线（不依赖 buffer）")
    if bufs is None or args.peak_chromatogram == "unsub":
        why = "unsub（明确指定）" if args.peak_chromatogram == "unsub" else "unsub（没有可用的全局 buffer）"
        return sp.edf_frames(profiles, *q_range), mask, why, None
    try:
        subs, *_ = raw.set_buffer_range(series, bufs, do_calcs=False)
    except Exception as exc:
        log(f"  ！全局 buffer {bufs} 扣减失败（{type(exc).__name__}: {exc}）→ 认峰改用未扣减曲线")
        return sp.edf_frames(profiles, *q_range), mask, "unsub（扣减失败）", None
    return sp.edf_frames(subs, *q_range), mask, f"sub（全局 buffer {sp.fmt_ranges(bufs)}）", subs


def detect_and_report(profiles, series, st, out, prefix, args):
    """认峰 + 把判定画成图/写成表（"代码判不了就看图"的那份依据）。返回 det。"""
    if sp is None:
        raise SystemExit("认峰需要同目录的 sec_peaks.py（没 import 成功）")
    log("")
    log("=== 洗脱峰检测（多峰识别）===")
    y, mask, src, _subs = chromatogram_and_mask(profiles, series, st, args)
    kw = peak_kwargs(args)
    if args.peak_ranges:
        det = sp.manual_peaks(y, sp.parse_ranges(args.peak_ranges), mask=mask,
                              buffer_mode=kw["buffer_mode"], buffer_min_len=kw["buffer_min_len"],
                              note="人工区间（--peak-ranges）")
    else:
        det = sp.detect_peaks(y, mask=mask, **kw)
    if args.peak_buffers and det["peaks"]:
        sp.apply_manual_buffers(det["peaks"], args.peak_buffers, len(y))
    det.update(timestamp=time.strftime("%Y-%m-%d %H:%M:%S"), chromatogram=src,
               input=os.path.abspath(os.path.expanduser(args.work_dir or args.series_dir)),
               prefix=prefix, mask_frames=int(np.sum(mask)) if mask is not None else 0)
    for line in sp.summary_lines(det):
        log("  " + line)
    reasons = sp.needs_visual_check(det)
    det["visual_check_reasons"] = reasons
    if reasons:
        log("  ⚠️ 必须看图确认（自动判定在这几处可能不对）：")
        for r in reasons:
            log("     · " + r)
    else:
        log("  ✅ 自动判定没有触发需要人眼复核的条件。")
    sp.write_peaks_csv(os.path.join(out, "tables", "sec_peaks.csv"), det)
    sp.write_peaks_json(os.path.join(out, "series", "sec_peaks.json"), det)
    if det["peaks"]:
        sp.render_overview(os.path.join(out, "series", "sec_peaks.png"), det)
        sp.render_zoom(os.path.join(out, "series", "sec_peaks_zoom.png"), det)
        log("  图：series/sec_peaks.png（总览：峰窗+buffer+逐帧参数）、"
            "series/sec_peaks_zoom.png（逐峰放大）")
    if args.peak_sweep:
        path = os.path.join(out, "series", "sec_peaks_sweep.png")
        sp.render_sweep(path, y, sp.sweep_combos(args.peak_min_prominence), mask=mask,
                        q_range=kw["q_range"])
        log(f"  灵敏度扫描图：series/sec_peaks_sweep.png（阈值不同 → 峰数不同，人眼挑一组）")
    return det


def run_downstream(profiles, series, sample_profile, out, prefix, args, st, steps,
                   det=None, peak=None):
    """一个样品（或**一个峰**）在 ``out`` 里跑完下游并落盘 meta + README。

    Guinier → IFT → MW → 形状重建 → RAW 报告 → run_meta.json → README.md。
    """
    rows_g, report_profiles, ifts, rows_i = [], [], [], []
    ift = dmax = denss_ift_obj = None
    if "guinier" in steps and sample_profile is not None:
        rows_g, report_profiles = step_guinier(sample_profile, st, out, args)
    if "ift" in steps and sample_profile is not None:
        ift, ifts, rows_i, dmax, denss_ift_obj = step_ift(sample_profile, st, out, prefix,
                                                          args.atsas_dir, args.ift_dmax)
    if "mw" in steps and sample_profile is not None:
        step_mw(sample_profile, st, out, ift, args.atsas_dir)
    shape_rows = {"denss": [], "dammif": [], "damaver": []}
    if "shape" in steps and sample_profile is not None:
        shape_rows = step_shape(ift, ifts, out, prefix, args, args.atsas_dir,
                                model_ift=denss_ift_obj) or shape_rows
    if "report" in steps and series is not None:
        step_report(profiles, sample_profile, report_profiles, ifts, series, out, prefix)
        try:
            raw.save_series(series, f"{prefix}_series.hdf5", os.path.join(out, "series"))
        except Exception as exc:
            log(f"  ！series 存盘失败：{type(exc).__name__}: {exc}")

    meta = dict(prefix=prefix, args={k: v for k, v in vars(args).items()},
                denss_rows=shape_rows.get("denss", []),
                dammif_rows=shape_rows.get("dammif", []),
                damaver_row=shape_rows.get("damaver", []),
                prefix0=prefix, input=os.path.abspath(os.path.expanduser(args.work_dir or args.series_dir)),
                n_frames=len(profiles), cfg=os.path.abspath(args.cfg),
                settings=dict(ImageHdrFormat=st.get('ImageHdrFormat'),
                              EnableNormalization=st.get('EnableNormalization'),
                              NormalizationList=st.get('NormalizationList'),
                              ATSASDir=st.get('ATSASDir')),
                steps=steps, model_engine=args.model_engine, n_models=args.n_models,
                symmetry=args.symmetry, guinier_rows=rows_g, ift_rows=rows_i, dmax=dmax,
                aatsas_dir=args.atsas_dir, timestamp=time.strftime("%Y-%m-%d %H:%M:%S"))
    if peak is not None:
        meta["peak"] = {k: v for k, v in peak.items()}
        meta["is_multi_peak_run"] = True
    if det is not None:
        meta["detected_peaks"] = [{k: p.get(k) for k in
                                   ("index", "apex", "window", "prominence", "snr", "width",
                                    "buffer_ranges", "subdir")} for p in det.get("peaks", [])]
        meta["chromatogram"] = det.get("chromatogram")
    with open(os.path.join(out, "run_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2, ensure_ascii=False, default=str)
    # 最后一步：把结果目录写成给人看的 README.md（由 write-saxs-results-readme 统一生成）
    try:
        p, n_keys = write_results_readme(out, "sec", getattr(args, "readme_script", None))
        log(f"  → 结果说明：{os.path.relpath(p, out)}（先看这份，再看其它文件；{n_keys} 条关键数字）")
    except Exception as exc:
        log(f"  ！README.md 未生成：{type(exc).__name__}: {exc}")
    return shape_rows, rows_g, rows_i, dmax, sample_profile


def combine_per_peak_params(per_peak, n):
    """把每个峰**用自己 buffer 算出来**的 Rg/I(0)/MW 拼成整条系列的数组（峰值外填 NaN）。

    这张图就是"平台判据"的可视化：每个峰在自己的 buffer 下应有一段 Rg/MW 平台。
    """
    if not per_peak:
        return None
    out = {"rg": {}, "i0": {}, "mw": {}}

    def place(vals, a, b):
        v = np.asarray(vals, dtype=float)[a:b + 1].copy()
        v[~np.isfinite(v) | (v <= 0)] = np.nan            # RAW 用 -1/0 当"没算出来"
        full = np.full(n, np.nan)
        full[a:b + 1] = v
        return full

    for item in per_peak:
        p = item["peak"]
        bufs = p.get("buffer_ranges") or []
        a = max(0, min([b[0] for b in bufs] + [p["window"][0]]))
        b = min(n - 1, max([b[1] for b in bufs] + [p["window"][1]]))
        lab = f"P{p['index']+1}"
        out["rg"][f"Rg {lab} (own buffer)"] = place(item["rg"], a, b)
        out["i0"][f"I(0) {lab}"] = place(item["i0"], a, b)
        out["mw"][f"Vc {lab}"] = place(item["vc"], a, b)
        out["mw"][f"Vp {lab}"] = place(item["vp"], a, b)
    return out


def read_mw_vc(pdir):
    """从峰子目录的 tables/mw.csv 里取 Vc 分子量（kDa）——只读产物，不重算。"""
    path = os.path.join(pdir, "tables", "mw.csv")
    if not os.path.exists(path):
        return None
    try:
        with open(path) as fh:
            for row in csv.DictReader(fh):
                if str(row.get("method", "")).lower() == "vc":
                    return float(row.get("mw"))
    except Exception:
        return None
    return None


def run_per_peak(profiles, series, st, out, prefix, det, args, steps):
    """多峰：**每个峰一个子目录**，各自在本峰的 buffer + 峰窗上重跑扣减与下游分析。"""
    if sp is None:
        raise SystemExit("多峰分析需要同目录的 sec_peaks.py（没 import 成功）")
    log("")
    log(f"=== 逐峰分析（{len(det['peaks'])} 个峰）===")
    peaks_root = os.path.join(out, "peaks")
    per_peak, summary = [], []
    for p in det["peaks"]:
        subdir = f"peak{p['index']+1:02d}_apex{p['apex']:05d}"
        pdir = os.path.join(peaks_root, subdir)
        for d in ("profiles", "tables", "ifts", "series", "models", "reports"):
            os.makedirs(os.path.join(pdir, d), exist_ok=True)
        p["subdir"] = os.path.join("peaks", subdir)
        # 01_integrated 是**共享**的（图像只积分一次）→ 相对符号链接过去，别重复占盘
        link = os.path.join(pdir, "profiles", "01_integrated")
        if not os.path.exists(link):
            try:
                os.symlink(os.path.relpath(os.path.join(out, "profiles", "01_integrated"),
                                           os.path.join(pdir, "profiles")), link)
            except OSError as exc:
                log(f"  ！共享 01_integrated 的符号链接建不了（{exc}）")
        pprefix = f"{prefix}_P{p['index']+1}"
        bufs = p.get("buffer_ranges") or None
        bl = None
        if args.baseline != "none" and bufs:
            # 多峰 + 基线校正：用**这个峰自己的** buffer 段当初末锚点（不是写死的 0,20;1800,1990）
            bl = f"{bufs[0][0]},{bufs[0][1]};{bufs[-1][0]},{bufs[-1][1]}"
        log(f"---- P{p['index']+1}：apex 帧 {p['apex']}，sample 窗 {p['window']}"
            f"（{p['window'][1]-p['window'][0]+1} 帧），buffer {sp.fmt_ranges(bufs or [])}"
            f" → {os.path.relpath(pdir, out)}")
        try:
            series_p, sample_p, ranges = step_ranges_and_subtraction(
                profiles, st, pdir, pprefix, args, buffers=bufs, sample=p["window"],
                baseline_ranges=bl, det=det,
                peak=dict(index=p["index"], apex=p["apex"], window=list(p["window"]),
                          buffer_ranges=bufs, subdir=p["subdir"], prominence=p["prominence"],
                          snr=p["snr"], width=p["width"]))
        except SystemExit as exc:
            log(f"  ！P{p['index']+1} 扣减/取样失败 → 跳过这个峰：{exc}")
            p["failed"] = str(exc)
            continue
        except Exception as exc:
            log(f"  ！P{p['index']+1} 扣减/取样出错 → 跳过这个峰：{type(exc).__name__}: {exc}")
            p["failed"] = f"{type(exc).__name__}: {exc}"
            continue
        try:
            shape_rows, rows_g, rows_i, dmax, _ = run_downstream(
                profiles, series_p, sample_p, pdir, pprefix, args, st, steps, det=det, peak=p)
            summary.append((p, shape_rows, rows_g, rows_i, dmax))
        except Exception as exc:
            log(f"  ！P{p['index']+1} 下游失败：{type(exc).__name__}: {exc}")
            summary.append((p, None, [], [], None))
        try:
            rg_a, _ = (np.asarray(x, float) for x in series_p.getRg())
            i0_a, _ = (np.asarray(x, float) for x in series_p.getI0())
            vc_a, _ = (np.asarray(x, float) for x in series_p.getVcMW())
            vp_a, _ = (np.asarray(x, float) for x in series_p.getVpMW())
            per_peak.append(dict(peak=p, rg=rg_a, i0=i0_a, vc=vc_a, vp=vp_a))
        except Exception as exc:
            log(f"  ！P{p['index']+1} 逐帧参数收集失败：{type(exc).__name__}: {exc}")

    # 顶层总览图重画一次：把每个峰自己的 Rg/I(0)/MW 拼在同一条系列上（平台判据一眼看）
    try:
        arrs = combine_per_peak_params(per_peak, len(det["y_raw"]))
        sp.render_overview(os.path.join(out, "series", "sec_peaks.png"), det, series_arrays=arrs,
                           title=f"SEC chromatogram + {len(det['peaks'])} peaks "
                                 f"(colour = peak; lower panels = each peak's own Rg / I(0) / MW)")
    except Exception as exc:
        log(f"  ！总览图重画失败：{type(exc).__name__}: {exc}")

    # 顶层：索引 README（每个峰的结果说明在它自己的子目录里）
    det["timestamp"] = time.strftime("%Y-%m-%d %H:%M:%S")
    extra = ["## 3. 产物布局", "",
             f"· **共享**：`profiles/01_integrated/`（{len(profiles)} 帧，图像只积分一次）、"
             f"`tables/sec_peaks.csv`（峰表）、`series/sec_peaks*.png`（判定图）、"
             f"`series/sec_peaks.json`（阈值与全部实测数字）",
             "· **每个峰一个子目录**（结构 = 该峰单独跑一遍的完整产物）："
             "`README.md`、`run_meta.json`、`profiles/`（02_buffer、03_subtracted、04_sample、"
             "06_guinier）、`tables/`、`ifts/`、`models/`、`reports/`、`series/`",
             "· 子目录里的 `profiles/01_integrated` 是**相对符号链接**（同一份积分结果，不重复占盘）。",
             "",
             "## 4. 每个峰一行速览", "",
             "| 峰 | 子目录 | 主拟合 Rg (Å) | MW (Vc, kDa) | 形状重建 |",
             "|---|---|---|---|---|"]
    for p, shape_rows, rows_g, rows_i, dmax in summary:
        rg = ''
        for r in rows_g or []:
            if r and r[0] == 'auto':
                rg = f"{float(r[1]):.1f}"
        mwv = read_mw_vc(os.path.join(out, p.get("subdir", "")))
        shape = []
        if (shape_rows or {}).get("denss"):
            shape.append("DENSS 电子云 ✓")
        if (shape_rows or {}).get("dammif"):
            shape.append(f"DAMMIF×{len(shape_rows['dammif'])} 珠模")
        if (shape_rows or {}).get("damaver"):
            shape.append("DAMAVER ✓")
        extra.append(f"| **P{p['index']+1}** | `{p['subdir']}` | {rg or '—'} | "
                     f"{('%.0f' % mwv) if mwv else '—'} | {'、'.join(shape) or '—'} |")
    if det.get("visual_check_reasons"):
        extra += ["", "## 5. 看图之前别引用数字", "",
                  "自动判定触发了下面这些条件（见 `series/sec_peaks.png` 与 "
                  "`series/sec_peaks_zoom.png`），确认或修正之后再引用各峰的拟合结果：", ""]
        extra += [f"· {r}" for r in det["visual_check_reasons"]]
    try:
        sp.write_peaks_readme(out, prefix, det,
                              png_overview=os.path.join(out, "series", "sec_peaks.png"),
                              png_zoom=os.path.join(out, "series", "sec_peaks_zoom.png"),
                              extra_lines=extra)
        log("  → 顶层索引：README.md（峰表 + 每个峰的产物入口）")
    except Exception as exc:
        log(f"  ！顶层 README 失败：{type(exc).__name__}: {exc}")

    meta = dict(prefix=prefix, args={k: v for k, v in vars(args).items()},
                multi_peak=True, n_peaks=len(det["peaks"]),
                input=os.path.abspath(os.path.expanduser(args.work_dir or args.series_dir)),
                n_frames=len(profiles), cfg=os.path.abspath(args.cfg), steps=steps,
                chromatogram=det.get("chromatogram"),
                peaks=[{k: p.get(k) for k in ("index", "apex", "window", "prominence", "snr",
                                              "width", "buffer_ranges", "subdir", "manual",
                                              "failed")} for p in det["peaks"]],
                visual_check_reasons=det.get("visual_check_reasons", []),
                mask_frames=det.get("mask_frames"), noise_sigma=det.get("noise_sigma"),
                timestamp=time.strftime("%Y-%m-%d %H:%M:%S"))
    with open(os.path.join(out, "run_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2, ensure_ascii=False, default=str)
    return det


# ----------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=FMT)
    ap.add_argument("--work-dir", default=None, help="emit-bl19u2-header-txt.py 产出的符号链接目录（tif+txt）")
    ap.add_argument("--series-dir", default=None, help="或直接给图像目录（要求每帧旁边就有同名 .txt）")
    ap.add_argument("--out-dir", required=True, help="产物根目录")
    ap.add_argument("--cfg", required=True, help="RAW 设置文件（线站下机的 <日期>.cfg）")
    ap.add_argument("--prefix", default=None, help="产物文件名前缀（默认取目录名）")
    ap.add_argument("--steps", default=",".join(STEPS), help=f"逗号分隔，可选 {STEPS}")
    ap.add_argument("--buffer-range", default=None,
                    help="手工指定 buffer 区 'start,end'，可给多段 's,e;s,e'（帧号 0 基；跨度大时给峰前+峰后两段）")
    ap.add_argument("--sample-range", default=None, help="手工指定样品区 'start,end'（帧号，0 基）")
    ap.add_argument("--baseline", choices=["none", "linear", "integral"], default="none",
                    help="扣减后是否再做基线校正（RAW 的 Linear/Integral）")
    ap.add_argument("--baseline-ranges", default="0,20;1800,1990",
                    help="--baseline linear 时的起止区间 's0,s1;e0,e1'")
    ap.add_argument("--ift-dmax", type=float, default=None,
                    help="显式 Dmax（A）→ 用 RAW 原生 DIFT 生成 DENSS 的输入 IFT。IFT/DENSS 被低 q 拖坏时用"
                         "（经验起手：Dmax ≈ 3×Guinier Rg）；不给则用 BIFT 自动定出的 IFT")
    ap.add_argument("--trim-qmax", type=float, default=None,
                    help="下游分析前丢掉 q 高于此值的点（1/A）——判断\"高 q 是模型表达不出来还是数据不好\"时用："
                         "同一条曲线降 q_max 后 χ² 若明显回落，说明是模型/信息量问题，不是样品问题")
    ap.add_argument("--trim-qmin", type=float, default=None,
                    help="下游分析（Guinier 表/IFT/MW）前丢掉 q 低于此值的点（1/A）；低 q 被寄生散射污染时用")
    ap.add_argument("--guinier-ranges", default=None,
                    help="手工指定多区间 'qlo:qhi,qlo:qhi'（1/A）；默认以 auto Rg 为锚铺 qRg 阶梯")
    ap.add_argument("--model-engine", choices=["auto", "denss", "dammif", "both", "none"], default="auto",
                    help="auto=电子云(DENSS)总跑 + 有 ATSAS 时再加珠模(DAMMIF)；二者可分开指定")
    ap.add_argument("--denss-mode", choices=["Fast", "Slow", "Custom"], default="Fast",
                    help="DENSS 模式（Fast 出得快、Slow 收敛更好、耗时更长）")
    ap.add_argument("--n-models", type=int, default=4, help="DAMMIF 模型数")
    ap.add_argument("--readme-script", default=None,
                    help="write-saxs-results-readme 的 scripts/write-readme.py 路径；"
                         "默认自动找（同类目下的同名技能目录），找不到只在日志里告警")
    ap.add_argument("--dammif-mode", choices=["Fast", "Slow", "Custom"], default="Fast",
                    help="DAMMIF 模式（Fast 出得快；Slow 更彻底、耗时显著更长）")
    ap.add_argument("--model-format", choices=["cif", "pdb"], default="cif",
                    help="DAMMIF/DAMAVER 输出的模型格式（pdb 便于直接看结构，cif 是 RAW 默认）")
    ap.add_argument("--symmetry", default="P1", help="DAMMIF 对称性")
    ap.add_argument("--atsas-dir", default=None, help="ATSAS bin 目录（装了就传，RAW 的 GNOM/DAMMIF 需要）")
    ap.add_argument("--no-header-normalization", action="store_true",
                    help="该系列没有逐帧 BL19U2 header txt 时用：不启用逐帧归一化（ImageHdrFormat=None）")
    ap.add_argument("--limit", type=int, default=None, help="只用前 N 帧（试跑）")
    # ---- 多峰（≥2 个洗脱峰）识别与逐峰分析。见 sec_peaks.py 与 SKILL.md「多峰 SEC 数据」
    ap.add_argument("--multi-peak", choices=["auto", "off", "always"], default=None,
                    help="auto=识别到 ≥2 个洗脱峰就逐峰分析（默认）；off=不做分峰（只出检测图/表）；"
                         "always=即使只有 1 个峰也用它自己的 buffer+峰窗（不走 RAW 自动选区间）")
    ap.add_argument("--peak-ranges", default=None,
                    help="**人工**指定峰窗口 'lo,hi;lo,hi'（0 基闭区间）→ 不看阈值，直接按它分峰"
                         "（视觉复核后的落点；先用 find-sec-peaks.py 看图定区间）")
    ap.add_argument("--peak-buffers", default=None,
                    help="**人工**指定每个峰的 buffer 's,e;s,e|s,e'（'|' 分峰、';' 分段）；不给则按峰自动取")
    ap.add_argument("--peak-q-range", default=None,
                    help="色谱图取的 q 窗口 'qlo,qhi'（1/A），默认 %g,%g——低 q 对组分敏感、"
                         "高 q 对噪声敏感" % sp.DEF_Q_RANGE if sp else "")
    ap.add_argument("--peak-chromatogram", choices=["sub", "unsub"], default="sub",
                    help="用哪条曲线认峰：sub=扣减后（组分分得开，但依赖 buffer 选得对）；"
                         "unsub=未扣减（不依赖 buffer，基线漂移更大）")
    ap.add_argument("--peak-min-prominence", type=float, default=None,
                    help="峰幅度下限（相对最高峰，0-1）")
    ap.add_argument("--peak-min-snr", type=float, default=None, help="峰高/基线噪声（1.4826*MAD）下限")
    ap.add_argument("--peak-min-height", type=float, default=None, help="峰高下限（相对最高峰，0-1）")
    ap.add_argument("--peak-min-width", type=int, default=None, help="半高宽下限（帧）")
    ap.add_argument("--peak-max-n", type=int, default=None, help="最多分析几个峰（多的按幅度排掉并告警）")
    ap.add_argument("--peak-baseline-window", type=int, default=None,
                    help="滚动中位基线窗口（帧）；默认按帧数自动（约 n/10，51–401）")
    ap.add_argument("--peak-smooth-window", type=int, default=None, help="savgol 平滑窗口（帧）；默认 min(51, n/2)")
    ap.add_argument("--peak-window", choices=["half-height", "valley"], default="half-height",
                    help="峰窗取法：half-height=半高宽（保守）；valley=到相邻峰谷底（整个峰）")
    ap.add_argument("--peak-buffer", choices=["local", "global"], default="local",
                    help="每个峰的 buffer 取法：local=紧邻的前后两段（默认，就是'峰前+峰后'）；"
                         "global=所有峰共用全部非峰区")
    ap.add_argument("--peak-buffer-min", type=int, default=None, help="单侧 buffer 段的帧数下限（帧）")
    ap.add_argument("--peak-buffer-max", type=int, default=0,
                    help="单侧 buffer 段的帧数上限（帧；0=不限制）——基线在漂移时用，把每段截到"
                         "紧邻峰的那一侧（水平匹配优先，实测 BSA 峰前 −0.10 / 峰后 +0.02）")
    ap.add_argument("--peak-sweep", action="store_true",
                    help="额外出阈值灵敏度扫描图 series/sec_peaks_sweep.png（供人眼挑参数）")
    args = ap.parse_args()
    if sp is not None:                       # 没给就用 sec_peaks 的默认值（单一来源）
        args.multi_peak = args.multi_peak or "auto"
        args.peak_q_range = args.peak_q_range or "%g,%g" % sp.DEF_Q_RANGE
        args.peak_min_prominence = sp.DEF_MIN_PROMINENCE if args.peak_min_prominence is None else args.peak_min_prominence
        args.peak_min_snr = sp.DEF_MIN_SNR if args.peak_min_snr is None else args.peak_min_snr
        args.peak_min_height = sp.DEF_MIN_HEIGHT if args.peak_min_height is None else args.peak_min_height
        args.peak_min_width = sp.DEF_MIN_WIDTH if args.peak_min_width is None else args.peak_min_width
        args.peak_max_n = sp.DEF_MAX_PEAKS if args.peak_max_n is None else args.peak_max_n
        args.peak_buffer_min = sp.DEF_BUFFER_MIN_LEN if args.peak_buffer_min is None else args.peak_buffer_min
    else:
        args.multi_peak = args.multi_peak or "off"

    src = args.work_dir or args.series_dir
    if not src:
        raise SystemExit("必须给 --work-dir 或 --series-dir")
    files = sorted(glob.glob(os.path.join(os.path.expanduser(src), "*.tif")))
    if args.limit:
        files = files[:args.limit]
    if not files:
        raise SystemExit(f"{src} 下没有 tif")
    out = os.path.abspath(os.path.expanduser(args.out_dir))
    prefix = args.prefix or os.path.basename(os.path.abspath(src).rstrip("/"))
    steps = [s.strip() for s in args.steps.split(",") if s.strip()]
    for d in ("profiles", "tables", "ifts", "series", "models", "reports"):
        os.makedirs(os.path.join(out, d), exist_ok=True)

    log(f"输入 {len(files)} 帧 ← {src}")
    log(f"输出 {out}（prefix={prefix}）| steps={steps}")
    st = make_settings(args.cfg, args.atsas_dir,
                       header_normalization=not args.no_header_normalization)
    log(f"settings: ImageHdrFormat={st.get('ImageHdrFormat')} EnableNormalization={st.get('EnableNormalization')} "
        f"NormalizationList={st.get('NormalizationList')} ATSASDir={st.get('ATSASDir')}")

    profiles = step_integrate(files, st, out, prefix, args.limit) if "integrate" in steps else None
    if profiles is None:
        log("！跳过积分则后续步骤无法进行（本脚本以 RAW API 在内存里传对象）")
        return

    series = raw.profiles_to_series(profiles, st)

    # ---- 认峰（--multi-peak off 时也出图/表：知道"其实只有一个峰/有两个峰"本身就有用）
    det = None
    if sp is not None and ("peaks" in steps or args.multi_peak != "off"):
        det = detect_and_report(profiles, series, st, out, prefix, args)
    n_peaks = len(det["peaks"]) if (det and det.get("ok")) else 0
    need_series = bool({"series", "guinier", "ift", "mw", "shape", "report"} & set(steps))
    if not need_series:
        if det is not None:
            with open(os.path.join(out, "run_meta.json"), "w") as fh:
                json.dump(dict(prefix=prefix, steps=steps, input=src, n_frames=len(files),
                               cfg=os.path.abspath(args.cfg), peaks_only=True,
                               n_peaks=n_peaks, chromatogram=det.get("chromatogram"),
                               visual_check_reasons=det.get("visual_check_reasons", []),
                               timestamp=time.strftime("%Y-%m-%d %H:%M:%S")),
                          fh, indent=2, ensure_ascii=False, default=str)
            try:                      # 只跑认峰时也给一份索引（人要能一眼看到峰表与判定图）
                for p in det["peaks"]:
                    p.setdefault("subdir", "（只认峰，未逐峰分析；去掉 --steps 限制重跑）")
                sp.write_peaks_readme(out, prefix, det,
                                      png_overview=os.path.join(out, "series", "sec_peaks.png"),
                                      png_zoom=os.path.join(out, "series", "sec_peaks_zoom.png"))
            except Exception as exc:
                log(f"  ！索引 README 失败：{type(exc).__name__}: {exc}")
            log(f"steps={steps} 里没有 series/下游 → 只出了峰检测产物（{out}/series/sec_peaks*.png、"
                f"{out}/tables/sec_peaks.csv）；确认峰划分后再跑全流程")
        return

    # ---- 分不分峰：--multi-peak auto（≥2 才分）/ always（1 个也分）/ off（不分）
    split = bool(det and det.get("ok") and n_peaks >= 1
                 and args.multi_peak != "off" and (n_peaks >= 2 or args.multi_peak == "always"))
    if args.peak_ranges:                      # 人工给了区间 → 一定按它分
        split = bool(det and det.get("ok") and n_peaks >= 1)
    if split:
        run_per_peak(profiles, series, st, out, prefix, det, args, steps)
    else:
        if args.multi_peak != "off":
            if n_peaks == 1:
                log("  （只认到 1 个洗脱峰 → 走原来的单峰流程：RAW 自动 buffer + 自动样品区；"
                    "想让它也用自己那一对 buffer/峰窗，加 --multi-peak always）")
            elif n_peaks == 0:
                log("  （没认到洗脱峰 → 走原来的单峰流程，手工 --buffer-range/--sample-range 仍然可用）")
        series_s, sample_profile, _ranges = step_ranges_and_subtraction(profiles, st, out, prefix,
                                                                       args, det=det)
        run_downstream(profiles, series_s, sample_profile, out, prefix, args, st, steps, det=det)
    log(f"完成。清单见 {out}/run_meta.json")


if __name__ == "__main__":
    main()
