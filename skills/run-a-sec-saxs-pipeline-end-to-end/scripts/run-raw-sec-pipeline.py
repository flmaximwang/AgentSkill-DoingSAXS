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


FMT = argparse.ArgumentDefaultsHelpFormatter
STEPS = ["integrate", "series", "guinier", "ift", "mw", "shape", "report"]


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


def step_ranges_and_subtraction(profiles, st, out, prefix, args):
    """buffer 区 → 扣减 →（可选）基线 → sample 区。返回样品平均曲线与 series。"""
    series = raw.profiles_to_series(profiles, st)

    # ---- buffer 区
    if args.buffer_range:
        b_rng = [[int(x) for x in part.split(",")] for part in args.buffer_range.split(";")]
        clip_ranges(b_rng, len(profiles), "buffer")
        b_ok = True
    else:
        ok, s, e = raw.find_buffer_range(series)
        if s is None or e is None:
            raise SystemExit(f"RAW 没找到 buffer 区（success={ok}）：系列里可能没有可识别的洗脱峰"
                             "（例如只取了峰前的帧）。用 --buffer-range start,end 手工指定（0 基帧号）。")
        b_ok, b_rng = bool(ok), [[int(s), int(e)]]
    log(f"buffer 区 {b_rng}（{'给定' if args.buffer_range else 'RAW 自动找'}, success={b_ok}）")
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
                start_range, end_range = args.baseline_ranges.split(";")
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
    if args.sample_range:
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
    log(f"样品区 {s_rng[0]}（{'给定' if args.sample_range else 'RAW 自动找'}, success={s_ok}, profile_type={pt}）")
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
        plot_series(out, frames, int_a, rg_a, i0_a, vc_a, vp_a, b_rng[0], s_rng[0])
    except Exception as exc:
        log(f"  ！series 图失败：{type(exc).__name__}: {exc}")

    ranges = dict(buffer=b_rng, sample=s_rng, profile_type=pt,
                  sample_range_success=s_ok, buffer_range_success=b_ok)
    with open(os.path.join(out, "series", "ranges.json"), "w") as fh:
        json.dump(ranges, fh, indent=2)
    return series, sample_profile, ranges


def plot_series(out, frames, int_a, rg_a, i0_a, vc_a, vp_a, b_rng, s_rng):
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
    args = ap.parse_args()

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

    series = sample_profile = None
    report_profiles, ifts, rows_g = [], [], []
    if "series" in steps or {"guinier", "ift", "mw", "shape", "report"} & set(steps):
        series, sample_profile, ranges = step_ranges_and_subtraction(profiles, st, out, prefix, args)
    if "guinier" in steps and sample_profile is not None:
        rows_g, report_profiles = step_guinier(sample_profile, st, out, args)
    ift, ifts, rows_i, dmax, denss_ift_obj = (None, [], [], None, None)
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
                prefix0=prefix, input=src, n_frames=len(files), cfg=os.path.abspath(args.cfg),
                settings=dict(ImageHdrFormat=st.get('ImageHdrFormat'),
                              EnableNormalization=st.get('EnableNormalization'),
                              NormalizationList=st.get('NormalizationList'),
                              ATSASDir=st.get('ATSASDir')),
                steps=steps, model_engine=args.model_engine, n_models=args.n_models,
                symmetry=args.symmetry, guinier_rows=rows_g, ift_rows=rows_i, dmax=dmax,
                aatsas_dir=args.atsas_dir, timestamp=time.strftime("%Y-%m-%d %H:%M:%S"))
    with open(os.path.join(out, "run_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2, ensure_ascii=False, default=str)
    # 最后一步：把结果目录写成给人看的 README.md（由 write-saxs-results-readme 统一生成）
    try:
        p, n_keys = write_results_readme(out, "sec", getattr(args, "readme_script", None))
        log(f"  → 结果说明：{os.path.relpath(p, out)}（先看这份，再看其它文件；{n_keys} 条关键数字）")
        log(f"     验收：python <{README_SKILL}>/scripts/verify-results-folder.py {out}")
    except Exception as exc:
        log(f"  ！README.md 未生成：{type(exc).__name__}: {exc}")
    log(f"完成。清单见 {out}/run_meta.json")


if __name__ == "__main__":
    main()
