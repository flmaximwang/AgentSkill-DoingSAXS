"""峰内梯度（把一条 SEC 洗脱峰当**稀释序列**用）。

为什么可以这么用：SEC 洗脱峰的**浓度沿帧变化**——峰顶最浓、两翼渐稀。所以同一条峰里
不同位置的帧天然构成一条 c 递减的序列，且**基质相同**（同一 buffer、同一批样品），
不需要另配样品。于是可以做三件在原流程里做不到的事：

1. **Rg / MW 的浓度依赖**：理想单分散体系 Rg 不随 c 变；Rg 随 c 降 → 排斥型相互作用
   （或浓度依赖的解离）；Rg 随 c 升 → 吸引/自缔合（低聚态）。
2. **c→0 外推**：把 Rg 外推到无限稀释（去掉相互作用），得到"真"Rg；SEC 峰内浓度未知，
   但**相对**浓度足够（外推只看相对标度，截距即为 c=0 时的值）。
3. **上升/下降两翼互检**：同一 c 在峰的两侧各出现一次（asc/desc），两者的 Rg/MW 必须在
   误差内一致——否则说明有基线漂移、组分变化或扣减残差，比单点判据更能暴露问题。

口径（都是为了"可比较"，不是自写拟合）：

- 所有切片都调 **RAW 自己的** ``set_sample_range`` + ``guinier_fit``；**用同一条 q 区间**
  （默认取主分析的 auto 区间）——换区间会让 Rg 差异变成"取点差异"，那就没法比了。
- 切片的"浓度"用**同一张色谱图**（扣减后低 q 窗口积分）在切片内的均值 / 峰顶值，是
  **相对**浓度（SEC 浓度未知，绝对刻度不适用）。

两种切片口径：``frames``（默认，峰窗内**连续**帧块，从一侧扫到另一侧；块内组分几乎不变）、
``height``（按相对峰高分带）。
"""

import csv
import json
import os

import numpy as np

DEF_N_SLICES = 5
DEF_MIN_FRAMES = 3
DEF_MIN_CONC = 0.08          # 相对浓度下限：更稀的空切片信噪撑不住 Rg
DEF_MODE = "frames"

CSV_FIELDS = ["idx", "side", "label", "frame_start", "frame_end", "n_frames", "conc_rel",
              "height_frac", "rg", "rg_err", "i0", "i0_err", "q_min", "q_max", "qrg_min",
              "qrg_max", "r2", "mw_vc", "mw_vp", "mw_vc_err", "verdict"]


# ----------------------------------------------------------------- 切片划分
def _peak_window(peak, n_frames_series):
    lo, hi = int(peak["window"][0]), int(peak["window"][1])
    return max(0, lo), min(n_frames_series - 1, hi)


def slices_for_peak(det, peak, n=DEF_N_SLICES, mode=DEF_MODE, min_frames=DEF_MIN_FRAMES,
                    min_conc=DEF_MIN_CONC, include_whole=True):
    """把峰窗切成若干"浓度切片"（返回 dict 列表，按浓度从高到低）。

    ``mode='frames'``：峰窗内 n 个**连续**帧块（等宽、覆盖整窗），块内组分几乎不变；
    ``mode='height'``：按相对峰高分带（带内帧可能不连续）。
    ``conc_rel`` = 该切片帧上（基线校正）色谱强度的均值 / 峰顶值 → **相对**浓度。
    """
    y = np.asarray(det["y_corr"], dtype=float)
    n_series = len(y)
    lo, hi = _peak_window(peak, n_series)
    apex = int(peak["apex"])
    top = float(y[apex]) if y[apex] > 0 else float(np.nanmax(y[lo:hi + 1]))
    out = []

    if include_whole and hi > lo:
        out.append(dict(label="whole-peak", side="all", window=[lo, hi], n_frames=hi - lo + 1,
                        conc_rel=float(np.nanmean(y[lo:hi + 1]) / top), height_frac=float(
                            np.nanmean(y[lo:hi + 1]) / top), note="峰窗整体（= 主分析用的样品区）"))

    if mode == "height":
        for k, (f_hi, f_lo) in enumerate(zip(np.linspace(1.0, min_conc, n + 1)[:-1],
                                             np.linspace(1.0, min_conc, n + 1)[1:])):
            idx = np.flatnonzero((y[lo:hi + 1] / top >= f_lo) & (y[lo:hi + 1] / top <= f_hi)) + lo
            if len(idx) < min_frames:
                continue
            out.append(dict(label=f"h{f_hi:.2f}-{f_lo:.2f}", side="band",
                            window=[int(idx.min()), int(idx.max())],
                            n_frames=int(len(idx)), conc_rel=float(np.nanmean(y[idx]) / top),
                            height_frac=float(np.nanmean(y[idx]) / top),
                            frames=[int(i) for i in idx]))
    else:
        width = hi - lo + 1
        if width <= 0:
            return out
        b = max(int(min_frames), int(round(width / max(1, n))))
        b = min(b, width)
        starts = np.linspace(lo, hi - b + 1, n).round().astype(int) if width > b else np.array([lo])
        for i, s in enumerate(sorted(set(int(x) for x in starts))):
            e = min(hi, s + b - 1)
            if e - s + 1 < min_frames:
                continue
            c = float(np.nanmean(y[s:e + 1]) / top)
            if c < min_conc:
                continue
            out.append(dict(label=f"f{s}-{e}", side="asc" if (s + e) / 2 <= apex else "desc",
                            window=[s, e], n_frames=e - s + 1, conc_rel=c, height_frac=c))

    out.sort(key=lambda d: -d["conc_rel"])
    for i, d in enumerate(out):
        d["idx"] = i
    return out


# ----------------------------------------------------------------- 逐切片拟合（全走 RAW）
def _qindex(profile, qval):
    q = np.asarray(profile.getQ())
    if len(q) == 0:
        return 0
    return int(np.argmin(np.abs(q - float(qval))))


def analyze_slices(series, st, slices, raw, qrange=None, mw=True):
    """对每个切片：RAW 取平均 → RAW 的 Guinier（**同一条 q 区间**）→ 可选 MW。只读结果。"""
    rows = []
    for sl in slices:
        rng = [[int(sl["window"][0]), int(sl["window"][1])]]
        prof = raw.set_sample_range(series, rng, profile_type="sub")
        row = dict(sl)
        try:
            if qrange:
                i0i, i1i = _qindex(prof, qrange[0]), _qindex(prof, qrange[1])
                if i1i <= i0i:
                    i1i = min(len(prof.getQ()) - 1, i0i + 2)
                fit = raw.guinier_fit(prof, i0i, i1i, settings=st)
            else:
                fit = raw.auto_guinier(prof, settings=st)
            (row["rg"], row["i0"], row["rg_err"], row["i0_err"], row["q_min"], row["q_max"],
             row["qrg_min"], row["qrg_max"], row["r2"]) = (float(x) for x in fit[:9])
        except Exception as exc:
            row["error"] = f"{type(exc).__name__}: {exc}"
        if mw and "error" not in row:
            for key, fn in (("mw_vc", raw.mw_vc), ("mw_vp", raw.mw_vp)):
                try:
                    res = fn(prof, settings=st)
                    row[key] = float(res[0])
                    if key == "mw_vc" and len(res) > 2:
                        row["mw_vc_err"] = float(res[2])
                except Exception:
                    pass
        rows.append(row)
    return rows


# ----------------------------------------------------------------- 浓度依赖判词
def _linfit(x, y, w=None):
    """加权最小二乘 y = a + b·x；返回 (a, b, se_a, se_b, r2)。x/y 已去过 NaN。"""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    if n < 2 or np.allclose(x, x[0]):
        return (float(np.nanmean(y)) if n else np.nan, 0.0, np.nan, np.nan, np.nan)
    w = np.ones(n) if w is None else np.asarray(w, float)
    w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
    if w.sum() <= 0:
        w = np.ones(n)
    sw, swx = w.sum(), (w * x).sum()
    swy, swxx, swxy = (w * y).sum(), (w * x * x).sum(), (w * x * y).sum()
    den = sw * swxx - swx ** 2
    if abs(den) < 1e-12:
        return float(swy / sw), 0.0, np.nan, np.nan, np.nan
    b = (sw * swxy - swx * swy) / den
    a = (swy - b * swx) / sw
    res = y - (a + b * x)
    dof = max(1, n - 2)
    s2 = float((w * res ** 2).sum() / dof)
    se_b = float(np.sqrt(s2 * sw / den))
    se_a = float(np.sqrt(s2 * swxx / den))
    ss_tot = float((w * (y - np.nanmean(y)) ** 2).sum())
    r2 = float(1 - (w * res ** 2).sum() / ss_tot) if ss_tot > 0 else np.nan
    return float(a), float(b), se_a, se_b, r2


def gradient_stats(rows):
    """Rg / MW 对相对浓度的最小二乘 + c→0 外推 + 判词。只统计"可用"切片。"""
    use = [r for r in rows if "error" not in r and np.isfinite(r.get("rg", np.nan))
           and r.get("rg", -1) > 0 and r.get("r2", 0) >= 0.9 and r.get("side") in ("asc", "desc", "band", "all")]
    c = np.array([r["conc_rel"] for r in use], float)
    rg = np.array([r["rg"] for r in use], float)
    rge = np.array([abs(r.get("rg_err", 0.0)) for r in use], float)
    i0 = np.array([r.get("i0", np.nan) for r in use], float)
    mw = np.array([r.get("mw_vc", np.nan) for r in use], float)
    st = dict(n_used=len(use), n_total=len(rows))
    if len(use) < 2:
        st.update(verdict="🔴 不可判", reason=f"可用切片只有 {len(use)} 个（要 ≥2，建议重跑 --gradient-slices 调大或放宽 min-conc）")
        return st
    w = 1.0 / np.where(rge > 0, rge, np.nanmedian(rge) if np.median(rge) > 0 else 1.0)
    a, b, se_a, se_b, r2 = _linfit(c, rg, w)
    st.update(rg0=a, rg0_err=se_a, slope=b, slope_err=se_b, fit_r2=r2,
              t_slope=float(abs(b) / se_b) if se_b and np.isfinite(se_b) and se_b > 0 else np.nan,
              span_c=float(np.nanmax(c) - np.nanmin(c)),
              d_rg=float(np.nanmax(rg) - np.nanmin(rg)),
              rg_med=float(np.nanmedian(rg)), med_rg_err=float(np.nanmedian(rge)))
    if np.isfinite(mw).sum() >= 2:
        am, bm, seam, sebm, r2m = _linfit(c[np.isfinite(mw)], mw[np.isfinite(mw)])
        st.update(mw_at_c0=am, mw_slope=bm, mw_fit_r2=r2m,
                  mw_t=float(abs(bm) / sebm) if sebm and sebm > 0 else np.nan)
    if np.isfinite(i0).sum() >= 2:
        lr_i = _linfit(np.log10(c[np.isfinite(i0)]), np.log10(i0[np.isfinite(i0)]))
        st.update(i0_loglog_slope=lr_i[1], i0_loglog_r2=lr_i[4])

    t = st.get("t_slope", np.nan)
    med_err = st["med_rg_err"] if st["med_rg_err"] > 0 else np.nan
    significant = (np.isfinite(t) and t >= 2.0) and (med_err != med_err or st["d_rg"] > 2 * med_err)
    if not significant:
        st.update(verdict="🟢 无浓度依赖",
                  reason=f"Rg 在相对浓度 {c.min():.2f}–{c.max():.2f} 内变化 {st['d_rg']:.2f} Å（< 2×中位误差 "
                         f"{st['med_rg_err']:.2f}）或斜率不显著（|t|={t:.1f}）→ 单分散、无可见相互作用；"
                         f"c→0 外推 Rg={a:.2f}±{se_a:.2f} Å 与直接值 {st['rg_med']:.2f} Å 一致")
    elif b < 0:
        st.update(verdict="🟡 Rg 随浓度下降",
                  reason=f"斜率 {b:+.1f}±{se_b:.1f} Å/单位相对浓度（|t|={t:.1f}），"
                         f"c→0 外推 Rg={a:.2f}±{se_a:.2f} Å > 直接值 {st['rg_med']:.2f} Å → "
                         f"排斥型相互作用（或浓度依赖解离）；比较 c→0 值与其他方法")
    else:
        st.update(verdict="🟡 Rg 随浓度上升",
                  reason=f"斜率 {b:+.1f}±{se_b:.1f} Å/单位相对浓度（|t|={t:.1f}），"
                         f"c→0 外推 Rg={a:.2f}±{se_a:.2f} Å < 直接值 {st['rg_med']:.2f} Å → "
                         f"吸引/自缔合（低聚态）——再看 MW 是否同向")
    if st.get("mw_t", 0) >= 2.0 and st.get("mw_slope", 0) > 0:
        st["verdict"] = "🟡 MW 与 Rg 同向上升（自缔合/低聚）"
        st["reason"] += "；MW 也随 c 上升（|t|=%.1f），支持自缔合" % st["mw_t"]
    if st.get("i0_loglog_slope") is not None:
        s_i = st["i0_loglog_slope"]
        if not (0.8 <= s_i <= 1.2):
            st["reason"] += (f"；注意 I(0)–c 双对数斜率 {s_i:.2f}（理想 ≈1）——先怀疑**浓度代理**"
                             f"（色谱窗口/滚动中位基线对宽峰峰顶的压缩）与扣减线性，"
                             f"它只影响 c 轴标度（即 c→0 外推的截距），不影响 Rg 随 c 的**方向**判断")
    return st


def flank_consistency(rows, tol=0.15):
    """上升翼 vs 下降翼在相近浓度下的 Rg 是否一致（同一 c 的重复测量）。"""
    asc = [r for r in rows if r.get("side") == "asc" and r.get("rg", -1) > 0]
    desc = [r for r in rows if r.get("side") == "desc" and r.get("rg", -1) > 0]
    pairs = []
    for a in asc:
        d = min(desc, key=lambda r: abs(r["conc_rel"] - a["conc_rel"])) if desc else None
        if d and abs(d["conc_rel"] - a["conc_rel"]) <= 0.15:
            pairs.append((a, d, abs(a["rg"] - d["rg"]) / max(1e-9, 0.5 * (a["rg"] + d["rg"]))))
    if not pairs:
        return None
    worst = max(pairs, key=lambda t: t[2])
    return dict(n_pairs=len(pairs),
                worst_rel_diff=worst[2],
                worst=dict(asc=worst[0]["label"], desc=worst[1]["label"],
                           rg_asc=worst[0]["rg"], rg_desc=worst[1]["rg"]),
                ok=bool(worst[2] <= tol),
                msg=("🟢 两翼在相近浓度下一致（最大相对差 %.1f%%）" % (100 * worst[2])) if worst[2] <= tol
                    else ("🟡 两翼在相近浓度下不一致（最大相对差 %.1f%%：%s %.2f Å vs %s %.2f Å）→ "
                          "查基线漂移/组分变化/扣减残差" % (100 * worst[2], worst[0]["label"],
                                                       worst[0]["rg"], worst[1]["label"], worst[1]["rg"])))


# ----------------------------------------------------------------- 落盘
def write_gradient(out, rows, stats, flank, prefix):
    """CSV + JSON + 3 张图 + 一段给人看的结论文字。返回结论行（给 README 用）。"""
    gdir = os.path.join(out, "gradient")
    os.makedirs(gdir, exist_ok=True)
    with open(os.path.join(gdir, "gradient_slices.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_FIELDS)
        for r in rows:
            w.writerow([r.get(k, "") for k in CSV_FIELDS])
    with open(os.path.join(gdir, "gradient.json"), "w") as fh:
        json.dump(dict(slices=[{k: r.get(k) for k in CSV_FIELDS} for r in rows],
                       stats=stats, flank=flank), fh, indent=2, ensure_ascii=False, default=str)
    try:
        plot_gradient(os.path.join(gdir, "gradient.png"), rows, stats, prefix)
    except Exception as exc:                                    # 图失败不影响数据
        pass
    return summary_lines(rows, stats, flank)


def summary_lines(rows, stats, flank):
    lines = ["**峰内梯度（峰内不同位置当稀释序列）**", ""]
    lines.append(f"· 切片 {stats.get('n_used')}/{stats.get('n_total')} 个可用（同一 q 区间、都走 RAW 的 "
                 f"`set_sample_range` + `guinier_fit`）")
    if stats.get("verdict"):
        lines.append(f"· 判定：{stats['verdict']} —— {stats.get('reason', '')}")
    if stats.get("rg0") is not None and np.isfinite(stats.get("rg0", np.nan)):
        lines.append(f"· c→0 外推 Rg={stats['rg0']:.2f}±{stats.get('rg0_err', float('nan')):.2f} Å"
                     f"（直接中位 {stats.get('rg_med', float('nan')):.2f} Å）；"
                     f"斜率 {stats.get('slope', float('nan')):+.1f} Å/单位相对浓度，拟合 r²="
                     f"{stats.get('fit_r2', float('nan')):.3f}")
    if flank:
        lines.append(f"· 两翼互检：{flank['msg']}")
    return lines


def plot_gradient(path, rows, stats, prefix):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    use = [r for r in rows if "error" not in r and r.get("rg", -1) > 0]
    c = np.array([r["conc_rel"] for r in use], float)
    rg = np.array([r["rg"] for r in use], float)
    rge = np.array([abs(r.get("rg_err", 0.0)) for r in use], float)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))

    ax = axs[0]
    for r in rows:
        ax.plot(r["window"], [r["conc_rel"]] * 2, lw=2, alpha=0.7,
                color={"asc": "tab:blue", "desc": "tab:orange", "all": "black", "band": "tab:green"}
                .get(r.get("side"), "grey"))
    ax.scatter(c, [r["conc_rel"] for r in use], s=8, color="k", zorder=3)
    ax.set_xlabel("frame"); ax.set_ylabel("relative concentration (slice mean / peak apex)")
    ax.set_title(f"{prefix}: slices along the peak (blue=asc, orange=desc, black=whole)")
    ax.grid(alpha=0.3)

    ax = axs[1]
    if use:
        ax.errorbar(c, rg, yerr=rge, fmt="o", ms=4, capsize=2, color="tab:blue")
    if stats.get("rg0") is not None and np.isfinite(stats.get("rg0", np.nan)) and use:
        xs = np.array([0.0, float(np.nanmax(c))])
        ax.plot(xs, stats["rg0"] + stats.get("slope", 0.0) * xs, ls="--", color="tab:red",
                label=f"fit: Rg0={stats['rg0']:.1f}±{stats.get('rg0_err', float('nan')):.1f} A")
        ax.legend(fontsize=8)
    ax.set_xlabel("relative concentration"); ax.set_ylabel("Rg (A)")
    ax.set_title("Rg vs concentration (c->0 = intercept)"); ax.grid(alpha=0.3)

    ax = axs[2]
    mw = np.array([r.get("mw_vc", np.nan) for r in use], float)
    i0 = np.array([r.get("i0", np.nan) for r in use], float)
    if np.isfinite(mw).any():
        ax.plot(c[np.isfinite(mw)], mw[np.isfinite(mw)], "o-", ms=4, color="tab:green", label="MW (Vc)")
    ax.set_xlabel("relative concentration"); ax.set_ylabel("MW (kDa)")
    ax.set_title("MW vs concentration"); ax.grid(alpha=0.3)
    if np.isfinite(i0).any() and np.nanmax(i0) > 0:
        ax2 = ax.twinx()
        ax2.plot(c[np.isfinite(i0)], i0[np.isfinite(i0)], "s--", ms=3, color="0.5",
                 label="I(0) (right axis)")
        ax2.set_ylabel("I(0) [a.u.]", fontsize=8)
    h1, l1 = ax.get_legend_handles_labels()
    if h1:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
