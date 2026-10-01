#!/usr/bin/env python3
"""SEC 系列的**多峰**识别与逐峰区间划分（run-raw-sec-pipeline.py / find-sec-peaks.py 共用）。

为什么要有这个模块
------------------
RAW 自己的 ``find_sample_range`` / ``find_buffer_range`` 只认**最大的那个峰**：

    SASCalc.py:4071  findSampleRange():  max_peak_idx = np.argmax(peak_params['peak_heights'])

所以一条有两个洗脱峰的 SEC 系列，RAW（和只用 RAW 的流水线）只会把最大的那个峰当样品，
**第二个峰被静默丢掉**——产物看起来完全正常，只是少了东西。本模块补上"认峰 + 逐峰圈区间"。

口径（哪些是 RAW 的、哪些是我们自己的）
--------------------------------------
* **峰形判据沿用 RAW 内部口径**：平滑用 ``SASCalc.smooth_data``（savgol），找峰用
  ``SASCalc.find_peaks``（= scipy.signal.find_peaks 的包装），峰值/半高宽/峰底来自它返回的
  ``peak_heights / prominences / widths / left_ips / right_ips / left_bases / right_bases``。
* **只有两处是我们加的**：① 归一化前先减掉一条**滚动中位基线**（原样除以全局最大值时，
  10% 量级的基线漂移会让整条曲线都贴着 1.0，噪声峰满天飞）；② 峰区间的**划分**
  （RAW 不给"哪个峰该配哪段 buffer"）。
* 不下任何拟合结论：本模块只回答"有几个峰、每个峰占哪几帧、扣减用哪几帧"。

色谱图用什么
------------
默认用**扣减后曲线在低 q 窗口 [0.01, 0.05] 1/Å 的积分强度**（``SASM.getIofQRange``，RAW 自己的
梯形积分）。理由：低 q 区对溶质最敏感、对噪声最不敏感，比总积分强度干净得多；而官方
``plot_series`` 的"某个 q 区间的强度"也是同一个思路。**总强度（``getIntI``）不要用来找峰**——
本机实测 BSA 2000 帧：总强度只有 ~10% 起伏且被束位漂移主导，任何阈值都会找出十几个假峰。

视觉检验
--------
自动判定一定会有边界情形（肩峰、聚集体前峰、气泡/快门野值帧）。``render_overview`` /
``render_zoom`` / ``render_sweep`` 三个出图函数就是留给"代码判不了就上眼睛"的那条路：
把检测结果画成图（峰位置、区间、buffer 段、逐帧 Rg/I0/MW 平台），人（或 agent 看图）确认后，
再用 ``--peak-ranges`` 手工覆盖。判据与操作序列见 SKILL.md 的「多峰 SEC 数据」一节。
"""
import csv
import json
import math
import os

import numpy as np

try:                                    # 只有真跑 RAW 时才需要
    import bioxtasraw.SASCalc as _sascalc
except Exception:                       # pragma: no cover - 纯区间运算不依赖 RAW
    _sascalc = None

# 默认参数（单位都标在名字里；CLI 侧保持一致，别只改一处）
DEF_Q_RANGE = (0.01, 0.05)      # 1/A，色谱图取这个 q 窗口
DEF_MIN_HEIGHT = 0.05           # 相对（基线校正后）最高峰的峰高下限
DEF_MIN_PROMINENCE = 0.2        # 相对最高峰的 prominence 下限（见下面注释：默认偏保守）
DEF_MIN_WIDTH = 5               # frame，半高宽下限（低于它的算尖刺/噪声）
DEF_MIN_SNR = 5.0               # 峰高 / 基线噪声（1.4826*MAD）下限
DEF_MAX_PEAKS = 6               # 最多分析几个峰（多的按 prominence 排掉并告警）
DEF_BUFFER_MIN_LEN = 10         # frame，单侧 buffer 段的帧数下限
# 为什么 prominence 默认 0.2 而不是 0.05：**扣减后的色谱图本身就带系统性纹波**
# （单个全局 buffer 平均 ± 束位漂移 → 剩余起伏可达主峰的 10–20%），本机实测 4EH2 那条
# 干净的单体/二聚体数据里，纹波在三个 q 窗口都能复现、峰值达 7–15σ、宽度 10–18 帧，
# 任何"相对噪声"或"多 q 窗口一致性"判据都杀不掉它们——只有"相对主峰的幅度"能。
# 想找更小的组分（<20% 主峰）：调低这个阈值，然后**必须**看 sec_peaks.png 逐个人工确认，
# 或改用 UV 痕/裁剪视频对照。见 SKILL.md「多峰 SEC 数据」一节。

PEAK_COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown',
               'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']


# ----------------------------------------------------------------- 小工具
def parse_ranges(text):
    """'150,320;460,520' → [[150,320],[460,520]]（0 基闭区间）。"""
    out = []
    for part in str(text).split(';'):
        part = part.strip()
        if not part:
            continue
        a, b = (int(x) for x in part.split(','))
        out.append([a, b])
    return out


def fmt_ranges(ranges):
    return ';'.join(f"{a},{b}" for a, b in ranges) if ranges else ''


def odd_at_most(value, upper):
    """取 ≤ upper 的最大奇数（≥1）。"""
    v = int(value)
    if v > upper:
        v = upper
    v = max(1, v)
    if v % 2 == 0:
        v -= 1
    return max(1, v)


def rolling_median(y, window):
    """滚动中位数（边缘用 edge padding，长度不变）。O(n·w)，n 是帧数很小，够用。"""
    y = np.asarray(y, dtype=float)
    n = len(y)
    w = odd_at_most(window, n if n % 2 == 1 else n - 1)
    if w <= 1:
        return y.copy()
    half = w // 2
    pad = np.pad(y, half, mode='edge')
    idx = np.lib.stride_tricks.sliding_window_view(pad, w)
    return np.median(idx, axis=1)


def smooth(y, window=None, order=5):
    """平滑：优先用 RAW 的 SASCalc.smooth_data（savgol），没有就本地 savgol。"""
    y = np.asarray(y, dtype=float)
    n = len(y)
    wl = odd_at_most(window if window else min(51, n // 2), n)
    if wl < 3:
        return y.copy()
    order = min(order, wl - 1)
    if _sascalc is not None:
        return _sascalc.smooth_data(y, window_length=wl, order=order)
    from scipy.signal import savgol_filter
    return savgol_filter(y, wl, order)


def find_peaks(y_norm, height, width, min_prominence):
    """找峰：优先 RAW 的 SASCalc.find_peaks（scipy 包装），再按 prominence 过滤。"""
    y_norm = np.asarray(y_norm, dtype=float)
    if _sascalc is not None:
        peaks, props = _sascalc.find_peaks(y_norm, height=height, width=width, rel_height=0.5)
    else:
        from scipy.signal import find_peaks as _fp
        peaks, props = _fp(y_norm, height=height, width=width, rel_height=0.5)
    keep = [i for i, p in enumerate(peaks) if props['prominences'][i] >= min_prominence]
    out = {}
    for key in ('peak_heights', 'prominences', 'widths', 'left_ips', 'right_ips',
                'left_bases', 'right_bases'):
        out[key] = np.asarray(props[key], dtype=float)[keep]
    return np.asarray(peaks, dtype=int)[keep], out


def edf_frames(sasms, q1, q2, label='q'):
    """逐帧积分强度（扣减后曲线的低 q 窗口）→ 色谱图。"""
    if not sasms:
        return np.zeros(0)
    q = np.asarray(sasms[0].getQ(), dtype=float)
    lo, hi = max(q1, q[0]), min(q2, q[-1])
    if not lo < hi:
        raise ValueError(f"色谱图的 q 窗口 {q1}–{q2} 与数据 q 范围 {q[0]:.4f}–{q[-1]:.4f} 不相交")
    return np.array([float(s.getIofQRange(lo, hi)) for s in sasms], dtype=float)


def no_beam_mask(int_unsub, frac=0.5):
    """无束流帧（快门没开/断束）：总强度低于全系列中位数 frac 倍的帧。

    本机实测 4DH1-KDPV-ZN 存在这种帧（总强度掉到 ~0），扣减后是一根巨大的负尖刺，
    会被 find_peaks 当峰。这类帧必须从色谱图判断与 buffer 里都剔掉。
    """
    v = np.asarray(int_unsub, dtype=float)
    if len(v) == 0:
        return np.zeros(0, dtype=bool)
    med = np.nanmedian(v)
    if not np.isfinite(med) or med <= 0:
        return np.zeros(len(v), dtype=bool)
    return ~np.isfinite(v) | (v < frac * med)


# ----------------------------------------------------------------- 主检测
def detect_peaks(y_raw, q_range=DEF_Q_RANGE, smooth_window=None, baseline_window=None,
                 min_height=DEF_MIN_HEIGHT, min_prominence=DEF_MIN_PROMINENCE,
                 min_width=DEF_MIN_WIDTH, min_snr=DEF_MIN_SNR, max_peaks=DEF_MAX_PEAKS,
                 mask=None, window_mode='half-height', buffer_mode='local',
                 buffer_min_len=DEF_BUFFER_MIN_LEN, buffer_exclude=None, buffer_max_len=0):
    """在色谱图 ``y_raw`` 上认峰。

    返回 dict：原始/平滑/基线/基线校正后的曲线、参数、以及 peaks 列表
    （每项：apex / height / prominence / width / half-height 边界 / 峰底边界 / sample 窗）。
    """
    y_raw = np.asarray(y_raw, dtype=float)
    n = len(y_raw)
    if n < 5:
        raise ValueError(f"帧数太少（{n}），做不了峰检测")
    y_use = y_raw.copy()
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.any():
            y_use[mask] = np.nan                      # 野值帧不参与平滑/基线
    y_fill = y_use.copy()
    if np.isnan(y_fill).any():                         # 线性插值补齐，仅用于平滑/基线
        good = ~np.isnan(y_fill)
        y_fill = np.interp(np.arange(n), np.flatnonzero(good), y_fill[good])

    y_sm = smooth(y_fill, smooth_window)
    wb = baseline_window or odd_at_most(max(21, min(401, n // 10)), n)
    base = rolling_median(y_sm, wb)
    y_corr = y_sm - base
    if mask is not None and mask.any():
        y_corr = y_corr.copy()
        y_corr[mask] = 0.0                             # 别让野值帧进峰判据
    top = float(np.max(y_corr)) if len(y_corr) else 0.0
    if not np.isfinite(top) or top <= 0:
        return dict(ok=False, reason="基线校正后整条色谱图没有正值（没有可识别的洗脱峰）",
                    y_raw=y_raw, y_smooth=y_sm, base=base, y_corr=y_corr,
                    y_norm=np.zeros(n), params=dict(), peaks=[])
    y_norm = np.clip(y_corr / top, 0.0, None)

    peaks, props = find_peaks(y_norm, min_height, min_width, min_prominence)

    # ---- 信噪比过滤：峰高 / 基线噪声（1.4826*MAD，用"不在候选峰内"的帧估）
    cand_mask = np.zeros(n, dtype=bool)
    for j in range(len(peaks)):
        lo = int(round(props['left_ips'][j]))
        hi = int(round(props['right_ips'][j]))
        cand_mask[max(0, lo):min(n - 1, hi) + 1] = True
    quiet = y_corr[~cand_mask]
    quiet = quiet[np.isfinite(quiet)]
    sigma = float(1.4826 * np.median(np.abs(quiet - np.median(quiet)))) if len(quiet) > 10 else np.nan
    keep = []
    for i in range(len(peaks)):
        h_abs = float(props['peak_heights'][i]) * top
        snr = h_abs / sigma if sigma and sigma > 0 else float('inf')
        props.setdefault('snrs', np.zeros(len(peaks)))
        props['snrs'][i] = snr
        if np.isfinite(sigma) and sigma > 0 and snr < min_snr:
            continue
        keep.append(i)
    keep = np.asarray(keep, dtype=int)

    # 按 prominence 排序取前 max_peaks 个，再按帧号排回去（保证输出顺序 = 洗脱顺序）
    order = keep[np.argsort(props['prominences'][keep])[::-1]] if len(keep) else keep
    dropped = []
    if len(order) > max_peaks:
        dropped = [int(peaks[i]) for i in order[max_peaks:]]
        order = order[:max_peaks]
    order = np.sort(order)

    plist = []
    for k, i in enumerate(order):
        apex = int(peaks[i])
        plist.append(dict(index=k, apex=apex,
                          height=float(props['peak_heights'][i]),
                          prominence=float(props['prominences'][i]),
                          width=float(props['widths'][i]),
                          snr=float(props['snrs'][i]) if len(props.get('snrs', [])) else float('nan'),
                          half_height_abs=float(props['peak_heights'][i]) * top,
                          left_ips=float(props['left_ips'][i]),
                          right_ips=float(props['right_ips'][i]),
                          left_base=int(props['left_bases'][i]),
                          right_base=int(props['right_bases'][i]),
                          on_no_beam_frame=bool(mask[apex]) if mask is not None and len(mask) == n else False,
                          near_no_beam_frame=bool(np.any(mask[max(0, apex - 3):min(n, apex + 4)]))
                          if mask is not None and len(mask) == n else False))
    resolve_windows(plist, n, y_corr, window_mode)
    buffer_ranges_for_peaks(plist, n, mode=buffer_mode, min_len=buffer_min_len,
                            exclude=buffer_exclude if buffer_exclude is not None else mask,
                            max_len=buffer_max_len)
    return dict(ok=bool(plist), reason='' if plist else
                "没有峰同时满足 峰高/prominence/半高宽/信噪比 四个下限（数据可能没有洗脱峰，"
                "或阈值偏严 → 用 find-sec-peaks.py --sweep 扫一遍灵敏度）",
                q_range=list(q_range), smooth_window=int(odd_at_most(smooth_window or min(51, n // 2), n)),
                baseline_window=int(wb), noise_sigma=sigma, top=top,
                params=dict(min_height=min_height, min_prominence=min_prominence,
                            min_width=min_width, min_snr=min_snr, max_peaks=max_peaks),
                dropped_apexes=dropped, window_mode=window_mode, buffer_mode=buffer_mode,
                mask=mask, y_raw=y_raw, y_smooth=y_sm, base=base, y_corr=y_corr,
                y_norm=y_norm, peaks=plist)


def resolve_windows(peaks, n, y_corr, mode='half-height'):
    """给每个峰定 sample 窗（0 基闭区间）。

    * ``half-height``：半高宽的两个交点（find_peaks 的 left_ips/right_ips）——峰中心，最保守；
    * ``valley``：到相邻峰的谷底（两峰顶点之间的最小值）——"整个峰都算"。
    """
    for k, p in enumerate(peaks):
        if mode == 'valley':
            left = int(math.floor(p['left_base']))
            right = int(math.ceil(p['right_base']))
            if k > 0:
                seg = y_corr[peaks[k - 1]['apex']:p['apex'] + 1]
                if len(seg):
                    left = peaks[k - 1]['apex'] + int(np.argmin(seg))
            if k < len(peaks) - 1:
                seg = y_corr[p['apex']:peaks[k + 1]['apex'] + 1]
                if len(seg):
                    right = p['apex'] + int(np.argmin(seg))
        else:
            left = int(round(p['left_ips']))
            right = int(round(p['right_ips']))
        left = max(0, min(left, n - 1))
        right = max(0, min(right, n - 1))
        if right < left:
            left, right = right, left
        if right - left < 2:                       # 至少 3 帧，否则给不出曲线
            left = max(0, p['apex'] - 1)
            right = min(n - 1, p['apex'] + 1)
        p['window'] = [int(left), int(right)]


def peak_windows(peaks):
    return [list(p['window']) for p in peaks]


def complement_blocks(n, windows):
    """帧 0..n-1 里不在任何峰窗内的连续段。"""
    mask = np.zeros(n, dtype=bool)
    for a, b in windows:
        mask[max(0, a):min(n - 1, b) + 1] = True
    blocks = []
    start = None
    for i in range(n):
        if not mask[i] and start is None:
            start = i
        elif mask[i] and start is not None:
            blocks.append([start, i - 1])
            start = None
    if start is not None:
        blocks.append([start, n - 1])
    return blocks


def buffer_ranges_for_peaks(peaks, n, mode='local', min_len=DEF_BUFFER_MIN_LEN, exclude=None,
                            max_len=0):
    """给每个峰挑扣减用的 buffer 段（0 基闭区间）。

    * ``local``（默认）：每个峰用它**紧邻的**前一段 + 后一段非峰区（就是"峰前+峰后"两段，
      这是官方教程与源A都推荐的双缓冲液区做法）；
    * ``global``：所有峰共用全部非峰区（等价于 RAW 的整条系列 buffer）。

    ``max_len`` 给了就把每段**截到紧邻峰的那一侧**的 max_len 帧内（0 = 不截）：
    基线在漂移时（实测 BSA 峰前 −0.10 / 峰后 +0.02），水平匹配比多平均更重要——
    段拉得越长，扣减后残留的斜漂移越大。

    同时给出每段的位置类别，便于人复核：
    ``edge-before``（第一个峰之前）/ ``between-peaks``（两峰之间，**可能含混合组分**）/ ``edge-after``。
    """
    windows = peak_windows(peaks)
    blocks = complement_blocks(n, windows)
    exclude = np.asarray(exclude, dtype=bool) if exclude is not None else None

    def usable(block):
        seg = block
        if exclude is not None:
            idx = [i for i in range(block[0], block[1] + 1) if not exclude[i]]
            if not idx:
                return None
            return [idx[0], idx[-1]]
        return seg

    def kind(block):
        if not peaks:
            return 'edge'
        if block[1] < peaks[0]['window'][0]:
            return 'edge-before'
        if block[0] > peaks[-1]['window'][1]:
            return 'edge-after'
        return 'between-peaks'

    for k, p in enumerate(peaks):
        if mode == 'global':
            chosen = [[b, kind(b)] for b in blocks]
        else:
            lo, hi = p['window']
            before = max((b for b in blocks if b[1] < lo), key=lambda b: b[1], default=None)
            after = min((b for b in blocks if b[0] > hi), key=lambda b: b[0], default=None)
            chosen = [[b, kind(b)] for b in (before, after) if b is not None]
        routed = []
        for b, kd in chosen:
            b2 = usable(b)
            if b2 is None:
                continue
            if max_len and b2[1] - b2[0] + 1 > max_len:
                if b2[1] < p['window'][0]:            # 峰前那段 → 留靠峰的右端
                    b2 = [b2[1] - max_len + 1, b2[1]]
                elif b2[0] > p['window'][1]:          # 峰后那段 → 留靠峰的左端
                    b2 = [b2[0], b2[0] + max_len - 1]
            if b2[1] - b2[0] + 1 < min_len:
                p.setdefault('buffer_warnings', []).append(
                    f"buffer 段 {b2} 只有 {b2[1]-b2[0]+1} 帧（< {min_len}）→ 丢弃该段")
                continue
            routed.append([int(b2[0]), int(b2[1]), kd])
        if not routed:
            p.setdefault('buffer_warnings', []).append(
                "没有可用的 buffer 段 → 回退到全局非峰区")
            routed = [[int(usable(b)[0]), int(usable(b)[1]), kind(b)]
                      for b in blocks if usable(b) is not None]
        p['buffers'] = routed
        p['buffer_ranges'] = [[a, b] for a, b, _ in routed]
    return peaks


def resolution(pk_a, pk_b, y_corr, top=None, smooth_window=15):
    """相邻两峰的分辨率判据：**谷底深度 / 较矮那个峰的高**。

    0 → 完全回到基线（真分开，各自可以独立扣减平均）；
    ≥0.5 → 部分重叠；≥0.9 → 基本没分开（这种不是"多峰"，是"未解析重叠"，
    该去做 SVD/EFA/REGALS 分解，不是这里切区间）。

    两个细节别省：① 峰高是"相对最高峰"的归一化值，谷底是绝对强度 → 必须一起除以
    ``top``（= 基线校正后的最大值）才能比；② 谷底要在**平滑过**的曲线上取，否则噪声
    让谷底随便跌到 0 以下，任何两峰都会被判成"已分开"。
    """
    lo, hi = pk_a['apex'], pk_b['apex']
    if hi <= lo:
        return 1.0
    y = np.asarray(y_corr, dtype=float)
    top = float(np.nanmax(y)) if top is None else float(top)
    if top <= 0:
        return 1.0
    v = smooth(y, smooth_window) if smooth_window else y
    valley = float(np.min(v[lo:hi + 1]))
    short = min(pk_a['height'], pk_b['height'])
    if short <= 0:
        return 1.0
    return max(0.0, valley / top) / short


def adjacent_valley_ratio(det, k):
    """det 里第 k 与第 k+1 个峰之间的"谷底/峰高"（归一化，判两峰分不分得开）。"""
    peaks = det['peaks']
    return resolution(peaks[k], peaks[k + 1], det['y_corr'], top=det.get('top'))


def manual_peaks(y_raw, ranges, mask=None, buffer_mode='local',
                 buffer_min_len=DEF_BUFFER_MIN_LEN, buffer_exclude=None, note='',
                 buffer_max_len=0):
    """用**人工给的区间**当峰（0 基闭区间）：曲线照旧算出来做参照，峰位/峰窗按给的来。

    这就是"代码判不了 → 眼睛判 → 用眼睛的结论覆盖"的那条路（见 SKILL.md 多峰章节）：
    ``--peak-ranges`` / ``find-sec-peaks.py --peak-ranges`` 都走这里。
    """
    det = detect_peaks(y_raw, mask=mask, min_height=0.0, min_prominence=0.0, min_snr=0.0,
                       min_width=1, max_peaks=1000, buffer_mode=buffer_mode,
                       buffer_min_len=buffer_min_len, buffer_exclude=buffer_exclude)
    y_corr, top, sigma, n = det['y_corr'], det['top'], det['noise_sigma'], len(det['y_corr'])
    plist = []
    for k, (lo, hi) in enumerate(ranges):
        lo = max(0, min(int(lo), n - 1))
        hi = max(lo, min(int(hi), n - 1))
        seg = y_corr[lo:hi + 1]
        apex = lo + int(np.argmax(seg))
        h_abs = float(y_corr[apex])
        h = h_abs / top if top > 0 else 0.0
        above = np.flatnonzero(seg >= 0.5 * h_abs)
        width = float(above[-1] - above[0] + 1) if len(above) else float(hi - lo + 1)
        plist.append(dict(index=k, apex=int(apex), height=float(h),
                          prominence=float(max(0.0, (h_abs - float(np.min(seg))) / top)) if top > 0 else 0.0,
                          width=width, snr=float(h_abs / sigma) if sigma and sigma > 0 else float('inf'),
                          half_height_abs=h_abs, left_ips=float(lo), right_ips=float(hi),
                          left_base=int(lo), right_base=int(hi), manual=True,
                          on_no_beam_frame=bool(mask[apex]) if mask is not None and len(mask) == n else False,
                          near_no_beam_frame=bool(np.any(mask[max(0, apex - 3):min(n, apex + 4)]))
                          if mask is not None and len(mask) == n else False,
                          window=[int(lo), int(hi)]))
    det['peaks'] = plist
    det['manual_ranges'] = [[int(a), int(b)] for a, b in ranges]
    det['ok'] = bool(plist)
    det['reason'] = '' if plist else '人工区间是空的'
    if note:
        det['manual_note'] = note
    buffer_ranges_for_peaks(plist, n, mode=buffer_mode, min_len=buffer_min_len,
                            exclude=buffer_exclude if buffer_exclude is not None else mask,
                            max_len=buffer_max_len)
    return det


def apply_manual_buffers(peaks, spec, n):
    """人工指定每个峰的 buffer：'s,e;s,e|s,e' —— ``|`` 分隔峰、``;`` 分隔段。

    给得比峰少的那些峰保持原来的 buffer；每段标 kind='manual'。
    """
    groups = [g for g in str(spec).split('|')]
    for k, p in enumerate(peaks):
        if k >= len(groups) or not groups[k].strip():
            continue
        rngs = [[int(x) for x in part.split(',')] for part in groups[k].split(';') if part.strip()]
        rngs = [[max(0, a), min(n - 1, b)] for a, b in rngs]
        p['buffers'] = [[a, b, 'manual'] for a, b in rngs]
        p['buffer_ranges'] = rngs
    return peaks


def resolution_R(pk_a, pk_b):
    """经典色谱分辨率 R = 2Δt/(w1+w2)（w = 半高宽 FWHM）。

    判据（色谱学惯例）：R ≥ 1.5 基线分离；1.0–1.5 部分重叠；< 1.0 基本没分开。
    和谷底判据一起看：两个判据一致才叫"确实分开了"。
    """
    w = 0.5 * (float(pk_a.get('width', 0)) + float(pk_b.get('width', 0)))
    if w <= 0:
        return float('nan')
    return abs(pk_b['apex'] - pk_a['apex']) / w


def label_R(R):
    if not np.isfinite(R):
        return '—'
    if R >= 1.5:
        return '🟢 基线分离'
    if R >= 1.0:
        return '🟡 部分重叠'
    return '🔴 基本没分开'


def buffer_level_mismatch(peaks, y_corr, top=None):
    """两侧 buffer 段的**水平差**（相对最高峰）——SEC 最常见的陷阱之一。

    峰前与峰后基线常不相等（束位漂移/柱温/损伤蛋白粘窗）。实测 BSA 2000 帧：峰前段在
    −0.10、峰后段在 +0.02，差到主峰的 ~60%（那次单靠一段 buffer 扣减，Rg 在 38–117 Å 乱跳）。
    返回 (差/峰高, 水平低的那段, 水平高的那段, 段号说明)；只有一个 buffer 段时返回 (0,...)。
    """
    top = float(np.nanmax(y_corr)) if top is None else float(top)
    worst = (0.0, None, None, '')
    if top <= 0:
        return worst
    for p in peaks:
        bufs = p.get('buffers', [])
        if len(bufs) < 2:
            continue
        levels = [(float(np.median(np.asarray(y_corr)[a:b + 1])), (a, b), kd) for a, b, kd in bufs]
        levels.sort()
        diff = abs(levels[-1][0] - levels[0][0]) / top
        if diff > worst[0]:
            worst = (diff, levels[0][1], levels[-1][1], f"P{p['index']+1}: {levels[0][2]} vs {levels[-1][2]}")
    return worst


def label_verdict(ratio):
    if ratio >= 0.9:
        return '🔴 基本没分开'
    if ratio >= 0.5:
        return '🟡 部分重叠'
    if ratio >= 0.2:
        return '🟢 基线附近分开'
    return '🟢 基线分离'


def label_verdict_en(ratio):
    """同一判据的英文版（图里用英文标签，控制台/README 用中文）。"""
    if ratio >= 0.9:
        return 'NOT resolved'
    if ratio >= 0.5:
        return 'partially overlapped'
    if ratio >= 0.2:
        return 'near baseline'
    return 'baseline separated'


def sweep_combos(base_prominence=DEF_MIN_PROMINENCE):
    """灵敏度扫描的网格：只看"幅度阈值 × 信噪比阈值"两个轴（变化最敏感的两项）。

    给"看图挑参数"用：同一张图里画几组阈值的检测结果，人眼挑与峰形相符的那一组。
    """
    out = []
    for prom in sorted({0.4, 0.25, round(base_prominence, 3), 0.12, 0.06, 0.03}, reverse=True):
        snr = 10 if prom >= 0.2 else (5 if prom >= 0.06 else 3)
        out.append((f"prom>={prom}, snr>={snr}", dict(min_prominence=prom, min_snr=snr)))
    return out


# ----------------------------------------------------------------- 出图
def _style_ax(ax, xlabel=False, ylabel=None, logy=False):
    ax.grid(alpha=0.3)
    if logy:
        ax.set_yscale('log')
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9)
    if xlabel:
        ax.set_xlabel('frame')


def render_overview(path, det, series_arrays=None, title=None, axs=None):
    """总览图：色谱图（峰/区间/buffer 都标出来）+ 逐帧 Rg / I(0) / MW。

    未扣减曲线（漂移大、量级大）会把基线校正后的曲线压扁 → 两者量级差 >5× 时
    给校正曲线开一条右轴（图例里写明 right axis），别把峰看没了。

    ``axs`` 给外部复用时是一个 4 个 axes 的数组（灵敏度扫描图用）。
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    own = axs is None
    if own:
        fig, axs = plt.subplots(4, 1, figsize=(13, 12), sharex=True)
    else:
        fig = axs[0].figure
    y_raw, y_sm, base, y_corr = det['y_raw'], det['y_smooth'], det['base'], det['y_corr']
    n = len(y_raw)
    x = np.arange(n)

    ax = axs[0]
    ax.plot(x, y_raw, lw=0.6, color='0.75', label='chromatogram (raw)')
    ax.plot(x, y_sm, lw=1.3, color='tab:blue', label='smoothed')
    ax.plot(x, base, lw=1.0, ls=':', color='tab:red', label='rolling-median baseline')
    raw_span = float(np.nanmax(y_raw) - np.nanmin(y_raw)) if n else 0.0
    corr_span = float(np.nanmax(y_corr) - np.nanmin(y_corr)) if n else 0.0
    twin = bool(corr_span > 0 and raw_span > 5 * corr_span)
    axc = ax.twinx() if twin else ax
    axc.plot(x, y_corr, lw=1.1, color='tab:green', alpha=0.9,
             label='baseline-corrected' + (' (right axis)' if twin else ''))
    for p in det['peaks']:
        c = PEAK_COLORS[p['index'] % len(PEAK_COLORS)]
        lo, hi = p['window']
        ax.axvspan(lo, hi, color=c, alpha=0.13)
        ax.axvline(p['apex'], color=c, lw=1.1, ls='-')
        axc.annotate(f"P{p['index']+1}", (p['apex'], y_corr[p['apex']]),
                     textcoords='offset points', xytext=(3, 8), fontsize=10, color=c, weight='bold')
        for a, b, kd in p.get('buffers', []):
            ax.axvspan(a, b, color=c, alpha=0.07, hatch='//', edgecolor=c, lw=0)
    if det.get('mask') is not None and np.asarray(det['mask']).any():
        for i in np.flatnonzero(det['mask']):
            ax.axvline(i, color='black', lw=1.2, alpha=0.5)
    ax.set_ylabel('I(q {}–{})  [a.u.]'.format(*det.get('q_range', DEF_Q_RANGE))
                  if not twin else 'raw / smoothed  [a.u.]', fontsize=9)
    if twin:
        axc.set_ylabel('baseline-corrected  [a.u.]', fontsize=9)
        axc.tick_params(labelsize=8)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = axc.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, ncol=2, loc='upper left')
    _style_ax(ax)
    ax.set_title(title or (f"SEC chromatogram + detected peaks ({len(det['peaks'])} peak(s); "
                           f"shaded = peak sample window, hatched = that peak's buffer segments, "
                           f"black lines = suspected no-beam frames)"))

    labels = ['Rg (A)', 'I(0)', 'MW (kDa)']
    keys = ['rg', 'i0', 'mw']
    for k in range(3):
        ax = axs[k + 1]
        arrs = (series_arrays or {}).get(keys[k])
        # 只有"真的有有限值"才画：全 NaN（逐峰参数没算出来 / RAW 全是 -1 哨兵）时画出来是空轴，
        # 读者会把它当成"值≈0"。宁可明说没有，也不要留一个看起来像 0 的空面板。
        drew = False
        if arrs:
            for lab, arr in arrs.items():
                a = np.asarray(arr, float)
                if np.isfinite(a).any():
                    ax.plot(x, a, lw=0.8, label=lab)
                    drew = True
            if drew:
                ax.legend(fontsize=8)
        if not drew:
            ax.text(0.02, 0.5, 'no per-frame values (RAW computed none for this peak: '
                               'all -1 / set_buffer_range not run)',
                    fontsize=8, transform=ax.transAxes, color='0.5')
        for p in det['peaks']:
            c = PEAK_COLORS[p['index'] % len(PEAK_COLORS)]
            ax.axvspan(p['window'][0], p['window'][1], color=c, alpha=0.13)
        ax.set_ylabel(labels[k], fontsize=9)
        if k == 2 and drew:
            ax.set_yscale('log')
        _style_ax(ax, xlabel=(k == 2))

    if own:
        fig.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path, dpi=130)
        plt.close(fig)
    return fig, axs


def render_zoom(path, det, pad=0.6):
    """逐峰放大图：每行一个峰，看"这个峰到没到基线、两侧 buffer 干不干净"。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    peaks = det['peaks']
    if not peaks:
        return None
    y_raw, y_sm, y_corr, base = det['y_raw'], det['y_smooth'], det['y_corr'], det['base']
    n = len(y_raw)
    fig, axs = plt.subplots(len(peaks), 1, figsize=(12, 2.8 * len(peaks)), squeeze=False)
    for k, p in enumerate(peaks):
        ax = axs[k][0]
        lo, hi = p['window']
        xspan = [lo, hi] + [x for a, b, _ in p.get('buffers', []) for x in (a, b)] + [p['apex']]
        w = max(5, hi - lo)
        a = max(0, int(min(xspan) - pad * w))
        b = min(n - 1, int(max(xspan) + pad * w))
        x = np.arange(a, b + 1)
        ax.plot(x, y_raw[a:b + 1], lw=0.6, color='0.75', label='raw')
        ax.plot(x, y_sm[a:b + 1], lw=1.3, color='tab:blue', label='smoothed')
        ax.plot(x, base[a:b + 1], lw=1.0, ls=':', color='tab:red', label='rolling-median baseline')
        ax.axhline(0, color='0.4', lw=0.6)
        c = PEAK_COLORS[p['index'] % len(PEAK_COLORS)]
        ax.axvspan(lo, hi, color=c, alpha=0.15)
        for aa, bb, kd in p.get('buffers', []):
            ax.axvspan(aa, bb, color='k', alpha=0.07)
            ax.annotate(kd, (aa + (bb - aa) / 2, ax.get_ylim()[0]), fontsize=7, color='0.35',
                        ha='center', va='bottom')
        ratio = adjacent_valley_ratio(det, p['index']) if p['index'] < len(peaks) - 1 else None
        ax.set_title(f"P{p['index']+1}: sample frames {lo}-{hi} ({hi-lo+1} frames, apex {p['apex']})"
                     f" | buffer {fmt_ranges([[aa,bb] for aa,bb,_ in p.get('buffers',[])])}"
                     f" | prominence {p['prominence']:.2f}x max, {p['snr']:.0f} sigma, "
                     f"FWHM {p['width']:.1f} frames"
                     + (f" | valley/height to next P{p['index']+2} = {ratio:.2f} {label_verdict_en(ratio)}"
                        if ratio is not None else ""),
                     fontsize=9, loc='left')
        _style_ax(ax, ylabel='I sub (a.u.)', xlabel=(k == len(peaks) - 1))
        if k == 0:
            ax.legend(fontsize=8)
        ax.set_xlim(a, b)          # 别让画到画面外的 buffer 阴影把 x 轴拉开
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def render_sweep(path, y_raw, combos, series_arrays=None, mask=None, q_range=DEF_Q_RANGE):
    """灵敏度扫描图：同一张图里画 grid 组参数的检测结果，供人眼挑一组。

    ``combos`` 是 (标签, kwargs) 列表；每个格子画色谱图 + 检测到的峰位（竖线）。
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ncol = min(3, len(combos))
    nrow = int(math.ceil(len(combos) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 2.6 * nrow), squeeze=False, sharex=True)
    x = np.arange(len(y_raw))
    for i, (label, kw) in enumerate(combos):
        ax = axs[i // ncol][i % ncol]
        try:
            det = detect_peaks(y_raw, q_range=q_range, mask=mask, **kw)
        except Exception as exc:                                   # pragma: no cover
            ax.text(0.5, 0.5, f"{label}\n失败：{exc}", fontsize=8, ha='center')
            continue
        ax.plot(x, det['y_corr'], lw=1.0, color='tab:green')
        nwin = len(det['peaks'])
        for p in det['peaks']:
            c = PEAK_COLORS[p['index'] % len(PEAK_COLORS)]
            ax.axvline(p['apex'], color=c, lw=1.0)
            ax.axvspan(p['window'][0], p['window'][1], color=c, alpha=0.13)
            ax.annotate(f"P{p['index']+1}", (p['apex'], det['y_corr'][p['apex']]), fontsize=8,
                        color=c, textcoords='offset points', xytext=(2, 4))
        ax.set_title(f"{label} -> {nwin} peak(s)", fontsize=9, loc='left')
        _style_ax(ax, ylabel='I_corr' if i % ncol == 0 else None,
                  xlabel=(i >= len(combos) - ncol))
    for j in range(len(combos), nrow * ncol):
        axs[j // ncol][j % ncol].axis('off')
    fig.suptitle("peak-detection sensitivity sweep (pick the row that matches the peaks you see "
                 "by eye, then pass the same thresholds as --peak-min-*)", fontsize=10)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


# ----------------------------------------------------------------- 落盘
def peaks_rows(det):
    """CSV/表用的行（每峰一行 + 相邻峰的谷底分辨与经典 R）。"""
    rows = []
    peaks = det['peaks']
    y_corr = det['y_corr']
    for k, p in enumerate(peaks):
        ratio = None
        if k < len(peaks) - 1:
            ratio = resolution(p, peaks[k + 1], y_corr, top=det.get('top'))
        keep = {kk: vv for kk, vv in p.items()
                if kk not in ('buffers', 'buffer_warnings')}
        keep['buffer_ranges'] = fmt_ranges(p.get('buffer_ranges', []))
        keep['buffer_kinds'] = ';'.join(kd for _, _, kd in p.get('buffers', []))
        keep['n_sample_frames'] = p['window'][1] - p['window'][0] + 1
        keep['valley_ratio_to_next'] = '' if ratio is None else round(ratio, 4)
        keep['resolution_to_next'] = '' if ratio is None else label_verdict(ratio)
        keep['chrom_R_to_next'] = ('%.2f' % resolution_R(p, peaks[k + 1])) if ratio is not None else ''
        keep['warnings'] = ' / '.join(p.get('buffer_warnings', []))
        rows.append(keep)
    return rows


CSV_FIELDS = ['index', 'apex', 'height', 'prominence', 'snr', 'width', 'half_height_abs',
              'left_ips', 'right_ips', 'left_base', 'right_base', 'window', 'n_sample_frames',
              'buffer_ranges', 'buffer_kinds', 'valley_ratio_to_next', 'resolution_to_next',
              'chrom_R_to_next', 'warnings']


def write_peaks_csv(path, det):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS, extrasaction='ignore')
        w.writeheader()
        for r in peaks_rows(det):
            r = dict(r)
            r['window'] = f"{r['window'][0]},{r['window'][1]}"
            w.writerow(r)
    return path


def write_peaks_json(path, det, extra=None):
    payload = dict(ok=det['ok'], reason=det.get('reason', ''), q_range=det.get('q_range'),
                   smooth_window=det.get('smooth_window'), baseline_window=det.get('baseline_window'),
                   params=det.get('params'), window_mode=det.get('window_mode'),
                   dropped_apexes=det.get('dropped_apexes', []),
                   peaks=[{kk: vv for kk, vv in p.items()} for p in det['peaks']])
    for p in payload['peaks']:
        p['buffer_ranges'] = p.get('buffer_ranges', [])
        p['buffers'] = [[a, b, k] for a, b, k in p.get('buffers', [])]
    if extra:
        payload.update(extra)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w') as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False, default=str)
    return path


def summary_lines(det):
    """给控制台/README 的几行结论（结论先行，再说依据）。"""
    peaks = det['peaks']
    if not det['ok']:
        return [f"❌ 没有检测到洗脱峰：{det.get('reason','')}"]
    p0 = det['params']
    if det.get('manual_ranges'):
        head = (f"峰窗口是**人工指定**的（--peak-ranges）：{fmt_ranges(det['manual_ranges'])}"
                f"（不经过阈值判定；下列幅度/信噪比只是事后算出来的参照值）")
    else:
        head = (f"检测到 **{len(peaks)} 个洗脱峰**（色谱图 = 扣减后 q "
                f"{det['q_range'][0]}–{det['q_range'][1]} 1/Å 的积分强度；判据 = RAW 的 "
                f"savgol 平滑 + find_peaks，阈值 prominence≥{p0['min_prominence']}（相对最高峰）、"
                f"信噪比≥{p0['min_snr']}σ、半高宽≥{p0['min_width']} 帧"
                f"（基线噪声 σ={det.get('noise_sigma', float('nan')):.3g}））")
    lines = [head]
    for p in peaks:
        r = ''
        if p['index'] < len(peaks) - 1:
            ratio = adjacent_valley_ratio(det, p['index'])
            R = resolution_R(p, peaks[p['index'] + 1])
            r = (f" | 与下一峰：谷底/峰高={ratio:.2f}，R={R:.2f} "
                 f"（{label_verdict(ratio)} / {label_R(R)}）")
        lines.append(f"  · P{p['index']+1}: apex 帧 {p['apex']}，sample 窗 {p['window'][0]}–{p['window'][1]}"
                     f"（{p['window'][1]-p['window'][0]+1} 帧），幅度 {p['prominence']:.2f}×主峰，"
                     f"{p['snr']:.0f}σ，半高宽 {p['width']:.1f} 帧，"
                     f"buffer {fmt_ranges(p.get('buffer_ranges', []))}{r}")
        for wmsg in p.get('buffer_warnings', []):
            lines.append(f"      ！{wmsg}")
        if p.get('on_no_beam_frame'):
            lines.append("      ！峰顶落在无束流帧上 → 这条数据本身有问题，先看 assess-saxs-raw-data-quality")
        elif p.get('near_no_beam_frame'):
            lines.append("      ！峰顶紧邻无束流/断束帧 → 这个峰可能是野值造成的，看 sec_peaks.png 确认")
    if det.get('dropped_apexes'):
        lines.append(f"  ！还有 {len(det['dropped_apexes'])} 个更小的峰没分析（apex {det['dropped_apexes']}）"
                     f"，需要就调低 --peak-min-prominence / --peak-min-snr，或调大 --peak-max-n")
    diff, lo_seg, hi_seg, who = buffer_level_mismatch(peaks, det['y_corr'], top=det.get('top'))
    if diff >= 0.15:
        lines.append(f"  ！buffer 两侧水平差 {diff:.2f}×主峰（{who}：{fmt_ranges([lo_seg])} vs "
                     f"{fmt_ranges([hi_seg])}）→ 单靠峰前+峰后两段扣减会把斜漂移带进基线；"
                     f"先看要不要 --baseline linear/integral，或把 buffer 换成同一水平的两段，"
                     f"或用 --peak-buffer-max N 把每段截到紧邻峰的 N 帧内（水平匹配优先）")
    return lines


def needs_visual_check(det):
    """哪些情况**必须**上眼睛（返回原因列表；空 = 自动结果可以直接用）。"""
    reasons = []
    peaks = det['peaks']
    p0 = det.get('params', {})
    min_prom = p0.get('min_prominence', DEF_MIN_PROMINENCE)
    if not det['ok']:
        reasons.append("没有自动认出峰（阈值/数据问题）")
    for p in peaks:
        marginal = []
        if p['prominence'] < 3 * min_prom:
            marginal.append(f"幅度只有 {p['prominence']:.2f}×主峰（阈值 {min_prom}）")
        if np.isfinite(p.get('snr', float('nan'))) and p['snr'] < 15:
            marginal.append(f"信噪比 {p['snr']:.0f}σ 偏低")
        if marginal:
            reasons.append(f"P{p['index']+1}（apex {p['apex']}）刚过阈值：{('；'.join(marginal))}"
                           f" → 可能是纹波/肩峰，看图确认")
    for p in peaks[:-1]:
        ratio = adjacent_valley_ratio(det, p['index'])
        R = resolution_R(p, peaks[p['index'] + 1])
        if ratio >= 0.2 or R < 1.5:
            reasons.append(f"P{p['index']+1}/P{p['index']+2}：谷底/峰高={ratio:.2f}，R={R:.2f}"
                           f"（{label_verdict(ratio)} / {label_R(R)}）→ 两峰可能没真分开；"
                           f"若确认是未解析重叠 → 转 deconvolve-overlapping-elution-peaks")
    for p in peaks:
        for _, _, kd in p.get('buffers', []):
            if kd == 'between-peaks':
                reasons.append(f"P{p['index']+1} 的 buffer 段落在两峰之间（谷底），可能含混合组分")
                break
        if p['window'][1] - p['window'][0] + 1 < 10:
            reasons.append(f"P{p['index']+1} 的 sample 窗只有 {p['window'][1]-p['window'][0]+1} 帧"
                           f"（Guinier/IFT 会很不可靠）")
        if p.get('on_no_beam_frame'):
            reasons.append(f"P{p['index']+1} 峰顶在无束流帧上")
        elif p.get('near_no_beam_frame'):
            reasons.append(f"P{p['index']+1}（apex {p['apex']}）峰顶紧邻无束流/断束帧 → 可能是野值造成的假峰")
    if det.get('mask') is not None and np.asarray(det['mask']).any():
        reasons.append(f"有 {int(np.sum(det['mask']))} 帧疑似无束流（已在峰检测与 buffer 里剔除，"
                       f"仍应回原始帧确认）")
    diff, lo_seg, hi_seg, who = buffer_level_mismatch(peaks, det['y_corr'], top=det.get('top'))
    if diff >= 0.15:
        reasons.append(f"两侧 buffer 段水平差 {diff:.2f}×主峰（{who}）→ 峰前/峰后基线不在同一水平，"
                       f"先决定 baseline 校正还是换 buffer 段")
    return reasons


def _rel(path, out_dir):
    """README 里写相对路径（产物目录相对自己），别写字面绝对路径。"""
    if not path:
        return ''
    try:
        p = os.path.relpath(os.path.abspath(path), os.path.abspath(out_dir))
    except Exception:
        p = os.path.basename(path)
    return p if not p.startswith('..') else os.path.basename(path)


def write_peaks_readme(out_dir, prefix, det, png_overview=None, png_zoom=None, extra_lines=None):
    """多峰运行的**顶层** README（索引：几个峰、每个峰在哪、产物在哪个子目录）。"""
    peaks = det['peaks']
    lines = [f"# {prefix} —— SEC-SAXS **多峰**处理结果", "",
             f"> 自动生成 {det.get('timestamp', '')} ｜ 脚本 `run-raw-sec-pipeline.py` ｜ "
             f"（这份是**索引**；每个峰自己的完整结果说明在它的子目录里）", "",
             "## 0. 结论先行", ""]
    lines += summary_lines(det)
    lines.append("")
    reasons = needs_visual_check(det)
    if reasons:
        lines += ["**必须先看图的理由**（自动判定在这几处可能不对）：", ""]
        lines += [f"· {r}" for r in reasons]
        lines += ["", f"看图：`{_rel(png_overview, out_dir)}`（总览）与 "
                      f"`{_rel(png_zoom, out_dir)}`（逐峰放大），"
                      "确认后如果不对，用 `--peak-ranges` 手工给区间重跑；"
                      "两条流水线之外的重叠峰分解 → `deconvolve-overlapping-elution-peaks`。", ""]
    else:
        lines += ["自动判定没有触发任何「需要看图」的条件（仍建议扫一眼下面的图）。", ""]

    lines += ["", "## 1. 峰表", "",
              "| 峰 | apex 帧 | sample 窗（帧） | 帧数 | 幅度（×主峰） | SNR | 半高宽（帧） | "
              "buffer 段（帧） | 位置 | 与下一峰（谷底/峰高、R） | 子目录 |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
    for p in peaks:
        ratio = adjacent_valley_ratio(det, p['index']) if p['index'] < len(peaks) - 1 else None
        kinds = '、'.join(kd for _, _, kd in p.get('buffers', []))
        res = ('—' if ratio is None else
               f"{ratio:.2f} {label_verdict(ratio)}；R={resolution_R(p, peaks[p['index']+1]):.2f} "
               f"{label_R(resolution_R(p, peaks[p['index']+1]))}")
        lines.append(f"| **P{p['index']+1}** | {p['apex']} | {p['window'][0]}–{p['window'][1]} | "
                     f"{p['window'][1]-p['window'][0]+1} | {p['prominence']:.2f} | {p['snr']:.0f}σ | "
                     f"{p['width']:.1f} | "
                     f"{fmt_ranges(p.get('buffer_ranges', [])) or '（无）'} | {kinds or '—'} | "
                     f"{res} | "
                     f"`{p.get('subdir', '')}` |")
    lines += ["", "位置一列的含义：`edge-before` = 第一个峰之前、`edge-after` = 最后一个峰之后、"
                  "`between-peaks` = 两峰之间的谷底（**扣减用它会带进另一个组分**，要人工判断）。", ""]

    lines += ["## 2. 看图与手工复核", "",
              f"· `{_rel(png_overview, out_dir)}`：色谱图 + 峰窗 + buffer 段 + 逐帧 Rg/I(0)/MW。",
              f"· `{_rel(png_zoom, out_dir)}`：每个峰的放大图（峰是否回到基线、buffer 是否平）。",
              "· 看图判四问：① 峰的个数/位置与裁剪区视频或 UV 痕一致；② 每个峰的 sample 窗框住峰中心；"
              "③ buffer 段落在真基线（不能落在肩部/谷底）；④ 逐帧 Rg/MW 在每个峰上有平台。",
              "· `find-sec-peaks.py` 可以**不重跑积分**地改参数重算（秒级），"
              "`--sweep` 出一张灵敏度扫描图供人眼挑。", ""]
    if extra_lines:
        lines += extra_lines + [""]
    path = os.path.join(out_dir, 'README.md')
    os.makedirs(out_dir, exist_ok=True)
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines))
    return path
