#!/usr/bin/env python3
"""在**已有产物**上重算 SEC 洗脱峰（不重跑图像积分）：出峰表 + 两张图，给人工/agent 视觉复核。

为什么要它：``run-raw-sec-pipeline.py`` 的峰判定有阈值，而"到底有几个峰、边界在哪"这件事
代码判不干净（肩峰、扣减纹波、野值帧都会骗过阈值）。所以把"改参数看一眼"做成一个**秒级**的
独立工具：不动积分、不动已扣减的数据，只重算色谱图与峰区间，并把结果画成图给人看。
眼睛定下来之后，两种用法：① 用同一组阈值（``--peak-min-*`` 等）跑主管线；
② 直接把人工区间写回主管线：``run-raw-sec-pipeline.py --peak-ranges "..." --peak-buffers "..."``。

色谱图的算法与 RAW **逐点一致**：扣减后曲线在 q 窗口内的**梯形积分**
（= ``SASM.getIofQRange``），无束流帧的判据用整条曲线的**梯形积分**（= ``SASM.getTotalI``）；
两处都是"同样的点、同样的梯形"，所以本工具与主管线（走 RAWAPI 对象）算出来的数相同。
唯一差别：主管线在内存里用 RAW 的对象，本工具直接读 ``.dat`` 三列文本，因此快得多。

典型用法::

    # 1) 先看自动判出几个峰（秒级）
    find-sec-peaks.py --products processed/SEC-SAXS/4DH2-676-apo-3
    # 2) 觉得阈值不对 → 出一张灵敏度扫描图，用眼睛挑一组
    find-sec-peaks.py --products <产物目录> --sweep
    # 3) 按眼睛定下来的区间重算（并把可复用的参数打印出来）
    find-sec-peaks.py --products <产物目录> --peak-ranges "684,712;796,830"
    # 4) 想试另一段 buffer（这一步要 RAW，只重算扣减、不重跑积分，所以仍比主管线快）
    find-sec-peaks.py --products <产物目录> --buffer-range "560,620;800,880" --cfg data/20261001.cfg
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sec_peaks as sp                                          # noqa: E402

FMT = argparse.ArgumentDefaultsHelpFormatter


def log(msg):
    print(msg, flush=True)


def read_dat_matrix(path):
    """读一份三列 .dat（q I err）→ (q, I, err)；跳过 '#' 注释行。"""
    try:
        with open(path) as fh:
            txt = fh.read()
    except Exception:
        return None
    rows = [ln for ln in txt.splitlines() if ln.strip() and not ln.lstrip().startswith('#')]
    if not rows:
        return None
    try:
        arr = np.array(' '.join(rows).split(), dtype=float)
    except ValueError:
        return None
    if arr.size % 3 != 0:
        return None
    arr = arr.reshape(-1, 3)
    q = arr[:, 0]
    if len(q) > 2 and np.any(np.diff(q) < 0):                   # 万一是倒序，排回来
        order = np.argsort(q)
        arr = arr[order]
    return arr[:, 0], arr[:, 1], arr[:, 2]


def load_dir_matrix(files):
    """一目录的 .dat → (q, I[n_frame, n_q], err)。q 网格取第一份；点数不一致的帧丢掉并报告。"""
    q = None
    Is, errs, dropped = [], [], 0
    for f in files:
        d = read_dat_matrix(f)
        if d is None:
            dropped += 1
            continue
        if q is None:
            q = d[0]
        if len(d[0]) != len(q):
            dropped += 1
            continue
        Is.append(d[1])
        errs.append(d[2])
    return q, np.asarray(Is), np.asarray(errs), dropped


def load_arrays_for_figure(products, n_frames):
    """产物里的逐帧参数（上一轮跑出来的 Rg/I0/MW）→ 给总览图的三个面板；行数不匹配就不用。"""
    path = os.path.join(products, 'tables', 'frame_params.csv')
    if not os.path.exists(path):
        return None
    try:
        import csv
        rows = list(csv.DictReader(open(path)))
        if len(rows) != n_frames:
            return None
        g = lambda k: np.array([float(r[k]) for r in rows])      # noqa: E731
        return dict(rg={'Rg (previous run)': g('rg')},
                    i0={'I(0) (previous run)': g('i0')},
                    mw={'Vc MW': g('vc_mw'), 'Vp MW': g('vp_mw')})
    except Exception:
        return None


def sweep_combos():                                          # noqa: F811
    """灵敏度扫描的网格（与主管线共用 sec_peaks 的默认值）。"""
    return sp.sweep_combos()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=FMT)
    ap.add_argument('--products', required=True, help='产物目录（要有 profiles/01_integrated 与 03_subtracted）')
    ap.add_argument('--out', default=None, help='图与表的输出目录（默认 = --products）')
    ap.add_argument('--chromatogram', choices=['sub', 'unsub'], default='sub',
                    help='用哪条曲线找峰：sub=扣减后（组分分得开，但依赖 buffer 选得对）；'
                         'unsub=未扣减（不依赖 buffer，但基线漂移大）')
    ap.add_argument('--buffer-range', default=None,
                    help="给了就**重新扣减**再找峰（'s,e;s,e'，0 基帧号）；这一步用 RAW，"
                         "但仍不重跑图像积分。不给则直接用产物里已扣减的曲线")
    ap.add_argument('--cfg', default=None, help='RAW 设置档（只在 --buffer-range 时用得上）')
    ap.add_argument('--peak-q-range', default='%g,%g' % sp.DEF_Q_RANGE,
                    help='色谱图取的 q 窗口（1/A），逗号分隔')
    ap.add_argument('--peak-min-prominence', type=float, default=sp.DEF_MIN_PROMINENCE,
                    help='峰幅度下限（相对最高峰，0-1）')
    ap.add_argument('--peak-min-snr', type=float, default=sp.DEF_MIN_SNR,
                    help='峰高/基线噪声（1.4826*MAD）下限')
    ap.add_argument('--peak-min-height', type=float, default=sp.DEF_MIN_HEIGHT,
                    help='峰高下限（相对最高峰，0-1）')
    ap.add_argument('--peak-min-width', type=int, default=sp.DEF_MIN_WIDTH, help='半高宽下限（帧）')
    ap.add_argument('--peak-max-n', type=int, default=sp.DEF_MAX_PEAKS, help='最多分析几个峰')
    ap.add_argument('--peak-baseline-window', type=int, default=None,
                    help='滚动中位基线的窗口（帧）；默认按帧数自动（约 n/10，51–401）')
    ap.add_argument('--peak-smooth-window', type=int, default=None, help='savgol 平滑窗口（帧）；默认 min(51, n/2)')
    ap.add_argument('--peak-window', choices=['half-height', 'valley'], default='half-height',
                    help='峰窗取法：half-height=半高宽（保守）；valley=到相邻峰谷底（整个峰）')
    ap.add_argument('--peak-buffer', choices=['local', 'global'], default='local',
                    help='每个峰的 buffer 取法：local=紧邻的前后两段（默认）；global=全部非峰区共用')
    ap.add_argument('--peak-buffer-min', type=int, default=sp.DEF_BUFFER_MIN_LEN,
                    help='单侧 buffer 段的帧数下限（帧）')
    ap.add_argument('--peak-buffer-max', type=int, default=0,
                    help='单侧 buffer 段的帧数上限（帧；0=不限制）——基线漂移时用它把每段截到紧邻峰的那一侧')
    ap.add_argument('--peak-ranges', default=None,
                    help="**人工**给峰区间 'lo,hi;lo,hi'（0 基闭区间）→ 跳过自动识别（视觉复核后的落点）")
    ap.add_argument('--peak-buffers', default=None,
                    help="人工给每峰的 buffer 's,e;s,e|s,e'（'|' 分峰，';' 分段）")
    ap.add_argument('--sweep', action='store_true', help='额外出一张阈值灵敏度扫描图（供人眼挑参数）')
    ap.add_argument('--no-plots', action='store_true', help='只出表与 json，不出图')
    args = ap.parse_args()

    products = os.path.abspath(os.path.expanduser(args.products))
    out = os.path.abspath(os.path.expanduser(args.out or products))
    q_range = tuple(float(x) for x in args.peak_q_range.split(','))
    integ = sorted(glob.glob(os.path.join(products, 'profiles/01_integrated/*.dat')))
    subd = sorted(glob.glob(os.path.join(products, 'profiles/03_subtracted/*.dat')))
    if not integ and not subd:
        raise SystemExit(f'{products} 里没有 profiles/01_integrated 或 03_subtracted，先跑主管线')

    series_arrays = None
    # ---- 色谱图 + 无束流帧掩码
    if args.buffer_range:
        import bioxtasraw.RAWAPI as raw
        log(f'重新扣减：buffer {args.buffer_range}（这一步走 RAW，不重跑积分）')
        if not integ:
            raise SystemExit('--buffer-range 需要 profiles/01_integrated/*.dat')
        st = raw.load_settings(args.cfg) if args.cfg else None
        profs = raw.load_profiles(integ, st)
        series = raw.profiles_to_series(profs, st)
        bufs = sp.parse_ranges(args.buffer_range)
        subs, rg, _rger, i0, _i0er, vc, _vcer, vp = raw.set_buffer_range(series, bufs)
        y = sp.edf_frames(subs, *q_range)
        n = len(y)
        series_arrays = dict(rg={'Rg': np.asarray(rg, float)}, i0={'I(0)': np.asarray(i0, float)},
                             mw={'Vc MW': np.asarray(vc, float), 'Vp MW': np.asarray(vp, float)})
        mask_files = integ
    else:
        src = subd if (args.chromatogram == 'sub' and subd) else integ
        if args.chromatogram == 'sub' and not subd:
            log(f'！没有 profiles/03_subtracted（未扣减）→ 退回用未扣减曲线找峰')
        q, I, _err, dropped = load_dir_matrix(src)
        if q is None:
            raise SystemExit(f'{src[0]} 所在的目录里没有能读的三列 .dat')
        if dropped:
            log(f'！{dropped} 帧读不了/点数与首帧不一致，已跳过')
        lo, hi = max(q_range[0], q[0]), min(q_range[1], q[-1])
        if not lo < hi:
            raise SystemExit(f'q 窗口 {q_range} 与数据 q 范围 {q[0]:.4f}–{q[-1]:.4f} 不相交')
        sel = (q >= lo) & (q <= hi)
        y = np.trapezoid(I[:, sel], q[sel], axis=1)               # = RAW 的 getIofQRange
        n = len(y)
        mask_files = integ
        series_arrays = load_arrays_for_figure(products, n)

    mask = None
    if mask_files:
        _q2, I2, _e2, _d2 = load_dir_matrix(mask_files)
        if I2 is not None and len(I2) == n:
            mask = sp.no_beam_mask(np.trapezoid(I2, _q2, axis=1))   # = RAW 的 getTotalI

    # ---- 找峰（自动 或 人工区间）
    common = dict(mask=mask, buffer_mode=args.peak_buffer, buffer_min_len=args.peak_buffer_min,
                  buffer_max_len=args.peak_buffer_max)
    if args.peak_ranges:
        det = sp.manual_peaks(y, sp.parse_ranges(args.peak_ranges), **common,
                              note='人工区间（--peak-ranges）')
    else:
        det = sp.detect_peaks(y, q_range=q_range, smooth_window=args.peak_smooth_window,
                              baseline_window=args.peak_baseline_window,
                              min_height=args.peak_min_height,
                              min_prominence=args.peak_min_prominence, min_snr=args.peak_min_snr,
                              min_width=args.peak_min_width, max_peaks=args.peak_max_n,
                              window_mode=args.peak_window, **common)
    if args.peak_buffers and det['peaks']:
        sp.apply_manual_buffers(det['peaks'], args.peak_buffers, n)
    det['timestamp'] = 'find-sec-peaks.py'
    det['source'] = products
    det['chromatogram'] = args.chromatogram if not args.buffer_range else 'sub(rebuffered)'

    log('')
    for line in sp.summary_lines(det):
        log(line)
    reasons = sp.needs_visual_check(det)
    log('')
    if reasons:
        log('⚠️ 必须看图确认（自动判定在这几处可能不对）：')
        for r in reasons:
            log('   · ' + r)
    else:
        log('✅ 自动判定没有触发需要人眼复核的条件。')

    # ---- 落盘
    png = os.path.join(out, 'series', 'sec_peaks.png')
    png_zoom = os.path.join(out, 'series', 'sec_peaks_zoom.png')
    if not args.no_plots:
        sp.render_overview(png, det, series_arrays=series_arrays)
        sp.render_zoom(png_zoom, det)
        log(f'\n图：{png}')
        log(f'    {png_zoom}')
        for i, p in enumerate(det['peaks']):
            p['subdir'] = ''
    csv_path = sp.write_peaks_csv(os.path.join(out, 'tables', 'sec_peaks.csv'), det)
    json_path = sp.write_peaks_json(os.path.join(out, 'series', 'sec_peaks.json'), det)
    log(f'表：{csv_path}')
    log(f'    {json_path}')
    if args.sweep and not args.peak_ranges:
        sweep = os.path.join(out, 'series', 'sec_peaks_sweep.png')
        sp.render_sweep(sweep, y, sweep_combos(), mask=mask, q_range=q_range)
        log(f'灵敏度扫描图：{sweep}')

    # ---- 打印"可以原样搬回主管线"的参数（人工定下来的区间/buffer）
    if det['peaks']:
        log('\n可搬回主管线的参数（复制即用；只在你**人工**确认过区间/buffer 时才需要）：')
        log('  --peak-ranges "%s"' % sp.fmt_ranges([p['window'] for p in det['peaks']]))
        log('  --peak-buffers "%s"' % '|'.join(
            sp.fmt_ranges(p.get('buffer_ranges', [])) for p in det['peaks']))
        log('  阈值（若非默认）：--peak-min-prominence %g --peak-min-snr %g --peak-q-range "%g,%g"'
            % (args.peak_min_prominence, args.peak_min_snr, q_range[0], q_range[1]))


if __name__ == '__main__':
    main()
