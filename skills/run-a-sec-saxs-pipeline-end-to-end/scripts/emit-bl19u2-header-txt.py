#!/usr/bin/env python3
"""给 SEC-SAXS 系列生成"每帧一份的 BL19U2 header .txt"——RAW 用它做逐帧归一化。

线站 batch/管式模式会给每帧写一份 `<帧名>.txt`（头里带 `Transmitted_Beam`），SEC 模式只给
`<系列>_1.Iochamber`（光强监视器时间序列）+ `<系列>_00001.log`（每帧 endTime）。RAW 自己不会
把这两样拼起来，但它是**认**这种 txt 的：`SASFileIO.parseBL19U2HeaderFile` 按 `<图像路径去掉扩展名>.txt`
读取（所以 txt 必须和 tif 同名同目录），配合 settings 里 `ImageHdrFormat='BL19U2, SSRF'` +
`EnableNormalization=True` + `NormalizationList=[['/','Transmitted_Beam']]`，积分时把每帧除以该帧
txt 里的 `Transmitted_Beam`（`SASImage.integrateCalibrateNormalize` → `calcExpression`）。所以归一化
完全在 RAW 里做，**不需要写归一化后的 tif**。

本脚本只做这一步：监视器 + 采集日志 → **每帧 txt 写进源数据目录（与 tif 并排）**；因子表与元数据
写到 --out-dir/norm/ 供复核。

用法：
  python emit-bl19u2-header-txt.py --series-dir <原始 tif 目录> \
      --monitor <系列>_1.Iochamber --log <系列>_00001.log --out-dir <产物目录> [--lag-auto]

若 tif 旁边**已经有线站给的逐帧 txt**（batch/管式模式会给，SEC 模式有时也给），别重生成——用
`--from-txt-dir <目录>` 只把它们读成 `normalization_factors.csv`（给裁剪视频做逐帧归一化用），
RAW 那边本来就直接读这些 txt。
"""
import argparse
import csv
import datetime as dt
import glob
import json
import os

import numpy as np

FMT = argparse.ArgumentDefaultsHelpFormatter

# ---- 监视器/日志解析（与 correct-sec-saxs-baseline/scripts/frame_flux_normalization.py 同一套） ----


def parse_monitor(path):
    """监视器文件：'#Time:YYYY-MM-DD hh:mm:ss' 头 + 'hh:mm:ss  value' 行。"""
    t, v, date = [], [], None
    for line in open(path, errors="ignore"):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            if s.startswith("#Time:"):
                date = dt.datetime.fromisoformat(s.split(":", 1)[1].strip()).date()
            continue
        p = s.split()
        try:
            hh, mm, ss = p[0].split(":")
            if date is None:
                raise SystemExit(f"监视器文件缺少 '#Time:YYYY-MM-DD ...' 头: {path}")
            t.append(dt.datetime.combine(date, dt.time(int(hh), int(mm), int(float(ss)))).timestamp())
            v.append(float(p[-1]))
        except (ValueError, IndexError):
            continue
    if not t:
        raise SystemExit(f"监视器文件没有解析出数据点: {path}")
    return np.array(t), np.array(v)


def parse_log(path):
    """采集日志：'t_req t_meas endTime  /path/to/file.tif'。返回每帧 endTime（秒）与文件名。"""
    end, names = [], []
    for ln in open(path, errors="ignore"):
        p = ln.split()
        if len(p) >= 4 and ln.rstrip().endswith(".tif"):
            try:
                end.append(dt.datetime.fromisoformat(p[2]).timestamp())
            except ValueError:
                continue
            names.append(os.path.basename(p[-1]))
    if not end:
        raise SystemExit(f"采集日志没有解析出帧: {path}")
    return np.array(end), names


def runmed(x, w):
    h = max(0, w // 2)
    return np.array([np.median(x[max(0, i - h):min(len(x), i + h + 1)]) for i in range(len(x))])


def monitor_per_frame(tmon, vmon, end, exposure, interval, lag_frames, n, min_val=None):
    """每帧曝光窗口 [te-exposure, te] 内监视器的中位数；te = endTime + lag*interval。

    min_val: 低于该值的采样点视为"无束流"野值丢弃（BL19U2 的 Iochamber 文件开头有几个 ~1e-13 的
    未开束流采样；不丢会让个别帧的窗口落进这段，因子爆掉）。
    """
    raw = np.full(n, np.nan)
    dropped = 0
    for i in range(n):
        te = end[i] + lag_frames * interval
        m = (tmon >= te - exposure) & (tmon <= te)
        vals = vmon[m]
        if min_val is not None and vals.size:
            keep = vals > min_val
            dropped += int((~keep).sum())
            vals = vals[keep]
        raw[i] = np.median(vals) if vals.size else np.nan
    idx = np.arange(n)
    ok = np.isfinite(raw)
    if ok.sum() < n:
        print(f"  警告: {n - ok.sum()} 帧的曝光窗口内没有可用监视器采样，已线性插值")
    if dropped:
        print(f"  注: 丢弃 {dropped} 个低于阈值的监视器采样（视为无束流野值）")
    return np.interp(idx, idx[ok], raw[ok])


def scan_lag(tot, tmon, vmon, end, exposure, interval, smooth, min_val=None,
             scan=60, step=2, subsample=4):
    """用"检测器总计数 vs 监视器"的相关系数扫 lag（两端时钟可能不同步，BL19U2 实测约 48 s）。"""
    n = len(tot)
    sub = np.arange(0, n, subsample)
    d_sub = runmed(tot, smooth)[sub]
    d_sub = d_sub / d_sub.mean()
    rows = []
    for lag in range(-scan, scan + 1, step):
        sm = runmed(monitor_per_frame(tmon, vmon, end, exposure, interval, lag, n, min_val), smooth)[sub]
        sm = sm / sm.mean()
        rows.append((lag, float(np.corrcoef(d_sub, sm)[0, 1])))
    best = max(rows, key=lambda r: r[1])
    print("=== lag 扫描（corr 越高越好）===")
    for lag, c in rows:
        if c > best[1] - 0.02:
            print(f"  lag {lag:>4} 帧 ({lag*interval:>6.0f} s)  corr {c:.4f}")
    print(f"→ 采用 lag {best[0]} 帧 (≈{best[0]*interval:.0f} s, corr {best[1]:.4f})")
    return best[0], rows


HDR = ("Description:                            {desc}\n"
       "Code:                                   \n"
       "Concentration [mg/ml]:                  \n"
       "Run Number:                             {run}\n"
       "Frame Number:                           {fn}\n"
       "Exposure time [s]:                      {expo}\n"
       "SAXS Position:                          {pos}\n"
       "Timestamp:              \t        {ts}\n"
       "SR Current:                             {sr}mA\n"
       "Transmitted_Beam:                       {tb}\n"
       "Wavelength [nm]:                        {wl}nm\n")


def write_header_txt(path, desc, run, frame, exposure, ts, sr, tb, wl, pos="saxs"):
    with open(path, "w") as fh:
        fh.write(HDR.format(desc=desc, run=run, fn=frame, expo=f"{exposure:.1f}", ts=ts,
                            sr=sr, tb=tb, wl=wl, pos=pos))


def read_tb_from_txt(path):
    """从一份 BL19U2 header txt 里取 Transmitted_Beam（找不到就报错，避免静默当成 1.0）。"""
    for ln in open(path, errors="ignore"):
        if ln.lower().startswith("transmitted_beam"):
            return float(ln.split(":", 1)[1].strip())
    raise SystemExit(f"{path} 里没有 Transmitted_Beam 行")


def factors_from_existing_txt(series, txt_dir, n_limit=None):
    """线站已给逐帧 txt 的情形：只读因子，不写任何 txt。"""
    tifs = sorted(glob.glob(os.path.join(series, "*.tif")))
    if not tifs:
        raise SystemExit(f"{series} 下没有 tif")
    if n_limit:
        tifs = tifs[:n_limit]
    missing = [f for f in tifs
               if not os.path.exists(os.path.join(txt_dir, os.path.splitext(os.path.basename(f))[0] + ".txt"))]
    if missing:
        raise SystemExit(f"{txt_dir} 里缺 {len(missing)} 份 <帧名>.txt（例如 "
                         f"{os.path.basename(missing[0])[:-4]}.txt）——线站 txt 没配齐，先补齐或用监控器模式生成")
    rows, tbs = [], []
    for i, tf in enumerate(tifs):
        stem = os.path.splitext(os.path.basename(tf))[0]
        tb = read_tb_from_txt(os.path.join(txt_dir, stem + ".txt"))
        tbs.append(tb)
        rows.append((i + 1, stem, "", "", f"{tb:.6e}"))
    return rows, np.array(tbs)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=FMT)
    ap.add_argument("--series-dir", required=True, help="原始 tif 目录（逐帧 txt 就写在这里，与 tif 并排）")
    ap.add_argument("--monitor", default=None, help="光强监视器文件（<系列>_1.Iochamber）")
    ap.add_argument("--log", default=None, help="采集日志（<系列>_00001.log，提供每帧 endTime）")
    ap.add_argument("--from-txt-dir", default=None,
                    help="线站已给逐帧 txt 时用：只从该目录读 Transmitted_Beam 生成因子表，不写 txt、不需要 monitor/log")
    ap.add_argument("--out-dir", required=True, help="产物目录（因子表/元数据写到 <out-dir>/norm/）")
    ap.add_argument("--txt-dir", default=None, help="逐帧 txt 的落盘目录（默认 = --series-dir，与 tif 并排）")
    ap.add_argument("--exposure", type=float, default=1.5, help="每帧曝光时长（s）；决定监视器取值窗口宽度")
    ap.add_argument("--interval", type=float, default=1.51, help="帧采样间隔（s，含曝光后空白）；时间轴推进与 lag 换算")
    ap.add_argument("--lag-frames", type=int, default=-32, help="监视器时间轴相对采集日志的偏移（帧）")
    ap.add_argument("--lag-auto", action="store_true", help="用检测器总计数扫 lag 并采用最优值（覆盖 --lag-frames）")
    ap.add_argument("--smooth", type=int, default=15, help="监视器滑动中位数窗口（帧）；抑制仪表逐点噪声")
    ap.add_argument("--monitor-min-frac", type=float, default=0.05,
                    help="丢弃低于（该比例×监视器全局中位）的采样点，视为无束流野值；0=不过滤")
    ap.add_argument("--tb-scale", type=float, default=1e8,
                    help="写进 txt 的 Transmitted_Beam = 监视器值×此常数（归一化只用到比值，尺度任意；"
                         "取 1e8 让强度与管式数据同量级，便于比较）")
    ap.add_argument("--description", default=None, help="写进 txt 的样品名（默认取目录名）")
    ap.add_argument("--run", type=int, default=1, help="写进 txt 的 Run Number")
    ap.add_argument("--sr-current", default="150.0", help="写进 txt 的 SR Current 占位值（mA）")
    ap.add_argument("--wavelength-nm", default="0.1033", help="写进 txt 的 Wavelength 占位值（nm）")
    ap.add_argument("--force", action="store_true", help="txt 已存在时覆盖（默认拒绝，避免盖掉线站原件）")
    ap.add_argument("--limit", type=int, default=None, help="只处理前 N 帧（试跑用）")
    args = ap.parse_args()

    series = os.path.abspath(os.path.expanduser(args.series_dir))
    out = os.path.abspath(os.path.expanduser(args.out_dir))
    txt_dir = os.path.abspath(os.path.expanduser(args.txt_dir)) if args.txt_dir else series
    norm_dir = os.path.join(out, "norm")
    os.makedirs(norm_dir, exist_ok=True)
    os.makedirs(txt_dir, exist_ok=True)

    if args.from_txt_dir:
        tdir = os.path.abspath(os.path.expanduser(args.from_txt_dir))
        rows, tbs = factors_from_existing_txt(series, tdir, args.limit)
        csv_path = os.path.join(norm_dir, "normalization_factors.csv")
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["frame", "stem", "endTime_epoch", "monitor_median", "Transmitted_Beam_x1"])
            w.writerows(rows)
        spread = float(100 * (np.percentile(tbs, 99) - np.percentile(tbs, 1)) / np.median(tbs))
        json.dump({"series_dir": series, "txt_dir": tdir, "n_frames": len(rows),
                   "source": "existing beamline txt (--from-txt-dir)", "tb_median": float(np.median(tbs)),
                   "tb_rel_spread_pct": spread,
                   "raw_settings_required": {"ImageHdrFormat": "BL19U2, SSRF", "EnableNormalization": True,
                                             "NormalizationList": [["/", "Transmitted_Beam"]]}},
                  open(os.path.join(norm_dir, "normalization_meta.json"), "w"), indent=2, ensure_ascii=False)
        print(f"=== 用线站已有 txt 生成因子表（不写 txt） ===")
        print(f"帧数 {len(rows)} | Transmitted_Beam 中位 {np.median(tbs):.4g} | 1–99 百分位相对展宽 {spread:.2f}%")
        print(f"因子表 : {csv_path}")
        print(f"元数据 : {norm_dir}/normalization_meta.json")
        return

    if not args.monitor or not args.log:
        raise SystemExit("需要 --monitor 与 --log（或改用 --from-txt-dir 读线站已有 txt）")

    tifs = sorted(glob.glob(os.path.join(series, "*.tif")))
    if not tifs:
        raise SystemExit(f"{series} 下没有 tif")
    tmon, vmon = parse_monitor(args.monitor)
    end, names = parse_log(args.log)
    n = min(len(tifs), len(end))
    if args.limit:
        n = min(n, args.limit)
    tifs, end = tifs[:n], end[:n]

    existing = [f for f in tifs
                if os.path.exists(os.path.join(txt_dir, os.path.splitext(os.path.basename(f))[0] + ".txt"))]
    if existing and not args.force:
        raise SystemExit(f"{txt_dir} 里已有 {len(existing)} 份同名 txt（例如 "
                         f"{os.path.basename(existing[0])[:-4]}.txt）——线站原件可能已在，"
                         f"确认要覆盖再加 --force")

    d = np.diff(end[:n]) if n > 1 else np.array([args.interval])
    med_d = float(np.median(d))
    print("=== 前置检查 ===")
    print(f"帧数 {n} | endTime 间隔中位 {med_d:.4f} s"
          f"{'' if abs(med_d-args.interval) < 0.02 else '  **与 --interval 不一致**'}"
          f" | 监视器 {len(tmon)} 点 / {(tmon[-1]-tmon[0]):.0f} s = {len(tmon)/(tmon[-1]-tmon[0]):.2f} 点/s")
    mon_med = float(np.median(vmon))
    min_val = args.monitor_min_frac * mon_med if args.monitor_min_frac > 0 else None
    print(f"监视器全局中位 {mon_med:.3e}；低于 {args.monitor_min_frac:.0%} 中位（无束流野值）的采样点 "
          f"{int((vmon <= min_val).sum()) if min_val is not None else 0} 个")

    lag = args.lag_frames
    lag_rows = None
    if args.lag_auto:
        try:
            import fabio
        except ImportError:
            raise SystemExit("--lag-auto 需要用装了 fabio 的 python（如 /Applications/BioXTASRAW/bin/python）")
        idx = np.arange(0, n, 4)
        tot = np.array([float(fabio.open(tifs[i]).data.sum()) for i in idx])
        lag, lag_rows = scan_lag(np.interp(np.arange(n), idx, tot), tmon, vmon, end,
                                 args.exposure, args.interval, args.smooth, min_val)

    tb_raw = monitor_per_frame(tmon, vmon, end, args.exposure, args.interval, lag, n, min_val)

    desc = args.description or os.path.basename(series.rstrip("/"))
    rows = []
    for i, tf in enumerate(tifs):
        stem = os.path.splitext(os.path.basename(tf))[0]
        tb = tb_raw[i] * args.tb_scale
        write_header_txt(os.path.join(txt_dir, stem + ".txt"), desc, args.run, i + 1, args.exposure,
                         dt.datetime.fromtimestamp(end[i]).strftime("%Y-%m-%d  %H:%M:%S.%f"),
                         args.sr_current, f"{tb:.6e}", args.wavelength_nm)
        rows.append((i + 1, stem, f"{end[i]:.4f}", f"{tb_raw[i]:.6e}", f"{tb:.6e}"))

    csv_path = os.path.join(norm_dir, "normalization_factors.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "stem", "endTime_epoch", "monitor_median", f"Transmitted_Beam_x{args.tb_scale:g}"])
        w.writerows(rows)

    meta = {
        "series_dir": series, "txt_dir": txt_dir, "monitor": args.monitor, "log": args.log, "n_frames": n,
        "exposure_s": args.exposure, "interval_s": args.interval, "lag_frames": lag,
        "lag_auto": bool(args.lag_auto), "smooth": args.smooth,
        "monitor_min_frac": args.monitor_min_frac, "tb_scale": args.tb_scale,
        "tb_first": float(tb_raw[0]), "tb_median": float(np.median(tb_raw)),
        "tb_rel_spread_pct": float(100 * (np.percentile(tb_raw, 99) - np.percentile(tb_raw, 1))
                                   / np.median(tb_raw)),
        "raw_settings_required": {"ImageHdrFormat": "BL19U2, SSRF", "EnableNormalization": True,
                                  "NormalizationList": [["/", "Transmitted_Beam"]]},
    }
    if lag_rows is not None:
        meta["lag_scan"] = [[int(l), round(c, 4)] for l, c in lag_rows]
    with open(os.path.join(norm_dir, "normalization_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2, ensure_ascii=False)

    print("\n=== 完成 ===")
    print(f"逐帧 txt : {txt_dir}/<帧名>.txt  ({n} 份，与 tif 并排，RAW 直接读)")
    print(f"因子表   : {csv_path}")
    print(f"元数据   : {norm_dir}/normalization_meta.json")
    print(f"Transmitted_Beam 中位 {np.median(tb_raw)*args.tb_scale:.4g}(x{args.tb_scale:g})，"
          f"1–99 百分位相对展宽 {meta['tb_rel_spread_pct']:.2f}%")
    print("RAW 端设置：ImageHdrFormat='BL19U2, SSRF'，EnableNormalization=True，"
          "NormalizationList=[['/','Transmitted_Beam']]")


if __name__ == "__main__":
    main()
