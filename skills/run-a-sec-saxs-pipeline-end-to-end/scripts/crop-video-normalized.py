#!/usr/bin/env python3
"""SEC-SAXS 帧序列 → "逐帧归一化后的裁剪区域"视频（不落归一化 tif）。

为什么这么做：SEC 的每帧曝光期间入射光强会漂，看原始 tif 的裁剪区会把"通量起伏"当成样品信号；
但把 2000 帧归一化后写成新 tif 又要几十 GB。所以在**读入内存的那一刻**乘上该帧因子
（factor_i = median(TB)/TB_i，TB 取自 emit-bl19u2-header-txt.py 产出的 normalization_factors.csv），
直接把像素流喂给 ffmpeg，只留 mp4 + 量化副产物。

画面约定（**两个轴都是从"大的那头"数**——BL19U2 视图给的那套坐标）：数组下标
  row = H-1-y       （y 从下往上数）
  col = W-1-x       （x **从右往左**数，即镜像）
实测标定（bsa_00061 对历史视频帧做相关：cols 724-794=W-1-680..W-1-750 → corr 0.66；若按 col=x 取 680-750 → corr 0.11）。
裁剪 → 固定灰阶窗口（抽样帧的 p0.5/p99.9 定死，逐帧自动拉伸会让漂移消失）
→ 8× 最近邻放大 → 左上角烧帧号 → rgb24 rawvideo 管道给 ffmpeg。

副产物：逐帧裁剪区 sum / 强度质心 CSV（= 束斑漂移监测）、uint16 堆栈 .npy（可选 --save-stack）、
整帧+红框+裁块放大的预览 PNG。
"""
import argparse
import csv
import glob
import os
import subprocess
import sys

import numpy as np

try:
    import fabio
except ImportError:
    fabio = None

from PIL import Image, ImageDraw, ImageFont

FMT = argparse.ArgumentDefaultsHelpFormatter
DEFAULT_FF = "/Users/maxim/.hermes/tools/ffmpeg-9.0.1-darwin-arm64/ffmpeg"


def read_norm_csv(path):
    """归一化因子表 → {stem: Transmitted_Beam}。"""
    tb = {}
    with open(path, newline="") as fh:
        r = csv.DictReader(fh)
        tb_col = [c for c in r.fieldnames if c.startswith("Transmitted_Beam")]
        if not tb_col:
            raise SystemExit(f"{path} 里没有 Transmitted_Beam 列（先跑 emit-bl19u2-header-txt.py）")
        for row in r:
            tb[row["stem"]] = float(row[tb_col[0]])
    if not tb:
        raise SystemExit(f"{path} 是空的")
    return tb


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=FMT)
    ap.add_argument("--series-dir", required=True, help="原始 tif 目录")
    ap.add_argument("--norm-csv", default=None,
                    help="normalization_factors.csv（提供每帧 Transmitted_Beam）；不给=不做逐帧归一化")
    ap.add_argument("--out", required=True, help="输出 mp4 路径")
    ap.add_argument("--x1", type=int, default=680, help="裁剪区 x 下界（横向）")
    ap.add_argument("--x2", type=int, default=750, help="裁剪区 x 上界（横向）")
    ap.add_argument("--y1", type=int, default=910, help="裁剪区 y 下界（纵向，**从下往上数**）")
    ap.add_argument("--y2", type=int, default=970, help="裁剪区 y 上界（纵向，**从下往上数**）")
    ap.add_argument("--x-origin", choices=["right", "left"], default="right",
                    help="x 从哪头数：right=col=W-1-x（BL19U2 视图的约定，默认）；left=col=x")
    ap.add_argument("--upscale", type=int, default=8, help="整数最近邻放大倍数（libx264+yuv420p 要求偶数边长）")
    ap.add_argument("--fps", type=int, default=20, help="播放帧率（2000 帧 @20fps ≈ 100 s）")
    ap.add_argument("--stride", type=int, default=1, help="抽帧步长；>1 时视频更快、更小")
    ap.add_argument("--limit", type=int, default=None, help="只处理前 N 帧（试跑用）")
    ap.add_argument("--save-stack", action="store_true", help="另存 uint16 堆栈 .npy（2000×61×71 ≈ 17 MB）")
    ap.add_argument("--csv-out", default=None, help="逐帧 sum/质心 CSV 路径（默认 <out>.centroid.csv）")
    ap.add_argument("--preview-out", default=None, help="预览 PNG 路径（默认 <out>.preview.png）")
    ap.add_argument("--ffmpeg", default=DEFAULT_FF, help="ffmpeg 可执行文件")
    ap.add_argument("--normalize", choices=["on", "off"], default="on",
                    help="on=按 TB 逐帧归一化（factor=median(TB)/TB_i）；off=原始像素")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(os.path.expanduser(args.series_dir), "*.tif")))
    if not files:
        raise SystemExit(f"{args.series_dir} 下没有 tif")
    if args.stride > 1:
        files = files[::args.stride]
    if args.limit:
        files = files[:args.limit]

    if args.norm_csv and args.normalize == "on":
        tb = read_norm_csv(args.norm_csv)
        tb_med = float(np.median([tb[os.path.splitext(os.path.basename(f))[0]]
                                  for f in files if os.path.splitext(os.path.basename(f))[0] in tb]))
        factors = []
        for f in files:
            stem = os.path.splitext(os.path.basename(f))[0]
            if stem not in tb:
                raise SystemExit(f"{stem} 不在 {args.norm_csv} 里")
            factors.append(tb_med / tb[stem])
        factors = np.array(factors)
    else:
        factors = np.ones(len(files))   # 没有监视器数据（或显式关掉）→ 不归一化
        print("（未提供 --norm-csv 或 --normalize off：不做逐帧归一化，因子全 1）")
    nrm_state = "on" if (args.norm_csv and args.normalize == "on") else "off"
    print(f"帧数 {len(files)} | 归一化 {nrm_state} | 因子 1–99%: "
          f"{np.percentile(factors, 1):.4f}–{np.percentile(factors, 99):.4f}（中位 1.0）")

    d0 = fabio.open(files[0]).data
    H, W = d0.shape
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    r1, r2 = H - 1 - args.y2, H - 1 - args.y1
    if args.x_origin == "right":
        c1, c2 = W - 1 - args.x2, W - 1 - args.x1      # x 从右数 → 镜像
    else:
        c1, c2 = args.x1, args.x2
    cw, ch = c2 - c1 + 1, r2 - r1 + 1
    print(f"图 {H}行×{W}列 | x {args.x1}-{args.x2}(从{'右' if args.x_origin=='right' else '左'}数) "
          f"y {args.y1}-{args.y2}(从下数) → 行 {r1}-{r2} 列 {c1}-{c2} | 裁块 {ch}行×{cw}列")

    # 第 1 遍：抽样定显示灰阶（在归一化后的像素上定，逐帧自动拉伸会让漂移消失）
    sample = []
    for i in range(0, len(files), max(1, len(files) // 80)):
        sample.append(fabio.open(files[i]).data[r1:r2 + 1, c1:c2 + 1].astype(np.float32) * factors[i])
    s = np.concatenate([x.ravel() for x in sample])
    lo, hi = np.percentile(s, 0.5), np.percentile(s, 99.9)
    print(f"抽样 {len(sample)} 帧: lo(p0.5)={lo:.1f} hi(p99.9)={hi:.1f}")

    cmd = [args.ffmpeg, "-y", "-f", "rawvideo", "-pix_fmt", "rgb24",
           "-s", f"{cw*args.upscale}x{ch*args.upscale}", "-r", str(args.fps), "-i", "-",
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", args.out]
    pipe = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 22)

    stack = np.zeros((len(files), ch, cw), dtype=np.uint16) if args.save_stack else None
    rows = []
    yy, xx = np.mgrid[0:ch, 0:cw]
    for i, p in enumerate(files):
        c = fabio.open(p).data[r1:r2 + 1, c1:c2 + 1].astype(np.float32) * factors[i]
        if stack is not None:
            stack[i] = np.clip(c, 0, 65535).astype(np.uint16)
        tot = float(c.sum())
        rows.append((i + 1, os.path.basename(p), factors[i], tot,
                     float((c * xx).sum() / tot) if tot else float("nan"),
                     float((c * yy).sum() / tot) if tot else float("nan")))
        v = np.clip((c - lo) / max(hi - lo, 1e-9) * 255.0, 0, 255).astype(np.uint8)
        im = Image.fromarray(v, "L").resize((cw * args.upscale, ch * args.upscale),
                                            Image.NEAREST).convert("RGB")
        dr = ImageDraw.Draw(im)
        dr.rectangle([0, 0, 210, 30], fill=(0, 0, 0))
        dr.text((6, 4), f"frame {i+1}  x{factors[i]:.3f}", font=font, fill=(255, 220, 0))
        pipe.stdin.write(np.asarray(im, dtype=np.uint8).tobytes())
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(files)}", flush=True)
    pipe.stdin.close()
    err = pipe.stderr.read().decode(errors="ignore")[-400:]
    rc = pipe.wait()
    if rc != 0:
        raise SystemExit(f"ffmpeg 失败 rc={rc}\n{err}")

    csv_out = args.csv_out or (os.path.splitext(args.out)[0] + ".centroid.csv")
    with open(csv_out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "file", "norm_factor", "crop_sum", "xcen", "ycen"])
        w.writerows(rows)

    if stack is not None:
        np.save(os.path.splitext(args.out)[0] + ".stack.npy", stack)

    preview = args.preview_out or (os.path.splitext(args.out)[0] + ".preview.png")
    full = Image.fromarray(np.clip(d0 / max(float(d0.max()), 1) * 255, 0, 255).astype(np.uint8), "L").convert("RGB")
    dr = ImageDraw.Draw(full)
    dr.rectangle([c1, r1, c2, r2], outline=(255, 0, 0), width=3)
    crop_v = np.clip((d0[r1:r2 + 1, c1:c2 + 1].astype(np.float32) - lo) / max(hi - lo, 1e-9) * 255,
                     0, 255).astype(np.uint8)
    crop_im = Image.fromarray(crop_v, "L").resize((cw * 6, ch * 6), Image.NEAREST).convert("RGB")
    canvas = Image.new("RGB", (full.width + crop_im.width + 20, max(full.height, crop_im.height)), (30, 30, 30))
    canvas.paste(full, (0, 0))
    canvas.paste(crop_im, (full.width + 20, 0))
    canvas.save(preview)

    n = len(files)
    print(f"\n视频: {args.out}  {os.path.getsize(args.out)/1e6:.1f} MB  "
          f"{n} 帧 / {args.fps} fps = {n/args.fps:.0f} s")
    print(f"副产物: {csv_out}" + (f" | {os.path.splitext(args.out)[0]}.stack.npy" if stack is not None else "")
          + f" | {preview}")


if __name__ == "__main__":
    main()
