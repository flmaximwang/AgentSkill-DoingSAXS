#!/usr/bin/env python
"""低 q 上翘归因（第三道检验，最直观）：看 2D 差分图。

把样品帧和 control 帧各自平均，按高 q 窗的因子缩放 control，做成差分图，然后回答一个问题：
低 q 的多余强度是**铺开的各向同性光晕**（真散射），还是**紧贴 beamstop 的一圈窄亮环**
（光束挡边缘/寄生散射，什么扣减都去不掉）——后者在 1D 曲线上照样表现为"低 q 上翘"。

同时印：
  * 最低 q 环带（r 6-30 px）里 样品/空白/差分 的均值与**方位角 RMS/均值**（各向同性判据），
  * 用 q_of_r 印出几个 q 对应的像素半径（判断"上翘点"在探测器上的位置），
  * 中心校验：手算径向平均必须复现 RAW 自己的 1D 形状（正确中心 RMS(log)~0.02，错的 ~0.22）。

    <RAW python> inspect-lowq-image.py <sample-dir> --cfg <day>.cfg [--out <png>]
"""
from __future__ import annotations

import argparse, glob, os, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qc_common as qc                                    # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("sample_dir", help="one sample directory (sample frames + control frames)")
    ap.add_argument("--cfg", required=True)
    ap.add_argument("--hdr-format", default=qc.HDR_FORMAT)
    ap.add_argument("--zoom", type=int, default=300, metavar="PIX",
                    help="half width of the zoom around the beam, in pixels")
    ap.add_argument("--out", default=None, help="default <sample-dir>/_lowq_image_check.png")
    args = ap.parse_args()
    d = os.path.abspath(os.path.expanduser(args.sample_dir))
    key = os.path.basename(d)
    s = qc.load_settings(args.cfg, args.hdr_format)

    files = sorted(glob.glob(os.path.join(d, "*.tif")))
    sk = qc.guess_sample_key(files, key)
    sam_files, ctl_runs, _ = qc.split_runs(files, sk)
    ctl_files = [f for v in ctl_runs.values() for f in v]

    imgs_s, _ = qc.raw.load_images(sam_files, settings=s)
    imgs_c, _ = qc.raw.load_images(ctl_files, settings=s)
    S = np.mean(np.asarray(imgs_s, float), axis=0)
    B = np.mean(np.asarray(imgs_c, float), axis=0)
    a, b = qc.average(sam_files, s), qc.average(ctl_files, s)
    f = qc.scale_factor(a, b)
    D = S - f * B

    cx, cy = qc.beam_center(S.shape, s)
    print("sample %s (%s) | %d sample + %d control frames | control scale=%.4f"
          % (key, sk, len(sam_files), len(ctl_files), f))
    print("beam at (%.1f, %.1f) px  [cfg Xcenter/Ycenter=%s/%s, row flipped]"
          % (cx, cy, s.get("Xcenter"), s.get("Ycenter")))
    for qq in (0.0068, 0.010, 0.02, 0.05):
        print("   q=%.4f -> r=%.1f px (%.2f mm)" % (qq, qc.r_of_q(qq, s),
                                                    qc.r_of_q(qq, s) * float(s.get("DetectorPixelSizeX")) / 1000.0))

    # centre check: manual radial average of the blank vs RAW's own 1D curve
    yy, xx = np.indices(S.shape)
    r_px = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)
    rr = np.arange(1, 300, 0.5)
    good = (np.asarray(s.get("Masks")["BeamStopMask"][0]) > 0) & np.isfinite(B) & (B > 0)
    man = np.array([B[(r_px >= r0 - .5) & (r_px < r0 + .5) & good].mean() for r0 in rr])
    qman = qc.q_of_r(rr, s)
    q1d, i1d = b.getQ(), b.getI()
    sel = (q1d > 0.03) & (q1d <= 0.4)
    mm = np.interp(q1d[sel], qman, man)
    mm, bb = mm / mm[0], i1d[sel] / i1d[sel][0]
    print("   centre check vs RAW 1D: RMS(log)=%.3f corr=%.3f"
          % (np.sqrt(np.mean((np.log(mm) - np.log(bb)) ** 2)), np.corrcoef(np.log(mm), np.log(bb))[0, 1]))

    for tag, img in (("sample", S), ("blank", B), ("diff", D)):
        ring = (r_px > 6) & (r_px < 30)
        v = img[ring]
        print("   %-6s ring(r6-30px): mean=%9.2f  azimuthal RMS/mean=%.3f"
              % (tag, v.mean(), v.std() / abs(v.mean()) if v.mean() else float("nan")))

    cxi, cyi = int(round(cx)), int(round(cy))
    z = args.zoom
    sl = (slice(max(cyi - z, 0), cyi + z), slice(max(cxi - z, 0), cxi + z))
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    for ax, (tag, img) in zip(axs, (("sample", S), ("blank x scale", f * B), ("difference", D))):
        sub = img[sl]
        vmax = np.percentile(sub, 99.5)
        ax.imshow(sub, cmap="inferno", vmin=0 if tag != "difference" else -vmax / 4, vmax=vmax)
        ax.set_title("%s | %s" % (key, tag), fontsize=10)
        ax.plot(cxi - sl[1].start, cyi - sl[0].start, "c+", ms=10, mew=0.8)
    plt.tight_layout()
    png = os.path.abspath(os.path.expanduser(args.out)) if args.out else os.path.join(
        d, "_lowq_image_check.png")
    plt.savefig(png, dpi=110)
    print("wrote %s" % png)


if __name__ == "__main__":
    main()
