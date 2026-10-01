#!/usr/bin/env python
"""低 q 上翘归因（第二道检验）：低 q 背景本身在不同 run 之间有多可复现？

对每一对 run（空白-空白、样品-空白）印不加任何缩放的比值 r(q)=A(q)/B(q) 在四个 q 上的值：

    r(q1)    q1 = 0.0068  最低可用点（上翘就住在这里）
    r(q2)    2*q1
    r(win)   0.30-0.44 的定标窗
    shape    r(q1)/r(win)：比值本身从定标窗到最低 q 变了多少

流水线只能扣掉一个**常数**因子（`scaleRelative`），所以 (shape-1) × (q1 处的背景电平)
会原封不动留在扣减曲线上，变成一个"看起来像样品的"低 q 上翘——哪怕里面一个蛋白都没有。

判读：把空白-空白 的 shape 当成背景可复现极限（本机 BL19U2 实测 ≤1.2%）。
      样品-空白 的 shape 若只比它大一点 -> 上翘是背景；大好几倍 -> 有真东西（再按浓度标度律判）。

    <RAW python> blank-reproducibility.py <root> --cfg <day>.cfg
"""
from __future__ import annotations

import argparse, glob, os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qc_common as qc                                    # noqa: E402

WIN = qc.SCALE_WINDOW


def ratio_profile(a, b, label):
    q = a.getQ()
    k1 = 1
    k2 = int(np.argmin(abs(q - 2 * q[k1])))
    m = (q >= WIN[0]) & (q <= WIN[1])
    r = a.getI() / b.getI()
    rq1, rq2, rw = float(r[k1]), float(r[k2]), float(np.mean(r[m]))
    level = float(a.getI()[k1])
    print("  %-30s r(q1)=%.4f r(q2)=%.4f r(win)=%.4f | shape=%.3f | 电平@q1=%.0f -> 残留 %+.0f counts"
          % (label, rq1, rq2, rw, rq1 / rw, level, (rq1 / rw - 1) * level))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--cfg", required=True)
    ap.add_argument("--hdr-format", default=qc.HDR_FORMAT)
    ap.add_argument("--out", default=None,
                    help="可选的 CSV，把每个 shape 收成一行（方便横向比样品）")
    args = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(args.root))
    s = qc.load_settings(args.cfg, args.hdr_format)
    rows = []

    for d in sorted(glob.glob(os.path.join(root, "*"))):
        if not os.path.isdir(d) or os.path.basename(d).startswith("_"):
            continue
        key = os.path.basename(d)
        files = glob.glob(os.path.join(d, "*.tif"))
        if not files:
            continue
        sk = qc.guess_sample_key(files, key)
        if sk is None:
            continue
        sam_files, ctl_runs, _ = qc.split_runs(files, sk)
        if not ctl_runs:
            continue
        print("== %s (%s | %d 帧)" % (key, sk, len(sam_files)))
        prof = {k: qc.average(v, s) for k, v in ctl_runs.items()}
        prof[sk] = qc.average(sam_files, s)
        ks = sorted(ctl_runs)
        for x, y in zip(ks, ks[1:]):
            ratio_profile(prof[x], prof[y], "%s / %s (空白-空白)" % (x, y))
            ratio_profile(prof[y], prof[x], "%s / %s (空白-空白 反)" % (y, x))
        for k in ks:
            ratio_profile(prof[sk], prof[k], "%s / %s (样品-空白)" % (sk, k))
        print("", flush=True)
        for k in ks:
            a, b = prof[sk], prof[k]
            q = a.getQ()
            m = (q >= WIN[0]) & (q <= WIN[1])
            rows.append((key, "%s/%s" % (sk, k), round(float((a.getI()[1] / b.getI()[1]) /
                                                             np.mean(a.getI()[m] / b.getI()[m])), 4)))
        for x, y in zip(ks, ks[1:]):
            a, b = prof[x], prof[y]
            q = a.getQ()
            m = (q >= WIN[0]) & (q <= WIN[1])
            rows.append((key, "%s/%s(空白)" % (x, y), round(float((a.getI()[1] / b.getI()[1]) /
                                                                  np.mean(a.getI()[m] / b.getI()[m])), 4)))

    if args.out:
        p = os.path.abspath(os.path.expanduser(args.out))
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            fh.write("sample,pair,shape\n")
            for r in rows:
                fh.write("%s,%s,%s\n" % r)
        print("wrote %s" % p)


if __name__ == "__main__":
    main()
