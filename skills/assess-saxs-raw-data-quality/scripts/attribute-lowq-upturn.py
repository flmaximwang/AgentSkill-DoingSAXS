#!/usr/bin/env python
"""低 q 上翘归因（第一道检验）：是样品的散射，还是背景对消不干净？

对每个样品目录，用和流水线完全相同的操作（RAW average / scaleRelative / subtract）算三样东西：

  (a) 样品 - 缩放后的全部 control（流水线报出来的那条）
  (b) control run1 - 缩放后的 control run2      <- 纯 buffer 对纯 buffer，里面没有蛋白
  (c) 样品 - run1 单独 / 样品 - run2 单独        <- 换一个空白，结果会不会变

对每条差分曲线印同一套低 q 量：
      I(q1)/I(2*q1)      比例 > 1 = 上翘
      dI = I(q1)-I(2 q1) 绝对过量（counts）
      n  (q<0.02 的 -dlogI/dlogq)  形状：n~4 像寄生散射，n<1 像相互作用式的浅增强
      dI/I(q1)

判读：如果 (b)（空白-空白，绝对没有一个蛋白）给出的 dI 与 (a) 同量级，
      那个"上翘"就是背景，不是样品；如果 (a) 显著大于 (b) 的可复现极限，
      再拿浓度标度律去判（见 references/upturn-attribution-protocol.md）。

    <RAW python> attribute-lowq-upturn.py <root> --cfg <day>.cfg [--qmin 0.010]
"""
from __future__ import annotations

import argparse, glob, json, os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qc_common as qc                                    # noqa: E402


def metrics(d, qmin=0.0):
    q, I = d.getQ(), d.getI()
    i0 = int(np.argmin(abs(q - qmin))) if qmin > 0 else 0
    q, I = q[i0:], I[i0:]
    k1 = 1
    k2 = int(np.argmin(abs(q - 2 * q[k1])))
    lo = q < 0.02
    pos = I[lo] > 0
    n = float(-np.polyfit(np.log(q[lo][pos]), np.log(I[lo][pos]), 1)[0]) if pos.sum() > 5 else float("nan")
    return dict(q1=float(q[k1]), q2=float(q[k2]), i_q1=float(I[k1]), i_q2=float(I[k2]),
                ratio=float(I[k1] / I[k2]), dI=float(I[k1] - I[k2]),
                dI_over_I=float((I[k1] - I[k2]) / I[k1]), exponent=n,
                i_lowq_mean=float(np.mean(I[lo])))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--cfg", required=True)
    ap.add_argument("--hdr-format", default=qc.HDR_FORMAT)
    ap.add_argument("--qmin", type=float, default=0.0, metavar="Q",
                    help="[1/A] 先丢掉 beam-stop 边缘的点再算（0.010 是本机 BL19U2 的经验值）")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(args.root))
    s = qc.load_settings(args.cfg, args.hdr_format)

    report = {}
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
        sam = qc.average(sam_files, s)
        bl = {k: qc.average(v, s) for k, v in ctl_runs.items()}
        pooled = qc.average([f for v in ctl_runs.values() for f in v], s)

        out = dict(sample=key, sample_run=sk, n_sample=len(sam_files),
                   blank_runs={k: len(v) for k, v in ctl_runs.items()})
        a, fa = qc.subtract(sam, pooled)
        out["sample_minus_pooled"] = dict(metrics(a, args.qmin), scale=round(fa, 4))
        out["per_blank"] = {}
        for k, b in bl.items():
            c, f = qc.subtract(sam, b)
            out["per_blank"][k] = dict(metrics(c, args.qmin), scale=round(f, 4))
        ks = sorted(ctl_runs)
        if len(ks) >= 2:
            for x, y in ((ks[0], ks[1]), (ks[1], ks[0])):
                c, f = qc.subtract(bl[x], bl[y])
                out.setdefault("blank_minus_blank", []).append(
                    dict(metrics(c, args.qmin), scale=round(f, 4), pair="%s - %s" % (x, y)))
        report[key] = out

        m = out["sample_minus_pooled"]
        print("%-14s 样品-全部空白 : ratio=%5.2f  dI=%+9.2f  n=%4.1f  dI/I=%.2f"
              % (key, m["ratio"], m["dI"], m["exponent"], m["dI_over_I"]))
        for bb in out.get("blank_minus_blank", []):
            print("%-14s 空白-空白      : ratio=%5.2f  dI=%+9.2f  n=%4.1f  dI/I=%.2f   (%s)"
                  % ("", bb["ratio"], bb["dI"], bb["exponent"], bb["dI_over_I"], bb["pair"]))
        for k, v in out["per_blank"].items():
            print("%-14s   只减 %-12s: ratio=%5.2f  dI=%+9.2f  n=%4.1f"
                  % ("", k, v["ratio"], v["dI"], v["exponent"]))
        print("", flush=True)

    p = os.path.abspath(os.path.expanduser(args.out)) if args.out else os.path.join(
        root, "_rawqc", "upturn_attribution.json")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh:
        json.dump(report, fh, indent=1, default=str)
    print("wrote %s" % p)


if __name__ == "__main__":
    main()
