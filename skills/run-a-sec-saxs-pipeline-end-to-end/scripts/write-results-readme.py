#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立入口：把已有结果目录写成给人看的 README.md（不重新处理数据）。

    /Applications/BioXTASRAW/bin/python write-results-readme.py --out <结果目录>

同样的实现也被 run-raw-sec-pipeline.py 在最后一步自动调用（见 scripts/results_readme.py）。
用在"已经跑完、但没有 README"的旧结果目录上最合适。
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from results_readme import write_results_readme  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--out", required=True, help="结果目录（含 run_meta.json 或 tables/）")
    ap.add_argument("--prefix", default=None, help="覆盖前缀（默认取 run_meta.json / 目录名）")
    a = ap.parse_args()
    if not os.path.isdir(a.out):
        raise SystemExit(f"不是目录：{a.out}")
    print("→", write_results_readme(os.path.abspath(a.out), a.prefix))


if __name__ == "__main__":
    main()
