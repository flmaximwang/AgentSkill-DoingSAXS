#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把一个 SAXS 结果目录写成给人看的 `README.md`（SEC / 管式共用同一套章节）。

    python write-readme.py <结果目录> [<结果目录> ...] [--mode auto|sec|tube] [--stdout]

只读产物（`run_meta.json` / `summary.json` / `tables/*.csv` / `models/` / 目录清单），
**不重新拟合、不臆造数字**。管线（`run-raw-sec-pipeline.py` / `run-raw-tube-pipeline.py`）
在最后一步自动调用同一个 `readme_common.write_readme()`；对已经跑完、但没 README 的旧目录，
直接跑本脚本补写即可。

写完想确认「README 里的关键数字确实和产物对得上」，再跑：
    python verify-results-folder.py <结果目录>
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import readme_common as rc  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("out_dirs", nargs="+", help="结果目录（含 run_meta.json 或 summary.json）")
    ap.add_argument("--mode", default="auto", choices=["auto", "sec", "tube"],
                    help="结果目录类型；auto = 按目录里的标记文件判断")
    ap.add_argument("--stdout", action="store_true", help="只打印，不写文件")
    a = ap.parse_args()
    rc_ok = True
    for d in a.out_dirs:
        d = os.path.abspath(os.path.expanduser(d))
        mode = None if a.mode == "auto" else a.mode
        if not os.path.isdir(d):
            print("  ✗ 不是目录：%s" % d)
            rc_ok = False
            continue
        if a.stdout:
            text, keys, f = rc.render(d, mode)
            if text is None:
                print("  ✗ 认不出这是哪一类结果目录（没有 run_meta.json / summary.json）：%s" % d)
                rc_ok = False
                continue
            print(text)
            continue
        path, keys, f = rc.write_readme(d, mode)
        if path is None:
            print("  ✗ 认不出这是哪一类结果目录（没有 run_meta.json / summary.json）：%s" % d)
            rc_ok = False
            continue
        avail, _ = rc._availability(f)
        print("  → %s  [%s，%d 条关键数字，%d 个文件]"
              % (os.path.relpath(path, d), f["mode_label"], len(keys), len(f["files"])))
        print("     数据可用性：%s" % avail)
        fl = rc.flags(f)
        if fl:
            print("     ⚠️ 需要警惕：%s" % "；".join(fl[:3]))
        print("     下一步：python verify-results-folder.py %s" % d)
    return 0 if rc_ok else 1


if __name__ == "__main__":
    sys.exit(main())
