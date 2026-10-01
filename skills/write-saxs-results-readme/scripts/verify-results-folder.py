#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""验收一个 SAXS 结果目录：**结构齐不齐** + **README 里的关键结果与产物对不对得上**。

    python verify-results-folder.py <结果目录> [<结果目录> ...] [--mode auto|sec|tube]

四道检查（每一道都可能单独红）：
  1. 结构：按结果目录契约（两种模式各自的 must / should / info 清单）逐项点名。
  2. README 章节：8 节是否齐全、顺序是否与规范一致（SEC 与管式必须一样）。
  3. README ↔ 产物：README 里该出现的每条关键数字，是否能在产物里找到原文。
  4. 产物 ↔ 产物：两条独立路径（json ↔ csv）读同一个数是否一致 —— 不一致就是产物自相矛盾。

退出码：0 = 没有 🔴；1 = 有 🔴（交付前必须处理）。
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import readme_common as rc  # noqa: E402

ICON = {"🟢": "🟢", "🟡": "🟡", "🔴": "🔴"}


def verify_one(d, mode=None, show_all=False):
    d = os.path.abspath(os.path.expanduser(d))
    print("=" * 78)
    print("验收 %s" % d)
    print("=" * 78)
    if not os.path.isdir(d):
        print("🔴 不是目录")
        return 1
    mode, rows = rc.structure_check(d, mode)
    if mode is None:
        print("🔴 认不出这是哪一类结果目录（既没有 run_meta.json 也没有 summary.json）")
        return 1
    print("\n[1] 结构（%s，契约 %d 项）" % (rc.MODE_LABEL[mode], len(rows)))
    n_red = 0
    for r in rows:
        if r["status"] == "🟢" and not show_all:
            continue
        if r["status"] == "🔴":
            n_red += 1
        print("  %s %-42s %s" % (r["status"], r["path"], r["note"]))
    print("  （🟢 %d 项全部就位%s）"
          % (sum(1 for r in rows if r["status"] == "🟢"),
             "，未逐条列出：加 --all 看" if not show_all else ""))

    f = rc.collect(d, mode)
    readme = os.path.join(d, "README.md")
    if not os.path.exists(readme):
        print("\n[2-3] README：🔴 不存在 —— 先跑 "
              "`write-readme.py %s`（管线跑完会自动写）" % d)
        return 1
    text = open(readme, errors="ignore").read()
    keys = rc.build(f)[1]
    print("\n[2-3] README（%s，%d 字节）" % ("README.md", len(text)))
    for name, st, note in rc.readme_check(d, text, keys, f):
        if st == "🔴":
            n_red += 1
        print("  %s %-34s %s" % (st, name, note))

    cs = rc.cross_source_checks(d, mode, f)
    if cs:
        print("\n[4] 产物 ↔ 产物（两条独立路径读同一个数）")
        for name, st, note in cs:
            if st == "🔴":
                n_red += 1
            print("  %s %-44s %s" % (st, name, note))

    print("\n结论：%s" % ("🔴 有 %d 处硬伤，交付前必须处理" % n_red if n_red
                        else "🟢 结构与关键结果都对得上"))
    return 1 if n_red else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("out_dirs", nargs="+", help="结果目录")
    ap.add_argument("--mode", default="auto", choices=["auto", "sec", "tube"],
                    help="结果目录类型；auto = 按标记文件判断")
    ap.add_argument("--all", action="store_true", help="结构清单里 🟢 的项也逐条列出")
    a = ap.parse_args()
    mode = None if a.mode == "auto" else a.mode
    rc_total = 0
    for d in a.out_dirs:
        rc_total |= verify_one(d, mode, a.all)
        print()
    return rc_total


if __name__ == "__main__":
    sys.exit(main())
