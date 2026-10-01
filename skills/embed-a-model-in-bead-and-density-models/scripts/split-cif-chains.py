#!/usr/bin/env python
"""按链拆分 mmCIF：输出 <stem>_<CHAIN>.cif（保留原 header，只留该链的 _atom_site 行）。

    python split-cif-chains.py <in.cif> [--out-dir DIR] [--chains A,B] [--list]

环境里没有 Biopython（RAW/ATSAS 的 python 都没有），所以自己扫 _atom_site loop。
"""
from __future__ import annotations

import argparse
import os
import sys


def split_row(line):
    out, cur, q = [], "", None
    for ch in line:
        if q:
            if ch == q:
                q = None
            else:
                cur += ch
        elif ch in "'\"":
            q = ch
        elif ch.isspace():
            if cur:
                out.append(cur)
                cur = ""
        else:
            cur += ch
    if cur:
        out.append(cur)
    return out


def read_cif(path):
    """-> (header_lines, hdr_cols, data_lines) of the first _atom_site loop."""
    lines = open(path, errors="ignore").read().splitlines()
    for i, ln in enumerate(lines):
        if ln.strip() != "loop_":
            continue
        j, cols = i + 1, []
        while j < len(lines) and lines[j].strip().startswith("_atom_site."):
            cols.append(lines[j].strip()[len("_atom_site."):].split()[0])
            j += 1
        if not cols:
            continue
        rows, k = [], j
        while k < len(lines):
            s = lines[k].strip()
            if not s or s.startswith("#") or s.startswith("loop_") or s.startswith("_"):
                break
            rows.append(lines[k])
            k += 1
        return lines[:i], cols, rows
    raise SystemExit("no _atom_site loop in %s" % path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("cif", help="input mmCIF")
    ap.add_argument("--out-dir", default=None, help="output directory (default: same as input)")
    ap.add_argument("--chains", default=None, help="逗号分隔，只导出这些链（默认全部）")
    ap.add_argument("--list", action="store_true", help="只列出链与原子数，不写文件")
    args = ap.parse_args()

    header, cols, rows = read_cif(args.cif)
    ci = None
    for name in ("label_asym_id", "auth_asym_id"):
        if name in cols:
            ci = cols.index(name)
            break
    if ci is None:
        raise SystemExit("no chain column (label_asym_id/auth_asym_id) in %s" % args.cif)

    groups = {}
    for r in rows:
        p = split_row(r)
        if len(p) <= ci:
            continue
        groups.setdefault(p[ci], []).append(r)

    if args.list:
        for ch in sorted(groups):
            print("  chain %-3s %6d atoms" % (ch, len(groups[ch])))
        return

    want = sorted(groups) if not args.chains else [c.strip() for c in args.chains.split(",")]
    out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.cif))
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(args.cif))[0]
    for ch in want:
        if ch not in groups:
            print("! chain %s not in %s (have %s)" % (ch, args.cif, ",".join(sorted(groups))),
                  file=sys.stderr)
            continue
        out = os.path.join(out_dir, "%s_%s.cif" % (stem, ch))
        with open(out, "w") as fh:
            fh.write("\n".join(header) + "\n")
            fh.write("loop_\n")
            for c in cols:
                fh.write("_atom_site.%s\n" % c)
            for r in groups[ch]:
                fh.write(r.rstrip() + "\n")
            fh.write("#\n")
        print("  %s  (%d atoms, chain %s)  <- %s" % (out, len(groups[ch]), ch, args.cif))


if __name__ == "__main__":
    main()
