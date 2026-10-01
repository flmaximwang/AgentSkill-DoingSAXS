#!/usr/bin/env python
"""Aggregate many pipeline runs into one table
(<processed>/<mode>/_summary/{summary.csv,summary.md}).

Reads, for every sub-directory of <root> that contains summary.json:
  summary.json, tables/guinier_results.json, tables/ift_summary.csv, tables/mw.csv
and writes one row per sample with the numbers a human judges the run by.
"""
import argparse, csv, glob, json, os, time

ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument("root", help="directory holding one sub-directory per sample run")
ap.add_argument("--out", default=None)
args = ap.parse_args()
root = os.path.abspath(os.path.expanduser(args.root))
out = os.path.abspath(os.path.expanduser(args.out)) if args.out else os.path.join(root, "_summary")
os.makedirs(out, exist_ok=True)

COLS = ["sample", "n_frames", "n_frame_outliers", "control_scale", "contrast_lowq",
        "auto_I0", "auto_Rg", "auto_qmin", "auto_qmax", "auto_R2",
        "rec_Rg", "rec_qmin", "rec_qmax", "rec_chi2red", "rec_passes", "range_Rg_spread",
        "ift_trusted", "ift_start", "ift_qmin", "ift_Dmax", "ift_Rg_real", "ift_chisq",
        "ift_rg_ref_soft", "MW_Vp_kDa", "MW_Vc_kDa", "denss_chi2", "denss_Rg_model",
        "denss_support_A3", "note"]

rows = []
for d in sorted(glob.glob(os.path.join(root, "*"))):
    sj = os.path.join(d, "summary.json")
    if not os.path.isdir(d) or not os.path.exists(sj):
        continue
    s = json.load(open(sj))
    r = {c: "" for c in COLS}
    r["sample"] = s.get("sample", os.path.basename(d))
    r["n_frames"] = s.get("n_frames", "")
    r["n_frame_outliers"] = s.get("n_frame_outliers", "")
    sub = s.get("subtraction") or {}
    r["control_scale"] = round(sub.get("control_scale_factor", 0), 4)
    r["contrast_lowq"] = round(sub.get("contrast_lowq") or 0, 4)
    ga = s.get("guinier_auto") or {}
    r["auto_Rg"] = round(ga.get("rg", 0), 2)
    r["auto_qmin"], r["auto_qmax"] = round(ga.get("qmin", 0), 5), round(ga.get("qmax", 0), 5)
    r["auto_R2"] = round(ga.get("r_sqr", 0), 4)
    gr = os.path.join(d, "tables", "guinier_results.json")
    if os.path.exists(gr):
        g = json.load(open(gr))
        r["auto_I0"] = round((g.get("auto") or {}).get("i0", 0), 2)
        rec = g.get("recommended") or {}
        if rec:
            r["rec_Rg"] = round(rec["rg"], 2)
            r["rec_qmin"], r["rec_qmax"] = round(rec["qmin"], 5), round(rec["qmax"], 5)
            r["rec_chi2red"] = round(rec["chi2_red"], 2)
            r["rec_passes"] = int(all(rec["gates"].values()))
        rr = [x["rg"] for x in g["ranges"] if "failed" not in x]
        if rr:
            r["range_Rg_spread"] = "%.1f-%.1f" % (min(rr), max(rr))
    it = s.get("ift") or {}
    runs = it.get("runs") or []
    chosen = it.get("chosen")
    row = next((x for x in runs if x.get("tag") == chosen and "failed" not in x), None)
    r["ift_trusted"] = int(bool(it.get("trusted")))
    if row:
        r["ift_start"] = row["tag"]
        r["ift_qmin"] = round(row["qmin"], 5)
        r["ift_Dmax"] = round(row["dmax"], 1)
        r["ift_Rg_real"] = round(row["rg_realspace"], 1)
        r["ift_chisq"] = round(row["chisq"], 3)
        r["ift_rg_ref_soft"] = int(row.get("rg_ref_soft", False))
    mw = s.get("mw") or {}
    for k, col in (("Vp_porod", "MW_Vp_kDa"), ("Vc", "MW_Vc_kDa")):
        v = mw.get(k)
        if isinstance(v, dict) and v.get("mw") is not None:
            r[col] = round(v["mw"], 1)
    sh = (s.get("shape") or {}).get("denss")
    if isinstance(sh, dict):
        r["denss_chi2"] = round(sh["chi2"], 3)
        r["denss_Rg_model"] = round(sh["rg_model"], 1)
        r["denss_support_A3"] = int(sh["support_volume"])
    elif isinstance(sh, str):
        r["note"] = sh[:60]
    rows.append(r)

with open(os.path.join(out, "summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS)
    w.writeheader()
    w.writerows(rows)

with open(os.path.join(out, "summary.md"), "w") as fh:
    fh.write("| " + " | ".join(COLS) + " |\n|" + "---|" * len(COLS) + "\n")
    for r in rows:
        fh.write("| " + " | ".join(str(r[c]) for c in COLS) + " |\n")

# ---- series check: does I0 scale monotonically across a dilution series? -------------
import re
groups = {}
for r in rows:
    m = re.match(r"^(.*?)[-_]?(\d+)$", r["sample"])
    if m and r["auto_Rg"] != "":
        groups.setdefault(m.group(1), []).append(r)
with open(os.path.join(out, "summary.md"), "a") as fh:
    fh.write("\n## 稀释/系列检查\n\n")
    for pref, g in sorted(groups.items()):
        if len(g) < 3:
            continue
        g = sorted(g, key=lambda r: r["sample"])
        base = float(g[0]["auto_I0"] or 0)
        fh.write("序列 `%s*`（%d 个）\n\n" % (pref, len(g)))
        fh.write("| 样品 | I0 | I0/I0(第一) | auto Rg | q_max | R2 | 区间 Rg 跨度 | 通过闸门 | IFT 可信 |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|\n")
        for r in g:
            i0 = float(r["auto_I0"] or 0)
            fh.write("| %s | %s | %s | %s | %s | %s | %s | %s | %s |\n"
                     % (r["sample"], r["auto_I0"], ("%.3f" % (i0 / base)) if base else "",
                        r["auto_Rg"], r["auto_qmax"], r["auto_R2"],
                        r["range_Rg_spread"], r["rec_passes"], r["ift_trusted"]))
        fh.write("\n")
# ---- 整个结果文件夹的 README（给人看的"这里有什么、先看哪个、结论是什么"） --------------
def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


small, big = [], []
for r in rows:
    rg, mwvc = _num(r["auto_Rg"]), _num(r["MW_Vc_kDa"])
    (big if ((rg and rg > 40) or (mwvc and mwvc > 60)) else small).append(r)

with open(os.path.join(root, "README.md"), "w") as fh:
    fh.write("# Tube-SAXS 结果总览\n\n")
    fh.write("> 自动生成 %s ｜ 每个样品一个文件夹，文件夹里各自的 `README.md` 说明**里面每个文件**；"
             "本文件说明**整个文件夹**怎么用\n\n" % time.strftime("%Y-%m-%d %H:%M"))
    fh.write("## 先看这三样\n\n")
    fh.write("1. `_summary/overview.png` —— 所有样品一张图（扣减曲线 / Kratky / I(0) / Rg）\n")
    fh.write("2. `_summary/summary.md` —— 一行一个样品的汇总表（31 列，含判据）\n")
    fh.write("3. `<样品>/README.md` —— 单个样品的结果、每个文件的用途与判据（要具体数字时看这个）\n\n")
    fh.write("## 这批数据分两类（按 Rg 与 Vc 分子量自动分）\n\n")
    for title, group, note in (
            ("小蛋白一类（Rg ≤ 40 Å 且 Vc-MW ≤ 60 kDa）", small,
             "形态上像单一蛋白；低 q 能不能当真还要看各样品 README 里的对比度与 IFT 闸门"),
            ("大颗粒/聚集一类（Rg > 40 Å 或 MW 明显偏大）", big,
             "Rg 与分子量都远超单体，别按单体解读；先查样品制备/离心，再谈形状")):
        fh.write("### %s\n\n" % title)
        if not group:
            fh.write("（无）\n\n")
            continue
        fh.write("| 样品 | auto Rg (Å) | Vc 分子量 (kDa) | 低 q 对比度 | IFT 可信 | 目录里有什么 |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for r in group:
            has3d = "有 3D" if r["denss_Rg_model"] != "" else "无 3D"
            fh.write("| [%s](%s/README.md) | %s | %s | %s | %s | %s |\n"
                     % (r["sample"], r["sample"], r["auto_Rg"] or "—", r["MW_Vc_kDa"] or "—",
                        (str(r["contrast_lowq"]) and "%.1f%%" % (float(r["contrast_lowq"]) * 100))
                        if r["contrast_lowq"] not in ("", None) else "—",
                        "是" if r["ift_trusted"] == 1 else "**否**", has3d))
        fh.write("\n%s\n\n" % note)
    fh.write("## 目录结构\n\n")
    fh.write("- `<样品>/` —— 该样品的一次完整处理（曲线、表、P(r)、3D、报告、workspace）\n"
             "  - 里面每个文件的含义见该目录的 `README.md`\n"
             "- `_summary/` —— 汇总：`summary.md` / `summary.csv`（一行一样品）、`overview.png`\n"
             "- `_logs/` —— 每个样品的**完整运行日志**（某样品结果不对劲时先看这里）\n"
             "- `_rawqc/` —— 原始帧质量评估（如果跑过 `assess-saxs-raw-data-quality`）\n"
             "\n")
    fh.write("## 怎么再生 / 换参数\n\n")
    fh.write("```bash\n"
             "for d in <项目>/data/Tube-SAXS/*/; do\n"
             "  <RAW python> run-raw-tube-pipeline.py --sample-dir \"$d\" --cfg <项目>/data/<日期>.cfg \\\n"
             "      --out-dir <项目>/processed/Tube-SAXS/$(basename \"$d\") --denss-mode Fast\n"
             "done\n"
             "<RAW python> summarize-tube-run.py <项目>/processed/Tube-SAXS\n"
             "```\n\n")
    fh.write("- 想跳过光束挡边缘的点：加 `--qmin 0.010`（本机 BL19U2 实测 q<0.010 不可用）。\n"
             "- 想换 3D 引擎/不建模型：`--model-engine none`，或 `--denss-mode Slow`。\n"
             "- 重跑会覆盖该样品的全部产物（包括它的 README.md）。\n")

print("wrote %s (%d samples)" % (os.path.join(out, "summary.csv"), len(rows)))
for r in rows:
    print("%-16s frames=%-4s  scale=%-7s autoRg=%-6s recRg=%-6s pass=%-2s spread=%-12s "
          "IFT=%-2s Dmax=%-6s RgI=%-6s MWvp=%-6s MWvc=%-6s denssRg=%s"
          % (r["sample"], r["n_frames"], r["control_scale"], r["auto_Rg"], r["rec_Rg"],
             r["rec_passes"], r["range_Rg_spread"], r["ift_trusted"], r["ift_Dmax"],
             r["ift_Rg_real"], r["MW_Vp_kDa"], r["MW_Vc_kDa"], r["denss_Rg_model"]))
