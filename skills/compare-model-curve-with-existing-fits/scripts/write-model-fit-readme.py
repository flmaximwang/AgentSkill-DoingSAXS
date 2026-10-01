#!/usr/bin/env python3
"""从 <产物目录>/comparison.json 生成人能看懂的 README.md（分点总结 + 明细表 + 🟢🟡🔴）。

只搬运产物里已有的数字，不重跑 crysol/datcmp、不重算 χ²（重算在 model-vs-fit.py 里）。
用法：write-model-fit-readme.py <产物目录>
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


def fmt(v, nd=3):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{nd}g}"
    return str(v)


def main() -> int:
    outdir = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    p = json.loads((outdir / "comparison.json").read_text())
    curves = p["curves"]
    used = [c for c in curves if c.get("ok", True)]
    good = [c for c in used if c.get("verdict") == "🟢"]
    warn = [c for c in used if c.get("verdict") == "🟡"]
    bad = [c for c in used if c.get("verdict") == "🔴"]
    failed = [c for c in curves if not c.get("ok", True)]
    models = [c for c in used if c.get("kind", "").startswith("本 skill")]
    others = [c for c in used if c not in models]
    rank = [c for c in used if c.get("primary_chi2") and not c.get("ranking_excluded")]
    best = min(rank, key=lambda c: c["primary_chi2"], default=None)
    others_rank = [c for c in rank if c in others]
    best_other = min(others_rank, key=lambda c: c["primary_chi2"], default=None)
    qmin = p.get("data_q_min") if p.get("data_q_min") is not None else p.get("qmin")
    qmax = p.get("data_q_max") if p.get("data_q_max") is not None else p.get("qmax")
    band = p.get("chi2_band_for_n_data") or [None, None]

    L: list[str] = []
    title = Path(p["processed"]).name if p.get("processed") else Path(p["data_original"]).stem
    L.append(f"# {title} —— model 曲线 vs 已有拟合\n")
    L.append(f"> 自动生成 {p['generated']} ｜ 脚本 `model-vs-fit.py` ｜ 数据 `{p['data_original']}` "
             f"｜ 拟合区间 q {fmt(qmin,4)}–{fmt(qmax,4)} Å⁻¹ ｜ ATSAS `{p['atsas_bin']}`\n")

    # ---- 结论速览
    L.append("## 结论速览\n")
    if bad:
        L.append(f"> 🟡 **需要警惕**：{len(bad)} 条曲线离实验曲线太远（{bad[0]['label']} "
                 f"red.χ²={fmt(bad[0].get('primary_chi2'))}）——它们不是「算错了」，"
                 f"是这套结构/重建不解释这条曲线。\n")
    if failed:
        L.append(f"> 🔴 **有 {len(failed)} 个模型没算出来**："
                 + "；".join(f"{c['label']}（{c.get('error')}）" for c in failed) + "\n")
    L.append("- **数据可用性**：用于比较的是 `data_used.dat`（"
             f"{fmt(qmin,4)}–{fmt(qmax,4)} Å⁻¹，{p.get('qmin') and '已用 datcrop 裁过' or '未裁剪'}）；"
             "本 README 的所有 χ² 都指**这一段 q**，不是整条曲线。\n")
    if models:
        L.append("- **本 skill 算出来的模型曲线**：")
        L.append("；".join(f"{c['label']} χ²(判据)={fmt(c.get('primary_chi2'))}{c.get('verdict','')}"
                          for c in models) + "\n")
    if best:
        L.append(f"- **判据 χ² 最小的一条**：{best['label']}（{fmt(best.get('primary_chi2'))}，"
                 f"口径 {best.get('primary_source')}，{best.get('verdict')}）。")
        if best_other and best.get("kind", "").startswith("本 skill"):
            L.append(f" 已有重建里最好的是 {best_other['label']}"
                     f"（{fmt(best_other.get('primary_chi2'))}，{best_other.get('primary_source')}）——"
                     + ("模型比它更贴。" if best["primary_chi2"] <= best_other["primary_chi2"]
                        else "模型**没有**比它更贴，说明高分辨结构还有问题（或它只是欠拟合）。") + "\n")
        else:
            L.append("\n")
    if good:
        L.append("- **表现较好的**：" + "；".join(
            f"🟢 {c['label']}（χ²={fmt(c.get('primary_chi2'))}）" for c in good) + "\n")
    if warn:
        L.append("- **表现较差的**：" + "；".join(
            f"🟡 {c['label']}（{'口径不同（.out 全 q）→ 不参与排序' if c.get('ranking_excluded') else c.get('reason')}）"
            for c in warn) + "\n")

    # ---- 明细表
    L.append("## 1. 对比明细表\n")
    L.append("> **怎么读**：三列 χ² 是三个口径，别混着引用 —— ①「程序自报」来自文件头/日志（DENSS 用 N 作分母，"
             "crysol 与 DAMMIF 用 N−1，实测同一个文件 2.665 / 4.233 / 4.334）；"
             "②「datcmp red.χ²」是 ATSAS 从**同一个文件**重算的统一口径（DENSS 那份 2.665 → 2.668）；"
             "③「同网格」是本 skill 把该曲线插到实验点网格、用实验误差重算 —— 网格本来就相同（标 `same_grid`）时它与 datcmp 一致，"
             "标「重采样」的说明那条曲线不在原始网格上（实测：原始 1150 点等距 3.8941e-4；DAMMIF 1149 点被重采样成 3.89e-4、与原始点最大差 4.95e-7 Å⁻¹；"
             "DENSS 845 点高 q 段更稀），插回原始网格的数只作参考。"
             "**判据一列（加粗）才是用来排序/判决的**：datcmp 可用就用它，datcmp 报哨兵时才退回同网格。「—（哨兵）」= datcmp 在这个文件上算不出有意义的值。"
             f"\n> 参考区间：χ²_n 的 95% 区间（n={p.get('chi2_band_for_n', '—')}）"
             f"= [{fmt(band[0])}, {fmt(band[1])}]。\n")
    L.append("| 曲线 | 来源 | 点数 | q 区间 Å⁻¹ | 程序自报 χ² | datcmp red.χ² | 同网格 red.χ² | "
             "**判据 χ²（口径）** | CorMap C / p | 模型 Rg Å | 能不能用 |\n"
             "|---|---|---|---|---|---|---|---|---|---|---|\n")
    for c in curves:
        dc = c.get("datcmp") or {}
        chi = (dc.get("chi-square") or {})
        cm = (dc.get("cormap") or {})
        cormap = (f"{fmt(cm.get('value'))} / {fmt(cm.get('p'))}" if cm else "—")
        if cm and not cm.get("usable"):
            cormap += "（哨兵）"
        reg = fmt(c.get("chi2_regrid"), 4) + ("" if c.get("same_grid") else "（重采样）")
        prim = (f"**{fmt(c.get('primary_chi2'), 4)}**（{c.get('primary_source')}）"
                if c.get("primary_chi2") is not None else "—")
        L.append("| {lab} | {kind} | {n} | {q0}–{q1} | {rep} | {dc} | {reg} | {prim} | {cm_} | {rg} | {v} {why} |\n".format(
            lab=c.get("label"), kind=c.get("kind"), n=fmt(c.get("n_points"), 5),
            q0=fmt(c.get("q_min"), 4),
            q1=fmt(c.get("q_max"), 4) if c.get("q_max") is not None else "—",
            rep=fmt(c.get("chi2_reported"), 4), dc=fmt(chi.get("value"), 4) if chi.get("usable") else "—",
            reg=reg, prim=prim, cm_=cormap,
            rg=fmt((c.get("stats") or {}).get("model_rg_A"), 3),
            v=c.get("verdict"), why=(c.get("reason") or "")))

    # ---- 参数
    L.append("\n## 2. 这次用的参数\n")
    L.append(f"- crysol：`--lm={p.get('lm')}`、`--smax={p.get('smax') or '默认 0.5'}`、"
             f"常数扣减 `{'开' if p.get('constant') else '关'}`、水化层 `{p.get('shell') or '默认 directional'}`、"
             f"角度单位 `{p.get('unit') or '自动猜'}`。\n")
    L.append(f"- 裁剪：{('datcrop --smin=%s --smax=%s' % (p.get('qmin'), p.get('qmax'))) if p.get('qmin') or p.get('qmax') else '未裁（整条曲线都参与拟合）'}。\n")
    L.append("- 复现命令：\n\n```bash\n" + "\n".join(p["commands"]) + "\n```\n")

    # ---- 目录
    L.append("\n## 3. 这个文件夹里有什么\n")
    L.append("| 路径 | 什么时候看 |\n|---|---|\n")
    L.append("| `README.md` | 现在这一份：结论 + 明细表 + 判据 |\n")
    L.append("| `comparison.csv` / `comparison.json` | 要数字/要喂给别的脚本时（json 里还有完整拟合参数） |\n")
    L.append("| `model-vs-fit.png` | 要看形状对不对、残差是不是随 q 成片偏移 |\n")
    L.append("| `fits/` | 从各程序目录**拷贝**过来的拟合文件（原始那份还在原处，没动过） |\n")
    L.append("| `crysol/<模型>/` | crysol 的原始产物：`.log`（参数与 χ²）、`.fit`（拟合四列）、`.int`/`.abs`（理论曲线）、`crysol.stdout.txt` |\n")
    L.append("| `data_used.dat` | 真正喂给 crysol 的那条曲线（裁过 q 就是裁过的版本） |\n")

    # ---- 判据
    L.append("\n## 4. 每个数字的判据\n")
    L.append("- **同网格 reduced χ²**：🟢 = 落在 χ²_n 的 95% 区间内且 CorMap p > 0.01（按 ATSAS 手册的 α=0.01 口径）；"
             "🟡 = 区间内但 CorMap 拒绝（残差有系统性形状偏差）、或超出区间不到 2 倍；🔴 = 超出区间 2 倍以上。\n")
    L.append("- **CorMap C / p**：检验残差的随机性，**不需要误差估计**（p 大 = 无法拒绝「拟合成立」）。"
             "曲线越长同样的 C 越显著，所以 n=1150 的曲线很容易被 CorMap 判成「有系统偏差」。\n")
    L.append("- **datcmp 的「哨兵」**：χ²≤0（如 BIFT 的 `.ift`）或 χ²=0 且 p=1.0（喂了普通 `.dat`）都表示"
             "**datcmp 在这个文件上算不出有意义的值**，绝不是「完美拟合」。\n")
    L.append("- **不同程序自报的 χ² 不可直接比**：DENSS 用 N 作分母（头里 2.665），crysol 与 DAMMIF 用 N−1 "
             "（本机实测 4.233 / 4.334），datcmp 又统一重算一遍（DENSS 那份变成 2.668）。\n")

    # ---- 没做的
    L.append("\n## 5. 这次没做的 / 不能信的\n")
    L.append("- **没有改模型**：本 skill 只回答「这个结构配不配这条曲线」，不做精修、不做系综、不改侧链。\n")
    L.append("- **没做绝对刻度归一**：χ² 与标度无关，但**不能**把 `.abs` 的绝对强度当浓度用（SEC 浓度未知）。\n")
    L.append("- **χ² 不是拟合好坏的唯一标准**：低 q 若有未被裁掉的背景/聚集，任何模型都会差 → 先回数据 QC 与 Guinier 判据。\n")
    L.append("- **datcmp 的 red.χ² 依赖文件里的误差估计**：SEC 扣减后的 σ 偏乐观（逐帧独立假设）→ 同一个模型在 SEC 数据上"
             "往往比在管式数据上更难达标，这是误差标定的问题，不是模型更差。\n")
    L.append("- 每个模型单独调过参数（lm/smax/constant）吗？没有 —— 同一次运行里所有模型用同一组参数，才可比。\n")

    # ---- 复核
    L.append("\n## 6. 想自己复核\n")
    L.append("- 重跑：把 `## 2` 里的命令按顺序跑一遍即可（模型、数据路径都是绝对的）。\n")
    L.append("- 看形状：`model-vs-fit.png` 上图叠加、中图标准化残差（残差成片同号 = 形状不对；逐点跳 = 噪声/误差被低估）。\n")
    L.append("- 换口径：加 `--no-constant` 再看一次 χ²（常数扣减＝允许吸收缓冲液失配，关掉它 χ² 会变差，"
             "这正是「失配有多大」的度量）。\n")

    # 行级后处理：段落/标题/列表之间补空行，压掉多余空行（markdown 才会分段显示）
    lines = "".join(L).split("\n")
    fixed: list[str] = []
    for ln in lines:
        starts_block = ln.startswith(("## ", "### ", "- ", "|", "> ", "```"))
        prev = fixed[-1] if fixed else ""
        cont = ((ln.startswith("- ") and prev.startswith("- "))
                or (ln.startswith("|") and prev.startswith("|"))
                or (ln.startswith("> ") and prev.startswith("> "))
                or (ln.startswith("```") and prev.startswith("```")))
        if starts_block and prev.strip() and not cont:
            fixed.append("")
        fixed.append(ln)
    text = re.sub(r"\n{3,}", "\n\n", "\n".join(fixed))
    (outdir / "README.md").write_text(text)
    print(f"{outdir / 'README.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
