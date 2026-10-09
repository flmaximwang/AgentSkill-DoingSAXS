---
name: triage-saxs-batch-quality
description: "SAXS 批次结果：样品信号排序、浓度该取多少。"
---

# 一批 SAXS 结果的取用判断：谁好、谁废、下次该上多少浓度

## When to Use

用户把**一批已经跑完的**静态/管式 SAXS 结果摆出来，问的是**取用**而不是**判因**：

- 「把信号比较好的样品按从好到差排出来」「这批里哪个能用、哪个别用」「这两个系列比一比」。
- 「X kDa 的蛋白，静态 SAXS 该上多少浓度」「这个浓度够不够」「上次那批的浓度还能不能照用」。

判据固定成两类问题，分别回答、不要混：**① 排序**（每个样品一份可比较的数字 + 判据 + 数）；
**② 浓度**（用不变量换算出来的区间 + 上样配方 + 复核判据）。

不属于本 skill：原始帧质量**归因**（低 q 上翘是聚集还是光束挡边缘）→ `assess-saxs-raw-data-quality`；
单条曲线的 Rg 取点是否可信 → `assess-guinier-fit-quality`；分子量该报哪个方法 → `choose-a-molecular-weight-method`；
把一批原始帧跑成结果 → `run-a-tube-saxs-pipeline-end-to-end` / `run-a-sec-saxs-pipeline-end-to-end`；
结果目录缺什么/写 README → `write-saxs-results-readme`。本 skill **只读**，不重跑、不改产物。

## Step 0 · 先找数字源，别信汇总目录

事实源是**逐样品目录**：`<项目>/processed/<模式>/<样品>/`

| 要什么 | 在哪 |
|---|---|
| 对比度、Guinier 自动区间、IFT/GNOM、MW 各路 | `<样品>/summary.json`（字段表见 `references/triage-saxs-batch-quality-quality-metrics-and-fields.md`） |
| 多区间 Guinier 的逐档闸门 | `<样品>/tables/guinier_multi_range.csv` |
| 逐帧质量（判定/上翘/单帧 q 上界/帧间残差） | `<模式>/_rawqc/raw_quality.md` + `<样品>/qc.json` |
| 单个样品的自述 | `<样品>/README.md` |

**顶层 `README.md`、`Rg-MW-table.md`、`_summary/` 都是生成物**：会陈旧，也可能被并发会话重排或整个搬走
（实测：顶层 README 引用的 `_summary/` 目录已经不在，而各样品目录里的数字完好）。⇒ **数字一律回样品目录取**，
引用汇总表前先确认它列的数量与样品目录对得上。

同理，`_rawqc/` 是**一次性跑出来的快照，可能不覆盖你关心的样品**（批中后加的样品没有）——
排序前先数一遍 `_rawqc/` 里有几个样品，缺的自己补跑（Step 3）。

## Step 1 · 抽一张可比较的表

```bash
# 逐样品读 summary.json（Rg / I0 / 对比度 / IFT / 分子量一次拿全）
python3 - <<'PY'
import json, os
for s in sorted(os.listdir('.')):
    p = os.path.join(s, 'summary.json')
    if not os.path.isfile(p): continue
    d = json.load(open(p)); sub = d.get('subtraction', {}); g = d.get('guinier_auto', {}) or {}
    ift = d.get('ift', {}) or {}; gk = ift.get(ift.get('chosen'), {}) or {}
    print('%-24s c=%.4f Rg=%7.2f qmax*Rg=%.2f IFT_Rg=%7.2f chisq=%7.3f %s' % (
        s, sub.get('contrast_lowq', float('nan')), g.get('rg', float('nan')),
        (g.get('qmax') or 0)*(g.get('rg') or 0), gk.get('rg_realspace', float('nan')),
        gk.get('chisq', float('nan')), gk.get('quality', '-')))
PY
```

## Step 2 · 排序判据（顺序固定，逐条给数字，别给印象分）

前三条决定「这条数据能不能用」，后三条是程度；拟合侧只做复核。

1. **原始帧判定**：`可用` > `勉强可用`（`_rawqc/raw_quality.md` 的结论列，阈值在它的 RULES 里）。
2. **扣减后低 q 上翘** `I(q1)/I(2q1)`：越接近 1 越好，>1.2 报警。**这是排在对比度之前的第一号质量判据。**
3. **单帧 SNR≥2 的 q 上界**（单帧 qmax）：越高越好，决定你能覆盖到 qRg 多少；<0.15 Å⁻¹ 报警。
4. **低 q 对比度**：同口径内比（口径不同不能同表，见黑名单）；<2% 报警、<0.5% 不可用。
5. **帧间起伏 / RAW 报的 σ**：≤4 正常，>4 说明通量归一化残差占主导。
6. 拟合侧复核：`summary.json` 的 `guinier_recommended` **非空** = 四条闸门全过（可直接引用 Rg）；
   否则读 `ift.<chosen>.quality`（GNOM 自己的评语）与 `chisq`。**`ift.trusted=True` 是流水线的闸门，
   不等于 GNOM 说「a REASONABLE solution」**——两者要分开报，别拿一个当另一个的证据。

**输出形态**：一行一个样品、从好到差，每行给「判定 + 上面三条硬数字 + 一句为什么排在这」；
末尾写清**排除了什么、还剩什么**（例如：所有样品都有低 q 上翘，是否聚集未归因 ⇒ 只能说到「存在低 q 过量成分」）。

## Step 3 · 缺失的原始帧 QC：补跑，但绝不写进共享目录

共享的 `processed/` 树可能有别的会话在写/重排，**不要把你的 QC 产物落回那里**。用 scratch 软链接根 + `--out`：

```bash
S=~/.hermes/cache/scratch/<批次名>; mkdir -p "$S/root" "$S/out"
for n in <样品1> <样品2> …; do ln -s "<项目>/data/<模式>/$n" "$S/root/$n"; done
<RAW python> raw-qc-frames.py "$S/root" --cfg "<项目>/data/<日期>.cfg" --qmin 0.010 --out "$S/out"
```

- 软链接足够（脚本只读帧、只往 `--out` 下写）；**`--out` 不给就会写进 `<root>/_rawqc`，所以必须显式给**。
- `--qmin 0.010`：本机 BL19U2（D=2680 mm / λ=1.033 Å / Pilatus 2M）实测 q<0.010 全是光束挡边缘，判低 q 前先裁掉。
- 要写回项目（共享产物）**先问用户**，并按仓库既有的目录口径落。

## Step 4 · 浓度：用不变量换算，不用感觉

1. **不变量**：`I(0) ∝ c·M_w` ⇒ 同一个 `I(0)` 电平对应的浓度满足 **`c·M` 守恒**；
   把推荐值写成 **`c·M`（mg·mL⁻¹·kDa⁻¹）** 就能跨蛋白搬运：`c₂ = c₁ · M₁ / M₂`。
2. **决定性指标是单帧 q 上界，不是对比度**：要覆盖 `qRg ≈ 6`（小蛋白 Rg≈25 Å ⇒ q 到 ~0.25 Å⁻¹）
   就得保住「单帧 SNR≥2 ≥0.23 Å⁻¹」这一档。
3. **上限由溶解度定，不由浓度效应定 —— 用 I(0)/c 证明**：把稀释序列的 `I(0)/c` 列出来，
   **恒定在 ±5% 内 ⇒ 测不到 A2/粒子间干涉**；随 c **下降** = 排斥（Rg/I(0) 被低估）、**上升** = 聚集。
   结论里要写上「排除了浓度效应」这一句。
4. **浓度从哪来**：帧 header 的 `Concentration [mg/ml]` 字段**常年为空**（实测），浓度要去
   beamtime 备忘/纯化记录里按「质量 ÷ 体积 + 梯度倍数」换算；**换算完用数据验证**：
   2 倍梯度必须给出**精确 2 倍的 `I(0)` 阶梯**。验证不了的因子 2 歧义（例如「几个稀释点含不含原液」）
   **要明说并给出两个区间**，不要挑一个当答案。
5. **上样配方**：3–4 档跨 8–10 倍（最浓那档按第 1 条定），**同一根毛细管 buffer→样品→buffer 夹心**，每档 20 帧；
   对照组用**样品自己的 buffer**（不是 ddh2o——水对照会把低 q 抬起来，看起来像上翘）。

## 检查点

| 触发 | 动作 |
|---|---|
| 🛑 要把 QC/排序结果写回 `processed/` 或 `_rawqc/` | STOP：先问用户；共享树可能有并发写者 |
| 🛑 想「顺手」重跑某些样品让排序好看 | STOP：本 skill 只读；重跑走对应流水线 skill，且要先问 |
| 🔴 报告里要引用某个样品的 Rg/MW | CHECKPOINT：连它的 q 区间 / 闸门状态一起报，别只给数 |

## 🚫 黑名单

- **不要把「对比度最高」当成「信号最好」**：对比度是相对低 q 背景电平算的，低 q 过量（聚集、光束挡边缘残留）
  会把它抬高。实测反例：一批里对比度最高（524%）的那支上翘 2.51、无一条 Guinier 区间可选、MW 970 kDa。
- **不要把两种对比度口径混进一张表**：原始帧口径（相对背景电平，百分之几百）与扣减后曲线口径（22–70%）
  是两套数，混用会得出相反的排序。表头要写清是哪一种。
- **不要把「上翘」直接说成聚集/相互作用**：未做归因（空白−空白、背景可复现性、2D 差分）之前，
  只能写「存在低 q 过量成分」，并把归因交给 `assess-saxs-raw-data-quality`。
- **不要用印象排序**：没有数字支撑的「看起来还行」不算结论；每一条排序理由都要能指回字段或文件名。
- **不要引用顶层的汇总表当事实**：它可能陈旧（生成时间早于后加的样品）或被并发会话搬走；逐样品目录才是事实源。

## 本机标定值（校准点，不是通用常量）

20 kDa 蛋白（A5-05）6 点 2 倍梯度、BL19U2 管式（D=2680 mm / λ=1.033 Å / Pilatus 2M / 1 s × 20 帧），
浓度按备忘换算 20 / 10 / 5 / 2.5 / 1.25 / 0.625 mg/mL：

- `I(0)` = 74.8 / 37.8 / 19.2 / 9.3 / 5.2 / 2.4 —— 精确 2× 阶梯（浓度标定自洽）
- `I(0)/c` = 3.74 / 3.78 / 3.84 / 3.73 / 4.13 / 3.76 —— 跨 32 倍恒定 ±5% ⇒ 无 A2
- 单帧 SNR≥2 的 q 上界 = 0.248 / 0.234 / 0.209 / 0.177 / 0.132 / 0.095 Å⁻¹
- 扣减后上翘 = 1.04 / 1.13 / 3.04 / 6.17 / 9.41 / 5.54 —— **越稀越翘**（背景/边缘假信号的特征），不是聚集

⇒ 可用区间：保住 `qmax ≥ 0.23 Å⁻¹` 需要 **`c·M ≳ 200`**（mg·mL⁻¹·kDa⁻¹），上翘 ≤1.2 也落在同一档。
换算到 30 kDa：**7–13 mg/mL，起点 10**；下限 ~5 mg/mL（再稀 q 上界掉到 0.18 以下）。

## Support files

| 文件 | 承担什么 |
|---|---|
| `references/triage-saxs-batch-quality-quality-metrics-and-fields.md` | `summary.json` 字段表、`raw_quality.md` 的列与 RULES 阈值、判据顺序的理由、浓度换算与稀释序列判据的完整推导、数据来源清单与已知坑 |

## Skill Structure

<!-- Generated by Scripts -->

```
triage-saxs-batch-quality/
├── SKILL.md  (155 lines)
├── test-prompts.json  (41 lines)
└── references/
    └── triage-saxs-batch-quality-quality-metrics-and-fields.md  (98 lines)
```

<!-- Generated by Scripts -->
