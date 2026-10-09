# 字段表、阈值与浓度换算（`triage-saxs-batch-quality` 的深度层）

## 1. `<样品>/summary.json` 的字段（管式流水线实测）

| 路径 | 含义 |
|---|---|
| `subtraction.contrast_lowq` | **扣减后曲线口径**的低 q 对比度；与原始帧 QC 的对比度**不是同一个数** |
| `subtraction.control_scale_factor` | control 缩放因子；`|f−1|>5%` 报警（通量/位置不匹配） |
| `subtraction.negative_kratky_frac` | Kratky 为负的比例 |
| `guinier_auto.{rg,qmin,qmax,r_sqr}` | 自动区间拟合（**趋势用**，不是可引用值）；`qmin*rg` 落在 0.5–1.0 才算取点合理 |
| `guinier_recommended` | **非空 = 四条闸门（n_pts / qRg / chi2 / stable）全过**，这一档 Rg 可直接引用；为 `None` 时它的 `reason` 写明卡在哪一条 |
| `ift.chosen` | 采用哪个 IFT 引擎（`gnom` / `bift` / `dift`） |
| `ift.<engine>.{rg_realspace,i0,chisq,dmax,quality,gates}` | 实空间 Rg、I(0)、拟合 χ²、Dmax、**GNOM 自己的评语**、三条闸门 |
| `ift.trusted` | 流水线的可信闸门 —— **不等于** GNOM 的 `quality` 文本；两者分开报 |
| `mw.{Vc,Vp_porod,datmw_bayes,datclass}` | 三类浓度无关 MW + 形状分类（`compact` / `extended` / `random-chain`） |
| `shape.denss.*` | 电子云重建（`chi2` / `rg_model` / `support_volume`） |

`<样品>/tables/`：`guinier_multi_range.csv`（逐档闸门，自动区间不可用时看这个）、`guinier_results.json`、
`ift_summary.{csv,json}`、`mw.csv`、`shape_results.json`。

## 2. `_rawqc/raw_quality.md` 的列与阈值

列：结论 / 对比度(低q) / 缩放因子 / 扣减后上翘 / `I(0)/I(0.2)` / SNR均值(lowq) / 单帧qmax / 平均qmax /
帧间与σ之比（低 q 与 q=0.2 两套）/ 离群帧 / 问题。脚本把阈值写成可辩论的 `RULES` 字典（照抄自实测输出）：

```
contrast_fatal 0.005   contrast_warn 0.02    ctl_drift_warn 0.03   ctl_spread_warn 0.05
sam_spread_warn 0.05   sam_drift_warn 0.10   outlier_warn 0.20     jump_warn 0.15
qmax_warn 0.15         err_ratio_warn 4.0    upturn_warn 1.20      factor_warn 0.05
run_ratio_warn 0.05
```

判读要点：

- **低 q 判据看不见样品侧的异常帧**（低 q 被与背景共模的寄生散射主宰）⇒ 结论按 **q=0.2 那一套**；两套都跑。
- `I(0)/I(0.2)` 是一条不依赖拟合的分型：紧凑小蛋白十几到几十，大颗粒/聚集体上百到上千。
- 对比度是**原始帧口径**（分母是背景电平），绝对值可以到几百 %；别把它和 `summary.json` 的
  `contrast_lowq` 放同一张表里排序。

## 3. 为什么排序把「上翘」排在「对比度」前面

对比度度量的是「信号比背景大多少」，低 q 过量（聚集、边缘残留）与真实散射**都能把它抬起来**；
而上翘度量的是**形状对不对**——一条上翘的曲线不可能给出可信的 Rg/Dmax。所以：形状指标（上翘、单帧 qmax）
先定「能不能用」，强度指标（对比度）只在**同一形状档内**比大小。

## 4. 浓度：从稀释序列反推

### 4.1 不变量

`I(0) = K · c · M_w · (Δρ)²`（同一蛋白/同一溶剂下 K、(Δρ)² 固定）⇒
**`c · M_w` 相同 ⇒ `I(0)` 相同 ⇒ SNR 与可用 q 上界相同**。于是：

```
推荐浓度（新蛋白） = 参考浓度 × (参考 M / 新 M)
```

把标定写成 **`c·M`（mg·mL⁻¹·kDa⁻¹）**，就能跨蛋白、跨分子量搬运，不必每次重测一整条梯度。

### 4.2 判据三件

1. **够不够**：单帧 SNR≥2 的 q 上界 ≥ 目标 q（`q_need ≈ 6 / Rg`，Rg 用 `≈ 0.4·M^(1/3)` nm 粗估）。
   小蛋白（20–30 kDa，Rg 20–25 Å）⇒ q 到 **0.25 Å⁻¹** 才够做像样的 Guinier + Kratky 尾部。
2. **过没过量**：`I(0)/c` 在同一梯度里**恒定 ±5%** ⇒ 无 A2；下降 = 排斥、上升 = 聚集。
3. **低 q 平不平**：`I(q1)/I(2q1) ≤ 1.2`；越稀越翘是背景/边缘的特征（固定背景被相对放大），
   不是「越稀越聚集」。

### 4.3 工作例子（本机 BL19U2 管式，20 kDa 蛋白，2× 梯度 20→0.625 mg/mL）

| c (mg/mL) | I(0) | I(0)/c | 单帧 qmax (Å⁻¹) | 上翘 | 原始帧判定 |
|---|---|---|---|---|---|
| 20 | 74.8 | 3.74 | 0.248 | 1.04 | 可用 |
| 10 | 37.8 | 3.78 | 0.234 | 1.13 | 可用 |
| 5 | 19.2 | 3.84 | 0.209 | 3.04 | 可用 |
| 2.5 | 9.3 | 3.73 | 0.177 | 6.17 | 勉强可用 |
| 1.25 | 5.2 | 4.13 | 0.132 | 9.41 | 勉强可用 |
| 0.625 | 2.4 | 3.76 | 0.095 | 5.54 | 勉强可用 |

读法：`I(0)` 精确 2× 阶梯 ⇒ 浓度标定自洽；`I(0)/c` 恒定 ⇒ 全浓度段无 A2；
`qmax ≥ 0.23` 只在前两档出现 ⇒ **`c·M ≳ 200` 是「数据够用」的经验闸**。
换算到 30 kDa ⇒ 6.7–13 mg/mL；到 60 kDa ⇒ 3.3–6.7 mg/mL。

## 5. 数据来源清单（按可信度）

1. `<样品>/summary.json` + `tables/` —— 逐样品事实源。
2. `_rawqc/raw_quality.md` + `<样品>/qc.json` —— 逐帧事实源（**可能不覆盖全部样品**）。
3. 帧 header `*.txt`：`Description` / `Concentration [mg/ml]` / 曝光 / 波长 / 透射；
   **Concentration 字段常年为空**，别它为空就断言样品浓度未知——去备忘/纯化记录换算。
4. beamtime 备忘（`Memos-*/SAXS 实验记录 *.md` 一类）：口述配液步骤 + 换算浓度 + 待确认项；
   它自己会标注哪些是「换算」、哪些是「口述未说」。
5. 顶层 `README.md` / `Rg-MW-table.md` / `_summary/` —— **生成物**，引用前先确认覆盖范围。

## 6. 已知坑

- `_summary/` 可能已被并发会话搬走，而顶层 README 还在引用它 ⇒ 数字回样品目录取。
- `_rawqc/` 是快照：批中后加的样品不在里面 ⇒ 排序前先数覆盖数，缺的用 scratch 软链接根补跑。
- `ift.trusted` 与 `ift.<engine>.quality` 是两件事；报告里不要把前者当成后者的证据。
- `guinier_auto.rg` 只是**趋势**；`guinier_recommended` 为空时，任何 Rg 都要连「卡在哪条闸门」一起报。
- 分馏/子样（`-Void` / `-PeakN`）是**本体管的子样**：排序时先确认比的是本体还是子样，别把两者混成一列。
