---
name: compare-model-curve-with-existing-fits
description: "用 CRYSOL（ATSAS）把 model 算成 SAXS 理论曲线并拟合到扣减后的 .dat，再用 datcmp 与 processed/ 里已有的拟合（DENSS/DAMMIF/IFT）同口径比 reduced χ²，落叠加图与 README。用于「模型配不配这条曲线」「模型和珠模型/电子云谁更贴」「哪个链/哪条组装态更贴」「crysol 命令行参数怎么给」「χ² 各家算出来不一样」「datcmp 怎么用在 .fit/.fir/.out 上」「.cif 模型能不能算」；不负责该不该信某个参数 / RAW 界面里的取舍（转 fit-a-high-resolution-model-to-data），不负责形状重建本身（转 evaluate-a-shape-reconstruction），不负责把模型嵌进重建（转 embed-a-model-in-bead-and-density-models）。"
source_book: ATSAS 4.1.4 官方手册 *CRYSOL*（Introduction / 命令行选项 / 输出文件 / 示例 / 常见报错）· *DATCMP*（标准化残差定义 / 三个统计检验的性质 / Assessing Model Fit 示例）· 文件格式页 *.dat*、*.fit*、*.fir*
source_chapter: 源A manuals/crysol.md · manuals/datcmp.md；源B file-formats/{dat,fit,fir}.md
tags: [saxs, atsas, crysol, datcmp, chi2, model-fitting, theory-curve, like-for-like]
related_skills:
  - slug: fit-a-high-resolution-model-to-data
    relation: contrasts-with
  - slug: evaluate-a-shape-reconstruction
    relation: contrasts-with
  - slug: embed-a-model-in-bead-and-density-models
    relation: contrasts-with
  - slug: run-a-tube-saxs-pipeline-end-to-end
    relation: depends-on
  - slug: run-a-sec-saxs-pipeline-end-to-end
    relation: depends-on
---

# 「model 算出来的曲线」和「已有的拟合」不是同一个数：先化成同口径，再谈谁更贴

一件事有两种问法，本 skill 处理第二种：

- 问法一（**能不能信这个拟合**）→ `fit-a-high-resolution-model-to-data`：计算器选谁、harmonics/N samples 怎么调、`minimal 曲线` 什么时候骗人。那一条是 RAW 界面里的判据。
- 问法二（**把这件事跑成可复现的产物，并和已经躺在 `processed/` 里的拟合比一比**）→ 本 skill：ATSAS 的 `crysol` 命令行 + `datcmp` + 一条能横向比的 χ² 列 + 图 + README。

**核心一句：同一份数据、同一个模型，χ² 会因为「谁算的」而变——程序自报的分母不同、各家写文件的 q 网格不同、误差列的口径也不同。所以本 skill 不引用各家的自报数字来排序，而是让 ATSAS 的 `datcmp` 从**同一个文件**重算一遍（这才是手册为「评拟合」准备的口径），并且只把 `datcmp` 列用于排序与判决。**

## R — 原文 (Reading)

> This fitting is done by varying three parameters: (i) total displaced solvent volume (ii) contrast of the hydration shell (iii) relative background
>
> 出处：源A crysol.md · Introduction

> `--constant` or `-cst`｜Enables constant subtraction. This operation accounts for possible systematic errors due to mismatched buffers in the experimental data.
>
> 出处：源A crysol.md · Arguments and Options

> `--lm <N>`｜Maximum order of harmonics … This value must be increased whenever the maximum scattering angle is increased (smax).
> `--smax <SM>`｜Maximum scattering angle in inverse angstroms … default: 0.5 Å⁻¹, maximum: 2.0 Å⁻¹
>
> 出处：源A crysol.md · Arguments and Options

> If no prefix is provided, output file names are generated based on the base name of the inputs … With prefix, the two .log files are merged into one.
> The first line of the .fit file contains the following fields: Dro｜Optimal hydration shell contrast. RGT｜Radius of gyration (in Å) estimated from the theoretical curve. Vol｜Optimal excluded volume (in Å³). Chi^2｜Discrepancy between theoretical and experimental curves.
>
> 出处：源A crysol.md · crysol Output Files

> Columns in the output file are: s, I_exp, Err and I_sim.
>
> 出处：源B fir.md（.fir 列定义）· 同理 fit.md：*Fit of the calculated scattering curve versus the experimental data.*

> Further, if a single regularised SAS data (.out) file is provided as an argument, DATCMP calculates the fit of the reconstructed scattering from the p(r) to the experimental data. If a single fit to experimental data (.fir) or fit to calculated data (.fit) file is provided, DATCMP compares the contained experimental data and the model fit.
>
> 出处：源A datcmp.md · Introduction

> r(s) = (I₁(s) − I₂(s)) / sqrt(σ₁(s)² + σ₂(s)²) … If the two data I₁(s) and I₂(s) are not significantly dissimilar, the standardized residuals should follow a Standard Normal Distribution.
>
> 出处：源A datcmp.md · Standardized Residuals

> CorMap test｜Purpose: Randomness of residuals … does not require error estimates and is robust when errors are missing or unreliable … longer curves make long runs more significant, so the same C can imply different p-values depending on n
> red. χ² test｜Requires accurate error estimates … the reduced χ² value heavily depends on the error estimates, incorrect errors will lead to incorrect conclusions … with correct error estimates, a "close to 1.0" reduced χ² is generally considered a good fit; "close to 1.0" depends on the number of data points n：n≈50 → 0.5–1.5；n≈2000 → 0.9–1.1；the exact ranges may be calculated from the χ²_n-distribution for any n
> Anderson-Darling test｜more sensitive to deviations in the tails than reduced χ²
>
> 出处：源A datcmp.md · Statistical Tests

> `% datcmp ly01.fir` … Pair-wise Correlation Map test … C = 12.000000, Pr(>C) = 0.044513 … Assuming a cutoff α of 0.01, the Pr(>C) = 0.044 is greater than α, hence the assumption that the model scattering fits the experimental data can not be rejected, a.k.a. “the model fits the data”.
>
> 出处：源A datcmp.md · Examples — Assessing Model Fit

> SASDATA｜Two or more experimental SAS data (.dat) files, or a single regularised SAS data (.out), fit to experimental data (.fir) or fit to calculated data (.fit) file.
>
> 出处：源A datcmp.md · Arguments and Options

> if the maximum scattering vector value satisfies s_max < 1.0, the data are assumed to be given in inverse ångströms (unit 1)
>
> 出处：源B dat.md · Heuristic unit detection

## I — 骨架 (Interpretation)

**先把"拟合"这个词拆开。** 本 skill 一次运行里有两种拟合：

| | 谁做的 | 输入 | 产物 | 它的 χ² 指什么 |
|---|---|---|---|---|
| **模型拟合**（本 skill 新算） | `crysol`（ATSAS） | 一个 `.pdb`/`.cif` + 一条 `.dat` | `<模型>.fit`/`.log`/`.int`/`.abs` | 理论曲线 vs **这条**实验曲线 |
| **已有拟合**（`processed/` 里的） | DENSS / DAMMIF / GNOM / BIFT（流水线跑的） | 同一批数据 | `*_denss_map.fit`、`*_dammif_XX.fir`、`ifts/*.out` | 各自重建 vs 同一批实验数据 |

两者**本来就不是同一个量**，原因有三条，全部在本机实跑核对过（数据：`DataProcess_2026.10.01`，SEC 样品 `4DH2-676-apo-3`）：

1. **分母不同。** DENSS 写在 `.fit` 头里的 `chi2=2.665` 是 `Σ((I_exp−I_fit)/σ)²/N`；`crysol` 与 `DAMMIF` 的 `Chi^2` 是同一个式子除以 `N−1`（实测：crysol 报 4.233 = 我们按 N−1 复算的 4.233；DAMMIF 报 4.334 = 复算 4.334；DENSS 按 N 复算 2.665 = 它头里的数，按 N−1 是 2.668）。
2. **q 网格不同。** 同一批原始数据是 1150 点、等距 `3.8941e-4 Å⁻¹`；DAMMIF 把它重采样成 1149 点等距 `3.89e-4`（与原始点最大差 `4.95e-7 Å⁻¹`，相对 4.8e-6）；DENSS 是 845 点（高 q 段更稀）；`crysol` 的 `.fit` 则**正好写在实验点网格上**。把重采样过的曲线插回原始网格会凭空引入偏差（实测 DENSS：自己口径 2.668 → 插值后 4.40，虚高 65%）。
3. **误差列不同。** 同一批数据在不同程序里带的 σ 可能被重采样/平滑过，而 red. χ² 对 σ 极其敏感（手册原话：*incorrect errors will lead to incorrect conclusions*）。

**结论：排序与判决只用 `datcmp` 从同一个文件重算的那一列。**手册为这件事专门准备了入口——给 `datcmp` **单个** `.fit`/`.fir`，它就比"文件里的实验数据 vs 文件里的模型曲线"（`.out` 则是"p(r) 重建 vs 实验数据"），三个检验各答一个问题：

| 检验 | 答什么问题 | 需要的输入 | 本机实测（4DH2-676-apo-3，crysol 拟合） |
|---|---|---|---|
| CorMap（C, p） | 残差**随机**吗 | 只要曲线，**不需要误差估计** | C=41、p=0 → 拒绝"拟合成立"（曲线 1150 点，长曲线更容易显著） |
| red. χ² | 残差**有多大**（相对误差棒） | **误差必须准确** | 1.2668（与 crysol 自报 1.267 一致） |
| Anderson-Darling（A², p） | 残差**尾部**像不像正态 | 误差必须准确 | 本 skill 没跑，需要时 `datcmp --test=anderson-darling` |

**判档用哪条线：** 手册说 red. χ² 的"接近 1.0"随 n 变（它给 n≈50 → 0.5–1.5、n≈2000 → 0.9–1.1 两个锚点，并说"精确区间可由 χ²_n 分布算出"）。本 skill 就按这句话算 **χ²_n 的 95% 区间**（n=50 → [0.65, 1.43]，与手册的 0.5–1.5 相容；n=1150 → [0.92, 1.08]；n=519 → [0.88, 1.13]），再加 CorMap 的 `α=0.01`：

- 🟢 = red. χ² 落在区间内 **且** CorMap 不拒绝（p > 0.01）；
- 🟡 = 落在区间内但 CorMap 拒绝（**大小合格、形状不合格** —— 两个检验答的不是同一个问题）；或超出区间但不到 2 倍；
- 🔴 = 超出区间 2 倍以上（残差远大于误差棒）。

**几条一眼可用的推论：**

1. **`--constant` 不是可选项，而是"缓冲液失配有多大"的度量。** 同一模型同一数据：不加常数 χ²=4.233，加 `--constant` 后 1.267（同一份 `.fit` 头里写着 `Chi^2: 1.267`，且**第 4 列已经含了这个常数**——把常数再加一次会得到 10.65）。差这么多说明失配量级可观；关掉它再看一次，是了解自己数据的好办法。
2. **裁剪 q 不是万灵药，但必须做。** crysol 会拿你给它的**整条** `.dat` 去拟合，低 q 的 beamstop 光晕/聚集区对任何模型都拟不上，会主导 χ²。用 ATSAS 的 `datcrop` 先裁到可用区再拟合；**裁完还是差，那才是模型/寡聚态的问题**（实测管式 `A5-05-1` + `A5.pdb`：未裁 161.6 → 裁到 0.0567–0.2592 后 155.0，几乎没变 → 不是 q 区间的问题，见 A1-2）。
3. **不带 `.dat` 跑 crysol 得到的是 minimal 理论曲线**：只有 `.int`/`.abs`/`.alm`/`.log`，**没有 `.fit`、没有 χ²**。那是"生成"不是"拟合"，不要拿它下结论。
4. **`.cif` 可以喂 crysol**（手册：`.pdb`/`.cif`/`.ent` 都行）；这跟 PDB2SAS 只吃 `.pdb` 不一样——两个计算器的差别归 `fit-a-high-resolution-model-to-data`。
5. **`datcmp` 的文件类型是硬约束**：它只认 `.dat`/`.out`/`.fir`/`.fit`。喂一个**普通 3/4 列 `.dat`**，它会把自己跟自己比并给出 `χ²=0.000000, p=1.000000` —— 一个看着完美的**假结果**；喂 BIFT 的 `.ift` 会得到哨兵 `-1.000000`。判据：**χ²≤0 或 (χ²=0 且 p=1) = 这个文件不可判**，不是"拟合完美"。
6. **IFP/IFT 那份拟合的数不能直接跟模型比。** GNOM 的 `.out` 里既有 p(r) 表也有拟合块，`datcmp <xxx.out>` 给的是**拟合块里那段全 q**（含被裁掉的低 q 与 q=0 点）上的残差：实测 4.2526，而 IFT 汇总表里 GNOM 自报的是 0.4631。两个数都对，只是量的不是同一段 q/同一套误差模型 → 标注"口径不同，只作参考"。

## A1 — 案例 (Past Application)

数据全部来自 2026-10-01 的 BL19U2 批次（`~/Repositories/DataProcess_2026.10.01`），命令都是 `model-vs-fit.py` 一条命令跑出来的。

**1. SEC `4DH2-676-apo-3` + `4DH2-676-apo-dimer.cif`：模型比两个从头算重建都贴**

- 命令：`--model models/4DH2-676-apo-dimer.cif --data processed/SEC-SAXS/4DH2-676-apo-3/profiles/04_sample/sample_avg.dat --processed processed/SEC-SAXS/4DH2-676-apo-3`（默认 `--constant`）。
- 结果：模型（二聚体，crysol 报 Rg 25.13 Å）`datcmp` red. χ²=**1.267**；同一批数据上 DENSS 电子云 **2.668**、DAMMIF×4 **4.25–4.33**、GNOM IFT **4.25**（口径不同，只作参考）。
- 判档：模型 🟡——red. χ²=1.267 落在 n=1150 的区间 [0.92, 1.08] **外**但不到 2 倍，且 CorMap p=0（残差有系统性形状偏差）。一句话结论：**这个高分辨结构与数据相容，且比 845/1149 点的从头算重建更贴**；剩下的系统偏差要么来自 SEC 扣减误差偏乐观，要么来自真实构象涨落。
- 附带发现：`--constant` 从 4.233 降到 1.267；`--lm=100` 与 `--lm=20` 结果完全一样（球状二聚体，20 阶谐波已够——这与 RAW 侧"高长径比才需要加 harmonics"的口径一致）。

**2. 管式 `A5-05-1` + `A5.pdb`：χ²=155 不是参数问题，是分子状态问题**

- 命令：`--model models/A5.pdb --qmin 0.05666 --qmax 0.25915 --data processed/Tube-SAXS/A5-05-1/profiles/03_subtracted/subtracted.dat --processed processed/Tube-SAXS/A5-05-1`。
- 结果：`A5.pdb` red. χ²=**155.0**（🔴）；同一批数据上 DAMMIF×4 **1.22–1.48**、DENSS **5.42**、GNOM IFT 0.943。
- 判读：`A5.pdb` 是 116 残基单链、crysol 给 Rg **15.5 Å**、MW 1.8×10⁴ Da；而这条曲线的 Guinier/IFT 给 Rg **20.8 Å**、Dmax 72.6 Å、Vc-MW 20.9 kDa。**残差是"曲线整体形状对不上"（残差图里成片 ±20σ），不是噪声** → 该问的是"这条曲线是不是二聚体/其它组装态"，而不是继续调 crysol 参数。顺手证伪一个偷懒解释：裁到干净区（0.0567–0.2592）后 χ² 只从 161.6 变到 155.0，q 区间不是主因。
- 一句可复用的判据：**χ² 大到几十以上时，先比 Rg（模型 vs Guinier/P(r)），再看残差是否成片同号**；Rg 差 30% 就别指望拟合。

**3. `4LI2-676` 三条链 A/B/C：用拟合给"哪条链/哪种组装"排序**

- 命令：三个 `--model 4LI2-676_{A,B,C}.cif` 一次跑（脚本对每个模型各跑一次 crysol，产物目录分开），数据 `profiles/04_sample/sample_avg.dat`。
- 结果：链 C **2.320** < 链 B 2.478 < 链 A 2.671（crysol Rg 分别 15.80/15.80/15.97 Å，与样品 Rg 16.2 Å 对得上）；而同一批数据上的 DENSS **5.053**、DAMMIF×4 **4.86–4.90**。
- 结论：**高分辨模型全面优于从头算重建**（2.3–2.7 vs 4.8–5.1），链条间差异只有 0.1–0.35，属于"同一蛋白的不同链"，不足以据此判断哪条链"正确"——要靠化学计量与界面证据。
- 反面例子（同一次运行的另外两个）：`4DH1-KDPV-ZN` + `4DH1-KDPV.cif` χ²=53.7（模型 Rg 15.8 Å vs 样品 Rg 25.2 Å）、`4EH2-KDPV-ZN` + `4EH2-KDPV.cif` χ²=13.5（20.4 Å vs 26.0 Å）——都是**拿单体去拟二聚体样品**的典型症状。

**4. 四个踩过的坑（都是实测，不是推测）**

- `datcmp` 喂普通 4 列 `.dat`（`ifts/ift_fit.dat`）→ `Reduced Chi^2 = 0.000000, p = 1.000000`，两个"数据集"还被标成同一文件。**这是假结果**，因为 `.dat` 不是拟合格式。BIFT 的 `.ift` → `-1.000000`（算不出）。
- DENSS `.fit` 只有 845 点（重采样网格），把它插回 1150 点原始网格会把 χ² 从 2.668 抬到 4.40。**别用插值后的数去说 DENSS 更差**。
- DAMMIF 写了两个文件：`dammif_01.fir`（4 列，含误差）与 `dammif_01.fit`（3 列，绘图用、无误差列，q 从 0 开始）。`datcmp` 对后者给哨兵 `-1`。**取 `.fir`，别取 `.fit`。**
- `4LI2-676_D.cif` 里只有水（HOH），crysol 直接报 `error: No such model or chain in structure: 4LI2-676_D.cif`。**多链体系先核对每条链里有没有原子**（拆链脚本见 `embed-a-model-in-bead-and-density-models` 的 `split-cif-chains.py`），否则"某个模型没结果"会被误读成"这个模型不配"。

## A2 — 触发场景 (Future Trigger)

**用户会在什么情境下遇到这类问题**

- 手上有一条扣减好的 `.dat` 和一个模型，想知道**这个结构配不配**这条曲线。
- `processed/` 里已经躺着 DENSS/DAMMIF 的结果，想知道**高分辨模型是不是比这些从头算重建更贴**。
- 想给几个候选（不同链、单体 vs 二聚体、AlphaFold vs 晶体）**排序**。
- 发现"同一个拟合，为什么 RAW 里看是 2.66、ATSAS 日志里是 2.67、我手算是 4.4"。
- 想用命令行（而不是 RAW 界面）把这件事跑成可复现的产物给别人看。
- 拿到一个 `.cif`，不确定 crysol 能不能吃（能）。

**语言信号**

- 「这个模型配不配这条曲线」「帮我算一下理论曲线跟数据叠一叠」
- 「模型和珠模型/电子云谁更贴」「哪个链更像」
- 「crysol 命令行怎么跑」「`--constant` 要不要加」「`--lm` 给多少」
- 「χ² 各家算出来不一样」「datcmp 怎么用在 `.fit` 上」「datcmp 给我 χ²=0 是不是完美拟合」
- 「`.out` 里 GNOM 的 χ² 是 0.46，datcmp 却给 4.25」

**与相邻 skill 的区别**

- 与 `fit-a-high-resolution-model-to-data`：那边是 RAW 界面里的**判据与陷阱**（minimal 曲线、默认计算器、harmonics/N samples 的症状—处方）；本 skill 是**命令行产物 + 与已有拟合的同口径比较**。两边都会说"要拟合不要 minimal"，判据细节以那边为准。
- 与 `evaluate-a-shape-reconstruction`：那边评的是**从头算重建本身**（DAMMIF/DENSS 的平均 NSD、聚类、χ² 能不能收）；本 skill 是把**已有高分辨结构**放到同一堆数据边上比。两边都用 χ²，但只有本 skill 的 `datcmp` 列可以跨程序横着比。
- 与 `embed-a-model-in-bead-and-density-models`：那边是**刚体叠合**（NSD / fitmap correlation），不回答"配不配数据"；本 skill 是**拟合**。两者结论可以不一致：叠得好不等于配得上数据。
- 与两条端到端流水线（`run-a-sec-`/`run-a-tube-saxs-pipeline-end-to-end`）：那两条**生产** `.dat` 与各种拟合文件；本 skill **消费**它们。本 skill 不生成曲线以外的东西。

## E — 执行步骤 (Execution)

1. **确认 ATSAS 可用。** `--atsas-dir` 指到 **bin 这一级**（例如 `/Applications/ATSAS-4.1.4-1/bin`），脚本会自己把 `ATSAS` 环境变量补成 bin 的父目录（crysol 靠它找资源文件）。
   完成标准：`crysol`/`datcmp`/`datcrop` 三个可执行都在。判停点：只有 crysol 没有 datcmp → 版本太老，装 ATSAS ≥ 4。
2. **挑一条可信的实验曲线。** 用 RAW 扣减后落盘的 `.dat`（SEC 走 `profiles/04_sample/sample_avg.dat`，管式走 `profiles/03_subtracted/subtracted.dat`）。
   完成标准：能说出这条曲线的可用 q 区间（从数据 QC / Guinier 判据来）。🔴 STOP-1：低 q 有未裁掉的光晕/聚集 → 下一步给 `--qmin`，别指望拟合掩盖它。
3. **定模型与分子状态。** 单体还是二聚体？多链 `.cif` 要不要只取一条链（`--chain`）？有没有只有水的链？
   完成标准：能说出"我要算的是哪个组装态、文件里是哪几条链"。🔴 STOP-2：crysol 报 `No such model or chain in structure` → 先检查该文件里到底有没有原子。
4. **跑一次（一条命令，可直接粘）。** 脚本路径 = 本 skill 目录下的 `scripts/model-vs-fit.py`（装好的副本在 `~/.hermes/skills/saxs/compare-model-curve-with-existing-fits/scripts/`，仓库副本在 `~/Repositories/AgentSkill-DoingSAXS/skills/compare-model-curve-with-existing-fits/scripts/`）；解释器必须用 RAW 自带那只（有 numpy/scipy/matplotlib）。
   ```bash
   /Applications/BioXTASRAW/bin/python <skill>/scripts/model-vs-fit.py \
     --model /abs/模型.cif --model /abs/第二个.pdb \
     --data /abs/曲线.dat --processed /abs/processed/<模式>/<样品> \
     --atsas-dir /Applications/ATSAS-4.1.4-1/bin --qmin 0.0067 --qmax 0.2592 --out /abs/产物目录
   ```
   实测的控制台输出（本机 SEC 样例，逐行就是判据）：
   ```
   [data] …/data_used.dat  N=1150  q=0.00643–0.45386 Å⁻¹
   [crysol] 4DH2-676-apo-dimer: χ²(程序)=1.267 χ²(datcmp)=1.267 → 🟡
   [已有] DENSS 电子云: …_denss_map: χ²(datcmp)=2.668 → 🔴
   [已有] DAMMIF 珠模: …_dammif_01: χ²(datcmp)=4.334 → 🔴
   [plot] OK …/model-vs-fit.png      [readme] OK …/README.md      [best] …
   ```
   `--help` 里有全部参数与默认值（默认 `--constant` 开、`--lm=20`、`--smax` 用 crysol 默认 0.5）。
   完成标准：`crysol/<模型>/` 下出现 `.fit`/`.log`/`.int`/`.abs`，控制台打印每个模型与每条已有拟合的判据 χ²。🔴 STOP-2：模型全部失败 → 看 `crysol.stdout.txt`（常见：只有水/未知元素/氢原子数）再决定换文件还是加 `--explicit-hydrogens`。
5. **读三列 χ²，别混用。** 程序自报（各家分母不同）／`datcmp`（统一口径，**判决用这列**）／同网格重算（只在 `regridded=false` 时才与 datcmp 同级；标"重采样"的只作参考）。
   完成标准：能说出"我判决用的是哪一列、为什么"。
6. **看判档与两个检验。** red. χ² 对照 χ²_n 95% 区间；CorMap 看残差随机性。
   完成标准：给出"🟢 相容 / 🟡 形状还有系统偏差 / 🔴 不解释这条曲线"三档之一，而不是"看起来差不多"。🔴 STOP-3 / STOP-4：χ² 大到几十以上先比 Rg 与残差形状（别继续调参数）；red. χ² 在区间内但 CorMap 拒绝时别改判成🟢，结论要写成"量级合格、形状仍有系统偏差"。
7. **读图做归因。** `model-vs-fit.png`：上图看整体形状与高 q 端；中图标准化残差——**成片同号 = 形状/组装态不对**（例：A1-2 的 ±20σ），**逐点乱跳 = 误差被低估或噪声**（例：4LI2 高 q 段）。
   完成标准：能指出残差偏在哪个 q 区间、对应哪种物理原因。
8. **换口径复核（可选但便宜）。** 加 `--no-constant` 再跑一次（看缓冲液失配量级）；给 `--qmin/--qmax` 裁到干净区再跑一次（看低 q 污染的影响）。
   完成标准：两次运行的判据 χ² 差值能解释（常数项贡献 / 低 q 污染贡献）。
9. **交付。** 产物目录里 `comparison.csv`（每行一条曲线：三列 χ²、CorMap、判据、判档）、`comparison.json`（含完整参数与复现命令）、`model-vs-fit.png`、`fits/`（拷贝过来的拟合文件）、`README.md`（分点总结 + 明细表 + 🟢🟡🔴 + 没做的）。两条生成规则：①顶部「需要警惕」只统计**参与排序**的曲线（口径不同而被排除的 IFT 行单独在明细表里标「只作参考」），不许把两类混在一个计数里；②每条曲线的「判据 χ²」必须写明口径来源（datcmp / 同网格 / 程序自报），读者据此才能知道这一列为什么可以横着比。
   **给用户的回答写成三段式**（结论先行，≤5 句）：① 判据 χ² 是多少、什么口径、落在哪一档；② 与已有拟合的横向位次（谁最小、差多少）；③ 排除了什么、还剩什么、下一步做什么 —— 不许只甩一张表。
   🔴 STOP-6：交付前逐项核对 README 的「这次没做的 / 不能信的」。完成标准：别人只看 `README.md` 就能复述结论、并知道每个数字的口径。

## 🔴 停点（CHECKPOINT / STOP）

> 六道闸门：命中就**停下来先修上游或问用户**，不许带着已知缺陷往下跑。每道都对应一次实跑踩到的坑。

| 停点 | 触发条件 | 停下来做什么 |
|---|---|---|
| 🔴 STOP-1 | 实验曲线的可用 q 区间说不清（低 q 有未裁掉的光晕/聚集，或扣减未确认） | 先 `datcrop --smin/--smax` 裁到可用区并写明依据；不许靠拟合掩盖 —— 实测裁前 161.6 / 裁后 155.0，低 q 污染会主导 χ² |
| 🔴 STOP-2 | 模型的分子状态没定（单体 vs 二聚体、多链取哪条、某条链里有没有原子） | 先定状态再算：实测拿单体（Rg 15.5 Å）拟二聚体样品（20.8 Å）→ χ²=155，调参数救不回来 |
| 🔴 STOP-3 | 判据 χ² 大到几十以上 | 先比 Rg（模型 vs Guinier/P(r)）与残差形状，再决定换模型/换组装态；**禁止**继续调 `--lm`/`--constant`/`--smax` 试图把它压下来 |
| 🔴 STOP-4 | red. χ² 在区间内但 CorMap 拒绝（p ≤ α） | 结论必须写成「量级合格、形状仍有系统偏差」，**不许**写成「相容/通过」 |
| 🔴 STOP-5 | 产物目录里出现你没写过的文件或 `_prev_*`，或 `--out` 落在别人的结果目录 | 停下问用户（并发会话会重排目录）；改到独立 `--out`，并在报告里点名落点 |
| 🔴 STOP-6 | 交付前 | 逐项核对 README 的「这次没做的 / 不能信的」：排除了什么、还剩什么、下一步做什么 —— 缺这节等于没交付 |

**E 节里的对应停点**：step 2 → STOP-1；step 3/4 → STOP-2；step 6 → STOP-3、STOP-4；step 9 → STOP-6。

## F — 故障与兜底（if-then 三段式）

> 每条触发条件都在本机实跑里出现过；"一线修复"是这一次有效的动作，"仍失败兜底"是它不管用时的下一步。**不许静默跳过**：任一条触发后，结论里要写明走了哪一路。

| 触发条件（怎么发现） | 一线修复 | 仍失败兜底 |
|---|---|---|
| crysol 报 `ATSAS resource files not found` | 显式给 `--atsas-dir /Applications/ATSAS-4.1.4-1/bin`，或 `export ATSAS=<安装根>`（脚本会按 bin 的父目录补 `$ATSAS`） | `crysol --version` 确认可执行；版本 < 4 就换安装（`datcmp` 是判决列，缺它等于本 skill 无法排序） |
| crysol 报 `unable to determine number of hydrogens for …` | 先 `--explicit-hydrogens`（文件里已带 H 时） | 再 `--alternative-names`（老 PDB 命名）；仍失败用 `--implicit-hydrogen 0` |
| crysol 报 `unable to determine element for …` | `--alternative-names` | `--sub-element C`（或该原子真实元素）覆盖 |
| crysol 报 `No such model or chain in structure: X.cif`（实测：`4LI2-676_D.cif` 里只有 HOH） | 打开文件确认有没有原子（只有水/离子的链不算模型） | 换文件；或用 `--chain <ID>` 只取有原子的链；多链体系先拆链核对（拆链脚本在 `embed-a-model-in-bead-and-density-models`） |
| `datcmp` 给 `χ²=0.000000, p=1.000000` | 你喂的是普通 `.dat`（它把自己跟自己比）→ 换成 `.fit`/`.fir`/`.out` | 手上只有 `.dat` 时，用本 skill 的同网格 χ² 列自算，并在报告里标「datcmp 不适用」 |
| `datcmp` 给哨兵 `-1.000000` | 该文件不是可识别的拟合格式（BIFT `.ift`；DAMMIF 的 3 列 `.fit` 绘图文件） | 找同名 4 列 `.fir`（脚本已自动优先取它）；找不到就把该行标「不可判」、不参与排序 |
| 同一条 IFT 拟合出现两个 χ²（GNOM `.out` 4.253 vs IFT 表里 0.4631） | 说明口径：`.out` 的拟合块覆盖全 q（含被裁掉的低 q 与 q=0），IFT 表里那个是 GNOM 自己拟合区间/误差模型下的值 | 该行标「口径不同、只作参考」，并从排序中排除 |
| χ² 大到几十~几百，且残差图整体成片同号 | 先比 Rg：模型 vs Guinier/P(r)。差 30% 就别指望拟合（实测 `A5.pdb`：15.5 Å vs 样品 20.8 Å → χ²=155）→ 怀疑组装态（单体/二聚体） | 换组装态/换模型重跑；**不要**靠调 `--lm`/`--constant` 硬凑 |
| χ² 大，且残差集中在低 q | 用 `datcrop --smin/--smax` 裁掉光束挡晕/聚集区再拟合 | 裁完仍差 → 回 `assess-saxs-raw-data-quality` / `assess-guinier-fit-quality`，而不是继续拟合（实测裁 q 只把 161.6 降到 155.0，主因不在 q 区间） |
| red. χ² 在区间内但 CorMap 拒绝（p ≤ α） | 结论写成「量级合格、形状仍有系统偏差」，并在残差图上指出偏差所在的 q 段 | 偏差集中在高 q → 怀疑模型侧链/寡聚态；集中在低 q → 回数据 QC |
| `scipy` 不可用（算不出 χ²_n 区间） | 脚本自动退化为「只报数、不下判」 | 换 `/Applications/BioXTASRAW/bin/python` 跑（含 numpy/scipy/matplotlib） |
| 出图失败（缺字形 / 后端 / 负误差条） | 图内文字保持英文（脚本默认），报错只告警；数字与 README 已先落盘 | 用 `--no-plot` 先交付数字与 README，事后再单独补图 |
| `--processed` 目录里出现你没写过的 `_prev_*` 或新布局 | 说明别的会话刚重排过该目录 → **停下问用户**，别写进别人的输出目录 | 一律用 `--out` 指到独立目录；报告里点名产物落点 |
| 长跑中途被打断（回合插话会杀后台进程） | 重跑同一命令（同名产物覆盖，参数一致即幂等） | 按样品/模型拆 `--out` 分块跑，完成一块核对一块 |

## B — 边界 (Boundary)

**不要用的场景**

- 实验曲线本身不可信（背景没扣好、聚集、低 q 未裁）→ 先回 `assess-saxs-raw-data-quality` / `assess-guinier-fit-quality`；本 skill 的结论会跟着错。
- 想要**重建形状**（"这个蛋白长什么样"）→ `evaluate-a-shape-reconstruction` 与两条流水线；本 skill 只判"这个已知结构配不配"。
- 想把模型**塞进**珠模型/电子云看贴合度 → `embed-a-model-in-bead-and-density-models`。
- 想**改模型**（精修、柔性拟合、系综、多构象加权）→ 超出本 skill（那类程序：`sreflex`/`nmator`/EOM/BilboMD 等；手册体系里的刚性体/系综建模另有专页）。
- 想把 `.abs` 的绝对强度当浓度用 → SEC 浓度未知，χ² 与标度无关但浓度是另一件事（→ `choose-a-molecular-weight-method`）。

**源里明确警告过的失败模式**

- **拿普通 `.dat` 喂 datcmp**（实测）：给出 `χ²=0, p=1.0` 的假"完美拟合"。只喂 `.fit`/`.fir`/`.out`。
- **拿插值后的 χ² 跨程序比**（实测）：重采样曲线插回原始网格会虚高（DENSS 2.668 → 4.40）。
- **把各程序自报的 χ² 当同一口径**（实测）：分母 N vs N−1、网格不同、误差列不同。
- **关掉 `--constant` 却不说**：缓冲液失配会被算进模型误差里（本例 1.267 → 4.233）。
- **拿 minimal 曲线（不带 `.dat`）的"没有 χ²"当成程序出错**：那是没做拟合，不是失败。
- **靠"调参数"救一个寡聚态不对的模型**：Rg 差 30% 时（A1-2）任何参数都救不回来。

**材料盲点**

- 手册**没有给"χ² 多大算好"的绝对阈值**，只给与 n 相关的经验区间并说"精确区间可由 χ²_n 分布算出"。本 skill 的 95% 区间就是这句话的实现，**不是手册原文的表格**；手册的 red. χ² 前提是"误差估计准确"，SEC 扣减后的 σ 偏乐观会让判档偏严——所以**排序（相对）比绝对判档更可靠**。
- 手册**没有定义 crysol 日志里的 `Probability of fit`**（只说它出现）。**不要把它当 p 值用**；要 p 值用 `datcmp` 的。本 skill 只记录它、不解释它。
- **datcmp 的 `--mode` 语义**（PAIRWISE/INDEPENDENT）只在多文件比较时才有意义；本 skill 每次只喂一个拟合文件，故不涉及。
- 手册给 `datcmp` 的例子（`ly01.fir`，C=12、p=0.0445）是**短曲线**；曲线越长同样的 C 越显著，跨数据集比 C 的绝对值没有意义。

## 参考文件

`references/crysol-datcmp-parameters-and-fits.md` —— crysol 命令行逐项参数（默认值/单位/何时改）、`.fit`/`.fir`/`.dat`/`.out` 四种文件的结构、datcmp 三个检验的性质与输出字段、本机实测的三列 χ² 对照表、以及判档阈值的确切来源。

## 相关 skills

- **fit-a-high-resolution-model-to-data** — `contrasts-with`：RAW 界面里的判据（计算器选择、minimal 曲线陷阱、harmonics/N samples 症状—处方）。本 skill 做的是同一件事的命令行可复现版本 + 与已有拟合的比较。
- **evaluate-a-shape-reconstruction** — `contrasts-with`：那边评从头算重建本身；本 skill 把已有高分辨结构放进同一堆数据里比。**两者不一致时，文档口径是"以拟合结论为准"**（拟合好而珠模型不吻合 → 错的是珠模型）。
- **embed-a-model-in-bead-and-density-models** — `contrasts-with`：刚体叠合 ≠ 拟合；叠得好不等于配得上数据。
- **run-a-sec-saxs-pipeline-end-to-end / run-a-tube-saxs-pipeline-end-to-end** — `depends-on`：它们生产 `.dat` 与 DENSS/DAMMIF/IFT 产物，本 skill 消费这些产物。
- **assess-saxs-raw-data-quality / assess-guinier-fit-quality** — 前置：曲线不可信时先回那边。
- **analyze-saxs-data-with-atsas** — ATSAS 程序的索引（要换工具时先查它）。
