# 管式数据的 control 配对、归一化与相对缩放（实测）

本文件是 `run-a-tube-saxs-pipeline-end-to-end` 的机制附录：三条"不做就会静默出错"的设置，以及各自的实测数字。
数据来源：`~/Repositories/DataProcess_2026.10.01/data/`（BL19U2，2026-10-01，Pilatus 2M，D=2680 mm，λ=1.0332 Å，
`20261001.cfg`），RAW 2.4.2，python `/Applications/BioXTASRAW/bin/python`。

## 1. 目录即实验记录：control 就在同一个目录里

BL19U2 的 batch（管式）模式，线站把**样品 run 和夹着它的空白 run 一起放进一个目录**：

| 目录 | 里面的 run | control 前缀 |
|---|---|---|
| `Tube-SAXS/A5-05-1` … `-6` | `ddh2o_00XX` + `A5-05-N_00YY`（+ 后一个 `ddh2o`） | `ddh2o` |
| `Tube-SAXS/BSA` | `pb7_0001` + `BSA_0002` + `pb7_0003` | `pb7` |
| `Tube-SAXS/5705` / `97df` / `877-apo-pb7` / `877-4zinc-pb7` | `pb7_00XX` 夹着样品 run | `pb7` |

采集顺序（按 run 号，`batch/` 汇总目录里可复原）：`pb7_0001 → BSA_0002 → pb7_0003 → 877-apo_0004 → pb7_0005 →
877-4zinc_0006 → pb7_0007 → 5705_0008 → pb7_0009 → 97df_0010 → pb7_0011 → ddh2o_0012 → A5-05-1_0013 → ddh2o_0014 →
A5-05-2_0015 → … → ddh2o_0022 → A5-05-6_0023`。即：**每组样品用同一种 blank 前后夹着**（前一组是 pb7，A5-05 那组是 ddH2O）。
判别法：`.txt` 里的 `Description` 字段 + run 号顺序。

脚本默认规则：目录名 = 样品前缀 → 同前缀帧是样品，其余全是 control；多前缀时 `--control-key ddh2o,pb7`。
一个目录里 60 帧通常是 3 个 run × 20 帧（control + 样品 + control）；`A5-05-6` 这类首尾样品只有一侧背景（40 帧 = 2 run）。

**两条来自上游 `organize-batch-saxs-dataset` 的判据，别自己发明**：

- 样品身份取**帧 1 的 `Description` 字段**，不是文件名前缀的子串——`877-apo-pb7` 是"在 pb7 里测的 877-apo 样品"，
  用子串判背景会把整个样品判成背景。
- 背景身份 = `Description` **精确等于**用户点名的背景名（典型 `pb7` / `ddh2o`）；目录还没归类时先跑那个 skill 拆目录。

## 2. 归一化：RAW 会读逐帧 txt，但 cfg 里两个开关是关着的

- `SASFileIO.parseBL19U2HeaderFile`（header 格式表里的 `'BL19U2, SSRF'`）按 **`<图像路径去扩展名>.txt`** 读，
  所以 txt 必须与 tif **同名同目录**（batch 模式线站本来就给了）。
- 归一化发生在积分时（`SASImage.integrateCalibrateNormalize` → `calcExpression`），由 settings 驱动：
  `ImageHdrFormat='BL19U2, SSRF'`、`EnableNormalization=True`、`NormalizationList=[['/','Transmitted_Beam']]`。
- 下机 `20261001.cfg` 的实测值：`ImageHdrFormat='None'`、`EnableNormalization=False`、`NormalizationList=[['/','Transmitted_Beam']]`
  —— 也就是说**归一化列表在，但开关是关的**，结果静默不归一化。
- 打开后的实测（A5-05-6 第一帧）：I(q=0.01) 从 181.7 → 403.5，比值 = 1/0.4447（该帧 `Transmitted_Beam`），逐帧精确对上。

推论：**不需要写归一化后的 tif**（几千帧 × 9 MB），归一化在积分时完成。

## 3. control 的相对缩放：1–3% 的高 q 常数残差

透射归一化之后，样品与 control 之间仍会差 ~2%（采集时间不同 → 束流/管位不同）。这个差值在扣减后留下一个
**近常数的高 q 残差**，症状是 `q²I` 在高 q 上翘；它会直接毁掉 IFT（BIFT 的 chisq 从 ~1 涨到几十）。

做法：因子 = `mean(I_sample(q) / I_control(q))`，取**粒子不散射的高 q 窗**（默认 q 0.30–0.44 Å⁻¹），
再用 RAW 自己的 `SASM.scaleRelative(f)` 作用在 control 上，然后 `subtract`。

本机 `A5-05-1` 实测：

| 取样窗口 | 因子 | 说明 |
|---|---|---|
| q 0.15–0.25 | **1.0937** | 窗太低：把蛋白自身的信号当背景了，明显偏大 |
| q 0.25–0.35 | 1.0200 | |
| q 0.30–0.44 | **1.0236** | 采用值 |
| q 0.35–0.44 | 1.0263 | 与 0.30–0.44 一致 → 说明窗口取对了（不一致就该往高 q 挪） |

效果：`auto_guinier` 的 R² 0.873 → **0.960**；扣减曲线高 q 归零。
判据：**同一份数据在两个相邻高 q 窗里给出的因子应彼此一致（差 <1%）**；因子偏离 1 超过 5% 时不要硬缩放，
先怀疑 control 配错了（换 run / 换 buffer）。

## 4. 多区间 Guinier 的闸门与列

脚本对每个区间存一份 `setQrange` 截断的 profile 副本，并用**RAW 返回的 Rg/I0** 画拟合线；
自己只算残差统计。`tables/guinier_multi_range.csv` 的列：

| 列 | 含义 / 闸门 |
|---|---|
| `n_pts` | 区间点数，闸门 ≥10（取点太少的"好 R²"没有意义） |
| `qmin` / `qmax` | 区间边界（报告里必须带） |
| `Rg` / `Rg_err` / `I0` | RAW `guinier_fit` 的返回值 |
| `qRg_min` / `qRg_max` | 闸门 `qRg_max ≤ 1.3`（球；棒 1.0 / 盘 1.7） |
| `R2` | 线性化的判定系数（仅供参考） |
| `chi2_red` | 误差加权残差的约化 χ²（用 RAW 的 Rg/I0 反推直线），闸门 ≤3，**看它是否≈1** |
| `curvature` | 前半段残差均值 − 后半段均值，>0 = smile（聚集/残留背景），<0 = frown（排斥） |
| `rg_split_rel` | 前后半段各自拟合的 Rg 相对差，闸门 ≤0.05（区间内是否真是一条直线） |
| `pass_*` | 四条闸门的逐条结果 |

本机 `A5-05-1` 实测（17 个区间，全表见产物）。**这是一个"没有可信 Rg"的样例**，正好说明为什么要多区间：

| q 区间 | Rg | qRg_max | R² | chi2_red |
|---|---|---|---|---|
| 0.0064–0.0465 | 11.4 | 0.53 | 0.608 | 12.8 |
| 0.0064–0.0754 | 17.2 | 1.30 | 0.960 | 130 |
| 0.0314–0.0695 | 17.8 | 1.24 | 0.982 | 49 |
| 0.0567–0.0754 | 19.5 | 1.47 | 0.9995 | **0.91** |

Rg 从 11.4 漂到 19.5 Å：越往高 q 挪越"像直线"（χ²→1），但那时 q·Rg 已 >1，早出了 Guinier 区。
结论写法：**"多区间扫描显示 Rg 随区间单调漂移 10.8–19.6 Å，没有任一区间同时满足 qRg≤1.3 与 χ²≈1 →
该曲线不服从单一 Guinier 定律（多分散 / 残留背景 / control 不匹配），不报 Rg"**。判据细节 → `assess-guinier-fit-quality`。

## 5. BIFT 的 Dmax 搜索域

RAW 默认 `minDmax/maxDmax = 10/400 Å`、`DmaxPoints=10`、`PrPoints=100`。对 Rg≈17 Å 的样品，
默认域会让 BIFT 在 400 Å 交出一个 chisq 20–70 的假解（不是报错，是安静地给一个错的 Dmax）。

脚本按 `0.7–1.6 × 3.1·Rg` 收窄（球状粒子的 Dmax≈3.1 Rg 起手），实测 `A5-05-1`：
收窄前 `Dmax=400.1 Å / Rg(realspace)=127.7 Å / chisq=32.8`；收窄后 `Dmax=48.4 Å / Rg=18.3 Å / chisq 仍高`。
`--ift-sweep N` 会扫 Dmax 并落 `tables/ift_dmax_sweep.csv`（看 chisq 谷底的平台）。

**低 q 的可用起点**：beamstop 挡到 r≤12 px，`q=0.0064` 对应的 r≈16 px 处只有 ~38% 的环未遮挡（完全未遮挡 ≈90%），
r≈25 px（q≈0.010）才到 ~75%。第一个数据点因此坐在光束挡边缘，携带一段扣不掉的边缘光晕——实测 `A5-05-5` 的首点 I 从 130 降到 37（砍到 q≥0.010）。
要做低 q 判断（上翘、Guinier、IFT）先把它砍掉：`--qmin 0.010`。

**这个域不是硬边界**：RAWAPI 原文写着 *"The value of Dmax can go beyond this bound in the optimization step"*。
实测 `A5-05-6`（最稀那条）在 48–109 Å 的网格下交出 **Dmax=741.4 Å、Rg(realspace)=269 Å、chisq=1.09**，
而把 IFT 起点从"分析窗起点"换成"Guinier 拟合自己的起点"（idx 78，q≥0.0368）后是 **Dmax=88.9 Å、Rg=25.3 Å、chisq=1.23**。
两条都"收敛"，只有闸门能分辨：本 skill 因此**跑两次 BIFT**并落 `tables/ift_summary.csv` 两行
（列：`tag / qmin / qmax / Dmax / Rg_realspace / chisq / rg_tol / rg_ref_soft / pass_rg / pass_dmax_over_rg / pass_dmax_within_grid / trusted`），
不可信的那条不进 3D、也不占用 `ifts/bift.ift` 这个名字。
注意 `chisq` 在低 q 信噪比很高的管式数据上会被系统误差主导（几十也可能代表"拟合形状是对的、误差被低估"），
要**和实空间 Rg 对 Guinier Rg 的一致性一起看**，不能只看一个数。

## 6. ATSAS 边界（RAW 里"珠模"那一支全是外壳）

- RAW 包里自带实现的求解器只有：`BIFT.py`、`DENSS.py`（带 `denss_resources/`）、`REGALS.py`。
- `RAWAPI.dammif` 的文档字符串：*"Creates a bead model (dummy atom) reconstruction using DAMMIF **from the ATSAS package** …
  Will raise `SASExceptions.NoATSASError` if ATSAS is not found."*，实现里把 `atsas_dir` 交给 `SASCalc`，最终执行
  `<ATSASDir>/dammif` 二进制。`gnom/dammin/datgnom/datmw/datclass/cifsup/crysol` 同理。
- GUI 闸门顺序（本机源码 RAW.py，v2.4.2）：`showDAMMIFFrame` 在 **1209 行**先判 `iftm.getParameter('algorithm') != 'GNOM'`
  → 弹 **"Wrong IFT type"**（*"DAMMIF can only process IFTs produced by GNOM. This was produced using BIFT."*）；
  到 **1289 行**才判 ATSAS 是否存在 → 弹 **"Can't find ATSAS"**。所以用 BIFT 的 IFT 去点 DAMMIF，**验不到 ATSAS 有没有**；
  最快的判定点是同一个 ATSAS 子菜单里的 **GNOM**（1063 行同一条消息）。
- 另外 DAMMIF 需要的是 **GNOM 的 `.out`**：`RAWAPI.dammif` 会先 `SASFileIO.writeOutFile(ift, …)`，
  对 BIFT 的 IFTM 这一步就返回 `None`（BIFT 没有 GNOM 的参数集），链条在更早一步就断。
- 结论：无 ATSAS 时本机可行组合 = **BIFT（P(r)/Dmax）+ DENSS（3D 电子密度）+ Vp/Vc（MW）**；
  装好 ATSAS（`https://biosaxs.com/download`，学术免费、需个性化 license）后用 `--atsas-dir <ATSAS>/bin`
  切到 GNOM/DAMMIF/DATMW/DATCLASS。

## 7. 本机实测一览（供回归对比）

`A5-05-x` = 一条 2 倍稀释序列（`Tube-SAXS/A5-05-N`，control = 夹着的 `ddh2o`）；数值为脚本 `summary.json` 的
`guinier_auto` 与 `mw` 节点：

| 样品 | 帧数 | I0 (auto-Guinier) | Rg (auto) | MW Vp / Vc (kDa) |
|---|---|---|---|---|
| A5-05-1 | 20 + 40 | 65.5 | 17.3 | 27.9 / 20.9 |
| A5-05-2 | 20 + 40 | 35.8 | 19.2 | |
| A5-05-3 | 20 + 40 | 20.2 | 21.7 | |
| A5-05-4 | 20 + 40 | 10.9 | 24.2 | |
| A5-05-5 | 20 + 40 | 5.3 | 23.8 | |
| A5-05-6 | 20 + 20 | 2.4 | 23.5 | |
| BSA | 20 + 40 | (Rg 79 Å、IFT 74 Å) | — | 114 / 75 |

`BSA` 那条用 pb7 做 control 时，Guinier 与 IFT 一致地给 Rg≈74–79 Å、Dmax≈460 Å，远大于单体 BSA（29 Å / 66 kDa）：
**先怀疑 control 不匹配或样品聚集，不要当成"拟合方法不对"**（两套独立方法给同一个数，说明数是从数据里来的）。
