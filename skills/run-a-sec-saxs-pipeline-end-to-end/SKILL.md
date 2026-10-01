---
name: run-a-sec-saxs-pipeline-end-to-end
description: "端到端跑一条 SEC-SAXS 系列（图像→报告，全在 RAW 里做）：补出每帧 BL19U2 header txt 让 RAW 逐帧归一化（不写归一化 tif）、裁剪区归一化视频、buffer/sample 区与扣减、多区间 Guinier 结果、IFT、分子量、DAMMIF/DENSS 形状重建，每个节点都落 .dat + 表格 + RAW PDF 报告。用于「把这条 SEC-SAXS 数据端到端跑一遍」「一千多帧怎么变成曲线和报告」「归一化视频怎么做」「珠模怎么出」「要各节点 profile 的 dat」；不负责单点判据（Guinier 取点→assess-guinier-fit-quality，P(r)/Dmax→compute-and-validate-p-of-r，MW 方法选择→choose-a-molecular-weight-method，重建评估→evaluate-a-shape-reconstruction）。"
version: 1.0.0
author: hermes
license: MIT
tags: [saxs, sec-saxs, bioxtas-raw, pipeline, normalization, bl19u2, video, guinier, ift, dammif, denss]
metadata:
  hermes:
    tags: [saxs, sec-saxs, bioxtas-raw, pipeline, normalization, bl19u2, video, guinier, ift, dammif, denss]
    related_skills:
      - slug: process-sec-saxs-series
        relation: composes-with
      - slug: script-raw-with-the-python-api
        relation: composes-with
      - slug: correct-sec-saxs-baseline
        relation: composes-with
      - slug: configure-bioxtas-raw-for-a-dataset
        relation: composes-with
      - slug: fit-a-high-resolution-model-to-data
        relation: composes-with
---

# 端到端跑一条 SEC-SAXS 系列（全程 RAW）

这条任务是**串流程**，不是教你判据：输入一个目录（SEC 的连续洗脱 tif + 监视器 + 采集日志），
输出「归一化视频 + 一条扣减曲线 + 多区间拟合结果 + IFT + 形状模型 + RAW 报告 + 每个节点的 .dat」。
中间**所有计算都由 RAW 完成**——脚本只负责拼参数与落盘。用户明确的要求：**不自己写积分/拟合/P(r)/重建算法，只调参数。**

一句话定位：**它回答"这条 SEC 系列怎么从图像一路跑到报告、每一步留什么文件"，不回答"这段 q 该不该取、这个 Rg 能不能信"（后者转判据类 skill）。**

## When to Use（什么时候用）

- 一整个 SEC-SAXS 系列要端到端跑一遍，还要视频 / 拟合图 / 珠模 / 各节点数据。
- 用户问「这批 SEC 数据怎么处理完」「一千多帧怎么变成一条曲线」「归一化视频怎么做」「珠模怎么出」「各节点的 dat 都要」。
- 线站只给了 `.Iochamber` + `.log`（SEC 模式没有逐帧 header txt），而 RAW 恰恰要靠这个 txt 归一化。

**不负责**：单点判据。Rg 取点与可信度 → `assess-guinier-fit-quality`；P(r)/Dmax/GNOM 选择 →
`compute-and-validate-p-of-r`；分子量方法选择 → `choose-a-molecular-weight-method`；重建结果评估 →
`evaluate-a-shape-reconstruction`；基线该用 Linear 还是 Integral → `correct-sec-saxs-baseline`。

## 数据形状（BL19U2 的两种模式，差别就在 txt）

| 模式 | 每帧图像 | 配套文件 |
|---|---|---|
| batch / 管式 | `<系列>_<run>_<帧号>.tif` | 每帧一份 `<同名>.txt`（头里 `Transmitted_Beam`、`SR Current`、`Exposure time`）+ `<系列>_<run>.Intensity` + `_00001.log` |
| **SEC（连续洗脱）** | `<系列>_<帧号>.tif` | **只有** `<系列>_1.Iochamber`（监视器时间序列）+ `<系列>_00001.log`（每帧 endTime） |

两个必须记住的实测事实：

- **监视器文件的"行"和"帧"不是一一对应**：本机实测 `bsa_1.Iochamber` 19651 行 / 2000 帧 ≈ **9.83 采样/帧**
  （6.5 点/s × 1.5 s 曝光）。所以只能**按时间窗口取中位数**，不能按行号对帧。
- **监视器开头有几行 `~1e-13` 的野值**（未开束流），正常值是 **`2.6e-08`** 量级；这些行必须丢掉，
  否则窗口落进这段的帧因子会爆掉。

## 核心机制：RAW 自己就会用逐帧 txt 归一化（别自己写）

1. RAW 的 header 格式清单里就有 **`BL19U2, SSRF`**（`SASFileIO.py:909` `parseBL19U2HeaderFile`），
   它按 **`<图像路径去掉扩展名>.txt`** 读取——**txt 必须和 tif 同名同目录**。
2. 归一化由 settings 驱动（`SASImage.integrateCalibrateNormalize` → `calcExpression`）：
   - `ImageHdrFormat = 'BL19U2, SSRF'`
   - `EnableNormalization = True`
   - `NormalizationList = [['/', 'Transmitted_Beam']]` —— 把每帧**除以该帧 txt 里的 `Transmitted_Beam`**
     （`SASImage.py:366-380` 多乘性项合并成一个 factor，`:520-545` 对非乘性项逐项 `scaleRawIntensity`）。
3. 线站下机的 `.cfg` 里通常**已经**是这套（本机 `20261001.cfg`：`NormalizationList=[["/","Transmitted_Beam"]]`），
   只有 `ImageHdrFormat` 是 `None`、`EnableNormalization` 是 `false` —— 两个值改掉即可。
4. 实测（bsa 前 5 帧）：开启归一化后 I(q) 恰好是原来的 **1/TB**（比率 3.894e7 = 1/2.568e-8，逐帧对上）。

**推论：不需要写归一化后的 tif。** 逐帧 txt 是 KB 级，归一化在积分时发生——省掉 2000 帧 × 9 MB ≈ 18 GB。

## 执行步骤

### Step 0 — 前置核对

- 用**装了 RAW 的 python**（本机 `/Applications/BioXTASRAW/bin/python`），且**不要站在 RAW 源码目录里**跑
  （源码树会遮蔽 site-packages 里编译好的 `sascalc_exts`，直接 ImportError）。
- 核对该会话的 `.cfg`（定心/距离/掩膜/标样）→ `configure-bioxtas-raw-for-a-dataset`。

### Step 1 — 逐帧 header txt：**先看线站给了什么**，再决定生成还是不生成

**1a. 线站已给逐帧 txt（`tif` 旁边就有 `<帧名>.txt`）→ 只读因子，不重生成：**

```bash
python emit-bl19u2-header-txt.py \
  --series-dir <原始 tif 目录> --from-txt-dir <同一个目录> --out-dir <产物根>
```

（只把每份 txt 里的 `Transmitted_Beam` 汇总成 `<out>/norm/normalization_factors.csv` 给 Step 2 用；
RAW 那边本来就直接读这些 txt，**不需要 monitor/log，也不要 `--force` 覆盖线站原件**。
本机 `4LI2-676` 第二次下机就是这种：1800 份 txt、TB 中位 0.521、1–99% 展宽 4.6%。）

**1b. 只给了监视器 + 日志 → 生成逐帧 txt：**

```bash
python emit-bl19u2-header-txt.py \
  --series-dir <原始 tif 目录> --monitor <系列>_1.Iochamber --log <系列>_00001.log \
  --out-dir <产物根> [--lag-auto]
```

- 每帧的 `Transmitted_Beam` = 该帧曝光窗口 `[endTime − 曝光, endTime]`（含时钟偏移 lag）内**监视器中位数**；
- 低于 `--monitor-min-frac`（默认 5%）× 全局中位的采样点按"无束流"丢弃，个别窗口没采样则线性插值（会打印条数）；
- **txt 直接写进源数据目录**（与 tif 并排，RAW 才找得到）；因子表与元数据写 `<out>/norm/`；
- 线站原件保护：目录里已有同名 txt 时默认**拒绝覆盖**，确认要盖再加 `--force`；
- 时钟偏移：`--lag-auto` 用"检测器总计数 vs 监视器"的相关系数扫（本机 bsa 实测 **−32 帧 ≈ −48 s，corr 0.869**）。

完成标准：`<系列目录>/<帧名>.txt` 份数 = 帧数；`normalization_factors.csv` 里 TB 的 1–99 百分位展宽是合理的
通量起伏量级（本机 5%），没有离群帧（>20% 偏离要回头看 lag 与野值过滤）。

### Step 2 — 归一化裁剪视频（`crop-video-normalized.py`）

```bash
python crop-video-normalized.py --series-dir <tif 目录> --norm-csv <out>/norm/normalization_factors.csv \
  --out <out>/video/<系列>_crop_x<x1>-<x2>_y<y1>-<y2>_cols<c1>-<c2>_rows<r1>-<r2>_8x_20fps.mp4 \
  --x1 .. --x2 .. --y1 .. --y2 ..       # x/y 都从"大的那头"数（见下）
```

**坐标约定（别猜，用旧帧标定）**：BL19U2 视图给的 x/y **两个轴都从大的那头数** → 数组下标
`row = H-1-y`、`col = W-1-x`（`--x-origin right` 为默认；从左数才用 `--x-origin left`）。
标定法：拿一张**已被接受过的旧视频帧**，把候选读法各裁一份、按同一条渲染链渲出来算相关系数——
本机实测 `col=x` **0.11** vs `col=W-1-x` **0.66**，行向扫描峰值落在 `rows 708-768 = H-1-y`。

- 因子 `median(TB)/TB_i` 在**读入内存时**乘上，直接喂 ffmpeg —— 全程不落归一化 tif；
- 灰阶窗口用抽样帧定死后**逐帧不再自动拉伸**（自动拉伸会把要看的漂移抹平）；
- 帧号与因子烧在左上角；副产物 = 逐帧 `sum`/质心 CSV（束斑漂移监测）+（可选）`.npy` 堆栈 + 预览 PNG；
- 裁剪坐标约定见下（**两个轴都从"大的那头"数**：`row = H−1−y`、`col = W−1−x`）——别按"左下原点 (x,y)"直读。

完成标准：`ffmpeg -i out.mp4 -f null -` 的 `frame=` 等于帧数（或帧数/步长）；抽一帧看确实落在束斑/束挡区。

### Step 3 — RAW 端到端（`run-raw-sec-pipeline.py`）

```bash
python run-raw-sec-pipeline.py --series-dir <tif 目录> --out-dir <产物根> --cfg <日期>.cfg \
  [--steps integrate,series,guinier,ift,mw,shape,report] \
  [--buffer-range "s,e[;s,e]"] [--sample-range s,e] [--baseline none|linear|integral] \
  [--trim-qmin q] [--trim-qmax q] \
  [--guinier-ranges "qlo:qhi,..."] \
  [--model-engine auto|denss|dammif|both|none] [--denss-mode Fast|Slow|Custom] \
  [--n-models 4] [--symmetry P1] [--atsas-dir <ATSAS>/bin]
```

内部只调用 `bioxtasraw.RAWAPI`（= GUI 面板背后的同一套实现）：

| 节点 | RAW 入口 | 产物 |
|---|---|---|
| 积分（含逐帧归一化） | `load_and_integrate_images(files, settings)` | `profiles/01_integrated/<帧名>.dat` × N + `tables/frames_integrated.csv` |
| series + buffer 区 + 扣减 | `profiles_to_series` → `find_buffer_range` → `set_buffer_range` | `series/<前缀>_series.hdf5`、`profiles/02_buffer/buffer_avg.dat`、`profiles/03_subtracted/<帧名>_sub.dat` × N |
| （可选）基线 | `find_baseline_range` / `set_baseline_correction` | `profiles/05_baseline/*.dat`；之后 `profile_type='baseline'` |
| 样品区 + 样品平均 | `find_sample_range` → `set_sample_range` | `profiles/04_sample/sample_avg.dat`、`series/ranges.json` |
| 逐帧参数 | `series.getFrames/getRg/getI0/getVcMW/getVpMW/getIntI` | `tables/frame_params.csv`、`series/series_plot.png` |
| 多区间 Guinier | `auto_guinier` + `guinier_fit(idx_min, idx_max)` | `profiles/06_guinier/guinier_<标签>.dat`、`tables/guinier_multi_range.csv` + `.png` |
| IFT | `bift`（原生）/ `auto_dmax`+`gnom`（需 ATSAS） | `ifts/<前缀>_bift.ift` / `_gnom.out`、`tables/ift_summary.csv` |
| 分子量 | `mw_vc` / `mw_vp` / `mw_bayes` | `tables/mw.csv` |
| 形状重建 | `denss`（**电子云，RAW 原生，默认总跑**）+ `dammif`/`damaver`（**珠模，有 ATSAS 才跑**） | `models/<前缀>_denss.mrc`（+`_support.mrc`/`_map.fit`/`_stats_by_step.dat`/`_denss.log`）；有 ATSAS 时另有 `models/<前缀>_dammif_0N-1.cif`（**ATSAS≥4.0 写 `.cif`，不是 `.pdb`**）+ DAMAVER 的 `<前缀>_damaver-{distances.txt,global-summary.txt,global-fsc.dat,cluster*-summary.txt,global-damaver.cif}` |
| 报告 | `save_report(pdf, dir, profiles, ifts, series)` | `reports/<前缀>_raw_report.pdf` |
| **结果说明（自动，最后一步）** | 纯读产物（`results_readme.py`，不重算） | **`<产物根>/README.md`** —— 给人看的目录导航 + 关键数字 + 判读红线 + 本次告警 |

**形状重建两个都要**（`--model-engine auto`）：**电子云（DENSS，`.mrc`）与珠模（DAMMIF，`.pdb`）是同一份 IFT 的两种重建，互不替代**——
DENSS 是 RAW 原生（numba），`Fast` 模式实测几秒到几十秒就出；DAMMIF 是 ATSAS 可执行文件的外壳，**没装 ATSAS 时自动只出电子云**（并在日志里明说），不要因此把 `models/` 留空。

**多区间 Guinier 是"让用户复核拟合过程"的主产物**：不传 `--guinier-ranges` 时，脚本以 `auto_guinier` 的 Rg 为锚
铺一条跨判据边界的阶梯（qRg 0.3–0.6 / 0.4–0.8 / 0.5–1.0 / 0.6–1.3），每个区间**一份 profile 副本**
（`setQrange` 截到该区间）单独存 `.dat` 并进 RAW 报告，同框图用 **RAW 返回的 Rg/I0** 画模型线（不是自写拟合）。
表里给全 `Rg / I0 / rg_err / i0_err / q_min / q_max / qRg_min / qRg_max / r²`，据此判断"哪段才合适"。

### Step 4 — 结果目录里的 `README.md`（**别再让人对着文件夹发呆**）

管线**最后一步自动**写 `<产物根>/README.md`：只读产物（`run_meta.json` / `tables/*.csv` /
`series/ranges.json` / `models/` 清单），**不重新拟合、不臆造数字**，每一节都能指回某个文件。

**开头必须是"最终拟合参数 + 误差 + 可靠性 + 误差怎么读"三张表**（用户明确要求：先看结论，且要表格）：

| 表 | 内容 | 关键写法 |
|---|---|---|
| **表 1 最终主拟合** | Rg ± `rg_err`、I(0) ± `i0_err`、拟合的 q 区间与 qRg 区间、`r²`、**跨区间一致性**、可靠性列 | 「采用区间」按 **收敛 + qRg 落在 0.3–1.35 + r² 最高**挑；一致性只统计**可用区间**（`rg>0` 且 `r²≥0.9`），并把被排除的区间**列出来说明为什么排除**（一票否决的错误做法：把 `rg=-1` 的失败区间也算进一致性） |
| **表 2 其它参数** | 每支 IFT 的 Dmax / Rg ± `rg_err`（含 `Dmax/Rg` 比与判词）、MW 四法（Vc 带 ± ，Bayesian 带 CI）、DENSS 的 χ²/Rg/support、DAMMIF×N 的 χ²/Rg/Dmax/MW、DAMAVER 平均 NSD ± std | 每行**必须给可靠性**：🟢 可用 / 🟡 有限 / 🔴 不可用，并写一句"为什么"（χ² ≤3 可靠、3–5 有限、>5 先修曲线；Dmax/Rg ≤3.2 正常、>4.5 判不可用；NSD ≤2 同一类形状） |
| **表 3 误差怎么读** | 每条误差**是什么 / 不包含什么 / 什么时候不能用** | 清楚区分**随机误差**（`rg_err` 只是该区间内最小二乘标准误差，区间收窄它必然变小）与**系统误差**（低 q 污染、背景失配、Dmax 选择、浓度假设）；MW 的 ± 只来自 Vc 一种方法；DAMAVER 的 NSD 是模型间差异、不是与数据的拟合误差 |

之后才是：输入与处理概要 → 只看三个文件 → 目录导航 → 关键数字明细 → 复现命令 → 判读红线 → 本次告警 → 缩写。

给**已跑完但没 README**的旧结果目录补写（不动数据）：

```bash
python write-results-readme.py --out <产物根>
```

## 复核点（跑完先看这几处，再看数字）

1. `tables/guinier_multi_range.csv`：区间之间 Rg 是否稳定；`qRg_max` 是否越过形状对应的上界（球≈1.3）；
   `r²` 与 `rg_err` 是否随区间缩窄而恶化。判据细节 → `assess-guinier-fit-quality`。
2. `series/series_plot.png`：峰上 `Rg/I0` 是否成平台（不成平台说明多组分/聚集）；
   buffer/sample 两条阴影就是脚本用的区间。
3. `tables/ift_summary.csv`：BIFT 的 `chi²`、`log_alpha`、`evidence`；GNOM 的 `dmax` 是自动定的，要复核。
4. `tables/mw.csv`：Vc/Vp（浓度无关）是否互相一致；SEC 峰内浓度未知，**不要**用 I0 标样/绝对刻度那几法
   → `choose-a-molecular-weight-method`。
5. `models/`：DAMMIF 至少 4 个模型看 a-score/平均 NSD；DENSS 看 FSC 分辨率 → `evaluate-a-shape-reconstruction`。

## 坑（都是实测撞出来的）

- **在 RAW 源码目录里跑会 ImportError `bioxtasraw.sascalc_exts`**：源码树遮蔽 site-packages 里编译好的包。
  要么换目录跑，要么在源码树里 `build_ext --inplace`。
- **`find_buffer_range` / `find_sample_range` 在"截断的系列"上会失败**（返回 `success=False`、区间为 `None`）：
  只取峰前的帧、或系列里没有可识别洗脱峰时必然如此 → 用 `--buffer-range/--sample-range` 手工给（0 基帧号）。
- **ATSAS 不在 = DAMMIF/GNOM/DATGNOM/CIFSUP 全不可用**：RAW 的 DAMMIF 是**外壳**，
  `RAW.py:1232-1242` 拼的是 `os.path.join(atsas_dir, 'dammif')`，`RAWAPI.dammif` 会抛 `NoATSASError`。
  没有 ATSAS 时跑 RAW 原生 `denss` + `bift`（都可用），装好 ATSAS（`https://biosaxs.com/download`，学术免费，
  需个性化 license）后加 `--atsas-dir <ATSAS>/bin` 即切到 DAMMIF/GNOM。
- **逐帧 txt 有三种来源，先看清楚再动手**：① 线站已经给了（本机 `4LI2-676` 第二次下机就带了 1800 份
  `<帧名>.txt`，和 tif 并排）→ **不要重新生成**，用 `--from-txt-dir <目录>` 只把它们的 `Transmitted_Beam`
  读成因子表（给视频用），RAW 那边本来就直接读这些 txt；② 只给了监视器 + 日志 → 用 Step 1 生成；
  ③ 两样都没有（本机 `4EH2-KDPV-ZN`：445 个 tif，无 `.Iochamber`/`.log`/`.txt`）→ 只能
  `--no-header-normalization` 跑（RAW `ImageHdrFormat=None`），视频省略 `--norm-csv`（因子全 1），
  **并在报告里写明"本系列未做通量归一化"**，别让读者以为做了。
- **逐帧 txt 写进源数据目录**是刻意的（RAW 只在同目录找 `<帧名>.txt`）；写之前目录里若已有同名 txt 会被拒绝，
  真要覆盖才 `--force`——线站原件优先。
- **`--buffer-range` 优先给峰前 + 峰后两段**（`"560,620;800,880"`）：SEC 峰后的基线常不等于峰前
  （窗口污染 / 束位漂移）。本机实测同一条 BSA 曲线，只用峰前一段时 Guinier 的 Rg 在 38–117 Å 之间随区间乱跳；
  换成峰前+峰后两段后落到 **23.7–27.4 Å、r² 0.91–0.99**（BSA 单体理论 ≈29 Å）。
- **`--trim-qmax`：判断"高 q 的失配是模型表达不出来，还是数据本身不好"**。同一条曲线把 q 上界压到
  0.25 / 0.15 再跑一遍，看 χ² 往哪边走：**本例 4LI2-676**（Rg 15.3 Å 的小蛋白，珠模 χ²=4.66）——
  IFT 的 χ² 反而**升**（1.98 → 2.6 → 3.1），说明失配**不是**高 q 独有，低/中 q 段本来就有平滑 P(r) 表达不了的结构；
  而珠模的 Σχ² 有 66% 来自 q 0.2–0.46 段（那里 σ 小、点数多，系统性偏 ~2σ）。两件事一起读才对：
  **"点多、误差小的高 q 段"会把模型表达能力的不足放大成 χ²**，但**不能**用它来解释整条曲线的失配。
  裁高 q 后 `auto_dmax` 常返回 −1 让 GNOM 报 `expected value ≥ 0 for option 'rmax', got '-1'` → 这时配 `--ift-dmax` 显式给 Dmax 即可。
- **高 q 该截到哪，由 I/σ 决定，别凭感觉也别人云亦云**：读 `profiles/04_sample/sample_avg.dat` 的误差列，
  按 q 分箱算中位 I/σ —— **中位 I/σ 掉到 ~2 以下的那个 q 就是 IFT/珠模该截的地方**；若高 q 段不只是弱、
  还**系统性为负**（负值比例过半），那叫**扣减过头**，跟"高 q 弱"是两回事，但一样得截。本机实测
  **4EH2-KDPV-ZN**（1500 帧、有监视器归一，色谱正常）：q ≤ 0.15 中位 I/σ 17–80，0.15–0.25 掉到 **0.91**
  （46% 负值），0.25–0.45 中位 **−4.28、98% 负值**。整段喂 IFT 的后果是**三条一起崩**：DAMMIF Slow
  χ² **15.03**（其它样品 1.1–4.9）、DENSS 撑不出支撑直接失败（`index -1 is out of bounds … size 0`）、
  MW **Vp 116 kDa**（同一批的贝叶斯/datclass 都是 ~29）——而且**换 Dmax 没用**（`--ift-dmax 79` 后 χ² 仍 15.03），
  说明病根在数据段不在 Dmax。换 `--trim-qmax 0.15` 后：DAMMIF χ² **1.40**、DENSS χ² 1.42 复活、
  Vp 33.8、DAMAVER NSD 0.590±0.070。同批对照 **4DH2-676-apo-3** 的信噪截点在 **0.217**（χ² 4.25）——
  **每条曲线各量一次**，不要照抄别的样品的数字。代价：高 q 截短后 DATCLASS 会以
  `insufficient data to integrate to s*Rg >= 5.0` 跳过、BIFT 的自动 Dmax 可能退化（实测 400 Å）→ 用 `--ift-dmax` 显式给。
- **低 q 被寄生散射污染会把 IFT 的 Dmax 拖到离谱值**：同一曲线未裁 q 时 BIFT 给 **Dmax=417 Å / Rg=149 Å**。
  用 `--trim-qmin 0.017` 显式裁掉低 q 段（内部就是 RAW 自己的 `setQrange`，不是自写拟合），再跑 IFT/MW。
- **ATSAS 接线三件事**（装好 ATSAS 后要一次对上）：① `--atsas-dir` 必须指到 **`bin` 这一级**
  （RAW 用 `os.path.split(atsas_dir)[0]` 反推 `ATSAS` 变量；本机实测 `/Applications/ATSAS-4.1.4-1/bin`）；
  ② `dammif` 写出的模型名是 **`<prefix>-1.<model_format>`**（默认 `cif`），而 `damaver` **只接受文件名、
  不接受路径**、且 `model_format` 必须与 dammif 一致 —— 对不上就报
  `FileNotFoundError: <prefix>-damaver-distances.txt`（看着像 DAMAVER 坏了，其实是输入清单错了）；
  ③ `mw_bayes` / `mw_datclass` **没有 `settings` 参数**（签名是 `profile, rg, i0, first, atsas_dir, ...`），
  传 `atsas_dir` 才生效，否则一直报 TypeError 被误当"没装 ATSAS"。
- **`damaver` 的返回值里第 2 个是标准差，不是均值**：`NSD=0.09` 这种日志看着很美，其实那是 stdev ——
  均值要读 `models/<前缀>_damaver-distances.txt` 的 `Mean value of nsd`（代表模型在
  `-global-summary.txt` 的 `Most representative`）。**别把 stdev 当一致性指标**。
- **Fast/Slow 与 NSD 没有单调关系，判质量还是看 χ² + 聚类数**：本机同一条 bsa，`Fast` 4 模型
  `χ²=2.08、平均 NSD 0.77±0.06`（2 簇）、`Slow` 4 模型 `χ²=2.08、平均 NSD 1.11±0.09`（2 簇）——
  Fast 反而"更雷同"（都收敛到相似的模糊解）。另两条 `Slow`：4LI2-676 `1.16±0.20`（3 簇）、
  4DH2-676-apo-3 `1.17±0.29`（3 簇）。代价：`Fast` 8–40 s/模型，`Slow` 约 3 min/模型（4 个 ≈12 min）——
  先用 Fast 看 χ² 能否建模，交付前用 Slow 复核。
- **χ² 大 = 这组数据还不配做从头建模**，别在重建参数里找答案。三条实测（同一套流程、同一天）：
  bsa `χ²=2.08 / Rg 28.0 Å / Dmax≈87 Å / MW 61 kDa / DAMAVER NSD 0.09`（可用）；4LI2-676
  `χ²=4.87 / Rg 16.3 Å / Dmax≈53 Å / MW 9 kDa / NSD 0.20`（可用）；4DH2-676-apo-3
  `χ²=4.27 / Rg 25.3 Å / Dmax≈96 Å / MW 28–30 kDa / NSD 0.29`（尚可，模型偏扁平，需与高分辨模型/其他证据对照）。
- **DENSS 报 `ValueError: The number of derivatives at boundaries does not match: expected 3, got 0+0`**（来自
  `DENSS.py:4269 regrid_Iq` 的三次样条）= **扣减后曲线已退化**（点数/单调性不足、噪声主导），本机 `4EH2-KDPV-ZN`
  就是这个报错。这一条与下面的 IFT 问题是一类：**先判数据，再谈重建**。
- **DENSS 报 `DENSS failed to run properly` / `IndexError: index -1 is out of bounds ... labeled_support == feature` = 它的输入 IFT 不可用**，不是 DENSS 参数没调好：那个下标来自 `DENSS.py:1931` —— 收缩包络把 support 压成**空**（`num_features == 0`）时 `sums` 长度为 0。
  根因几乎总是 IFT 被低 q 拖出**离谱的 Dmax**（本机 4LI2-676 实测：Guinier Rg=15.9 Å，BIFT 却给 Dmax=189–357 Å，
  DENSS 盒子 side>1000 Å、密度摊薄 → 塌陷）。处置顺序：① 先把 q 裁到 qRg_min ≈ 0.4–0.5（`--trim-qmin`）；
  ② 仍不对就用**显式 Dmax** 走 RAW 原生 DIFT 喂 DENSS（`--ift-dmax 55`，经验起手 Dmax ≈ 3×Guinier Rg）；
  ③ 只在 ① ② 之后才考虑 DENSS 自己的步骤/盒子参数。**不要在 IFT 坏的时候去调 DENSS 参数。**
- **`models/` 空着 = 没跑 `shape` 步，或缺 ATSAS 又被当成"珠模出不来就什么都不出"**：正确姿势是
  `--model-engine auto` —— 电子云（DENSS）总出，珠模（DAMMIF）有 ATSAS 才出；两者都缺才叫失败。
- **帧区间是 0 基闭区间，`end` 给到帧数就会 `IndexError`**：`445` 帧的系列写 `--buffer-range "60,220;400,445"`
  → `SECM.averageFrames` 里 `list index out of range`（本机实测）。脚本现在有 `clip_ranges` 自动收到 `n−1` 并告警，
  但**区间里的数字仍应自己核对**（峰位置看副产物色谱图最保险）。
- **判"这条系列能不能用"的硬指标**（任一命中就是数据问题，不是参数问题）：Guinier `Rg < 5 Å` 或 `r² < 0`；
  `Vp` 为负；逐帧 `Rg` 在几十到几埃之间乱跳；DENSS 建不起样条。本机 `4EH2-KDPV-ZN`（445 帧、无监视器、
  洗脱事件只抬高 ~3.7%、两种区间都给出 Rg 1.27–1.29 Å）就是这种：**在产物目录写一份
  `结果不可用-README.md` 说明判定依据**（数据事实 + 试过的区间对照 + 三重证据 + 建议），别让表格里的数字被当结果引用。
- **脚本被并发编辑时，先 `git status` 再跑，必要时用快照跑**：本机同一仓库有别的会话在改同一个脚本，
  撞到过两次半成品（`step_ift` 改返回 5 值而 `main` 还解 4；`step_shape` 里 `itf`/`ift` 拼写），
  两次都是"跑到一半 ValueError/NameError"。稳妥做法：`git log -1` 记下 revision → `cp` 到 scratch 当快照 →
  跑快照 → 报告里写清用的哪个 sha。
- **"结果看不懂"是缺交付物，不是缺解释**：产物目录必须有一份 `<产物根>/README.md`（管线最后一步自动写，
  旧目录用 `write-results-readme.py --out <目录>` 补）。它只读产物、不重算，所以写完要和 `run_meta.json`
  的时间戳对一下（避免"README 比数据旧"）。数据本身不可用的系列，另写 `结果不可用-README.md`（见下条）。
- **表里 `rg = -1` 不是负数，是 RAW 的失败哨兵值**（该区间的 Guinier 拟合没收敛/点数不够）——按"此区间不可用"读，
  不要当成数值；同理 `r² < 0` 表示拟合比取平均还差。
- **裁剪坐标：两个轴都从"大的那头"数** → `row = H−1−y`、**`col = W−1−x`**（即 180° 旋转；本机 BL19U2
  实测标定：`col=x` 时与旧视频帧相关 0.11，`col=W−1−x` 时 0.66，行向扫描峰值恰为 `H−1−y`）。脚本
  `crop-video-normalized.py` 默认 `--x-origin right`。裁出来一片均匀背景/位置偏移说明约定错了，别微调数字，
  回 `frame-sequence-to-video` 的"用已接受的旧帧算相关系数"重标法。
- **ffmpeg 的 libx264 + yuv420p 要求偶数边长**：61×71 的裁块靠整数放大（8×）顺带解决。
- 视频的 `-pix_fmt rgb24` 必须与写进管道的字节一致（写 RGB 就声明 rgb24），否则帧数会变 3 倍。

## 相关 skills

- **process-sec-saxs-series** — 本 skill 的"判据版"：区间怎么选、峰上 Rg 平台怎么看、SEC 为什么不能用绝对刻度。
- **script-raw-with-the-python-api** — RAWAPI 的骨架与返回元组顺序（本 skill 的脚本就是它的落点）。
- **correct-sec-saxs-baseline** — 扣减后仍漂移时选 Linear/Integral，以及本 skill 用的监视器归一套路。
- **configure-bioxtas-raw-for-a-dataset** — 换实验日/仪器时先核对 `.cfg`。
- **evaluate-a-shape-reconstruction** / **fit-a-high-resolution-model-to-data** — 重建与高分辨模型对照。
