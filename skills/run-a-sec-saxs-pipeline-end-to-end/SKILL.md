---
name: run-a-sec-saxs-pipeline-end-to-end
description: "端到端跑一条 SEC-SAXS 系列：多峰逐峰分析、逐帧归一化、裁剪区视频、多区间 Guinier 表、IFT、分子量、DAMMIF/DENSS 形状重建（图像→报告全程 RAW；无逐帧 txt 时补 BL19U2 header txt），每个节点落 .dat + 表格 + 图 + RAW PDF 报告 + 结果 README。用于「把这条 SEC-SAXS 数据端到端跑一遍」「这条系列有两个洗脱峰、怎么分别分析」「一千多帧怎么变成曲线和报告」「归一化视频怎么做」；不负责单点判据（Guinier 取点→assess-guinier-fit-quality，P(r)/Dmax→compute-and-validate-p-of-r，MW 方法选择→choose-a-molecular-weight-method，重建评估→evaluate-a-shape-reconstruction，未解析重叠峰分解→deconvolve-overlapping-elution-peaks）。"
version: 1.3.0
author: hermes
license: MIT
tags: [saxs, sec-saxs, bioxtas-raw, pipeline, normalization, bl19u2, video, guinier, ift, dammif, denss, multi-peak, peak-detection]
metadata:
  hermes:
    tags: [saxs, sec-saxs, bioxtas-raw, pipeline, normalization, bl19u2, video, guinier, ift, dammif, denss, multi-peak, peak-detection]
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
      - slug: write-saxs-results-readme
        relation: composes-with
      - slug: deconvolve-overlapping-elution-peaks
        relation: contrasts-with
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

**同名样品可能有两个下机批次，先数帧数与配套文件再决定走哪条路**：本机 `4EH2-KDPV-ZN` 既有**445 帧、无 `.txt`/`.Iochamber`/`.log`** 的旧批次（属于下表第 ③ 种，只能 `--no-header-normalization`），也有 **1500 帧 + 1500 份线站逐帧 txt + 1 份 `.Iochamber` + 1 份 `.log`** 的新批次（第 ① 种，**不要**去重生成 txt、更不要用 `--no-header-normalization` 白丢逐帧归一化）。

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
  [--steps integrate,peaks,series,guinier,ift,mw,shape,report] \
  [--multi-peak auto|off|always] [--peak-ranges "lo,hi;lo,hi"] [--peak-buffers "s,e;s,e|s,e"] \
  [--peak-q-range "qlo,qhi"] [--frame-flag-q "qlo,qhi"|none] [--peak-min-prominence 0.2] [--peak-min-snr 5] \
  [--peak-buffer local|global] [--peak-buffer-max N] [--peak-sweep] \
  [--buffer-range "s,e[;s,e]"] [--sample-range s,e] [--baseline none|linear|integral] \
  [--trim-qmin q] [--trim-qmax q] \
  [--guinier-ranges "qlo:qhi,..."] \
  [--model-engine auto|denss|dammif|both|none] [--denss-mode Fast|Slow|Custom] \
  [--n-models 4] [--symmetry P1] [--atsas-dir <ATSAS>/bin]
```

（认峰与逐峰分析的全部参数见下面「多峰 SEC 数据」一节；只跑 `--steps integrate,peaks` 就是
"先看有几个峰"。）

内部只调用 `bioxtasraw.RAWAPI`（= GUI 面板背后的同一套实现）：

| 节点 | RAW 入口 | 产物 |
|---|---|---|
| **认峰（多峰识别）** | `sec_peaks.detect_peaks`（归一化/平滑/找峰沿用 RAW 的 `SASCalc.smooth_data` + `SASCalc.find_peaks`） | `series/sec_peaks.png`（判定图）、`series/sec_peaks_zoom.png`、`tables/sec_peaks.csv`、`series/sec_peaks.json`；**≥2 个峰时**另有 `peaks/peakNN_apexNNNNN/`（每峰一棵完整产物树） |
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
| **结果说明（自动，最后一步）** | 纯读产物（`write-saxs-results-readme` 的 `readme_common.py`，不重算） | **`<产物根>/README.md`** —— 结论速览 + 参数明细表 + 关键结果 + 判读红线 + 本次告警（与管式那条逐条对齐） |

**形状重建两个都要**（`--model-engine auto`）：**电子云（DENSS，`.mrc`）与珠模（DAMMIF，`.pdb`）是同一份 IFT 的两种重建，互不替代**——
DENSS 是 RAW 原生（numba），`Fast` 模式实测几秒到几十秒就出；DAMMIF 是 ATSAS 可执行文件的外壳，**没装 ATSAS 时自动只出电子云**（并在日志里明说），不要因此把 `models/` 留空。

**多区间 Guinier 是"让用户复核拟合过程"的主产物**：不传 `--guinier-ranges` 时，脚本以 `auto_guinier` 的 Rg 为锚
铺一条跨判据边界的阶梯（qRg 0.3–0.6 / 0.4–0.8 / 0.5–1.0 / 0.6–1.3），每个区间**一份 profile 副本**
（`setQrange` 截到该区间）单独存 `.dat` 并进 RAW 报告，同框图用 **RAW 返回的 Rg/I0** 画模型线（不是自写拟合）。
表里给全 `Rg / I0 / rg_err / i0_err / q_min / q_max / qRg_min / qRg_max / r²`，据此判断"哪段才合适"。

### Step 4 — 结果目录里的 `README.md`（交给 `write-saxs-results-readme`）

管线**最后一步自动**调用隔壁技能 `write-saxs-results-readme` 的 `readme_common.write_readme(out, "sec")`
写 `<产物根>/README.md`：只读产物（`run_meta.json` / `tables/*.csv` / `series/ranges.json` / `models/` 清单），
**不重新拟合、不臆造数字**。本 skill **不再自带 README 模板** —— 管式那条流水线用的是同一份实现，
所以两份 README 的章节、表格列、判据与「采用区间/采用那一支 P(r)」口径逐条对齐。
要改格式/加一节，去 `write-saxs-results-readme` 改，别在这里就地改。

README 的 8 节（固定顺序）：结论速览（含「先看这几个文件就够」）→ 最终拟合参数（明细表：环节/参数/值/
程序误差/能不能用/误差怎么读）→ 关键结果 → 这个文件夹里有什么（按建议阅读顺序）→ 每个数字的判据 →
这次没做的 / 不能信的（含各 IFT 引擎全表 + 三条闸门）→ 本次用的参数（复现用）→ 想自己复核（+缩写+目录树）。
开头那两张表必须带**可靠性列**与**误差怎么读**列：程序误差只是该区间内的最小二乘标准误差，
真不确定度看跨区间一致性、换起跑点后参数怎么变。

管线找不到那份实现时**只告警不报错**（README 是交付物，不是计算步骤）。手工补写/验收（旧目录、不动数据）：

```bash
python <write-saxs-results-readme>/scripts/write-readme.py <产物根>       # 补写 README.md
python <write-saxs-results-readme>/scripts/verify-results-folder.py <产物根>   # 验收：结构与关键结果对齐
```

## 多峰 SEC 数据（≥2 个洗脱峰）：逐峰识别 + 逐峰分析

**为什么原来会丢峰**：RAW 自己只认一个峰 —— `SASCalc.findSampleRange()` 里
`max_peak_idx = np.argmax(peak_params['peak_heights'])`（`SASCalc.py:4071`），
`find_buffer_range` 也一样只用最大峰的位置去定搜索窗。所以一条有两个组分的 SEC 系列，
用"RAW 自动"跑出来**只有一条曲线，而且只是最大那个峰的**，第二个峰被静默丢掉
（产物看起来完全正常，只是少了东西）。本 skill 的管线补上这一步：
**先认峰，≥2 个就逐峰建子目录、各用自己的 buffer 与峰窗分别扣减与分析**。

### 多峰 0 — 先看有几个峰（几十秒，不跑下游）

```bash
python run-raw-sec-pipeline.py --series-dir <tif 目录> --cfg <日期>.cfg \
  --out-dir <产物根> --steps integrate,peaks
```

产物：`series/sec_peaks.png`（总览）、`series/sec_peaks_zoom.png`（逐峰放大）、
`tables/sec_peaks.csv`（峰表）、`series/sec_peaks.json`（阈值与全部实测数字）。
**先按"多峰 2"看图**，再决定要不要跑全流程。

### 多峰 1 — 怎么认峰（口径与阈值）

- 色谱图 = **扣减后曲线在 q 窗口 [0.01, 0.05] 1/Å 的积分强度**（`--peak-q-range`，
  `SASM.getIofQRange` 的梯形积分）。**不要用总强度**：实测 BSA 2000 帧的总强度只有 ~10% 起伏、
  且被束位漂移主导，任何阈值都会找出十几个假峰。
- 基线：**滚动中位数**（窗口默认 ≈ n/10、51–401 帧，`--peak-baseline-window`）——
  先把漂移当基线减掉，峰才浮出来。
- 平滑与找峰**沿用 RAW 内部口径**：`SASCalc.smooth_data`（savgol，窗口 min(51, n/2)）+
  `SASCalc.find_peaks`（= `scipy.signal.find_peaks`）。
- 四个下限（都可用 `--peak-min-*` 改）：

  | 参数 | 默认 | 含义 |
  |---|---|---|
  | `--peak-min-prominence` | 0.2 | 峰幅度 / 最高峰（**主判据**）|
  | `--peak-min-snr` | 5 | 峰高 / 基线噪声（1.4826×MAD）|
  | `--peak-min-width` | 5 帧 | 半高宽下限（滤尖刺）|
  | `--peak-min-height` | 0.05 | 峰高 / 最高峰（兜底）|

- **为什么幅度默认卡在 0.2**：扣减后的色谱图本身带**系统性纹波**（单个全局 buffer 平均 ×
  束位漂移 → 剩余起伏可达主峰的 10–20%）。实测 4EH2 那条数据里，纹波在三个 q 窗口都能复现、
  峰值 7–15σ、半高宽 10–18 帧 —— "相对噪声"与"多 q 窗口一致性"这类判据**都杀不掉它们**，
  只有"相对主峰的幅度"能。要找 <20% 主峰的小组分：调低阈值，然后**必须**逐个人工看图。
- 无束流帧（总强度 < 0.5×全系列中位）在峰检测与 buffer 挑选里都被剔除：实测 4DH1-KDPV-ZN
  有 16 帧这种帧，扣减后是一根巨大的负尖刺，会被 find_peaks 当成峰。

### 多峰 2 — 视觉复核是流程的一步，不是可选项

管线自己算一份"必须看图"的条件清单（`sec_peaks.needs_visual_check`），命中任一就写进日志、
`README.md` 与 `sec_peaks.json` 的 `visual_check_reasons`，并在 README 里明写
**"看图之前别引用数字"**：

- 峰刚过阈值（幅度 < 3×阈值 或 SNR < 15σ）→ 可能是纹波/肩峰；
- 相邻两峰：谷底/峰高 ≥ 0.2 **或** 经典分辨率 R = 2Δt/(w1+w2) < 1.5 → 两峰可能没真分开；
- 某个峰的 buffer 段落在**两峰之间**（谷底）→ 会带进另一个组分；
- 某个峰的 sample 窗 < 10 帧 → Guinier/IFT 不可靠；
- 两侧 buffer 段**水平差 ≥ 0.15×主峰** → 峰前/峰后基线不在同一水平；
- 有疑似无束流帧，或峰顶紧邻这种帧。

**看图判四件事**：
1. 峰的个数/位置与**裁剪区视频**或 UV 痕对得上吗？
2. 每个峰的 sample 窗（彩带）框住峰中心了吗？两侧 buffer 段（斜纹）落在**真基线**上，
   而不是肩部或谷底？
3. 逐峰放大图里峰是否回到基线？两峰之间的谷底是**降到基线**（真分开）还是**只降一点点**（未解析重叠）？
4. 总览图下面三栏里，每个峰在**自己的 buffer** 下有没有一段 Rg/MW **平台**？

### 多峰 3 — 改参数只看峰：`find-sec-peaks.py`（秒级，不重跑积分）

```bash
python find-sec-peaks.py --products <产物目录>                    # 看自动结果
python find-sec-peaks.py --products <目录> --sweep                # 出一张阈值扫描图，人眼挑一组
python find-sec-peaks.py --products <目录> --peak-min-prominence 0.1 --peak-q-range 0.02,0.06
python find-sec-peaks.py --products <目录> --peak-ranges "690,720;800,840"    # 人工给区间
python find-sec-peaks.py --products <目录> --buffer-range "560,620;800,880"   # 换 buffer 重算扣减（走 RAW）
```

- 它直接读 `profiles/*.dat` 做梯形积分（与 RAW 的 `getTotalI` / `getIofQRange` 逐点一致），
  **不重跑图像积分**；
- `--sweep` 把 6 组阈值画在同一张图上（每组标注峰数），人眼挑最符合峰形的一组；
- 每次都会打印**可直接搬回主管线的参数**（`--peak-ranges` / `--peak-buffers` / 阈值）。

### 多峰 4 — 把视觉结论写回主管线

```bash
# 自动结果可用 → 什么都不用给（默认 --multi-peak auto）
# 人工确认了峰窗与 buffer → 原样搬回去
python run-raw-sec-pipeline.py ... \
  --peak-ranges "801,832;943,977" \
  --peak-buffers "0,800;833,942|833,942;978,1499"
```

- `--multi-peak auto`（默认）：≥2 个峰 → 逐峰分析；只有 1 个峰 → 走原来的单峰流程
  （RAW 自动 buffer + 自动样品区，行为与从前一致）。
- `--multi-peak always`：即使只有 1 个峰，也用它**自己的**前后两段 buffer + 峰窗。
- `--multi-peak off`：不分峰，只出检测图/表（"其实只有/有两个峰"这个信息本身就有用）。

### 多峰 5 — 产物布局

```
<产物根>/
  README.md                  ← 索引：峰表 + 每个峰的入口 + "看图之前别引用数字"清单
  tables/sec_peaks.csv       ← 峰表（apex/窗/幅度/SNR/半宽/buffer/谷底比/R/告警）
  series/sec_peaks.png       ← 判定图（峰窗 + buffer + 每峰自己的 Rg/I(0)/MW）
  series/sec_peaks_zoom.png  ← 每峰一张放大图
  series/sec_peaks.json      ← 全部阈值与实测数字（可复现）
  run_meta.json              ← n_peaks / peaks[] / visual_check_reasons
  profiles/01_integrated/    ← **共享**（图像只积分一次）
  peaks/peak01_apex00812/    ← 每个峰一个完整产物目录
  peaks/peak02_apex00957/    ←   README.md / run_meta.json / profiles/ / tables/ / ifts/ / models/ / reports/ / series/
```

每个峰的子目录 = 那一峰单独跑一遍的结果，只是 `profiles/01_integrated` 是**相对符号链接**
（不重复占盘）。**每个峰的扣减是用本峰的 buffer 段重算的**，所以子目录里的 Rg/I(0)/MW
是该峰自己的口径；顶层 `sec_peaks.png` 把各峰曲线拼在一起，方便比"哪个峰的 Rg 平台平"。
峰子目录里的 `README.md` 由共享实现 `write-saxs-results-readme` 写（与单峰跑一遍同款 8 节模板）；
顶层 `README.md` 是本 skill 写的**索引**（峰表 + 每峰入口）。

### 多峰 6 — buffer 与边界的判断规则

- 每个峰的 buffer = **紧邻它的前一段 + 后一段非峰区**（`--peak-buffer local`，默认）——
  就是"峰前+峰后两段"，因为 SEC 峰后的基线常不等于峰前。`global` = 所有峰共用全部非峰区。
- **两峰之间的谷底不是纯 buffer**（是两组分的混合）：管线把它标成 `between-peaks` 并列入
  看图理由。若谷底明显没回到基线，正确做法不是"再抠一段 buffer"，而是转
  `deconvolve-overlapping-elution-peaks`（SVD/EFA/REGALS，判据是组分数与浓度曲线，不是平台）。
- 峰窗取法：`--peak-window half-height`（默认，半高宽，保守）或 `valley`（到相邻峰谷底，整个峰）。
- `--peak-ranges` 与所有区间都是 **0 基闭区间**。
- 多峰模式下 `--baseline linear` **自动用该峰自己的 buffer 段当初末锚点**
  （不再用写死的 `0,20;1800,1990`）。

### 多峰 7 — 实测锚点（本机 5 条数据，2026-10-01）

| 样品 | 帧数 | 检测结果 | 结论 |
|---|---|---|---|
| 4LI2-676 | 1800 | **1 峰**（apex 904，窗 888–933，46 帧） | 单峰，走原流程，无看图理由 |
| 4EH2-KDPV-ZN | 1500 | **2 峰**：812（1.00×）+ 957（0.71×），R=4.4 | 峰间谷底回到基线 → 真分开；主峰外确有第二个组分 |
| 4DH2-676-apo-3 | 1500 | **2 峰**：698（0.46×）+ 810（1.00×） | 前峰 0.46× 且谷底平 → **必须看图**确认是独立组分还是肩 |
| bsa | 2000 | 2 峰：616（0.23×）+ 689（1.00×），但**两侧 buffer 水平差 ≈0.6×主峰** | 漂移主导 → 先 `--baseline` 或换 buffer，再谈分不分峰 |
| 4DH1-KDPV-ZN | 1200 | 16 帧无束流；遮掉后仍有帧 28 处尖刺（7σ、0.55×） | 野值帧；看图确认后再决定要不要分析这个"峰" |

（"×"= 该峰 prominence / 最高峰 prominence。）

### 多峰 8 — 什么时候**不要**用多峰这一套

- 谷底只降到很小（谷底/峰高 ≥ 0.5）→ 那是**未解析重叠**，不是"两个峰"：转
  `deconvolve-overlapping-elution-peaks`；硬切区间只会把两条曲线一起弄脏。
- 数据本身不干净（无束流帧、扣减过头、低 q 上翘）→ 先 `assess-saxs-raw-data-quality`；
  多峰识别的全部前提是"扣减后的曲线可信"。
- 洗脱事件只抬高几个百分点、噪声又大的数据 → 阈值怎么调都分不出峰，别硬分；
  按判据写一份"数据不可用/不可分辨"的说明（见下面的"结果不可用-README"那条）。

## 峰内梯度：把**一条峰的不同位置**当稀释序列用

SEC 洗脱峰的**浓度沿帧变化**（峰顶最浓、两翼渐稀），所以同一条峰天然是一条 **c 递减序列**，
而且同一 buffer、同一批样品，**不需要另配稀释样品**。流水线把它做成一个节点（`--gradient`，默认开），
每峰产出 `gradient/gradient_slices.csv`、`gradient.json`、`gradient.png` 与 README 第 8 节。

**能回答三个原流程答不了的问题**：

1. **Rg / MW 的浓度依赖**：理想单分散体系 Rg 不随 c 变；Rg 随 c **下降** → 排斥型相互作用
   （或浓度依赖解离）；Rg 随 c **上升** → 吸引/自缔合（低聚态，再看 MW 是否同向）。
2. **c→0 外推**：把 Rg 外推到无限稀释（去掉相互作用）→ `Rg0 ± err`（c→0 截距）与直接中位值对比。
3. **两翼互检**：同一个 c 在峰的上升翼与下降翼各出现一次，两者的 Rg/MW 必须在误差内一致——
   比任何单点判据都更容易暴露基线漂移、组分变化、扣减残差。

**三条铁律（不遵守这个节点就没意义）**：

- **所有切片共用同一条 q 区间**（默认沿用主分析的 `auto` 区间，`--gradient-q` 可显式给）。
  换区间 → Rg 差异里混进"取点差异"，就不能归因于浓度。
- **浓度是相对的**（扣减后低 q 窗口积分 / 峰顶值）。SEC 的绝对浓度未知——所以结论只能写
  "随（相对）浓度如何变化"、以及 c→0 外推值，**不能**写"在 x mg/mL 时"。
- **切片要连续**（`--gradient-mode frames`，默认）：一个切片内的帧是**连续**帧块，块内组分几乎不变。
  `height` 模式按相对峰高分带（帧可能不连续），适合峰形不规则时。

```bash
python run-raw-sec-pipeline.py --series-dir <tif 目录> --out-dir <产物根> --cfg <日期>.cfg \
  --gradient auto|off [--gradient-slices 5] [--gradient-mode frames|height] \
  [--gradient-min-frames 3] [--gradient-min-conc 0.08] [--gradient-q "qlo,qhi"|auto]
```

**产物与怎么读**：`gradient_slices.csv` 每行一个切片（`side` 分 asc/desc/band/all、`conc_rel` 相对浓度、
`window` 帧区间、Rg±err、I(0)、qRg 边界、r²、MW）；`gradient.png` 三格 =
峰上切片位置（蓝=上升翼、橙=下降翼、黑=整窗）→ **Rg vs c**（带 c→0 外推线）→ MW 与 I(0) vs c。
判词按"斜率是否显著（|t| = |slope|/SE ≥ 2 且变化 > 2×中位 Rg 误差）"三档给：
🟢 无浓度依赖 / 🟡 随浓度升或降（并写清方向与可能机制）/ 🔴 不可判（可用切片不足）。

**自检三条**：

- `I(0)` 对 c 的**双对数斜率应 ≈1**（理想稀释线性）。偏离时先怀疑**浓度代理本身**（色谱窗口、滚动中位
  基线对宽峰峰顶的压缩）与扣减线性——它只影响 c 轴的标度（也就是 c→0 外推的那个截距），
  **不影响**"Rg 随 c 升/降"这个方向判断。本机实测两峰都是 0.43–0.45，方向结论不受影响，
  但外推值当"量级参考"而不是精密值。
- Rg 误差**必然随 c 降低变大**（越稀越差）——所以看趋势要带误差棒，不要拿最稀的点下结论。
- 峰顶附近的切片不一定最可信（上样/柱头效应、聚集峰常在最前沿）。

**实测锚点（本机 4EH2-KDPV-ZN，两峰 1500 帧，同一 q 区间）**：主峰 P1（apex 812，窗 801-832）
5 个切片 + 整窗，相对浓度 0.28–0.88；判词与数字见该结果目录的 `README.md` 第 8 节与
`gradient/gradient_slices.csv`。

**什么时候别用它**：峰窗只有几帧（切片凑不出 ≥2 个可用）；洗脱事件只抬高几个百分点（浓度梯度太小，
c 范围不够）；峰内有明显肩/双组分（那先分峰或去卷积，别把混合物的"浓度依赖"当相互作用）。

## 复核点（跑完先看这几处，再看数字）

0. **`series/sec_peaks.png` + `README.md` 的"看图之前别引用数字"清单（多峰时先看这个）**：
   几个峰、每个峰的窗与 buffer 落在哪、`visual_check_reasons` 里有没有条目。
   有 → 按「多峰 2」的眼看四问复核，或 `--peak-ranges/--peak-buffers` 覆盖后重跑。
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
  ③ 两样都没有（本机 **`4EH2-KDPV-ZN` 的 445 帧旧批次**：无 `.Iochamber`/`.log`/`.txt`）→ 只能
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
- **帧区间是 0 基闭区间，`end` 给到帧数就会 `IndexError`**：445 帧的系列（4EH2 旧批次）写 `--buffer-range "60,220;400,445"`
  → `SECM.averageFrames` 里 `list index out of range`（本机实测）。脚本现在有 `clip_ranges` 自动收到 `n−1` 并告警，
  但**区间里的数字仍应自己核对**（峰位置看副产物色谱图最保险）。
- **判"这条系列能不能用"的硬指标**（任一命中就是数据问题，不是参数问题）：Guinier `Rg < 5 Å` 或 `r² < 0`；
  `Vp` 为负；逐帧 `Rg` 在几十到几埃之间乱跳；DENSS 建不起样条。本机 `4EH2-KDPV-ZN` 的 **445 帧旧批次**（无监视器、
  洗脱事件只抬高 ~3.7%、两种区间都给出 Rg 1.27–1.29 Å）就是这种：**在产物目录写一份
  `结果不可用-README.md` 说明判定依据**（数据事实 + 试过的区间对照 + 三重证据 + 建议），别让表格里的数字被当结果引用。
- **脚本被并发编辑时，先 `git status` 再跑，必要时用快照跑**：本机同一仓库有别的会话在改同一个脚本，
  撞到过两次半成品（`step_ift` 改返回 5 值而 `main` 还解 4；`step_shape` 里 `itf`/`ift` 拼写），
  两次都是"跑到一半 ValueError/NameError"。稳妥做法：`git log -1` 记下 revision → `cp` 到 scratch 当快照 →
  跑快照 → 报告里写清用的哪个 sha。
- **"结果看不懂"是缺交付物，不是缺解释**：产物目录必须有一份 `<产物根>/README.md`（管线最后一步由
  `write-saxs-results-readme` 自动写，旧目录用它的 `write-readme.py <目录>` 补）。它只读产物、不重算，
  **必须最后写**（README 比 `tables/*.csv` 旧 = 读者拿到旧数字；`verify-results-folder.py` 会把这条标红）。
  数据本身不可用的系列，另写 `结果不可用-README.md`（见下条）。
- **表里 `rg = -1` 不是负数，是 RAW 的失败哨兵值**（该区间的 Guinier 拟合没收敛/点数不够）——按"此区间不可用"读，
  不要当成数值；同理 `r² < 0` 表示拟合比取平均还差。
- **裁剪坐标：两个轴都从"大的那头"数** → `row = H−1−y`、**`col = W−1−x`**（即 180° 旋转；本机 BL19U2
  实测标定：`col=x` 时与旧视频帧相关 0.11，`col=W−1−x` 时 0.66，行向扫描峰值恰为 `H−1−y`）。脚本
  `crop-video-normalized.py` 默认 `--x-origin right`。裁出来一片均匀背景/位置偏移说明约定错了，别微调数字，
  回 `frame-sequence-to-video` 的"用已接受的旧帧算相关系数"重标法。
- **ffmpeg 的 libx264 + yuv420p 要求偶数边长**：61×71 的裁块靠整数放大（8×）顺带解决。
- 视频的 `-pix_fmt rgb24` 必须与写进管道的字节一致（写 RGB 就声明 rgb24），否则帧数会变 3 倍。

### 多峰相关的坑（都是实测撞出来的）

- **逐峰 `tables/frame_params.csv` 里 rg/i0/vc/vp 可能**整列都是 `-1`**，而且**不报错**：
  RAW 只对"可用帧"算这些参数（`SECM.subtractAllSASMs`：该帧强度 / buffer 平均强度 > `calc_thresh`(1.02)
  → 标记可用；`SASCalc.run_secm_calcs` 只对**连续 `window_size`(默认 5) 帧都被标记**的窗做计算）。
  默认按**总强度**判定，而 SEC 数据的总强度被束位/通量漂移主导——本机实测 `4EH2-KDPV-ZN`：1500 帧里
  只有 **69 帧**被标记、且全是散落的噪声帧 → 一个 5 帧窗都凑不齐 → `run_secm_calcs` 把 rg/i0/vc/vp
  **全部写成 -1**。表里看着像"算了，只是值是负"，其实是**根本没算**（`SASCalc.py:2766` 的兜底），
  而 `-1` 与"拟合失败"的哨兵值**长得一模一样**。处置：用 `--frame-flag-q "0.01,0.05"`（默认已开）
  把判定换成**低 q 窗口积分强度**（与认峰同一窗口）→ 同一数据 P1 从 0 帧变 **112 帧有值**、
  其中 34 帧落在 P1 自己的峰窗内，Rg 中位 24.0 Å（与它的 Guinier 24.3 Å 对得上）。
  要回到 RAW 默认口径传 `--frame-flag-q none`。**看到逐帧参数整列 -1 先查这个，别去调 Guinier 区间。**

- **RAW 自己只认最大的那个峰**：`SASCalc.findSampleRange()` 里 `max_peak_idx =
  np.argmax(peak_params['peak_heights'])`，`find_buffer_range` 也只用最大峰定搜索窗。
  → 一条有两个组分的系列，用"RAW 自动"只会给你**最大峰那一条**曲线，第二个峰**静默消失**。
  本管线默认 `--multi-peak auto` 处理；即使 `--multi-peak off`，`tables/sec_peaks.csv` 也会
  告诉你实际认到几个峰。
- **认峰别用总强度**（`series.getIntI`）：实测 BSA 2000 帧总强度只有 ~10% 起伏、被束位漂移主导，
  照 RAW 内部口径（归一化后 height=0.4）能找出 22 个"峰"。用**低 q 窗口 [0.01,0.05] 的积分强度**
  （本 skill 默认，`--peak-q-range`）——低 q 对溶质敏感、对噪声不敏感。
- **扣减后的色谱图带系统性纹波**：单个全局 buffer 平均 × 束位漂移 → 剩余起伏可达主峰的 10–20%。
  实测 4EH2 那条数据里，纹波在三个 q 窗口都复现、峰值 7–15σ、半高宽 10–18 帧，
  **"相对噪声"与"多 q 窗口一致性"这两类判据都杀不掉它**。这就是 `--peak-min-prominence`
  默认卡 0.2×主峰的原因；要挖更小的组分，调低阈值后**必须**逐个人工看图。
- **谷底判据有两个必踩的坑**（写代码时实测撞到）：① 谷底要在**平滑过**的曲线上取，
  否则噪声让谷底随便跌到 0 以下，任何两峰都会被判成"已分开"；② 峰高是"相对最高峰"的
  归一化值、谷底是绝对强度，**得同除一个 top 才能比**。另外给一个经典色谱分辨率
  `R = 2Δt/(w1+w2)`（≥1.5 基线分离、1.0–1.5 部分重叠、<1.0 未分开），**两个判据一起看**。
- **无束流/断束帧**（整帧总强度 < 0.5×全系列中位）：实测 4DH1-KDPV-ZN 有 16 帧，扣减后是一根
  巨大的负尖刺，会被 find_peaks 当峰。现在检测与 buffer 挑选都剔除它们，但**峰顶紧邻这种帧**
  的峰仍会被标出来要人工看图（不要拿它当结果）。
- **多峰模式下每个峰都要重跑一次 `set_buffer_range`**（因为要用本峰自己的 buffer 段），
  所以逐帧 Rg/I(0)/MW 是**逐峰口径**；顶层 `sec_peaks.png` 把各峰曲线拼在一起才看得出
  "哪个峰的 Rg 平台平"。
- **`--steps integrate,peaks` 是最省的第一步**：本机 1500 帧 ≈ **85 s**（含积分）出峰表与判定图，
  先看几个峰再决定要不要跑下游，别一上来就跑全套。
- **多峰时顶层 `README.md` 是索引**（由 `sec_peaks.write_peaks_readme` 写），每个峰自己的
  README 由共享实现 `write-saxs-results-readme` 写在峰子目录里。**别**用
  `write-saxs-results-readme/scripts/write-readme.py <产物根>` 覆盖顶层（它写的是单样品版，
  会把索引换掉；要补写就补到峰子目录上）。
- **峰子目录里的 `profiles/01_integrated` 是相对符号链接**（图像只积分一次，不重复占盘）：
  拷贝/打包产物时要跟随链接（`cp -rL`、`tar -h`），否则那 1500 份 .dat 不在包里。
- **`find_buffer_range` 会在整条系列上失败**（返回 `success=False`）：实测 4EH2-KDPV-ZN
  1500 帧就是这种。认峰因此退回**未扣减**曲线（不依赖 buffer），并照常出峰表——
  不要因为"RAW 自动 buffer 失败"就以为这条系列没法处理。

## 相关 skills

- **process-sec-saxs-series** — 本 skill 的"判据版"：区间怎么选、峰上 Rg 平台怎么看、SEC 为什么不能用绝对刻度。
- **deconvolve-overlapping-elution-peaks** — `contrasts-with`：两峰之间的**谷底没回到基线**时，
  正确做法不是切区间，而是做 SVD/EFA/REGALS 分解（判据是组分数与浓度曲线）——本 skill 的
  多峰流程只处理"峰之间有基线"的情况。
- **assess-saxs-raw-data-quality** — 多峰识别的全部前提是"扣减后曲线可信"；无束流帧、扣减过头、
  低 q 上翘这些先判数据，再谈分几个峰。
- **script-raw-with-the-python-api** — RAWAPI 的骨架与返回元组顺序（本 skill 的脚本就是它的落点）。
- **correct-sec-saxs-baseline** — 扣减后仍漂移时选 Linear/Integral，以及本 skill 用的监视器归一套路。
- **configure-bioxtas-raw-for-a-dataset** — 换实验日/仪器时先核对 `.cfg`。
- **write-saxs-results-readme** — 结果目录的 README 与验收（本 skill 第 4 步调它；管式那条也调它，
  所以两份 README 逐节对齐）。
- **evaluate-a-shape-reconstruction** / **fit-a-high-resolution-model-to-data** — 重建与高分辨模型对照。
