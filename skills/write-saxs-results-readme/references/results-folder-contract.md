# 结果目录契约（结构 + 章节 + 字段别名）

本文件是 `scripts/readme_common.py` 里 `STRUCTURE` / `SECTION_ORDER` / `FILE_DOC` / `RULES` / `ALIAS`
的文字版。**改代码必须同步改这里**（反之亦然），两份不一致 = 验收与实际行为不一致。

## 1. 结构契约（`STRUCTURE`）

级别：**must** 缺 = 交付不完整（验收 🔴）；**should** 缺 = 第 5 节必须写明原因（🟡）；
**info** 有就说（缺了不扣分）。目录存在但 0 文件，非 must 记 🟡。

### SEC-SAXS（标记文件 `run_meta.json`）

| 级别 | 路径 | 说明 |
|---|---|---|
| must | `run_meta.json` | 本次运行的参数与关键结果（机器可读） |
| must | `profiles/01_integrated` | 逐帧积分曲线（含逐帧归一化） |
| must | `profiles/03_subtracted` | 逐帧扣减曲线（系列图的数据源） |
| must | `profiles/04_sample` | 样品区平均曲线（下游分析用的就是它） |
| must | `tables/guinier_multi_range.csv` | 多区间 Guinier 全表 |
| must | `tables/ift_summary.csv` | 各 IFT 引擎的 Dmax/Rg/χ² |
| must | `tables/mw.csv` | 分子量（Vc / Vp / 贝叶斯 / DATCLASS） |
| should | `series` | 序列对象、洗脱曲线与区间记录（`ranges.json`） |
| should | `profiles/02_buffer` | buffer 平均曲线 |
| should | `profiles/06_guinier` | 每个 q 区间的 Guinier 拟合副本 |
| should | `reports` | RAW 自己出的 PDF 报告 |
| should | `ifts` | IFT 解（P(r) 与 RAW/ATSAS 原文件） |
| info | `tables/frame_params.csv` | 逐帧 Rg/I0/Vc-MW |
| info | `tables/frames_integrated.csv` | 逐帧透射与低 q 强度 |
| info | `models` | 3D 重建（电子云 `.mrc` / 珠模 `.cif` / DAMAVER） |
| info | `video` | 归一化裁剪视频与束斑质心轨迹 |
| info | `norm` | 逐帧归一化因子表 |

### 管式 / 静态（标记文件 `summary.json`）

| 级别 | 路径 | 说明 |
|---|---|---|
| must | `summary.json` | 参数、关键结果与逐节点状态（机器可读） |
| must | `tables/guinier_multi_range.csv` | 多区间 Guinier 全表 + 每条闸门 |
| must | `tables/ift_summary.csv` | 各引擎/起跑点的 Dmax/Rg/χ² + 三条闸门 |
| must | `tables/mw.csv` | 分子量（Vp / Vc / 贝叶斯 / DATCLASS） |
| must | `profiles/01_control` | 对照（buffer）平均曲线 |
| must | `profiles/02_sample` | 样品平均曲线 |
| must | `profiles/03_subtracted` | 扣减曲线（最终交付的那条） |
| must | `profiles/04_guinier` | 每个候选区间一份截断后的曲线 |
| must | `ifts` | IFT 解与 P(r) |
| should | `tables/shape_results.json` | 形状重建全表（电子云 / 每个珠模 / DAMAVER） |
| should | `qc.png` | 四联诊断图（扣减曲线 / Kratky / 多区间 Rg 与 χ² / P(r)） |
| should | `reports` | RAW 自己出的 PDF 报告 |
| should | `models` | 3D 重建（电子云 `.mrc` / 珠模 `.cif` / DAMAVER） |
| info | `norm/frame_qc.csv` | 逐帧质检（低 q 电平 / 透射 / 离群） |
| info | `frames` | 逐帧 1D 曲线（仅在 `--save-frames` 时） |

## 2. 章节规范（`SECTION_ORDER` + 各节内容）

**两种模式必须有这 8 节，且顺序一致**（验收会逐节比对标题）：

```
## 0. 结论速览
## 1. 最终拟合参数（明细表）
## 2. 关键结果
## 3. 这个文件夹里有什么（按建议阅读顺序）
## 4. 每个数字的判据
## 5. 这次没做的 / 不能信的
## 6. 本次用的参数（复现用）
## 7. 想自己复核          ← 内含 ### 这些缩写、### 目录树（两层）
```

### §1 的固定行（`环节 / 参数 / 值 / 程序误差 / 能不能用 / 误差·不确定度怎么读`）

| 环节 | 参数 | 出现条件 |
|---|---|---|
| Guinier 拟合 | Rg | 有采用区间 |
| Guinier 拟合 | I(0) | 同上 |
| Guinier 拟合 | R² / χ²_red / qRg | 同上（SEC 的表没有 `chi2_red` 一列时显示 —） |
| Guinier 拟合 | Rg（auto，仅参考） | 没有采用区间但有 auto |
| IFT / P(r) | Dmax | 有采用的那一支 |
| IFT / P(r) | Rg（实空间） | 该支给了 Rg |
| IFT / P(r) | χ²（拟合优度） | 该支给了 chisq |
| GNOM | Total estimate / α | 该支是 GNOM 且带这两个字段（管式的 `summary.json.ift.gnom`） |
| 分子量 | Vc（浓度无关） | `mw.csv` 有 Vc |
| 分子量 | Vp（Porod 体积） | 有 Vp（体积取**修正后**的 `pvol_cor`） |
| 分子量 | DATMW（贝叶斯） | 有（CI 与区间外概率一起给） |
| 分子量 | DATCLASS | 有（形状分类 + 第二返回值） |
| 珠模 DAMMIF（N 个） | χ² | 有珠模 |
| 珠模 DAMMIF | Rg / Dmax / MW | 有珠模 |
| 珠模一致性 DAMAVER | 平均 NSD | 有 DAMAVER |
| 电子云 DENSS | χ² / 模型 Rg | 有 DENSS |

「能不能用」列的颜色规则：珠模 χ² ≤3 绿 / 3–5 黄 / >5 红；NSD ≤2 且单 cluster 绿、分簇黄；
DENSS 在 IFT 不可信时一律红（无意义）；Dmax/Rg ≤4.5 才算过闸门。

### §2 的固定行（`项目 / 数值 / 怎么看`）

对照与缩放（SEC 换成 buffer/样品帧区间）→ 归一化 → 低 q 对比度 → 分析窗 →
auto-Guinier → 采用的 Guinier 区间 → 跨区间一致性 → P(r)/IFT 采用那一支 →
分子量 Vc / Vp / DATMW / DATCLASS → 电子云 DENSS → 珠模 DAMMIF → DAMAVER 一致性。

### §3 的条目字典（`FILE_DOC`，只列真实存在的）

条目 = `(适用模式, 通配, 是什么, 什么时候看)`，顺序就是建议阅读顺序：README → `qc.png`（管式）→
`reports/*` → 扣减曲线 → 逐帧曲线 → 平均曲线 → `tables/guinier_multi_range.csv` → 多区间副本 →
`tables/ift_summary.csv` → `ifts/pr*.dat` / `*.ift` / `*.out` / `ift_fit*.dat` → `tables/mw.csv` →
`models/*`（电子云 → 珠模 → DAMAVER → distances）→ `tables/shape_results.json` → `norm/*` →
`*.hdf5` → `series/*`（SEC）→ `video/*`（SEC）→ `tables/frame_params.csv`（SEC）→ `summary.json` /
`run_meta.json` → `frames/*`（管式）。清单里没提到的文件在末尾汇总成一行「还有 N 个文件」。

### §4 的判据表（`RULES`，两种模式完全相同）

高 q 缩放因子 / 低 q 对比度 / Guinier R² / 多区间 Rg 一致性 / IFT 三条闸门 / Dmax÷Rg / Vc 分子量 /
DENSS χ² / DENSS support 体积 / 珠模 χ² / 珠模 Rg·Dmax / DAMAVER 平均 NSD —— 每行三档（好 / 勉强 / 差）。

### §5 的结构

跳过的节点（如「没有可信的 P(r)」）+ 失败原因（DENSS / 珠模 / 归一化缺失）+ 告警 +
**IFT 各起跑点/引擎全表**：`起跑点 / 引擎 / q 范围 / Dmax / Rg_real / chisq / Rg 闸门 / Dmax÷Rg 闸门 /
网格闸门 / 可信`；SEC 没有网格信息时该列显示 —。

## 3. 字段别名表（`ALIAS`：canonical → 产物里的列名）

| canonical | SEC | 管式 | 备注 |
|---|---|---|---|
| `label` | `range_label` | `tag` / `idx_min-idx_max` | 管式没有区间名时用 `idx a-b` |
| `rg` / `rg_err` | `rg` / `rg_err` | `Rg` / `Rg_err` | Guinier 表 |
| `i0` / `i0_err` | `i0` / `i0_err` | `I0` / （无） | |
| `qmin` / `qmax` | `q_min` / `q_max` | `qmin` / `qmax` | |
| `qrg_min` / `qrg_max` | `qRg_min` / `qRg_max` | 同左 | |
| `r2` | `r_sqr` | `R2` | |
| `chi2_red` | （无） | `chi2_red` | SEC 的表没这一列 |
| `dmax` / `dmax_err` | `dmax` / `dmax_err` | `Dmax` / `Dmax_err` | |
| `rg_real` | `rg`（IFT 表里的 rg 就是实空间） | `Rg_realspace` | **注意 guinier 表的 `rg` 是另一回事** |
| `chisq` | `chi_sq` | `chisq` | |
| `method` | `method` | `engine` | |
| `mw` | `mw` | `MW_kDa` | |
| 闸门 | （无，现算） | `pass_rg` / `pass_dmax_over_rg` / `pass_dmax_within_grid` / `trusted` | 管式产物自带的优先 |

### 分子量四法的返回值顺序（`mw.csv` 的 `detail1..4` = 管式的 `aux[0..3]`）

| RAW 入口 | 返回值 | 本 README 怎么用 |
|---|---|---|
| `mw_vc` | `(mw, vcor, mw_err, qmax)` | MW + ± 误差（`vcor` 是体积、`qmax` 是上限） |
| `mw_vp` | `(mw, pvol_cor, pvol, qmax)` | MW；体积报**修正后**的 `pvol_cor` |
| `mw_bayes` | `(mw, mw_prob, ci_lower, ci_upper, ci_prob)` | MW + CI + 区间外概率 |
| `mw_datclass` | `(mw, 形状分类, …)` | MW + 形状分类（第二返回值，`-1` = 没给） |

`rg = -1`、`-1.0 kDa`、`r² < 0`、`inf`、`nan` 是**哨兵/失败值**，显示成「未收敛（哨兵 -1）」
「inf（发散）」「—」，不当数字用。

## 4. 统一的选值口径（代码位置）

| 口径 | 函数 | 说明 |
|---|---|---|
| 采用区间 | `pick_adopted_guinier()` + `declared_verdict()` | 产物闸门全过才采信；标了 NOT recommended 就退回通用规则 |
| 跨区间一致性 | `guinier_spread()` | 只统计收敛且 R² ≥ 0.9；被排除的列出来 |
| 采用那一支 P(r) | `pick_adopted_ift()` | chosen / trusted → 闸门 → χ² |
| 三档可用性 | `_availability()` | 见下表 |
| 缺什么 / 不能信什么 | `structure_check()` / `missing_items()` / `flags()` | 结构契约 / 第 5 节 / §0 警示块 |

| 数据可用性 | 条件 |
|---|---|
| 🟢 可用 | 采用区间 R² ≥ 0.99 且区间间 Rg 漂 ≤ 5% 且 IFT 有可信的一支（没有 IFT 节点时不算这条） |
| 🟡 有限可用 | 有可引用的 Rg，但：跨区间漂 5–15% / R² 0.95–0.99 / 产物标了 NOT recommended |
| 🔴 不可用/存疑 | 没有可引用的 Rg（全部未收敛、没有区间过闸门），或采用区间 R² < 0.95 |
