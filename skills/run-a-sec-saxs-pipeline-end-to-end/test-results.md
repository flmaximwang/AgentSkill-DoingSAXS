# test-results.md — `run-a-sec-saxs-pipeline-end-to-end`

本轮（2026-10-02）为**多峰 SEC-SAXS 缺陷修复**，version 1.1.0 → 1.2.0。
记录两类凭据：**① 路由盲测（三次，含一次回退判断）② 真实数据实跑证据（4 次运行，全部 exit 0）**。

---

## 0. 这次修的是什么

**缺陷**：一条 SEC-SAXS 系列如果有 ≥2 个洗脱峰，原流程只出**最大那个峰**的一条曲线，其余峰**静默消失**
（产物看起来完全正常，只是少了东西）。用户报告现象，根因在 RAW 自己身上：

```
/Applications/BioXTASRAW/lib/python3.12/site-packages/bioxtasraw/SASCalc.py:4071
    max_peak_idx = np.argmax(peak_params['peak_heights'])      # findSampleRange() 只认最高峰
    ...  find_buffer_range() 的搜索窗也是绕着同一个最大峰定的
```

**修法**：认峰 → ≥2 个峰就**逐峰建子目录**（每峰用"紧邻本峰的前后两段 buffer + 本峰峰窗"各跑一遍
扣减与下游）→ 顶层 README 作索引。**代码判不干净的部分做成视觉复核的一步**（用户明确要求：
"不一定完全用代码解决，可以在代码难以判断时，fallback 到使用视觉检验"）。

---

## 1. 路由盲测：方法

| 项 | 值 |
|---|---|
| 候选 skill | **21 个**（default profile `saxs` 类目全部） |
| 可见信息 | 每个 candidate 只有 **description 前 57 字符**（`agent/skill_utils.py:761` + `...`），**看不到 SKILL.md 正文** |
| 题目 | 本 skill `test-prompts.json` 的 **9 条 = 5 正面 + 4 诱饵**，交错排列，评测者不知正负 |
| 评测者 | **2 个独立子代理**（互不可见，各自独立判 9 题） |
| 判据 | 正面应判本 skill；诱饵应判**指定的那条**兄弟 skill（不是"只要不是本 skill 就算对"） |

题目与金标（编号 = 盲测时的交错顺序）：

| # | 题面（摘要） | 金标 |
|---|---|---|
| 1 | 把一条 SEC-SAXS 系列端到端跑一遍（图像→归一化→扣减→Rg/MW→IFT→珠模→报告） | 本 skill |
| 2 | 扣减完的曲线 Rg 取哪一段 / q_max 到多少 | `assess-guinier-fit-quality` |
| 3 | 只有 `.Iochamber` + 采集日志、没有每帧 txt，怎么在 RAW 里做逐帧归一化 | 本 skill |
| 4 | 扣减完基线还在抬，Linear 还是 Integral | `correct-sec-saxs-baseline` |
| 5 | 想自己判断拟合过程合不合适，能不能把不同数据范围都拟合一遍 | 本 skill |
| 6 | **色谱图上有两个洗脱峰，第二个峰也要一套完整结果**（新增） | 本 skill |
| 7 | 两峰太近、谷底没回基线，怎么分成两条曲线（新增诱饵） | `deconvolve-overlapping-elution-peaks` |
| 8 | 看束斑/样品有没有漂做个小视频，别写归一化 tif | 本 skill |
| 9 | 珠模一堆模型被剃掉、平均 NSD 0.8 算好吗 | `evaluate-a-shape-reconstruction` |

---

## 2. 三轮结果

| 轮次 | 本 skill 的 57 字符头 | A | B | 合计 |
|---|---|---|---|---|
| 1 | `端到端跑一条 SEC-SAXS 系列（图像→报告，全程 RAW）：先认洗脱峰，≥2 个峰就逐峰建子目录、各用本峰的` | 6/9（正面 2/5 · 诱饵 4/4） | 6/9（正面 2/5 · 诱饵 4/4） | 12/18 |
| 2 | `端到端跑一条 SEC-SAXS 系列（图像→报告）：多峰逐峰分析、逐帧归一化（BL19U2 header txt）` | **8/9**（正面 4/5 · 诱饵 4/4） | 6/9（正面 2/5 · 诱饵 4/4） | 14/18 |
| **3（定稿）** | `端到端跑一条 SEC-SAXS 系列：多峰逐峰分析、逐帧归一化、裁剪区视频、多区间 Guinier 表、IFT、分` | 7/9（正面 3/5 · 诱饵 4/4） | 7/9（正面 3/5 · 诱饵 4/4） | 14/18 |

逐题判定（✓ = 与金标一致）：

| # | 金标 | 轮 1 A/B | 轮 2 A/B | 轮 3 A/B |
|---|---|---|---|---|
| 1 端到端 | 本 skill | ✓/✓ | ✓/✓ | ✓/✓ |
| 2 Guinier 取点 | assess-guinier | ✓/✓ | ✓/✓ | ✓/✓ |
| 3 逐帧归一化 | 本 skill | ✗pro/✗pro | ✓/✗pro | ✗pro/✗pro |
| 4 基线 | correct-baseline | ✓/✓ | ✓/✓ | ✓/✓ |
| 5 多区间拟合 | 本 skill | ✗g/✗g | ✗g/✗p | ✗g/✗g |
| 6 **多峰逐峰** | 本 skill | ✓/✓ | ✓/✓ | ✓/✓ |
| 7 未解析重叠 | deconvolve | ✓/✓ | ✓/✓ | ✓/✓ |
| 8 视频 | 本 skill | ✗raw/✗none | ✓/✗raw | ✓/✓ |
| 9 珠模评估 | evaluate | ✓/✓ | ✓/✓ | ✓/✓ |

（`pro` = `process-sec-saxs-series`，`g` = `assess-guinier-fit-quality`，`p` = `compute-and-validate-p-of-r`，
`raw` = `assess-saxs-raw-data-quality`。）

### 定稿理由与**接受的代价**

- **新功能题 #6（多峰）三轮两票全中**，新诱饵 #7（未解析重叠峰）三轮两票全对 → 这次改造的**路由意图成立**。
- **4 条诱饵在 3 轮 × 2 评测者 = 6 次判定中 24/24 全对**：本 skill 从不抢兄弟 skill 的活（零误召）。
- 第 3 轮与第 2 轮同分（14/18），但第 3 轮**两评测者 9 题逐题完全一致**（第 2 轮有 2 处分歧）——
  按仓库既有口径**可复现优先**，取第 3 轮头部。第 3 轮同时把六类钩子塞进 57 字符：
  `多峰 / 逐峰 / 逐帧归一化 / 裁剪区视频 / 多区间 Guinier / IFT`。
- **代价（不再调参，按"两三轮即停"记录）**：正面 #3（逐帧归一化）与 #5（多区间拟合）在这一版稳定漏判，
  两票都落给**最像的兄弟 skill**（`process-sec-saxs-series` / `assess-guinier-fit-quality`）。
  归因：57 字符窗口里"端到端"与具体机制词互相挤占；这两题本身也确实是**判据 skill 合理可答**的问法
  （#5 问的是"合不合适"，判据版有权接）。这是明确的取舍，不是未发现的缺陷。

---

## 3. 真实数据实跑证据（本机，4 次运行全部 exit 0）

数据：`~/Repositories/DataProcess_2026.10.01`（`data/SEC-SAXS/<样品>`，cfg `data/20261001.cfg`），
解释器 `/Applications/BioXTASRAW/bin/python`，ATSAS `/Applications/ATSAS-4.1.4-1/bin`。
产物根一律在 `processed/SEC-SAXS/_multipeak-test/`。

| 运行 | 命令要点 | 结果 |
|---|---|---|
| **A** 全流程/双峰 | `4EH2-KDPV-ZN`（1500 帧），`--n-models 2` | 认到 **2 峰**（apex 812 / 957）→ `peaks/peak01_apex00812/`、`peaks/peak02_apex00957/` 各一套完整产物（profiles/tables/ifts/models/reports/series + 8 节模板 README）；顶层 README 为索引 |
| **B** 只认峰 | `--steps integrate,peaks` | 只出 `series/sec_peaks.png`、`tables/sec_peaks.csv` + 索引 README；`run_meta.json` 带 `peaks_only=true` |
| **C** 单峰 + `--multi-peak always` | `4LI2-676`（1800 帧） | 认到 **1 峰**（apex 904、窗 888–933、45 帧半宽、59σ）→ 仍建 `peaks/peak01_apex00904/` |
| **D** 人工区间（视觉回退路径） | `4EH2-KDPV-ZN`，`--peak-ranges "801,832;943,977"` + `--peak-buffers` | 峰窗与 buffer **逐字按人工给的**落地；每峰一份 138 行共享模板 README + 顶层索引 |

**四条系列的判定锚点**（2026-10-02）：

- `4LI2-676` → **单峰**（走原流程，行为不变）；
- `4EH2-KDPV-ZN` → **两个峰**：主峰 Rg 24.3 Å / MW(Vc) 37 kDa、二峰 Rg 14.0 Å / MW 7 kDa，谷底比 0.00、经典 R=4.45（🟢 基线分离）—— 主峰+二聚体那一类；
- `4DH2-676-apo-3` → 前肩 0.46× 主峰（标为**必须看图**）；
- `bsa` → 两侧 buffer 水平差 ≈0.6× 主峰（漂移主导，先基线校正再谈分不分峰）；
- `4DH1-KDPV-ZN` → 16 帧无束流帧（已剔除，仍要求人工看图）。

| **E** 逐帧参数的可用帧判定（`--frame-flag-q`） | `4EH2`（同 A 的两峰） | 修前 0/1500 帧有值（整列 -1 哨兵）→ 修后 P1 **112 帧**、P2 **110 帧**有值，且各 34/35 帧落在本峰窗内；P1 Rg 中位 24.0 Å 与其 Guinier 24.3 Å 相符 |

---

## 4. 本轮改代码踩到 / 修掉的坑

1. `SECM.averageFrames` 之外，逐峰子目录只建了 `profiles/` → `tables/frame_params.csv` 写入时
   `FileNotFoundError`（**已修**：先建齐 `profiles/tables/ifts/series/models/reports`）。
2. 谷底比判据两处单位不一致（峰高做了归一化、谷底没有；且谷底取的是未平滑曲线）→ 现在同除 `top`
   且在平滑曲线上取谷底（**已修**）。
3. `find-sec-peaks` 的 zoom 图 x 轴被 buffer 阴影撑开（**已修**）。
4. 图内中文/emoji 在本机 matplotlib 字体下缺字形 → 图标签一律英文、控制台/README 中文（**已修**）。
5. 逐峰流程中某个峰失败会整批中断 → 逐峰 `try/except` 隔离，失败的峰写进告警（**已修**）。
6. **人工区间（`--peak-ranges`）模式下仍打印"阈值 prominence≥0.0 / 半高宽≥1 帧"**（那是内部算曲线用的
   参数），会误导读者。已改成"峰窗口是**人工指定**的：…（不经过阈值判定；下列幅度/信噪比只是事后参照值）"。
   这条是运行 **D** 时看到的（**已修**，`d3829b7`）。

7. **逐峰 `frame_params.csv` 整列 `-1`（静默算不出来）**：RAW 的可用帧判定（`SECM.subtractAllSASMs`
   按**总强度** > 1.02×buffer 打标记）+ `window_size=5` 的连续窗要求 → 本机 `4EH2` 只有 69/1500 帧被标记
   且都是散落噪声帧 → 一个窗都凑不齐 → `run_secm_calcs` 把 rg/i0/vc/vp 全写成 `-1`（"没算"与"算失败"
   同形）。**已修**：新增 `--frame-flag-q`（默认 `0.01,0.05`，用**低 q 窗口积分强度**当判定口径，
   与认峰同窗）→ 同数据 P1 0 → **112 帧**有值、34 帧落在本峰窗内，Rg 中位 24.0 Å（Guinier 24.3 Å 对得上）；
   并在日志/`ranges.json` 里显式报告"有值 N/总帧数"，全 -1 时直接告警。
8. **总览图下三格（Rg/I(0)/MW）画成"空轴"**：逐峰参数全 `-1` 被掩成 NaN 后，`plot` 什么都没画、
   自动量程给 ±0.05 → 面板看着像"值 ≈ 0"。**已修**：只在"真的有有限值"时才画，否则在面板里写明
   `no per-frame values (RAW computed none for this peak: all -1)`。这条是**看图**看出来的
   （跑完流程截图逐格核），不是读代码读出来的。

另记两条**没改**但有据的判断：

- **总强度（`getIntI`）不能用来认峰**：BSA 2000 帧总强度只有 ~10% 起伏、被束位漂移主导，按 RAW 口径能
  找出 **22 个假峰**（`probe_detect.py` 实测）。改用**扣减后低 q 窗口 [0.01,0.05] 1/Å 的积分强度**。
- **阈值只能靠"相对主峰幅度"卡**：扣减后色谱图带**系统性纹波**（10–20% 主峰、7–15σ、宽 10–18 帧，
  且多个 q 窗口一致 → 任何"相对噪声/多窗口一致性"判据都杀不掉）。故默认 `prominence≥0.2×(主峰)`
  而非更松的 0.05（松阈值下 `4EH2-KDPV-ZN` 会报 6 个峰：4 真峰 + 2 纹波）。

---

## 5. 遗留 / 不做

- 正面 #3/#5 的漏判（见 §2 代价）——按"两三轮即停"记录，不再调参。
- 谷底**未回基线**的未解析重叠峰**不**在流水线里硬切区间 → 转 `deconvolve-overlapping-elution-peaks`。
  这是设计边界，不是遗留缺陷。
