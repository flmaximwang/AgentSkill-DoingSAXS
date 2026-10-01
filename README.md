# AgentSkill-DoingSAXS

"做 SAXS 这件事本身"的 skill 包：**拟合之前先判数据**，出现怪现象时怎么把责任判给数据、几何还是模型，
**把一条 SEC-SAXS 系列或一批管式/静态帧端到端跑完**（图像 → 报告）并交付一份人能看懂的结果说明，
**把高分辨模型嵌进从头算重建**（珠模型 + DENSS 电子云）并报出可解释的分数，
以及**用 ATSAS 的 crysol 从 model 算曲线、与 `processed/` 里已有的拟合同口径比较**（datcmp + 同网格 χ² + 图 + README）。
（单点判据 skill 在 `AgentSkill-UsingBioXTASRAW`；这里放**评估/归因** + **端到端流水线** + **交付物/对齐** + **嵌入** + **模型↔已有拟合的比较**。）

| 来源 | 内容 | 沉淀成 |
|---|---|---|
| 2026-10-01 管式批次（BL19U2，11 个样品：A5-05-1…6 稀释序列 + 5705/877-apo/877-4zinc/97df/BSA）的实跑 | 逐帧质量（漂移/离群/对比度/SNR/误差诚实度）、低 q 上翘的三道归因检验、beamstop 几何与 q_min 依据 | [`skills/assess-saxs-raw-data-quality/`](skills/assess-saxs-raw-data-quality/) |
| 同批 A5-05-1…6 的 `embed/` 实跑（GNOM→DAMMIF×15→DAMAVER + CIFSUP；ChimeraX fitmap 的 R 阶梯）+ ATSAS CIFSUP/SUPCOMB 与 ChimeraX fitmap 官方文档 | 珠模型嵌入（CIFSUP/NSD，含分数只写在输出文件里、selection 口径、系综波动）与电子云嵌入（fitmap 的 resolution/search/seed 与姿态散布判据），以及 PyMOL 出图（ChimeraX `--nogui` 存不了图） | [`skills/embed-a-model-in-bead-and-density-models/`](skills/embed-a-model-in-bead-and-density-models/) |
| ATSAS 4.1.4 官方手册 *CRYSOL*（命令行/输出/常见报错）+ *DATCMP*（标准化残差 / 三个统计检验 / Assessing Model Fit）+ `.dat`/`.fit`/`.fir` 格式页；以及 2026-10-01 批次上 `crysol`/`datcmp`/`datcrop` 的实跑（4DH2-676-apo-3、A5-05-1、4LI2-676 三条链） | 模型→理论曲线→拟合数据这条链上「各家 χ² 不是同一个数」的全套口径：分母（N vs N−1）、q 网格（重采样）、误差列、`datcmp` 单文件用法与哨兵值、判档阈值的确切来源 | [`skills/compare-model-curve-with-existing-fits/`](skills/compare-model-curve-with-existing-fits/) |
| 两条流水线的 README 各写各的（章节不同、同一个数位置不同、采用口径不一致）这一实际问题 | 抽出**唯一**一份结果目录契约与实现：8 节固定模板 + 字段别名表 + 「采用区间 / 采用那一支 P(r) / 三档可用性」统一口径 + 四道验收（结构 / 章节 / README↔产物 / 产物↔产物） | [`skills/write-saxs-results-readme/`](skills/write-saxs-results-readme/) |
| **多峰 SEC 数据被静默丢峰**（RAW 的 `findSampleRange` 只取 `argmax` 那个峰）+ 本机 5 条 SEC 系列的实跑（4EH2-KDPV-ZN 确为两个组分：主峰 Rg 24 Å / MW 37 kDa 与二峰 Rg 14 Å / MW 7 kDa） | SEC 流水线补上**认峰 → 逐峰建子目录 → 各用本峰的 buffer 与峰窗分别扣减分析**，并把「代码判不干净」的部分做成**视觉复核的一步**（判定图 + 阈值扫描图 + `--peak-ranges` 覆盖） | [`skills/run-a-sec-saxs-pipeline-end-to-end/`](skills/run-a-sec-saxs-pipeline-end-to-end/)（`scripts/sec_peaks.py`、`scripts/find-sec-peaks.py`） |

**5 个 skill**：1 个评估/归因原子技能（24 个候选判据 → 通过 1 个）+ 2 条端到端流水线（SEC 连续洗脱帧 / 管式静态帧，
**共用同一份结果目录契约与 README 实现**）+ 1 个"模型嵌入重建"的可执行流程 + 1 条"模型↔已有拟合"的同口径比较流程。
评估/归因那条的流程：真实数据上磨出来的一条流程——
**先量几何（中心/掩膜/q_min）→ 逐帧 QC → 扣减后形状 → 低 q 上翘归因（空白-空白 / 背景形状失配 / 2D 差分）→ 完整结论（排除了什么、还剩什么、下一步做什么实验）**。
它**不做任何拟合**：拟合归 `AgentSkill-UsingBioXTASRAW` 的 15 个 skill，这里只判"那些拟合肥不肥"。

> **状态**：已 push 到 <https://github.com/flmaximwang/AgentSkill-DoingSAXS>（**public**，远程用 SSH，default=main）；
> 6 个 skill 均已安装进 default profile 的 **`saxs`** 类目（三段式标识符 + skills.sh/community，pin 见下方表格）。
> `compare-model-curve-with-existing-fits`（2026-10-02）的 scan verdict 也是 **CAUTION**：`scripts/model-vs-fit.py` 命中
> skills-guard 的 `python_os_environ`（读 `$ATSAS` 定位 ATSAS，2×HIGH）与 `python_subprocess`（调本机 crysol/datcmp/datcrop，1×MEDIUM），
> 无网络、无外发，按先例用 `--force` 安装。
> `embed-a-model-in-bead-and-density-models` 的 scan verdict 是 **CAUTION**：`scripts/embed-model.py` 命中
> skills-guard 的 `python_os_environ` 与 `python_subprocess` 两条规则（2×HIGH exfiltration on `os.environ`、
> 3×MEDIUM execution on `subprocess`），因此按既有先例用 `--force` 安装。命中项的实质是**本机可执行文件调用**
> （ATSAS `cifsup`、ChimeraX、PyMOL），无网络、无外发；安装后 `diff -rq` 与仓库逐字节一致。

安装记录（pin = 安装时仓库的 commit，`hermes skills check` 用它对账）：

| skill | 安装 pin | scan verdict | 安装时内容哈希 |
|---|---|---|---|
| `assess-saxs-raw-data-quality` | `eff1219` | safe | `sha256:a9d3d37700dc7b9f` |
| `run-a-sec-saxs-pipeline-end-to-end` | `cf55322dcc03d05ce465ea7e97646e6a47d05d46` | safe（`sec_peaks.py` 命中 `string_reversal` 的 LOW 提示，实为 `[::-1]` 排序，判 ALLOWED） | —（本轮新增 `sec_peaks.py`/`find-sec-peaks.py`、删 `results_readme.py`/`write-results-readme.py`，哈希口径见 `hermes skills check`） |
| `run-a-tube-saxs-pipeline-end-to-end` | `4fd536d` | safe | `sha256:8724528379778d13` |
| `write-saxs-results-readme` | `f45c923` | safe | —（本轮把 `test-results.md` 加进技能目录，哈希口径见 `hermes skills check`） |
| `embed-a-model-in-bead-and-density-models` | `ac9ad19` | caution | `sha256:d73f8356a2578300` |
| `compare-model-curve-with-existing-fits` | `17feeba` | caution | `sha256:fef50a862f8dfd97` |

## 索引

| skill | 用途 | 可执行入口 |
|---|---|---|
| [run-a-sec-saxs-pipeline-end-to-end](skills/run-a-sec-saxs-pipeline-end-to-end/SKILL.md) | **端到端跑一条 SEC-SAXS 系列**（全程 RAW，不自写拟合）：**先认洗脱峰（≥2 个峰就逐峰建子目录、各用本峰的 buffer 与峰窗分别扣减分析；RAW 只认最大那个峰，原流程会静默丢峰）** → 逐帧归一化（补 BL19U2 header txt / 读线站已有 txt）→ 归一化裁剪视频 → buffer/sample 区与扣减 → 多区间 Guinier → IFT → MW → 形状重建（DENSS 电子云 + ATSAS DAMMIF 珠模 + DAMAVER）→ RAW PDF 报告 → **按 `write-saxs-results-readme` 的模板写结果目录 `README.md`**（多峰时顶层是索引、每峰子目录各一份），每个节点都落 `.dat`/表格 | `references/bl19u2-header-normalization.md`、`scripts/`（5 个脚本：emit-bl19u2-header-txt / crop-video-normalized / run-raw-sec-pipeline / **sec_peaks** / **find-sec-peaks**） |
| [run-a-tube-saxs-pipeline-end-to-end](skills/run-a-tube-saxs-pipeline-end-to-end/SKILL.md) | **端到端跑一批管式/静态帧**（目录里样品 run + 夹着它的 control run）：逐帧归一化（cfg 两处开关）→ control 相对缩放（高 q 窗）→ 扣减 → 多区间 Guinier（14 列判据 + 四条闸门）→ BIFT 两次起跑 + GNOM 的 P(r)（三条闸门）→ MW → 电子云 + 珠模 + DAMAVER → 报告/workspace → 四联诊断图 `qc.png` + 同一套 README；批处理另给 `summarize-tube-run.py` / `plot-tube-overview.py` | `references/tube-control-pairing-and-scaling.md`、`scripts/`（3 个脚本：run-raw-tube-pipeline / summarize-tube-run / plot-tube-overview） |
| [write-saxs-results-readme](skills/write-saxs-results-readme/SKILL.md) | **结果目录的验收 + README**：两条流水线共用的唯一实现（`readme_common.py`）——8 节固定模板、字段别名表（`rg`↔`Rg`、`q_min`↔`qmin`、`mw`↔`MW_kDa`、`detail1..4`↔`aux`）、统一选值口径（采用区间 / 采用那一支 P(r) / 三档可用性）、哨兵值翻译（`rg=-1`/`inf`/`nan`）；`verify-results-folder.py` 四道验收（结构 / 章节 / README↔产物关键数字 / 产物↔产物去重矛盾 + 时间戳） | `references/results-folder-contract.md`、`references/verification-rules.md`、`scripts/write-readme.py`、`scripts/verify-results-folder.py` |
| [assess-saxs-raw-data-quality](skills/assess-saxs-raw-data-quality/SKILL.md) | 一批 SAXS 原始帧的质量评估 + 低 q 上翘归因：先量几何与掩膜（含 RAW 读 Pilatus 的 **y 翻转**、beamstop 边缘 → `--qmin`），再逐帧判据表（对比度/漂移/离群/误差诚实度/SNR-qmax），再做三道归因（**空白-空白可复现极限**、**shape×电平**、**2D 差分：各向同性光晕 vs 紧贴 beamstop 的窄亮环**），最后按"浓度标度律"写出合格结论；含"单条曲线上翘 ≠ 相互作用"的判据与下一步实验设计 | `references/bl19u2-geometry-and-mask.md`、`references/upturn-attribution-protocol.md`、`scripts/`（5 个可执行脚本 + 1 个共用件） |
| [compare-model-curve-with-existing-fits](skills/compare-model-curve-with-existing-fits/SKILL.md) | **把 model 用 CRYSOL 算成 SAXS 曲线并拟合到扣减后的 `.dat`，再与 `processed/` 里已有的拟合（DENSS `.fit` / DAMMIF `.fir` / GNOM `.out`）同口径比**：`crysol --constant`（常数扣减＝缓冲液失配量级）+ `datcrop` 先裁可用 q 区间 + `datcmp <单个 .fit/.fir/.out>`（ATSAS 统一的 red.χ² 与 CorMap，**判决只用这列**）+ 同网格 χ² 交叉核对（标出哪条是重采样网格）+ 三列 χ² 为什么不是一个数（分母 N vs N−1、网格、误差列）+ 判档阈值（χ²_n 95% 区间 + CorMap α=0.01）+ 叠加图/残差图/明细表/README | `references/crysol-datcmp-parameters-and-fits.md`、`scripts/model-vs-fit.py`、`scripts/plot-model-vs-fit.py`、`scripts/write-model-fit-readme.py` |
| [embed-a-model-in-bead-and-density-models](skills/embed-a-model-in-bead-and-density-models/SKILL.md) | **把高分辨模型（PDB）嵌进珠模型与 DENSS 电子云**：珠模型分支 = GNOM→DAMMIF×N→DAMAVER 取 `damfilt` 共识模型 + CIFSUP（`--method=NSD`，REGRID＝SUPCOMB fast 口径）给 **NSD**（分数只在输出文件里，附 selection 对照与系综波动）；电子云分支 = ChimeraX fitmap（**必须 `resolution` 才有 correlation、必须 `search+seed` 才是全局解**）跑 R 阶梯并报**姿态散布**（实测 20 Å 下存在相距 28.9 Å 的等分解）；另含 CA 几何复核（CA→最近 bead 距离）、两张分坐标系图（PyMOL，ChimeraX `--nogui` 存不了图）、每个样品一份 `embed/README.md` | `references/bead-and-density-embedding-parameters.md`、`scripts/embed-model.py`、`scripts/render-embed-figure.py`、`scripts/write-embed-readme.py` |

## 多峰 SEC-SAXS（≥2 个洗脱峰）

RAW 的 `SASCalc.findSampleRange()` 只取 `argmax(peak_heights)` 那个峰，所以一条有两个组分的 SEC 系列
用原流程跑出来**只有最大峰的那一条曲线**，第二个峰**静默消失**（产物看起来完全正常，只是少了东西）。
SEC 流水线现在默认 `--multi-peak auto`：认峰（扣减后低 q 窗口积分的色谱图 + 滚动中位基线 +
RAW 的 `savgol`/`find_peaks` 口径）→ ≥2 个峰就 `peaks/peakNN_apexNNNNN/` **每峰一套完整产物**
（各自用"紧邻本峰的前后两段 buffer + 本峰峰窗"重跑扣减与下游分析），顶层 `README.md` 是索引。

"代码判不干净"的部分做成**视觉复核的一步**：

- 判定图 `series/sec_peaks.png`（峰窗 + buffer 段 + 每峰自己的 Rg/I(0)/MW）与 `series/sec_peaks_zoom.png`（逐峰放大）；
- `needs_visual_check` 把"刚过阈值 / 两峰没分开（谷底比与经典 R 两个判据）/ buffer 落在谷底 /
  两侧 buffer 水平差 ≥0.15×主峰 / 无束流帧"写成**必须看图**的条件，命中就写进 README 的
  「看图之前别引用数字」一节，并从中给出下一步（`--baseline`、`--peak-buffer-max`、换 buffer 段）；
- `find-sec-peaks.py`（秒级、不重跑图像积分）改参数重看、`--sweep` 出阈值扫描图、
  `--peak-ranges/--peak-buffers` 把眼睛的结论写回主管线。
- **逐帧参数（Rg/I0/MW）的"可用帧判定"另有开关**：RAW 默认按总强度判定，SEC 数据上会被束位/通量漂移
  主导（实测 1500 帧只标记 69 个散落噪声帧 → 逐峰 `frame_params.csv` **整列 -1**，且不报错）。
  流水线默认改用低 q 窗口积分强度（`--frame-flag-q "0.01,0.05"`，与认峰同一窗口），同一数据 P1
  从 0 帧变 112 帧有值；要回到 RAW 默认口径传 `--frame-flag-q none`。

边界：两峰之间**谷底没回到基线** = 未解析重叠，转 `deconvolve-overlapping-elution-peaks`（SVD/EFA/REGALS），
不要在流水线里硬切区间。

实测锚点（本机 5 条 SEC 系列，2026-10-02）：`4LI2-676` 单峰；`4EH2-KDPV-ZN` **两个峰**
（主峰 Rg 24.3 Å / MW 37 kDa、二峰 Rg 14.0 Å / MW 7 kDa，R=4.45 基线分离 —— 主峰+二聚体那一类）；
`4DH2-676-apo-3` 前肩 0.46× 主峰（标为必须看图）；`bsa` 两侧 buffer 水平差 ≈0.6×主峰（漂移主导，
先做基线校正再谈分不分峰）；`4DH1-KDPV-ZN` 有 16 帧无束流帧。

## 结果目录的交付契约（`write-saxs-results-readme`）

两条流水线**不再各自写 README**：它们最后一步调同一份实现，产出同一套 8 节模板
（结论速览 → 最终拟合参数明细表 → 关键结果 → 文件清单 → 判据表 → 没做的/不能信的 →
复现参数 → 想自己复核），并用 `verify-results-folder.py` 验收。**判据文本、章节顺序、阈值只有这一份定义**，
改格式去那个 skill 改（`references/results-folder-contract.md` 是它的文字版）。

```bash
python <write-saxs-results-readme>/scripts/write-readme.py <结果目录>          # 补写/重写 README.md（不重算）
python <write-saxs-results-readme>/scripts/verify-results-folder.py <结果目录>  # 退出码 1 = 有 🔴
```

## 质量凭据（盲测）

按 Hermes 路由时的**真实 57 字符截断**（`agent/skill_utils.py:761`）把候选 skill 的 description 交给**独立评测者**逐条判路由
（题目 10 条 = 6 正面 + 4 诱饵 · 交错排列 · 评测者看不到 SKILL.md 正文）：

| 轮次 | 对象 | 评测者 A | 评测者 B | 处置 |
|---|---|---|---|---|
| 1（首轮） | `embed-a-model-in-bead-and-density-models`，17 个候选，10 条 | **9/10**（正面 5/6 · 诱饵 **4/4**） | **9/10**（正面 5/6 · 诱饵 **4/4**） | 两评测者 10 题**逐题完全一致**；唯一错误 #6（ChimeraX `--nogui` 存图报 OpenGL）两票都判 *none*——归因为**可见窗口放不下**（57 字符已给主触发词），按既有口径记代价、不雕题面 |
| 1 | `compare-model-curve-with-existing-fits`，21 个候选，10 条 | **9/10**（正面 5/6 · 诱饵 **4/4**） | **9/10**（正面 **6/6** · 诱饵 **4/4**） | 两条错误都落在与 `fit-a-high-resolution-model-to-data` 的边界上（诱饵 #7「RAW 界面的 pdb 没有 Chi² / 要不要加 harmonics」被抢走；A 另把正面 #2 判给对方）。**试改一版头部（把「理论曲线」换成「命令行」）后诱饵 #7 被抢回、正面却掉 3–4 条（6/10、7/10）→ 回退用首轮头部**；取舍与代价见该 skill 的 `test-results.md` |
| 1（首轮） | `write-saxs-results-readme`，13 个候选，9 条 | **9/9**（正面 **5/5** · 诱饵 **4/4**） | **9/9**（正面 **5/5** · 诱饵 **4/4**） | 两评测者**逐题完全一致、零错误**：两条最像的诱饵（#9 也是「写 README」但属批处理根目录 → `result-folder-readme`；#8 README 里也写 NSD 但判据归 `evaluate-a-shape-reconstruction`）都没被误收；无事可调，按既有口径记录即止 |
| 2 | 同上，改了头部重测 | 6/10（正面 2/6 · 诱饵 4/4） | 7/10（正面 3/6 · 诱饵 4/4） | 已回退（见上一行处置） |
| 1（首轮） | `assess-saxs-raw-data-quality`，17 个候选，10 条 | **7/10**（正面 4/6 · 诱饵 **4/4**） | **7/10**（正面 4/6 · 诱饵 **4/4**） | 4 条诱饵 100% 未误收（流水线/Guinier 判据/P(r)/MW 各归各位）；三条一致错误里 #10 认定为**金标偏严**（修订后 8/10 · 8/10），#5（砍 q_min 砍到哪）与 #9（上机前排 control）属**可见窗口放不下** → 按既有口径记代价、不再调参 |
| 1 | `run-a-sec-saxs-pipeline-end-to-end`（**多峰改造后**，21 个候选，9 条） | 6/9（正面 2/5 · 诱饵 **4/4**） | 6/9（正面 2/5 · 诱饵 **4/4**） | 新功能题 #6（两个洗脱峰各要一套）与新增诱饵 #7（未解析重叠峰）**两票都对**；但新头部把「逐帧归一化」挤出了 57 字符窗口 → 正面 #3（逐帧归一化）、#8（视频）两票丢给兄弟 skill → **改头重测** |
| 2 | 同上，头部改为「（图像→报告）：多峰逐峰分析、逐帧归一化（BL19U2 header txt）」 | **8/9**（正面 4/5 · 诱饵 **4/4**） | 6/9（正面 2/5 · 诱饵 **4/4**） | #3、#8 在 A 侧修回，但两评测者出现 2 处分歧（不稳定）、#5 两票仍丢 → 再试一版把 Guinier/视频也塞进窗口 |
| 3（定稿） | 同上，头部改为「端到端跑一条 SEC-SAXS 系列：多峰逐峰分析、逐帧归一化、裁剪区视频、多区间 Guinier 表、IFT、分」 | **7/9**（正面 3/5 · 诱饵 **4/4**） | **7/9**（正面 3/5 · 诱饵 **4/4**） | 与第 2 轮同分（14/18）但**两评测者 9 题逐题完全一致**（可复现优先）；六类钩子都在窗口内；两条漏判（#3、#5）都落给**最像的兄弟 skill**（判据版有权接这两问）；**诱饵 3 轮×2 评测者 24/24 全对**（零误召）→ 按「两三轮即停」定稿，取舍与代价见该 skill 的 `test-results.md` |

darwin 优化：`compare-model-curve-with-existing-fits` 2026-10-02 走完 9 维基线（84.3）→ 三轮 paired 3-0 keep（dim3 if-then 故障表 / dim4 六道 🔴 STOP 闸门 / dim8 可粘命令+三段式回答，另修一处复审发现的示例不一致）→ 92.7，记录见该 skill 的 `test-results.md` 第 5 节。

细节（含每条错误归因与金标修订理由）见 [`test-results.md`](test-results.md)；
`embed-a-model-in-bead-and-density-models` 的首轮盲测见 [`skills/embed-a-model-in-bead-and-density-models/test-results.md`](skills/embed-a-model-in-bead-and-density-models/test-results.md)。
`write-saxs-results-readme` 的首轮盲测（含 13 个候选的截断文本与逐题判定）见 [`skills/write-saxs-results-readme/test-results.md`](skills/write-saxs-results-readme/test-results.md)。
`run-a-sec-saxs-pipeline-end-to-end` 的多峰改造（含 4 次真实数据实跑 + 三轮盲测逐题矩阵）见 [`skills/run-a-sec-saxs-pipeline-end-to-end/test-results.md`](skills/run-a-sec-saxs-pipeline-end-to-end/test-results.md)。

## 安装（三段式标识符，按仓库内路径，不需要 tap；`--category` 只决定落点）

```bash
for s in assess-saxs-raw-data-quality run-a-sec-saxs-pipeline-end-to-end \
         run-a-tube-saxs-pipeline-end-to-end write-saxs-results-readme \
         embed-a-model-in-bead-and-density-models \
         compare-model-curve-with-existing-fits; do
  hermes skills install <owner>/AgentSkill-DoingSAXS/skills/$s --category saxs -y
done
```

更新：`hermes skills update <skill 名>`（改了本仓库并 push main 之后）。

## 运行前提

- 脚本用 **RAW 自带的解释器**跑：`/Applications/BioXTASRAW/bin/python`（含 pyFAI/numba/matplotlib）。
  管式/嵌入那两条另需 **ATSAS ≥ 4**（本机 `/Applications/ATSAS-4.1.4-1`，`--atsas-dir` 指到 `bin` 级；
  建议 `ln -s /Applications/ATSAS-4.1.4-1 ~/ATSAS-4.1.4-1` 让 RAW 自己找得到）。
- 不要在 RAW 源码目录里跑（`sascalc_exts` 会被遮蔽）。
- `compare-model-curve-with-existing-fits` 另需：**ATSAS ≥ 4**（`crysol` 算/拟合、`datcmp` 评拟合、`datcrop` 裁 q；
  脚本按 参数 → `$ATSAS` → 常见路径 自动找 bin 目录）。
- `embed-a-model-in-bead-and-density-models` 另需：**ChimeraX**（fitmap，本机 1.9）、**PyMOL**（出图，
  `/Applications/PyMOL.app/Contents/bin/python3.10`）。
- SEC 输入 = 一条连续洗脱的 tif 系列（可能只有 `.Iochamber` + `.log`）；管式输入 = "每样品一个目录、
  目录里同时有样品帧与夹着它的 control 帧"的布局（BL19U2 批处理风格，命名形如 `A5-05-1_0013_00001.tif`，
  run 键 = 前两段）。
- `write-saxs-results-readme` 与两条流水线是**兄弟技能**（同类目下的同名目录）：管线按这个相对位置找它，
  也可以 `--readme-script` 显式指定；找不到时管线只告警、不写 README。

## 与其它包的边界

| 问题 | 去哪 |
|---|---|
| 这批原始帧好不好 / 上翘是真还是假 / q_min 取多少 | **本包**（`assess-saxs-raw-data-quality`） |
| **SEC 连续洗脱帧**图像→曲线→Rg/P(r)/MW/3D 端到端跑一遍 | **本包**（`run-a-sec-saxs-pipeline-end-to-end`） |
| **管式/静态帧**端到端跑一遍 | **本包**（`run-a-tube-saxs-pipeline-end-to-end`） |
| 结果目录缺不缺东西 / README 怎么写、怎么验收、两条流水线怎么对齐 | **本包**（`write-saxs-results-readme`） |
| **把高分辨模型嵌进珠模型 / DENSS 电子云并报分**（NSD / fitmap correlation） | **本包**（`embed-a-model-in-bead-and-density-models`） |
| **用 model 算理论曲线、跟 `processed/` 里已有的拟合比 χ²** | **本包**（`compare-model-curve-with-existing-fits`） |
| RAW 界面里的判据：默认计算器选谁、minimal 曲线陷阱、harmonics / N samples 怎么调、某个 χ² 该不该信 | `AgentSkill-UsingBioXTASRAW` 的 `fit-a-high-resolution-model-to-data`（本包的命令行版本会引用它的判据） |
| Guinier 怎么取点、P(r) 用哪个程序、MW 用哪一法、重建怎么评、模型配不配数据 | `AgentSkill-UsingBioXTASRAW` 的 13 个判据 skill |
| ATSAS 命令行（GNOM/DAMMIF/DATMW…） | `AgentSkill-UsingATSAS` |
