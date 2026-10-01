# 质量记录（compare-model-curve-with-existing-fits）

## 1. 盲测（路由）

**方法**：候选集 21 个 = profile `saxs/` 类目现存的 19 个（含 13 个判据 skill、2 条端到端流水线、`assess-saxs-raw-data-quality`、`embed-…`）+ `mals/` 的 `analyze-saxs-data-with-atsas` + 本 skill。
**可见信息**：每个候选只给 `name :: description 的前 57 个字符`（Hermes 路由时真实的截断：`agent/skill_utils.py` `SKILL_PROMPT_DESC_LIMIT=60`，渲染为 `desc[:57]+"..."`）；评测者**看不到 SKILL.md 正文**。
**评测者**：每轮 2 位独立子代理（互不可见、同 prompt、同候选表、不看答案），每条请求必须给唯一 skill 名 + 置信度。
**题目**：10 条（**6 正面 + 4 诱饵**，交错排列，不告知类型）→ 见 `test-prompts.json`。

| # | 类型 | 请求（摘要） | 金标 | 轮1-A | 轮1-B | 轮2-A | 轮2-B |
|---|---|---|---|---|---|---|---|
| 1 | 正面 | 有扣减好的 .dat + 二聚体 cif，算曲线拟合看配不配 | 本 skill | **本 skill** (h) | **本 skill** (h) | fit-a-high-res (h) ✗ | fit-a-high-res (h) ✗ |
| 2 | 正面 | processed 里已有 DENSS/DAMMIF，想知道模型是不是更贴 | 本 skill | fit-a-high-res (m) ✗ | **本 skill** (m) | **本 skill** (h) | **本 skill** (m) |
| 3 | 正面 | crysol 命令行参数怎么给（--constant/--lm/smax） | 本 skill | **本 skill** (m) | **本 skill** (h) | analyze-saxs-atsas (m) ✗ | **本 skill** (h) |
| 4 | 正面 | RAW 2.66 / ATSAS 日志 2.67 / 自己算 4.4，信哪个 | 本 skill | **本 skill** (m) | **本 skill** (l) | fit-a-high-res (m) ✗ | **本 skill** (m) |
| 5 | 正面 | datcmp 怎么用在 .fit/.fir/.out 上；χ²=0、p=1 是完美拟合吗 | 本 skill | **本 skill** (h) | **本 skill** (h) | **本 skill** (h) | **本 skill** (h) |
| 6 | 正面 | 4LI2 四条链拆开，哪条链最贴 | 本 skill | **本 skill** (m) | **本 skill** (h) | fit-a-high-res (m) ✗ | fit-a-high-res (h) ✗ |
| 7 | 诱饵 | RAW 里 Plot 的 pdb 没有 Chi²、低 q 差很多，要不要加 harmonics | fit-a-high-res | **本 skill** (l) ✗ | **本 skill** (m) ✗ | fit-a-high-res (m) ✓ | fit-a-high-res (m) ✓ |
| 8 | 诱饵 | A5.pdb 嵌进珠模型/DENSS 电子云，NSD 0.87 算好吗 | embed-a-model | embed (h) ✓ | embed (h) ✓ | embed (h) ✓ | embed (h) ✓ |
| 9 | 诱饵 | 15 个 DAMMIF 分 3 簇、平均 NSD 0.89，重建能用吗 | evaluate-a-shape-recon | evaluate (h) ✓ | evaluate (h) ✓ | evaluate (h) ✓ | evaluate (h) ✓ |
| 10 | 诱饵 | 管式原始帧一路跑到 Rg/P(r)/MW/DENSS + README | tube-pipeline | tube-pipeline (h) ✓ | tube-pipeline (h) ✓ | tube-pipeline (h) ✓ | tube-pipeline (h) ✓ |

**得分**

| 轮次 | description 头（可见 57 字） | 评测者 A | 评测者 B |
|---|---|---|---|
| 1 | `用 CRYSOL（ATSAS）把 model 算成 SAXS 理论曲线并拟合到扣减后的 .dat，再用 datcm…` | **9/10**（正面 5/6 · 诱饵 4/4） | **9/10**（正面 6/6 · 诱饵 **4/4**） |
| 2 | `用 ATSAS 的 crysol 命令行把 model 拟合到扣减后的 .dat，再用 datcmp 与 proc…` | 6/10（正面 2/6 · 诱饵 4/4） | 7/10（正面 3/6 · 诱饵 4/4） |

**采用轮 1 的头部**（总分 18/20 vs 13/20）。轮 2 的改法（把 `理论曲线` 换掉、前置 `命令行`）确实把诱饵 #7 抢了回来（两位评测者都改判给 `fit-a-high-resolution-model-to-data`），代价是丢掉了 3–4 条正面 —— 因为触发本 skill 的用户话术里正是「算理论曲线」「拟合到 .dat」，而 `命令行` 这个词把 #3 引到了 ATSAS 手册索引 skill。**两轮的取舍已记在下面，不再继续调题面。**

## 2. 与 `fit-a-high-resolution-model-to-data` 的互抢（已知代价，不修）

两条 skill 天生在同一个话题上：一条是 **RAW 界面里的判据与参数陷阱**，一条是 **命令行产物 + 与已有拟合的同口径比较**。可见 57 字里放不下这个区分：

- 轮 1 的错误方向：诱饵 #7（RAW 界面、没有 Chi²、要不要加 harmonics）被本 skill 抢走（A 低置信、B 中置信）；评测者 A 另把正面 #2（"模型是不是比 DENSS/DAMMIF 更贴"）判给了对方。
- 判据依据（与既有先例一致）：#7 的正确归属需要看到 `minimal 曲线陷阱 / harmonics / N samples` 这些词，而对方 skill 的可见头也只写到「计算即拟合：不要先生成 mini…」，同样不含 harmonics —— 属于**可见窗口放不下**，不是题面缺陷。
- **代价**：主题词相邻的两条 skill 之间可能多一跳（用户会被引到"计算即拟合"的判据 skill，再从那里被引回来）。风险低：两条都在同一包族里、都指向 ATSAS/CRYSOL，且两边正文互相点名。

## 3. 实跑验证（产物级，不是"跑通了"）

环境：ATSAS 4.1.4-1（`/Applications/ATSAS-4.1.4-1/bin`），解释器 `/Applications/BioXTASRAW/bin/python`，数据 `~/Repositories/DataProcess_2026.10.01`。四个场景全部实跑过、产物落盘。

| 场景 | 命令要点 | 实测结果 |
|---|---|---|
| SEC `4DH2-676-apo-3` + `4DH2-676-apo-dimer.cif` | `--processed processed/SEC-SAXS/4DH2-676-apo-3`（默认 `--constant`） | 模型 `datcmp` red.χ²=**1.267**（crysol 自报同为 1.267，CorMap C=41/p=0）；同批 DENSS **2.668**、DAMMIF×4 **4.25–4.33**、GNOM `.out` 4.253（口径不同，标"不参与排序"）。`--lm=100` 与 20 结果一致；`--no-constant` 4.233 |
| 管式 `A5-05-1` + `A5.pdb` | `--qmin 0.05666 --qmax 0.25915`（`datcrop` 先裁） | `A5.pdb` red.χ²=**155.0**（🔴，模型 Rg 15.5 Å vs 样品 20.8 Å）；未裁时为 161.6 → 裁 q 只降 6.6，主因不是 q 区间。同批 DAMMIF 1.22–1.48、DENSS 5.42、GNOM 0.943 |
| SEC `4LI2-676` + 三条链 A/B/C | 三个 `--model` 一次跑 | **2.320 / 2.478 / 2.671**（Rg 均 15.8–16.0 Å，与样品 16.2 Å 相符）vs 同批 DENSS 5.053、DAMMIF 4.86–4.90 → 高分辨模型全面更贴 |
| 失败路径 | `--model 4LI2-676_D.cif`（该文件里只有水） | crysol 报 `error: No such model or chain in structure`，脚本记 `ok=false`、README 顶部写「有 1 个模型没算出来」并附原始报错；同时 `--fit` 显式指定（不用 `--processed`）路径也已验证 |

**工具级坑（都已复核）**：

- `datcmp` 喂普通 4 列 `.dat`（`ifts/ift_fit.dat`）→ `χ²=0.000000, p=1.000000`（自己跟自己比，假"完美拟合"）；喂 BIFT `.ift` → `-1.000000`。脚本按"χ²≤0 或 (χ²=0 且 p=1)"判为不可判。
- χ² 分母实测：crysol/DAMMIF 用 `N−1`（复算 1.2668/4.334 与自报一致），DENSS 头里的 2.665 用 `N`（用 `N−1` 是 2.668 = datcmp 的值）。
- q 网格实测：原始 1150 点等距 3.8941e-4 Å⁻¹；DAMMIF 1149 点重采样成等距 3.89e-4（与原始点最大差 4.95e-7 Å⁻¹）；DENSS 845 点。把 DENSS 曲线插回原始网格会把 χ² 从 2.668 抬到 **4.40** → 因此"同网格"列只在与原始网格一致（`regridded=false`）时才与 datcmp 同级。
- `dammif_01.fit`（3 列、无误差列、q 从 0 开始）不是可判文件：脚本在同名 4 列 `.fir` 存在时自动跳过它。

## 4. 待办缺口（交给下一次）

- 手册没给"χ² 多大算好"的绝对阈值，本 skill 的 95%（χ²_n 分布）区间是手册"可由 χ²_n 算"那句话的实现，不是原文表格；若以后拿到官方或文献的绝对判据，应替换并在 `references/` 里记出处。
- `Anderson-Darling` 检验本 skill 只透传参数、未纳入判档（手册说它对尾部比 red.χ² 更敏感）——若以后出现"χ² 合格但尾部差"的实例，再决定是否加进判档。
- 系综/多构象的判据（同一条曲线对多个构象加权）不在本 skill 范围，未验证。

## 5. darwin-skill 优化（2026-10-02）

**基线（Phase 1，triage 用）**：结构 63.6 + 实测 20.7 = **84.3**。dim8=9 —— 2 prompt × 2 臂实跑：
带本 skill → 输出 χ²=1.267（datcmp 口径）、点出 CorMap 拒绝、并警告 DENSS 重采样网格插值会虚高；
不带 skill（baseline 臂）→ **漏掉 `--constant`**，得到 χ²=4.27，据此把结论写成「高分辨模型不比 DAMMIF 珠模更贴、差异在噪声内」（方向性错误）。
加权短板：dim3 3.6 / dim4 3.0 / dim8 2.3。

**三轮 paired 优化（每位 judge 在一次 call 内读改前+改后，3 位多数决）**：

| 轮 | 维度 | 改了什么 | paired 结果 | 处置 | commit |
|---|---|---|---|---|---|
| 1 | dim3 失败模式编码 | 新增 F 节 if-then 三段式故障表 14 条（触发条件/一线修复/仍失败兜底；全部来自实跑：datcmp 哨兵、只有水的链、重采样虚高、裁 q 只降 6.6 等） | 3-0 better（3 clear） | keep | `01293fb` |
| 2 | dim4 检查点设计 | 新增六道 🔴 STOP 闸门表（区间未定/分子状态未定/χ² 几十以上/区间内但 CorMap 拒绝/目录被并发重排/交付前），并把 E 节「判停点」升级为显性 🔴 STOP 标记 | 3-0 better（2 clear/1 slight） | keep | `c66edff` |
| 3 | dim8 产物与落地 | Step 4 换成可直接粘的绝对路径命令 + 实测控制台输出样例；Step 9 加「给用户的三段式回答」；生成器顶部告警只统计**参与排序**的曲线 | 3-0 better（1 clear/1 slight/1 未注） | keep | `756f6be` |
| — | 复审修正 | 第 3 轮 slight 指出「示例 `--qmin` 与样例输出不一致」→ 注明白未裁/裁后只影响 `[data]` 行 | 不计分（复审缺陷） | 修 | `fd2ed6d` |

**结果**：84.3 → **92.7**（绝对分只作 triage、且本轮为 dry_run 复评；keep 依据是三轮 paired 3-0）。
维度变化：dim3 7→9、dim4 5→9、dim2 9→10、dim5 9→10、dim9 8→9；体积 25.9 KB → 33.0 KB（仍在 150% 上限内）。
**收手理由**：MAX_ROUNDS=3 用尽，且第 3 轮已出现 margin=slight 的声音（F 表与停点表有部分重复）——按 HL-4 见好就收，不再硬凑。
**已验证**：优化后重新实跑生成产物（SEC 样例），README 的告警行与「判据 χ²（口径）」列按新规则输出；仓库与 profile 副本 `diff -rq` 逐字节一致。
