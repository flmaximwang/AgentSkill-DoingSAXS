# AgentSkill-DoingSAXS

"做 SAXS 这件事本身"的 skill 包：**拟合之前先判数据**，出现怪现象时怎么把责任判给数据、几何还是模型，
**把一条 SEC-SAXS 系列端到端跑完**（图像 → 报告）并交付一份人能看懂的结果说明，
以及**把高分辨模型嵌进从头算重建**（珠模型 + DENSS 电子云）并报出可解释的分数。
（单点判据 skill 在 `AgentSkill-UsingBioXTASRAW`；这里放**评估/归因** + **端到端流水线** + **嵌入/对齐**。）

| 来源 | 内容 | 沉淀成 |
|---|---|---|
| 2026-10-01 管式批次（BL19U2，11 个样品：A5-05-1…6 稀释序列 + 5705/877-apo/877-4zinc/97df/BSA）的实跑 | 逐帧质量（漂移/离群/对比度/SNR/误差诚实度）、低 q 上翘的三道归因检验、beamstop 几何与 q_min 依据 | [`skills/assess-saxs-raw-data-quality/`](skills/assess-saxs-raw-data-quality/) |
| 同批 A5-05-1…6 的 `embed/` 实跑（GNOM→DAMMIF×15→DAMAVER + CIFSUP；ChimeraX fitmap 的 R 阶梯）+ ATSAS CIFSUP/SUPCOMB 与 ChimeraX fitmap 官方文档 | 珠模型嵌入（CIFSUP/NSD，含分数只写在输出文件里、selection 口径、系综波动）与电子云嵌入（fitmap 的 resolution/search/seed 与姿态散布判据），以及 PyMOL 出图（ChimeraX `--nogui` 存不了图） | [`skills/embed-a-model-in-bead-and-density-models/`](skills/embed-a-model-in-bead-and-density-models/) |

**3 个 skill**：1 个评估/归因原子技能（24 个候选判据 → 通过 1 个）+ 1 条端到端 SEC 流水线（工程产物）+ 1 个"模型嵌入重建"的可执行流程。评估/归因那条的流程：：真实数据上磨出来的一条流程——
**先量几何（中心/掩膜/q_min）→ 逐帧 QC → 扣减后形状 → 低 q 上翘归因（空白-空白 / 背景形状失配 / 2D 差分）→ 完整结论（排除了什么、还剩什么、下一步做什么实验）**。
它**不做任何拟合**：拟合归 `AgentSkill-UsingBioXTASRAW` 的 15 个 skill，这里只判"那些拟合肥不肥"。

> **状态**：已 push 到 <https://github.com/flmaximwang/AgentSkill-DoingSAXS>（**public**，远程用 SSH，default=main）；
> 已安装进 default profile 的 **`saxs`** 类目（三段式标识符 + skills.sh/community，scan verdict **SAFE**）：
> `assess-saxs-raw-data-quality`（pin `eff1219`）、`run-a-sec-saxs-pipeline-end-to-end`（pin `fc148f6`）。
> `embed-a-model-in-bead-and-density-models` 同期并入同一 `saxs` 类目（pin `f6c48ad`，2026-10-01）。
> 注意它的 scan verdict 是 **CAUTION**：`scripts/embed-model.py` 命中 skills-guard 的 `python_os_environ` 与
> `python_subprocess` 两条规则（2×HIGH exfiltration on `os.environ`、3×MEDIUM execution on `subprocess`），
> 因此按既有先例用 `--force` 安装。命中项的实质是**本机可执行文件调用**（ATSAS `cifsup`、ChimeraX、PyMOL），
> 无网络、无外发；安装后 `diff -rq` 与仓库逐字节一致。
> 后者并入时另加了一项交付物：**结果目录里自动写 `README.md`**（只读产物、不重算；见该 skill 的 Step 4）。

## 索引

| skill | 用途 | 可执行入口 |
|---|---|---|
| [run-a-sec-saxs-pipeline-end-to-end](skills/run-a-sec-saxs-pipeline-end-to-end/SKILL.md) | **端到端跑一条 SEC-SAXS 系列**（全程 RAW，不自写拟合）：逐帧归一化（补 BL19U2 header txt / 读线站已有 txt）→ 归一化裁剪视频 → buffer/sample 区与扣减 → 多区间 Guinier → IFT → MW → 形状重建（DENSS 电子云 + ATSAS DAMMIF 珠模 + DAMAVER）→ RAW PDF 报告 → **结果目录自动写 `README.md`**（目录导航 + 关键数字 + 判读红线 + 本次告警），每个节点都落 `.dat`/表格 | `references/bl19u2-header-normalization.md`、`scripts/`（5 个脚本：emit-bl19u2-header-txt / crop-video-normalized / run-raw-sec-pipeline / results_readme / write-results-readme） |
| [assess-saxs-raw-data-quality](skills/assess-saxs-raw-data-quality/SKILL.md) | 一批 SAXS 原始帧的质量评估 + 低 q 上翘归因：先量几何与掩膜（含 RAW 读 Pilatus 的 **y 翻转**、beamstop 边缘 → `--qmin`），再逐帧判据表（对比度/漂移/离群/误差诚实度/SNR-qmax），再做三道归因（**空白-空白可复现极限**、**shape×电平**、**2D 差分：各向同性光晕 vs 紧贴 beamstop 的窄亮环**），最后按"浓度标度律"写出合格结论；含"单条曲线上翘 ≠ 相互作用"的判据与下一步实验设计 | `references/bl19u2-geometry-and-mask.md`、`references/upturn-attribution-protocol.md`、`scripts/`（5 个可执行脚本 + 1 个共用件） |
| [embed-a-model-in-bead-and-density-models](skills/embed-a-model-in-bead-and-density-models/SKILL.md) | **把高分辨模型（PDB）嵌进珠模型与 DENSS 电子云**：珠模型分支 = GNOM→DAMMIF×N→DAMAVER 取 `damfilt` 共识模型 + CIFSUP（`--method=NSD`，REGRID＝SUPCOMB fast 口径）给 **NSD**（分数只在输出文件里，附 selection 对照与系综波动）；电子云分支 = ChimeraX fitmap（**必须 `resolution` 才有 correlation、必须 `search+seed` 才是全局解**）跑 R 阶梯并报**姿态散布**（实测 20 Å 下存在相距 28.9 Å 的等分解）；另含 CA 几何复核（CA→最近 bead 距离）、两张分坐标系图（PyMOL，ChimeraX `--nogui` 存不了图）、每个样品一份 `embed/README.md` | `references/bead-and-density-embedding-parameters.md`、`scripts/embed-model.py`、`scripts/render-embed-figure.py`、`scripts/write-embed-readme.py` |

## 质量凭据（盲测）

按 Hermes 路由时的**真实 57 字符截断**（`agent/skill_utils.py:761`）把候选 skill 的 description 交给**独立评测者**逐条判路由
（候选 17 个 · 题目 10 条 = 6 正面 + 4 诱饵 · 交错排列 · 评测者看不到 SKILL.md 正文）：

| 轮次 | 对象 | 评测者 A | 评测者 B | 处置 |
|---|---|---|---|---|
| 1（首轮） | `embed-a-model-in-bead-and-density-models`，17 个候选，10 条 | **9/10**（正面 5/6 · 诱饵 **4/4**） | **9/10**（正面 5/6 · 诱饵 **4/4**） | 两评测者 10 题**逐题完全一致**；唯一错误 #6（ChimeraX `--nogui` 存图报 OpenGL）两票都判 *none*——归因为**可见窗口放不下**（57 字符已给主触发词），按既有口径记代价、不雕题面 |
| 1（首轮） | 本 skill，17 个候选，10 条 | **7/10**（正面 4/6 · 诱饵 **4/4**） | **7/10**（正面 4/6 · 诱饵 **4/4**） | 4 条诱饵 100% 未误收（流水线/Guinier 判据/P(r)/MW 各归各位）；三条一致错误里 #10 认定为**金标偏严**（修订后 8/10 · 8/10），#5（砍 q_min 砍到哪）与 #9（上机前排 control）属**可见窗口放不下** → 按既有口径记代价、不再调参 |

细节（含每条错误归因与金标修订理由）见 [`test-results.md`](test-results.md)；
本 skill 的首轮盲测见 [`skills/embed-a-model-in-bead-and-density-models/test-results.md`](skills/embed-a-model-in-bead-and-density-models/test-results.md)。

## 安装（三段式标识符，按仓库内路径，不需要 tap；`--category` 只决定落点）

```bash
hermes skills install <owner>/AgentSkill-DoingSAXS/skills/assess-saxs-raw-data-quality --category saxs -y
hermes skills install <owner>/AgentSkill-DoingSAXS/skills/run-a-sec-saxs-pipeline-end-to-end --category saxs -y
hermes skills install <owner>/AgentSkill-DoingSAXS/skills/embed-a-model-in-bead-and-density-models --category saxs -y
```

更新：`hermes skills update <skill 名>`（改了本仓库并 push main 之后）。

## 运行前提

- 脚本用 **RAW 自带的解释器**跑：`/Applications/BioXTASRAW/bin/python`（含 pyFAI/numba/matplotlib）。
- 不要在 RAW 源码目录里跑（`sascalc_exts` 会被遮蔽）。
- `embed-a-model-in-bead-and-density-models` 另需：**ATSAS ≥ 4**（`cifsup` = SUPCOMB 后继；GNOM/DAMMIF/DAMAVER 也来自它，本机 `/Applications/ATSAS-4.1.4-1`）、**ChimeraX**（fitmap，本机 1.9）、**PyMOL**（出图，`/Applications/PyMOL.app/Contents/bin/python3.10`）。
- 输入是"每样品一个目录、目录里同时有样品帧与夹着它的 control 帧"的布局（BL19U2 批处理风格）；
  文件命名形如 `A5-05-1_0013_00001.tif`（run 键 = 前两段）。

## 与其它包的边界

| 问题 | 去哪 |
|---|---|
| 这批原始帧好不好 / 上翘是真还是假 / q_min 取多少 | **本包** |
| **SEC 连续洗脱帧**图像→曲线→Rg/P(r)/MW/3D 端到端跑一遍 | **本包**（`run-a-sec-saxs-pipeline-end-to-end`） |
| **把高分辨模型嵌进珠模型 / DENSS 电子云并报分**（NSD / fitmap correlation） | **本包**（`embed-a-model-in-bead-and-density-models`） |
| 管式/静态帧端到端跑一遍 | `AgentSkill-UsingBioXTASRAW` 的 `run-a-tube-saxs-pipeline-end-to-end` |
| Guinier 怎么取点、P(r) 用哪个程序、MW 用哪一法、重建怎么评、模型配不配数据 | `AgentSkill-UsingBioXTASRAW` 的 13 个判据 skill |
| ATSAS 命令行（GNOM/DAMMIF/DATMW…） | `AgentSkill-UsingATSAS` |
