# AgentSkill-DoingSAXS

"做 SAXS 这件事本身"的 skill 包：**拟合之前先判数据**，以及出现怪现象时怎么把责任判给数据、几何还是模型。
（下游的判据/流水线 skill 在 `AgentSkill-UsingBioXTASRAW`；这里放的是**评估与归因**。）

| 来源 | 内容 | 沉淀成 |
|---|---|---|
| 2026-10-01 管式批次（BL19U2，11 个样品：A5-05-1…6 稀释序列 + 5705/877-apo/877-4zinc/97df/BSA）的实跑 | 逐帧质量（漂移/离群/对比度/SNR/误差诚实度）、低 q 上翘的三道归因检验、beamstop 几何与 q_min 依据 | [`skills/assess-saxs-raw-data-quality/`](skills/assess-saxs-raw-data-quality/) |

**1 个 skill**（24 个候选判据 → 通过 1 个原子技能）：真实数据上磨出来的一条流程——
**先量几何（中心/掩膜/q_min）→ 逐帧 QC → 扣减后形状 → 低 q 上翘归因（空白-空白 / 背景形状失配 / 2D 差分）→ 完整结论（排除了什么、还剩什么、下一步做什么实验）**。
它**不做任何拟合**：拟合归 `AgentSkill-UsingBioXTASRAW` 的 15 个 skill，这里只判"那些拟合肥不肥"。

> **状态**：已 push 到 <https://github.com/flmaximwang/AgentSkill-DoingSAXS>（**public**，远程用 SSH，default=main）；
> 已安装进 default profile 的 **`saxs`** 类目（三段式标识符 + skills.sh/community，scan verdict **SAFE**，pin `eff1219`）。

## 索引

| skill | 用途 | 可执行入口 |
|---|---|---|
| [assess-saxs-raw-data-quality](skills/assess-saxs-raw-data-quality/SKILL.md) | 一批 SAXS 原始帧的质量评估 + 低 q 上翘归因：先量几何与掩膜（含 RAW 读 Pilatus 的 **y 翻转**、beamstop 边缘 → `--qmin`），再逐帧判据表（对比度/漂移/离群/误差诚实度/SNR-qmax），再做三道归因（**空白-空白可复现极限**、**shape×电平**、**2D 差分：各向同性光晕 vs 紧贴 beamstop 的窄亮环**），最后按"浓度标度律"写出合格结论；含"单条曲线上翘 ≠ 相互作用"的判据与下一步实验设计 | `references/bl19u2-geometry-and-mask.md`、`references/upturn-attribution-protocol.md`、`scripts/`（5 个可执行脚本 + 1 个共用件） |

## 质量凭据（盲测）

按 Hermes 路由时的**真实 57 字符截断**（`agent/skill_utils.py:761`）把候选 skill 的 description 交给**独立评测者**逐条判路由
（候选 17 个 · 题目 10 条 = 6 正面 + 4 诱饵 · 交错排列 · 评测者看不到 SKILL.md 正文）：

| 轮次 | 对象 | 评测者 A | 评测者 B | 处置 |
|---|---|---|---|---|
| 1（首轮） | 本 skill，17 个候选，10 条 | **7/10**（正面 4/6 · 诱饵 **4/4**） | **7/10**（正面 4/6 · 诱饵 **4/4**） | 4 条诱饵 100% 未误收（流水线/Guinier 判据/P(r)/MW 各归各位）；三条一致错误里 #10 认定为**金标偏严**（修订后 8/10 · 8/10），#5（砍 q_min 砍到哪）与 #9（上机前排 control）属**可见窗口放不下** → 按既有口径记代价、不再调参 |

细节（含每条错误归因与金标修订理由）见 [`test-results.md`](test-results.md)。

## 安装（三段式标识符，按仓库内路径，不需要 tap；`--category` 只决定落点）

```bash
hermes skills install <owner>/AgentSkill-DoingSAXS/skills/assess-saxs-raw-data-quality --category saxs -y
```

更新：`hermes skills update assess-saxs-raw-data-quality`（改了本仓库并 push main 之后）。

## 运行前提

- 脚本用 **RAW 自带的解释器**跑：`/Applications/BioXTASRAW/bin/python`（含 pyFAI/numba/matplotlib）。
- 不要在 RAW 源码目录里跑（`sascalc_exts` 会被遮蔽）。
- 输入是"每样品一个目录、目录里同时有样品帧与夹着它的 control 帧"的布局（BL19U2 批处理风格）；
  文件命名形如 `A5-05-1_0013_00001.tif`（run 键 = 前两段）。

## 与其它包的边界

| 问题 | 去哪 |
|---|---|
| 这批原始帧好不好 / 上翘是真还是假 / q_min 取多少 | **本包** |
| 图像→曲线→Rg/P(r)/MW/3D 端到端跑一遍 | `AgentSkill-UsingBioXTASRAW` 的两条 pipeline skill |
| Guinier 怎么取点、P(r) 用哪个程序、MW 用哪一法、重建怎么评 | `AgentSkill-UsingBioXTASRAW` 的 13 个判据 skill |
| ATSAS 命令行（GNOM/DAMMIF/DATMW…） | `AgentSkill-UsingATSAS` |
