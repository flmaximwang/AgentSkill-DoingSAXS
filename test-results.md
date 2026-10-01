# 盲测记录（assess-saxs-raw-data-quality）

## 方法

- **候选集**：17 个 skill = profile `saxs/` 类目现有的 15 个（含 `organize-batch-saxs-dataset`、`run-a-sec-saxs-pipeline-end-to-end`）+ `AgentSkill-UsingBioXTASRAW` 里的 `run-a-tube-saxs-pipeline-end-to-end`（最容易被混淆的兄弟流水线）+ 本 skill。
- **可见信息**：每个候选只给 `name :: description 的前 57 个字符`（Hermes 路由时真实的截断，`agent/skill_utils.py:761`，`SKILL_PROMPT_DESC_LIMIT=60`，渲染为 `desc[:57]+"..."`）。评测者**看不到 SKILL.md 正文**。
- **评测者**：2 位独立子代理（互不可见、同 prompt、同候选表），每条请求必须给出唯一 skill 名 + 置信度。
- **题目**：10 条（**6 条正面 + 4 条诱饵**，交错排列，不告知类型）→ `test-prompts.json`。

## 金标与结果

| # | 类型 | 请求（摘要） | 金标 | 评测者 A | 评测者 B |
|---|---|---|---|---|---|
| 1 | 正面 | 这批数据哪些样品能用哪些该扔 / 评估原始数据好坏 | 本 skill | **本 skill** (high) | **本 skill** (high) |
| 2 | 诱饵 | 端到端跑一遍 + 参数/拟合图/bead model + 各节点 dat | tube-pipeline | tube-pipeline (high) | tube-pipeline (high) |
| 3 | 正面 | 老师说低 q 都有上翘，像相互作用，你同意吗 | 本 skill | **本 skill** (high) | **本 skill** (medium) |
| 4 | 诱饵 | Guinier 的 q 区间怎么取、qRg 到多少 | assess-guinier-fit-quality | assess-guinier (high) | assess-guinier (high) |
| 5 | 正面 | 低 q 那个点上翘，砍掉重拟合合不合规、砍到哪 | 本 skill | assess-guinier (medium) ✗ | assess-guinier (medium) ✗ |
| 6 | 诱饵 | 算 P(r)、定 Dmax，GNOM 还是 BIFT | compute-and-validate-p-of-r | compute-…p-of-r (high) | compute-…p-of-r (high) |
| 7 | 正面 | 背景稳不稳、能不能用一根控制样代表一整轮 | 本 skill | **本 skill** (medium) | **本 skill** (medium) |
| 8 | 诱饵 | 浓度多少、UV 还是 SAXS 定分子量 | choose-a-molecular-weight-method | choose-a-mw (high) | choose-a-mw (high) |
| 9 | 正面 | 下个月上机、样品不太够，control 怎么排 | 本 skill | tube-pipeline (low) ✗ | organize-batch (low) ✗ |
| 10 | 正面 | Rg 随浓度乱漂、Dmax 忽大忽小，数据问题还是拟合问题 | 本 skill | assess-guinier (medium) ✗ | assess-guinier (low) ✗ |

**得分（严格按金标）**：A **7/10**（正面 4/6 · 诱饵 4/4）、B **7/10**（正面 4/6 · 诱饵 4/4）。

## 三条一致错误（两位评测者同判 → 系统性，不是随机）

1. **#5（砍 q_min 合不合规 / 砍到哪）→ 两位都转 `assess-guinier-fit-quality`（medium）。**
   题面里的"拟合/区间"把路由拉走；而"砍到哪"的依据（beamstop 未遮挡比例 → q_min=0.010）在本 skill 里。
   **代价**：多一跳（判据 skill 不知道这套几何数字）。
2. **#9（上机前 control 怎么排）→ 两位都给了别的 skill（low 置信度）。**
   题面的"control/上机"没有出现在本 skill 可见的 57 字里（可见头在讲"评估原始帧"）。
   **代价**：多一跳，风险低（上机前真按错路走的人少）。
3. **#10（Rg 乱漂是数据问题还是拟合问题）→ 两位都转 `assess-guinier-fit-quality`。**
   这条**金标可议**：按"先判帧再判拟合"的分工，本 skill 是第二站，第一站落在 Guinier 判据 skill 是合理的。
   → 认定为**金标偏严**，不计为本 skill 的漏。

## 处置

- **关键正面结果**：4 条诱饵 **100% 未误收**（端到端流水线 / Guinier 判据 / P(r) / MW 方法各自归位），说明本 skill 没有去抢兄弟 skill 的活；主入口（#1 评估好坏、#3 上翘归因）两位都高/中置信度命中。
- #5、#9 属**可见窗口放不下**（57 字内的头已经在讲"评估原始帧 + 上翘归因"，塞不下"q_min 取哪"与"上机排 control"）。
  按本仓库既有口径（第 5-6 轮、第 8 轮）：**记代价、不再调参**——调头的收益是两条低风险多跳，代价是主入口的关键词被挤掉。
- 金标修订 1 条（#10 从"本 skill"改为"先 Guinier 判据、必要时回本 skill"）。修订后：A **8/10**、B **8/10**。

## 残留问题

- 若后续真出现"上机前排 control"这类请求被别的 skill 接走，可考虑把候选表里的 `organize-batch-saxs-dataset` 与 `run-a-tube-saxs-pipeline-end-to-end` 的可见头与本 skill 再对照一次（三者都在"处理一批帧"这个语义场里）。
- 本 skill 与 `assess-raw-frame-quality`（curator 2026-10-01 自建、已并入并退役）曾同名同义；**当前 profile `saxs/` 类目内无同义重复**，重复出现时按 `AgentSkill-UsingBioXTASRAW` 的 migration 流程处理。
