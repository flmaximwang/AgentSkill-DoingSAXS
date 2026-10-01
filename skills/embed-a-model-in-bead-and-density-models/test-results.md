# 路由盲测结果 —— embed-a-model-in-bead-and-density-models

方法（与包内其它 skill 同口径）：按 Hermes 路由的**真实 57 字符截断**（`agent/skill_utils.py:761`
里 description 只渲染 `desc[:57] + '...'`）把候选 description 交给**两个独立评测者**逐条判路由。
候选 **17 个**（本包 3 个 + `AgentSkill-UsingBioXTASRAW` 的 14 个 saxs 家族 skill，即最容易互相抢词的那一圈），
题目 **10 条** = 6 正面 + 4 诱饵，交错排列，评测者**看不到任何 SKILL.md 正文**、只能看截断后的 description。
评测者由 `delegate_task` 起两个隔离上下文的子代理（互不可见），各自独立作答。

## 结果

| 轮次 | 对象 | 评测者 A | 评测者 B | 一致？ |
|---|---|---|---|---|
| 1（首轮） | 本 skill，17 候选，10 题 | **9/10**（正面 5/6 · 诱饵 4/4） | **9/10**（正面 5/6 · 诱饵 4/4） | 10/10 逐题完全一致 |

逐题（A/B 的 choice 完全相同）：

| # | 类型 | 题 | 判到 | 金标 | 对否 |
|---|---|---|---|---|---|
| 1 | 正面 | 把 A5.pdb 嵌到珠模型和 DENSS 电子云里怎么叠 | embed-a-model-in-bead-and-density-models | 同 | ✅ |
| 2 | 正面 | CIFSUP 不打印，NSD 从哪读、0.87 算好吗 | embed-… | 同 | ✅ |
| 3 | 正面 | fitmap 只有 overlap，correlation 是 None | embed-… | 同 | ✅ |
| 4 | 正面 | 两个分数差不多、摆法差很远的解报哪个 | embed-… | 同 | ✅ |
| 5 | 正面 | 只有 denss.mrc 没有珠模型，要摆模型进去 | embed-… | 同 | ✅ |
| 6 | 正面 | ChimeraX 无头模式存图报 OpenGL 不可用 | **none** | embed-… | ❌（见下） |
| 7 | 诱饵 | 15 个 DAMMIF 模型、3 个簇、平均 NSD 0.89 能用吗 | evaluate-a-shape-reconstruction | 同 | ✅ |
| 8 | 诱饵 | CRYSOL 理论曲线低 q 差很多，模型错还是数据错 | fit-a-high-resolution-model-to-data | 同 | ✅ |
| 9 | 诱饵 | P(r) 的 Dmax 取多少、GNOM 还是 BIFT | compute-and-validate-p-of-r | 同 | ✅ |
| 10 | 诱饵 | 原始帧低 q 上翘是真散射还是 beamstop 边缘 | assess-saxs-raw-data-quality | 同 | ✅ |

**4 条诱饵 100% 未被误收**（重建评估 / 直接拟合 / P(r) / 原始帧 QC 各归各位）——这是本 skill 最要紧的一项，
因为"嵌入"和"拟合"最容易被混为一谈。两个评测者对 10 道题给出**完全相同**的答案，说明这个分界不是运气。

## 唯一那条错误的归因（#6）

- 现象：两个评测者都把"ChimeraX `--nogui` 存图报 OpenGL rendering is not available"判成 **none**，而不是本 skill。
- 归因：**可见窗口放不下**。本 skill 的触发词名额（前 57 字符）已经给了"把 model 嵌进珠模型
  （DAMMIF/DAMAVER）与 DENSS 电子云时用：CIFSUP 叠珠模型给 N…"，渲染工具（ChimeraX/PyMOL）在 57 字符之外；
  而这条题问的是**渲染环境**这个形态，脱离上下文确实可以判成"通用工具问题"。
- 处置：**不改 description**（改了就挤掉主触发词），按包内既有口径记代价：这条属于**会话内追问**
  （用户已经在本 skill 的流程里，Agent 手上有 SKILL.md 全文，B 段明确写了 ChimeraX 存不了图 → 用 PyMOL），
  不属于"冷启动要靠 description 独立命中"的场景。同一取舍在 `assess-saxs-raw-data-quality` 的首轮里
  也出现过（#5 砍 q_min、#9 上机前排 control），当时的处置也是记代价、不再调参。
- 若要进一步压缩这个代价，可把题面改成"这批嵌入图我怎么渲出来"（带上"嵌入图"这个本 skill 的词），
  但按既有约定不为了盲测分数去雕题面；如实记在案。

## 结论

- 正面命中 5/6、诱饵 0 误收、两评测者完全一致 → 本 skill 的 description 在**真实 57 字符窗口**下可用。
- 已知代价：渲染/出图类追问在冷启动下会被判成"无对应 skill"（会话内不受影响）。

---

# darwin 实测（dim8，2026-10-01）

方法：3 道题 × 2 臂 = 6 个隔离上下文的子代理。**带技能臂**读本 SKILL.md 全文后作答，
**不带技能臂**禁止读任何 skill 文件、只靠自身知识作答；两臂都禁止写文件、禁跑流水线。
按「输出是否完成用户意图 / 相对 baseline 的提升 / 有无 skill 引入的负面效果」打分。

| 题 | 带技能 | 不带技能 | 判定 |
|---|---|---|---|
| A：ChimeraX `--nogui` 存图失败 + 交付必须有 `.pse`（黄色 cartoon / 半透明白包络） | **9** | 5 | clear：baseline 建议灰底、把两个坐标系存进**一个** pse、给 `sphere_transparency 0.4` 且完全没提 `transparency_mode`（照做珠子仍不透） |
| B：同一张电子云两个分数相近、摆法差很远的解报哪个 | **9** | 6.5 | clear：baseline 给的是通用方法学（bootstrap/手性/聚类），没有「报法」与姿态散布口径，也没指向 CRYSOL 那条出口 |
| C：CRYSOL 低 q 差很多（诱饵题） | 7 | 7 | tie：两臂都答得像样，也都越到 `fit-a-high-resolution-model-to-data` 的地盘 |

带技能臂额外产出 3 条 skill 缺口（已评估、本轮未改，留档）：无 Linux/`xvfb-run`、ChimeraX `--offscreen` 的对照路径；
没有「两个摆法的分数差算不算并列」的阈值口径（脚本只在姿态散布 > 5 Å 时写 `ambiguity`）；
珠模型分支（CIFSUP NSD）没有与电子云分支对称的多解诊断。
