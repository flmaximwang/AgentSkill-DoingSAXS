# 路由盲测（首轮）：`write-saxs-results-readme`

做法与包内其它技能一致：按 Hermes 路由时的**真实 57 字符截断**（`agent/skill_utils.py:761`
渲染 `desc[:57] + "…"`）把候选技能交给**两个独立评测者**逐条判路由。
评测者只看到「技能名 | 截断后的 description」，看不到 SKILL.md 正文，也看不到金标；
题目 = **5 条正面**（应走本技能）+ **4 条诱饵**（同一批数据/同一批产物，但应走别的技能）。

## 候选（13 个，含本技能；截断形态）

| 候选 | 可见的 description 前 57 字 |
|---|---|
| 1 run-a-sec-saxs-pipeline-end-to-end | 端到端跑一条 SEC-SAXS 系列（图像→报告，全在 RAW 里做）：补出每帧 BL19U2 header tx… |
| 2 run-a-tube-saxs-pipeline-end-to-end | 管式/静态 SAXS 一批帧端到端出全套结果（目录内样品+control→扣减→多区间 Guinier→BIFT/… |
| 3 assess-saxs-raw-data-quality | 评估一批 SAXS 原始帧质量并归因低 q 上翘（真散射 / beamstop 边缘 / 背景残留）：先量几何与掩… |
| 4 assess-guinier-fit-quality | 从一条扣减过的 SAXS 曲线读 Rg/I0，并判断这次 Guinier 拟合信不信得过：n_min 低 q 取点… |
| 5 compute-and-validate-p-of-r | 要算 P(r)/距离分布、定 Dmax、选 GNOM/DIFT/BIFT 时用：三法按下游重建程序选，Dmax 八… |
| 6 choose-a-molecular-weight-method | 分子量六法怎么选、能不能信（SEC 用它判低聚态/单体二聚体）：两轴定位（RAW 原生 4 法 vs ATSAS … |
| 7 evaluate-a-shape-reconstruction | DAMMIF/DENSS 形状重建做完后评估能不能用：a-score/平均 NSD/被剔模型数/聚类数/各模型 χ… |
| 8 result-folder-readme | 批处理跑完要给结果目录写人能看懂的 README 时用：分点总结+参数明细表+🟢🟡🔴。 |
| 9 embed-a-model-in-bead-and-density-models | 把 model 嵌进珠模型（DAMMIF/DAMAVER）与 DENSS 电子云时用：CIFSUP 叠珠模型给 N… |
| 10 process-sec-saxs-series | 处理 SEC-SAXS 系列（连续洗脱帧）；SEC 浓度未知→禁用 I0 标样/绝对刻度；色谱图→buffer/s… |
| 11 reduce-saxs-frames-to-curves | 把 SAXS 帧还原成 1D 曲线（积分→平均→扣减→存 .dat）；SEC 连续洗脱帧不归它（走 process… |
| 12 correct-sec-saxs-baseline | 扣减后强度-帧号仍漂移时按性质选基线校正：束流/仪器漂移→Linear，毛细管污垢→Integral；含过校正识别… |
| **13 write-saxs-results-readme（本技能）** | **SAXS 结果目录（SEC / 管式）验收结构 + 生成给人看的 README.md：两种流水线共用同一套 8 节…** |

## 题目与两位评测者的判定（逐条一致）

| # | 题面（用户原话摘要） | 金标 | A 判 | B 判 |
|---|---|---|---|---|
| 1 | 结果目录要交付，看结构齐不齐、README 对不对 | **13 本技能** | 13 ✓ | 13 ✓ |
| 2 | 两条流水线的 README 格式不一样，能统一吗 | **13 本技能** | 13 ✓ | 13 ✓ |
| 3 | 半年前的结果没有 README，补一份 | **13 本技能** | 13 ✓ | 13 ✓ |
| 4 | README 里的 Rg 怎么确认不是手抄错的 | **13 本技能** | 13 ✓ | 13 ✓ |
| 5 | 文件夹里有哪些文件、先看哪个、哪些是最终产物 | **13 本技能** | 13 ✓ | 13 ✓ |
| 6 | 扣减曲线 Rg 取哪一段、q_max 到多少 | 4 `assess-guinier-fit-quality` | 4 ✓ | 4 ✓ |
| 7 | 这批管式 tif 端到端跑一遍，各节点 dat 都要 | 2 `run-a-tube-…-pipeline-end-to-end` | 2 ✓ | 2 ✓ |
| 8 | 珠模 NSD 0.8、2 个 cluster，能用吗 | 7 `evaluate-a-shape-reconstruction` | 7 ✓ | 7 ✓ |
| 9 | 批处理 22 个样品，根目录写一份总览 README | 8 `result-folder-readme` | 8 ✓ | 8 ✓ |

| 评测者 | 正面（应走本技能） | 诱饵（不该走本技能） | 与金标不一致的题 |
|---|---|---|---|
| A | **5/5** | **4/4** | 无 |
| B | **5/5** | **4/4** | 无 |

## 归因

- **诱饵 4/4 全未误收**，包括最像的两条：#9（也是「写 README」——但那是**批处理根目录**的总览，
  归 `result-folder-readme`）与 #8（README 里恰好也写 NSD，但「能不能用」的判据归
  `evaluate-a-shape-reconstruction`）。说明 description 里的**排除条款**（「不负责数字本身的判据…也不负责
  跑流水线」）在 57 字符之外仍然通过主触发词把边界划住了。
- **#7 未被误收**：管式端到端那条的可见头以「管式/静态 SAXS 一批帧端到端出全套结果」开头，
  与「结果目录怎么写」不抢词。
- 与包内前两轮不同，本轮**两评测者逐题完全一致、零错误**，因此没有金标修订与调参空间；
  按既有口径记录即止，不再雕题面。

## 复现

```bash
# 候选截断文本与题目：skills/write-saxs-results-readme/test-prompts.json
# 用 delegate_task 发 2 个独立评测者（只看截断 description，不看 SKILL.md）
```
