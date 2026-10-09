---
name: write-saxs-results-readme
description: "SAXS 结果目录（SEC / 管式）验收结构 + 生成给人看的 README.md：两种流水线共用同一套 8 节模板与「采用区间/采用那一支 P(r)/三档可用性」口径，README 里每个数字都能追回产物原文。用于「这个结果目录缺什么、能不能交付」「README 看不懂 / 两份 README 格式不一样」「跑完要一份说明」「旧目录补写 README」；不负责数字本身的判据（Rg 取点→assess-guinier-fit-quality，P(r)/Dmax→compute-and-validate-p-of-r，MW→choose-a-molecular-weight-method，重建→evaluate-a-shape-reconstruction），也不负责跑流水线（SEC→run-a-sec-saxs-pipeline-end-to-end，管式→run-a-tube-saxs-pipeline-end-to-end）。再含 README 的写作口味（分点总结在前、🟢🟡🔴 三档、明细表带「误差怎么读」、生成器补块间空行、顶部警惕计数只算参与判决的行），折自通用 result-folder-readme，本包一律以 8 节为准。"
version: 1.0.0
author: hermes
license: MIT
tags: [saxs, readme, deliverable, results-folder, verification, sec-saxs, tube-saxs, bioxtas-raw]
metadata:
  hermes:
    tags: [saxs, readme, deliverable, results-folder, verification, sec-saxs, tube-saxs, bioxtas-raw]
    related_skills:
      - slug: run-a-sec-saxs-pipeline-end-to-end
        relation: composes-with
      - slug: run-a-tube-saxs-pipeline-end-to-end
        relation: composes-with
      - slug: assess-saxs-raw-data-quality
        relation: composes-with
---

# SAXS 结果目录：验收结构 + 写一份人看得懂的 README

**这条 skill 是「结果目录怎么交付」的唯一文本来源**：SEC 与管式两条端到端流水线**都不再各自写 README**，
它们产出的目录由本 skill 的同一份实现（`scripts/readme_common.py`）套同一个 8 节模板生成 README，
并由 `verify-results-folder.py` 验收「结构齐不齐 + README 里的关键数字与产物对不对得上」。

一句话定位：**它回答"这个结果文件夹缺什么、合不合规范、README 里的数是不是真来自产物"，不回答
"这个 Rg 能不能信"（后者是判据类 skill 的事）。**

## When to Use（什么时候用）

- 一条流水线刚跑完 / 或跑完了一阵子，要一份**给人看的结果说明**（结论、关键数字、判据、怎么复核）。
- 要交付或归档一个结果目录，先问**结构齐不齐**（哪些 must 产物没落盘、哪些节点被跳过、为什么）。
- 用户说「我看不懂」「两份 README 不一样」「怎么知道 README 里的数字没写错」。
- 两条流水线的 README 要对齐时：改模板只改这里（**不要**在各流水线里就地改格式）。
- 已经跑完但没 README 的旧目录：补写（不动数据）。

**不负责**：

| 问题 | 去哪 |
|---|---|
| 这批原始帧好不好 / 低 q 上翘是真还是假 / q_min 取多少 | `assess-saxs-raw-data-quality` |
| Rg 取哪一段、P(r)/Dmax 用哪个引擎、MW 用哪一法、重建怎么评 | `assess-guinier-fit-quality` / `compute-and-validate-p-of-r` / `choose-a-molecular-weight-method` / `evaluate-a-shape-reconstruction` |
| SEC 系列端到端跑一遍 / 管式一批帧端到端跑一遍 | `run-a-sec-saxs-pipeline-end-to-end` / `run-a-tube-saxs-pipeline-end-to-end` |
| 高分辨模型嵌进重建 | `embed-a-model-in-bead-and-density-models` |

## 契约：一个结果目录长什么样

两种模式各有 must / should / info 三档（缺 must = 交付不完整，缺 should = 第 5 节必须写明原因）：

| 级别 | SEC（`run_meta.json` 为标记） | 管式（`summary.json` 为标记） |
|---|---|---|
| must | `run_meta.json`、`profiles/01_integrated`、`profiles/03_subtracted`、`profiles/04_sample`、`tables/{guinier_multi_range,ift_summary,mw}.csv` | `summary.json`、`profiles/01_control`、`profiles/02_sample`、`profiles/03_subtracted`、`profiles/04_guinier`、`ifts/`、`tables/{guinier_multi_range,ift_summary,mw}.csv` |
| should | `series/`、`profiles/02_buffer`、`profiles/06_guinier`、`reports/`、`ifts/` | `tables/shape_results.json`、`qc.png`、`reports/`、`models/` |
| info | `tables/frame_params.csv`、`tables/frames_integrated.csv`、`models/`、`video/`、`norm/` | `norm/frame_qc.csv`、`frames/` |

完整的逐项说明与「章节规范」见 `references/results-folder-contract.md`。

### 本机项目里的落点（BL19U2 / DataProcess_2026.10.01）

产物根不是随便一个目录，本机固定是
`/Users/maxim/Repositories/DataProcess_2026.10.01/processed/<模式>/<样品>/`（**该仓库不是 git 仓库** ——
别再往上一层 `git status`，也没有版本可回滚；这个结果目录的"可追溯"全靠本 skill 写的 README 与表）：

- `data/<模式>/<样品>/` 放**原始帧**（线站逐帧 `.txt` 与 `.tif` 并排），`data/*.cfg` 是当天的 BL19U2 配置；
- `processed/<模式>/<样品>/` 放**结果**：`profiles/01_integrated|02_buffer|03_subtracted|04_sample`、
  `tables/ift_summary.csv`、`models/`、`reports/`、`video/`，以及本 skill 写的 **`README.md`**；
- **embed 类产物在"同级"的 `processed/SEC-SAXS/_embed/<样品>/`**（不塞进样品目录里）；
- 本机解释器 `/Applications/BioXTASRAW/bin/python`，ATSAS 在 `/Applications/ATSAS-4.1.4-1/bin`。

「一个结果目录」= 上面那个 `processed/<模式>/<样品>/`；`_embed/`、`_summary/`、`_logs/` 都是**它的兄弟**目录，
不在本 skill 契约的 must/should 清单里。

## 8 节模板（两种模式**必须一样、顺序一样**）

| 节 | 内容 | 两种模式的差别 |
|---|---|---|
| `## 0. 结论速览` | 🟡 警示引用块 + 数据可用性 + 表现较好/较差的参数（**与第 1 节同一张表**）+ 「只看这几个文件就够」 | 只看文件那 3 条按模式给（管式先 `qc.png`，SEC 先报告 + `series_plot.png`） |
| `## 1. 最终拟合参数（明细表）` | `环节 / 参数 / 值 / 程序误差 / 能不能用 / 误差·不确定度怎么读`，Guinier→IFT→GNOM→MW→珠模→DAMAVER→DENSS 固定顺序 | 同一套**行**，有数据就出、没有就整行不出现 |
| `## 2. 关键结果` | `项目 / 数值 / 怎么看`：对照与缩放、归一化、对比度、分析窗、auto-Guinier、采用区间、跨区间一致性、P(r)、MW、3D | SEC 给 buffer/样品帧区间，管式给对照 run 名与缩放因子 |
| `## 3. 这个文件夹里有什么（按建议阅读顺序）` | 无序列表，每条 `路径 —— 是什么（什么时候看）`；只列真实存在的文件，末尾报「还有 N 个文件没列出」 | 条目按模式筛选（同一份字典，见 reference） |
| `## 4. 每个数字的判据` | `数字 / 好 / 勉强 / 差` 一张表 | **完全相同** |
| `## 5. 这次没做的 / 不能信的` | 跳过的节点 + 失败原因 + 告警；附 IFT 各起跑点/引擎全表（含三条闸门） | 告警按各自产物给（SEC 有「没有 models/」类；管式有「没 --save-frames」类） |
| `## 6. 本次用的参数（复现用）` | bash 代码块（管式直接引用 `summary.json` 记的原始命令行）+ RAW 设置/ATSAS 一行 | 命令行来源不同 |
| `## 7. 想自己复核` | workspace 怎么开、表在哪、怎么重跑、低 q 该找谁、判据找谁；`### 这些缩写`；`### 目录树（两层）` | 同一套 |

## 统一的选值口径（这是「两边对齐」的机制，不许各写各的）

1. **采用区间（Rg）**：① 产物自己声明了推荐区间**且闸门全过** → 用它；② 产物声明了但**有闸门没过**
   （它的 `reason` 里写着 NOT recommended）→ **不采用**，退回通用规则，并在 README 里写明这句话；
   ③ 通用规则 = 收敛（`rg > 0`）+ `qRg ∈ [0.28, 1.35]` + R² 最高；④ 都不满足 → **没有可引用的 Rg**，不报数。
2. **跨区间一致性**：只统计「收敛且 R² ≥ 0.9」的区间，其余**列出来说明为什么排除**
   （把 `rg = -1` 的失败区间算进一致性 = 一票否决，是错的）。
3. **采用那一支 P(r)**：先看产物声明的 `chosen` / `trusted`，再看三条闸门（Rg 对得上 / Dmax ≤ 4.5·Rg /
   Dmax 未越出搜索网格），再看 χ²；SEC 没有产物级 `trusted`，用同一套闸门现算。
4. **数据可用性三档**：🟢 = 采用区间 R² ≥ 0.99、区间间 Rg 漂 ≤ 5%、IFT 有可信的一支；
   🟡 = 有可引用的 Rg 但上面打折（漂移 / R² 略低 / 产物标了 NOT recommended）；
   🔴 = 没有可引用的 Rg，或采用区间 R² < 0.95。**三档的定义写在 `_availability()` 里，两种模式同一条代码路径。**
5. **字段别名**：两种流水线的表头不一样（`rg` vs `Rg`、`q_min` vs `qmin`、`mw` vs `MW_kDa`、
   `chi_sq` vs `chisq`、`detail1..4` vs `aux`…），全部由 `ALIAS` 表归一 —— 见 reference 的别名表。

## 执行步骤

### Step 1 — 写 README（只读产物，不重算）

```bash
python write-readme.py <结果目录> [<结果目录> ...]      # 写 <目录>/README.md
python write-readme.py --stdout <结果目录>            # 只打印（不落盘，适合先看一眼）
```

管线（`run-raw-sec-pipeline.py` / `run-raw-tube-pipeline.py`）在**最后一步自动调同一份实现**：
它们按「同类目下的同名技能目录」找 `write-readme.py`，找不到就只告警不报错（README 是交付物，不是计算步骤）。
两侧都有 `--readme-script <路径>` 可显式指定（脚本被拷走单独用时用得上）。

完成标准：`<结果目录>/README.md` 存在，且 8 节齐全、顺序与上面的表一致。

### Step 2 — 验收（交付前必跑）

```bash
python verify-results-folder.py <结果目录> [--all]     # 退出码 1 = 有 🔴
```

四道检查：

1. **结构**：按契约逐项点名（🔴 缺 must / 🟡 缺 should 或空目录 / 🟢 就位）。
2. **README 章节**：8 节齐不齐、顺序对不对。
3. **README ↔ 产物**：README 里该出现的每条关键数字是否都能在产物里找到原文
   （被人手改过、或产物重算过而 README 没重写，都会在这里红）。
4. **产物 ↔ 产物**：两条独立路径（json ↔ csv）读同一个数是否一致；CSV 是四舍五入后的文本，
   比较按**它的精度**（`round(a, nd) == round(b, nd)`），不当成矛盾。
   另有时间戳检查：关键表比 README 新 ⇒ README 旧了（🔴）。

完成标准：退出码 0（只有 🟡 可以带进交付，但要在报告里点名）。

两种模式各跑一遍，**关键结果必须逐条对得上**（同一套章节、同一套表格列、同一批闸门）——
这就是"确保两条流水线交付物对齐"的可执行判据。

## 面向人的写法约定（改模板时必须遵守）

- **memo 式的条目一律写成无序列表**（`- `）。不要把若干短句用 `·`、`；` 串成一大段"看着像段落"的碎片。
- **要成段的地方必须是有逻辑的完整句子**（README 顶部出处说明、每张表上面的「怎么读」引用块）。
- 表格**每一行都要带「怎么看 / 什么时候看」一列**：只给数字、不给判据的表不算交付物。
- 颜色只有三档：🟢 可用 / 🟡 需要警惕 / 🔴 不可用；**结论先行**，明细在后。
- 数字格式统一走 `_fmt` / `_fmt_err` / `_fmt_mw`：误差两位有效数字、`rg = -1` / `inf` / `nan`
  一律翻译成「未收敛（哨兵 -1）」「inf（发散）」「—」，**绝不把哨兵值当数字印出来**。
- 上面这些口味的完整推导、以及折自通用 `result-folder-readme` 的对照与 **6 节 vs 8 节**的取舍，见 `references/results-folder-contract.md` §5。

## 坑（都是实测撞出来的）

- **`summary.json` 里的布尔可能被 `default=str` 写成字符串 `"False"`**，而 `"False"` 在 Python 里是**真**值
  —— 读闸门必须显式归一化（`_gate()`），否则「没过闸门」会被读成「全过」。
- **管式 `guinier_recommended` 不为空 ≠ 可用**：它可能是**明确标了 NOT recommended** 的（`gates` 里有 false、
  `reason` 里写着NOT recommended）。要读 `gates`，不能只看有没有这个键（本机 A5-05-1 就是这样：
  Rg 在区间之间从 10.8 漂到 19.6 Å，产物自己写了 NOT recommended）。
- **管式旧版 `tables/mw.csv` 的 `aux` 是直接 `str(list)` 写进去的**（内含逗号 → 列被撑开、表头与数据错位），
  `csv.DictReader` 会读歪。本 skill 用按位置读的 `_read_mw_csv()`；新写入用 `;` 连接（见流水线的 `node_mw`）。
- **`inf` / `nan` 不是合法 JSON**：`_parse_aux` 走的是「逐 token 试 float」的回退路径
  （5705 的 `datmw_bayes` 就是 `inf`，Datclass 是 `'random-chain', -1.0` 这种混合类型）。
- **SEC 的 `run_meta.json` 把表存成位置列表**（`guinier_rows` / `ift_rows`），列序与同目录 CSV 相同
  —— 顺序变了要同步改 `_SEC_G_KEYS` / `_SEC_I_KEYS`。
- **SEC 的 `denss_rows[0][0]` 是整幅密度数组的字符串**（不是 chi²）：chi²/Rg/体积/盒子在 `[1..4]`。
- **README 必须**是管线最后写的那一个文件**：否则时间戳检查会红（写 README 之后又生成图/表 = 读者拿到的是旧数字）。
  管式的图要放在 README **之前**（图炸了也要写出 README），SEC 的 README 本来就在最后。
- **判据文本不要在两处各写一份**：`RULES`、`ABBREV`、`FILE_DOC` 都只有这一份；发现格式不一致时改这里，
  改完对两条流水线各跑一次验收。
- 空目录（存在但 0 文件）不算「就位」：非 must 的记 🟡（`frames/` 没加 `--save-frames` 就是这种）。

## 相关 skills

- `run-a-sec-saxs-pipeline-end-to-end` / `run-a-tube-saxs-pipeline-end-to-end` — 产出结果目录的两条流水线；
  它们调本 skill 的 `readme_common.write_readme()`，自己不写 README。
- `assess-saxs-raw-data-quality` — 低 q 能不能信该找它（README 第 7 节直接把读者导过去）。
- `assess-guinier-fit-quality` / `compute-and-validate-p-of-r` / `choose-a-molecular-weight-method` /
  `evaluate-a-shape-reconstruction` — 第 1/4 节那张表背后的判据细节。
- `result-folder-readme` — 批处理**根目录**汇总 README 的写法（每个样品一份 README 的上一层）。本 skill 已把它的净值（分点总结在前、生成器补块间空行、顶部警惕计数口径、一次性小任务边界）折进 `references/results-folder-contract.md` §5，且本包一律以 8 节为准。

## Skill Structure

<!-- Generated by Scripts -->

```
write-saxs-results-readme/
├── SKILL.md  (194 lines)
├── test-prompts.json  (69 lines)
├── test-results.md  (61 lines)
├── references/
│   ├── results-folder-contract.md  (200 lines)
│   └── verification-rules.md  (61 lines)
└── scripts/
    ├── readme_common.py  (1732 lines)
    ├── verify-results-folder.py  (95 lines)
    └── write-readme.py  (66 lines)
```

<!-- Generated by Scripts -->
