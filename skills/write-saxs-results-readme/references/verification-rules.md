# 验收规则（`verify-results-folder.py` 到底在验什么）

验收的目的只有一个：**交付前证明「README 里写的关键结果」确实来自这个目录里的产物，且两条流水线的
交付物逐条对得上**。所以四道检查都围绕"能不能追回原文"，而不是"数字好不好看"（那是判据类 skill 的事）。

## 四道检查

| # | 检查 | 判据 | 红色意味着 |
|---|---|---|---|
| 1 | 结构 | 按契约逐项：must 缺 = 🔴；should 缺 / 目录存在但空 = 🟡；🟢 = 就位 | 交付不完整（比如 `tables/ift_summary.csv` 没落盘） |
| 2 | README 章节 | 8 节齐全且顺序与 `SECTION_ORDER` 一致 | 模板被改过 / 手工删节 |
| 3 | README ↔ 产物 | `build(facts)` 算出的**每条关键数字**都要能在 README 文本里找到原文 | 有人手改了 README 的数字，或产物重算过而 README 没重写 |
| 4 | 产物 ↔ 产物 | 两条独立路径（`summary.json`/`run_meta.json` ↔ `tables/*.csv`）读同一个数必须一致 | 产物自相矛盾（其中一张表被重算/被改过） |

外加**时间戳检查**（第 3 组的最后一行）：管线最后一步才写 README，所以

- `tables/*.csv`、`tables/*.json`、`summary.json`、`run_meta.json` 比 README 新 ⇒ 🔴（README 旧了，先重写）；
- 其它产物（`models/*`、`qc.png`、`reports/*`、视频…）比 README 新 ⇒ 🟡（数字不受影响，但可重跑一次同步）。

## 容差：CSV 是四舍五入后的文本

`tables/*.csv` 里写的是**按位四舍五入**的数（`%.2f` / `%.4f` / `%.5f` 不等），`summary.json` 里是全精度。
所以第 4 道检查**按 CSV 自己的精度**比：`_same(a, b, b_text)` → `round(a, nd) == round(b, nd)`，
`nd` = CSV 文本的小数位数。例：`48.41485`（json）与 `48.41`（csv）**一致**，不是矛盾。
两个都带小数点的数按更"粗"的那一个比；没有小数位时退回 `1e-9` 相对容差。

同理，管式推荐区间的行**身份**用 `(qmin, qmax)` 这一对（都按 CSV 精度比）来匹配，
多个行落在同一对四舍五入值上时取 Rg 最接近的那个 —— 只按 `qmin` 匹配会匹错行
（实测 `idx 129-162` 与 `idx 129-177` 的 `qmin` 都写成 `0.05666`）。

## 档位与退出码

- 🔴 = 硬伤（缺 must 产物 / 章节缺失或乱序 / 关键数字追不回原文 / 产物互相矛盾 / README 比关键表旧）；
- 🟡 = 提示（缺 should 产物、空目录、非关键产物比 README 新）；
- 退出码：**0 = 没有 🔴**，1 = 有 🔴（交付前必须处理）。

## 负例实测（确保验收真的会红）

在本机 `processed/` 的真实结果目录副本上逐条打过：

| 负例 | 做法 | 验收结果 |
|---|---|---|
| README 里的数字被手改 | 把 `Rg 18.97 Å` 改成 `19.97 Å` | 🔴 `README ↔ 产物`：`Guinier 拟合 Rg = 18.97 Å` 这条找不到原文 |
| 缺 must 产物 | 删 `tables/ift_summary.csv` | 🔴 结构：`tables/ift_summary.csv 缺失（交付不完整）`，退出码 1 |
| 产物自相矛盾 / README 旧 | 重算/改关键表而没重写 README | 🔴 `README 与数据的时间戳`：`关键表比 README 新 N 分钟`，退出码 1 |
| README 不存在 | 删 `README.md` | 🔴 提示先跑 `write-readme.py` |
| 空目录 | `frames/` 存在但 0 文件 | 🟡（非 must） |

正例：`processed/{SEC-SAXS/bsa,Tube-SAXS/A5-05-1,Tube-SAXS/07PR0111-02-6}` 三份副本
（含一个 DENSS 失败 + 珠模失败的坏样品）全部 🟢，退出码 0，两条模式的章节与表格逐节一致。

## 与逐条目检的关系

`verify-results-folder.py` 只回答"这个目录自洽吗"。数字**本身**可不可信（这段 q 该不该取、
Dmax 定得对不对、重建能不能用）仍然走判据类 skill：

- Rg 取点与可信度 → `assess-guinier-fit-quality`
- P(r) / Dmax / GNOM 选择 → `compute-and-validate-p-of-r`
- MW 用哪一法 → `choose-a-molecular-weight-method`
- 重建评估（a-score / NSD / 聚类）→ `evaluate-a-shape-reconstruction`
- 低 q 上翘到底是真散射还是光束挡边缘 → `assess-saxs-raw-data-quality`
