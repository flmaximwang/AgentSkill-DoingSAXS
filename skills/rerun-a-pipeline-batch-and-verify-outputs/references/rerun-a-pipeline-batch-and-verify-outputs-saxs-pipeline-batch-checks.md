# SAXS 流水线整批重跑：具体路径与核对清单

## 布局

- 原始帧：`<项目>/data/<模式>/<样品>/*.tif`（+ 每帧 BL19U2 header `.txt`）
- 结果：`<项目>/processed/<模式>/<样品>/{frames,ifts,models,norm,profiles,reports,tables}`
- 旁支：`_summary/{summary.csv,summary.md,overview.png}`、`_logs/<样品>.log`、`_prev_<日期>/`、`_embed/<样品>/`
- 解释器固定用 `/Applications/BioXTASRAW/bin/python`（系统 python 没有 RAW/ATSAS 依赖）
- ATSAS 的 `bin` 目录用 `/Applications/ATSAS-<版本>/bin`；RAW 自己找 ATSAS 走 `~/ATSAS-<版本>` 符号链接（官方找法：`ln -s /Applications/ATSAS-<版本> ~/ATSAS-<版本>`）

## 管式（tube）流水线

脚本：`~/.hermes/skills/saxs/run-a-tube-saxs-pipeline-end-to-end/scripts/`

- 逐样品：`run-raw-tube-pipeline.py --sample-dir <data/…/<样品>> --cfg <data/…/*.cfg> --out-dir <processed/…/<样品>> --model-engine auto --denss-mode Fast --n-models 4`
- 汇总：`summarize-tube-run.py <结果根>` → `_summary/summary.{csv,md}`（31 列）
- 总览图：`plot-tube-overview.py <结果根>` → `_summary/overview.png`
- 根 README（自动分「小蛋白 / 大颗粒」两类）由汇总这步产出；每个样品目录另有自己的 `README.md`

跑完必做的核对：

1. 逐日志找失败：`for f in _logs/*.log; do grep -nE 'failed|失败|not written|Traceback' $f; done`
2. 逐样品核对 README 是否齐：`for d in */; do [ -f "$d/README.md" ] || echo "缺 README: $d"; done`
3. 汇总样品数 = 数据目录样品数（代码数，不是眼睛数）

已知失败信号（内部节点失败，退出码仍 0）：

- `3D DENSS failed: DENSS failed to run properly`，紧邻 scipy `The number of derivatives at boundaries does not match` → 曲线本身坏（低 q 上翘 / IFT 未过闸门）
- `珠模 DAMMIF #N 失败: FileNotFoundError …/dammif_0N-1.cif` → DAMMIF 没写出模型，是同一坏曲线的下游
- `WARNING: README.md not written: 'str' object has no attribute 'get'` → 生成器把「电子云失败」的错误字符串当 dict 用；该样品没有 README

这三条同源：**数据没过闸门，不是漏跑**。汇报时把这类样品的「哪些节点没出、为什么」单独列一节。

## 嵌入（embed）落点与口径

- 约定落点 `<结果根>/_embed/<样品>/`（**不是**样品目录里——样品目录会被下次流水线重跑整批替换）；脚本 argparse 默认却是 `<样品>/embed`，所以必须显式 `--out-dir`。
- 人看的 README 由**单独脚本** `write-embed-readme.py` 生成，入参是**目录**（给结果根的 `_embed` 一次刷全部；给 `embed_results.json` 路径会静默什么都不做）。
- 珠模型默认复用流水线已产出的 DAMAVER 共识（`models/` 下），不重打；`--rebuild-bead` 才会自己重新 GNOM→DAMMIF×N→DAMAVER。
- 若跑的时候样品还没有 `models/denss.mrc` 或 IFT 没过闸门，脚本仍会写出「no DENSS map」/「IFT was NOT trusted」的**空壳 README + embed_results.json**——结果树里看到这种空壳 = 上一次跑在流水线之前，直接覆盖即可，别当有效结果。
- 同一份数据在结果树里可能有**两套**（一次是第一遍按默认落点跑进样品目录的、一次是后来按约定落点跑进 `_embed/` 的）。认身份靠 `embed_results.json` 里的 `out_dir` 字段（写明是哪次调用写的）+ 文件 mtime；只看目录名会以为是同一份东西，进而把两套数字混着报。
- **报 NSD 必须写明用的是哪个珠模型文件**：同一 `models/` 下同时有 `*-global-damfilt.cif`（滤过平均≈最可能模型，bead 数约为平滑平均的一半）、`*-damaver.cif`（平滑平均）与 `*-cluster00N-*` 变体；同一份数据同一曲线换文件，NSD 可差 ~0.2（实测 0.903 vs 1.143，bead 数 589 vs 1177）。口径跟着文件走，不写文件名的 NSD 没有可比性。
