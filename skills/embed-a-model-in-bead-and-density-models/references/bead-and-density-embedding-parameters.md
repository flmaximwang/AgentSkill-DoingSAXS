# 参数、调用式与实跑数字（bead + density 嵌入）

本文件是 `SKILL.md` 的参考件：完整的命令行/API 调用式、参数表、结果字段、以及 2026-10-01
A5-05-1…6 的实跑数字。**所有数字都是本机实跑得到的**，不是推算。

---

## 1. 输入契约：样品目录里什么是靶子、什么是输入

上游流水线（`run-a-tube-saxs-pipeline-end-to-end` / `run-a-sec-saxs-pipeline-end-to-end`）每个样品写：

| 路径 | 用途 |
|---|---|
| `<样品>/profiles/03_subtracted/subtracted.dat` | 扣减曲线（q, I, err 三列）→ GNOM 的输入 |
| `<样品>/tables/ift_summary.json` | 各 IFT run 的 `idx_min/qmin/qmax/dmax/rg_realspace/chisq/gates/trusted` + `chosen`（代建珠模型就取这里） |
| `<样品>/summary.json` | 同一批数字的上层汇总（`ift.runs` / `ift.chosen` / `ift.trusted` / `shape.denss`） |
| `<样品>/models/denss.mrc` | **电子云靶子**（`denss_current.mrc` 是过程量，`denss_support.mrc` 是支撑掩膜，别拿错） |
| `<样品>/models/*dammif*` `*damaver*` `*damfilt*` | 若上游已跑过珠模型，**珠模型靶子**在这里；本次 A5 批没有（`shape.dammif` 当年写的是 "not run"），所以是代建 |
| `models/*.pdb`（**项目根**，不是样品目录） | 高分辨模型（本次：`/…/DataProcess_2026.10.01/models/A5.pdb`） |

**两条与"产物寿命"有关的实测**（2026-10-01）：

- **样品目录会被上游流水线重跑清空**：同一批数据常有并发会话在重跑 `run-raw-tube-pipeline.py`，它重建整个
  `<样品>/` 目录——把结果写进 `<样品>/embed/` 会被下一次重跑删掉。**默认输出目录因此改成同级
  `<processed>/_embed/<样品>/`**（`--out-dir` 可覆盖）。
- **上游自己也会建珠模型**：流水线带 `--model-engine auto` 时会在 `<样品>/models/` 写下
  `damaver-cluster001-damfilt.cif` 等；本 skill 检测到就直接**复用**（优先 `*damfilt*`＝滤过平均＝最可能模型，
  其次 `*damaver*`），不再重建。注意别把两个 glob 合并成一次 `sorted()`——那样 `damaver` 会排在 `damfilt` 前面。

`read_sample_meta()` 的取用顺序：`tables/ift_summary.json` → `summary.json`；只取 `chosen` 那条 run，
优先 `trusted=true`。**不要**自己另选 q 窗口——那会把上游判据废掉（见 §5 的对照实测）。

---

## 2. 代建珠模型：RAWAPI 调用式（本机实测通过）

用 RAW 自带解释器跑（`pyFAI/numba/matplotlib` 都在；**不要在 RAW 源码目录里跑**）：

```python
import sys; sys.path.insert(0, "/Applications/BioXTASRAW/lib/python3.12/site-packages")
import bioxtasraw.RAWAPI as raw

s = raw.load_settings(<项目>/data/<日期>.cfg)     # 或 raw.load_settings(None)
s.set("ATSASDir", "/Applications/ATSAS-4.1.4-1/bin")   # ★ 必须是 bin 目录
prof = raw.load_profiles([<样品>/profiles/03_subtracted/subtracted.dat], settings=s)[0]

ift, dmax, rg, i0, rg_err, i0_err, total_est, chi2, alpha, quality = raw.gnom(
        prof, dmax=<P(r) 的 Dmax>, rg=<P(r) 的 Rg>,
        idx_min=<trusted run 的 idx_min>, idx_max=<q≤qmax 的最后一个点>,
        settings=s, atsas_dir="/Applications/ATSAS-4.1.4-1/bin",
        save_ift=True, savename="<prefix>.out", datadir=<bead 输出目录>)

for i in range(1, 16):
    r = raw.dammif(ift, "<prefix>_dammif_%02d" % i, <bead 目录>, mode="Fast", symmetry="P1",
                   model_format="cif", write_ift=False, ift_name="<prefix>.out",
                   settings=s, atsas_dir="/Applications/ATSAS-4.1.4-1/bin")
    # r = (chi_sq, rg, dmax, mw, excluded_volume)

raw.damaver(files=<纯文件名列表>, prefix="<prefix>", datadir=<bead 目录>,
            symmetry="P1", model_format="cif", settings=s, atsas_dir=…/bin)
```

要点 / 坑（都实测过）：

- **ATSAS 路径有两副面孔**：RAWAPI 要 `…/bin`（它拿父目录当 `ATSAS` 变量）；直接跑 ATSAS 二进制要
  `export ATSAS=/Applications/ATSAS-4.1.4-1`（安装根）。喂错就是 `NoATSASError: 'Cannot find gnom.'`。
- **`datgnom` 没有 Dmax 参数**：命令行帮助里只有 `-r/--rg`（Rg，必填），Dmax 由它自己定
  （实测 A5-05-1：喂 `-r 20.77` → 自动 Dmax **55.6 Å**，而 P(r) 是 **72.6 Å**）。要控 Dmax 就用
  `raw.gnom(dmax=…)`。
- **`model_format='pdb'` 在 ATSAS 4 下不可靠**：DAMMIF 实际仍写出 `…-1.cif`，而 RAWAPI 按 pdb 名回读 →
  `FileNotFoundError: …/xxx_dammif_01-1.pdb`。**统一用 `cif`**（RAWAPI 也按 cif 回读），要 PDB 再用
  `cif2pdb`/CIFSUP 转。
- **`dammif` 的返回**：`(chi_sq, rg, dmax, mw, excluded_volume)`；ATSAS 4 下 `mw` 实测为 0.0（别当结果用）。
- **DAMAVER 产物命名**（ATSAS 4，`--model-format=cif` 时）：
  `<prefix>-global-damfilt.cif`（滤过平均＝最可能模型，**推荐做嵌入靶子**）、`-global-damaver.cif`
  （平均）、`-global-damstart.cif`（给 DAMMIN 精修的起点）、`-cluster00N-*.cif`（按簇分）、
  `<prefix>-distances.txt`（两两 NSD 矩阵 + `Mean value of nsd` / `Standard deviation of nsd`）、
  `<prefix>-global-summary.txt`（哪些模型 Included、哪个是 Most representative）。
  取 consensus 的优先级：`global-damfilt` > `cluster001-damfilt` > `global-damaver` > `cluster001-damaver`。
- **bead 半径与数目**：写在模型 cif 头 `_atsas_dummy_atom_model.value 2.300`（本次 2.30 Å）；
  数目随系综变（A5-05-1：669 beads@15 模型；两次 3 模型的试跑分别是 589 / 617）。
- 时间参考（本机 M 系列）：GNOM ~0.1 s；DAMMIF `Fast` **~7–8 s/个**；DAMAVER 15 个模型 ~30–60 s。

---

## 3. CIFSUP：参数表、分数在哪、实测 NSD

```
cifsup [OPTIONS] <STATIC> <MOVABLE>          # STATIC=珠模型（或不动的那个）, MOVABLE=原子模型
```

| 选项 | 说明 | 本 skill 用法 |
|---|---|---|
| `--method=` | `NSD` / `NCC` / `ICP` / `RMSD`，默认 **ICP** | 一定要显式写 `--method=NSD`（ICP 给的是另一种量纲，实测 8–9 那种数） |
| `--selection=` | `ALL`（全原子）/ `BACKBONE`（主链 CA）/ `REGRID`（两边都转成 DAM）/ `SHELL`（只取最外层原子），默认 `ALL` | 默认用 `REGRID`（官方＝SUPCOMB fast mode）；四个都跑一遍记对照 |
| `--beads=N` | 仅 `REGRID` 用，默认 **2000** | 记下实际值，因为 NSD 随它变（实测 2000→0.744，1200→0.686） |
| `-e/--enantiomorphs` | 搜对映体，**默认 YES** | 保持默认（dummy atom 模型有手性歧义） |
| `--lm/--ns/--smax` | 只在 `method=NCC` 时有效 | 不用 |
| `-o/--output=` | 输出被移动/旋转后的 **MOVABLE** 模型 | 同时写一份 `.cif`（读分用）和一份 `.pdb`（给人看） |
| `-h/--help` | — | — |

**分数在哪：** CIFSUP **没有运行时输出**（手册原文）。NSD 写在输出文件里：

```bash
cifsup --method=NSD --selection=REGRID -o out.cif <bead.cif> model.pdb   # stdout 是空的
grep -m1 '^score' out.cif        # → score                  0.874620
```
PDB 输出时同一条信息在 `REMARK 265` 段（`Template file / Movable file / Enantiomorphs searched /
Superposition method / Input atom selection` + `_atsas_superposition` 的旋转矩阵与 `score`）。

**实测（A5-05-1，与同一份 bead 模型/同一份模型对比）：**

| `--selection` | `--method=NSD` | `--method=ICP` |
|---|---|---|
| ALL | 2.202 | 8.176 |
| SHELL | 2.202 | 8.176 |
| BACKBONE | 1.141 | 8.553 |
| **REGRID** | **0.875** | 9.422 |

- 方向对称性：把 template/movable 对调，REGRID 0.744407 vs 0.744313、ALL 1.658905 vs 1.658646（差在数值误差级别）。
- **同一数据不同 DAMMIF 系综**：3 模型系综 → NSD 0.744（617 beads）；另一批 3 模型 → 0.933（589 beads）；
  15 模型 → 0.875（669 beads）。所以 NSD 只有配上"哪批珠模型"才可解释。

**判读（只有官方这一条）：** NSD→0 = 理想叠合；**>1 = 两个物体系统性地不同**（SUPCOMB 手册原文）。
本 skill 不编"合格线"。

---

## 4. ChimeraX fitmap：命令、字段、实测

```
open <denss.mrc>                      # #1 参考图（不动）
open <model.pdb>                      # #2 被嵌的模型
fitmap #2 in #1 resolution 15 metric correlation search 200 seed 42 logFits fit_r15.csv listFits true
save <out.pdb> models #2
```

| 参数 | 作用 | 坑 |
|---|---|---|
| `resolution R`（Å） | 把原子模型高斯化成图，做 **map-in-map** 拟合 | **不写就退化成 atoms-in-map**，只最大化"原子位置上的平均图值"，`correlation` 会是 **None** |
| `metric` | `overlap`（默认）/ `correlation` / `cam` | 本次用 `correlation`；报数时 `correlation` 与 `correlation_about_mean` 都给 |
| `search N` + `seed M` | N 个随机初始位置 + 局部优化 | fitmap 是**局部**优化器；`search 0` 时结果取决于模型初始位置；可复现要固定 seed（官方还要求先 `view initial`） |
| `logFits <csv>` | 把每个唯一解放进 CSV | 列名固定：`Rxx…Tz`（两次，共 24 个），`correlation correlation_about_mean overlap average_map_value points atoms_outside_contour clash contour_level steps shift angle` |
| `listFits true` | 在 Log 里打印 Top 相关值列表 | 本次输出形如 `Correlations and times found: 0.8812 (18), 0.8777 (…), …` |
| `envelope` / `zeros` / `levelInside` | 控制用哪些格点、多少比例必须在轮廓内 | 本次保持默认（`envelope` 默认 true、`levelInside` 0.1） |

**实测 R 阶梯（A5-05-1，`search 200 seed 42`）：**

| resolution | best correlation | correlation_about_mean | n_unique_fits | top3 | 相对最高分（25 Å）的姿态 CA-RMSD |
|---|---|---|---|---|---|
| 10 Å | 0.8678 | 0.5548 | 14 | 0.8678 / 0.8649 / 0.8646 | 1.70 Å |
| 15 Å | 0.8812 | 0.6997 | 18 | 0.8812 / 0.8777 / 0.8767 | 2.61 Å |
| 20 Å | 0.8782 | 0.7596 | 24 | 0.8782 / 0.8772 / 0.8762 | **28.89 Å** |
| **25 Å** | **0.8820** | 0.8046 | 21 | 0.8820 / 0.8790 / 0.8790 | 0（基准） |

- 分布很平：相关值 0.868–0.882，**但姿态不是**——20 Å 找到另一个摆法，分数只差 0.004。
- 即使同一个 R，前两名也能是不同摆法（15 Å：0.8812 vs 0.8777，差 0.0035）。
- **可复现的前提是 `view initial` + 固定 seed**：只用固定 seed 时，上面这张表在两次会话间会整体换一套值（见 §4 实测 0.8973→0.9122 vs 0.8678→0.8812）。
- 单模型 DENSS 没有 FSC → **没有分辨率可依据**，所以才用阶梯 + 姿态散布当判据（这是替代口径）。

---

## 5. 对照实测：q 窗口决定 GNOM 给不给你"可疑解"

同一条 A5-05-1 曲线、同一个 Dmax=72.58 Å：

| q 窗口 | 点数 | GNOM 实空间 Rg | χ² | 总估计 | 结论串 |
|---|---|---|---|---|---|
| 整条曲线（idx_min=0，含低 q 上翘） | 1150 | 20.52 | **338.3** | 0.4663 | **a SUSPICIOUS solution**（DISCRP 18.28，理想 0.7） |
| 上游可信窗口（idx_min=129, q≥0.0567） | 521 | 20.79 | **0.94** | 0.613 | **a REASONABLE solution** |

→ 代建珠模型时必须按上游 `chosen` 的 q 窗口喂 GNOM。A5-05-2 是一个"照做但仍然差"的例子：
它被选中的 run 是 `window-start`（Dmax 75.5 Å），GNOM χ² = **19.85**，DAMMIF 15 个模型的 χ² 也都在
**19.5–20.0**（模型再好也救不回 P(r) 本身的拟合差）——这种时候要在结果里点出来，而不是只报 NSD。

---

## 6. 出图与 `.pse`：PyMOL（ChimeraX 无头模式渲染不了）

ChimeraX `--nogui --script x.cxc` 里 `save out.png` 会抛
`LimitationError: Unable to save images because OpenGL rendering is not available`
（脚本不退出，只是**没有 png**）。所以图用 PyMOL：

```bash
env -u PYTHONPATH /Applications/PyMOL.app/Contents/bin/python3.10 render-embed-figure.py \
    <fitted.pdb> <bead.cif|denss.mrc> <out.png> <beads|density> [选项]
# 选项（默认值即交付口径）：--overlay-transparency 0.9|0.7  --level auto  --smooth auto
#   --bead-radius <auto=CIF 头里的半径>  --transparency-mode 1  --sample-color yellow
#   --overlay-color white  --frame-margin 1.1  --size 1200  --pse <路径> | --no-pse
```

**样式是交付口径，不是审美偏好**：**黑底**（`--bg black`）+ `sample` 黄色 cartoon + `beads`/`density` 半透明**白色**；
每张 png 旁边必须有一个同名 `.pse`（`cmd.save()`），打开就是同一张图。

- **背景必须黑**：白底上白色半透明包络「白对白」，珠模型与电子云都看不清（本机实测：白底渲出来的珠球几乎与背景同色）。
  改黑底后同一组参数立刻可读；两个面板同一口径，`.pse` 里存的是 `bg_color black`。
- 黑底下的透明度：珠子 **0.85**（0.92 时包络发灰发暗）、等值面 **0.7**（0.5 时包络在画面里糊成一大片灰）。

- 必须 `pymol.finish_launching(['pymol','-cq'])`（本机 CLI 二进制被 `biorazer_pymol` 挡着；会打印一条
  `ModuleNotFoundError: No module named 'biorazer_pymol'`，无害）。
- 对象名**不能用 `model`**（PyMOL 保留字，会被改名成 `model_`，后续选择器全部失效）。
- **珠球透明度有两个前置条件**（实测，2026-10-01，`sphere_transparency 0.8` 与全不透明渲染逐像素比）：
  1. 珠子认的是 `sphere_transparency`（`transparency` 只管 surface，老脚本对珠球设 `transparency` 等于没设）；
  2. 还要把 `transparency_mode` 从 PyMOL 默认的 **2** 改成 **1**——mode 2 在 ray trace 里把珠球透明度吃掉：

     | transparency_mode | 与不透明渲染的平均像素差 | 结论 |
     |---|---|---|
     | 2（默认） | **2.8 / 255** | 看着不透 |
     | 0 / 1 / 3 | **8.4 / 255** | 真的透了 |

  脚本现取 `--transparency-mode 1`。
- **珠球半径用模型自己的值**：`.cif` 头里的 `_atsas_dummy_atom_model.value`（本批 2.2 / 2.3 Å），
  实测 `sphere_scale` 0.5 的 vdW 画法会让珠子沿 4.4 Å 晶格互相吞并成米粒状（视觉上像另一种重建）。
- **透明度的实际观感**（同一张 A5-05-1 电子云图，测\"黄色像素占比 / 黄色像素的 R−B 均值\"）：

  | 等值面 transparency | 0.35 | 0.45 | 0.65 | 0.75 | 0.85 |
  |---|---|---|---|---|---|
  | 黄色像素占比 | 15.4% | 16.3% | 17.0% | 16.8% | 16.4% |
  | 黄色像素 R−B | 69.7 | 87.3 | **123.0** | 143.2 | 164.3 |

  → 包络越透，样品越饱和；取 **0.7** 时白包络仍清楚可读、黄色已经饱和（0.55 时偏灰）。
- **DENSS 等值面**：A5 批的 `.mrc` 最大 ~0.58、均值 ~0.0015、体素 **6.8 Å**（32³ / 218 Å 盒子）。
  - 直接按 6.8 Å 体素渲会出现大片平面刻面；出图前用 `--smooth auto`（= 一个体素的 Gaussian）再求等值面。
  - **平滑后必须重解 level**：同一个 level 0.02 在平滑图上包围体积会涨到 **279%**；正确口径是让包围体积
    回到 `denss.log` 的 `Final Support Volume`（本次 56392 Å³），实测 level 0.09539 → 56077 Å³ = **99.4%**。
  - 该体积由 `mrcmap.py` 直接从 `.mrc` 数格点算出（`volume_above`），不是估计值。
- **框架**（`--zoom-buffer`，默认 8 Å；包围盒按**一个体素**外扩——等值面在体素之间插值，实测 18 个样品里有 1 个
  （909-TC-APO，体素 11.6 Å、包络 1.49e5 Å³）把等值面顶到画面角上），余下同句：：珠球面板直接 `zoom(sample or beads)`；电子云面板不能用 `zoom()`
  （**surface 对象不是原子选择，PyMOL 会报 `Invalid selection name`**），改用「把包络包围盒做成 8 个临时 pseudoatom
  再 zoom」——包络体积的 voxel 下标 → XYZ 用 `cmd.get_extent(map)` 换算（实测该 extent 覆盖的是**体素中心**：32 体素 /
  217.7 Å 的盒子报 210.9 Å = 31 个间隔，所以映射是 `lo + idx/(N-1)*(hi-lo)`，不加半格）。
  **不要**用手动拉远相机（把 view 的 9/10/11 乘系数）：PyMOL 的裁剪面不会跟着走，实测 ×1.35 整幅画面变暗糊掉。
- **两个靶子两张图**：珠模型与 DENSS 图各自居中在各自原点，相对取向无约束；混画会得到"模型戳出珠球"的假象。
- 补 `.pse` 不必重跑科学步骤：`embed-model.py <样品> --out-dir <结果目录> --figures-only`
  （读现成 `embed_results.json`，只重渲两个面板，几秒一个样品）。


---

## 7. 几何复核：CA 到最近 bead 的距离

NSD 是一个归一化距离，看不住"模型是不是真在格子里"，补一条直白度量（`embed-model.py` 里没默认算，
需要时用下面的片段）：

```python
d = min_distance_each_CA_to_any_bead      # 粗判：中位数、>8 Å 的比例
```

A5-05-1 实测（669 bead 的 `damfilt` 共识模型，bead 半径 2.30 Å）：CA→最近 bead **中位数 3.35 Å**、
最大 24.6 Å、**65.5% 在 8 Å 内**。注意 `damfilt`（滤过平均，只保留多数模型都有的 bead）比 `damaver`
（未滤过平均）稀疏，换 damfilt↔damaver 或换一批 DAMMIF 模型这个比例都会动 —— 它是复核，不是判据。

---

## 8. 实跑总表（2026-10-01，A5-05-1…6，`models/A5.pdb`，DAMMIF Fast ×15）

（下表由 `embed_results.json` 汇总；`trusted` 指上游 IFT 闸门，`—` 表示该样品没有 DENSS 图/不可信。）

| 样品 | IFT run | trusted | Dmax (Å) | Rg(P(r)) | GNOM χ² | DAMMIF χ² | DAMAVER NSD | beads | CIFSUP NSD (REGRID) | NSD (ALL/BACKBONE) | fitmap 最好解 (R / correlation) | 姿态散布 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A5-05-1 | guinier-fit-start | ✅ | 72.6 | 20.79 | 0.94 | 1.37–1.97 | 0.73±0.29 (2 簇) | 669 | **0.875** | 2.20 / 1.14 | 25 Å / **0.8820** | 28.9 Å |
| A5-05-2 | window-start | ✅ | 75.5 | 20.42 | **19.85** | 19.5–20.0 | — | — | **0.775** | 1.83 / 1.03 | 20 Å / 0.8980 | 27.9 Å |
| A5-05-3 | guinier-fit-start | ❌ | 108.9 | 23.45 | — | — | 0.60±0.05 (**6 簇**) | 414 | 1.231 | 3.51 / 1.54 | 无 DENSS 图 | — |
| A5-05-4 | guinier-fit-start | ❌ | 121.1 | 28.20 | — | — | — | — | 1.973 | 6.27 / 2.52 | 无 DENSS 图 | — |
| A5-05-5 | guinier-fit-start | ❌ | 174.2 | 41.84 | — | — | — | — | **4.426** | 12.17 / 4.65 | 无 DENSS 图 | — |
| A5-05-6 | guinier-fit-start | ✅ | 88.9 | 25.28 | — | 0.88–0.89 | 0.76±0.25 (2 簇) | 492 | 1.931 | 5.44 / 2.26 | 25 Å / 0.9061 | 26.2 Å |

读这张表的三个要点：

1. **同一个模型、同一批稀释序列，NSD 从 0.775 到 4.426**：可信的三个样品里 A5-05-1/2 落在 <1（0.875 / 0.775），
   A5-05-6 是 1.93（>1 = 系统性不同）；四个未过闸门的样品 NSD 全部 ≤1 之外（1.23 / 1.97 / 4.43）。
   NSD 在这里表现为"这条曲线给出的包络里能不能塞进这个模型"，而不是"模型对不对"。
2. **A5-05-2 是一个"照做但仍然差"的例子**：它的 GNOM χ² = 19.85，DAMMIF 15 个模型也全在 19.5–20.0 ——
   模型再好也救不回 P(r) 本身的拟合差；这种样品要在报告里点出来。
3. **电子云分支的相关值全在 0.89–0.91**（三个有图的样品），而**姿态散布全在 26–29 Å**：
   说明这套包络对"摆在哪"的约束很弱，报分时必须连姿态散布一起报。

> 每个样品的完整数字在 `<样品>/embed/embed_results.json`；本表的 A5-05-1 行是本次实跑的一手数据，
> 其余行以 JSON 为准（同一批跑的）。
