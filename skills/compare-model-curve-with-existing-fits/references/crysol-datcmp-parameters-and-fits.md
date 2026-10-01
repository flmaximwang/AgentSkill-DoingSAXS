# CRYSOL / DATCMP / datcrop 参数与文件格式速查（ATSAS 4.1.4）

> 全部来自 ATSAS 4.1.4 官方手册（`biosaxs-com/atsas` 4.1.4：manuals/crysol、manuals/datcmp、manuals/datcrop、file-formats/{dat,fit,fir,abs,int}），
> 加本机（ATSAS 4.1.4-1，2026-10-01）在 `DataProcess_2026.10.01` 上的实测数字。**实测值与手册值分开标注**。

## 1. crysol 命令行

用法：`crysol [OPTIONS] [FILE(S)]`，FILE(S) = ≥1 个坐标文件（`.pdb`/`.cif`/`.ent`）+ 可选 0..n 个实验数据（`.dat`）。

| 选项 | 默认 | 单位 / 范围 | 什么时候改 |
|---|---|---|---|
| `--lm N` | **20** | 1–100，球谐最大阶数 | 高长径比物体调大；**`--smax` 调大时必须跟着调大**（手册原话）；球状体系改它没用（实测 4DH2-676-apo-dimer：`--lm=100` 与 20 结果完全一致） |
| `--fb N` | 17 | 10–18，Fibonacci 网格阶 | 只在 `--shell=directional` 下有效；提精度、代价是 CPU |
| `--ns N` | 101 | ≤10001 | 计算数据点个数（理论曲线的 q 点数） |
| `--smax S` / `-sm` | **0.5 Å⁻¹** | ≤2.0 | 拟合/计算的最大 q；数据上界更高时要显式给，否则理论曲线"够不到" |
| `-u/--unit` | 猜 | `u`/1/2/3/4 | 1 = `s=4πsinθ/λ` in Å⁻¹；2 = nm⁻¹；3 = `2sinθ/λ` in Å⁻¹；4 = `2sinθ/λ` in nm⁻¹。**启发式只认 1 和 2**：`s_max<1.0` 判 Å⁻¹，否则 nm⁻¹；WAXS 数据必须显式指定 |
| `--dns V` | 0.334 e/Å³ | — | 溶剂电子密度（纯水）；高盐可略高 |
| `--dro V` | 0.03 e/Å³ | — | 水化层对比度；**打开最小化时这两个会被拟合**，要手改得先 `--skip-minimization` |
| `--constant` / `-cst` | 关 | — | **常数扣减**：吸收"缓冲液失配"这类系统偏差。实测同一模型同一数据：关 4.233 → 开 **1.267** |
| `--skip-minimization` | 关 | — | 用给定的 dns/dro 直接算，不拟合 |
| `--shell` | directional | `directional`(经典 CRYSOL) / `water`(旧 CRYSOL3) | 换水化层模型；CRYSOL3 已并入 CRYSOL，故用 `--shell=water` |
| `--model ID` / `--chain ID` | all | — | NMR 多构象选 model；多链 `.cif` 选链（也可只算某条链来比） |
| `--explicit-hydrogens` | 关 | — | 文件里已含 H 时用，避免去查氢原子数 |
| `--implicit-hydrogen N` | — | ≥0 | 覆盖 `unable to determine number of hydrogens for …` 报错 |
| `--sub-element X` | — | — | 覆盖 `unable to determine element for …` |
| `--alternative-names` | 关 | — | 老 PDB 的原子命名 |
| `--energy eV` | 无 | — | 只在反常 SAXS 的能量校正用 |
| `-p/--prefix` | 输入基名 | — | **给了前缀后，两次日志（算曲线 + 拟合）合并成一个 `<prefix>.log`**，产物名也变成 `<prefix>.fit` |

**拟合都拟合了什么**（手册 Introduction）：三个参数 —— ①总置换溶剂体积 ②水化层对比度 ③相对本底。所以"生成理论曲线"和"拟合数据"是同一件事。

## 2. crysol 的输出与文件结构

| 后缀 | 类型 | 内容 |
|---|---|---|
| `.log` | ASCII | 屏幕输出副本（参数 + `Chi-square of fit` + `Probability of fit`） |
| `.int` | ASCII | 理论曲线，任意刻度 |
| `.abs` | ASCII | 理论曲线，绝对刻度 1/(cm·(mg/ml)) |
| `.fit` | ASCII | **拟合**：理论 vs 实验。第 1 行字段 = `Dro`（最优水化层对比度）、`RGT`（理论曲线 Rg, Å）、`Vol`（最优排除体积, Å³）、`Chi^2` |
| `.alm` | **二进制** | 球谐振幅（中间量，可删） |

- `.fit` 的数据列 = ①q ②I_exp ③Err ④I_fit(**已含常数项**，实测：把常数再加一次 χ² 从 1.267 变成 10.65)。
- `.fir`（DAMMIF/dammin 等写的）= 同四列，表头行是 `sExp | iExp | Err | iFit(+Const) | Chi^2=`。
- **不带 `.dat` 调用**：只有 `.int`/`.abs`/`.alm`/`.log`，**没有 `.fit`、没有 Chi^2** —— 这是 minimal 理论曲线，不是拟合。
- 一个模型一次调用 → 产物名不打架、日志无歧义（本 skill 的脚本就是这么做的；手册的多模型示例是一次调用多个模型，效果等价）。

## 3. datcrop（裁 q 区间）

`datcrop [OPTIONS] [SASDATA]`：`--smin V` / `--smax V`（**注意是 s 不是 q**）、`--first N` / `--last N`、`-o FILE`（默认 stdout）。
输出文件顶部写 `Sample description:`、`Parent(s): <源文件>`、`datcrop-smin/smax` —— 保留来源痕迹，别把它们当数据行。

## 4. datcmp

`datcmp [OPTIONS] <SASDATA(S)>`：`--mode PAIRWISE|INDEPENDENT`、`--test CORMAP|CHI-SQUARE|ANDERSON-DARLING`、`--adjust FWER|FDR`、`--alpha 0.01`、`-f FULL|CSV`、`--prefix`。

**单文件用法（本 skill 用这条）**：

| 喂进去的文件 | datcmp 算什么 |
|---|---|
| `.fit` / `.fir` | 文件里的 **实验数据 vs 模型拟合** ← 评拟合 |
| `.out`（GNOM 正则化） | p(r) 重建散射 vs 实验数据 |
| `.dat`（普通 3/4 列） | **没有意义**：实测 `χ²=0.000000, p=1.000000`（自己跟自己比），且把同一文件同时标成 1 和 2 |
| `.ift`（BIFT） | 算不出：哨兵 `-1.000000` |

**三个检验**（手册 Statistical Tests）：

| 检验 | 用途 | 判据值 | 性质 |
|---|---|---|---|
| CorMap | 残差的**随机性**（用符号） | `C`（p） | 不需要误差估计；对成片的形状偏差敏感、对逐点噪声不敏感；**曲线越长同样的 C 越显著** |
| red. χ² | 标准化残差的**大小** | `Chi^2`（p） | **要求误差准确**；"接近 1"随 n 变（n≈50 → 0.5–1.5；n≈2000 → 0.9–1.1）；精确区间可由 χ²_n 分布算 |
| Anderson-Darling | 标准化残差的**尾部**是否正态 | `A²`（p） | 比 red. χ² 更敏感于尾部偏差；同样要求误差准确 |

标准化残差定义：`r(s) = (I₁(s) − I₂(s)) / sqrt(σ₁² + σ₂²)`；相似时 `r` 应近似标准正态、落在 ±3 内。

CSV 输出一行：`<测试名>, <i>, <j>, <值>, <p>, <p_adj>`（`--format=csv`），本 skill 取第 4/5/6 个字段。

**手册的判读示例**（Advancing Model Fit）：`datcmp ly01.fir` → `C=12.000000, Pr(>C)=0.044513`；取 α=0.01 时 p>α → **不能拒绝"模型拟合数据"**。

本 skill 用的判档（**依据写在这里**）：

- 区间 = χ²_n 的 95% 区间 = `[chi2.ppf(0.025,n)/n, chi2.ppf(0.975,n)/n]`（手册说"精确区间可由 χ²_n 分布算出"，这是那句话的实现；手册自己的两个锚点 n≈50→0.5–1.5、n≈2000→0.9–1.1 与之相容：我们算的 n=50 → [0.65,1.43]、n=2000 → [0.94,1.06]，比手册的经验区间窄，因为手册那把尺子更松）。
- 🟢 = red. χ² 落在区间内 **且** CorMap `p > alpha`（默认 α=0.01）。
- 🟡 = 落在区间内但 CorMap 拒绝；或超出区间不到 2 倍。
- 🔴 = 超出区间 2 倍以上。
- 实测参考值：n=519 → [0.88, 1.13]；n=845 → [0.91, 1.10]；n=1149/1150 → [0.92, 1.08]。
- **绝对判档依赖误差标定**：SEC 扣减后的 σ 偏乐观（逐帧独立假设）→ 同一个模型在 SEC 曲线上更难达标；**排序（相对）比绝对判档可靠**。

## 5. 三列 χ² 的实测对照（同一批数据 `SEC-SAXS/4DH2-676-apo-3`）

| 曲线 | 文件里自报 | datcmp 重算 | 本 skill 插值到实验网格 | 点数 | q 网格 |
|---|---|---|---|---|---|
| crysol 拟合（`4DH2-676-apo-dimer`，`--constant`） | 1.267 | **1.2668** | 1.2678（无需插值，网格相同） | 1150 | 实验点网格 |
| DENSS（`*_denss_map.fit`） | 2.665（= 除以 N） | **2.668** | 4.3996（**虚高**，重采样网格） | 845 | 高 q 段更稀 |
| DAMMIF（`*_dammif_01.fir`） | 4.334 | **4.3336** | 4.3336 | 1149 | 重采样成等距 3.89e-4 |
| GNOM（`*_gnom.out`） | IFT 表里 0.4631 | 4.2526 | —（`.out` 不按四列读） | — | 拟合块含被裁掉的低 q 与 q=0 |
| BIFT（`*_bift.ift`） | 1.247（IFT 表） | −1（哨兵） | — | — | — |

网格实测：原始 1150 点等距 `3.8941e-4 Å⁻¹`；DAMMIF 1149 点等距 `3.89e-4`（与原始点最大差 `4.95e-7 Å⁻¹`）；DENSS 845 点。

## 6. 运行环境（本机）

- 解释器：`/Applications/BioXTASRAW/bin/python`（含 numpy/scipy/matplotlib；**不要在 RAW 源码目录里跑**，`sascalc_exts` 会被遮蔽）。
- ATSAS：`/Applications/ATSAS-4.1.4-1/bin`。脚本按 参数 → `$ATSAS` → 常见路径 找；找到后若无 `$ATSAS` 就用 bin 的父目录补上（crysol 靠它找 CDD 等资源）。
- 图内文字用英文：matplotlib 默认字体没有中文字形，中文会画成方框（控制台消息仍可中文）。

## 7. 已知报错与处置（手册 Common Issues + 实测）

| 报错 / 现象 | 原因 | 处置 |
|---|---|---|
| `unable to determine number of hydrogens for …` | 残基名对不上 CDD | `--alternative-names` / `--explicit-hydrogens` / `--implicit-hydrogen N` |
| `unable to determine element for …` | 老文件原子名/元素缺失 | `--alternative-names` / `--sub-element X` |
| `No such model or chain in structure: xxx.cif`（**实测**：`4LI2-676_D.cif` 里只有水） | 该文件/链里没有原子 | 换文件或用 `--chain`；多链体系先拆链核对（拆链脚本在 `embed-a-model-in-bead-and-density-models`） |
| `datcmp` 给 `χ²=0, p=1` | 喂了普通 `.dat` | 只喂 `.fit`/`.fir`/`.out` |
| `datcmp` 给 `-1` | 文件不是可识别的拟合格式（BIFT `.ift`；DAMMIF 的 3 列 `.fit` 绘图文件） | 换 4 列的 `.fir`；把该行标为"不可判" |
| 同一拟合两个不同的 χ²（如 GNOM 0.46 vs datcmp 4.25） | `.out` 的拟合块覆盖全 q（含被裁掉的低 q 与 q=0），GNOM 自报的是它自己拟合区间/误差模型下的值 | 说明口径；该行不参与跨程序排序 |
