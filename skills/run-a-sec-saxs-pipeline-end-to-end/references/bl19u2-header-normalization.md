# BL19U2 逐帧 header txt → RAW 归一化：机制、证据、核对法

本文件是 `run-a-sec-saxs-pipeline-end-to-end` 的支撑材料：为什么"补一份 txt"就能让 RAW 做逐帧归一化、
怎么证明它真的生效、以及监视器 → txt 这一步的每个参数凭什么这么取。**行号对应 BioXTAS RAW 2.4.2 源码**
（`~/Repositories/bioxtasraw`，装在本机 `/Applications/BioXTASRAW` 环境里）。

## 1. 机制链（源码级）

| 环节 | 位置 | 做什么 |
|---|---|---|
| header 格式注册 | `SASFileIO.py:968` `'BL19U2, SSRF' : parseBL19U2HeaderFile` | RAW 的 header 格式下拉里本来就有这个选项（还有 I711/I911-4/MAXLab、F2/G1/CHESS、BioCAT/APS、P12/Petra III） |
| 找哪个文件 | `SASFileIO.py:909` `parseBL19U2HeaderFile` | `fname, ext = os.path.splitext(filename); countFilename = fname + '.txt'` —— **图像路径去掉扩展名 + `.txt`**，即必须和 tif 同名同目录 |
| 怎么解析 | 同上 | 逐行 `line.split(':')`，第一段当 key、其余拼回当 value（所以 `Description:` 里的冒号不会切错），空值也是合法 value |
| 什么时候解析 | `SASFileIO.py:1284` `loadImageFile` → `:1334/:1360` `loadHeader(filename, new_filename, hdr_fmt)` | 载入图像时**自动**按 settings 的 `ImageHdrFormat` 找 header，不需要显式调用（API 里也有 `RAWAPI.load_counter_values`，见 `RAWAPI.py:489`） |
| 挂到数据上 | `RAWAPI.integrate_image` → `parameters['counters']` | counters 字典挂到 profile 参数里（`p.getParameter('counters')` 可见）——本机实测 `counters.TB=2.568050e-08` 与 txt 里写的完全一致 |
| 归一化求值 | `SASImage.py:165-166, 366-380` | `EnableNormalization` + `NormalizationList`：每一项是 `(op, expr)`；expr 交给 `calcExpression(expr, img_hdr, file_hdr)`（`SASImage.py:54`，用 `SASParser.PyMathParser` 把 header 的 key 当变量注入）求值 |
| 合并成因子 | `SASImage.py:366-380` | 全是 `/` `*` 时合并成单个 `norm_factor` 一次性乘；含 `+`/`-` 时走 `:520-545` 对每条曲线 `scaleRawIntensity`/`offsetRawIntensity` |

线站下机的 `.cfg` 里通常**已经写好**这一项（本机 `20261001.cfg`）：

```json
"EnableNormalization": false,
"ImageHdrFormat": "None",
"NormalizationList": [["/", "Transmitted_Beam"]],
```

→ 只有 `ImageHdrFormat` 与 `EnableNormalization` 两个值需要改（本 skill 的 `make_settings()` 干这件事）。

## 2. 怎么证明归一化真的生效（可复现的对照法）

同一批帧，开/关 `EnableNormalization` 各积分一次，比较 `I(q0)`：

```
EnableNormalization=False: bsa_00003.tif  I(q0)=5.851972e+02
EnableNormalization=True : bsa_00003.tif  I(q0)=2.284594e+10
比率 = 3.904e7 = 1 / 2.561493e-08 = 1 / (该帧 txt 里的 Transmitted_Beam)
```

逐帧比对：5 帧的比率分别等于 1/TB_i（相对偏差 <1e-6）。**判据：比率严格等于 1/TB 就是对的；
只差一个常数倍说明你比较的是"归一化到不同参考值"的两份数据，不是归一化没生效。**

## 3. 监视器 → 每帧 Transmitted_Beam 的每个参数

| 参数 | 取值依据（本机 BL19U2 bsa 实测） |
|---|---|
| 采样密度 | `bsa_1.Iochamber` 19651 行 / 3021 s = **6.50 点/s**；帧间隔 1.51 s（曝光 1.5 s + 空白 0.01 s）→ **每帧窗口约 9.8 个采样点** ⇒ 行 ≠ 帧，必须按时间窗口取中位数 |
| 窗口 | `[endTime + lag·interval − exposure, endTime + lag·interval]`；`exposure` 定窗口宽度（1.5 s），`interval` 推进时间轴与换算 lag（1.51 s）——两者不可混用 |
| lag（时钟偏移） | 两端时钟不同步：本机扫出 **−32 帧 ≈ −48 s，corr 0.869**（判据 = 检测器总计数 vs 监视器的相关系数，逐 lag 扫 ±60 帧） |
| 野值过滤 | 开头 5 行是未开束流的 `-1.35e-13`、`1.43e-13`×4，正常值 `2.576e-08`；按"低于 5%×全局中位即丢弃"处理（1e-10~1e-9 区间实测 0 个采样点，不存在"半开束流"的中间态） |
| 无采样窗口 | 丢完野值后个别窗口会空（本机 1 帧）→ 线性插值，并打印条数 |
| 平滑 | 监视器逐点抖动是仪表噪声，用 15 帧滑动中位数（`--smooth`）再取窗口值；直接把逐点值写进 txt 会把噪声灌进强度 |
| `--tb-scale` | 写进 txt 的值 = 监视器值 × 常数（默认 1e8）。归一化只用**比值**，常数不影响任何相对结果；取 1e8 是为了让强度落在与管式数据同量级，避免 profile 出现 1e10 这种数 |

**尺度提醒**：不做绝对刻度时，归一化常数会被吸收进任意刻度；若哪天要做绝对刻度，
应按 `put-saxs-data-on-an-absolute-scale` 的两条顺序约束重算标准，而不是沿用这里的 1e8。

## 4. 为什么不能"按行号对齐帧"

监视器 6.5 点/s、帧 0.66 帧/s，两者比例不是整数且会随采集抖动漂移；再加上两端时钟有 48 s 的系统偏移。
按行号取（例如第 i 帧取第 i 行）等于把监视器当成"每帧一个数"，误差随帧号线性增长。本机实测：
按行号对齐得到的"因子"与按窗口对齐的相关系数只有 ~0.2（异常），按窗口对齐是 0.87。

## 5. 与 `correct-sec-saxs-baseline` 的关系

那边的 `frame_flux_normalization.py` 走的是**另一条路线**：算完因子后**写归一化后的 int32 tif**
（适合已经决定要落一份归一化数据的场合，代价是几十 GB）。本 skill 的路线是**只写 txt、让 RAW 在积分时归一化**，
一个归一化 tif 都不落。两条路线的因子算法是同一套（同一份监视器解析 + 窗口中位数 + lag 标定），
差别只在"归一化发生在磁盘上还是发生在 RAW 内"。
