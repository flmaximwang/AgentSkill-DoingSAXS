#!/usr/bin/env python3
"""用 CRYSOL(ATSAS) 从 model 算 SAXS 理论曲线、拟合到扣减后的 .dat，并与 processed/ 里已有的拟合同口径比较。

三件事：① 跑 ATSAS 的 crysol（可选先用 datcrop 裁 q 区间）得到 <模型>.fit；
② 把每个拟合文件（本 skill 刚算的 + processed/ 里已有的 DENSS/DAMMIF/IFT）交给 ATSAS 的 datcmp
做同一套统计检验，并在**同一条 q 网格、同一列误差**上重算 reduced χ²；
③ 落 comparison.csv/json + 叠加图 + README.md（分点总结 + 明细表 + 🟢🟡🔴）。

参数、默认值与单位由 argparse 生成，见 --help；不在 docstring 里重抄。
产物布局见 --out 目录下的 README.md。
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------- ATSAS 定位

ATSAS_FALLBACKS = [
    "/Applications/ATSAS-4.1.4-1/bin",   # 本机安装位置
    os.path.expanduser("~/ATSAS-4.1.4-1/bin"),
    "/Applications/ATSAS/bin",
]

# 已有拟合文件在 processed/<模式>/<样品>/ 下的常见落点（旧名 → 标签）
FIT_PATTERNS = [
    ("models/*denss*.fit", "DENSS 电子云"),
    ("models/*denss*.fir", "DENSS 电子云"),
    ("models/*dammif*.fir", "DAMMIF 珠模"),
    ("models/*damaver*.fir", "DAMAVER 共识珠模"),
    ("models/*.fir", "珠模（未识别程序）"),
    ("ifts/*.out", "GNOM IFT"),
    ("ifts/*.ift", "BIFT IFT"),
]

# datcmp 只认 .out/.fit/.fir（见 ATSAS 手册）；喂普通 .dat 会静默给出 χ²=0/p=1.0 的假"完美拟合"
DATCMP_OK_SUFFIXES = {".out", ".fit", ".fir"}

# 这些后缀不是四列拟合文件，列不能按 SASDATA 读
NON_COLUMN_SUFFIXES = {".out", ".ift"}


def find_atsas(explicit: str | None) -> Path:
    """按 参数 → $ATSAS → 常见安装位置 的顺序找一个含 crysol+datcmp+datcrop 的 bin 目录。"""
    cands: list[str] = []
    if explicit:
        cands.append(explicit)
    if os.environ.get("ATSAS"):
        env = Path(os.environ["ATSAS"])
        cands += [str(env / "bin"), str(env)]
    cands += ATSAS_FALLBACKS
    for c in cands:
        p = Path(c)
        if (p / "crysol").exists() and (p / "datcmp").exists():
            return p
    raise SystemExit(
        "找不到 ATSAS（要 bin 目录里同时有 crysol 与 datcmp）。"
        "用 --atsas-dir 指定，或 export ATSAS=<安装根>。试过：" + ", ".join(cands)
    )


def atsas_env(bindir: Path) -> dict:
    """ATSAS 的二进制靠 $ATSAS 找资源文件（CDD 等）；没设就按 bin 的父目录补上。"""
    env = dict(os.environ)
    root = bindir.parent
    if not env.get("ATSAS") or not Path(env["ATSAS"]).exists():
        env["ATSAS"] = str(root)
    return env


def run(cmd: list[str], cwd: Path | None, env: dict) -> tuple[int, str]:
    proc = subprocess.run(
        cmd, cwd=str(cwd) if cwd else None, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    return proc.returncode, proc.stdout


# ---------------------------------------------------------------- 曲线读写


def load_sasdata(path: Path) -> np.ndarray:
    """读 3 列（q, I, err）或 4 列（q, I, err, I_fit）SASDATA；忽略注释/元数据行。"""
    rows: list[list[float]] = []
    for line in Path(path).read_text(errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith(("#", "!", ";", ">")):
            continue
        parts = s.split()
        if len(parts) < 3:
            continue
        try:
            vals = [float(x) for x in parts[:4]]
        except ValueError:
            continue  # 头部文字行（sExp | iExp …）、datcrop 的 "Parent(s): …"
        rows.append(vals)
    if not rows:
        raise SystemExit(f"没读到数值行：{path}")
    width = min(len(r) for r in rows)
    arr = np.array([r[:width] for r in rows], dtype=float)
    return arr


def fit_header_value(path: Path, key: str = r"chi\s*\^?\s*2") -> float | None:
    """从文件头几行里抓程序自己写的 χ²（crysol/dammif `Chi^2=`、DENSS `chi2=`）。"""
    try:
        head = Path(path).read_text(errors="replace").splitlines()[:3]
    except OSError:
        return None
    pat = re.compile(key + r"\s*[:=]\s*([-+0-9.eE]+)", re.I)
    for line in head:
        m = pat.search(line)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                return None
    return None


def ncols(path: Path) -> int:
    return int(load_sasdata(path).shape[1])


def align(qf: np.ndarray, yf: np.ndarray, qt: np.ndarray) -> tuple[np.ndarray, bool]:
    """在 qf 网格上取 qt 处的强度。qt 的每个点都能在 qf 里找到（1e-6 相对容差，ATSAS 只写 7 位有效数字）时
    直接取值；否则插值。返回 (值, 是否需要插值)。crysol/DAMMIF 的拟合写在实验点网格（或其子集）上，
    实测：不插值时 χ² 与 datcmp 完全一致，插值会凭空抬高 56%。"""
    qf = np.asarray(qf, float)
    yf = np.asarray(yf, float)
    qt = np.asarray(qt, float)
    if len(qf) >= 2:
        j = np.clip(np.searchsorted(qf, qt), 1, len(qf) - 1)
        left, right = qf[j - 1], qf[j]
        k = np.where(np.abs(qt - left) <= np.abs(right - qt), j - 1, j)
        tol = 1e-6 * np.maximum(np.abs(qt), 1e-6)
        if np.all(np.abs(qf[k] - qt) <= tol):
            return yf[k], False
    out = np.interp(qt, qf, yf)
    pos = yf > 0
    if pos.sum() >= 2:
        qp, yp = qf[pos], yf[pos]
        m = (qt >= qp[0]) & (qt <= qp[-1])
        if m.any():
            out[m] = np.exp(np.interp(qt[m], qp, np.log(yp)))
    return out, True


def resample(qf: np.ndarray, yf: np.ndarray, qt: np.ndarray) -> np.ndarray:
    """把 (qf, yf) 采到 qt 上。网格相同就原样返回 —— crysol/DAMMIF 的 .fit/.fir 本来就写在实验点网格上，
    插值只会凭空引入误差（实测：crysol 同一文件不插值时 χ² 与 datcmp 完全一致，插值后偏高 56%）。
    网格不同时在 log 强度上插值（曲线跨几个数量级），只在正强度区间做，其余回退线性。"""
    yf = np.asarray(yf, float)
    if same_grid(qf, qt):
        return yf
    out = np.interp(qt, qf, yf)
    pos = yf > 0
    if pos.sum() >= 2:
        qp, yp = qf[pos], yf[pos]
        m = (qt >= qp[0]) & (qt <= qp[-1])
        if m.any():
            out[m] = np.exp(np.interp(qt[m], qp, np.log(yp)))
    return out


def chi2_regrid(data: np.ndarray, fit: np.ndarray) -> dict:
    """在实验点网格 + 实验误差列上重算 reduced χ²（交叉核对用的「同口径」数）。"""
    q, I, err = data[:, 0], data[:, 1], data[:, 2]
    qf, If = fit[:, 0], fit[:, 3]
    lo, hi = max(float(q.min()), float(qf.min())), min(float(q.max()), float(qf.max()))
    m = (q >= lo) & (q <= hi) & (err > 0)
    n = int(m.sum())
    res = {"n_common": n, "q_common_min": lo, "q_common_max": hi, "regridded": None,
           "chi2_regrid": None}
    if n < 3:
        return res
    vals, regridded = align(qf, If, q[m])
    res["regridded"] = bool(regridded)
    r = (I[m] - vals) / err[m]
    ss = float((r * r).sum())
    res["chi2_regrid"] = ss / (n - 1)
    res["chi2_regrid_over_N"] = ss / n
    return res


def chi2_band(n: int) -> tuple[float | None, float | None]:
    """手册说 red.χ² 的合理区间可由 χ²_n 分布算出（n=50→0.5–1.5；n=2000→0.9–1.1 是它的两个锚点）。

    这里用 χ²_n 的 95% 区间：比手册的"经验区间"窄，但可以给任意 n 一个一致的答案。
    """
    if n is None or n < 2:
        return None, None
    try:
        from scipy.stats import chi2 as _chi2
    except ImportError:
        return None, None
    return float(_chi2.ppf(0.025, n) / n), float(_chi2.ppf(0.975, n) / n)


# ---------------------------------------------------------------- crysol


LOG_LINE = re.compile(r"^(?P<key>[^:]+?)\s*\.{2,}\s*:\s*(?P<val>.*?)\s*$")


def parse_crysol_log(text: str) -> dict:
    """crysol 的 .log 里同一批字段会出现两次（先算曲线、再拟合数据）；取最后一次的值。"""
    kv: dict[str, str] = {}
    for line in text.splitlines():
        m = LOG_LINE.match(line)
        if m:
            kv[m.group("key").strip()] = m.group("val")
    def num(key: str) -> float | None:
        v = kv.get(key)
        if v is None:
            return None
        try:
            return float(v)
        except ValueError:
            return None
    return {
        "raw": kv,
        "chi2_fit": num("Chi-square of fit"),
        "prob_fit": num("Probability of fit"),
        "model_rg_A": num("Rg (Atoms - Excluded volume + Shell) [A]"),
        "model_rg_slope_A": num("Rg from the slope of net intensity [A]"),
        "model_mw_Da": num("Molecular weight [Da]"),
        "model_volume_A3": num("Excluded Volume [A^3]"),
        "harmonics": num("Number of spherical harmonics"),
        "smax_fit": num("Maximum alm angle of experimental data [A^-1]"),
        "n_exp": num("Number of experimental data points"),
        "scale": num("Scaling factor of calculated intensities"),
        "constant": num("Added constant to calculated intensities"),
        "n_atoms": num("Total number of atoms read"),
        "n_residues": num("Number of residues read"),
        "constant_enabled": (kv.get("Enabled constant for fit", "") or "").lower().startswith("y"),
    }


def run_crysol(bindir: Path, env: dict, model: Path, data: Path, outdir: Path, a: argparse.Namespace) -> dict:
    """一个模型一次调用：crysol 会把算的曲线与拟合日志合并成 <prefix>.log，产物留在 outdir。"""
    stem = model.stem
    work = outdir / "crysol" / stem
    work.mkdir(parents=True, exist_ok=True)
    cmd = [str(bindir / "crysol")]
    if a.lm is not None:
        cmd += [f"--lm={a.lm}"]
    if a.smax is not None:
        cmd += [f"--smax={a.smax}"]
    if a.constant:
        cmd += ["--constant"]
    if a.shell:
        cmd += [f"--shell={a.shell}"]
    if a.unit:
        cmd += [f"--unit={a.unit}"]
    if a.chain:
        cmd += [f"--chain={a.chain}"]
    if a.model_id:
        cmd += [f"--model={a.model_id}"]
    cmd += [f"--prefix={stem}", str(model.resolve()), str(data.resolve())]
    rc, log = run(cmd, work, env)
    (work / "crysol.stdout.txt").write_text(log)
    logfile = work / f"{stem}.log"
    stats = parse_crysol_log(logfile.read_text(errors="replace") if logfile.exists() else log)
    fitfile = work / f"{stem}.fit"
    rec = {
        "label": stem,
        "kind": "本 skill 用 CRYSOL 算的模型曲线",
        "model": str(model),
        "command": " ".join(cmd),
        "workdir": str(work),
        "returncode": rc,
        "ok": rc == 0 and fitfile.exists(),
        "fit_file": str(fitfile) if fitfile.exists() else None,
        "int_file": str(work / f"{stem}.int") if (work / f"{stem}.int").exists() else None,
        "abs_file": str(work / f"{stem}.abs") if (work / f"{stem}.abs").exists() else None,
        "log_file": str(logfile) if logfile.exists() else None,
        "stats": stats,
    }
    if not rec["ok"]:
        tail = [ln for ln in log.strip().splitlines() if ln.strip()][-3:]
        rec["error"] = " | ".join(tail) or f"crysol 退出码 {rc}"
    return rec


# ---------------------------------------------------------------- datcmp


def run_datcmp(bindir: Path, env: dict, fitfile: Path, test: str) -> dict:
    """datcmp 给单个 .out/.fit/.fir 时，比的就是该文件里的「实验数据 vs 模型拟合」。"""
    cmd = [str(bindir / "datcmp"), f"--test={test}", "--format=csv", str(fitfile)]
    rc, out = run(cmd, None, env)
    rec = {"test": test, "command": " ".join(cmd), "returncode": rc, "raw": out.strip()}
    row = None
    for line in out.splitlines():
        parts = [p.strip() for p in line.split(",")]
        if len(parts) >= 5 and parts[0] and not parts[0].startswith(("Hypothesis", "Alternative")):
            row = parts
    if row:
        try:
            rec["value"] = float(row[-3])
            rec["p"] = float(row[-2])
            rec["p_adj"] = float(row[-1])
        except ValueError:
            pass
    # 哨兵：χ² ≤ 0（datcmp 算不出，如 BIFT 的 .ift）或 χ²=0 & p=1.0（喂了普通 .dat → 自己跟自己比）
    v = rec.get("value")
    rec["usable"] = bool(v is not None and v > 0 and not (v == 0.0 and rec.get("p") == 1.0))
    if not rec["usable"]:
        rec["sentinel"] = ("χ²≤0：这个文件不是 datcmp 认得的拟合格式"
                           if v is not None and v <= 0 else
                           "χ²=0 且 p=1：datcmp 把普通 .dat 当成了「自己 vs 自己」，此处的『完美拟合』无意义")
    return rec


# ---------------------------------------------------------------- 组装


def add_curve(rec: dict, data: np.ndarray, bindir: Path, env: dict, a: argparse.Namespace) -> dict:
    ff = rec.get("fit_file")
    if ff and Path(ff).exists():
        suffix = Path(ff).suffix.lower()
        if suffix in NON_COLUMN_SUFFIXES:
            # IFT 的 .out/.ift 不是四列拟合文件：.out 里是 p(r) 表 + 拟合块，.ift 是 BIFT 自己的格式
            rec["ift_only"] = True
            rec["ranking_excluded"] = True   # 口径不同（含被裁掉的低 q），不参与横向排序
            rec["note"] = ("IFT 文件不是四列拟合（列是 p(r)/多段文本）→ 只能交给 datcmp；"
                           "且 .out 的拟合块含被裁掉的低 q 与 q=0 点，datcmp 的数与 IFT 表里的 χ² 不是一回事")
        else:
            arr = load_sasdata(Path(ff))
            rec["n_points"] = int(arr.shape[0])
            rec["q_min"] = float(arr[:, 0].min())
            rec["q_max"] = float(arr[:, 0].max())
            rec["n_columns"] = int(arr.shape[1])
            rec["chi2_reported"] = fit_header_value(Path(ff))
            if arr.shape[1] >= 4:
                rec.update(chi2_regrid(data, arr))
            else:
                rec["note"] = "只有 3 列（无误差列）→ 算不了 reduced χ²，只能当形状参考"
        if suffix in DATCMP_OK_SUFFIXES and not a.skip_datcmp:
            rec["datcmp"] = {
                t: run_datcmp(bindir, env, Path(ff), t)
                for t in ("chi-square", "cormap")
            }
    return rec


def discover_existing(processed: Path, outdir: Path) -> list[tuple[Path, str]]:
    found: dict[str, str] = {}
    for pat, label in FIT_PATTERNS:
        for p in sorted(processed.glob(pat)):
            if outdir.resolve() in p.resolve().parents:
                continue
            found.setdefault(str(p.resolve()), label)
    by_stem: dict[str, list[Path]] = {}
    for k in found:
        by_stem.setdefault(Path(k).stem, []).append(Path(k))
    keep: list[tuple[Path, str]] = []
    for path, label in [(Path(k), v) for k, v in found.items()]:
        sibs = by_stem.get(path.stem, [])
        if len(sibs) > 1:
            try:
                widths = {p: ncols(p) for p in sibs}
            except SystemExit:
                widths = {}
            if widths and widths.get(path, 99) < 4 and any(w >= 4 for w in widths.values()):
                continue          # 同名 4 列文件（.fir）已经覆盖这个 3 列的绘图文件
        keep.append((path, label))
    return keep


def primary_chi2(rec: dict) -> tuple[float | None, str, int | None]:
    """选出「用来判」的那个 χ²：datcmp（ATSAS 统一口径）优先，其次同网格重算，最后程序自报。"""
    dc = (rec.get("datcmp") or {}).get("chi-square") or {}
    if dc.get("usable"):
        return dc["value"], "datcmp", rec.get("n_points")
    if rec.get("chi2_regrid") is not None and not rec.get("regridded"):
        return rec["chi2_regrid"], "同网格重算", rec.get("n_common")
    if rec.get("chi2_reported") is not None:
        return rec["chi2_reported"], "程序自报", rec.get("n_points")
    return None, "无", None


def verdict(rec: dict, a: argparse.Namespace, n_data: int | None = None) -> tuple[str, str]:
    """🟢/🟡/🔴 + 一句判据。阈值出处写进报告：χ²_n 的 95% 区间（手册说该区间可由 χ²_n 分布算）+ CorMap α。"""
    if not rec.get("ok", True):
        return "🔴", f"crysol 没算出来：{rec.get('error', '见 crysol.stdout.txt')}"
    c, source, n = primary_chi2(rec)
    if c is None:
        return "🟡", rec.get("note") or "没有可用于判定的 χ²（列不全或解析失败）"
    if rec.get("ift_only"):
        return "🟡", (f"IFT 文件：datcmp 给的 reduced χ²={c:.3f} 是**全 q 范围**（含被裁掉的低 q 与 q=0 点）上的残差，"
                      f"与上面各行口径不同 → 只作参考，不参与排序")
    approx = ""
    if n is None:                      # .out 这类文件里没有点数：用实验数据点数近似，并在报告里说明
        n, approx = n_data, "≈实验点数"
    lo, hi = chi2_band(n)
    dc = (rec.get("datcmp") or {}).get("cormap") or {}
    p = dc.get("p") if dc.get("usable") else None
    if lo is None:
        return "🟡", f"{source} reduced χ²={c:.3f}；没算出参考区间（无 scipy / 点数为 0），只报数不下判"
    band = f"χ²_n(n={n}{approx}) 95% 区间 [{lo:.2f}, {hi:.2f}]"
    if lo is not None and lo <= c <= hi and (p is None or p > a.alpha):
        return "🟢", f"{source} reduced χ²={c:.3f} 落在 {band}" + (
            "" if p is None else f"，CorMap p={p:.3g}>{a.alpha}")
    if lo is not None and lo <= c <= hi:
        return "🟡", (f"{source} reduced χ²={c:.3f} 在 {band} 内，但 CorMap p={p:.3g}≤{a.alpha} "
                      f"→ 残差有系统性形状偏差（曲线越长 CorMap 越容易显著；绝对判档依赖误差标定）")
    if hi is not None and c <= 2 * hi:
        return "🟡", f"{source} reduced χ²={c:.3f} 超出 {band} 但不到 2 倍"
    return "🔴", (f"{source} reduced χ²={c:.3f} 远在 {band} 之外 → 残差远大于误差棒（绝对判档；"
                  f"误差被高估/低调都会移动这个数，横向排序不受影响）")


def write_outputs(outdir: Path, a: argparse.Namespace, curves: list[dict], data: np.ndarray,
                  data_used: Path, atsas: Path, cmds: list[str], processed: Path | None) -> None:
    # 把每条参与比较的拟合文件拷到 fits/，一处就能看全（原始文件不动）
    fitdir = outdir / "fits"
    fitdir.mkdir(exist_ok=True)
    taken: set[str] = set()
    for c in curves:
        ff = c.get("fit_file")
        if not ff or not Path(ff).exists():
            continue
        src = Path(ff)
        name = f"{src.stem}{src.suffix}"        # 用原名（纯 ASCII），标签里有冒号/中文不适合当文件名
        i = 2
        while name in taken or ((fitdir / name).exists()
                                and (fitdir / name).resolve() != src.resolve()):
            name = f"{src.stem}_{i}{src.suffix}"
            i += 1
        taken.add(name)
        dst = fitdir / name
        if src.resolve() != dst.resolve():
            shutil.copy2(src, dst)
        c["fit_copy"] = str(dst)
    for c in curves:                      # 把「用来判的那个 χ²」也写进 json，画图/下游直接用
        c["primary_chi2"], c["primary_source"], _ = primary_chi2(c)
    n_ref = int(curves[0].get("n_common") or (data[:, 2] > 0).sum())
    lo, hi = chi2_band(n_ref)
    payload = {
        "generated": __import__("time").strftime("%Y-%m-%d %H:%M:%S"),
        "atsas_bin": str(atsas),
        "data_used": str(data_used),
        "data_original": str(Path(a.data).resolve()),
        "processed": str(processed) if processed else None,
        "qmin": a.qmin, "qmax": a.qmax, "constant": a.constant, "lm": a.lm,
        "data_q_min": float(data[:, 0].min()), "data_q_max": float(data[:, 0].max()),
        "smax": a.smax, "shell": a.shell, "unit": a.unit, "alpha": a.alpha,
        "chi2_band_for_n_data": [lo, hi],
        "chi2_band_for_n": n_ref,
        "commands": cmds,
        "curves": curves,
    }
    (outdir / "comparison.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2))
    cols = ["label", "kind", "n_points", "q_min", "q_max", "chi2_reported",
            "chi2_datcmp", "chi2_regrid", "regridded", "n_common", "cormap_C", "cormap_p",
            "primary_chi2", "primary_source", "model_rg_A", "verdict", "reason", "fit_file"]
    with open(outdir / "comparison.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for c in curves:
            dc = (c.get("datcmp") or {})
            prim, src, _ = primary_chi2(c)
            w.writerow({
                "label": c.get("label"), "kind": c.get("kind"),
                "n_points": c.get("n_points"), "q_min": c.get("q_min"), "q_max": c.get("q_max"),
                "chi2_reported": c.get("chi2_reported"),
                "chi2_datcmp": (dc.get("chi-square") or {}).get("value"),
                "chi2_regrid": c.get("chi2_regrid"), "regridded": c.get("regridded"),
                "n_common": c.get("n_common"),
                "primary_chi2": prim, "primary_source": src,
                "cormap_C": (dc.get("cormap") or {}).get("value"),
                "cormap_p": (dc.get("cormap") or {}).get("p"),
                "model_rg_A": (c.get("stats") or {}).get("model_rg_A"),
                "verdict": c.get("verdict"), "reason": c.get("reason"),
                "fit_file": c.get("fit_file"),
            })


# ---------------------------------------------------------------- 主流程


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="用 CRYSOL(ATSAS) 从 model 算 SAXS 曲线并拟合到 .dat，再与已有拟合同口径比较。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--model", action="append", required=True, metavar="FILE",
                   help="高分辨模型（.pdb/.cif/.ent，可重复给多个；珠模型 .pdb 也能喂）")
    p.add_argument("--data", required=True, metavar="FILE",
                   help="扣减后的实验曲线 .dat（q 单位 Å⁻¹）")
    p.add_argument("--processed", metavar="DIR",
                   help="样品结果目录（processed/<模式>/<样品>/）；给了就自动找里面的已有拟合")
    p.add_argument("--fit", action="append", default=[], metavar="FILE",
                   help="额外指定的已有拟合文件（.out/.fit/.fir），可重复")
    p.add_argument("--out", metavar="DIR", help="产物目录（默认 <processed>/model-fit，否则 ./model-fit）")
    p.add_argument("--qmin", type=float, metavar="Q", help="拟合下界 Å⁻¹；给了就用 datcrop 先裁（None=不裁）")
    p.add_argument("--qmax", type=float, metavar="Q", help="拟合上界 Å⁻¹；给了就用 datcrop 先裁")
    p.add_argument("--constant", dest="constant", action="store_true", default=True,
                   help="允许常数扣减（吸收缓冲液失配；手册的单模型拟合示例就用它）")
    p.add_argument("--no-constant", dest="constant", action="store_false",
                   help="关掉常数扣减（想跟『参数更少』的那次比时用）")
    p.add_argument("--lm", type=int, metavar="N", default=20,
                   help="球谐最大阶数（crysol 默认 20，上限 100；细长物体加大，smax 调大时必须加大）")
    p.add_argument("--smax", type=float, metavar="Q", default=None,
                   help="拟合用的最大 q，Å⁻¹（None=用 crysol 默认 0.5）")
    p.add_argument("--shell", choices=["directional", "water"], default=None,
                   help="水化层模型（None=用 crysol 默认 directional）")
    p.add_argument("--unit", choices=["u", "1", "2", "3", "4"], default=None,
                   help="实验数据的角度单位（1=4πsinθ/λ in Å⁻¹；None=让 ATSAS 猜）")
    p.add_argument("--chain", default=None, metavar="ID", help="只算该链（多链 .cif 用）")
    p.add_argument("--model-id", dest="model_id", default=None, metavar="ID",
                   help="只算该 model 号（NMR 多构象 .pdb 用）")
    p.add_argument("--atsas-dir", default=None, metavar="DIR",
                   help="ATSAS 的 bin 目录（默认依次找 $ATSAS、/Applications/ATSAS-4.1.4-1/bin）")
    p.add_argument("--alpha", type=float, default=0.01,
                   help="显著性水平（datcmp 判「能不能拒绝拟合」用）")
    p.add_argument("--skip-datcmp", action="store_true", help="跳过 datcmp（只算 crysol + 同网格 χ²）")
    p.add_argument("--no-plot", action="store_true", help="跳过出图")
    a = p.parse_args(argv)

    atsas = find_atsas(a.atsas_dir)
    env = atsas_env(atsas)
    data_orig = Path(a.data).expanduser().resolve()
    processed = Path(a.processed).expanduser().resolve() if a.processed else None
    outdir = Path(a.out).expanduser().resolve() if a.out else (
        processed / "model-fit" if processed else Path.cwd() / "model-fit")
    outdir.mkdir(parents=True, exist_ok=True)
    cmds: list[str] = []

    # 1) 裁 q 区间（要给 crysol 一条"可用区"上的曲线，不然低 q artifact 会主导 χ²）
    data_used = outdir / "data_used.dat"
    if data_orig.resolve() != data_used.resolve():
        shutil.copy2(data_orig, data_used)
    if a.qmin is not None or a.qmax is not None:
        cmd = [str(atsas / "datcrop")]
        if a.qmin is not None:
            cmd.append(f"--smin={a.qmin}")
        if a.qmax is not None:
            cmd.append(f"--smax={a.qmax}")
        cmd += [str(data_orig), "-o", str(data_used)]
        rc, out = run(cmd, None, env)
        cmds.append(" ".join(cmd))
        if rc != 0:
            raise SystemExit(f"datcrop 失败（{rc}）：{out.strip()[:400]}")
    data = load_sasdata(data_used)
    print(f"[data] {data_used}  N={data.shape[0]}  q={data[0,0]:.5f}–{data[-1,0]:.5f} Å⁻¹")

    # 2) 每个模型跑一次 crysol（一个模型一次调用 → 产物名不打架、日志无歧义）
    curves: list[dict] = []
    for m in a.model:
        mp = Path(m).expanduser().resolve()
        if not mp.exists():
            raise SystemExit(f"模型不存在：{mp}")
        rec = run_crysol(atsas, env, mp, data_used, outdir, a)
        cmds.append(rec["command"])
        rec = add_curve(rec, data, atsas, env, a)
        rec["verdict"], rec["reason"] = verdict(rec, a, len(data))
        st = rec.get("stats") or {}
        if rec["ok"]:
            prim, src, _ = primary_chi2(rec)
            print(f"[crysol] {rec['label']}: χ²(程序)={st.get('chi2_fit')} "
                  f"χ²({src})={prim if prim is None else round(prim, 3)} → {rec['verdict']}")
        else:
            print(f"[crysol] {rec['label']}: 失败 —— {rec.get('error')}")
        curves.append(rec)

    # 3) 已有的拟合（DENSS/DAMMIF/IFT…）走同一条流水线
    existing: list[tuple[Path, str]] = []
    if processed:
        existing += discover_existing(processed, outdir)
    existing += [(Path(f).expanduser().resolve(), "指定的已有拟合") for f in a.fit]
    seen = {c.get("fit_file") for c in curves}
    for path, label in existing:
        if str(path) in seen or not path.exists():
            continue
        rec = {"label": f"{label}: {path.stem}", "kind": label, "fit_file": str(path)}
        rec = add_curve(rec, data, atsas, env, a)
        rec["verdict"], rec["reason"] = verdict(rec, a, len(data))
        prim, src, _ = primary_chi2(rec)
        print(f"[已有] {rec['label']}: χ²({src})={prim if prim is None else round(prim, 3)} → {rec['verdict']}")
        curves.append(rec)

    if not curves:
        raise SystemExit("一条曲线都没有：检查 --model/--data/--processed")

    write_outputs(outdir, a, curves, data, data_used, atsas, cmds, processed)

    # 4) 出图（崩了不能丢数字：数字先落盘，出图放最后并吞掉异常）
    if not a.no_plot:
        try:
            import subprocess as _sp
            proc = _sp.run([sys.executable, str(Path(__file__).with_name("plot-model-vs-fit.py")),
                            str(outdir)], text=True, capture_output=True, env=env)
            print(f"[plot] {'OK' if proc.returncode == 0 else 'FAILED'} "
                  f"{(proc.stdout + proc.stderr).strip()[-300:]}")
        except Exception as exc:  # noqa: BLE001
            print(f"[plot] 出图失败（数字已落盘）：{exc}", file=sys.stderr)

    # 5) README（分点总结 + 明细表 + 🟢🟡🔴）
    try:
        rc, out = run([sys.executable, str(Path(__file__).with_name("write-model-fit-readme.py")), str(outdir)],
                      None, env)
        print(f"[readme] {'OK' if rc == 0 else 'FAILED'} {out.strip()[-300:]}")
    except Exception as exc:  # noqa: BLE001
        print(f"[readme] 写 README 失败：{exc}", file=sys.stderr)

    ranked = [(c["label"], primary_chi2(c)[0], primary_chi2(c)[1]) for c in curves
              if primary_chi2(c)[0] is not None and not c.get("ranking_excluded")]
    if ranked:
        best = min(ranked, key=lambda t: t[1])
        print(f"[best] 判据 χ² 最小：{best[0]}（{best[2]} {best[1]:.3f}）")
    print(f"[out] {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
