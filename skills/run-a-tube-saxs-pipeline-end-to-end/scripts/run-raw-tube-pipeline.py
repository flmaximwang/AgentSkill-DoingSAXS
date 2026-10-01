#!/usr/bin/env python
"""End-to-end SAXS pipeline for TUBE / static samples (BL19U2 batch mode), driven
entirely by BioXTAS RAW's own algorithms through bioxtasraw.RAWAPI.

    frames -> integrate -> average -> control scaling -> subtract
           -> multi-range Guinier (gated) -> IFT/P(r) (gated) -> MW -> 3D -> report

Every node is written to disk; every profile is written as .dat.

    python run-raw-tube-pipeline.py --sample-dir <dir> --cfg <day>.cfg --out-dir <out> \
        [--steps integrate,average,subtract,guinier,ift,mw,shape,report,workspace]

The sibling skill `run-a-sec-saxs-pipeline-end-to-end` covers SEC elution series;
this one covers static/tube datasets where each sample run sits next to its own
control runs in the same directory.

Nothing here re-implements integration / Guinier / IFT / reconstruction: the script
only sets parameters and writes files.  RAW does the maths.
"""
from __future__ import annotations

import argparse, glob, json, os, re, shutil, sys, time
from collections import Counter

import numpy as np
import bioxtasraw.RAWAPI as raw


# --------------------------------------------------------------------------- helpers
def log(*a):
    print("[%s]" % time.strftime("%H:%M:%S"), *a, flush=True)


def jdump(obj, path):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=str)


def mkdirs(root, *subs):
    for s in subs:
        os.makedirs(os.path.join(root, s), exist_ok=True)


ATSAS_GLOBS = ("/Applications/ATSAS*/bin", "~/ATSAS*/bin", "/usr/local/ATSAS*/bin",
               "/opt/ATSAS*/bin", "~/atsas/*/bin")


def find_atsas_dir():
    """Locate the ATSAS bin directory (GNOM/DAMMIF/DAMAVER/DATCLASS live there).

    RAW's own SASUtils.findATSASDirectory() only globs ``~/ATSAS*`` and otherwise falls back
    to ``/Applications/ATSAS/bin`` -- an install at ``/Applications/ATSAS-4.1.4-1`` is
    invisible to it (returns ''), which is why RAW's GNOM/DAMMIF menus claimed "Can't find
    ATSAS" while ATSAS was in fact installed.  Probe the usual layouts first and only then
    ask RAW.
    """
    for pat in ATSAS_GLOBS:
        for d in sorted(glob.glob(os.path.expanduser(pat)), reverse=True):
            if any(os.path.exists(os.path.join(d, x)) for x in ("gnom", "gnom.exe")):
                return d
    try:
        from bioxtasraw.SASUtils import findATSASDirectory
        return findATSASDirectory() or ""
    except Exception:
        return ""


class Ctx:
    """RAW settings + output layout shared by all nodes."""

    def __init__(self, cfg, out, hdr_format, atsas_dir=None, save_frames=False):
        self.out = os.path.abspath(os.path.expanduser(out))
        self.save_frames = save_frames
        mkdirs(self.out, "frames", "norm", "profiles", "profiles/01_control",
               "profiles/02_sample", "profiles/03_subtracted", "profiles/04_guinier",
               "profiles/06_ift", "tables", "ifts", "models", "reports")
        self.s = raw.load_settings(cfg)
        # A .cfg saved at the beamline usually has ImageHdrFormat=None and
        # EnableNormalization=False, so the per-frame .txt is never parsed and the
        # transmission normalization silently does nothing.  Force both.
        self.s.set("ImageHdrFormat", hdr_format)
        self.s.set("EnableNormalization", True)
        self.s.set("NormalizationList", [["/", "Transmitted_Beam"]])
        atsas_dir = atsas_dir or self.s.get("ATSASDir") or find_atsas_dir()
        if atsas_dir:
            self.s.set("ATSASDir", atsas_dir)
        self.atsas_dir = atsas_dir or ""

        def _have(*names):
            return bool(self.atsas_dir) and all(
                shutil.which(x) or os.path.exists(os.path.join(self.atsas_dir, x))
                for x in names)

        self.has_atsas = _have("gnom", "dammif")        # what the bead-model path needs
        # ATSAS >= 3.1 folded damclust into damaver, so damaver is the clustering entry
        self.has_damaver = _have("damaver")
        self.has_datclass = _have("datclass")

    def set(self, key, value):
        self.s.set(key, value)

    def save(self, profile, name, *sub):
        raw.save_profile(profile, name, os.path.join(self.out, *sub), settings=self.s)


def classify(files, sample_key, control_keys):
    """Sample frames = prefix matching the sample key; everything else = control."""
    key = re.sub(r"[^0-9A-Za-z]", "", sample_key).lower()
    sam, ctl = [], []
    for f in files:
        b = re.sub(r"[^0-9A-Za-z]", "", os.path.basename(f).split("_")[0]).lower()
        if control_keys:
            (ctl if b in control_keys else sam).append(f)
        else:
            (sam if b == key else ctl).append(f)
    return sorted(sam), sorted(ctl)


# ------------------------------------------------------------------- node 1: frames
def node_integrate(ctx, sam_files, ctl_files):
    log("integrate %d sample + %d control frames" % (len(sam_files), len(ctl_files)))
    t = time.time()
    sam, _ = raw.load_and_integrate_images(sam_files, settings=ctx.s)
    ctl, _ = raw.load_and_integrate_images(ctl_files, settings=ctx.s)
    rows = []
    for i, p in enumerate(list(sam) + list(ctl)):
        h = p.getParameter("counters") or {}
        q, I = p.getQ(), p.getI()
        rows.append(dict(file=p.getParameter("filename"),
                         kind="sample" if i < len(sam) else "control",
                         description=h.get("Description", ""),
                         transmission=h.get("Transmitted_Beam", ""),
                         sr_current=h.get("SR_Current", ""),
                         exposure=h.get("Exposure_time_s", ""),
                         i_lowq=float(I[0]),
                         i_midq=float(I[int(np.argmin(abs(q - 0.2)))])))
        if ctx.save_frames:
            kind = "sample" if i < len(sam) else "control"
            ctx.save(p, "%s_%s.dat" % (kind, os.path.splitext(p.getParameter("filename"))[0]),
                     "frames")
    log("  integrated in %.1fs" % (time.time() - t))
    for kind in ("sample", "control"):
        vals = np.array([r["i_lowq"] for r in rows if r["kind"] == kind])
        if vals.size:
            med = float(np.median(vals))
            for r in rows:
                if r["kind"] == kind:
                    r["dev_from_median"] = r["i_lowq"] / med - 1.0
                    r["outlier"] = bool(abs(r["i_lowq"] / med - 1.0) > 0.15)
    with open(os.path.join(ctx.out, "norm", "frame_qc.csv"), "w") as fh:
        fh.write("file,kind,description,transmission,sr_current,exposure,i_lowq,i_midq,"
                 "dev_from_median,outlier\n")
        for r in rows:
            fh.write("%s,%s,%s,%s,%s,%s,%.6g,%.6g,%.4f,%d\n" % (
                r["file"], r["kind"], r["description"], r["transmission"], r["sr_current"],
                r["exposure"], r["i_lowq"], r["i_midq"], r.get("dev_from_median", float("nan")),
                int(r.get("outlier", False))))
    jdump(rows, os.path.join(ctx.out, "norm", "frame_qc.json"))
    bad = [r["file"] for r in rows if r.get("outlier")]
    log("  frames deviating >15%% from their run median: %s" % (bad if bad else "none"))
    return sam, ctl, rows


# ------------------------------------------------- node 2-4: average + subtract
def node_average_subtract(ctx, sam, ctl, scale_window=(0.30, 0.44), do_scale=True):
    sam_avg = raw.average(sam)
    ctl_avg = raw.average(ctl)
    ctx.save(ctl_avg, "control_avg.dat", "profiles", "01_control")
    ctx.save(sam_avg, "sample_avg.dat", "profiles", "02_sample")

    # The control runs are not measured at exactly the same flux/position as the
    # sample runs, so a 1-3 % level difference survives transmission normalization
    # and leaves a constant residual at high q which wrecks the IFT.  Read the
    # factor off the q window where the particle does not scatter and apply it with
    # RAW's own relative-scale operation.
    factor = 1.0
    q = ctl_avg.getQ()
    m = (q >= scale_window[0]) & (q <= scale_window[1])
    if do_scale and m.sum() >= 5:
        factor = float(np.mean(sam_avg.getI()[m] / ctl_avg.getI()[m]))
        if abs(factor - 1.0) > 1e-4:
            ctl_avg = ctl_avg.copy()
            ctl_avg.scaleRelative(factor)
            ctx.save(ctl_avg, "control_avg_scaled.dat", "profiles", "01_control")

    sub = raw.subtract([sam_avg], ctl_avg)[0]
    ctx.save(sub, "subtracted.dat", "profiles", "03_subtracted")
    q, I = sub.getQ(), sub.getI()
    i_sam = sam_avg.getI()
    hi = slice(int(len(q) * 0.6), int(len(q) * 0.8))
    info = dict(n_points=len(q), qmin=float(q[0]), qmax=float(q[-1]),
                control_scale_factor=factor, scale_window=list(scale_window),
                contrast_lowq=float(I[0] / i_sam[0]) if i_sam[0] else None,
                sample_over_control_midq=float(np.mean(i_sam[hi] / ctl_avg.getI()[hi])),
                negative_kratky_frac=float(np.mean((q ** 2 * I)[hi] < 0)))
    log("subtraction: control scale %.4f (q %.2f-%.2f) | contrast %.3f of sample I(q_min)"
        % (factor, scale_window[0], scale_window[1], info["contrast_lowq"]))
    if abs(factor - 1.0) > 0.05:
        log("  WARNING: scale factor %.3f deviates >5%% from 1 - check the control run"
            % factor)
    return sub, info


def pick_analysis_window(sub, snr_min=2.0, qmax_cap=0.35, run=20, qmin=0.0):
    """Largest q ending a contiguous run of points with I>0 and I/err >= snr_min.

    qmin (A^-1) drops the first points: on BL19U2 the beam stop covers r<=12 px and the
    rings up to ~r 20 px are only partly unmasked, so q <~0.008 sits on the beam-stop edge
    and carries a beam-stop halo that no background subtraction removes (measured: cutting
    there removes ~75-80 % of the apparent low-q upturn of the dilute samples).
    """
    q, I, E = sub.getQ(), sub.getI(), sub.getErr()
    snr = np.divide(I, E, out=np.zeros_like(I), where=E > 0)
    good = (snr >= snr_min) & (I > 0)
    start = int(np.argmin(abs(q - qmin))) if qmin > 0 else 0
    cap = int(np.argmin(abs(q - qmax_cap)))
    for i in range(cap, -1, -1):
        if good[max(0, i - run + 1):i + 1].all():
            return (start if i - start >= 10 else 0), i
    return start, start


# --------------------------------------------------- node 5: multi-range Guinier
def _chi2_red(q, I, E, rg, i0, i0idx, i1idx):
    """Error-weighted reduced chi^2 of RAW's own fit line (residual statistic)."""
    m = slice(i0idx, i1idx + 1)
    fit = np.log(i0) - (rg ** 2 / 3.0) * q[m] ** 2
    res = np.log(I[m]) - fit
    sig = np.divide(E[m], I[m], out=np.full_like(res, np.nan), where=I[m] > 0)
    ok = np.isfinite(sig) & (sig > 0)
    if ok.sum() < 5:
        return np.nan, np.nan
    chi2 = float(np.sum((res[ok] / sig[ok]) ** 2) / (ok.sum() - 2))
    half = ok.sum() // 2
    curv = float(np.mean(res[ok][:half]) - np.mean(res[ok][half:]))   # >0 smile
    return chi2, curv


def node_guinier(ctx, sub, qmax_idx, r_gate_qrg, min_pts=10, chi2_max=3.0,
                 qrg_rungs=(0.8, 1.0, 1.2, 1.3), start_fracs=(0.05, 0.10, 0.15, 0.20),
                 qmin_idx=0):
    q, I, E = sub.getQ(), sub.getI(), sub.getErr()
    auto = raw.auto_guinier(sub, settings=ctx.s)
    rg0 = auto[0] if auto[0] > 0 else 30.0

    starts = sorted({max(qmin_idx, s) for s in (
        0,
        int(np.argmin(abs(q - auto[4]))) if auto[4] > 0 else 0,
        *[int(qmax_idx * f) for f in start_fracs])})
    ends = sorted({int(np.argmin(abs(q - min(qrg / rg0, q[qmax_idx])))) for qrg in qrg_rungs})

    rows = []
    for i0 in starts:
        for i1 in ends:
            if i1 <= i0 or i1 - i0 + 1 < min_pts or i0 < qmin_idx:
                continue
            try:
                r = raw.guinier_fit(sub, i0, i1, settings=ctx.s)
            except Exception as e:
                rows.append(dict(idx_min=i0, idx_max=i1, failed=str(e)))
                continue
            rg, i0v, rg_e, i0_e, qmin, qmax, qrg_min, qrg_max, r2 = (
                r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8])
            chi2r, curv = _chi2_red(q, I, E, rg, i0v, i0, i1)
            half = (i0 + i1) // 2
            try:
                rg_split = abs(raw.guinier_fit(sub, i0, half, settings=ctx.s)[0] -
                               raw.guinier_fit(sub, half, i1, settings=ctx.s)[0]) / rg
            except Exception:
                rg_split = np.nan
            rows.append(dict(idx_min=i0, idx_max=i1, n_pts=i1 - i0 + 1, qmin=float(qmin),
                             qmax=float(qmax), rg=float(rg), rg_err=float(rg_e),
                             i0=float(i0v), i0_err=float(i0_e), qrg_min=float(qrg_min),
                             qrg_max=float(qrg_max), r_sqr=float(r2), chi2_red=chi2r,
                             curvature=curv, rg_split_rel=float(rg_split),
                             gates=dict(n_pts=i1 - i0 + 1 >= min_pts,
                                        qrg=(qrg_max <= r_gate_qrg),
                                        chi2=(np.isfinite(chi2r) and chi2r <= chi2_max),
                                        stable=(np.isfinite(rg_split) and rg_split <= 0.05))))

    ok_rows = [r for r in rows if "failed" not in r]
    good = [r for r in ok_rows if all(r["gates"].values())]
    rec = None
    if good:
        med = float(np.median([r["rg"] for r in good]))
        rec = sorted(good, key=lambda r: (-r["n_pts"], abs(r["rg"] - med)))[0]
        rec["reason"] = ("passes every gate; widest passing range with Rg at the median "
                         "of the passing set (Rg %.2f vs median %.2f)" % (rec["rg"], med))
    elif ok_rows:
        rec = sorted(ok_rows, key=lambda r: (r["chi2_red"], -r["n_pts"]))[0]
        fails = [k for k, v in rec["gates"].items() if not v]
        rec["reason"] = ("NOT recommended - no range passes all gates (failed: %s). Best "
                         "chi2_red=%.2f at q %.4f-%.4f, but Rg drifts %.2f-%.2f A across "
                         "ranges => the curve does not follow a single Guinier law "
                         "(polydispersity / residual background / wrong control)."
                         % (",".join(fails), rec["chi2_red"], rec["qmin"], rec["qmax"],
                            min(r["rg"] for r in ok_rows), max(r["rg"] for r in ok_rows)))

    # one profile copy per range, trimmed with RAW's own setQrange
    off = sub.getQrange()[0]
    for r in ok_rows:
        p = sub.copy()
        p.setQrange((off + r["idx_min"], off + r["idx_max"] + 1))
        ctx.save(p, "guinier_%.5f-%.5f.dat" % (r["qmin"], r["qmax"]), "profiles", "04_guinier")
    if rec:
        p = sub.copy()
        p.setQrange((off + rec["idx_min"], off + rec["idx_max"] + 1))
        ctx.save(p, "recommended_range.dat", "profiles", "04_guinier")

    jdump(dict(auto=dict(rg=float(auto[0]), i0=float(auto[1]), rg_err=float(auto[2]),
                         i0_err=float(auto[3]), qmin=float(auto[4]), qmax=float(auto[5]),
                         qrg_max=float(auto[7]), r_sqr=float(auto[10])),
               ranges=rows, recommended=rec, n_ranges=len(rows), n_passing=len(good)),
          os.path.join(ctx.out, "tables", "guinier_results.json"))
    with open(os.path.join(ctx.out, "tables", "guinier_multi_range.csv"), "w") as fh:
        fh.write("idx_min,idx_max,n_pts,qmin,qmax,Rg,Rg_err,I0,qRg_min,qRg_max,R2,chi2_red,"
                 "curvature,rg_split_rel,pass_n,pass_qRg,pass_chi2,pass_stable\n")
        for r in rows:
            if "failed" in r:
                fh.write("%d,%d,-,-,-,-,-,-,-,-,-,-,-,-,,,,\n" % (r["idx_min"], r["idx_max"]))
                continue
            g = r["gates"]
            fh.write("%d,%d,%d,%.5f,%.5f,%.3f,%.3f,%.2f,%.3f,%.3f,%.5f,%.3f,%.4f,%.4f,"
                     "%d,%d,%d,%d\n" % (r["idx_min"], r["idx_max"], r["n_pts"], r["qmin"],
                                        r["qmax"], r["rg"], r["rg_err"], r["i0"], r["qrg_min"],
                                        r["qrg_max"], r["r_sqr"], r["chi2_red"], r["curvature"],
                                        r["rg_split_rel"], g["n_pts"], g["qrg"], g["chi2"],
                                        g["stable"]))
    log("Guinier: auto Rg=%.2f (q %.4f-%.4f, R2=%.3f) | %d ranges, %d pass gates"
        % (auto[0], auto[4], auto[5], auto[10], len(rows), len(good)))
    log("  %s" % (rec["reason"] if rec else "no usable range"))
    return auto, rows, rec


# --------------------------------------------------------------- node 6: IFT / P(r)
def node_ift(ctx, sub, rec, i0, i1, dmax_scale=3.1, dmax_pts=10, engines=("bift", "gnom"),
             gnom_dmax=None):
    """Run RAW's BIFT twice - once from the start of the analysis window, once from the
    Guinier fit's own start index - then keep the trusted one.

    Why both: BIFT's minDmax/maxDmax only bound the *grid search*; the joint optimiser may
    return a much larger Dmax (RAWAPI.bift docstring: "The value of Dmax can go beyond this
    bound in the optimization step").  Feeding the contaminated low-q points in is what makes
    it run away (measured: A5-05-6 gave Dmax=741 A from a 48-109 A window, chisq 1.09 - a
    perfect-looking fit to a meaningless P(r)).  RAW's own hint is `use_guinier_start`, i.e.
    start the IFT where the Guinier fit says the data become usable.
    """
    rg = rec["rg"] if rec else raw.auto_guinier(sub, settings=ctx.s)[0]
    lo, hi = max(10.0, 0.7 * dmax_scale * rg), max(30.0, 1.6 * dmax_scale * rg)
    ctx.set("minDmax", lo)
    ctx.set("maxDmax", hi)
    ctx.set("DmaxPoints", dmax_pts)
    ctx.set("PrPoints", 100)

    starts = [(i0, "window-start")]
    if rec and rec["idx_min"] > i0:
        starts.append((rec["idx_min"], "guinier-fit-start"))

    rows, ifts = [], {}
    for idx_min, tag in starts:
        log("IFT [%s]: q idx %d-%d, Dmax grid %.0f-%.0f A" % (tag, idx_min, i1, lo, hi))
        t = time.time()
        res = raw.bift(sub, idx_min=idx_min, idx_max=i1, settings=ctx.s, single_proc=True)
        ift = res[0]
        if ift is None or not (res[1] > 0):
            rows.append(dict(tag=tag, idx_min=idx_min, failed="BIFT returned no solution"))
            log("  failed")
            continue
        dmax, rgr, i0r, dmax_e, rgr_e = res[1], res[2], res[3], res[4], res[5]
        chisq = res[7]
        # The Rg the IFT is compared against is only a real reference when the multi-range
        # Guinier scan actually passed its gates; if that reference itself failed, the
        # agreement test gets a soft (25 %) tolerance instead of 10 % and is flagged.
        soft = not (rec and all(rec["gates"].values()))
        rg_tol = 0.25 if soft else 0.10
        gates = dict(rg_vs_guinier=abs(rgr - rg) / rg <= rg_tol,
                     dmax_over_rg=(dmax / rg) <= 4.5,
                     dmax_within_grid=(dmax <= 1.5 * hi))
        row = dict(tag=tag, idx_min=idx_min, qmin=float(sub.getQ()[idx_min]), qmax=float(sub.getQ()[i1]),
                   dmax=float(dmax), dmax_err=float(dmax_e), rg_realspace=float(rgr),
                   rg_err=float(rgr_e), i0=float(i0r), chisq=float(chisq),
                   rg_tol=rg_tol, rg_ref_soft=bool(soft),
                   seconds=time.time() - t, gates=gates, trusted=all(gates.values()))
        rows.append(row)
        ifts[tag] = ift
        log("  Dmax=%.1f+-%.1f  Rg(real)=%.1f  chisq=%.2f  trusted=%s  %s  rg_tol=%.0f%%%s"
            "  (%.0fs)" % (dmax, dmax_e, rgr, chisq, row["trusted"], gates,
                          rg_tol * 100, " [soft ref]" if soft else "", row["seconds"]))

    def _rank(r):                      # fewest failed gates first, then chi^2
        return (sum(1 for v in r["gates"].values() if not v), r["chisq"])

    solved = [r for r in rows if "failed" not in r]
    trusted = [r for r in solved if r["trusted"]]
    best = min(trusted, key=lambda r: r["chisq"]) if trusted else (
        min(solved, key=_rank) if solved else None)
    out = dict(runs=rows, chosen=best["tag"] if best else None,
               trusted=bool(best and best.get("trusted")))
    if not best:
        log("IFT: no usable solution - 3D will be skipped")
    elif not best.get("trusted"):
        log("IFT: no *trusted* solution (see gates); best is %s - 3D will be skipped and the "
            "P(r) must not be reported" % best["tag"])
    ift = ifts.get(best["tag"]) if best else None
    # ---------------------------------------------------------------- GNOM (needs ATSAS)
    # DAMMIF only accepts a GNOM IFT, so when ATSAS is present GNOM is what the bead model
    # eats.  It must run on the *clean* window: the points below the Guinier start carry the
    # beam-stop halo, and DAMMIF's chi^2 (computed against the data on the whole q range it
    # was given) blows up on them even though GNOM's own chi^2 still looks fine.  Measured on
    # A5-05-1 with identical settings: window 0.0064-0.259 -> DAMMIF chi^2 556 (Rg 20.5,
    # Dmax 72), window 0.0567-0.259 -> chi^2 1.4 (Rg 20.8, Dmax 81).
    gnom_ift, gnom_row = None, None
    if "gnom" in engines:
        if not ctx.has_atsas:
            log("IFT [gnom]: skipped (会用到 GNOM/DAMMIF 时才需要 ATSAS，未找到)")
        else:
            try:
                rg_ref = float(rec["rg"]) if rec else float(raw.auto_guinier(sub, settings=ctx.s)[0])
                gi0 = int(rec["idx_min"]) if (rec and rec.get("idx_min") is not None
                                              and rec["idx_min"] > i0) else i0
                dmax_g = float(gnom_dmax) if gnom_dmax else (
                    float(best["dmax"]) if (best and best.get("trusted"))
                    else round(dmax_scale * rg_ref, 1))
                log("IFT [gnom]: q idx %d-%d (q %.4f-%.4f), Dmax=%.1f A (RG ref %.2f)"
                    % (gi0, i1, sub.getQ()[gi0], sub.getQ()[i1], dmax_g, rg_ref))
                t = time.time()
                # rg / idx_min / idx_max are all passed explicitly: with cut_dam=True and
                # rg=None the wrapper reads analysis_dict['guinier'], which does not exist
                # on a profile loaded from a .dat -> KeyError (RAWAPI.py:3146).
                g = raw.gnom(sub, dmax=dmax_g, rg=rg_ref, idx_min=gi0, idx_max=i1,
                             atsas_dir=ctx.atsas_dir, settings=ctx.s, cut_dam=False,
                             save_ift=True, datadir=os.path.join(ctx.out, "ifts"),
                             savename="gnom.out")
                gnom_ift = g[0]
                dmax, rgr, i0r = (float(g[1]), float(g[2]), float(g[3]))
                rg_e, i0_e = float(g[4]), float(g[5])
                total_est, chisq, alpha = float(g[6]), float(g[7]), float(g[8])
                quality = str(g[9]) if len(g) > 9 else ""
                gates = dict(rg_vs_guinier=abs(rgr - rg_ref) / rg_ref <= (0.25 if soft else 0.10),
                             dmax_over_rg=(dmax / rg_ref) <= 4.5,
                             dmax_within_grid=True)      # GNOM fits Dmax directly, no grid
                gnom_row = dict(tag="gnom", engine="GNOM", idx_min=gi0,
                                qmin=float(sub.getQ()[gi0]), qmax=float(sub.getQ()[i1]),
                                dmax=dmax, dmax_err=float("nan"), rg_realspace=rgr, rg_err=rg_e,
                                i0=i0r, chisq=chisq, rg_tol=0.25 if soft else 0.10,
                                rg_ref_soft=bool(soft), total_est=total_est, alpha=alpha,
                                quality=quality, seconds=time.time() - t,
                                gates=gates, trusted=all(gates.values()))
                rows.append(gnom_row)
                ift_summary_note = np.column_stack([gnom_ift.r, gnom_ift.p, gnom_ift.err])
                np.savetxt(os.path.join(ctx.out, "ifts", "pr_gnom.dat"), ift_summary_note,
                           header="r(A)\tP(r)\terr\t[GNOM IFT, q %.4f-%.4f, Dmax %.1f]"
                                  % (sub.getQ()[gi0], sub.getQ()[i1], dmax_g))
                np.savetxt(os.path.join(ctx.out, "ifts", "ift_fit_gnom.dat"),
                           np.column_stack([gnom_ift.q_orig, gnom_ift.i_orig,
                                            gnom_ift.err_orig, gnom_ift.i_fit]),
                           header="q\tI_measured\terr\tI_fit(GNOM)")
                log("  Dmax=%.1f  Rg(real)=%.1f  chisq=%.3f  total_est=%.3f  alpha=%.2f  %s"
                    "  trusted=%s  (%.0fs)"
                    % (dmax, rgr, chisq, total_est, alpha, quality, gnom_row["trusted"],
                       gnom_row["seconds"]))
                if gnom_row["trusted"] and (best is None or not best.get("trusted")
                                            or chisq < best["chisq"]):
                    best, ift = gnom_row, gnom_ift
                    out["chosen"], out["trusted"] = "gnom", True
                    log("  GNOM 取代 BIFT 成为本次报告的 P(r)（闸门更干净、chisq 更低）")
            except Exception as e:
                log("IFT [gnom] failed: %s: %s" % (type(e).__name__, e))
    out["runs"] = rows
    out["gnom"] = gnom_row
    if gnom_row is not None and gnom_row["trusted"] and out.get("chosen") == "gnom":
        ift = gnom_ift                       # GNOM 赢：报告的 P(r) 换成 GNOM 的
    elif gnom_ift is not None and not (gnom_row and gnom_row["trusted"]):
        raw.save_ift(gnom_ift, "gnom_untrusted.ift", os.path.join(ctx.out, "ifts"))

    # --------------------------- 报告的 P(r)：哪个引擎赢就报哪个（GNOM 优先于 BIFT）----
    if ift is not None and best is not None and best.get("trusted"):
        engine = best.get("engine", "BIFT")
        raw.save_ift(ift, "gnom.out" if engine == "GNOM" else "bift.ift",
                     os.path.join(ctx.out, "ifts"))
        np.savetxt(os.path.join(ctx.out, "ifts", "pr.dat"),
                   np.column_stack([ift.r, ift.p, ift.err]),
                   header="r(A)\tP(r)\terr\t[%s, start=%s - 这就是本次报告采用的 P(r)]"
                          % (engine, best["tag"]))
        np.savetxt(os.path.join(ctx.out, "ifts", "ift_fit.dat"),
                   np.column_stack([ift.q_orig, ift.i_orig, ift.err_orig, ift.i_fit]),
                   header="q\tI_measured\terr\tI_fit(%s)" % engine)
        log("IFT: reported P(r) = %s (%s), Dmax=%.1f Rg=%.1f chisq=%.3f"
            % (engine, best["tag"], best["dmax"], best["rg_realspace"], best["chisq"]))
    elif ift is not None:
        # nothing passed: keep the best-effort solution under a name that says so, so it can
        # be inspected in the RAW GUI without being mistaken for a result
        tag = best["tag"]
        raw.save_ift(ift, "%s_untrusted.ift" % tag, os.path.join(ctx.out, "ifts"))
        np.savetxt(os.path.join(ctx.out, "ifts", "pr_untrusted_%s.dat" % tag),
                   np.column_stack([ift.r, ift.p, ift.err]),
                   header="r(A)\tP(r)\terr\t[%s, start=%s - UNTRUSTED, do not report]"
                          % (best.get("engine", "BIFT"), tag))

    jdump(out, os.path.join(ctx.out, "tables", "ift_summary.json"))
    with open(os.path.join(ctx.out, "tables", "ift_summary.csv"), "w") as fh:
        fh.write("tag,engine,idx_min,qmin,qmax,Dmax,Dmax_err,Rg_realspace,Rg_err,I0,chisq,"
                 "rg_tol,rg_ref_soft,pass_rg,pass_dmax_over_rg,pass_dmax_within_grid,"
                 "trusted,seconds\n")
        for r in rows:
            if "failed" in r:
                fh.write("%s,%s,%d,-,-,-,-,-,-,-,-,,-,,,,0\n"
                         % (r["tag"], r.get("engine", "BIFT"), r["idx_min"]))
                continue
            g = r["gates"]
            fh.write("%s,%s,%d,%.5f,%.5f,%.2f,%.2f,%.2f,%.2f,%.2f,%.4f,%.2f,%d,%d,%d,%d,%d,"
                     "%.0f\n"
                     % (r["tag"], r.get("engine", "BIFT"), r["idx_min"], r["qmin"], r["qmax"],
                        r["dmax"], r["dmax_err"], r["rg_realspace"], r["rg_err"], r["i0"],
                        r["chisq"], r["rg_tol"], int(r["rg_ref_soft"]), g["rg_vs_guinier"],
                        g["dmax_over_rg"], g["dmax_within_grid"], r["trusted"], r["seconds"]))
    return (ift if (best and best.get("trusted")) else None), out, gnom_ift


def node_ift_sweep(ctx, sub, i0, i1, rg, n=5, span=0.35):
    center = 3.1 * rg
    rows = []
    for f in np.linspace(1.0 - span, 1.0 + span, n):
        d = center * f
        ctx.set("minDmax", max(10.0, d - 5.0))
        ctx.set("maxDmax", d + 5.0)
        try:
            r = raw.bift(sub, idx_min=i0, idx_max=i1, settings=ctx.s, single_proc=True)
            rows.append(dict(dmax_target=float(d), dmax=float(r[1]), rg=float(r[2]),
                             chisq=float(r[7])))
        except Exception as e:
            rows.append(dict(dmax_target=float(d), failed=str(e)))
        log("  Dmax sweep %.0f A -> %s" % (d, rows[-1].get("chisq", rows[-1].get("failed"))))
    with open(os.path.join(ctx.out, "tables", "ift_dmax_sweep.csv"), "w") as fh:
        fh.write("dmax_target,dmax,Rg_realspace,chisq\n")
        for r in rows:
            if "failed" in r:
                fh.write("%.1f,-,-,- (%s)\n" % (r["dmax_target"], r["failed"]))
            else:
                fh.write("%.1f,%.1f,%.2f,%.4f\n" % (r["dmax_target"], r["dmax"], r["rg"],
                                                    r["chisq"]))
    return rows


# ---------------------------------------------------------------------- node 7: MW
def node_mw(ctx, sub, vp_mw=None):
    mw = {}
    for name, fn in (("Vp_porod", raw.mw_vp), ("Vc", raw.mw_vc)):
        try:
            v = fn(sub, settings=ctx.s)
            mw[name] = dict(mw=float(v[0]), aux=[float(x) for x in v[1:]])
        except Exception as e:
            mw[name] = dict(failed=str(e))
    if ctx.has_atsas:
        for name, fn in (("datmw_bayes", raw.mw_bayes), ("datclass", raw.mw_datclass)):
            try:
                v = fn(sub, atsas_dir=ctx.atsas_dir)
                mw[name] = dict(mw=float(v[0]), aux=[float(x) for x in v[1:]])
            except Exception as e:
                mw[name] = dict(failed=str(e))
    else:
        mw["note"] = ("datmw/datclass/GNOM/DAMMIF/DAMAVER 需要 ATSAS（没找到）——"
                      "本次只有 RAW 原生的 Vp/Vc/BIFT/DENSS")
    with open(os.path.join(ctx.out, "tables", "mw.csv"), "w") as fh:
        fh.write("method,MW_kDa,aux,note\n")
        for k, v in mw.items():
            if isinstance(v, dict):
                fh.write("%s,%s,%s,%s\n" % (k, v.get("mw", ""), v.get("aux", ""),
                                            v.get("failed", "")))
            else:
                fh.write("%s,,,,%s\n" % (k, v))
    log("MW: %s" % json.dumps({k: (v.get("mw") if isinstance(v, dict) else v)
                               for k, v in mw.items()}))
    return mw


# -------------------------------------------------------------------- node 8: 3D
def node_shape(ctx, ift, engine="auto", mode="Fast", n_electrons=None, symmetry=0,
               n_models=4, dammif_mode="Fast", model_format="cif", dammif_symmetry="P1",
               gnom_ift=None):
    """形状重建。 电子云（RAW 原生 DENSS）与珠模（ATSAS DAMMIF）是**同一份 IFT 的两种重建**，
    互不替代：DENSS 给电子密度图（.mrc），DAMMIF 给 dummy-atom 珠模（.cif/.pdb）。
    默认两个都出；没有 ATSAS 时退化成"只有电子云"，并在结果里写明原因。

    珠模只能用 GNOM 的 IFTM（DAMMIF 吃 GNOM 的 .out）——所以这里要 `gnom_ift`，不是 BIFT 那个。
    DAMAVER 只对"同一份数据的多个 DAMMIF 模型"做一致性平均（mean NSD + cluster），
    也是 ATSAS ≥3.1 里 damclust 的替代入口。
    """
    out = {"engine": engine, "n_models": int(n_models)}
    if engine == "none":
        log("3D: skipped (engine=none)")
        return out
    do_denss = engine in ("auto", "denss", "both")
    do_dammif = engine in ("auto", "dammif", "both")
    mdir = os.path.join(ctx.out, "models")
    os.makedirs(mdir, exist_ok=True)

    # ---------------------------------------------------------------- 电子云（DENSS）
    if do_denss:
        if ift is None:
            out["denss"] = "skipped: no IFT"
            log("3D DENSS: skipped (no IFT)")
        else:
            t = time.time()
            try:
                res = raw.denss(ift, "denss", mdir, mode=mode,
                                n_electrons=n_electrons or 10000, symmetry=symmetry,
                                settings=ctx.s)
                out["denss"] = dict(mode=mode, chi2=float(res[1]), rg_model=float(res[2]),
                                    support_volume=float(res[3]), side=int(res[4]),
                                    n_electrons=int(n_electrons or 10000),
                                    seconds=time.time() - t,
                                    map=os.path.join(mdir, "denss.mrc"))
                log("3D DENSS (%s): chi2=%.3f Rg_model=%.1f A support=%.0f A^3 (%.0fs)"
                    % (mode, res[1], res[2], res[3], time.time() - t))
            except Exception as e:
                out["denss"] = "failed: %s" % e
                log("3D DENSS failed: %s" % e)

    # ---------------------------------------------------------------- 珠模（DAMMIF）
    if do_dammif:
        if not ctx.has_atsas:
            out["dammif"] = "skipped: ATSAS not found"
            log("珠模: skipped（没找到 ATSAS，只有电子云）")
        elif gnom_ift is None:
            out["dammif"] = "skipped: no GNOM IFT (DAMMIF only eats GNOM .out)"
            log("珠模: skipped（没有 GNOM 的 IFTM，DAMMIF 吃不了 BIFT）")
        else:
            rows, files = [], []
            for i in range(1, int(n_models) + 1):
                prefix = "dammif_%02d" % i
                t = time.time()
                try:
                    r = raw.dammif(gnom_ift, prefix, mdir, mode=dammif_mode,
                                   symmetry=dammif_symmetry, atsas_dir=ctx.atsas_dir,
                                   model_format=model_format, settings=ctx.s)
                    row = dict(model=prefix, chisq=float(r[0]), rg=float(r[1]),
                               dmax=float(r[2]), mw=float(r[3]),
                               excluded_volume=float(r[4]), seconds=time.time() - t,
                               file="%s-1.%s" % (prefix, model_format))
                    rows.append(row)
                    files.append(row["file"])
                    log("  珠模 DAMMIF #%d: chi2=%.1f Rg=%.1f Dmax=%.1f MW=%.0f kDa (%.0fs)"
                        % (i, row["chisq"], row["rg"], row["dmax"], row["mw"],
                           row["seconds"]))
                except Exception as e:
                    rows.append(dict(model=prefix, failed="%s: %s" % (type(e).__name__, e)))
                    log("  珠模 DAMMIF #%d 失败: %s: %s" % (i, type(e).__name__, e))
                    break
            out["dammif"] = rows
            # DAMAVER：模型间一致性（mean NSD）+ 聚类（ATSAS>=3.1 里替代 damclust）
            if len(files) >= 2:
                if not ctx.has_damaver:
                    out["damaver"] = "skipped: damaver not found"
                else:
                    try:
                        a = raw.damaver(files, "damaver", mdir,
                                        symmetry=dammif_symmetry,
                                        model_format=model_format,
                                        atsas_dir=ctx.atsas_dir)
                        clusters = [dict(num=int(c.num), rep=str(c.rep_model))
                                    for c in (a[7] or [])]
                        res_val = a[4] if len(a) > 4 else None
                        out["damaver"] = dict(
                            mean_nsd=float(a[0]), stdev_nsd=float(a[1]),
                            rep_model=str(a[2]), n_models=len(files),
                            resolution=(float(res_val) if isinstance(res_val, (int, float))
                                        and float(res_val) > 0 else None),
                            res_unit=(str(a[6]) if len(a) > 6 else ""),
                            clusters=clusters,
                            files=[str(k) for k in (a[3] or {})])
                        log("  DAMAVER: 平均 NSD=%.3f ± %.3f，代表模型 %s，%d 个 cluster"
                            % (a[0], a[1], a[2], len(clusters)))
                    except Exception as e:
                        out["damaver"] = "failed: %s: %s" % (type(e).__name__, e)
                        log("  DAMAVER 失败: %s: %s" % (type(e).__name__, e))
            else:
                out["damaver"] = "skipped: 需要 ≥2 个珠模（这次只有 %d 个）" % len(files)

    jdump(out, os.path.join(ctx.out, "tables", "shape_results.json"))
    return out


# ------------------------------------------------------------------ node 9: report
def node_report(ctx, sub, ift, denss_paths=None):
    try:
        raw.save_report("raw_report.pdf", os.path.join(ctx.out, "reports"),
                        profiles=[sub], ifts=[ift] if ift else [],
                        denss_data=denss_paths or [])
        log("report -> reports/raw_report.pdf")
    except Exception as e:
        log("report failed: %s" % e)


def node_workspace(ctx, sub, ift, key):
    """One .hdf5 workspace so the whole chain can be inspected in the RAW GUI."""
    try:
        raw.save_workspace("%s_workspace" % key, ctx.out, profiles=[sub],
                           ifts=[ift] if ift else [])
        p = os.path.join(ctx.out, "%s_workspace.hdf5" % key)
        log("workspace -> %s (%d bytes)" % (p, os.path.getsize(p)))
    except Exception as e:
        log("workspace failed: %s" % e)


# ------------------------------------------------------------------------- plots
def plots(ctx, sub, rows, rec, ift, sweep=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    q, I, E = sub.getQ(), sub.getI(), sub.getErr()

    def _m(i, e, positive=False):
        """matplotlib refuses negative yerr and can't log-plot non-positive I."""
        i, e = np.asarray(i, float), np.asarray(e, float)
        ok = np.isfinite(i) & np.isfinite(e) & (e >= 0)
        return ok & (i > 0) if positive else ok

    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(2, 3, 1)
    m = _m(I, E, positive=True)
    ax.errorbar(q[m], I[m], yerr=E[m], fmt=".-", ms=2)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q (1/A)"); ax.set_ylabel("I")
    ax.set_title("subtracted (log-log)")
    ax = fig.add_subplot(2, 3, 2)
    m = _m(I, E)
    ax.errorbar(q[m], (q ** 2 * I)[m], yerr=(q ** 2 * E)[m], fmt=".-", ms=2)
    ax.set_xlabel("q"); ax.set_ylabel("q^2 I"); ax.set_title("Kratky")
    ax = fig.add_subplot(2, 3, 3)
    for r in rows:
        if "failed" in r:
            continue
        m = (q >= r["qmin"]) & (q <= r["qmax"]) & _m(I, E, positive=True)
        ax.errorbar(q[m] ** 2, np.log(I[m]), fmt=".", ms=2, alpha=0.25)
    if rec:
        m = (q >= rec["qmin"]) & (q <= rec["qmax"]) & _m(I, E, positive=True)
        ax.errorbar(q[m] ** 2, np.log(I[m]), yerr=(E / I)[m], fmt="o", ms=3,
                    label="recommended")
        xx = np.linspace(0, rec["qmax"] ** 2, 10)
        ax.plot(xx, np.log(rec["i0"]) - xx * rec["rg"] ** 2 / 3, "r-")
        ax.legend(fontsize=8)
    ax.set_xlabel("q^2"); ax.set_ylabel("ln I"); ax.set_title("Guinier fan (all ranges)")
    if ift is not None:
        ax = fig.add_subplot(2, 3, 4)
        m = _m(ift.p, ift.err)
        ax.errorbar(np.asarray(ift.r)[m], np.asarray(ift.p)[m], yerr=np.asarray(ift.err)[m],
                    fmt=".-", ms=2)
        ax.set_xlabel("r (A)"); ax.set_ylabel("P(r)"); ax.set_title("P(r) [RAW BIFT]")
        ax = fig.add_subplot(2, 3, 5)
        m = _m(ift.i_orig, ift.err_orig, positive=True)
        ax.errorbar(np.asarray(ift.q_orig)[m], np.asarray(ift.i_orig)[m],
                    yerr=np.asarray(ift.err_orig)[m], fmt=".", ms=2, label="data")
        mf = np.isfinite(np.asarray(ift.i_fit, float)) & (np.asarray(ift.i_fit, float) > 0)
        ax.plot(np.asarray(ift.q_orig)[mf], np.asarray(ift.i_fit)[mf], "-", label="IFT fit")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q"); ax.set_ylabel("I")
        ax.legend(fontsize=8); ax.set_title("IFT fit")
    if sweep:
        ax = fig.add_subplot(2, 3, 6)
        ax.plot([r["dmax"] for r in sweep if "failed" not in r],
                [r["chisq"] for r in sweep if "failed" not in r], "o-")
        ax.set_yscale("log"); ax.set_xlabel("Dmax (A)"); ax.set_ylabel("chi^2")
        ax.set_title("IFT Dmax sweep")
    plt.tight_layout()
    plt.savefig(os.path.join(ctx.out, "qc.png"), dpi=110)
    log("plots -> qc.png")


# -------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--sample-dir", required=True,
                    help="directory holding the sample frames and the control frames")
    ap.add_argument("--cfg", required=True, help="RAW .cfg saved for this dataset")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--sample-key", default=None,
                    help="file-name prefix of the sample (default: directory name)")
    ap.add_argument("--control-key", default=None,
                    help="comma separated prefixes treated as control (default: everything "
                         "that is not the sample prefix)")
    ap.add_argument("--hdr-format", default="BL19U2, SSRF",
                    help="RAW image-header reader (per-frame .txt must sit next to the tif)")
    ap.add_argument("--scale-window", type=float, nargs=2, default=[0.30, 0.44],
                    metavar=("QMIN", "QMAX"), help="q window for the control scale factor")
    ap.add_argument("--no-scale", action="store_true", help="skip the control scale factor")
    ap.add_argument("--qrg-max", type=float, default=1.3, help="Guinier gate on q*Rg")
    ap.add_argument("--snr-min", type=float, default=2.0, help="I/err cut for the q window")
    ap.add_argument("--qmin", type=float, default=0.0,
                    help="lowest q used (A^-1); >0 drops the beam-stop-edge points — on this "
                         "beamline q<~0.009 lies on the beam stop halo (see SKILL.md)")
    ap.add_argument("--ift-sweep", type=int, default=0, help="sweep Dmax over N values")
    ap.add_argument("--ift-engine", default="both", choices=["bift", "gnom", "both"],
                    help="IFT back-ends.  gnom needs ATSAS and is silently skipped without "
                         "it; it is also the only IFT DAMMIF accepts, so the bead model "
                         "needs it")
    ap.add_argument("--gnom-dmax", type=float, default=None,
                    help="Dmax handed to GNOM (default: the trusted BIFT Dmax, else 3.1*Rg)")
    ap.add_argument("--model-engine", default="auto",
                    choices=["auto", "none", "denss", "dammif", "both"],
                    help="形状重建后端：denss=电子云(.mrc)，dammif=珠模(.cif/.pdb)+DAMAVER，"
                         "auto=有 ATSAS 时两个都出、没有就只出电子云")
    ap.add_argument("--denss-mode", default="Fast", choices=["Fast", "Slow", "Custom"])
    ap.add_argument("--symmetry", type=int, default=0, help="n-fold symmetry for DENSS")
    ap.add_argument("--n-models", type=int, default=4,
                    help="珠模个数（DAMMIF 每次随机起跑；≥2 个才能做 DAMAVER 一致性）")
    ap.add_argument("--dammif-mode", default="Fast", choices=["Fast", "Slow", "Custom"],
                    help="DAMMIF 模式（Fast 在本机约 7 s/个）")
    ap.add_argument("--model-format", default="cif", choices=["cif", "pdb"],
                    help="珠模文件格式。注意 DAMMIF 的 Fast/Slow 模式只写 cif（RAW 只在 "
                         "Custom 模式传 --model-format），所以 pdb 实际拿不到文件")
    ap.add_argument("--dammif-symmetry", default="P1", help="DAMMIF/DAMAVER 对称性（P1 最保守）")
    ap.add_argument("--atsas-dir", default=None,
                    help="ATSAS bin 目录；默认自动探测（/Applications/ATSAS*/bin、~/ATSAS*/bin…）")
    ap.add_argument("--save-frames", action="store_true", help="write every frame as .dat")
    ap.add_argument("--steps", default="ift,mw,shape,report,workspace",
                    help="optional nodes to run; the backbone (integrate, average, "
                         "control scaling, subtract, multi-range Guinier) always runs")
    args = ap.parse_args()
    steps = {s.strip() for s in args.steps.split(",") if s.strip()}
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # for readme_for_run

    sample_dir = os.path.abspath(os.path.expanduser(args.sample_dir))
    key = args.sample_key or os.path.basename(sample_dir.rstrip("/"))
    files = sorted(glob.glob(os.path.join(sample_dir, "*.tif")))
    control_keys = None
    if args.control_key:
        control_keys = {re.sub(r"[^0-9A-Za-z]", "", c).lower()
                        for c in args.control_key.split(",") if c.strip()}
    sam_files, ctl_files = classify(files, key, control_keys)
    if not sam_files or not ctl_files:
        sys.exit("need both sample and control frames in %s" % sample_dir)

    ctx = Ctx(args.cfg, args.out_dir, args.hdr_format, args.atsas_dir, args.save_frames)
    log("sample '%s': %d frames | controls: %d from %s | ATSAS: %s"
        % (key, len(sam_files), len(ctl_files),
           dict(Counter(os.path.basename(f).split("_")[0] for f in ctl_files)),
           (ctx.atsas_dir or "not found")))

    sam, ctl, frames = node_integrate(ctx, sam_files, ctl_files)
    sub, sub_info = node_average_subtract(ctx, sam, ctl, tuple(args.scale_window),
                                          do_scale=not args.no_scale)
    i0, i1 = pick_analysis_window(sub, snr_min=args.snr_min, qmin=args.qmin)
    qmin_idx = int(np.argmin(abs(sub.getQ() - args.qmin))) if args.qmin > 0 else 0
    log("analysis q window: idx %d-%d (q %.4f-%.4f)"
        % (i0, i1, sub.getQ()[i0], sub.getQ()[i1]))
    auto, rows, rec = node_guinier(ctx, sub, i1, args.qrg_max, qmin_idx=qmin_idx)
    ift, ift_res, sweep, denss_res, gnom_ift = None, None, None, None, None
    if "ift" in steps:
        engines = ("bift",) if args.ift_engine == "bift" else (
            ("gnom",) if args.ift_engine == "gnom" else ("bift", "gnom"))
        ift, ift_res, gnom_ift = node_ift(ctx, sub, rec, i0, i1, engines=engines,
                                          gnom_dmax=args.gnom_dmax)
        if args.ift_sweep > 0:
            sweep = node_ift_sweep(ctx, sub, i0, i1, (rec["rg"] if rec else float(auto[0])),
                                   n=args.ift_sweep)
    mw = node_mw(ctx, sub) if "mw" in steps else {}
    if "shape" in steps and args.model_engine != "none":
        ne = None
        vp = (mw.get("Vp_porod") or {}).get("mw") if isinstance(mw.get("Vp_porod"), dict) else None
        if vp:
            ne = int(round(vp * 537))       # ~0.537 electrons per Da (protein)
            log("  DENSS n_electrons from Vp MW %.1f kDa -> %d" % (vp, ne))
        denss_res = node_shape(ctx, ift, engine=args.model_engine, mode=args.denss_mode,
                               n_electrons=ne, symmetry=args.symmetry,
                               n_models=args.n_models, dammif_mode=args.dammif_mode,
                               model_format=args.model_format,
                               dammif_symmetry=args.dammif_symmetry, gnom_ift=gnom_ift)
    if "report" in steps:
        node_report(ctx, sub, ift)
    if "workspace" in steps:
        node_workspace(ctx, sub, ift, key)
    ctl_runs = sorted({"_".join(os.path.basename(f).split("_")[:2]) for f in ctl_files})
    jdump(dict(sample=key, sample_dir=sample_dir, cfg=args.cfg, command=" ".join(sys.argv),
               atsas_dir=ctx.atsas_dir, has_atsas=bool(ctx.has_atsas),
               control_runs={r: 0 for r in ctl_runs}, subtraction=sub_info,
               guinier_auto=dict(rg=float(auto[0]), qmin=float(auto[4]), qmax=float(auto[5]),
                                 r_sqr=float(auto[10])),
               guinier_recommended=rec, ift=ift_res, ift_dmax_sweep=sweep, mw=mw,
               shape=denss_res, n_frames=len(sam_files) + len(ctl_files),
               n_frame_outliers=sum(1 for r in frames if r.get("outlier"))),
          os.path.join(ctx.out, "summary.json"))
    # human-readable guide to this result folder (批处理时每个样品一份)
    try:
        import readme_for_run
        readme_for_run.write_readme(ctx.out)
        log("README.md -> %s (结果说明与判据)" % os.path.join(ctx.out, "README.md"))
    except Exception as e:                                   # never fail a run over this
        log("WARNING: README.md not written: %s" % e)
    try:                            # figures come last: a plotting error must not lose numbers
        plots(ctx, sub, rows, rec, ift, sweep)
    except Exception as e:
        log("WARNING: plots failed (%s) - summary.json / README.md are already on disk" % e)
    log("done -> %s" % ctx.out)


if __name__ == "__main__":
    main()
