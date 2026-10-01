#!/usr/bin/env python
"""
Embed a high-resolution model (PDB) into the two ab-initio representations of one
SAXS sample:

  bead model       dummy-atom DAM, built here with GNOM -> DAMMIF xN -> DAMAVER
                   (or taken from <sample>/models/ if the sample already has one)
  electron density DENSS .mrc produced by the upstream pipeline

Embedding is rigid-body superposition, not refinement:

  bead model       ATSAS CIFSUP (the ATSAS >= 4 successor of SUPCOMB):
                   template = bead model, movable = atomic model
                   score    = NSD (normalised spatial discrepancy; > 1 means the two
                              objects differ systematically, -> 0 means identical)
  electron density ChimeraX `fitmap` with a global search of the atomic model in the map
                   score    = map-in-map correlation / correlation-about-mean / overlap

Both branches report their score AND their ambiguity diagnostics (CIFSUP selection
sensitivity; fitmap pose spread across resolutions).  The scores say whether the model
sits inside the envelope; they do NOT say whether the model fits the data - that is
`fit-a-high-resolution-model-to-data` (CRYSOL / PDB2SAS).

Run with RAW's own interpreter (pyFAI/numba/matplotlib present):

    /Applications/BioXTASRAW/bin/python embed-model.py <sample-dir> --model A5.pdb \
        --out-dir <sample-dir>/embed

Outputs (in --out-dir): bead/ (GNOM .out, DAMMIF models, DAMAVER consensus),
<stem>_in_beads.pdb + _in_beads.cif, <stem>_in_density_r<res>.pdb, fitmap_r*.csv,
cifsup_*.log, embed_results.json, and two PyMOL panels (one per frame):
<stem>_in_beads.png + <stem>_in_density.png
"""
from __future__ import annotations

import argparse
import csv as csvmod
import glob
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

ATSAS_DEFAULT = "/Applications/ATSAS-4.1.4-1"
CHIMERA_DEFAULT = "/Applications/ChimeraX-1.9.app/Contents/bin/ChimeraX"
PYMOL_DEFAULT = "/Applications/PyMOL.app/Contents/bin/python3.10"


def log(*a):
    print("[%s]" % time.strftime("%H:%M:%S"), *a, flush=True)


def jdump(obj, path):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=str)


def atsas_bin(atsas_dir):
    """RAWAPI wants the ATSAS *bin* directory; the ATSAS environment variable used by the
    ATSAS binaries themselves wants the install root."""
    p = os.path.expanduser(atsas_dir).rstrip("/")
    return p if os.path.basename(p) == "bin" else os.path.join(p, "bin")


def atsas_root(atsas_dir):
    p = os.path.expanduser(atsas_dir).rstrip("/")
    return os.path.dirname(p) if os.path.basename(p) == "bin" else p


# --------------------------------------------------------------- sample metadata
def read_sample_meta(sample_dir):
    """Pull the IFT / P(r) values the upstream pipeline decided on.

    Preference: tables/ift_summary.json -> summary.json.  Only a *trusted* run is used:
    GNOM must not be fed the q-range that the Guinier/IFT gates rejected - with the low-q
    up-turn included GNOM reports a bad discrepancy measure and a chi^2 in the hundreds
    (measured on A5-05-1: chi2 338, DISCRP 18.3, "a SUSPICIOUS solution" vs chi2 0.9,
    "a REASONABLE solution" from the trusted window).
    """
    meta = {"sample_dir": os.path.abspath(sample_dir), "source": None}
    runs, chosen_tag = None, None
    for cand in ("tables/ift_summary.json", "summary.json"):
        p = os.path.join(sample_dir, cand)
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        blk = d.get("ift", d) or {}
        r = blk.get("runs")
        if not r:
            continue
        runs, chosen_tag, meta["source"] = r, blk.get("chosen"), cand
        break
    meta["runs"] = runs or []
    if runs:
        chosen = next((r for r in runs if r.get("tag") == chosen_tag and r.get("trusted")), None)
        chosen = chosen or next((r for r in runs if r.get("tag") == chosen_tag), None) or runs[-1]

        def g(*ks):
            return next((chosen[k] for k in ks if chosen.get(k) is not None), None)

        meta.update(dmax=g("dmax", "Dmax"), rg=g("rg_realspace", "Rg_realspace"),
                    qmin=g("qmin"), qmax=g("qmax"), idx_min=g("idx_min"),
                    trusted=bool(chosen.get("trusted")), tag=chosen.get("tag"))
    meta["dat"] = os.path.join(sample_dir, "profiles", "03_subtracted", "subtracted.dat")
    meta["density_map"] = None
    for cand in ("models/denss.mrc", "models/denss_current.mrc"):   # the final map, not _support
        if os.path.exists(os.path.join(sample_dir, cand)):
            meta["density_map"] = os.path.join(sample_dir, cand)
            break
    meta["existing_bead_models"] = sorted(
        glob.glob(os.path.join(sample_dir, "models", "*damfilt*"))
        + glob.glob(os.path.join(sample_dir, "models", "*damaver*"))
        + glob.glob(os.path.join(sample_dir, "models", "*bead*")))
    return meta


# --------------------------------------------------------------- bead model branch
def build_bead_model(meta, out_bead, args):
    """GNOM (RAWAPI, trusted q-window) -> DAMMIF xN -> DAMAVER.  Returns (bead_model, info)."""
    sys.path.insert(0, "/Applications/BioXTASRAW/lib/python3.12/site-packages")
    import bioxtasraw.RAWAPI as raw

    abin = atsas_bin(args.atsas_dir)
    prefix = os.path.basename(meta["sample_dir"].rstrip("/"))
    os.makedirs(out_bead, exist_ok=True)
    settings = raw.load_settings(meta["cfg"]) if meta.get("cfg") else raw.load_settings(None)
    settings.set("ATSASDir", abin)

    prof = raw.load_profiles([meta["dat"]], settings=settings)[0]
    q = np.asarray(prof.getQ())
    dmax = float(args.dmax if args.dmax is not None else meta["dmax"])
    rg = float(args.rg if args.rg is not None else meta["rg"])
    i0 = int(meta.get("idx_min") or 0)
    i1 = int(np.searchsorted(q, meta["qmax"]) - 1) if meta.get("qmax") else len(q) - 1
    info = {"gnom": dict(dmax=dmax, rg=rg, idx_min=i0, qmin=float(q[i0]), qmax=float(q[i1]),
                         n_points=i1 - i0 + 1, trusted=meta.get("trusted"), tag=meta.get("tag"),
                         source=meta.get("source"))}

    out_name = "%s.out" % prefix
    t = time.time()
    (ift, g_dmax, g_rg, g_i0, rg_err, i0_err, total_est,
     chi_sq, alpha, quality) = raw.gnom(prof, dmax, rg=rg, idx_min=i0, idx_max=i1,
                                        settings=settings, atsas_dir=abin, save_ift=True,
                                        savename=out_name, datadir=out_bead)
    info["gnom"].update(dmax_out=g_dmax, rg_out=g_rg, i0_out=g_i0, rg_err=rg_err, i0_err=i0_err,
                        total_estimate=total_est, chi_sq=chi_sq, alpha=alpha, quality=quality,
                        seconds=time.time() - t)
    log("GNOM: Dmax=%.1f A  Rg(real)=%.2f A  chi2=%.2f  total estimate=%.2f (%s)"
        % (g_dmax, g_rg, chi_sq, total_est, quality))

    runs, t = [], time.time()
    for i in range(1, args.n_models + 1):
        r = raw.dammif(ift, "%s_dammif_%02d" % (prefix, i), out_bead, mode=args.dammif_mode,
                       symmetry=args.symmetry, model_format=args.model_format,
                       write_ift=False, ift_name=out_name, settings=settings, atsas_dir=abin)
        runs.append(dict(chi2=r[0], rg=r[1], dmax=r[2], mw=r[3], excluded_volume=r[4]))
        log("DAMMIF %2d/%d: chi2=%.3f  Rg=%.2f A  Dmax=%.1f A" % (i, args.n_models, r[0], r[1], r[2]))
    info["dammif"] = dict(mode=args.dammif_mode, symmetry=args.symmetry, n_models=args.n_models,
                          model_format=args.model_format, seconds=time.time() - t, runs=runs)

    ext = "." + args.model_format
    files = sorted(os.path.basename(f) for f in
                   glob.glob(os.path.join(out_bead, "%s_dammif_*%s" % (prefix, ext))))
    t = time.time()
    raw.damaver(files, prefix, out_bead, symmetry=args.symmetry,
                model_format=args.model_format, settings=settings, atsas_dir=abin)
    dist = os.path.join(out_bead, "%s-distances.txt" % prefix)
    nsd_mean = nsd_std = None
    if os.path.exists(dist):
        for line in open(dist):
            if "Mean value of nsd" in line:
                nsd_mean = float(line.split(":")[1])
            if "Standard deviation of nsd" in line:
                nsd_std = float(line.split(":")[1])
    info["damaver"] = dict(inputs=files, seconds=time.time() - t, mean_nsd=nsd_mean,
                           std_nsd=nsd_std, distances=dist,
                           summary=os.path.join(out_bead, "%s-global-summary.txt" % prefix),
                           clusters=sorted(os.path.basename(f) for f in
                                           glob.glob(os.path.join(out_bead, "%s-cluster*summary.txt" % prefix))))
    # the consensus bead model to embed into: DAMAVER "filtered" average (most probable model)
    pick = None
    for pat in ("%s-global-damfilt%s", "%s-cluster001-damfilt%s", "%s-global-damaver%s",
                "%s-cluster001-damaver%s"):
        hits = sorted(glob.glob(os.path.join(out_bead, pat % (prefix, ext))))
        if hits:
            pick = hits[0]
            break
    if pick is None:
        pick = sorted(glob.glob(os.path.join(out_bead, "%s_dammif_*%s" % (prefix, ext))))[0]
    info["bead_model"] = pick
    info["bead_model_beads"] = sum(1 for ln in open(pick) if ln[:4] in ("ATOM", "HETA"))
    log("DAMAVER: mean NSD=%.2f (sd %.2f), %d cluster summaries -> consensus bead model %s (%d beads)"
        % (nsd_mean if nsd_mean is not None else -1, nsd_std if nsd_std is not None else -1,
           len(info["damaver"]["clusters"]), os.path.basename(pick), info["bead_model_beads"]))
    return pick, info



def bead_provenance(bead_dir, prefix, bead_model):
    """Recover what is knowable about a REUSED bead model straight from its files, so a re-run
    does not silently drop the provenance of the model it embeds into."""
    info = {"reused": True}
    models = [f for f in sorted(glob.glob(os.path.join(bead_dir, "%s_dammif_*" % prefix)))
              if f.endswith((".cif", ".pdb"))]
    info["n_dammif_models_on_disk"] = len(models)
    dist = os.path.join(bead_dir, "%s-distances.txt" % prefix)
    if os.path.exists(dist):
        for line in open(dist):
            if "Mean value of nsd" in line:
                info["mean_nsd"] = float(line.split(":")[1])
            if "Standard deviation of nsd" in line:
                info["std_nsd"] = float(line.split(":")[1])
    info["clusters"] = sorted(os.path.basename(f) for f in
                              glob.glob(os.path.join(bead_dir, "%s-cluster*summary.txt" % prefix)))
    if bead_model and os.path.exists(bead_model):
        info["bead_model_beads"] = sum(1 for ln in open(bead_model) if ln[:4] in ("ATOM", "HETA"))
        if bead_model.lower().endswith(".cif"):
            for line in open(bead_model):
                if "Dummy atom radius" in line:
                    try:
                        info["dummy_atom_radius_A"] = float(line.split()[-1])
                    except ValueError:
                        pass
                if line.startswith("ATOM"):
                    break
    info["bead_model"] = bead_model
    return info


# --------------------------------------------------------------- CIFSUP branch
def cifsup(atsas_dir, template, movable, out_path, selection, beads=None):
    """Runs CIFSUP; returns the score written into the output file.

    CIFSUP has NO runtime output (manual: "CIFSUP does not have any runtime output"), so the
    NSD has to be read back from the header of the file it writes:
        _atsas_superposition  ...  score <NSD>
    """
    exe = shutil.which("cifsup") or os.path.join(atsas_bin(atsas_dir), "cifsup")
    cmd = [exe, "--method=NSD", "--selection=%s" % selection]
    if beads:
        cmd += ["--beads=%d" % beads]
    cmd += ["-o", out_path, template, movable]
    p = subprocess.run(cmd, capture_output=True, text=True,
                       env=dict(os.environ, ATSAS=atsas_root(atsas_dir)))
    txt = (p.stdout or "") + (p.stderr or "")
    score = None
    if os.path.exists(out_path):
        for line in open(out_path):
            if line.startswith("score"):
                score = float(line.split()[1])
                break
    return score, txt


def align_to_beads(bead_model, model, out_dir, stem, args):
    res = {"bead_model": bead_model, "selection": args.cifsup_selection, "method": "NSD"}
    pdb_out = os.path.join(out_dir, "%s_in_beads.pdb" % stem)
    cif_out = os.path.join(out_dir, "%s_in_beads.cif" % stem)
    nsd, txt = cifsup(args.atsas_dir, bead_model, model, cif_out, args.cifsup_selection,
                      args.cifsup_beads)
    if nsd is None:
        raise RuntimeError("CIFSUP produced no output for selection=%s: %s"
                           % (args.cifsup_selection, " ".join(txt.split())[:300]))
    res["nsd"] = nsd
    res["aligned_model"] = pdb_out
    res["aligned_model_cif"] = cif_out
    open(os.path.join(out_dir, "cifsup_%s.log" % args.cifsup_selection.lower()), "w").write(txt)
    # a PDB copy of the same superposition, for viewers that want PDB
    cifsup(args.atsas_dir, bead_model, model, pdb_out, args.cifsup_selection, args.cifsup_beads)
    # sensitivity of the score to the CIFSUP selection (documented as info, not as a verdict)
    table = {}
    for sel in ("ALL", "BACKBONE", "REGRID", "SHELL"):
        tmp = os.path.join(out_dir, "cifsup_%s.cif" % sel.lower())
        s, _ = cifsup(args.atsas_dir, bead_model, model, tmp, sel, args.cifsup_beads)
        table[sel] = s
        if sel != args.cifsup_selection and os.path.exists(tmp):
            os.remove(tmp)
    res["nsd_by_selection"] = table
    res["beads_in_template"] = sum(1 for ln in open(bead_model) if ln[:4] in ("ATOM", "HETA"))
    res["regrid_beads"] = args.cifsup_beads
    # geometric sanity check NSD cannot give: is the model actually inside the bead lattice?
    try:
        mca = ca_coords(model)                       # works on the input PDB (same frame as the output)
        bcd = coords_any(bead_model)
        d = np.linalg.norm(mca[:, None, :] - bcd[None, :, :], axis=-1).min(1)
        res["ca_to_nearest_bead"] = dict(
            median_A=float(np.median(d)), max_A=float(d.max()),
            frac_within_8A=float((d <= 8).mean()),
            n_beads=int(bcd.shape[0]), n_ca=int(mca.shape[0]))
    except Exception as e:
        res["ca_to_nearest_bead"] = "failed: %s" % e
    log("CIFSUP NSD (%s) = %s | by selection: %s | CA->bead median %s A"
        % (args.cifsup_selection, nsd, {k: (round(v, 3) if v else v) for k, v in table.items()},
           res["ca_to_nearest_bead"].get("median_A") if isinstance(res["ca_to_nearest_bead"], dict) else "?"))
    return res


# --------------------------------------------------------------- density branch
def ca_coords(pdb):
    out = []
    for line in open(pdb):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(out)


def coords_any(path):
    """Coordinates of any ATOM record, for both PDB (fixed columns) and ATSAS mmCIF
    (whitespace tokens; dummy atoms are written as `ATOM 1 C CA . ASP A 1 ? x y z ...`)."""
    cif = path.lower().endswith(".cif")
    out = []
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        if cif:
            t = line.split()
            out.append([float(t[9]), float(t[10]), float(t[11])])
        else:
            out.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(out)


def chimera_run(chimera, lines, workdir, name):
    cxc = os.path.join(workdir, name)
    with open(cxc, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    p = subprocess.run([chimera, "--nogui", "--exit", "--script", cxc],
                       capture_output=True, text=True, cwd=workdir)
    return (p.stdout or "") + (p.stderr or "")


def dock_in_density(map_path, model, out_dir, stem, args):
    """ChimeraX fitmap with a global search, repeated over a resolution ladder.

    Why a global search: fitmap is a *local* optimiser, so the result depends on where the
    model happens to sit when the fit starts (and fitting atoms without `resolution` gives
    map-average metrics only - the correlation column stays None).  Why a ladder: the model
    is Gaussian-smeared at `resolution`, and at SAXS envelope resolutions a wrong pose can
    score as well as the right one (measured on A5-05-1: a pose 28 A away from the best one
    tied it, 0.9096 vs 0.9122, at 20 A).  The pose spread across the ladder is therefore
    part of the result, not an extra.
    """
    res_list = [float(x) for x in str(args.fitmap_resolutions).split(",")]
    fits = []
    for res in res_list:
        tag = "r%02d" % int(res)
        pdb = os.path.join(out_dir, "%s_in_density_%s.pdb" % (stem, tag))
        csvf = os.path.join(out_dir, "fitmap_%s.csv" % tag)
        txt = chimera_run(args.chimera, [
            "open %s" % map_path,
            "open %s" % model,
            "view initial",             # doc: reset the model/view state before a global search
            "fitmap #2 in #1 resolution %g metric %s search %d seed %d logFits %s listFits true"
            % (res, args.fitmap_metric, args.fitmap_search, args.fitmap_seed, csvf),
            "save %s models #2" % pdb,
        ], out_dir, "fitmap_%s.cxc" % tag)
        open(os.path.join(out_dir, "fitmap_%s.log" % tag), "w").write(txt)
        rows = list(csvmod.DictReader(open(csvf), delimiter=" ")) if os.path.exists(csvf) else []
        best = max(rows, key=lambda r: float(r.get("correlation") or -9)) if rows else None
        entry = dict(resolution=res, n_unique_fits=len(rows), csv=csvf,
                     log=os.path.join(out_dir, "fitmap_%s.log" % tag))
        if best:
            entry.update(aligned_model=pdb,
                         correlation=float(best["correlation"]),
                         correlation_about_mean=float(best["correlation_about_mean"]),
                         overlap=float(best["overlap"]),
                         average_map_value=float(best["average_map_value"]),
                         points=int(best["points"]), shift=float(best["shift"]),
                         angle=float(best["angle"]), steps=int(best["steps"]),
                         top3_correlations=[round(float(r["correlation"]), 4)
                                            for r in sorted(rows, key=lambda r: -float(r["correlation"]))[:3]])
        fits.append(entry)
        log("fitmap r=%4.1f A: best correlation=%s (%d unique fits) top3=%s"
            % (res, entry.get("correlation"), entry["n_unique_fits"], entry.get("top3_correlations")))
    scored = [f for f in fits if f.get("correlation") is not None]
    out = {"map": map_path, "metric": args.fitmap_metric, "search": args.fitmap_search,
           "seed": args.fitmap_seed, "fits": fits}
    if scored:
        best = max(scored, key=lambda f: f["correlation"])
        out.update(best_resolution=best["resolution"], best_correlation=best["correlation"],
                   best_correlation_about_mean=best["correlation_about_mean"],
                   best_aligned_model=best["aligned_model"], n_unique_fits=best["n_unique_fits"])
        ref = ca_coords(best["aligned_model"])
        spread = {}
        for f in scored:
            if f.get("aligned_model") and os.path.exists(f["aligned_model"]):
                c = ca_coords(f["aligned_model"])
                if c.shape == ref.shape:
                    spread["r%02d" % int(f["resolution"])] = float(
                        np.sqrt(((c - ref) ** 2).sum(1).mean()))
        out["pose_rmsd_vs_best_A"] = spread
        out["pose_spread_max_A"] = max(spread.values()) if spread else None
        if out["pose_spread_max_A"] and out["pose_spread_max_A"] > 5:
            out["ambiguity"] = ("poses from different resolutions differ by %.1f A CA-RMSD: the "
                                "envelope supports more than one placement - do not quote a single "
                                "pose as the answer" % out["pose_spread_max_A"])
    return out


def render_panel(pymol_python, render_script, fitted_pdb, overlay, out_png, kind, param, workdir):
    """Render ONE panel with PyMOL.

    ChimeraX --nogui cannot save images ("Unable to save images because OpenGL rendering is
    not available"), so figures go through PyMOL.  Two panels, never one: the bead model and
    the DENSS map each live in their own frame (centred on their own origin, relative
    orientation unconstrained), so drawing them together fakes a "model sticking out of the
    beads" that is not a measurement.
    """
    cmd = [pymol_python, render_script, fitted_pdb, overlay, out_png, kind, str(param)]
    p = subprocess.run(cmd, capture_output=True, text=True,
                       cwd=workdir, env={k: v for k, v in os.environ.items() if k != "PYTHONPATH"})
    if not os.path.exists(out_png):
        raise RuntimeError("PyMOL produced no image: %s" % " ".join(
            ((p.stdout or "") + (p.stderr or "")).split())[:300])
    return out_png


# --------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(
        description="Embed a high-resolution model into a bead model and into a DENSS "
                    "electron-density map (one SAXS sample directory).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("sample_dir", help="processed sample directory, e.g. .../Tube-SAXS/A5-05-1")
    ap.add_argument("--model", required=True, help="high-resolution model to embed (.pdb)")
    ap.add_argument("--out-dir", default=None, help="output directory (default: <sample_dir>/embed)")
    ap.add_argument("--bead-model", default=None,
                    help="use this existing bead model instead of building one")
    ap.add_argument("--dmax", type=float, default=None,
                    help="Dmax (A) for GNOM/DAMMIF (default: trusted IFT value from the sample tables)")
    ap.add_argument("--rg", type=float, default=None,
                    help="Rg (A) for GNOM (default: trusted IFT value from the sample tables)")
    ap.add_argument("--n-models", type=int, default=15,
                    help="DAMMIF reconstructions to average (10-20; 15 is the ATSAS recommendation)")
    ap.add_argument("--dammif-mode", default="Fast", choices=["Fast", "Slow"],
                    help="DAMMIF mode (Fast for triage, Slow for figures)")
    ap.add_argument("--symmetry", default="P1", help="DAMMIF/DAMAVER symmetry")
    ap.add_argument("--model-format", default="cif", choices=["cif", "pdb"],
                    help="bead-model file format (ATSAS 4 DAMMIF writes cif; RAWAPI reads cif back)")
    ap.add_argument("--cifsup-selection", default="REGRID",
                    choices=["ALL", "BACKBONE", "REGRID", "SHELL"],
                    help="CIFSUP selection (REGRID with --method=NSD is the documented SUPCOMB fast-mode equivalent)")
    ap.add_argument("--cifsup-beads", type=int, default=None,
                    help="beads for the REGRID representation (CIFSUP default: 2000)")
    ap.add_argument("--fitmap-resolutions", default="10,15,20,25",
                    help="comma-separated resolutions (A) for the fitmap ladder")
    ap.add_argument("--fitmap-metric", default="correlation", choices=["correlation", "overlap"],
                    help="ChimeraX fitmap metric used for map-in-map fitting")
    ap.add_argument("--fitmap-search", type=int, default=200,
                    help="random initial placements for the fitmap global search (0 = local only)")
    ap.add_argument("--fitmap-seed", type=int, default=42,
                    help="random seed of the fitmap global search (fix it for reproducibility)")
    ap.add_argument("--atsas-dir", default=ATSAS_DEFAULT, help="ATSAS installation directory")
    ap.add_argument("--chimera", default=CHIMERA_DEFAULT, help="ChimeraX executable (fitmap)")
    ap.add_argument("--pymol-python", default=PYMOL_DEFAULT,
                    help="PyMOL's python (figures: ChimeraX --nogui cannot save images)")
    ap.add_argument("--render-script", default=None,
                    help="render-embed-figure.py path (default: next to this script)")
    ap.add_argument("--rebuild-bead", action="store_true",
                    help="rebuild the bead model even if <out>/bead/ already has one")
    ap.add_argument("--skip-bead", action="store_true", help="skip the bead-model branch")
    ap.add_argument("--skip-density", action="store_true", help="skip the electron-density branch")
    ap.add_argument("--no-figure", action="store_true", help="do not render the two PyMOL panels")
    args = ap.parse_args()

    sample_dir = os.path.abspath(os.path.expanduser(args.sample_dir))
    out = os.path.abspath(os.path.expanduser(args.out_dir or os.path.join(sample_dir, "embed")))
    os.makedirs(out, exist_ok=True)
    model = os.path.abspath(os.path.expanduser(args.model))
    if not os.path.exists(model):
        raise SystemExit("model not found: %s" % model)
    stem = os.path.splitext(os.path.basename(model))[0]
    results = {"sample": os.path.basename(sample_dir.rstrip("/")), "model": model, "out_dir": out}

    meta = read_sample_meta(sample_dir)
    cfgs = glob.glob(os.path.join(sample_dir, "..", "..", "..", "data", "*.cfg"))
    meta["cfg"] = sorted(cfgs)[-1] if cfgs else None
    results["input"] = {k: meta.get(k) for k in ("dat", "density_map", "dmax", "rg", "qmin", "qmax",
                                                 "idx_min", "tag", "trusted", "source")}
    if not meta.get("trusted", False):
        results["warning"] = ("the sample's IFT was NOT trusted (gates failed): the bead model is "
                              "built from a P(r) the upstream pipeline rejected - treat the NSD as "
                              "indicative only")
    log("%s | IFT %s (trusted=%s) Dmax=%.1f A Rg=%.2f A | density=%s"
        % (results["sample"], meta.get("tag"), meta.get("trusted"), meta.get("dmax") or -1,
           meta.get("rg") or -1, os.path.basename(meta.get("density_map") or "-")))

    if not args.skip_bead:
        try:
            bead_dir = os.path.join(out, "bead")
            cached = [] if args.rebuild_bead else sorted(
                glob.glob(os.path.join(bead_dir, "*global-damfilt*"))
                + glob.glob(os.path.join(bead_dir, "*damaver*damfilt*")))
            if args.bead_model:
                bead, binfo = args.bead_model, {"bead_model": args.bead_model, "reused": True}
            elif cached:
                bead = cached[0]
                binfo = bead_provenance(bead_dir, os.path.basename(meta["sample_dir"].rstrip("/")), bead)
                log("reusing bead model %s (use --rebuild-bead to rebuild)" % os.path.basename(bead))
            elif meta["existing_bead_models"]:
                bead = meta["existing_bead_models"][0]
                binfo = {"bead_model": bead, "reused": True}
            else:
                bead, binfo = build_bead_model(meta, bead_dir, args)
            results["bead_model"] = binfo
            results["beads_fit"] = align_to_beads(bead, model, out, stem, args)
        except Exception as e:
            results["beads_fit"] = {"status": "failed: %s" % e}
            log("bead branch failed: %s" % e)

    if not args.skip_density:
        if not meta["density_map"]:
            results["density_fit"] = {"status": "no DENSS map under <sample>/models/"}
            log("density branch skipped: no DENSS map")
        else:
            try:
                results["density_fit"] = dock_in_density(meta["density_map"], model, out, stem, args)
            except Exception as e:
                results["density_fit"] = {"status": "failed: %s" % e}
                log("density branch failed: %s" % e)

    if (not args.no_figure and results.get("bead_model")
            and results.get("density_fit", {}).get("best_aligned_model")):
        rs = args.render_script or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                "render-embed-figure.py")
        figs = {}
        try:
            figs["beads"] = render_panel(
                args.pymol_python, rs, results["beads_fit"]["aligned_model"],
                results["bead_model"]["bead_model"],
                os.path.join(out, "%s_in_beads.png" % stem), "beads", 0.5, out)
            figs["density"] = render_panel(
                args.pymol_python, rs, results["density_fit"]["best_aligned_model"],
                meta["density_map"],
                os.path.join(out, "%s_in_density.png" % stem), "density", 0.02, out)
            results["figures"] = figs
            log("figures: %s , %s" % (figs["beads"], figs["density"]))
        except Exception as e:
            results["figures"] = "failed: %s" % e
            log("figures failed: %s" % e)

    jdump(results, os.path.join(out, "embed_results.json"))
    log("wrote %s" % os.path.join(out, "embed_results.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
