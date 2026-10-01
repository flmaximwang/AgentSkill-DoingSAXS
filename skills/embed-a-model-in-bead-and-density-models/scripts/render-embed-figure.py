#!/usr/bin/env python
"""Render ONE panel of the embed figure with PyMOL (headless) and save a .pse session.

The panel style is fixed by the skill's delivery spec, so every delivered picture and
session reads the same way:

    high-resolution model        -> yellow cartoon    (object `sample`)
    bead model / DENSS density   -> translucent white (object `beads` / `density`)

The .pse is written next to the .png (same basename) unless --pse says otherwise; open it
with `pymol <file>.pse`.  ChimeraX --nogui cannot save images ("Unable to save images
because OpenGL rendering is not available"), which is why figures go through PyMOL.
"""
from __future__ import annotations

import argparse
import os
import sys

import pymol

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mrcmap  # noqa: E402  (ships next to this script)

# --overlay-transparency default depends on the representation: overlapping spheres pile up
# (every layer multiplies the transmitted light), so the bead envelope needs a much larger
# value than a single isosurface shell before the sample stays visible through it.
DEFAULT_TRANSPARENCY = {"beads": 0.85, "density": 0.7}


def bead_radius_from_cif(path):
    """ATSAS writes the dummy-atom radius into the model's own header
    (`_atsas_dummy_atom_model.value`, e.g. 2.200).  Drawing the beads at that radius makes
    the envelope the model's actual bead volume instead of a scaled-up vdW guess."""
    for line in open(path):
        if line.startswith("_atsas_dummy_atom_model.value"):
            return float(line.split()[1].strip("'"))
        if line.startswith("ATOM"):
            break
    return None


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Render one embed panel (sample = yellow cartoon, envelope = translucent "
                    "white) and save the matching .pse session.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("fitted", help="model after embedding (.pdb), in the overlay's frame")
    ap.add_argument("overlay", help="bead model (.cif/.pdb) or DENSS map (.mrc)")
    ap.add_argument("out_png", help="output image (.png)")
    ap.add_argument("kind", choices=["beads", "density"],
                    help="what the overlay is: dummy-atom spheres, or a map isosurface")
    ap.add_argument("--level", default="auto",
                    help="isosurface level for kind=density, in map units; 'auto' takes the "
                         "level whose enclosed volume matches DENSS's own support volume "
                         "(denss.log 'Final Support Volume'), falling back to 0.02")
    ap.add_argument("--sphere-scale", type=float, default=0.5,
                    help="PyMOL sphere_scale for kind=beads, used only when the bead radius "
                         "cannot be taken from the model")
    ap.add_argument("--bead-radius", type=float, default=None,
                    help="draw dummy atoms at this radius (A) instead of scaling vdW radii "
                         "(default: the model's own _atsas_dummy_atom_model.value)")
    ap.add_argument("--overlay-transparency", type=float, default=None,
                    help="0=opaque, 1=invisible (default: %s for beads, %s for density)"
                         % (DEFAULT_TRANSPARENCY["beads"], DEFAULT_TRANSPARENCY["density"]))
    ap.add_argument("--transparency-mode", type=int, default=1,
                    help="PyMOL transparency_mode; mode 2 (PyMOL's default) silently drops "
                         "sphere transparency in ray tracing, so bead spheres come out opaque")
    ap.add_argument("--smooth", default="auto",
                    help="Gaussian sigma to apply to the map before isosurfacing, in Angstrom "
                         "(0 = off); 'auto' uses one voxel of the map, which removes the voxel "
                         "facets of a coarse DENSS map without moving the envelope's volume "
                         "(the level is re-derived from the smoothed map)")
    ap.add_argument("--sample-color", default="yellow", help="colour of the atomic model")
    ap.add_argument("--overlay-color", default="white", help="colour of the envelope")
    ap.add_argument("--bg", default="black",
                    help="background colour; black, because a translucent white envelope on a "
                         "white background is nearly invisible")
    ap.add_argument("--size", type=int, default=1200, help="ray-traced image size (px)")
    ap.add_argument("--zoom-buffer", type=float, default=8.0,
                    help="margin (A) left around the drawn objects when framing the shot")
    ap.add_argument("--ray-trace-mode", type=int, default=0,
                    help="PyMOL ray_trace_mode (1 = black outlines, 0 = off)")
    ap.add_argument("--ss-ref", default=None,
                    help="model whose secondary structure is transferred onto the fitted PDB")
    ap.add_argument("--pse", default=None,
                    help="session file to write (default: <out_png> with a .pse suffix)")
    ap.add_argument("--no-pse", action="store_true", help="do not write the session file")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    t = args.overlay_transparency
    if t is None:
        t = DEFAULT_TRANSPARENCY[args.kind]
    pse = args.pse or os.path.splitext(args.out_png)[0] + ".pse"

    overlay, level, m = args.overlay, None, None
    if args.kind == "density":
        m = mrcmap.read(args.overlay)
        support = mrcmap.denss_support_volume(
            os.path.join(os.path.dirname(os.path.abspath(args.overlay)), "denss.log"))
        sigma = float(m["voxel"].mean()) if args.smooth == "auto" else float(args.smooth)
        if sigma > 0:
            m["data"] = mrcmap.smooth(m, sigma)
            # keep the smoothed copy beside the figure, so the rendered surface is reproducible
            overlay = os.path.join(
                os.path.dirname(os.path.abspath(args.out_png)),
                "%s_map_smooth%gA.mrc" % (os.path.splitext(os.path.basename(args.out_png))[0],
                                          round(sigma, 2)))
            mrcmap.write_like(m, overlay, m["data"])
        # the level is read off the map that actually gets contoured (smoothing moves the
        # volume-level curve, so a level derived from the raw map over-inflates the surface)
        if args.level != "auto":
            level = float(args.level)
        else:
            level = mrcmap.level_for_volume(m, support) if support else 0.02
        vol = mrcmap.volume_above(m, level)
        print("density: level=%.5f | enclosed %.0f A^3 | denss.log support %s A^3 (%s) | "
              "smooth %g A | voxel %.3g A"
              % (level, vol, ("%.0f" % support) if support else "n/a",
                 ("%.1f%%" % (100.0 * vol / support)) if support else "no denss.log",
                 sigma, float(m["voxel"].mean())))

    pymol.finish_launching(["pymol", "-cq"])
    from pymol import cmd

    cmd.load(args.fitted, "sample")
    cmd.load(overlay, "density_map" if args.kind == "density" else "beads")
    if args.kind == "beads":
        cmd.hide("everything", "beads")
    try:
        cmd.dss("sample")
    except Exception as e:
        print("dss failed:", e)
    # CIFSUP's PDB output can defeat PyMOL's DSS (0 helices found, measured); take the
    # assignment from the original model, which PyMOL reads correctly.
    try:
        nH = []
        cmd.iterate("sample and name CA and ss H", "nH.append(resi)", space={"nH": nH})
        if args.ss_ref and os.path.exists(args.ss_ref) and len(nH) < 3:
            cmd.load(args.ss_ref, "ssref")
            cmd.dss("ssref")
            pairs = []
            cmd.iterate("ssref and name CA", "pairs.append((chain,resi,ss))",
                        space={"pairs": pairs})
            for chain, resi, s in pairs:
                if s in ("H", "E"):
                    cmd.alter("sample and chain %s and resi %s" % (chain, resi), "ss='%s'" % s)
            cmd.delete("ssref")
            print("secondary structure transferred from", args.ss_ref)
    except Exception as e:
        print("ss transfer failed:", e)

    # mode 2 (PyMOL's default) drops sphere transparency in ray tracing, so bead panels come
    # out opaque and hide the sample (measured 2026-10-01, mean |pixel diff| against an opaque
    # render at sphere_transparency 0.8: 2.8/255 in mode 2 vs 8.4/255 in modes 0/1/3).
    cmd.set("transparency_mode", args.transparency_mode)

    if args.kind == "beads":
        cmd.show("spheres", "beads")
        radius = args.bead_radius or bead_radius_from_cif(overlay)
        if radius:
            # draw the beads at the model's own dummy-atom radius (spheres then touch along the
            # 4.4 A lattice instead of overlapping into capsules, which is what a vdW-based
            # sphere_scale produced: measured 2026-10-01 on damaver-cluster001-damaver.cif)
            cmd.alter("beads", "vdw=%f" % radius)
            cmd.set("sphere_scale", 1.0, "beads")
            print("beads: drawn at the model's dummy-atom radius %.3f A" % radius)
        else:
            cmd.set("sphere_scale", args.sphere_scale, "beads")
        cmd.set("sphere_quality", 2, "beads")
        cmd.color(args.overlay_color, "beads")
        # spheres need sphere_transparency: the generic `transparency` setting is surfaces-only
        cmd.set("sphere_transparency", t, "beads")
    else:
        cmd.isosurface("density", "density_map", level)
        cmd.color(args.overlay_color, "density")
        cmd.set("transparency", t, "density")
        cmd.disable("density_map")

    cmd.show("cartoon", "sample")
    cmd.color(args.sample_color, "sample")
    cmd.set("cartoon_transparency", 0.0, "sample")
    cmd.bg_color(args.bg)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_trace_mode", args.ray_trace_mode)
    cmd.set("antialias", 2)
    cmd.set("cartoon_side_chain_helper", 1)
    # framing: zoom() is an atom selection, so it cannot see an isosurface object (a surface
    # has no atoms).  The density panel therefore gets a temporary 8-pseudoatom box around the
    # envelope - the map's own extent (in PyMOL units) plus the voxels above the level give the
    # box - because zoom() on atoms is also what sets the clip planes correctly: a hand-rolled
    # "dolly out" by scaling the camera position leaves the clip planes behind and the shot
    # comes out fogged/clipped (measured at 1.35x: the whole frame went dark).
    framed = False
    if args.kind == "density" and m is not None:
        # get_extent returns [[xlo,ylo,zlo],[xhi,yhi,phi]] spanning the voxel *centres*
        # (measured: a 32-voxel / 217.7 A brick reports 210.9 A = 31 intervals), so the
        # voxel-index -> XYZ map is lo + idx/(N-1)*(hi-lo), no half-voxel shift.
        ext = np.array(cmd.get_extent("density_map"), dtype=float)
        idx = np.argwhere(m["data"] > level)
        if ext.shape == (2, 3) and len(idx):
            frac = idx / (np.array(m["data"].shape)[None, :] - 1.0)
            xyz = ext[0] + frac * (ext[1] - ext[0])
            for p in np.vstack([xyz.min(0), xyz.max(0),
                                xyz[[xyz[:, 0].argmax(), xyz[:, 0].argmin()]],
                                xyz[[xyz[:, 1].argmax(), xyz[:, 1].argmin()]],
                                xyz[[xyz[:, 2].argmax(), xyz[:, 2].argmin()]]]):
                cmd.pseudoatom("envelope_box", pos=[float(v) for v in p])
            cmd.hide("everything", "envelope_box")
            cmd.zoom("envelope_box or sample", args.zoom_buffer)
            cmd.delete("envelope_box")
            framed = True
    if not framed:
        have = set(cmd.get_names("objects"))
        drawn = [n for n in ("sample", "beads") if n in have]
        cmd.zoom(" or ".join(drawn) if drawn else "all", args.zoom_buffer)
    cmd.ray(args.size, args.size)
    cmd.png(args.out_png, dpi=150)
    if not args.no_pse:
        cmd.save(pse)
    for p in (args.out_png,) + (() if args.no_pse else (pse,)):
        if not os.path.exists(p):
            raise SystemExit("not written: %s" % p)
    print("wrote %s (%d bytes)" % (args.out_png, os.path.getsize(args.out_png)))
    if not args.no_pse:
        print("wrote %s (%d bytes)" % (pse, os.path.getsize(pse)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
