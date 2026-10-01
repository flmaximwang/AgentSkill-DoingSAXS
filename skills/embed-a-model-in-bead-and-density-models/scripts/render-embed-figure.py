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

import pymol

# --overlay-transparency default depends on the representation: overlapping spheres pile up
# (every layer multiplies the transmitted light), so the bead envelope needs a much larger
# value than a single isosurface shell before the sample stays visible through it.
DEFAULT_TRANSPARENCY = {"beads": 0.85, "density": 0.55}


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
    ap.add_argument("--level", type=float, default=0.02,
                    help="isosurface level for kind=density (in map units)")
    ap.add_argument("--sphere-scale", type=float, default=0.5,
                    help="PyMOL sphere_scale for kind=beads")
    ap.add_argument("--overlay-transparency", type=float, default=None,
                    help="0=opaque, 1=invisible (default: %s for beads, %s for density)"
                         % (DEFAULT_TRANSPARENCY["beads"], DEFAULT_TRANSPARENCY["density"]))
    ap.add_argument("--smooth", type=float, default=0.0,
                    help="Gaussian sigma (voxels) applied to the map before isosurfacing, "
                         "0 = off; coarse maps (6-7 A voxels) look faceted without it")
    ap.add_argument("--sample-color", default="yellow", help="colour of the atomic model")
    ap.add_argument("--overlay-color", default="white", help="colour of the envelope")
    ap.add_argument("--size", type=int, default=1200, help="ray-traced image size (px)")
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

    pymol.finish_launching(["pymol", "-cq"])
    from pymol import cmd

    cmd.load(args.fitted, "sample")
    cmd.load(args.overlay, "density_map" if args.kind == "density" else "beads")
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

    if args.kind == "beads":
        cmd.show("spheres", "beads")
        cmd.set("sphere_scale", args.sphere_scale, "beads")
        cmd.set("sphere_quality", 2, "beads")
        cmd.color(args.overlay_color, "beads")
        # spheres need sphere_transparency: the generic `transparency` setting is surfaces-only,
        # which is why bead panels used to come out opaque (measured 2026-10-01).
        cmd.set("sphere_transparency", t, "beads")
    else:
        src = "density_map"
        if args.smooth and args.smooth > 0:
            cmd.map_new("density_smooth", "gaussian", args.smooth, "density_map")
            src = "density_smooth"
        cmd.isosurface("density", src, args.level)
        cmd.color(args.overlay_color, "density")
        cmd.set("transparency", t, "density")
        cmd.disable("density_map")

    cmd.show("cartoon", "sample")
    cmd.color(args.sample_color, "sample")
    cmd.set("cartoon_transparency", 0.0, "sample")
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    cmd.set("ray_trace_mode", args.ray_trace_mode)
    cmd.set("antialias", 2)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.orient()
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
