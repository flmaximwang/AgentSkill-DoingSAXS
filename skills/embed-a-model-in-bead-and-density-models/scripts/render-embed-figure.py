#!/usr/bin/env python
"""Render one panel of the embed figure with PyMOL (headless).

ChimeraX --nogui cannot save images ("Unable to save images because OpenGL rendering is
not available"), so the figure is rendered with PyMOL.

    /Applications/PyMOL.app/Contents/bin/python3.10 render-embed-figure.py \
        <fitted_model.pdb> <overlay> <out.png> <beads|density> [level_or_output_dir]

overlay:
    beads    a dummy-atom bead model (.cif/.pdb) -> drawn as spheres
    density  a DENSS map (.mrc)                   -> drawn as an isosurface at `level`
"""
import os
import sys

import pymol

pymol.finish_launching(["pymol", "-cq"])
from pymol import cmd

fitted, overlay, out, kind = sys.argv[1:5]
param = float(sys.argv[5]) if len(sys.argv) > 5 else (0.02 if kind == "density" else 0.35)
ss_ref = sys.argv[6] if len(sys.argv) > 6 else None

cmd.load(fitted, "prot")
cmd.load(overlay, "ovl")
cmd.hide("everything", "ovl")
try:
    cmd.dss("prot")
except Exception as e:
    print("dss failed:", e)
# CIFSUP's PDB output can defeat PyMOL's DSS (0 helices found, measured); take the assignment
# from the original model, which PyMOL reads correctly.
try:
    nH = []
    cmd.iterate("prot and name CA and ss H", "nH.append(resi)", space={"nH": nH})
    if ss_ref and os.path.exists(ss_ref) and len(nH) < 3:
        cmd.load(ss_ref, "ssref")
        cmd.dss("ssref")
        pairs = []
        cmd.iterate("ssref and name CA", "pairs.append((chain,resi,ss))", space={"pairs": pairs})
        for chain, resi, s in pairs:
            if s in ("H", "E"):
                cmd.alter("prot and chain %s and resi %s" % (chain, resi), "ss='%s'" % s)
        cmd.delete("ssref")
        print("secondary structure transferred from", ss_ref)
except Exception as e:
    print("ss transfer failed:", e)
if kind == "beads":
    cmd.show("spheres", "ovl")
    cmd.set("sphere_scale", param, "ovl")
    cmd.color("gray70", "ovl")
    cmd.set("transparency", 0.45, "ovl")        # otherwise the spheres hide the model
else:
    cmd.isosurface("ovl_surf", "ovl", param)
    cmd.color("palecyan", "ovl_surf")
    cmd.set("transparency", 0.55, "ovl_surf")
cmd.show("cartoon", "prot")
cmd.color("marine", "prot")
cmd.bg_color("white")
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("cartoon_side_chain_helper", 1)
cmd.orient()
cmd.ray(1200, 1200)
cmd.png(out, dpi=150)
print("wrote", out, os.path.exists(out))
