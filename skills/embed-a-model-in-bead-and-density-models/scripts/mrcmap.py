#!/usr/bin/env python
"""Minimal MRC/CCP4 map I/O for the embed scripts (numpy only, no mrcfile).

Needed for two things the figure step must be able to state:

  * the *volume* enclosed by the displayed isosurface, so the contour level can be checked
    against the support volume DENSS itself reports in `denss.log`;
  * a smoothed copy of a coarse map (DENSS maps from the BL19U2 pipeline are 32^3 on a
    218 A box = 6.8 A voxels, which ray-traces into visible facets).

`write_like` copies the source header verbatim and replaces only the data block, so the
copy stays readable by PyMOL/ChimeraX.
"""
from __future__ import annotations

import struct

import numpy as np

MODES = {0: "i1", 1: "i2", 2: "f4", 6: "u2"}


def read(path):
    raw = open(path, "rb").read()
    nx, ny, nz, mode = struct.unpack_from("<4i", raw, 0)
    cell = struct.unpack_from("<3f", raw, 40)
    nsymbt = struct.unpack_from("<i", raw, 92)[0]
    off = 1024 + max(nsymbt, 0)
    if mode not in MODES:
        raise ValueError("unsupported MRC mode %r in %s" % (mode, path))
    dt = np.dtype("<" + MODES[mode])
    n = nx * ny * nz
    data = np.frombuffer(raw[off:off + n * dt.itemsize], dtype=dt).reshape(nz, ny, nx)
    voxel = np.array([abs(c) / m for c, m in zip(cell, (nx, ny, nz))], dtype=float)
    return {"data": data.astype("f4"), "header": raw[:off], "dims": (nx, ny, nz),
            "voxel": voxel, "mode": mode, "path": path,
            "voxel_volume": float(abs(voxel.prod()))}


def volume_above(m, level):
    """Volume (A^3) enclosed by the isosurface at `level` on map `m`."""
    return float((m["data"] > level).sum()) * m["voxel_volume"]


def level_for_volume(m, target_volume):
    """Largest level whose enclosed volume is still >= target (the volume-level curve is a
    step function on a coarse map, so pick the level of the first voxel below the target
    count instead of interpolating)."""
    n_target = int(round(target_volume / m["voxel_volume"]))
    flat = np.sort(m["data"].ravel())[::-1]
    if n_target <= 0 or n_target >= flat.size:
        return float(flat[min(max(n_target, 0), flat.size - 1)])
    return float(flat[n_target - 1])


def smooth(m, sigma_A):
    """Separable Gaussian smoothing of the data, sigma given in Angstrom."""
    sigma_vox = [max(sigma_A / v, 1e-6) for v in m["voxel"]]
    out = m["data"].astype("f4")
    for axis, s in enumerate(sigma_vox):
        if s <= 0.05:
            continue
        r = int(max(1, round(3 * s)))
        x = np.arange(-r, r + 1, dtype="f4")
        k = np.exp(-0.5 * (x / s) ** 2)
        k /= k.sum()
        shape = [1, 1, 1]
        shape[axis] = k.size
        out = _convolve_axis(out, k.reshape(shape), axis)
    return out


def _convolve_axis(a, k, axis):
    pad = [(0, 0)] * a.ndim
    pad[axis] = (k.size // 2, k.size // 2)
    ap = np.pad(a, pad, mode="edge")
    out = np.zeros_like(a, dtype="f4")
    for i, w in enumerate(k.ravel()):
        sl = [slice(None)] * a.ndim
        sl[axis] = slice(i, i + a.shape[axis])
        out += w * ap[tuple(sl)]
    return out


def write_like(m, path, data):
    with open(path, "wb") as fh:
        fh.write(m["header"])
        fh.write(np.asarray(data, dtype="<f4").tobytes())
    return path


def denss_support_volume(denss_log):
    """The 'Final Support Volume' DENSS printed, in A^3 (None if the log is absent/odd)."""
    try:
        for line in open(denss_log):
            if "Final Support Volume" in line:
                return float(line.split()[-1])
    except OSError:
        return None
    return None
