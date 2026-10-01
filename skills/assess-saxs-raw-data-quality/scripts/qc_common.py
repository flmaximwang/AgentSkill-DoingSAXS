#!/usr/bin/env python
"""Shared helpers for the raw-data assessment scripts (RAW does all the maths).

Everything here is a thin wrapper around bioxtasraw.RAWAPI plus the two conventions this
beamline/dataset layout needs:

  * frame file names are  <run>_<NNNNN>.tif  where <run> = "<sample>_<number>" (e.g.
    A5-05-1_0013_00001.tif -> run key "A5-05-1_0013").  One run = one block of frames
    measured on one filling (usually 20 frames).
  * the sample is the run whose key starts with the directory name; every other run in the
    same directory is a control (buffer) run - prefix match is EXACT after stripping
    punctuation ("877-apo-pb7" contains "pb7" as a buffer suffix, not as a background).

The settings builder forces the two switches a beamline .cfg usually leaves off; without
them the per-frame .txt is never parsed and transmission normalization silently does
nothing (see the pipeline skill's reference for the 2.25x error this caused).

Run these scripts with the RAW interpreter, e.g. /Applications/BioXTASRAW/bin/python.
"""
from __future__ import annotations

import os
import re

import numpy as np
import bioxtasraw.RAWAPI as raw

HDR_FORMAT = "BL19U2, SSRF"
SCALE_WINDOW = (0.30, 0.44)          # q window where the particle does not scatter


def load_settings(cfg, hdr_format=HDR_FORMAT):
    s = raw.load_settings(os.path.expanduser(cfg))
    s.set("ImageHdrFormat", hdr_format)
    s.set("EnableNormalization", True)
    s.set("NormalizationList", [["/", "Transmitted_Beam"]])
    return s


def runkey(path):
    """A5-05-1_0013_00001.tif -> 'A5-05-1_0013'."""
    m = re.match(r"^(.*)_(\d{5})\.tif$", os.path.basename(path))
    return m.group(1) if m else os.path.basename(path)


def _san(s):
    return re.sub(r"[^0-9A-Za-z]", "", s).lower()


def split_runs(files, sample_key):
    """Return (sample_files, {control_run_key: files}, sample_run_key)."""
    runs = {}
    for f in files:
        runs.setdefault(runkey(f), []).append(f)
    if sample_key:
        sam = next((k for k in sorted(runs) if _san(k) == _san(sample_key)), None)
    else:
        sam = None
    if sam is None:                      # no explicit key: the sample is the directory name
        raise ValueError("no run matches %r; runs seen: %s" % (sample_key, sorted(runs)))
    return sorted(runs[sam]), {k: sorted(v) for k, v in runs.items() if k != sam}, sam


def guess_sample_key(files, dirname):
    """sample run key for a sample directory (prefix match on the directory name)."""
    runs = sorted({runkey(f) for f in files})
    want = _san(dirname)
    hit = [k for k in runs if _san(k).startswith(want)]
    return hit[0] if hit else None


def average(files, s):
    profs, _ = raw.load_and_integrate_images(sorted(files), settings=s)
    return raw.average(profs)


def scale_factor(a, b, window=SCALE_WINDOW):
    """mean(a/b) in the flat high-q window: the constant the pipeline removes."""
    q = a.getQ()
    m = (q >= window[0]) & (q <= window[1])
    return float(np.mean(a.getI()[m] / b.getI()[m])) if m.sum() >= 5 else 1.0


def subtract(a, b, window=SCALE_WINDOW):
    """a - scaled(b), RAW's own subtract; returns (difference_profile, factor)."""
    f = scale_factor(a, b, window)
    bc = b.copy()
    if abs(f - 1.0) > 1e-4:
        bc.scaleRelative(f)
    return raw.subtract([a], bc)[0], f


def beam_center(image_shape, s):
    """(x, y) pixel of the beam in a loaded image.

    RAW's Pilatus reader hands back the array flipped in y relative to the .cfg's Ycenter,
    so the row index is (n_rows - 1 - Ycenter).  Check it: the manual radial average with
    this centre must reproduce RAW's own 1D I(q) shape (RMS(log) ~0.02 for a right centre
    against ~0.22 for the naive one).
    """
    return float(s.get("Xcenter")), image_shape[0] - 1 - float(s.get("Ycenter"))


def q_of_r(r_px, s):
    """q (1/A) at a pixel radius, using the same geometry RAW uses."""
    px_mm = float(s.get("DetectorPixelSizeX")) / 1000.0
    dist = float(s.get("SampleDistance"))
    lam = float(s.get("WaveLength"))
    return 4 * np.pi * np.sin(np.arctan(np.asarray(r_px, float) * px_mm / dist) / 2) / lam


def r_of_q(q, s):
    """inverse of q_of_r (bilinear search on a 0.1 px grid)."""
    rr = np.arange(0, 3000, 0.1)
    return float(rr[int(np.argmin(abs(q_of_r(rr, s) - q)))])


def unmasked_fraction_by_radius(s, rmax=80, dr=2):
    """radial profile of the beam-stop mask: how much of each ring is usable."""
    mask = np.asarray(s.get("Masks")["BeamStopMask"][0]) > 0
    cx, cy = beam_center(mask.shape, s)
    yy, xx = np.indices(mask.shape)
    r = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)
    out = []
    for r0 in np.arange(0, rmax, dr):
        m = (r >= r0) & (r < r0 + dr)
        out.append((float(r0), float(mask[m].mean())))
    return out
