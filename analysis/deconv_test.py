"""
@file deconv_test.py
@brief Can the blur be deconvolved? Each blurred lone-atom kernel against the unblurred one convolved with the blur the
physics predicts, and the atom's total phase before and after.

Deconvolution assumes  reconstruction = clean atom (*) blur kernel.  The physics fixes the kernel for each source:
  source size   sideways Gaussian, sigma = FWHM / 2.355
  focal spread  Gaussian along the beam, sigma = the focal spread (rms)
  vibration     isotropic 3-D Gaussian, sigma = Pb's per-axis rms displacement at T, from the simulator's own Debye
                model (sim/simulate_4dstem.py, B_at_temperature): 0.0385 / 0.0635 / 0.1068 A at 20 / 100 / 295 K
  combinations  the sources' Gaussians compounded (variances add per axis)
For every blurred kernel A with base B (same probe, alpha, scan, box) and prediction P = B (*) kernel, inside a cylinder
r <= R_CMP around the atom through the whole box (background = each layer's median, as everywhere in zblur):
  resid_none   |A - B| / |A|: doing nothing
  resid_conv   |A - P| / |A|: the convolution model, as the physics gives it
  scale, resid_shape   the best single factor a on P, and |A - aP| / |A|: shape alone, contrast forgiven
  kept         sum(A) / sum(B) within r <= R_SUM: the atom's total phase kept. A convolution with a normalised kernel
               keeps it (1, up to what spills out of the box). Meaningless on the ARM base, whose ringing sums to ~0
  kept_pred    the same for P (the spill a convolution would have)
  resid_proj   |proj A - proj P| / |proj A|, proj = the sum along the beam (the 2-D projected image) inside r <= R_CMP:
               the part a sideways deconvolution could act on
  peak_A / peak_P, fwhm_A / fwhm_P   the centre column's peak (rad) and FWHM (A), reconstruction and prediction
  amplify      the share of the base's Fourier support (|T_base| > 5 % of its max) where |T_A| > 1.2 |T_base|: frequencies
               a blur would have to amplify. Any blur that averages over positions (source points, focus values, thermal
               displacements) has a kernel that is a probability distribution, |H(k)| <= 1 at every k, so this is zero
               for a real convolution of any shape
P is aligned to A before the residuals (whole-slice shifts along the beam within +-1 slice, a sub-pixel sideways shift
from the projections' centroids), so a residual is never an offset: the test is as kind to the convolution as it can be.

    ~/hyperspy-bundle/bin/python analysis/deconv_test.py ~/Desktop/zblur/g0d ~/Desktop/zblur/g1 ~/Desktop/zblur/g1b
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter, shift as nd_shift

sys.path.insert(0, os.path.dirname(__file__))
import blobs  # noqa: E402
import zprofile as zp  # noqa: E402

PB_U = {20: 0.0385, 100: 0.0635, 295: 0.1068}     # A, per-axis rms (B / 8 pi^2)^0.5 from the simulator's Debye model
R_CMP = 1.5                                       # A, the comparison cylinder
R_SUM = 2.0                                       # A, the cylinder for the total phase


def bg(K):
    return K - np.median(K, axis=(1, 2), keepdims=True)


def kernel_sigmas(r):
    """(sigma along the beam, sigma sideways) in A of the convolution the physics predicts for r's blur sources."""
    sz, sxy = r["fs"], r["src"] / 2.3548
    if r["T"] is not None:
        u = PB_U[int(round(r["T"]))]
        sz, sxy = np.hypot(sz, u), np.hypot(sxy, u)
    return sz, sxy


def predict(b, r):
    """The unblurred kernel convolved with r's predicted blur (vacuum, zero, beyond the box)."""
    sz, sxy = kernel_sigmas(r)
    return gaussian_filter(bg(b["K"]), sigma=(sz / b["dz"], sxy / b["dx"], sxy / b["dx"]), mode="constant", cval=0.0)


def cylinder(K, dx, R):
    nz, ny, nx = K.shape
    yy, xx = np.mgrid[:ny, :nx]
    return np.broadcast_to(np.hypot(yy - ny // 2, xx - nx // 2) * dx <= R, K.shape)


def centroid(img, dx):
    ny, nx = img.shape
    yy, xx = np.mgrid[:ny, :nx]
    w = np.clip(img, 0, None) * (np.hypot(yy - ny // 2, xx - nx // 2) * dx <= 0.5)
    return float((w * yy).sum() / w.sum()), float((w * xx).sum() / w.sum())


def align(P, A, dx):
    """P moved onto A: the sideways offset of the projections' centroids, then the whole-slice shift along the beam
    (within +-1 slice) that fits A best."""
    cyA, cxA = centroid(A.sum(0), dx)
    cyP, cxP = centroid(P.sum(0), dx)
    P = nd_shift(P, (0, cyA - cyP, cxA - cxP), order=1, mode="constant", cval=0.0)
    m = cylinder(A, dx, R_CMP)
    best = None
    for k in (-1, 0, 1):
        Q = np.roll(P, k, axis=0)
        if k > 0:
            Q[:k] = 0
        elif k < 0:
            Q[k:] = 0
        s = float(np.dot(A[m], Q[m]) / np.dot(Q[m], Q[m]))
        e = float(np.linalg.norm(A[m] - s * Q[m]))
        if best is None or e < best[0]:
            best = (e, Q)
    return best[1]


def centre(K, r):
    m = zp.psf3d(K, r["dz"], z0=r["z0"], dx=r["dx"], core_A=3.0)
    return float(m["peak"]), float(m["fwhm_z"])


def compare(r, b):
    A, B = bg(r["K"]), bg(b["K"])
    P0 = predict(b, r)
    P = align(P0, A, r["dx"])
    m = cylinder(A, r["dx"], R_CMP)
    a, bb, p = A[m], align(B, A, r["dx"])[m], P[m]
    na = np.linalg.norm(a)
    s = float(np.dot(a, p) / np.dot(p, p))
    ms = cylinder(A, r["dx"], R_SUM)
    m2 = m[0]
    pa, pp = A.sum(0)[m2], P.sum(0)[m2]
    sz, sxy = kernel_sigmas(r)
    E, _, _, ok = blobs.envelope(r, b)
    pkA, fwA = centre(r["K"], r)
    pkP, fwP = centre(P0 + np.median(b["K"], axis=(1, 2), keepdims=True), r)
    return dict(sigma_z=sz, sigma_xy=sxy, resid_none=float(np.linalg.norm(a - bb) / na),
                resid_conv=float(np.linalg.norm(a - p) / na), scale=s, resid_shape=float(np.linalg.norm(a - s * p) / na),
                kept=float(A[ms].sum() / B[ms].sum()), kept_pred=float(P0[ms].sum() / B[ms].sum()),
                resid_proj=float(np.linalg.norm(pa - pp) / np.linalg.norm(pa)), peak_A=pkA, peak_P=pkP, fwhm_A=fwA,
                fwhm_P=fwP, amplify=float(np.mean(np.nan_to_num(E[ok]) > 1.2)))


def main(*packs):
    rs = blobs.load(packs)
    bases = {blobs.base_key(r): r for r in rs if not r["fs"] and not r["src"] and r["T"] is None}
    rows = []
    for r in rs:
        b = bases.get(blobs.base_key(r))
        if b is None or b is r:
            continue
        c = blobs.core(r)
        clean = bool(np.isfinite(c["fwhm_z"]) and c["dip"] > -0.9)
        rows.append(dict(name=blobs.describe(r), clean=clean, kernel=r["kernel"], **compare(r, b)))
    tag = "_".join(os.path.basename(os.path.normpath(p)) for p in packs)
    out = os.path.join(blobs.HERE, "results", f"deconv_{tag}.tsv")
    cols = ["name", "clean", "sigma_z", "sigma_xy", "resid_none", "resid_conv", "scale", "resid_shape", "kept",
            "kept_pred", "resid_proj", "peak_A", "peak_P", "fwhm_A", "fwhm_P", "amplify", "kernel"]
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for row in rows:
            w.writerow([f"{row[k]:.4f}" if isinstance(row[k], float) else row[k] for k in cols])
    print(f"{'kernel':40s} {'none':>5s} {'conv':>5s} {'scale':>5s} {'shape':>5s} {'proj':>5s} {'kept':>5s} "
          f"{'pkA':>5s} {'pkP':>5s} {'fwA':>5s} {'fwP':>5s} {'ampl':>5s}")
    for row in rows:
        print(f"{row['name'] + ('' if row['clean'] else ' *'):40s} {row['resid_none']:5.2f} {row['resid_conv']:5.2f} "
              f"{row['scale']:5.2f} {row['resid_shape']:5.2f} {row['resid_proj']:5.2f} {row['kept']:5.2f} "
              f"{row['peak_A']:5.2f} {row['peak_P']:5.2f} {row['fwhm_A']:5.2f} {row['fwhm_P']:5.2f} {row['amplify']:5.2f}")
    print("wrote", out)


if __name__ == "__main__":
    main(*(sys.argv[1:] or [os.path.expanduser(f"~/Desktop/zblur/{p}") for p in ("g0d", "g1", "g1b")]))
