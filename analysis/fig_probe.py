"""
@file fig_probe.py
@brief "Ptychography fitted the phase": the true, start and fitted probe phase across the aperture, and what is left.

Reads analysis/probe_pack.py's files (one .npz per released run, the probes in the aperture plane). For each run: the true probe's
phase, the start's (the search's errors), the fitted probe's, and fitted minus true. The fitted probe is free in a real-space shift
and a global phase (the object takes both up), so before differencing it is moved back by the shift that maximises its overlap with
the truth and its mean phase against the truth is removed. Phases are shown inside the aperture only (|true| above half its peak).

    ~/hyperspy-bundle/bin/python analysis/fig_probe.py <dir of npz> [name-pattern ...] [--out <file stem>]
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def align(F, T, n):
    """F moved by the real-space shift that best matches T (found on the full n x n grid), then its mean phase against T removed"""
    h = F.shape[0] // 2
    pad = np.zeros((n, n), complex)
    c = n // 2
    pad[c - h:c + h, c - h:c + h] = F * np.conj(T)
    xc = np.fft.ifft2(np.fft.ifftshift(pad))                 # the cross-correlation of the two probes in real space
    iy, ix = np.unravel_index(np.argmax(np.abs(xc)), xc.shape)
    sy = iy if iy < n // 2 else iy - n
    sx = ix if ix < n // 2 else ix - n
    ky, kx = np.meshgrid(np.arange(-h, h), np.arange(-h, h), indexing="ij")
    G = F * np.exp(+2j * np.pi * (ky * sy + kx * sx) / n)     # undo the shift (s) the correlation found (checked on a rolled truth)
    G = G * np.exp(-1j * np.angle(np.sum(G * np.conj(T))))
    return G, (sy, sx)


def panel(axs, z, title):
    T, F = z["true_k"].astype(complex), z["fitted_k"].astype(complex)
    S = z["start_k"].astype(complex) if "start_k" in z.files else None
    n, da = int(z["n"]), float(z["d_alpha_mrad"])
    m = np.abs(T) > 0.5 * np.abs(T).max()
    G, shift = align(F, T, n)
    h = T.shape[0] // 2
    ext = [-h * da, h * da, -h * da, h * da]
    if S is not None:
        S, _ = align(S, T, n)
    cols = (("true", T), ("start", S), ("fitted", G), ("fitted - true", G * np.conj(T) if G is not None else None))
    for ax, (lab, A) in zip(axs, cols):
        ax.set_xticks([]); ax.set_yticks([])
        if A is None:
            ax.text(0.5, 0.5, "no start\nwritten", ha="center", va="center", transform=ax.transAxes, fontsize=6); continue
        ph = np.where(m, np.angle(A), np.nan)
        ax.imshow(ph, extent=ext, origin="lower", cmap="twilight", vmin=-np.pi, vmax=np.pi, interpolation="nearest")
        ax.set_title(lab, fontsize=7, pad=2)
        r = float(z["convergence_mrad"]) * 1.06
        ax.set_xlim(-r, r); ax.set_ylim(-r, r)
    d = np.angle(G * np.conj(T))[m]
    axs[0].set_ylabel(title, fontsize=7)
    ov = f"overlap {float(z['ov_fitted']):.3f}" + (f" (start {float(z['ov_start']):.3f})" if "ov_start" in z.files else "")
    axs[-1].text(1.04, 0.5, f"rms {np.sqrt(np.mean(d ** 2)):.2f} rad\n{ov}\nshift {shift} px", transform=axs[-1].transAxes,
                 fontsize=6, va="center")
    return float(np.sqrt(np.mean(d ** 2)))


def main(*args):
    a = list(args)
    out = None
    if "--out" in a:
        i = a.index("--out"); out = a[i + 1]; del a[i:i + 2]
    d, pats = a[0], (a[1:] or ["*"])
    files = sorted({f for p in pats for f in glob.glob(os.path.join(d, f"{p}.npz"))})
    if not files:
        raise SystemExit(f"no npz in {d} for {pats}")
    fig, axes = plt.subplots(len(files), 4, figsize=(6.0, 1.55 * len(files)), squeeze=False)
    for row, f in zip(axes, files):
        name = os.path.basename(f)[:-4].replace("recon_", "")
        rms = panel(row, np.load(f), name[:28])
        print(f"{name}: rms fitted - true {rms:.3f} rad")
    fig.tight_layout()
    out = out or os.path.join(d, "fig_probe")
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=250, bbox_inches="tight")
    print(f"wrote {out}.pdf and .png")


if __name__ == "__main__":
    main(*sys.argv[1:])
