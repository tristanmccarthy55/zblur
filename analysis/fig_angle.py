"""
@file fig_angle.py
@brief Figure 1 of the single-atom paper (draft, 2026-11-19): what opening past the flat region does to one atom, probe released.

Top two rows: an 80 mrad aperture over a corrector phase-flat only to x (x = 30, 40, 60; 80 is the perfect probe), each with the
phase across its aperture (wrapped, the flat disc dashed) above the atom it reconstructs (x-z section through the atom's peak).
Next two rows: the perfect probe with an aperture of x, the same way. Bottom: the depth extent against x for both lines, as FWHM
(c) and as W80, the depth holding 80 % of the atom's phase (d), the probe released (solid) and known (dotted). Every reconstruction
here starts from the crystal search's C1 and C3 errors with A1 right (`_pu_noa1`, P9/P10); no noise. Phases from the simulated
tableaux (runs/flat_tableaux.tsv, runs/window_loss.py chi, checked against abTEM); kernels from ~/Desktop/zblur/gflat.

    ~/hyperspy-bundle/bin/python analysis/fig_angle.py [out_dir]
"""
from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "analysis"))
sys.path.insert(0, os.path.join(HERE, "runs"))
import fig_flat as ff  # noqa: E402
import fig_probe as fpr  # noqa: E402
import flat_plot as fp  # noqa: E402
import window_loss as wl  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

PSF = os.path.expanduser("~/Desktop/zblur/gflat/psf")
TAG = "_afb"                                  # the C1 + A1 fit at gain 100 from the full search errors (P15)
PROBES = os.path.expanduser("~/Desktop/zblur/probes/probes_afb")
XS = (30, 40, 60, 80)
DX, DZ = 0.04922, 1.0
FLAT = {30: "Pb_a080_flat30_single_st0p1_w8_s0p25_bz40p000_rs90_NL40", 40: "Pb_a080_flat40_single_st0p1_w8_s0p25_bz40p000_NL40",
        60: "Pb_a080_flat60_single_st0p1_w8_s0p25_bz20p000_NL20", 80: "Pb_a080_ideal_df0_single_st0p1_w8_s0p25_bz20p000_NL20"}
PERF = {30: "Pb_a030_ideal_df0_single_st0p1_w8_s0p25_bz40p000_NL40", 40: "Pb_a040_ideal_df0_single_st0p1_w8_s0p25_bz40p000_NL40",
        60: "Pb_a060_ideal_df0_single_st0p1_w8_s0p25_bz20p000_NL20", 80: "Pb_a080_ideal_df0_single_st0p1_w8_s0p25_bz20p000_NL20"}


def tableaux():
    out = {}
    for r in csv.DictReader((l for l in open(os.path.join(HERE, "runs", "flat_tableaux.tsv")) if not l.startswith("#")), delimiter="\t"):
        out[int(r["label"][-3:])] = (float(r["c1"]), json.loads(r["aber_json"]))
    return out


def phase_map(ax, x, aperture, tab, lim=85.0, n=301):
    """the probe's phase across its aperture (mrad axes), wrapped, piston of the flat disc removed; the flat disc dashed"""
    t = np.linspace(-lim, lim, n)
    TX, TY = np.meshgrid(t, t)
    k = 1e-3 / wl.LAM
    inside = np.hypot(TX, TY) <= aperture
    if x in tab and aperture == 80:
        c1, ab = tab[x]
        ph = wl.chi(TX * k, TY * k, c1, ab)
        ph = ph - ph[np.hypot(TX, TY) <= x].mean()
    else:
        ph = np.zeros_like(TX)
    ph = np.angle(np.exp(1j * ph))
    ax.imshow(np.where(inside, ph, np.nan), extent=[-lim, lim, -lim, lim], origin="lower", cmap="twilight", vmin=-np.pi, vmax=np.pi,
              interpolation="nearest")
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(aperture * np.cos(th), aperture * np.sin(th), color="0.3", lw=0.6)
    if aperture == 80 and x < 80:
        ax.plot(x * np.cos(th), x * np.sin(th), "--", color="k", lw=0.8)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")


def atom(ax, kernel, zhalf=6.0, xhalf=1.5):
    K = fp.blobs.zp.load_phase(os.path.join(PSF, f"psf_{kernel}{TAG}_vol.npy"))
    col = fp.pg.core_column(K, DX)
    k0 = int(np.argmax(col))
    _, ky, kx = np.unravel_index(np.argmax(K), K.shape)
    xz = K[:, ky, :]
    z = (np.arange(K.shape[0]) - k0) * DZ
    xx = (np.arange(K.shape[2]) - kx) * DX
    m = np.abs(xx) <= xhalf
    ax.imshow(xz[:, m] / xz.max(), extent=[xx[m][0] - DX / 2, xx[m][-1] + DX / 2, z[-1] + DZ / 2, z[0] - DZ / 2], aspect="auto",
              cmap="inferno", vmin=0, vmax=1, interpolation="nearest")
    ax.set_ylim(zhalf, -zhalf); ax.set_xlim(-xhalf, xhalf); ax.set_facecolor("k")   # past the box's end: black, not white
    ax.text(0.04, 0.97, f"{xz.max():.2f} rad", transform=ax.transAxes, color="w", fontsize=6, va="top")


def main(*args):
    out_dir = args[0] if args else os.path.join(os.path.dirname(HERE), "ptychoshelves-clean", "report", "aberrations_paper", "figs")
    tab = tableaux()
    rs = fp.rows([os.path.expanduser(f"~/Desktop/zblur/{d}") for d in ("gflat", "g0d", "g1")])
    fig = plt.figure(figsize=(7.2, 9.6))
    gs = GridSpec(6, 4, figure=fig, height_ratios=[0.85, 1.45, 0.85, 1.45, 0.25, 2.0], hspace=0.42, wspace=0.12)
    for row, (label, kernels, aperture_of) in enumerate(((("a", "80 mrad aperture, phase flat only to x"), FLAT, lambda x: 80),
                                                         (("b", "perfect probe, aperture x"), PERF, lambda x: x))):
        for i, x in enumerate(XS):
            axp = fig.add_subplot(gs[2 * row, i])
            phase_map(axp, x, aperture_of(x), tab)
            axp.set_title(f"x = {x} mrad", fontsize=7, pad=2)
            axa = fig.add_subplot(gs[2 * row + 1, i])
            atom(axa, kernels[x])
            if row == 1:
                axa.set_xlabel("x (Å)", fontsize=7, labelpad=1)
            else:
                axa.set_xticklabels([])
            if i == 0:
                axp.text(-0.45, 1.32, label[0], transform=axp.transAxes, fontweight="bold", fontsize=8)
                axp.text(-0.30, 1.32, label[1] + " (phase across the aperture; the atom below)", transform=axp.transAxes, fontsize=7.5)
                axa.set_ylabel("z (Å)", fontsize=7)
            else:
                axa.set_yticklabels([])
            axa.tick_params(labelsize=6)
    sub = gs[5, :].subgridspec(1, 2, wspace=0.25, width_ratios=[1, 1.1])
    ax = fig.add_subplot(sub[0, 0])
    a0, b0 = ff.pick(rs, "", 0.0)
    ff.lines(ax, a0, b0, "w80c", faint=True)
    a, b = ff.pick(rs, TAG, 0.0)
    ff.lines(ax, a, b, "w80c")
    ax.plot([], [], ":", color="0.4", lw=1.2, label="probe known")
    ax.legend(loc="upper right", frameon=False, fontsize=6.5)
    ax.set_xlim(25, 85); ax.set_ylim(0, 12)
    ax.set_xlabel("x (mrad)"); ax.set_ylabel("$W_{80}$ (Å)")
    ax.set_title("c   depth holding 80 % of the phase", loc="left", fontweight="bold", fontsize=8)
    # (d) the probe the fit found: start (the search's errors), fitted, fitted minus true, across the aperture
    pr = sub[0, 1].subgridspec(2, 3, wspace=0.06, hspace=0.25)
    for i, (lab, kern) in enumerate((("perfect 80", PERF[80]), ("flat 30, 80", FLAT[30]))):
        z = np.load(os.path.join(PROBES, f"recon_{kern}{TAG}.npz"))
        T, F, S = (z[k].astype(complex) for k in ("true_k", "fitted_k", "start_k"))
        n, da = int(z["n"]), float(z["d_alpha_mrad"])
        m = np.abs(T) > 0.5 * np.abs(T).max()
        G, _ = fpr.align(F, T, n); Sa, _ = fpr.align(S, T, n)
        h = T.shape[0] // 2
        for j, (title, A) in enumerate((("start", Sa * np.conj(T)), ("fitted", G * np.conj(T)), ("fitted, defocus removed", None))):
            axq = fig.add_subplot(pr[i, j])
            if A is None:                                  # the fitted error with its defocus part removed: what is left besides C1
                ky, kx = np.mgrid[-h:h, -h:h] * da / 1000
                d = np.angle(G * np.conj(T))
                th2 = (kx ** 2 + ky ** 2)[m]
                M = np.c_[np.ones_like(th2), th2]
                c = np.linalg.lstsq(M, d[m], rcond=None)[0]
                A = np.exp(1j * (d - c[0] - c[1] * (kx ** 2 + ky ** 2)))
            axq.imshow(np.where(m, np.angle(A), np.nan), extent=[-h * da, h * da, -h * da, h * da], origin="lower", cmap="twilight",
                       vmin=-np.pi, vmax=np.pi, interpolation="nearest")
            r_ = float(z["convergence_mrad"]) * 1.05
            axq.set_xlim(-r_, r_); axq.set_ylim(-r_, r_); axq.set_xticks([]); axq.set_yticks([])
            if i == 0:
                axq.set_title(title, fontsize=6.5, pad=2)
            if j == 0:
                axq.set_ylabel(lab, fontsize=6.5)
    fig.text(0.535, 0.305, "d   phase error of the probe, start and fitted", fontweight="bold", fontsize=8)
    os.makedirs(out_dir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out_dir, f"fig_angle.{ext}"), dpi=300, bbox_inches="tight")
    print(f"wrote {os.path.join(out_dir, 'fig_angle')}.pdf and .png")


if __name__ == "__main__":
    main(*sys.argv[1:])
