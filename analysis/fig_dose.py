"""
@file fig_dose.py
@brief Figure 2 of the single-atom paper (DRAFT, 2026-11-22): what dose does to one atom's depth extent, three tracks from Figure 1.

(a) The continuous W80 (analysis/flat_plot.py width_cont) against dose for the perfect 30 mrad probe, the 80 mrad aperture flat only to
30, and the perfect 80 mrad probe, with the true probe held (`gflat_dosefix`, `gflat_d1e6fix`), on the simulation's 712 px patterns
(914 at flat to 30; solid) and summed 4 x 4 into a camera (CAMERA_BIN=4, hollow, dashed); "no noise" at the right. (b) Their atoms
(x-z sections through the peak, each scaled to its own peak) at each dose. Draft: the released, fitted probe under noise and the
choice between fine pixels with the Poisson model and a coarser camera are still to come (P17, P20); a missing point is an empty panel.

    ~/hyperspy-bundle/bin/python analysis/fig_dose.py [out_dir]
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "analysis"))
import flat_plot as fp  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

PACKS = [os.path.expanduser(f"~/Desktop/zblur/{d}") for d in ("gflat", "g0d", "g1")]
TRACKS = (("perfect 30", "A", 30, "#1f5fa8", "--"), ("flat 30, 80", "B", 30, "#c8501e", "-"), ("perfect 80", "A", 80, "#1f5fa8", "-"))
DOSES = (1e6, 1e7, 1e8, None)                     # None = no noise
XPOS = {1e6: 1e6, 1e7: 1e7, 1e8: 1e8, None: 1e9}


def find(rs, line, x, dose, cam=1, tag=""):
    c = sorted([r for r in rs if r["geom"] == "single" and r["line"] == line and r["x"] == x and r["tag"] == tag and r["dose"] == dose
                and r["camera"] == cam and r["fs"] == 0], key=lambda r: r["early"])
    return c[0] if c else None


def kernel_file(r):
    for pk in PACKS:
        f = os.path.join(pk, "psf", f"psf_{r['kernel']}_vol.npy")
        if os.path.exists(f):
            return f


def section(ax, r, zhalf=8.0, xhalf=1.5, dx=0.04922, dz=1.0):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor("k")
    if r is None:
        ax.set_facecolor("0.92"); ax.text(0.5, 0.5, "not run", ha="center", va="center", transform=ax.transAxes, fontsize=6, color="0.4")
        return
    K = fp.blobs.zp.load_phase(kernel_file(r))
    col = fp.pg.core_column(K, dx); k0 = int(np.argmax(col))
    _, ky, kx = np.unravel_index(np.argmax(K), K.shape)
    xz = K[:, ky, :]
    z = (np.arange(K.shape[0]) - k0) * dz; x = (np.arange(K.shape[2]) - kx) * dx; m = np.abs(x) <= xhalf
    ax.imshow(xz[:, m] / xz.max(), extent=[x[m][0] - dx / 2, x[m][-1] + dx / 2, z[-1] + dz / 2, z[0] - dz / 2], aspect="auto",
              cmap="inferno", vmin=0, vmax=1, interpolation="nearest")
    ax.set_ylim(zhalf, -zhalf); ax.set_xlim(-xhalf, xhalf)
    ax.text(0.04, 0.96, f"{r['w80c']:.1f} Å", transform=ax.transAxes, color="w", fontsize=5.5, va="top")


def main(*args):
    out_dir = args[0] if args else os.path.join(os.path.dirname(HERE), "ptychoshelves-clean", "report", "aberrations_paper", "figs")
    rs = fp.rows([p for p in PACKS if os.path.isdir(p)])
    fig = plt.figure(figsize=(7.2, 4.4))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1, 1.15], wspace=0.22)
    ax = fig.add_subplot(gs[0, 0])
    for lab, line, x, col, ls in TRACKS:
        for cam, mfc, lsc in ((1, col, ls), (4, "white", ":")):
            pts = [(XPOS[d], r["w80c"]) for d in DOSES if (r := find(rs, line, x, d, cam)) is not None]
            if len(pts) < 1:
                continue
            ax.plot([p[0] for p in pts], [p[1] for p in pts], lsc, color=col, lw=1.3, marker="o", ms=4, mfc=mfc,
                    label=lab + (" (712 px)" if cam == 1 else ", camera"))
    ax.set_xscale("log"); ax.set_xlim(5e5, 2e9)
    ax.set_xticks([1e6, 1e7, 1e8, 1e9]); ax.set_xticklabels(["$10^6$", "$10^7$", "$10^8$", "no\nnoise"])
    ax.set_ylim(0, 24); ax.set_xlabel("dose (e/Å$^2$)"); ax.set_ylabel("$W_{80}$ (Å)")
    ax.set_title("a   depth extent against dose, probe known (draft)", loc="left", fontweight="bold", fontsize=8)
    ax.legend(loc="upper right", frameon=False, fontsize=5.8, ncol=1)
    sub = gs[0, 1].subgridspec(3, 4, wspace=0.05, hspace=0.12)
    for i, (lab, line, x, col, ls) in enumerate(TRACKS):
        for j, d in enumerate(DOSES):
            axs = fig.add_subplot(sub[i, j])
            section(axs, find(rs, line, x, d))
            if i == 0:
                axs.set_title("no noise" if d is None else f"{d:.0e}".replace("e+0", "e"), fontsize=6.5, pad=2)
            if j == 0:
                axs.set_ylabel(lab, fontsize=6.5)
    fig.text(0.53, 0.93, "b   the atoms, probe known, 712 px (z +-8 Å, x +-1.5 Å)", fontweight="bold", fontsize=8)
    os.makedirs(out_dir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out_dir, f"fig_dose_draft.{ext}"), dpi=300, bbox_inches="tight")
    print(f"wrote {os.path.join(out_dir, 'fig_dose_draft')}.pdf and .png")


if __name__ == "__main__":
    main(*sys.argv[1:])
