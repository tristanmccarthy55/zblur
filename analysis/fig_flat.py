"""
@file fig_flat.py
@brief The flat experiment's figure for the write-up (2026-11-13): one lone Pb atom's depth extent against angle.

(a) W80 against x, the true probe fixed: the perfect probe with an aperture of x mrad, and an 80 mrad aperture flat only to x.
(b) x-z sections of three of those atoms. Held back until their runs are complete (--found, --extra): the probe released from a
quarter of the misses, dashed in (a); (c) the two lines with 2.5 A of focal spread; (d) z-share against how far the start was.
Numbers come from analysis/flat_plot.py's rows (W80, z-share); kernels from the packs in ~/Desktop/zblur (gflat, g0d, g1).

    ~/hyperspy-bundle/bin/python analysis/fig_flat.py [out_dir]               # fig_flat.pdf / .png: a and b
    ~/hyperspy-bundle/bin/python analysis/fig_flat.py --found --extra         # fig_flat_found_extra: all four, for us
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import flat_plot as fp  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402

PACKS = [os.path.expanduser(f"~/Desktop/zblur/{d}") for d in ("gflat", "g0d", "g1")]
BLUE, ORANGE, GREY = "#1f5fa8", "#c8501e", "0.55"
START = {"_pu0": 0.0, "_pu025": 0.25, "_pu050": 0.5, "_pu": 1.0}       # fraction of the crystal search's misses (A1 0.6 waves at 80 mrad)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "font.family": "sans-serif"})


def pick(rs, tag, fs):
    """the line-A and line-B points of one treatment; the earlier runs fill line A where this treatment has none."""
    pts = [r for r in rs if r["geom"] == "single" and r["tag"] == tag and r["fs"] == fs and r["dose"] is None]
    have = {(r["line"], r["x"]) for r in pts if not r["early"]}
    pts = [r for r in pts if not r["early"] or (r["line"], r["x"]) not in have]
    a = sorted([r for r in pts if r["line"] == "A" and r["x"] <= 80], key=lambda r: r["x"])
    b = sorted([r for r in pts if r["line"] == "B"], key=lambda r: r["x"])
    a80 = [r for r in a if r["x"] == 80]
    if a80:
        b = sorted(b + [a80[0]], key=lambda r: r["x"])
    return a, b


def lines(ax, a, b, key, faint=False, labels=True, join_b=True, dashed=False):
    for pts, col, lab, join in ((a, BLUE, "perfect probe, aperture $x$", True), (b, ORANGE, "80 mrad aperture, flat only to $x$", join_b)):
        if not pts:
            continue
        xs, ys = [p["x"] for p in pts], [p[key] for p in pts]
        if faint:
            ax.plot(xs, ys, ":", color=col, lw=1.2, alpha=0.45)
        elif dashed:
            ax.plot(xs, ys, "--", color=col, lw=1.1, marker="s", ms=3,
                    label="released from a quarter of the misses" if (labels and col == BLUE) else None)
        else:
            ax.plot(xs, ys, linestyle="-" if join else "none", color=col, lw=1.4, marker="o", ms=4, label=lab if labels else None)


def section(ax, file, title, dz, dx, zhalf=6.0, xhalf=1.5):
    """x-z section through the atom's own peak (atoms sit up to a pixel off the crop centre: the G0d 80 mrad one peaks one row
    off, and a cut through the centre row reads its shoulder, 1.70 rad, not its peak, 1.86). Returns its depth profile."""
    K = fp.blobs.zp.load_phase(file)
    col = fp.pg.core_column(K, dx)
    k0 = int(np.argmax(col))
    nz, ny, nx = K.shape
    _, ky, kx = np.unravel_index(np.argmax(K), K.shape)
    xz = K[:, ky, :]
    z = (np.arange(nz) - k0) * dz
    x = (np.arange(nx) - kx) * dx
    m = np.abs(x) <= xhalf
    ax.imshow(xz[:, m] / xz.max(), extent=[x[m][0] - dx / 2, x[m][-1] + dx / 2, z[-1] + dz / 2, z[0] - dz / 2], aspect="auto",
              cmap="inferno", vmin=0, vmax=1, interpolation="nearest")
    ax.set_ylim(zhalf, -zhalf)
    ax.set_title(f"{title}\n{xz.max():.2f} rad", fontsize=7)
    ax.set_xlabel("x (Å)")
    return z, col


def profile(ax, z, col, dz, **kw):
    """the depth profile drawn as the slices it is made of: one flat step per 1 A reconstruction slice"""
    zs = np.repeat(z, 2) + np.tile([-dz / 2, dz / 2], len(z))
    ax.plot(np.repeat(col, 2), zs, **kw)


def main(*args):
    """a and b by default: what is solid. --found adds the probe released from a quarter of the misses to a (dashed); --extra adds
    c (focal spread) and d (the start ladder), held back from the write-up until those runs are complete."""
    found, extra = "--found" in args, "--extra" in args
    paths = [a for a in args if not a.startswith("--")]
    out_dir = paths[0] if paths else os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                                                  "ptychoshelves-clean", "report", "aberrations_paper", "figs")
    rs = fp.rows([p for p in PACKS if os.path.isdir(p)])
    fig = plt.figure(figsize=(7.2, 5.4 if extra else 2.9))
    gs = GridSpec(2 if extra else 1, 2, figure=fig, width_ratios=[1, 1.25], hspace=0.55, wspace=0.28)

    # (a) the graph: true probe fixed (and, with --found, released from a quarter of the misses)
    ax = fig.add_subplot(gs[0, 0])
    a, b = pick(rs, "", 0.0)
    lines(ax, a, b, "w80")
    if found:
        aq, bq = pick(rs, "_pu025", 0.0)
        lines(ax, aq, bq, "w80", dashed=True)
    ax.axhspan(2.6, 3.4, color="0.88", zorder=0, label="perfect probe at 80-100 mrad")
    ax.set_xlim(25, 85); ax.set_ylim(0, 16 if found else 12); ax.set_yticks(range(0, 13, 2) if found else range(0, 11, 2))
    ax.set_xlabel("$x$ (mrad)"); ax.set_ylabel("$W_{80}$ (Å)")
    ax.set_title("a   probe known" + (" (solid), probe found (dashed)" if found else ""), loc="left", fontweight="bold")
    ax.legend(loc="upper right", frameon=False)

    # (b) x-z sections: the perfect 30 mrad probe, 80 mrad flat only to 30, the perfect 80 mrad probe; then their depth profiles
    # on one axis in rad, so strength and length compare directly (each section is scaled to its own peak)
    sub = GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[0, 1], wspace=0.15, width_ratios=[1, 1, 1, 1.35])
    files = {"perfect 30": ("gflat", "psf_Pb_a030_ideal_df0_single_st0p1_w8_s0p25_bz40p000_NL40_vol.npy", BLUE, "--"),
             "flat 30, 80": ("gflat", "psf_Pb_a080_flat30_single_st0p1_w8_s0p25_bz40p000_rs90_NL40_vol.npy", ORANGE, "-"),
             "perfect 80": ("g0d", "psf_Pb_a080_ideal_df0_single_st0p1_w8_s0p25_bz20p000_NL20_vol.npy", BLUE, "-")}
    axp = fig.add_subplot(sub[0, 3])
    for i, (title, (pack, f, col_, ls)) in enumerate(files.items()):
        axs = fig.add_subplot(sub[0, i])
        z, col = section(axs, os.path.expanduser(f"~/Desktop/zblur/{pack}/psf/{f}"), title, dz=1.0, dx=0.04922)
        profile(axp, z, col, 1.0, color=col_, ls=ls, lw=1.2, label=title)
        if i == 0:
            axs.set_ylabel("z along the beam (Å)")
            axs.text(-0.55, 1.22, "b", transform=axs.transAxes, fontweight="bold", fontsize=8)
        else:
            axs.set_yticklabels([])
    axp.axvline(0, color="0.7", lw=0.6)
    axp.set_ylim(6, -6); axp.set_yticklabels([])
    axp.set_xlim(-0.05, 0.62); axp.set_xticks([0, 0.25, 0.5])
    axp.set_xlabel("phase (rad)")
    axp.set_title("along the beam\n(core, 0.15 Å)", fontsize=7)
    axp.legend(loc="upper right", frameon=False, fontsize=6, handlelength=1.6)

    if extra:
        # (c) focal spread 2.5 A, probe known (the line without it faint)
        ax = fig.add_subplot(gs[1, 0])
        lines(ax, a, b, "w80", faint=True)
        a2, b2 = pick(rs, "", 2.5)
        lines(ax, a2, b2, "w80", join_b=False)
        ax.set_xlim(25, 85); ax.set_ylim(0, 11)
        ax.set_xlabel("$x$ (mrad)"); ax.set_ylabel("$W_{80}$ (Å)")
        ax.set_title("c   probe known, focal spread 2.5 Å", loc="left", fontweight="bold")
        # (d) the probe released, against how far the start was from the truth
        ax = fig.add_subplot(gs[1, 1])
        for line, x, col, lab in (("A", 80, BLUE, "perfect probe, 80 mrad"), ("B", 30, ORANGE, "80 mrad aperture, flat to 30")):
            pts = sorted([(START[r["tag"]], r["zshare"], r["clean"]) for r in rs if r["geom"] == "single" and r["fs"] == 0
                          and r["tag"] in START and r["line"] == line and r["x"] == x])
            if pts:
                ax.plot([p[0] for p in pts], [p[1] for p in pts], "-", color=col, lw=1.4, label=lab)
                for f0, zsh, clean in pts:
                    ax.plot(f0, zsh, "o", ms=4, color=col, mfc=col if clean else "white")
            fixed = [r for r in rs if r["geom"] == "single" and r["fs"] == 0 and r["tag"] == "" and r["line"] == line and r["x"] == x]
            if fixed:
                ax.axhline(fixed[0]["zshare"], color=col, ls=":", lw=1.0)
        ax.set_xticks([0, 0.25, 0.5, 1.0]); ax.set_xticklabels(["truth", "1/4", "1/2", "all"])
        ax.set_xlabel("start probe: share of the search's misses")
        ax.set_ylabel("z-share, phase beyond half a cell"); ax.set_ylim(0, 0.88)
        ax.set_title("d   probe released", loc="left", fontweight="bold")
        ax.legend(loc="upper left", frameon=False)

    name = "fig_flat" + ("_found" if found else "") + ("_extra" if extra else "")
    os.makedirs(out_dir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out_dir, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    print(f"wrote {os.path.join(out_dir, name)}.pdf and .png")


if __name__ == "__main__":
    main(*sys.argv[1:])
