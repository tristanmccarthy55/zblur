"""
@file fig_npj.py
@brief The paper's figures drawn one topic per figure in npj style (2026-11-26): Arial 5-7 pt, bold lowercase panel letters,
       88 or 180 mm wide, no titles inside panels (the legend carries them), every atom labelled with W80 and its peak phase.

    ~/hyperspy-bundle/bin/python analysis/fig_npj.py [out_dir]
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "analysis"))
sys.path.insert(0, os.path.join(HERE, "runs"))
import fig_deck as d  # noqa: E402  (kernel names, loading and the measures)
import window_loss as wl  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

MM = 1 / 25.4
BLUE, ORANGE, GREY, INK = "#2a78d6", "#eb6834", "#8a8984", "#1a1a1a"
plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "lines.linewidth": 1.0, "lines.markersize": 3.5, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK,
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "legend.frameon": False, "axes.titlesize": 7, "axes.titleweight": "normal", "pdf.fonttype": 42})


def label(fig, ax, s, dx=-0.0, dy=0.0):
    """npj panel letter: 8 pt bold lowercase, at the axes' top-left corner."""
    b = ax.get_position()
    fig.text(b.x0 + dx, b.y1 + dy, s, fontsize=8, fontweight="bold", va="bottom", ha="left")


def atom(ax, k, zhalf=5.0, xhalf=1.5, left=False, bottom=False, fail=None, text=True):
    """x-z section through the atom's peak, scaled to its own peak; labelled W80 and peak phase."""
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(0.5)
    ax.set_facecolor("black")              # outside the cut-out volume
    if d.kfile(k) is None:
        ax.set_xticks([]); ax.set_yticks([]); ax.text(0.5, 0.5, "not run", ha="center", va="center", transform=ax.transAxes, fontsize=6)
        return None
    m = d.measure(k); K, dz = m["K"], m["dz"]
    k0 = int(np.argmax(m["col"])); _, ky, kx = np.unravel_index(np.argmax(K), K.shape)
    xz = K[:, ky, :]; z = (np.arange(K.shape[0]) - k0) * dz; x = (np.arange(K.shape[2]) - kx) * d.DX; sel = np.abs(x) <= xhalf
    im = ax.imshow(xz[:, sel] / max(xz.max(), 1e-9), extent=[x[sel][0] - d.DX / 2, x[sel][-1] + d.DX / 2, z[-1] + dz / 2, z[0] - dz / 2],
                   aspect="auto", cmap="inferno", vmin=0, vmax=1, interpolation="nearest")
    ax.set_ylim(zhalf, -zhalf); ax.set_xlim(-xhalf, xhalf)
    ax.set_xticks([-1, 0, 1] if bottom else []); ax.set_yticks(([-4, -2, 0, 2, 4] if zhalf <= 6 else [-10, -5, 0, 5, 10]) if left else [])
    if bottom:
        ax.set_xlabel("x (Å)")
    if left:
        ax.set_ylabel("z (Å)")
    if text:
        ax.text(0.05, 0.97, fail if fail else f"$W_{{80}}$ {m['w80']:.2f} Å\npeak {m['peak']:.2f} rad", transform=ax.transAxes, color="white",
                fontsize=5.5, va="top", ha="left", linespacing=1.15)
    return im


def cbar(fig, im, rect, lab="phase / own peak", phase=False):
    cax = fig.add_axes(rect)
    cb = fig.colorbar(im, cax=cax); cb.set_label(lab, fontsize=6); cb.ax.tick_params(labelsize=5.5, width=0.5, length=2)
    cb.outline.set_linewidth(0.5)
    if phase:
        cb.set_ticks([-np.pi, 0, np.pi]); cb.set_ticklabels(["−π", "0", "π"])
    else:
        cb.set_ticks([0, 0.5, 1])


def save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=600 if ext == "png" else None)
    plt.close(fig); print("wrote", name)


def tableau(x):
    r = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "runs", "flat_tableaux.tsv")) if l.startswith(f"flat_a080_f0{x:02d}")][0]
    return float(r[4]), json.loads(r[8])


def phase_map(ax, x):
    defocus, ab = tableau(x)
    k = np.linspace(-0.08, 0.08, 401) / wl.LAM
    KX, KY = np.meshgrid(k, k)
    ch = wl.chi(KX, KY, defocus, ab)
    r = wl.LAM * np.hypot(KX, KY) * 1e3
    ph = np.angle(np.exp(1j * (ch - np.median(ch[r <= x]))))
    ph[r > 80] = np.nan
    im = ax.imshow(ph, cmap="twilight", vmin=-np.pi, vmax=np.pi, extent=[-80, 80, -80, 80], origin="lower", interpolation="nearest")
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(x * np.cos(t), x * np.sin(t), "--", color="white", lw=0.6)
    ax.plot(80 * np.cos(t), 80 * np.sin(t), "-", color=INK, lw=0.5)
    ax.set_xlim(-84, 84); ax.set_ylim(-84, 84); ax.set_aspect("equal"); ax.axis("off")
    return im


# ------------------------------------------------------------------------------------------------------------------- figures
def fig_atoms_aperture(out):
    """a: perfect probe, aperture x = 30, 40, 60, 80. b: 80 mrad aperture flat only to x = 30, 40, 60 (phase across the aperture,
    then the atom), each column the same x as above it. Two-pass refined probe."""
    fig = plt.figure(figsize=(180 * MM, 118 * MM))
    gs = GridSpec(3, 4, figure=fig, left=0.07, right=0.9, top=0.95, bottom=0.08, hspace=0.3, wspace=0.08, height_ratios=[1, 0.95, 1])
    for j, x in enumerate((30, 40, 60, 80)):
        ax = fig.add_subplot(gs[0, j]); im = atom(ax, d.kname({30: "p30", 40: "p40", 60: "p60", 80: "p80"}[x], "_afb2"), left=(j == 0))
        ax.set_title(f"aperture {x} mrad", pad=3)
        if j == 0:
            label(fig, ax, "a", dx=-0.055, dy=0.01)
    cbar(fig, im, [0.915, 0.67, 0.012, 0.28])
    for j, x in enumerate((30, 40, 60)):
        ax = fig.add_subplot(gs[1, j]); pim = phase_map(ax, x)
        ax.set_title(f"flat to {x} mrad", pad=2)
        if j == 0:
            label(fig, ax, "b", dx=-0.055, dy=0.01)
        ax2 = fig.add_subplot(gs[2, j]); im = atom(ax2, d.kname({30: "f30", 40: "f40", 60: "f60"}[x], "_afb2"), left=(j == 0), bottom=True)
    cbar(fig, pim, [0.915, 0.40, 0.012, 0.2], lab="phase (rad)", phase=True)
    cbar(fig, im, [0.915, 0.07, 0.012, 0.25])
    save(fig, out, "npj_atoms_aperture")


def fig_w80_vs_x(out):
    rs = d.tsv_rows()
    def pts(line, tagend, known=False):
        sel = [r for r in rs if r["line"] == line and float(r["fs_A"]) == 0 and (r["treatment"] == "true probe, fixed" if known else r["kernel"].endswith(tagend))]
        return sorted({float(r["x_mrad"]): float(r["w80_cont_A"]) for r in sel if float(r["x_mrad"]) <= 80}.items())
    fig, ax = plt.subplots(figsize=(88 * MM, 62 * MM)); fig.subplots_adjust(left=0.14, right=0.97, bottom=0.17, top=0.96)
    for line, col in (("A", BLUE), ("B", ORANGE)):
        f, k = pts(line, "_afb2"), pts(line, None, known=True)
        if line == "B":
            p80f = [p for p in pts("A", "_afb2") if p[0] == 80]; p80k = [p for p in pts("A", None, True) if p[0] == 80]
            f, k = sorted(f + p80f), sorted(k + p80k)
        ax.plot([p[0] for p in f], [p[1] for p in f], "-o", color=col)
        ax.plot([p[0] for p in k], [p[1] for p in k], ":", color=col, lw=1.0)
    ax.set_xlabel("x (mrad)"); ax.set_ylabel("$W_{80}$ (Å)"); ax.set_xticks([30, 40, 60, 80]); ax.set_ylim(0, 10); ax.set_xlim(26, 84)
    h = [Line2D([], [], color=BLUE, marker="o"), Line2D([], [], color=ORANGE, marker="o"), Line2D([], [], color=GREY, ls="-"),
         Line2D([], [], color=GREY, ls=":")]
    ax.legend(h, ["perfect probe, aperture x", "80 mrad aperture, flat only to x", "probe refined (two passes)", "probe known"], loc="upper right")
    save(fig, out, "npj_w80_vs_x")


def fig_a1_atoms(out):
    keys = ["p30", "p40", "p60", "p80", "f30", "f40", "f60"]
    waves = {"p30": 0.08, "p40": 0.15, "p60": 0.33, "p80": 0.59, "f30": 0.59, "f40": 0.59, "f60": 0.59}
    fig = plt.figure(figsize=(180 * MM, 74 * MM))
    gs = GridSpec(2, 7, figure=fig, left=0.07, right=0.92, top=0.86, bottom=0.12, hspace=0.12, wspace=0.07)
    for j, key in enumerate(keys):
        a = fig.add_subplot(gs[0, j]); im = atom(a, d.kname(key, "_pu_noa1"), left=(j == 0))
        a.set_title(f"{d.NAME[key]}\n{waves[key]:.2f} waves", pad=3, linespacing=1.2)
        b = fig.add_subplot(gs[1, j]); atom(b, d.kname(key, "_pu"), left=(j == 0), bottom=True)
        if j == 0:
            label(fig, a, "a", dx=-0.06, dy=0.06); label(fig, b, "b", dx=-0.06, dy=0.0)
    fig.text(0.008, 0.66, "A1 right", rotation=90, fontsize=7, va="center"); fig.text(0.008, 0.31, "A1 3.6 Å off", rotation=90, fontsize=7, va="center")
    cbar(fig, im, [0.935, 0.12, 0.01, 0.74])
    save(fig, out, "npj_a1_atoms")


def fig_a1_ladder(out):
    errs = [0, 0.25, 0.5, 1.0, 3.6]; tags = ["_pu_noa1", "_pu_a1x0p25", "_pu_a1x0p5", "_pu_a1x1p0", "_pu"]
    fig, axs = plt.subplots(1, 2, figsize=(120 * MM, 50 * MM)); fig.subplots_adjust(left=0.1, right=0.98, bottom=0.2, top=0.92, wspace=0.35)
    for ax, what, lab, s in ((axs[0], "w80", "$W_{80}$ (Å)", "a"), (axs[1], "peak", "peak phase (rad)", "b")):
        for key, col in (("p80", BLUE), ("f30", ORANGE)):
            ax.plot(errs, [d.measure(d.kname(key, t))[what] for t in tags], "-o", color=col)
        ax.axvline(0.77, color=GREY, ls=":", lw=0.8); ax.set_xlabel("A1 error in the start probe (Å)"); ax.set_ylabel(lab); ax.set_xlim(-0.15, 3.8)
        ax.text(0.83, 0.97, "π/4", transform=ax.get_xaxis_transform(), fontsize=6, color=GREY, va="top")
        label(fig, ax, s, dx=-0.08, dy=0.0)
    axs[0].legend([Line2D([], [], color=BLUE, marker="o"), Line2D([], [], color=ORANGE, marker="o")],
                  ["perfect 80 mrad", "flat to 30, opened to 80"], loc="center right", bbox_to_anchor=(1.0, 0.36))
    save(fig, out, "npj_a1_ladder")


DOSE_FAIL = {("p80", "_a1fitp", "1e6"): "no atom", ("f30", "_poisson", "1e6"): "no atom", ("f30", "_a1fitp", "1e7"): "broken",
             ("f30", "_a1fitp", "1e6"): "no atom"}
DOSE_ROWS = [("p80", "", "_poisson", "perfect 80, probe known", BLUE, "-"), ("p80", "_a1fit", "_a1fitp", "perfect 80, A1 refined", BLUE, "--"),
             ("f30", "", "_poisson", "flat to 30 at 80, probe known", ORANGE, "-"), ("f30", "_a1fit", "_a1fitp", "flat to 30 at 80, A1 refined", ORANGE, "--")]


def fig_dose_line(out):
    X = {None: 1e9, "1e8": 1e8, "1e7": 1e7, "1e6": 1e6}
    fig, ax = plt.subplots(figsize=(88 * MM, 66 * MM)); fig.subplots_adjust(left=0.13, right=0.97, bottom=0.2, top=0.96)
    for key, t0, tn, lab, col, ls in DOSE_ROWS:
        pts, bad = [], []
        for dose in (None, "1e8", "1e7", "1e6"):
            if (key, tn, dose) in DOSE_FAIL:
                bad.append(X[dose]); continue
            pts.append((X[dose], d.measure(d.kname(key, t0) if dose is None else d.kname(key, tn, dose))["w80"]))
        pts.sort(); ax.plot([p[0] for p in pts], [p[1] for p in pts], ls, color=col, marker="o", label=lab)
        for b in bad:
            ax.plot(b, 9.0 if col == BLUE else 9.6, "x", color=col, ms=4, mew=1)
    ax.set_xscale("log"); ax.set_xlim(6e5, 2e9); ax.set_ylim(0, 14)
    ax.set_xticks([1e6, 1e7, 1e8, 1e9]); ax.set_xticklabels(["$10^6$", "$10^7$", "$10^8$", "no noise"]); ax.minorticks_off()
    ax.set_xlabel("dose (e Å$^{-2}$)"); ax.set_ylabel("$W_{80}$ (Å)")
    ax.legend(loc="upper right"); ax.text(7e5, 10.6, "× no atom or broken", fontsize=6, color=GREY, va="center")
    save(fig, out, "npj_dose_line")


def fig_dose_atoms(out):
    fig = plt.figure(figsize=(120 * MM, 118 * MM))
    gs = GridSpec(4, 4, figure=fig, left=0.17, right=0.89, top=0.94, bottom=0.08, hspace=0.08, wspace=0.07)
    for i, (key, t0, tn, lab, col, ls) in enumerate(DOSE_ROWS):
        for j, dose in enumerate((None, "1e8", "1e7", "1e6")):
            ax = fig.add_subplot(gs[i, j])
            im = atom(ax, d.kname(key, t0) if dose is None else d.kname(key, tn, dose), left=(j == 0), bottom=(i == 3),
                      fail=DOSE_FAIL.get((key, tn, dose)))
            if i == 0:
                ax.set_title("no noise" if dose is None else f"$10^{dose[-1]}$ e Å$^{{-2}}$", pad=3)
            if j == 0:
                ax.text(-0.55, 0.5, lab.replace(", ", "\n"), transform=ax.transAxes, rotation=90, va="center", ha="center", fontsize=6)
    cbar(fig, im, [0.905, 0.08, 0.015, 0.86])
    save(fig, out, "npj_dose_atoms")


def fig_focal_spread(out):
    rs = d.tsv_rows()
    def line(fs, ln):
        return sorted((float(r["x_mrad"]), float(r["w80_cont_A"])) for r in rs if r["kernel"].endswith("_pu_noa1") and float(r["fs_A"]) == fs and r["line"] == ln)
    fig, axs = plt.subplots(1, 4, figsize=(180 * MM, 52 * MM), sharey=True)
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.2, top=0.86, wspace=0.12)
    for j, (ax, fs) in enumerate(zip(axs, (0.0, 2.5, 5.0, 7.29))):
        a = line(fs, "A"); b = sorted(line(fs, "B") + [p for p in a if p[0] == 80])
        ax.plot([p[0] for p in a], [p[1] for p in a], "-o", color=BLUE); ax.plot([p[0] for p in b], [p[1] for p in b], "-o", color=ORANGE)
        ax.set_xticks([30, 40, 60, 80]); ax.set_xlabel("x (mrad)"); ax.set_ylim(0, 21)
        ev = {0.0: "", 2.5: " (0.10 eV)", 5.0: " (0.21 eV)", 7.29: " (0.30 eV)"}[fs]
        ax.set_title(f"focal spread {fs:g} Å{ev}", pad=3)
        label(fig, ax, "abcd"[j], dx=-0.035, dy=0.04)
    axs[0].set_ylabel("$W_{80}$ (Å)")
    axs[0].legend([Line2D([], [], color=BLUE, marker="o"), Line2D([], [], color=ORANGE, marker="o")],
                  ["perfect probe, aperture x", "80 mrad, flat only to x"], loc="upper right")
    save(fig, out, "npj_focal_spread")


SRC = [(0, "_NL20"), (0.05, "_src0p05_NL20"), (0.1, "_src0p1_NL20"), (0.2, "_src0p2_NL20"), (0.4, "_src0p4_NL20")]


def fig_source(out):
    ms = [d.measure(d.P["p80"] + t) for s, t in SRC]; src = [s for s, t in SRC]
    fig, axs = plt.subplots(1, 2, figsize=(120 * MM, 50 * MM)); fig.subplots_adjust(left=0.1, right=0.98, bottom=0.2, top=0.92, wspace=0.35)
    for ax, what, lab, s in ((axs[0], "w80", "$W_{80}$ (Å)", "a"), (axs[1], "peak", "peak phase (rad)", "b")):
        ax.plot(src, [m[what] for m in ms], "-o", color=BLUE); ax.set_xlabel("source size (Å FWHM)"); ax.set_ylabel(lab); ax.set_ylim(0, None)
        label(fig, ax, s, dx=-0.08, dy=0.0)
    save(fig, out, "npj_source")


def fig_source_atoms(out):
    fig = plt.figure(figsize=(120 * MM, 52 * MM))
    gs = GridSpec(1, 4, figure=fig, left=0.08, right=0.89, top=0.86, bottom=0.18, wspace=0.07)
    for j, (s, t) in enumerate([SRC[0], SRC[2], SRC[3], SRC[4]]):
        ax = fig.add_subplot(gs[0, j]); im = atom(ax, d.P["p80"] + t, left=(j == 0), bottom=True, zhalf=6)
        ax.set_title(f"source {s:g} Å", pad=3)
    cbar(fig, im, [0.905, 0.18, 0.015, 0.68])
    save(fig, out, "npj_source_atoms")


def fig_ghost(out):
    ks = [(d.kname("f30"), "probe known", GREY), (d.kname("f30", "_afb2"), "C1 + A1 refined, two passes", ORANGE), (d.kname("f30", "_a1fit"), "A1 refined alone", BLUE)]
    fig = plt.figure(figsize=(180 * MM, 62 * MM))
    gs = GridSpec(1, 5, figure=fig, left=0.06, right=0.89, top=0.9, bottom=0.18, wspace=0.08, width_ratios=[2.2, 0.35, 1, 1, 1])
    ax = fig.add_subplot(gs[0, 0])
    for k, lab, col in ks:
        m = d.measure(k); c = m["col"]; k0 = int(np.argmax(c)); z = (np.arange(len(c)) - k0) * m["dz"]
        ax.plot(z, c, "-o", color=col, ms=2.5, label=lab)
    ax.axhline(0, color=GREY, lw=0.4); ax.set_xlim(-12, 12); ax.set_xlabel("depth from the peak (Å)"); ax.set_ylabel("core phase (rad)")
    ax.legend(loc="upper right"); label(fig, ax, "a", dx=-0.045, dy=0.01)
    for j, (k, lab, col) in enumerate(ks):
        a = fig.add_subplot(gs[0, 2 + j]); im = atom(a, k, left=(j == 0), bottom=True); a.set_title(lab.replace(", ", ",\n"), pad=3, linespacing=1.2)
        if j == 0:
            label(fig, a, "b", dx=-0.05, dy=0.05)
    cbar(fig, im, [0.905, 0.18, 0.01, 0.72])
    save(fig, out, "npj_ghost")


def fig_camera(out):
    pairs = [(d.kname("p80"), "perfect 80\n712 px"), (d.P["p80"] + "_detb4_NL20", "perfect 80\n178 px camera"),
             (d.kname("f30"), "flat to 30\n914 px"), (d.P["f30"] + "_detb4_NL40", "flat to 30\n228 px camera")]
    fig = plt.figure(figsize=(120 * MM, 56 * MM))
    gs = GridSpec(1, 4, figure=fig, left=0.08, right=0.89, top=0.84, bottom=0.17, wspace=0.07)
    for j, (k, lab) in enumerate(pairs):
        ax = fig.add_subplot(gs[0, j]); im = atom(ax, k, left=(j == 0), bottom=True); ax.set_title(lab, pad=3, linespacing=1.2)
    cbar(fig, im, [0.905, 0.17, 0.015, 0.67])
    save(fig, out, "npj_camera")


def fig_kde(out):
    def w80_kde(p, dz, h, frac=0.8, step=0.005):
        q = np.clip(p, 0, None); z = np.arange(-10, len(q) * dz + 10, step); dens = np.zeros_like(z)
        for i, v in enumerate(q):
            dens[(z >= i * dz) & (z < (i + 1) * dz)] += v / dz
        if h > 0:
            g = np.arange(-5 * h, 5 * h + step, step); g = np.exp(-0.5 * (g / h) ** 2); dens = np.convolve(dens, g / g.sum(), mode="same")
        F = np.cumsum(dens) * step; F /= F[-1]; j = np.searchsorted(F, F + frac); ok = j < len(z)
        return float((z[j[ok]] - z[ok]).min())
    hs = [0, 0.25, 0.5, 1.0]
    fig, ax = plt.subplots(figsize=(88 * MM, 62 * MM)); fig.subplots_adjust(left=0.14, right=0.78, bottom=0.17, top=0.96)
    ends = []
    for key in ("p30", "p40", "p60", "p80", "f30", "f40", "f60"):
        c = d.measure(d.kname(key))["col"]; v = [w80_kde(c, 1.0, h) for h in hs]; rel = [100 * (x / v[0] - 1) for x in v]
        hi = key in ("p80", "p30")
        ax.plot(hs, rel, "-o", color=(BLUE if key == "p80" else ORANGE) if hi else "#bdbcb6", ms=3 if hi else 2, lw=1.0 if hi else 0.7, zorder=3 if hi else 2)
        ends.append([rel[-1], key])
    ends.sort()
    for i in range(1, len(ends)):
        ends[i][0] = max(ends[i][0], ends[i - 1][0] + 3.0)
    for y, key in ends:
        ax.text(1.04, y, d.NAME[key], fontsize=6, va="center", transform=ax.get_yaxis_transform() if False else ax.transData)
    ax.set_xlabel("KDE bandwidth σ (Å)"); ax.set_ylabel("$W_{80}$ change (%)"); ax.set_xlim(-0.04, 1.0); ax.set_xticks(hs)
    ax.set_clip_on(False)
    save(fig, out, "npj_kde")


def fig_p30noise(out):
    ks = [(d.kname("p30"), "no noise\n1 Å slices", None), (d.kname("p30", "_poisson", "1e8"), "$10^8$ e Å$^{-2}$\n1 Å slices", "diverged"),
          (d.kname("p30", "_poisson", "1e7"), "$10^7$ e Å$^{-2}$\n1 Å slices", "diverged"), (d.kname("p30", "_poisson", "1e6"), "$10^6$ e Å$^{-2}$\n1 Å slices", "diverged"),
          (d.kname("p30", nl=20), "no noise\n2 Å slices", None)]
    ZH = 10.0
    fig = plt.figure(figsize=(150 * MM, 60 * MM))
    gs = GridSpec(1, 5, figure=fig, left=0.07, right=0.9, top=0.84, bottom=0.16, wspace=0.07)
    for j, (k, lab, fail) in enumerate(ks):
        ax = fig.add_subplot(gs[0, j]); im = atom(ax, k, left=(j == 0), bottom=True, fail=fail, zhalf=ZH); ax.set_title(lab, pad=3, linespacing=1.2)
    cbar(fig, im, [0.915, 0.16, 0.012, 0.68])
    save(fig, out, "npj_p30noise")


def main(*args):
    out = args[0] if args else os.path.join(os.path.dirname(HERE), "ptychoshelves-clean", "report", "aberrations_paper", "figs", "npj")
    for f in (fig_atoms_aperture, fig_w80_vs_x, fig_a1_atoms, fig_a1_ladder, fig_dose_line, fig_dose_atoms, fig_focal_spread, fig_source,
              fig_source_atoms, fig_ghost, fig_camera, fig_kde, fig_p30noise):
        f(out)


if __name__ == "__main__":
    main(*sys.argv[1:])
