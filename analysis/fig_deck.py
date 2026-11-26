"""
@file fig_deck.py
@brief First-attempt figures for the paper's figure list (2026-11-26): one PNG per figure, drawn from the lone-atom kernels.

fig2_a1        which of the search's errors matters: the seven atoms with A1 right and with A1 3.6 A off; the A1 ladder.
fig3_dose      W80 against dose (Poisson, probe known and A1 fitted) with the atoms at each dose.
fig4_coherence W80 against x at four focal spreads (one coherent mode); source size on the perfect 80 mrad probe.
s_ghost        flat to 30: true probe, two passes, A1 alone (depth profiles and sections).
s_camera       the sim's patterns against a 4 x 4 camera, no noise.
s_kde          W80 of the histogram profile against a Gaussian KDE of growing bandwidth.
s_p30noise     perfect 30 mrad under noise on 1 A slices, and its noise-free atom on 2 A slices.

    ~/hyperspy-bundle/bin/python analysis/fig_deck.py [out_dir]
"""
from __future__ import annotations

import csv
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "analysis"))
import flat_plot as fp  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

C = dict(bg="#fcfcfb", ink="#0b0b0b", ink2="#52514e", grid="#e4e3df", perfect="#2a78d6", flat="#eb6834", third="#1baf7a")
plt.rcParams.update({"font.size": 12, "axes.edgecolor": C["ink2"], "axes.labelcolor": C["ink"], "xtick.color": C["ink2"],
                     "ytick.color": C["ink2"], "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": C["bg"],
                     "axes.facecolor": C["bg"], "savefig.facecolor": C["bg"], "lines.linewidth": 2, "axes.titlesize": 13,
                     "axes.titleweight": "bold", "axes.titlelocation": "left", "legend.frameon": False})
DX = 0.04922
P = {"p30": "Pb_a030_ideal_df0_single_st0p1_w8_s0p25_bz40p000", "p40": "Pb_a040_ideal_df0_single_st0p1_w8_s0p25_bz40p000",
     "p60": "Pb_a060_ideal_df0_single_st0p1_w8_s0p25_bz20p000", "p80": "Pb_a080_ideal_df0_single_st0p1_w8_s0p25_bz20p000",
     "f30": "Pb_a080_flat30_single_st0p1_w8_s0p25_bz40p000_rs90", "f40": "Pb_a080_flat40_single_st0p1_w8_s0p25_bz40p000",
     "f60": "Pb_a080_flat60_single_st0p1_w8_s0p25_bz20p000"}
NAME = {"p30": "perfect 30", "p40": "perfect 40", "p60": "perfect 60", "p80": "perfect 80", "f30": "flat to 30", "f40": "flat to 40",
        "f60": "flat to 60"}
NL = {"p30": 40, "p40": 40, "p60": 20, "p80": 20, "f30": 40, "f40": 40, "f60": 20}


def kname(key, tag="", dose=None, nl=None):
    return f"{P[key]}{'' if dose is None else '_dose' + dose}_NL{nl or NL[key]}{tag}"


def kfile(k):
    fs = sorted(glob.glob(os.path.expanduser(f"~/Desktop/zblur/*/psf/psf_{k}_vol.npy")))
    return fs[0] if fs else None


def dz_of(k):
    m = re.search(r"_bz(\d+)p(\d+).*_NL(\d+)", k)
    return float(f"{m.group(1)}.{m.group(2)}") / int(m.group(3))


def measure(k):
    K = fp.blobs.zp.load_phase(kfile(k))
    col = fp.pg.core_column(K, DX)
    return dict(K=K, col=col, w80=fp.width_cont(col, dz_of(k)), peak=float(K.max()), dz=dz_of(k))


def section(ax, k, note=None, zhalf=8.0, xhalf=1.5, title=None, fail=None):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor("black")
    for s in ax.spines.values():
        s.set_visible(False)
    if title:
        ax.set_title(title, fontsize=11, fontweight="normal", loc="center", color=C["ink"], pad=4)
    if kfile(k) is None:
        ax.set_facecolor("#e9e8e4"); ax.text(0.5, 0.5, "not run", ha="center", va="center", transform=ax.transAxes, fontsize=10, color=C["ink2"])
        return None
    m = measure(k)
    K, dz = m["K"], m["dz"]
    k0 = int(np.argmax(m["col"]))
    _, ky, kx = np.unravel_index(np.argmax(K), K.shape)
    xz = K[:, ky, :]
    z = (np.arange(K.shape[0]) - k0) * dz; x = (np.arange(K.shape[2]) - kx) * DX; sel = np.abs(x) <= xhalf
    ax.imshow(xz[:, sel] / max(xz.max(), 1e-9), extent=[x[sel][0] - DX / 2, x[sel][-1] + DX / 2, z[-1] + dz / 2, z[0] - dz / 2],
              aspect="auto", cmap="inferno", vmin=0, vmax=1, interpolation="nearest")
    ax.set_ylim(zhalf, -zhalf); ax.set_xlim(-xhalf, xhalf)
    txt = fail if fail else (note if note is not None else f"{m['w80']:.2f} Å\n{m['peak']:.2f} rad")
    ax.text(0.05, 0.96, txt, transform=ax.transAxes, color="white", fontsize=10, va="top")
    return m


def tsv_rows():
    rs = [{k.strip(): v for k, v in r.items()} for r in csv.DictReader(open(os.path.join(HERE, "results", "flat_plot.tsv")), delimiter="\t")]
    return [r for r in rs if r["geom"] == "single" and r["camera_bin"] == "1" and r["dose_e_per_A2"] == "-"]


def save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, f"{name}.png"), dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def fig2_a1(out):
    keys = ["p30", "p40", "p60", "p80", "f30", "f40", "f60"]
    waves = {"p30": 0.08, "p40": 0.15, "p60": 0.33, "p80": 0.59, "f30": 0.59, "f40": 0.59, "f60": 0.59}
    fig = plt.figure(figsize=(14, 9.2))
    gs = GridSpec(3, 7, figure=fig, height_ratios=[1, 1, 1.25], hspace=0.35, wspace=0.08)
    for j, key in enumerate(keys):
        a = fig.add_subplot(gs[0, j]); section(a, kname(key, "_pu_noa1"), title=f"{NAME[key]}\n{waves[key]:.2f} waves")
        b = fig.add_subplot(gs[1, j]); section(b, kname(key, "_pu"))
        if j == 0:
            a.text(-0.12, 0.5, "A1 right", transform=a.transAxes, rotation=90, va="center", ha="right", fontsize=12, color=C["ink"])
            b.text(-0.12, 0.5, "A1 3.6 Å off", transform=b.transAxes, rotation=90, va="center", ha="right", fontsize=12, color=C["ink"])
    fig.text(0.125, 0.975, "a   the seven atoms, A1 right and A1 left 3.6 Å off (x–z through the peak, ±8 Å deep, ±1.5 Å across; "
             "waves = that error at each probe's own edge)", fontsize=13, fontweight="bold", color=C["ink"])
    errs = [0, 0.25, 0.5, 1.0, 3.6]
    tags = ["_pu_noa1", "_pu_a1x0p25", "_pu_a1x0p5", "_pu_a1x1p0", "_pu"]
    sub = gs[2, :].subgridspec(1, 2, wspace=0.28)
    for i, (what, lab) in enumerate((("w80", "W80 (Å)"), ("peak", "peak phase (rad)"))):
        ax = fig.add_subplot(sub[0, i])
        for key, col in (("p80", C["perfect"]), ("f30", C["flat"])):
            v = [measure(kname(key, t))[what] for t in tags]
            ax.plot(errs, v, "-o", color=col, ms=8, label=NAME[key] + (" mrad" if key == "p80" else ", opened to 80"))
            ax.text(errs[-1] + 0.08, v[-1], NAME[key], color=C["ink2"], va="center", fontsize=11)
        ax.axvline(0.77, color=C["ink2"], ls=":", lw=1.5)
        ax.text(0.82, 0.5, "π/4 criterion, 0.77 Å", transform=ax.get_xaxis_transform(), rotation=90, fontsize=10, color=C["ink2"], va="center")
        ax.set_xlabel("A1 error in the start probe (Å)"); ax.set_ylabel(lab); ax.set_xlim(-0.1, 4.3)
        ax.grid(axis="y", color=C["grid"], lw=0.8); ax.set_axisbelow(True)
        ax.set_title(("b   " if i == 0 else "c   ") + ("A1 ladder: depth extent" if i == 0 else "A1 ladder: peak phase"))
        if i == 0:
            ax.legend(loc="upper left", fontsize=10, bbox_to_anchor=(0.25, 1.0))
    save(fig, out, "fig2_a1")


def fig3_dose(out):
    doses = [None, "1e8", "1e7", "1e6"]
    X = {None: 1e9, "1e8": 1e8, "1e7": 1e7, "1e6": 1e6}
    rows = [("p80", "", "_poisson", "perfect 80, probe known", C["perfect"], "-"),
            ("p80", "_a1fit", "_a1fitp", "perfect 80, A1 fitted", C["perfect"], "--"),
            ("f30", "", "_poisson", "flat to 30 at 80, probe known", C["flat"], "-"),
            ("f30", "_a1fit", "_a1fitp", "flat to 30 at 80, A1 fitted", C["flat"], "--")]
    failed = {("p80", "_a1fitp", "1e6"): "no atom", ("f30", "_poisson", "1e6"): "no atom", ("f30", "_a1fitp", "1e7"): "broken",
              ("f30", "_a1fitp", "1e6"): "no atom"}
    fig = plt.figure(figsize=(15, 7.6))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1, 1.05], wspace=0.12)
    ax = fig.add_subplot(gs[0, 0])
    for key, t0, tn, lab, col, ls in rows:
        pts, bad = [], []
        for d in doses:
            k = kname(key, t0) if d is None else kname(key, tn, d)
            if (key, tn, d) in failed:
                bad.append(X[d]); continue
            pts.append((X[d], measure(k)["w80"]))
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], ls, color=col, marker="o", ms=8, label=lab)
        for b in bad:
            ax.plot(b, 9.0 if col == C["perfect"] else 9.5, "x", color=col, ms=10, mew=2.5)
    ax.text(1.4e6, 9.95, "× = no atom or broken", fontsize=10, color=C["ink2"])
    ax.set_xscale("log"); ax.set_xlim(6e5, 2e9); ax.set_ylim(0, 10.2)
    ax.set_xticks([1e6, 1e7, 1e8, 1e9]); ax.set_xticklabels(["$10^6$", "$10^7$", "$10^8$", "no\nnoise"])
    ax.set_xlabel("dose (e Å$^{-2}$, over the scanned area)"); ax.set_ylabel("W80 (Å)")
    ax.grid(axis="y", color=C["grid"], lw=0.8); ax.set_axisbelow(True)
    ax.set_title("a   depth extent against dose (Poisson likelihood, 712 px)")
    ax.legend(loc="upper right", fontsize=10, bbox_to_anchor=(1.0, 0.9))
    ax.text(6.5e5, 0.35, "perfect 30: diverged at every dose on 1 Å slices", fontsize=10, color=C["ink2"])
    sub = gs[0, 1].subgridspec(4, 4, wspace=0.05, hspace=0.12)
    for i, (key, t0, tn, lab, col, ls) in enumerate(rows):
        for j, d in enumerate(doses):
            a = fig.add_subplot(sub[i, j])
            k = kname(key, t0) if d is None else kname(key, tn, d)
            section(a, k, title=("no noise" if d is None else f"$10^{d[-1]}$") if i == 0 else None, fail=failed.get((key, tn, d)))
            if j == 0:
                a.text(-0.08, 0.5, lab.replace(", ", "\n"), transform=a.transAxes, rotation=90, va="center", ha="right", fontsize=9.5,
                       color=C["ink"])
    fig.text(0.53, 0.93, "b   the atoms at each dose (x–z, ±8 Å deep, ±1.5 Å across; each scaled to its own peak)", fontsize=13,
             fontweight="bold", color=C["ink"])
    save(fig, out, "fig3_dose")


def fig4_coherence(out):
    rs = tsv_rows()
    def line(fs, line):
        pts = sorted((float(r["x_mrad"]), float(r["w80_cont_A"])) for r in rs if r["kernel"].endswith("_pu_noa1")
                     and float(r["fs_A"]) == fs and r["line"] == line)
        return pts
    fig = plt.figure(figsize=(15, 9.4))
    gs = GridSpec(2, 4, figure=fig, height_ratios=[1, 1.05], hspace=0.42, wspace=0.18)
    for j, fs in enumerate((0.0, 2.5, 5.0, 7.29)):
        ax = fig.add_subplot(gs[0, j])
        a = line(fs, "A"); b = line(fs, "B")
        p80 = [p for p in a if p[0] == 80]
        b = sorted(b + p80)
        ax.plot([p[0] for p in a], [p[1] for p in a], "-o", color=C["perfect"], ms=7, label="perfect probe, aperture x")
        ax.plot([p[0] for p in b], [p[1] for p in b], "-o", color=C["flat"], ms=7, label="80 mrad, flat only to x")
        ax.set_ylim(0, 21); ax.set_xticks([30, 40, 60, 80]); ax.set_xlabel("x (mrad)")
        ax.grid(axis="y", color=C["grid"], lw=0.8); ax.set_axisbelow(True)
        ev = {0.0: "", 2.5: " (0.10 eV)", 5.0: " (0.21 eV)", 7.29: " (0.30 eV)"}[fs]
        ax.set_title(("a   " if j == 0 else "") + f"focal spread {fs:g} Å{ev}", fontsize=12)
        if j == 0:
            ax.set_ylabel("W80 (Å)"); ax.legend(loc="upper right", fontsize=9.5)
        else:
            ax.set_yticklabels([])
    src = [0, 0.05, 0.1, 0.2, 0.4]
    tags = ["_NL20", "_src0p05_NL20", "_src0p1_NL20", "_src0p2_NL20", "_src0p4_NL20"]
    ms = [measure(P["p80"] + t) for t in tags]
    for i, (what, lab) in enumerate((("w80", "W80 (Å)"), ("peak", "peak phase (rad)"))):
        ax = fig.add_subplot(gs[1, i])
        ax.plot(src, [m[what] for m in ms], "-o", color=C["perfect"], ms=8)
        ax.set_xlabel("source size (Å FWHM)"); ax.set_ylabel(lab)
        ax.grid(axis="y", color=C["grid"], lw=0.8); ax.set_axisbelow(True)
        ax.set_title(("b   source size: depth extent" if i == 0 else "c   source size: peak phase"), fontsize=12)
        ax.set_ylim(0, None)
    sub = gs[1, 2:].subgridspec(1, 4, wspace=0.06)
    for j, (s, t) in enumerate(zip((0, 0.1, 0.2, 0.4), ("_NL20", "_src0p1_NL20", "_src0p2_NL20", "_src0p4_NL20"))):
        a = fig.add_subplot(sub[0, j]); section(a, P["p80"] + t, title=f"source {s:g} Å", zhalf=6)
    fig.text(0.53, 0.47, "d   perfect 80 mrad, true probe held: the atom as the source grows (x–z, ±6 Å deep)", fontsize=12, fontweight="bold", color=C["ink"])
    save(fig, out, "fig4_coherence")


def s_ghost(out):
    ks = [(kname("f30"), "true probe held", C["ink2"]), (kname("f30", "_afb2"), "C1 + A1 fitted, two passes", C["flat"]),
          (kname("f30", "_a1fit"), "A1 fitted alone", C["perfect"])]
    fig = plt.figure(figsize=(14, 6.2))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.1, 1], wspace=0.15)
    ax = fig.add_subplot(gs[0, 0])
    for k, lab, col in ks:
        m = measure(k); c = m["col"]; k0 = int(np.argmax(c)); z = (np.arange(len(c)) - k0) * m["dz"]
        ax.plot(z, c, "-o", color=col, ms=5, label=f"{lab}: W80 {m['w80']:.2f} Å, peak {m['peak']:.2f} rad")
    ax.axhline(0, color=C["grid"], lw=1); ax.set_xlim(-12, 12)
    ax.set_xlabel("depth from the peak (Å)"); ax.set_ylabel("core phase minus ring (rad)")
    ax.set_title("a   flat to 30 at 80 mrad: depth profiles"); ax.legend(loc="upper right", fontsize=10)
    sub = gs[0, 1].subgridspec(1, 3, wspace=0.06)
    for j, (k, lab, col) in enumerate(ks):
        section(fig.add_subplot(sub[0, j]), k, title=lab.replace(", ", "\n"))
    fig.text(0.56, 1.0, "b   the atoms (x–z, ±8 Å)", fontsize=13, fontweight="bold", color=C["ink"])
    save(fig, out, "s_ghost")


def s_camera(out):
    pairs = [(kname("p80"), "perfect 80\n712 px"), (P["p80"] + "_detb4_NL20", "perfect 80\n178 px camera"),
             (kname("f30"), "flat to 30\n914 px"), (P["f30"] + "_detb4_NL40", "flat to 30\n228 px camera")]
    fig, axs = plt.subplots(1, 4, figsize=(11, 4.6), gridspec_kw=dict(wspace=0.06))
    for a, (k, lab) in zip(axs, pairs):
        section(a, k, title=lab)
    save(fig, out, "s_camera")


def s_kde(out):
    def w80_kde(p, dz, h, frac=0.8, step=0.005):
        q = np.clip(p, 0, None); z = np.arange(-10, len(q) * dz + 10, step); dens = np.zeros_like(z)
        for i, v in enumerate(q):
            dens[(z >= i * dz) & (z < (i + 1) * dz)] += v / dz
        if h > 0:
            g = np.arange(-5 * h, 5 * h + step, step); g = np.exp(-0.5 * (g / h) ** 2); dens = np.convolve(dens, g / g.sum(), mode="same")
        F = np.cumsum(dens) * step; F /= F[-1]; j = np.searchsorted(F, F + frac); ok = j < len(z)
        return float((z[j[ok]] - z[ok]).min())
    hs = [0, 0.25, 0.5, 1.0]
    fig, ax = plt.subplots(figsize=(9, 5.6))
    ends = []
    for key in ("p30", "p40", "p60", "p80", "f30", "f40", "f60"):
        c = measure(kname(key))["col"]; v = [w80_kde(c, 1.0, h) for h in hs]; rel = [100 * (x / v[0] - 1) for x in v]
        hi = key in ("p80", "p30")
        ax.plot(hs, rel, "-o", color=(C["perfect"] if key == "p80" else C["flat"]) if hi else "#b9b8b2", ms=7 if hi else 5,
                lw=2 if hi else 1.2, zorder=3 if hi else 2)
        ends.append([rel[-1], key, hi])
    ends.sort()
    for i in range(1, len(ends)):                      # keep the end labels at least 2.6 % apart
        ends[i][0] = max(ends[i][0], ends[i - 1][0] + 2.6)
    for y, key, hi in ends:
        ax.text(1.03, y, NAME[key], fontsize=10, color=C["ink"] if hi else C["ink2"], va="center")
    ax.set_xlabel("KDE bandwidth σ (Å); 0 = the histogram used"); ax.set_ylabel("W80 change (%)"); ax.set_xlim(-0.05, 1.25)
    ax.grid(axis="y", color=C["grid"], lw=0.8); ax.set_axisbelow(True)
    ax.set_title("A KDE lengthens every atom, the compact ones most (true probe, no noise)")
    save(fig, out, "s_kde")


def s_p30noise(out):
    ks = [(kname("p30"), "no noise\n1 Å slices", None), (kname("p30", "_poisson", "1e8"), "$10^8$\n1 Å slices", "diverged"),
          (kname("p30", "_poisson", "1e7"), "$10^7$\n1 Å slices", "diverged"), (kname("p30", "_poisson", "1e6"), "$10^6$\n1 Å slices", "diverged"),
          (kname("p30", nl=20), "no noise\n2 Å slices", None)]
    fig, axs = plt.subplots(1, 5, figsize=(13, 4.8), gridspec_kw=dict(wspace=0.06))
    for a, (k, lab, fail) in zip(axs, ks):
        section(a, k, title=lab, fail=fail, zhalf=12)
    save(fig, out, "s_p30noise")


def main(*args):
    out = args[0] if args else os.path.join(os.path.dirname(HERE), "ptychoshelves-clean", "report", "aberrations_paper", "figs", "deck")
    for f in (fig2_a1, fig3_dose, fig4_coherence, s_ghost, s_camera, s_kde, s_p30noise):
        f(out)


if __name__ == "__main__":
    main(*sys.argv[1:])
