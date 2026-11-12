"""
@file flat_plot.py
@brief The flat experiment's figure: one lone Pb atom's depth FWHM against angle, the perfect probe against an 80 mrad
aperture flat only to x (runs/gflat*.tsv, runs/flat_family.py, prediction P6 in runs/PREDICTIONS.md).

Line A, the perfect probe with an aperture of x mrad; line B, an 80 mrad aperture over a corrector flat only to x (at x = 80 it
is line A's point). One panel per probe treatment found: the true probe held fixed, and the probe released from a start
(_pu: the crystal search's misses; _pu0: the truth; _pu_noho: the misses with the higher orders left out). Depth FWHM is
analysis/blobs.py's (read off the centre column); a kernel blobs.py calls NOT CLEAN is drawn hollow. Points from the earlier
runs of the same recipe (G0d, G1: grouping 32 against gflat's "32;16") are drawn as squares. With focal-spread rows (runs with
fs), a second figure: the same lines at that focal spread, with sqrt(lambda / (pi fs)), the angle past which the spread scrambles
the probe's phase (my estimate, untested), marked.

    ~/hyperspy-bundle/bin/python analysis/flat_plot.py                 # ~/Desktop/zblur/gflat, g0d, g1
    ~/hyperspy-bundle/bin/python analysis/flat_plot.py <pack dir> ...
"""
from __future__ import annotations

import csv
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import blobs  # noqa: E402
import predict_groups as pg  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAM = 0.0196877
# the lone atom (gflat*) and the grid of atoms on one plane (ggrid: the experimentalist version), each with its own scan and boxes
RECIPE = {"single": dict(step_A="0.1", win_A="8", det_mrad="200", slice_A="0.25"),
          "grid": dict(step_A="0.5", win_A="20", det_mrad="200", slice_A="0.25")}
BOX_HIGH = {"single": 20.0, "grid": 27.525}     # the box above 40 mrad; 30 and 40 mrad sit in a 40 A box in both
TREATMENT = {"": "true probe, fixed", "_pu": "released from the search's misses", "_pu0": "released from the truth",
             "_pu_noho": "released, misses, higher orders not given"}
EARLY = ("g0d", "g1")


def rtag(kernel):
    m = re.search(r"_NL\d+(_\w+)?$", kernel)
    return (m.group(1) or "") if m else ""


def load(packs):
    """blobs.load for both geometries (blobs.load keeps lone atoms only)."""
    out = []
    for pack in packs:
        for man in sorted(glob.glob(os.path.join(pack, "manifest_*.tsv"))):
            for r in csv.DictReader(open(man), delimiter="\t"):
                f = os.path.join(pack, "psf", f"psf_{r['kernel']}_vol.npy")
                if r.get("geom") not in RECIPE or not os.path.exists(f):
                    continue
                r = dict(r, file=f, K=blobs.zp.load_phase(f), dz=float(r["dz_A"]), dx=float(r["dx_A"]), fs=float(r["fs_A"]),
                         src=float(r["src_A"]), T=(None if r["T_K"] in ("-", "") else float(r["T_K"])))
                r["z0"] = int(r["zdrop"]) * r["dz"]
                out.append(r)
    return out


def rows(packs):
    seen, out = set(), []
    for r in load(packs):
        if any(r.get(k) != v for k, v in RECIPE[r["geom"]].items()) or r["T"] is not None or r["src"]:
            continue
        if r["kernel"] in seen:
            continue
        seen.add(r["kernel"])
        c = blobs.core(r)
        aber = r["aber"]
        x = float(r["alpha"]) if aber == "ideal" else (float(aber[4:]) if aber.startswith("flat") else None)
        if x is None:
            continue
        # gflat's own slicing and boxes: 1 A recon slices, a 40 A box at 30 and 40 mrad, 20 A above (earlier runs hold others)
        if abs(r["dz"] - 1.0) > 0.05 or abs(float(r["boxz_A"]) - (40.0 if x <= 40 else BOX_HIGH[r["geom"]])) > 0.01:
            continue
        early = os.path.basename(os.path.dirname(os.path.dirname(r["file"]))) in EARLY
        col = pg.core_column(r["K"], r["dx"])
        out.append(dict(kernel=r["kernel"], geom=r["geom"], line="A" if aber == "ideal" else "B", x=x, fs=r["fs"], tag=rtag(r["kernel"]),
                        fwhm=c["fwhm_z"], dip=c["dip"], side=c["side"], peak=c["peak"], w80=width(col, r["dz"], 0.8),
                        zshare=pg.zshare_of(col, r["dz"]), clean=bool(np.isfinite(c["fwhm_z"]) and c["dip"] > -0.9), early=early))
    return out


def width(p, dz, frac):
    """the shortest run of slices that holds frac of the centre column's positive phase (counts tails, which a FWHM does not)"""
    q = np.clip(p, 0, None)
    c = np.concatenate([[0.0], np.cumsum(q)])
    n, best = len(q), len(q)
    for i in range(n):
        j = int(np.searchsorted(c, c[i] + frac * c[-1], side="left"))
        if j <= n:
            best = min(best, j - i)
    return best * dz


METRIC = {"w80": ("W80: depth holding 80 % of the atom's phase, A", 3.0),
          "fwhm": ("depth FWHM (core only, misses tails), A", 1.25),
          "zshare": ("z-share: phase beyond half a cell along z", 0.012)}


def panel(ax, rs, title, fs=0.0, metric="w80"):
    a = sorted([r for r in rs if r["line"] == "A"], key=lambda r: r["x"])
    b = sorted([r for r in rs if r["line"] == "B"], key=lambda r: r["x"])
    a80 = [r for r in a if r["x"] == 80]
    if a80:
        b = b + [dict(a80[0], line="B")]          # flat to 80 is the perfect 80 mrad probe
        b = sorted(b, key=lambda r: r["x"])
    for pts, col, lab in ((a, "#1f5fa8", "perfect probe, aperture x"), (b, "#c8501e", "80 mrad aperture, flat only to x")):
        if not pts:
            continue
        ax.plot([p["x"] for p in pts], [p[metric] for p in pts], "-", color=col, lw=1.6, label=lab, zorder=1)
        for p in pts:
            ax.plot(p["x"], p[metric], "s" if p["early"] else "o", ms=7, color=col, mfc=col if p["clean"] else "white", mew=1.6, zorder=2)
    label, floor = METRIC[metric]
    ax.axhline(floor, color="0.6", lw=4, alpha=0.4, zorder=0, label="the perfect probe at 80-100 mrad (1 A recon slices)")
    if fs:
        tc = 1000 * np.sqrt(LAM / (np.pi * fs))
        ax.axvline(tc, color="0.3", ls="--", lw=1.1, label=f"sqrt(lambda / pi fs) = {tc:.0f} mrad (estimate)")
    ax.set_xlabel("x: aperture (perfect probe) or flat angle (80 mrad aperture), mrad")
    ax.set_ylabel(label, fontsize=9)
    vals = [p[metric] for p in rs if np.isfinite(p[metric])]
    ax.set_xlim(25, 85); ax.set_ylim(0, 1.12 * max(vals + [floor * 2]))
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)


def main(*packs):
    packs = packs or tuple(os.path.expanduser(f"~/Desktop/zblur/{d}") for d in ("gflat",) + EARLY)
    rs = rows([p for p in packs if os.path.isdir(p)])
    if not rs:
        raise SystemExit("no kernels of the flat recipe found")
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True); os.makedirs(os.path.join(HERE, "figs"), exist_ok=True)
    with open(os.path.join(HERE, "results", "flat_plot.tsv"), "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["geom", "line", "x_mrad", "fs_A", "treatment", "w80_A", "fwhm_z_A", "zshare", "dip", "side", "peak_rad", "clean", "earlier_run", "kernel"])
        for r in sorted(rs, key=lambda r: (r["geom"], r["fs"], r["tag"], r["line"], r["x"])):
            w.writerow([r["geom"], r["line"], f"{r['x']:g}", f"{r['fs']:g}", TREATMENT.get(r["tag"], r["tag"]), f"{r['w80']:.1f}", f"{r['fwhm']:.2f}",
                        f"{r['zshare']:.3f}", f"{r['dip']:.2f}", f"{r['side']:.2f}", f"{r['peak']:.3f}", int(r["clean"]), int(r["early"]), r["kernel"]])
    for geom, fs in sorted({(r["geom"], r["fs"]) for r in rs}):
        sub = [r for r in rs if r["fs"] == fs and r["geom"] == geom]
        tags = [t for t in TREATMENT if any(r["tag"] == t for r in sub)]
        fig, axes = plt.subplots(len(tags), len(METRIC), figsize=(5.2 * len(METRIC), 4.3 * len(tags)), squeeze=False)
        for row, t in zip(axes, tags):
            # the earlier fixed-probe runs fill line A wherever this treatment has no point of its own
            pts = [r for r in sub if r["tag"] == t]
            if t == "":
                have = {(r["line"], r["x"]) for r in pts if not r["early"]}
                pts = [r for r in pts if not r["early"] or (r["line"], r["x"]) not in have]
            for ax, m in zip(row, METRIC):
                panel(ax, pts, f"{TREATMENT[t]}" + (f", focal spread {fs:g} A" if fs else ""), fs, m)
        axes[0][0].legend(fontsize=7, loc="upper right")
        fig.suptitle(("One lone Pb atom" if geom == "single" else "Pb atoms on one plane, 3 A apart (kernel per atom)") + ": does opening past the flat region cost depth?  (filled: clean kernel; hollow: NOT CLEAN; "
                     "square: earlier run of the recipe)", fontsize=9)
        fig.tight_layout()
        out = os.path.join(HERE, "figs", "flat_plot" + ("" if geom == "single" else "_grid") + (f"_fs{fs:g}".replace(".", "p") if fs else "") + ".png")
        fig.savefig(out, dpi=150); plt.close(fig)
        print(f"wrote {out}")
    for r in sorted(rs, key=lambda r: (r["geom"], r["fs"], r["tag"], r["line"], r["x"])):
        print(f"  {r['geom']:6s} fs {r['fs']:<4g} {TREATMENT.get(r['tag'], r['tag']):42s} line {r['line']} x {r['x']:>3g}: FWHM {r['fwhm']:5.2f} A, "
              f"dip {r['dip']:+.2f}, side {r['side']:+.2f}{'' if r['clean'] else '  NOT CLEAN'}{'  (earlier run)' if r['early'] else ''}")


if __name__ == "__main__":
    main(*sys.argv[1:])
