"""
@file g0_report.py
@brief G0: does the single-atom depth PSF converge as the slices get finer, and do we reproduce Lei & Wang?

Reads an unpacked zblur pack (its manifest and psf/ kernels), measures every kernel's 3-D PSF (zprofile.psf3d), and
writes results/<pack>_table.tsv, figs/<pack>_profiles.png (the depth profiles through the atom's centre, normalised,
each run's slice counts overlaid, with Lei & Wang's 80 mrad FWHM of 2.8 A as a reference Gaussian) and
figs/<pack>_xz.png (x-z sections through the atom, as their SI Fig. 6).

    mkdir -p ~/Desktop/zblur/g0 && tar xzf ~/Desktop/zblur/zblur_g0_<ts>.tgz -C ~/Desktop/zblur/g0
    ~/hyperspy-bundle/bin/python analysis/g0_report.py ~/Desktop/zblur/g0 [~/Desktop/zblur/g0c ...]

Several packs are merged (a rerun's slices join the sim they came from). Per simulation with three or more slice
counts it fits FWHM_z^2 = W^2 + (c dz)^2: W is the core's own depth width, c dz what the slicing adds. The depth
response rings rather than falling like a Gaussian, so the side lobes are reported too: the deepest dip within 3 A
of the peak and the highest secondary lobe 2-6 A from it, both relative to the peak.
"""
from __future__ import annotations

import csv
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import zprofile as zp  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEI_WANG_FWHM_80 = 2.8     # A, a single Pr atom, 300 kV, 80 mrad, full pattern (arXiv 2506.22352 SI Note 5)


def lobes(z, p, zpk, peak):
    after = (np.abs(z - zpk) <= 3.0) & (np.abs(z - zpk) > 0.2)
    side = (np.abs(z - zpk) >= 2.0) & (np.abs(z - zpk) <= 6.0)
    return (float(p[after].min() / peak) if after.any() else np.nan, float(p[side].max() / peak) if side.any() else np.nan)


def load(pack):
    out = []
    for man in sorted(glob.glob(os.path.join(pack, "manifest_*.tsv"))):
        for r in csv.DictReader(open(man), delimiter="\t"):
            kern = r.get("kernel") or f"{r['sim']}_NL{r['nl']}"
            f = os.path.join(pack, "psf", f"psf_{kern}_vol.npy")
            if not os.path.exists(f):
                print(f"  missing kernel {os.path.basename(f)} (recon or extraction failed?)")
                continue
            K = zp.load_phase(f)
            dz, zdrop = float(r["dz_A"]), int(r["zdrop"])
            dx = float(r.get("dx_A") or zp.DX)
            m = zp.psf3d(K, dz, z0=zdrop * dz, dx=dx)
            dip, side = lobes(m["z"], m["p"], m["z_peak"], m["peak"])
            out.append(dict(r, **m, dz=dz, dx=dx, dip=dip, side=side, kern=kern, file=f))
    if not out and not glob.glob(os.path.join(pack, "manifest_*.tsv")):
        raise SystemExit(f"no manifest_*.tsv in {pack}")
    return out


def slicing_fit(rs):
    """Per sim: FWHM_z^2 = W^2 + (c dz)^2, a straight line in dz^2. -> {sim: (W, c, n)}."""
    fits = {}
    for s in dict.fromkeys(r["sim"] for r in rs):
        g = [r for r in rs if r["sim"] == s and np.isfinite(r["fwhm_z"])]
        if len(g) < 3:
            continue
        x = np.array([r["dz"] ** 2 for r in g]); y = np.array([r["fwhm_z"] ** 2 for r in g])
        c2, w2 = np.polyfit(x, y, 1)
        fits[s] = (float(np.sqrt(max(w2, 0))), float(np.sqrt(max(c2, 0))), len(g))
    return fits


def main(*packs):
    rs = [r for p in packs for r in load(p)]
    tag = "_".join(os.path.basename(os.path.normpath(p)) for p in packs)
    if not rs:
        raise SystemExit("no kernels")
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", f"{tag}_table.tsv"), "w") as fh:
        fh.write("kernel\tsim\tnl\tdz_A\tdx_A\tpeak_rad\tz_peak_A\tfwhm_z_A\tsigma_z_A\tcore_share\tdip\tside_lobe\tfwhm_xy_A\n")
        for r in rs:
            fh.write(f"{r['kern']}\t{r['sim']}\t{r['nl']}\t{r['dz']:.3f}\t{r['dx']:.4f}\t{r['peak']:.3f}\t{r['z_peak']:.2f}\t"
                     f"{r['fwhm_z']:.2f}\t{r['sigma_z']:.2f}\t{r['core_share']:.2f}\t{r['dip']:.2f}\t{r['side']:.2f}\t{r['fwhm_xy']:.3f}\n")
    print(f"{'kernel':<66} {'dz':>5} {'FWHMz':>6} {'core':>5} {'dip':>6} {'side':>5} {'FWHMxy':>7}")
    for r in rs:
        print(f"{r['kern']:<66} {r['dz']:5.2f} {r['fwhm_z']:6.2f} {r['core_share']:5.2f} {r['dip']:6.2f} {r['side']:5.2f} {r['fwhm_xy']:7.3f}")
    for s, (w, c, n) in slicing_fit(rs).items():
        print(f"slicing fit {s}: FWHM_z^2 = W^2 + (c dz)^2 over {n} slicings -> W = {w:.2f} A, c = {c:.2f}")
    plot(rs, tag)


def plot(rs, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sims = list(dict.fromkeys(r["sim"] for r in rs))
    fig, axs = plt.subplots(1, len(sims), figsize=(5.2 * len(sims), 3.8), squeeze=False)
    for ax, s in zip(axs[0], sims):
        group = sorted([r for r in rs if r["sim"] == s], key=lambda r: -r["dz"])
        cm = plt.cm.viridis(np.linspace(0.15, 0.85, len(group)))
        for c, r in zip(cm, group):
            ax.plot(r["z"] - r["z_peak"], r["p"] / r["peak"], "-o", ms=2.5, lw=1.4, color=c,
                    label=f"dz {r['dz']:.2f} Å: FWHM {r['fwhm_z']:.2f} Å")
        zz = np.linspace(-10, 10, 400); sg = LEI_WANG_FWHM_80 / 2.3548
        if group[0]["alpha"] == "80":
            ax.plot(zz, np.exp(-zz ** 2 / (2 * sg ** 2)), "--", color="0.5", lw=1, label="Lei & Wang 80 mrad, 2.8 Å")
        ax.axhline(0.5, color="0.8", lw=0.8)
        ax.set_xlim(-10, 10); ax.set_xlabel("depth from the peak (Å)"); ax.set_ylabel("phase / peak")
        ax.set_title(s.replace("_", " "), fontsize=8, loc="left"); ax.legend(fontsize=7, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figs", f"{tag}_profiles.png"), dpi=140)
    fig, axs = plt.subplots(1, len(rs), figsize=(2.2 * len(rs), 4.2), squeeze=False)
    for ax, r in zip(axs[0], rs):
        ext = [r["x"][0], r["x"][-1], r["z"][-1] + r["dz"] / 2, r["z"][0] - r["dz"] / 2]
        ax.imshow(r["xz"] / r["peak"], extent=ext, aspect="auto", cmap="viridis", vmin=0, vmax=1, interpolation="nearest")
        ax.set_title(f"NL {r['nl']}, dz {r['dz']:.2f}", fontsize=8); ax.set_xlabel("x (Å)")
    axs[0][0].set_ylabel("depth z (Å)")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figs", f"{tag}_xz.png"), dpi=140)
    print("wrote figs", f"{tag}_profiles.png", f"{tag}_xz.png")


if __name__ == "__main__":
    main(*(sys.argv[1:] or [os.path.expanduser("~/Desktop/zblur/g0")]))
