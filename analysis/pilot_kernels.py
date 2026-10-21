"""
@file pilot_kernels.py
@brief Phase 0: the depth response of a single Pb / Ti atom under each blur source, from the campaign's EXISTING kernels.

Every run of the ARM200F and 300 kV mono campaigns was scored with Pb and Ti kernels: grids of one element on one depth
plane, simulated with the same physics as the run (aberrations, partial coherence, phonons) and reconstructed the same
way (true probe fixed, or probe modes). Those kernels are single-atom point-spread functions, already made, at the
campaign's 18 slices of 1.53 A (80 mrad). This reads them all and reports, per condition:

  peak   the core phase in the atom's own layer (rad): how much contrast the blur leaves
  conc   the share of the positive core phase within +-1.5 A of the true depth: how much stays at the atom
  sig_m  the second moment of the positive core phase about its centroid (A): a model-free depth width
  sig_g  a slice-averaged Gaussian fit (A); only meaningful where the blur spans several slices

At 1.53 A slices the clean response sits in one or two layers, so this is a coarse first look; the fine-slice runs
(runs/PLAN.md) are the measurement. Writes figs/pilot_profiles.png and results/pilot_table.tsv.

    ~/hyperspy-bundle/bin/python analysis/pilot_kernels.py
"""
from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import zprofile as zp  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEARCH = ["~/Desktop/mono_phase1", "~/Desktop/arm_block2", "~/Desktop/arm_split"]
ATOM_Z = zp.BOXZ / 2
NL_OF_ALPHA = {40: 4, 50: 7, 60: 10, 65: 12, 70: 14, 75: 16, 80: 18, 90: 23}

# tag -> (group, label, chromatic focal spread A rms, source A FWHM, temperature K or None, probe)
COH_FEG = (7.29, 0.4)


def describe(tag):
    m = re.match(r"armf_a(\d+)(?:_(.*))?$", tag)
    alpha, rest = int(m.group(1)), (m.group(2) or "")
    fs = src = 0.0
    T = None
    probe = "fixed"
    g = re.search(r"fs(\d+(?:p\d+)?)", rest)
    if g:
        fs = float(g.group(1).replace("p", "."))
    g = re.search(r"src(\d+(?:p\d+)?)", rest)
    if g:
        src = float(g.group(1).replace("p", "."))
    if re.search(r"(^|_)coh(_|$)", rest):
        fs, src = COH_FEG
    if "ph16" in rest:
        g = re.search(r"_T(\d+)", rest)
        T = int(g.group(1)) if g else 295
    if re.search(r"_mix", rest) and not rest.endswith("_fixed"):
        probe = "modes"
    if rest.endswith("mixF_fixed"):
        probe = "fixed"
    return dict(alpha=alpha, fs=fs, src=src, T=T, probe=probe)


def find_kernels(el="Pb"):
    out = {}
    for d in SEARCH:
        for f in glob.glob(os.path.join(os.path.expanduser(d), "analysis", "*", "psf", f"psf_{el}_*_vol.npy")):
            tag = os.path.basename(f)[len(f"psf_{el}_"):-len("_vol.npy")]
            if re.match(r"armf_a\d+(_|$)", tag) and "w140" not in tag:
                out.setdefault(tag, f)
    return out


def metrics(path, nl):
    dz = zp.BOXZ / nl
    zdrop = int(round(zp.ZVAC / dz))
    K = zp.load_phase(path)
    p = zp.profile(K)
    z = (zdrop + np.arange(len(p)) + 0.5) * dz
    pp = np.clip(p, 0, None)
    near = np.abs(z - ATOM_Z) <= 1.5
    mu = float((pp * z).sum() / pp.sum())
    sig_m = float(np.sqrt((pp * (z - mu) ** 2).sum() / pp.sum()))
    g = zp.fit(p, dz, z0=zdrop * dz)
    return dict(dz=dz, z=z, p=p, peak=float(p.max()), conc=float(pp[near].sum() / pp.sum()), sig_m=sig_m, mu=mu,
                sig_g=(g["sigma"] if g else np.nan), sig_g_err=(g["sigma_err"] if g else np.nan))


def name(d):
    bits = []
    if d["fs"]:
        bits.append(f"fs {d['fs']:g}")
    if d["src"]:
        bits.append(f"src {d['src']:g}")
    if d["T"]:
        bits.append(f"{d['T']} K")
    s = " + ".join(bits) or "coherent"
    return s + (" (modes)" if d["probe"] == "modes" else "")


def main():
    rows = []
    for el in ("Pb", "Ti"):
        for tag, f in sorted(find_kernels(el).items()):
            d = describe(tag)
            nl = NL_OF_ALPHA.get(d["alpha"])
            if nl is None or nl < 14:
                continue                                    # too few layers to say anything about depth
            m = metrics(f, nl)
            rows.append(dict(el=el, tag=tag, name=name(d), **d, **m))
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", "pilot_table.tsv"), "w") as fh:
        fh.write("el\talpha\tname\tprobe\tfs_A\tsrc_A\tT_K\tdz_A\tpeak_rad\tconc\tsig_m_A\tsig_g_A\ttag\n")
        for r in rows:
            fh.write(f"{r['el']}\t{r['alpha']}\t{r['name']}\t{r['probe']}\t{r['fs']:g}\t{r['src']:g}\t{r['T'] or ''}\t"
                     f"{r['dz']:.3f}\t{r['peak']:.3f}\t{r['conc']:.2f}\t{r['sig_m']:.2f}\t{r['sig_g']:.2f}\t{r['tag']}\n")
    for el in ("Pb", "Ti"):
        print(f"\n{el}: alpha  condition                          peak   conc  sig_m  sig_g")
        for r in sorted([r for r in rows if r["el"] == el], key=lambda r: (r["alpha"], r["T"] or 0, r["fs"], r["src"], r["probe"])):
            print(f"   {r['alpha']:>3}  {r['name']:<34} {r['peak']:6.3f} {r['conc']:5.2f} {r['sig_m']:6.2f} {r['sig_g']:6.2f}")
    plot(rows)


def plot(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    groups = [("chromatic focal spread", lambda r: r["T"] is None and r["src"] == 0 and r["probe"] == "fixed" and r["alpha"] == 80),
              ("source size", lambda r: r["T"] is None and r["fs"] == 0 and r["probe"] == "fixed" and r["alpha"] == 80),
              ("both at once", lambda r: r["T"] is None and r["fs"] > 0 and r["src"] > 0 and r["probe"] == "fixed" and r["alpha"] == 80),
              ("phonons", lambda r: r["T"] is not None and r["probe"] == "fixed"),
              ("probe modes", lambda r: r["probe"] == "modes")]
    fig, axs = plt.subplots(2, len(groups), figsize=(3.2 * len(groups), 6.4), sharey="row", squeeze=False)
    ref = {el: next(r for r in rows if r["el"] == el and r["tag"] == "armf_a080") for el in ("Pb", "Ti")}
    for i, el in enumerate(("Pb", "Ti")):
        for j, (t, sel) in enumerate(groups):
            ax = axs[i, j]
            ax.step(ref[el]["z"], ref[el]["p"], where="mid", color="k", lw=1.6, label="coherent")
            rs = sorted([r for r in rows if r["el"] == el and sel(r) and r["tag"] != "armf_a080"],
                        key=lambda r: (r["T"] or 0, r["fs"], r["src"]))
            cm = plt.cm.viridis(np.linspace(0.1, 0.9, max(len(rs), 1)))
            for c, r in zip(cm, rs):
                ax.step(r["z"], r["p"], where="mid", color=c, lw=1.2, label=r["name"])
            ax.axvline(ATOM_Z, color="0.6", ls=":", lw=1)
            ax.set_title(f"{el}: {t}", fontsize=10, loc="left")
            ax.set_xlabel("depth z (Å)")
            if j == 0:
                ax.set_ylabel("core phase (rad)")
            ax.legend(fontsize=6.5, frameon=False)
            ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = os.path.join(HERE, "figs", "pilot_profiles.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
