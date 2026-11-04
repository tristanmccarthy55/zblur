"""
@file predict_groups.py
@brief Does a measured single-atom blur say when a crystal reconstruction succeeds, element by element?

One atom's numbers against the crystal's score under the same physics. The atoms: the crystal run's own kernels (a sparse
grid of single atoms through the same simulation and reconstruction, 18 slices of 1.53 A) and zblur's lone atoms (G1, G1b:
fine scan, 1 A slices). The score: phase1_scores.json, the campaign's official right-place numbers, as success = (right -
guess) / (1 - guess): 0 at guessing, 1 at perfect. Vibration runs are as run (sqrt(3) too warm, a8b9264): kernel and crystal
share the error, so the relation holds and no temperature may be read from it.

Measures of one atom, taken on the core column (centre disc minus a ring, per layer; radii in A so every kernel is cut alike):
  zshare    the share of the squared core phase further than c / 2 from the atom's own layer: how much has left its cell
  fwhm      full width at half maximum of the core column along z
  contrast  the core's peak over the same element's peak with nothing added
  S         sideways blur (A): the kz = 0 row of the transfer against the base's, fitted as a Gaussian (blobs.py's form)
Also: the lone atoms block-averaged to 1.53 A slices, in case the gap between lone atoms and grids is only the slicing; and the
zshare along the line a tilted column takes through the blob (superposition assumed, untested).

zblur's ARM atoms are NOT a measure of the ARM probe: zb_submit.sh sets no REGION_SIDE, so they ran in a 70 A periodic box with a
35 A solver window, which leaves 35.7 % of the 80 mrad probe outside (0.7 % at 60 mrad, 40.8 % at 90). They ring; the campaign's own
kernels (105 A window, 210 A region) do not (zshare 0.31 against 0.011, same grid, same 18 slices). Printed as lone_arm_window35 for
the record only. The ideal-probe atoms are small enough to be unaffected.

The threshold on zshare was picked after seeing these runs. It is fixed in runs/PREDICTIONS.md so that runs made after that
date can test it; do not change it here.

    ~/hyperspy-bundle/bin/python analysis/predict_groups.py
Writes results/predict_groups.tsv and results/predict_groups_tilt.tsv, and prints the tables.
"""
from __future__ import annotations

import csv
import glob
import itertools
import json
import os
import re
import sys

import numpy as np
from scipy.ndimage import shift as nd_shift
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))
import zprofile as zp  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCORES = os.path.join(HERE, "..", "ptychoshelves-clean", "aberration_experiment", "mono300", "phase1_scores.json")
GRID_DIRS = ["~/Desktop/mono_phase1", "~/Desktop/arm_block2", "~/Desktop/arm_split"]
LONE_DIRS = ["~/Desktop/zblur/g0d", "~/Desktop/zblur/g1", "~/Desktop/zblur/g1b"]
BLOBS = os.path.join(HERE, "results", "blobs_g0d_g1_g1b.tsv")

C_A = 3.905                 # the spacing along a column
DZ_GRID, DX_GRID = 27.525 / 18, 0.0739
R_CORE_A, R_RING_A = 0.148, (0.74, 1.48)      # the grid kernels' 2 px and 10-20 px at 0.0739 A
ZSHARE_WALL = 0.038         # fixed 2026-11-04, runs/PREDICTIONS.md
TILTS = (0.0, 0.0175, 0.035, 0.052, 0.07, 0.105)

# crystal run -> kernel tag. Left out: the old-engine mode legs (_mix, before 689b6d6) and the failed 2-stage search (_ps).
GRID_TAG = {"known_a080": "armf_a080", "split_coh": "armf_a080_coh", "split_ph": "armf_a080_ph16", "split_ph100": "armf_a080_ph16_T100",
            "T100_fixed": "armf_a080_ph16_T100_fs2p5_src0p1_mixF_fixed", "T20_fixed": "armf_a080_ph16_T20_fs2p5_src0p1_mixF_fixed",
            "T100_mixF": "armf_a080_ph16_T100_fs2p5_src0p1_mixF", "T20_mixF": "armf_a080_ph16_T20_fs2p5_src0p1_mixF",
            "fs2p5_src0p1_s3_ps": "armf_a080_fs2p5_src0p1_s3_ps", "T20_dose1e8_ps": "armf_a080_ph16_T20_fs2p5_src0p1_dose1e8_ps"}
# lone-atom physics -> the crystal run with the same physics, probe fixed
LONE_RUN = {"fs0p75": "fs0p75", "fs1p5": "fs1p5", "fs2p5": "fs2p5", "fs3p5": "fs3p5", "fs5": "fs5", "src0p05": "src0p05", "src0p1": "src0p1",
            "src0p2": "src0p2", "src0p4": "src0p4", "fs2p5_src0p1": "fs2p5_src0p1", "fs2p5_src0p2": "fs2p5_src0p2", "T295": "split_ph",
            "T100": "split_ph100", "fs2p5_src0p1_T20": "T20_fixed"}


def discs(K, dx, cx=None, cy=None):
    nz, ny, nx = K.shape
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(yy - ny // 2, xx - nx // 2) * dx
    return r <= R_CORE_A, (r >= R_RING_A[0]) & (r <= R_RING_A[1])


def core_column(K, dx):
    core, ring = discs(K, dx)
    return np.array([L[core].mean() - L[ring].mean() for L in K])


def zshare_of(p, dz, cut=C_A / 2):
    z = (np.arange(len(p)) + 0.5) * dz
    d = np.abs(z - z[int(np.argmax(p))])
    return float((p[d > cut] ** 2).sum() / (p ** 2).sum())


def measures(K, dx, dz):
    p = core_column(K, dx)
    z = (np.arange(len(p)) + 0.5) * dz
    return dict(peak=float(p.max()), zshare=zshare_of(p, dz), fwhm=zp.fwhm_interp(z, p))


def transfer(K, dz, dx):
    """As blobs.transfer: |FFT| of the Hann-windowed kernel, azimuthally averaged."""
    nz, ny, nx = K.shape
    w = np.hanning(nz)[:, None, None] * np.outer(np.hanning(ny), np.hanning(nx))[None]
    F = np.abs(np.fft.fftshift(np.fft.fftn((K - np.median(K, axis=(1, 2), keepdims=True)) * w)))
    kz = np.fft.fftshift(np.fft.fftfreq(nz, dz))
    ky, kx = np.meshgrid(np.fft.fftshift(np.fft.fftfreq(ny, dx)), np.fft.fftshift(np.fft.fftfreq(nx, dx)), indexing="ij")
    q = np.hypot(ky, kx)
    dq = 1.0 / (nx * dx)
    qb = np.arange(0, q.max(), dq)
    idx = np.digitize(q.ravel(), qb)
    n = np.maximum(np.bincount(idx, minlength=len(qb) + 1)[1:len(qb)], 1)
    T = np.array([np.bincount(idx, weights=F[k].ravel(), minlength=len(qb) + 1)[1:len(qb)] / n for k in range(nz)])
    return T, kz, qb[:T.shape[1]] + dq / 2


def sideways(K, base, dz, dx):
    """S (A): the kz = 0 row of K's transfer over the base's, log E = -2 pi^2 S^2 q^2, through the origin (blobs.fit_separable)."""
    T, kz, q = transfer(K, dz, dx)
    Tb, _, _ = transfer(base, dz, dx)
    k0 = int(np.argmin(np.abs(kz)))
    ok = Tb[k0] > 0.05 * Tb.max()
    e = np.where(ok, T[k0] / np.where(ok, Tb[k0], 1), np.nan)
    m = ok & (e > 0.02)
    if m.sum() < 3:
        return np.nan
    x, y = q[m] ** 2, np.log(e[m])
    return float(np.sqrt(max(-np.sum(x * y) / np.sum(x * x), 0) / (2 * np.pi ** 2)))


def tilted_zshare(K, dx, dz, theta):
    """zshare along the line a column tilted by theta takes through the blob, averaged over four directions."""
    core, ring = discs(K, dx)
    p0 = core_column(K, dx)
    z = (np.arange(len(p0)) - int(np.argmax(p0))) * dz
    out = []
    for a in (0, np.pi / 2, np.pi, 3 * np.pi / 2):
        p = []
        for k, L in enumerate(K):
            s = z[k] * np.tan(theta) / dx
            L2 = L if abs(s) < 1e-9 else nd_shift(L, (-s * np.sin(a), -s * np.cos(a)), order=1, mode="nearest")
            p.append(L2[core].mean() - L2[ring].mean())
        p = np.array(p)
        out.append((p[np.abs(z) > C_A / 2] ** 2).sum() / (p ** 2).sum())
    return float(np.mean(out))


def reslice(K, dz, dz_new):
    """Block-average a kernel along z to thicker layers (each new layer the mean of the old ones it covers)."""
    nz = K.shape[0]
    n_new = int(np.floor(nz * dz / dz_new))
    edges = np.arange(n_new + 1) * dz_new
    cum = np.concatenate([np.zeros((1,) + K.shape[1:]), np.cumsum(K, axis=0) * dz])
    z = np.arange(nz + 1) * dz
    i = np.clip(np.searchsorted(z, edges, side="right") - 1, 0, nz - 1)
    f = ((edges - z[i]) / dz)[:, None, None]
    C = cum[i] + f * (cum[i + 1] - cum[i])
    return np.diff(C, axis=0) / dz_new


def spear(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    return float(spearmanr(x[ok], y[ok]).correlation)


def best_cut(x, labels):
    """How many runs one threshold on x gets right against labels (True below the threshold), at its best."""
    return max(sum((a <= v) == b for a, b in zip(x, labels)) for v in sorted(set(x)))


def lone_key(name):
    bits = []
    for pat, lab in ((r"_fs(\d+(?:p\d+)?)", "fs"), (r"_src(\d+(?:p\d+)?)", "src"), (r"_T(\d+)", "T")):
        g = re.search(pat, name)
        if g:
            bits.append(lab + g.group(1))
    return "_".join(bits) or "base"


def main():
    S = {r["key"]: r for r in json.load(open(SCORES)) if r.get("found") and r.get("right")}
    succ = lambda k, e: (S[k]["right"][e] - S[k]["guess"][e]) / (1 - S[k]["guess"][e])
    for k in S:
        if k not in GRID_TAG and not k.endswith("_mix") and k != "fs2p5_src0p1_ps":
            GRID_TAG[k] = "armf_a080_" + k + ("" if k.endswith("_mixF") else "_fixed")
    files = {}
    for d in GRID_DIRS:
        for f in glob.glob(os.path.expanduser(d) + "/**/psf_*_vol.npy", recursive=True):
            files.setdefault(os.path.basename(f)[4:-8], f)
    base = {e: zp.load_phase(files[f"{e}_armf_a080"]) for e in ("Pb", "Ti")}
    b0 = {e: core_column(base[e], DX_GRID).max() for e in base}

    # ---------------------------------------------------------------- the crystal's own kernels
    rows = []
    for key, tag in GRID_TAG.items():
        if key not in S or f"Pb_{tag}" not in files or f"Ti_{tag}" not in files:
            continue
        r = dict(set="grid", run=key, modes=("mix" in key or key.endswith("_ps")), verdict=S[key].get("verdict"))
        for e in ("Pb", "Ti"):
            K = zp.load_phase(files[f"{e}_{tag}"])
            m = measures(K, DX_GRID, DZ_GRID)
            r.update({f"{e}_zshare": m["zshare"], f"{e}_fwhm": m["fwhm"], f"{e}_contrast": m["peak"] / b0[e],
                      f"{e}_S": sideways(K, base[e], DZ_GRID, DX_GRID), f"{e}_right": S[key]["right"][e], f"{e}_guess": S[key]["guess"][e],
                      f"{e}_success": succ(key, e)})
        rows.append(r)
    fx = [r for r in rows if not r["modes"]]
    md = [r for r in rows if r["modes"]]

    print(f"crystal's own kernels: {len(fx)} fixed-probe runs, {len(md)} mode runs")
    print("rank correlation with success, measure against the same element's own score:")
    for name, L in (("fixed", fx), ("modes", md), ("all", rows)):
        for e in ("Pb", "Ti"):
            print(f"  {name:6s} {e}: " + "  ".join(f"{m} {spear([r[f'{e}_{m}'] for r in L], [r[f'{e}_success'] for r in L]):+.2f}"
                                                  for m in ("zshare", "fwhm", "contrast", "S")))
    lab = [r["Pb_right"] >= 0.5 for r in fx]
    print(f"\nPb >= 50 % against one threshold, fixed probe: zshare {best_cut([r['Pb_zshare'] for r in fx], lab)} of {len(fx)}, "
          f"S {best_cut([r['Pb_S'] for r in fx], lab)}, fwhm {best_cut([r['Pb_fwhm'] for r in fx], lab)}, "
          f"1 - contrast {best_cut([1 - r['Pb_contrast'] for r in fx], lab)}")
    print(f"the fixed threshold zshare <= {ZSHARE_WALL} on fixed: "
          f"{sum((r['Pb_zshare'] <= ZSHARE_WALL) == (r['Pb_right'] >= .5) for r in fx)} of {len(fx)}; on modes: "
          f"{sum((r['Pb_zshare'] <= ZSHARE_WALL) == (r['Pb_right'] >= .5) for r in md)} of {len(md)}")
    for cut in (1.5, 1.95, 2.5, 3.0):
        v = [zshare_of(core_column(zp.load_phase(files["Pb_" + GRID_TAG[r["run"]]]), DX_GRID), DZ_GRID, cut) for r in fx]
        print(f"  zshare with the cut at {cut} A: rank correlation, fixed Pb, {spear(v, [r['Pb_success'] for r in fx]):+.2f}")
    # chance: the same threshold search on shuffled labels (one number), and a two-number box
    rng = np.random.default_rng(0)
    z_, s_ = [r["Pb_zshare"] for r in fx], [r["Pb_S"] for r in fx]
    n1 = best_cut(z_, lab)
    hit1 = sum(best_cut(z_, list(rng.permutation(lab))) >= n1 for _ in range(2000))
    print(f"  a threshold this good on shuffled labels: {hit1} of 2000 (shuffling keeps the 10 / 10 split only)")
    print("\nfixed probe, Pb, sorted by zshare:")
    for r in sorted(fx, key=lambda r: r["Pb_zshare"]):
        print(f"  {r['run']:16s} zshare {r['Pb_zshare']:.3f}  S {r['Pb_S']:.3f}  contrast {r['Pb_contrast']:.2f}  fwhm {r['Pb_fwhm']:.2f}"
              f"  Pb right {100 * r['Pb_right']:3.0f} guess {100 * r['Pb_guess']:3.0f}  | Ti zshare {r['Ti_zshare']:.3f} right {100 * r['Ti_right']:3.0f} guess {100 * r['Ti_guess']:3.0f}")
    print("modes, Pb, sorted by zshare:")
    for r in sorted(md, key=lambda r: r["Pb_zshare"]):
        print(f"  {r['run']:20s} zshare {r['Pb_zshare']:.3f}  S {r['Pb_S']:.3f}  Pb right {100 * r['Pb_right']:3.0f}  | Ti zshare {r['Ti_zshare']:.3f} right {100 * r['Ti_right']:3.0f}")

    # ---------------------------------------------------------------- lone atoms
    lone = {}
    dxs = {}
    for d in LONE_DIRS:
        d = os.path.expanduser(d)
        for man in glob.glob(d + "/manifest_*.tsv"):
            for m in csv.DictReader(open(man), delimiter="\t"):
                dxs[m["kernel"]] = float(m["dx_A"])
        for f in glob.glob(d + "/**/psf_Pb_a080_*NL20_vol.npy", recursive=True):
            n = os.path.basename(f)[4:-8]
            lone.setdefault(("arm" if "_arm" in n else "ideal", lone_key(n)), (f, n))
    blobs = {r["name"].split("  [")[0].strip(): r for r in csv.DictReader(open(BLOBS), delimiter="\t")}
    blob_name = {"fs0p75": "fs 0.75", "fs1p5": "fs 1.5", "fs2p5": "fs 2.5", "fs3p5": "fs 3.5", "fs5": "fs 5", "src0p05": "src 0.05", "src0p1": "src 0.1",
                 "src0p2": "src 0.2", "src0p4": "src 0.4", "fs2p5_src0p1": "fs 2.5, src 0.1", "fs2p5_src0p2": "fs 2.5, src 0.2", "T295": "295 K",
                 "T100": "100 K", "fs2p5_src0p1_T20": "fs 2.5, src 0.1, 20 K"}
    print("\nlone Pb atoms against the crystal's Pb success, fixed probe (radii in A as the grid kernels):")
    for probe in ("ideal", "arm"):
        L = []
        for phys, run in LONE_RUN.items():
            if (probe, phys) not in lone or run not in S:
                continue
            f, n = lone[(probe, phys)]
            K = zp.load_phase(f)
            dx = dxs.get(n, zp.DX)
            m = measures(K, dx, 1.0)
            K15 = reslice(K, 1.0, DZ_GRID)
            b = blobs.get(f"{probe} 80 mrad NL20, {blob_name[phys]}", {})
            r = dict(set="lone_ideal" if probe == "ideal" else "lone_arm_window35", run=phys, modes=False, verdict=S[run].get("verdict"), Pb_zshare=m["zshare"], Pb_fwhm=m["fwhm"],
                     Pb_contrast=float(b.get("contrast", "nan")), Pb_S=float(b.get("S", "nan")), Pb_right=S[run]["right"]["Pb"], Pb_guess=S[run]["guess"]["Pb"],
                     Pb_success=succ(run, "Pb"), zshare_1p53=zshare_of(core_column(K15, dx), DZ_GRID), dx_A=dx)
            L.append(r)
        rows += L
        y = [r["Pb_success"] for r in L]
        lab_l = [r["Pb_right"] >= 0.5 for r in L]
        print(f"  {'ideal' if probe == 'ideal' else 'arm, 35 A window, NOT a measure of the probe':5s} n={len(L)}: " + "  ".join(f"{m} {spear([r[f'Pb_{m}'] for r in L], y):+.2f}" for m in ("zshare", "fwhm", "contrast", "S"))
              + f"  zshare at 1.53 A slices {spear([r['zshare_1p53'] for r in L], y):+.2f}")
        print(f"        Pb >= 50 % against one threshold: zshare {best_cut([r['Pb_zshare'] for r in L], lab_l)} of {len(L)}, at 1.53 A slices "
              f"{best_cut([r['zshare_1p53'] for r in L], lab_l)}, S {best_cut([r['Pb_S'] for r in L], lab_l)}, 1 - contrast {best_cut([1 - r['Pb_contrast'] for r in L], lab_l)}")
        for r in sorted(L, key=lambda r: r["zshare_1p53"]):
            print(f"        {r['run']:18s} zshare {r['Pb_zshare']:.3f}  at 1.53 A {r['zshare_1p53']:.3f}  S {r['Pb_S']:.3f}  contrast {r['Pb_contrast']:.2f}  crystal Pb {100 * r['Pb_right']:3.0f}")

    # ---------------------------------------------------------------- tilt
    tilt_rows = []
    print("\nzshare along a tilted column, from the crystal's kernels (superposition assumed):")
    print("  run               element right% | " + " ".join(f"{1e3 * t:6.1f}" for t in TILTS) + "  mrad")
    for key in ("known_a080", "fs2p5_src0p1", "fs3p5_src0p1", "src0p2", "fs2p5_src0p2", "fs3p5_src0p2", "fs5", "fs7p29", "src0p4", "T20_fixed", "T100_fixed", "split_ph100"):
        for e in ("Pb", "Ti"):
            K = zp.load_phase(files[f"{e}_{GRID_TAG[key]}"])
            v = [tilted_zshare(K, DX_GRID, DZ_GRID, t) for t in TILTS]
            tilt_rows.append(dict(run=key, element=e, right=S[key]["right"][e], **{f"t{int(round(1e3 * t))}": x for t, x in zip(TILTS, v)}))
            print(f"  {key:17s} {e:2s} {100 * S[key]['right'][e]:5.0f} | " + " ".join(f"{x:6.3f}" for x in v))

    cols = ["set", "run", "modes", "verdict", "Pb_right", "Pb_guess", "Pb_success", "Pb_zshare", "Pb_fwhm", "Pb_contrast", "Pb_S", "zshare_1p53",
            "Ti_right", "Ti_guess", "Ti_success", "Ti_zshare", "Ti_fwhm", "Ti_contrast", "Ti_S"]
    with open(os.path.join(HERE, "results", "predict_groups.tsv"), "w") as f:
        w = csv.DictWriter(f, cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "results", "predict_groups_tilt.tsv"), "w") as f:
        w = csv.DictWriter(f, list(tilt_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(tilt_rows)


if __name__ == "__main__":
    main()
