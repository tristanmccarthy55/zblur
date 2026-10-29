"""
@file blobs.py
@brief The blob model: each single-atom kernel as a few numbers, each blur source as an envelope, and the tilt at which
one atom's blob stops reaching the next atom's core.

For every kernel in the packs given (G0d / G1 / G1b, lone atoms):
  core     contrast (the core's peak phase, relative to its base), sigma_xy (2-D Gaussian in the peak layer), sigma_z
           (slice-averaged Gaussian on the centre column, +-3 A)
  lobes    the deepest dip within 3 A of the peak and the highest side lobe 2-6 A from it, relative to the peak
  cone     the hourglass: its opening (the radius of the strongest ring against |dz|, a straight-line fit) and its height
           one lattice spacing away (the ring's phase at |dz| = SPACING, relative to the core)
  envelope E(q, kz) = |T_source| / |T_base| (T = FFT of the windowed kernel, averaged over azimuth), fitted to the
           predicted form of its source: exp(-2 pi^2 d^2 kz^2) (focal spread d), exp(-2 pi^2 s^2 q^2) (source
           s = FWHM / 2.355), exp(-2 pi^2 u^2 (q^2 + kz^2)) (phonons) -- each returns an effective width to compare
           with the true one. A combination's envelope is compared with the product of its single sources'.
  tilt     neighbours along a column spacing SPACING, tilted by theta, sit at d = SPACING (sin theta, 0, cos theta)
           from each other. C(theta) = the blob's mean positive phase over the neighbour's core (r <= 2 px, both
           neighbours) over the blob's own core: how much of one atom lands on the next. theta_10 / theta_5 = where it
           falls below 10 / 5 %. This treats the reconstruction as a sum of blobs; G2 tests that on real columns.

The base of a kernel is the run with the same probe, alpha, scan and box and no blur source (G0d's two at 80 mrad).

    ~/hyperspy-bundle/bin/python analysis/blobs.py ~/Desktop/zblur/g0d ~/Desktop/zblur/g1 ~/Desktop/zblur/g1b
"""
from __future__ import annotations

import csv
import glob
import os
import sys

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(__file__))
import zprofile as zp  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPACING = 3.905           # A, the perovskite cell along the column (PbTiO3 / SrTiO3 a)
THETAS = np.linspace(0, 0.15, 61)


# ------------------------------------------------------------------------------------------------------ loading
def load(packs):
    out = []
    for pack in packs:
        for man in sorted(glob.glob(os.path.join(pack, "manifest_*.tsv"))):
            for r in csv.DictReader(open(man), delimiter="\t"):
                if r.get("geom") != "single":
                    continue
                f = os.path.join(pack, "psf", f"psf_{r['kernel']}_vol.npy")
                if not os.path.exists(f):
                    print(f"  missing {os.path.basename(f)}")
                    continue
                r = dict(r, file=f, K=zp.load_phase(f), dz=float(r["dz_A"]), dx=float(r["dx_A"]),
                         fs=float(r["fs_A"]), src=float(r["src_A"]), T=(None if r["T_K"] in ("-", "") else float(r["T_K"])))
                r["z0"] = int(r["zdrop"]) * r["dz"]
                out.append(r)
    return out


def base_key(r):
    return (r["aber"], r["alpha"], r["step_A"], r["win_A"], r["boxz_A"], r["nl"])


# ------------------------------------------------------------------------------------------------------ core and lobes
def core(r):
    K, dz, dx = r["K"], r["dz"], r["dx"]
    m = zp.psf3d(K, dz, z0=r["z0"], dx=dx, core_A=3.0)
    k = int(np.argmax(m["p"]))
    nz, ny, nx = K.shape
    W = int(round(0.4 / dx))
    yy, xx = np.mgrid[-W:W + 1, -W:W + 1] * dx
    patch = K[k, ny // 2 - W:ny // 2 + W + 1, nx // 2 - W:nx // 2 + W + 1] - np.median(K[k])
    g2 = lambda X, A, s, c: (A * np.exp(-(X[0] ** 2 + X[1] ** 2) / (2 * s ** 2)) + c).ravel()
    try:
        (A, sxy, c), _ = curve_fit(g2, (xx, yy), patch.ravel(), p0=[patch.max(), 0.06, 0.0], bounds=([0, 0.005, -np.inf], [np.inf, 1.0, np.inf]))
    except RuntimeError:
        sxy = np.nan
    z = m["z"]; p = m["p"]; zpk = m["z_peak"]; pk = m["peak"]
    after = (np.abs(z - zpk) <= 3.0) & (np.abs(z - zpk) > 0.2)
    side = (np.abs(z - zpk) >= 2.0) & (np.abs(z - zpk) <= 6.0)
    # sigma_z from the FWHM read straight off the centre column (/2.355): a Gaussian fit collapses wherever a dip sits
    # inside its window (the ARM base: -0.22 at +1 A), and the FWHM is what the tilt and spacing comparisons need
    return dict(peak=pk, z_peak=zpk, k_peak=k, sigma_z=m["fwhm_z"] / 2.3548, fwhm_z=m["fwhm_z"], sigma_xy=float(sxy),
                dip=float(p[after].min() / pk) if after.any() else np.nan, side=float(p[side].max() / pk) if side.any() else np.nan)


# ------------------------------------------------------------------------------------------------------ the cone
def cone(r, c):
    """Azimuthally averaged phase per layer against radius; the ring's radius against |dz| -> the opening; its height at
    |dz| = SPACING (interpolated between layers) relative to the core."""
    K, dz, dx = r["K"], r["dz"], r["dx"]
    nz, ny, nx = K.shape
    yy, xx = np.mgrid[:ny, :nx]
    rad = np.hypot(yy - ny // 2, xx - nx // 2) * dx
    edges = np.arange(0, rad.max(), dx)
    idx = np.digitize(rad.ravel(), edges)
    prof = np.array([np.bincount(idx, weights=(K[k] - np.median(K[k])).ravel(), minlength=len(edges) + 1)[1:len(edges)]
                     / np.maximum(np.bincount(idx, minlength=len(edges) + 1)[1:len(edges)], 1) for k in range(nz)])
    rc = edges[:prof.shape[1]] + dx / 2
    z = r["z0"] + (np.arange(nz) + 0.5) * dz
    dzs = z - c["z_peak"]
    ring_r, ring_h, use = [], [], []
    for k in range(nz):
        if abs(dzs[k]) < 1.5:
            continue
        sel = rc >= 0.15                                  # outside the core
        j = int(np.argmax(prof[k][sel]))
        ring_r.append(rc[sel][j]); ring_h.append(prof[k][sel][j] / c["peak"]); use.append(dzs[k])
    use, ring_r, ring_h = np.array(use), np.array(ring_r), np.array(ring_h)
    opening = float(np.sum(np.abs(use) * ring_r) / np.sum(use ** 2)) if len(use) else np.nan     # r = tan(theta) |dz|
    h_at = {}
    for sgn, lab in ((1, "below"), (-1, "above")):
        s = use * sgn > 0
        h_at[lab] = float(np.interp(SPACING, np.abs(use[s]), ring_h[s])) if s.sum() >= 2 else np.nan
    return dict(opening=opening, cone_below=h_at["below"], cone_above=h_at["above"])


# ------------------------------------------------------------------------------------------------------ envelopes
def transfer(K, dz, dx):
    """|FFT| of the kernel, windowed (Hann in z and in-plane), azimuthally averaged -> T[kz, q], kz (1/A), q (1/A)."""
    nz, ny, nx = K.shape
    w = np.hanning(nz)[:, None, None] * np.outer(np.hanning(ny), np.hanning(nx))[None]
    F = np.abs(np.fft.fftshift(np.fft.fftn((K - np.median(K, axis=(1, 2), keepdims=True)) * w)))
    kz = np.fft.fftshift(np.fft.fftfreq(nz, dz))
    ky, kx = np.meshgrid(np.fft.fftshift(np.fft.fftfreq(ny, dx)), np.fft.fftshift(np.fft.fftfreq(nx, dx)), indexing="ij")
    q = np.hypot(ky, kx)
    dq = 1.0 / (nx * dx)
    qb = np.arange(0, q.max(), dq)
    idx = np.digitize(q.ravel(), qb)
    T = np.array([np.bincount(idx, weights=F[k].ravel(), minlength=len(qb) + 1)[1:len(qb)]
                  / np.maximum(np.bincount(idx, minlength=len(qb) + 1)[1:len(qb)], 1) for k in range(nz)])
    return T, kz, qb[:T.shape[1]] + dq / 2


def envelope(r, b):
    T, kz, q = transfer(r["K"], r["dz"], r["dx"])
    Tb, _, _ = transfer(b["K"], b["dz"], b["dx"])
    ok = Tb > 0.05 * Tb.max()
    E = np.where(ok, T / np.where(ok, Tb, 1), np.nan)
    KZ, Q = np.meshgrid(kz, q, indexing="ij")
    return E, KZ, Q, ok


def fit_envelope(E, KZ, Q, ok, kind):
    """The effective width of the predicted form for this kind of source (A), by least squares on log E."""
    y = np.log(np.clip(E[ok], 1e-6, None))
    x = {"fs": KZ[ok] ** 2, "src": Q[ok] ** 2, "T": KZ[ok] ** 2 + Q[ok] ** 2}[kind]
    a = float(np.sum(x * y) / np.sum(x * x))           # log E = -2 pi^2 w^2 x, through the origin
    return float(np.sqrt(max(-a, 0) / (2 * np.pi ** 2)))


# ------------------------------------------------------------------------------------------------------ tilt
def tilt_curve(r, c):
    """C(theta): the blob's positive phase averaged over a neighbour's core (r <= 2 px, one layer), both neighbours,
    relative to the blob's own core value."""
    K = np.clip(r["K"] - np.median(r["K"], axis=(1, 2), keepdims=True), 0, None)
    nz, ny, nx = K.shape
    dz, dx = r["dz"], r["dx"]
    k0, y0, x0 = c["k_peak"], ny // 2, nx // 2
    disc = [(a, b) for a in range(-2, 3) for b in range(-2, 3) if a * a + b * b <= 4]
    own = np.mean([K[k0, y0 + a, x0 + b] for a, b in disc])
    out = []
    for th in THETAS:
        vals = []
        for sgn in (1, -1):
            kz = k0 + sgn * SPACING * np.cos(th) / dz
            xs = x0 + sgn * SPACING * np.sin(th) / dx
            pts = np.array([[kz] * len(disc), [y0 + a for a, _ in disc], [xs + b for _, b in disc]])
            vals.append(map_coordinates(K, pts, order=1, mode="constant", cval=0.0).mean())
        out.append(np.mean(vals) / own)
    return np.array(out)


def crossing(c, level):
    below = np.where(c < level)[0]
    return float(THETAS[below[0]]) if len(below) else np.nan


# ------------------------------------------------------------------------------------------------------ main
def describe(r):
    bits = [f"{r['aber']} {int(r['alpha'])} mrad NL{r['nl']}"]
    if r["fs"]:
        bits.append(f"fs {r['fs']:g}")
    if r["src"]:
        bits.append(f"src {r['src']:g}")
    if r["T"]:
        bits.append(f"{r['T']:g} K")
    return ", ".join(bits)


def main(*packs):
    rs = load(packs)
    if not rs:
        raise SystemExit("no lone-atom kernels")
    bases = {base_key(r): r for r in rs if not r["fs"] and not r["src"] and r["T"] is None}
    rows = []
    for r in rs:
        c = core(r); cn = cone(r, c); tc = tilt_curve(r, c)
        b = bases.get(base_key(r))
        # a kernel that is noise (0.5 A slices: the engine's floor) has a peak no stronger than its own lobes
        clean = bool(np.isfinite(c["fwhm_z"]) and c["side"] < 0.6 and c["dip"] > -0.6)
        row = dict(name=describe(r) + ("" if clean else "  [NOT CLEAN]"), kernel=r["kernel"], clean=clean, **c, **cn,
                   contrast=np.nan, w_fs=np.nan, w_src=np.nan, w_T=np.nan,
                   theta10=crossing(tc, 0.10), theta5=crossing(tc, 0.05), C0=float(tc[0]), tilt=tc, r=r)
        if b is not None and b is not r:
            row["contrast"] = c["peak"] / core(b)["peak"]
            E, KZ, Q, ok = envelope(r, b)
            for kind in ("fs", "src", "T"):
                row[f"w_{kind}"] = fit_envelope(E, KZ, Q, ok, kind)
        rows.append(row)
    write(rows, packs)


def write(rows, packs):
    tag = "_".join(os.path.basename(os.path.normpath(p)) for p in packs)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    cols = ["name", "contrast", "sigma_xy", "sigma_z", "fwhm_z", "dip", "side", "opening", "cone_below", "cone_above",
            "w_fs", "w_src", "w_T", "C0", "theta10", "theta5", "kernel"]
    with open(os.path.join(HERE, "results", f"blobs_{tag}.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for w in rows:
            fh.write("\t".join(f"{w[k]:.4f}" if isinstance(w[k], float) else str(w[k]) for k in cols) + "\n")
    print(f"{'blob':<34} {'contr':>5} {'s_xy':>5} {'s_z':>5} {'dip':>5} {'side':>5} {'open':>5} {'cone@c':>6} "
          f"{'w_fs':>5} {'w_src':>5} {'w_T':>5} {'C(0)':>5} {'th10':>5} {'th5':>5}")
    for w in sorted(rows, key=lambda w: w["name"]):
        f = lambda v, d=2: "  -  " if not np.isfinite(v) else f"{v:5.{d}f}"
        print(f"{w['name']:<34} {f(w['contrast'])} {f(w['sigma_xy'], 3)} {f(w['sigma_z'])} {f(w['dip'])} {f(w['side'])} "
              f"{f(w['opening'], 3)} {f(np.nanmean([w['cone_below'], w['cone_above']])):>6} {f(w['w_fs'])} {f(w['w_src'], 3)} "
              f"{f(w['w_T'], 3)} {f(w['C0'])} {f(w['theta10'], 3)} {f(w['theta5'], 3)}")
    plot_tilt(rows, tag)


def plot_tilt(rows, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for w in sorted([w for w in rows if w["clean"]], key=lambda w: w["name"]):
        ax.plot(THETAS * 1000, w["tilt"] * 100, lw=1.4, label=w["name"])
    ax.axhline(10, color="0.6", ls=":", lw=1); ax.axhline(5, color="0.6", ls=":", lw=1)
    ax.set_xlabel("tilt off the zone axis (mrad)")
    ax.set_ylabel("one blob over the next atom's core (%)")
    ax.set_title(f"neighbours {SPACING} A apart along the column", fontsize=10, loc="left")
    ax.legend(fontsize=7, frameon=False); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figs", f"blobs_tilt_{tag}.png"), dpi=140)
    print("wrote", f"results/blobs_{tag}.tsv", f"figs/blobs_tilt_{tag}.png")


if __name__ == "__main__":
    main(*(sys.argv[1:] or [os.path.expanduser("~/Desktop/zblur/g0d")]))
