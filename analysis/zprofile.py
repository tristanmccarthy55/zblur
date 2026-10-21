"""
@file zprofile.py
@brief The depth (z) profile of a reconstructed single atom, and its width from a slice-averaged Gaussian fit.

A kernel volume is what ptychoshelves-clean's analysis/atomfind/extract_psf.py writes: exp(i * K), K the averaged phase
crop of the grid atoms, layers [zdrop, NL - zdrop) of an NL-layer reconstruction of a box BOXZ deep (the vacuum bands
dropped), +-W px in-plane around the atom. The depth profile here is the phase integrated over a disc around the atom,
minus the per-layer background of a ring outside it, so a lateral blur (source size) does not change it, only a depth
blur does. Each reconstruction layer holds the average of the response over its thickness dz, so the model is a
Gaussian AVERAGED OVER EACH SLICE plus a constant: its sigma is the underlying depth width, not inflated by the slicing.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf

BOXZ = 27.525       # the campaign's box: 5 cells + 2 x 4 A vacuum
ZVAC = 4.0


def load_phase(path):
    v = np.load(path)
    return np.angle(v) if np.iscomplexobj(v) else np.asarray(v, float)


def profile(K, r_core=2, r_bg=(10, 20)):
    """The atom's CORE phase per layer: the mean over r <= r_core px at the crop centre minus the mean of a ring r_bg.

    Not the laterally integrated phase: that is the zero-frequency part of each layer, which multislice ptychography does
    not resolve in depth at all (the missing cone), so it comes out spread over the whole box. The depth information is
    carried by the atom's sharp core, so the core minus its local surroundings is what is fitted."""
    nz, ny, nx = K.shape
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(yy - ny // 2, xx - nx // 2)
    core, ring = r <= r_core, (r >= r_bg[0]) & (r <= r_bg[1])
    return np.array([K[k][core].mean() - K[k][ring].mean() for k in range(nz)])


def slice_gauss(z, A, mu, sig, c, dz):
    """A Gaussian of area A, centre mu, width sig, averaged over slices of thickness dz centred on z, plus c."""
    s2 = np.sqrt(2) * sig
    return A * (erf((z + dz / 2 - mu) / s2) - erf((z - dz / 2 - mu) / s2)) / (2 * dz) + c


def fit(p, dz, z0=0.0):
    """Fit slice_gauss to a profile p sampled at layer centres z0 + (k + 0.5) dz. -> dict(sigma, fwhm, mu, err...)."""
    z = z0 + (np.arange(len(p)) + 0.5) * dz
    k = int(np.argmax(p))
    A0 = max(p[k] - np.median(p), 1e-6) * dz * 1.5
    f = lambda zz, A, mu, sig, c: slice_gauss(zz, A, mu, sig, c, dz)
    try:
        popt, pcov = curve_fit(f, z, p, p0=[A0, z[k], max(dz, 1.0), np.median(p)],
                               bounds=([0, z[0] - dz, 0.05, -np.inf], [np.inf, z[-1] + dz, 30.0, np.inf]), maxfev=20000)
    except RuntimeError:
        return None
    err = np.sqrt(np.diag(pcov))
    resid = p - f(z, *popt)
    return dict(A=popt[0], mu=popt[1], sigma=popt[2], c=popt[3], sigma_err=err[2], fwhm=2.3548 * popt[2],
                rms_resid=float(np.sqrt(np.mean(resid ** 2)) / max(p.max() - popt[3], 1e-9)), z=z, p=p,
                model=lambda zz: f(zz, *popt))


def kernel_width(path, nl, boxz=BOXZ, zvac=ZVAC):
    """sigma_z of one kernel file, given the NL of the reconstruction it came from."""
    dz = boxz / nl
    zdrop = int(round(zvac / dz))
    K = load_phase(path)
    if K.shape[0] != nl - 2 * zdrop:
        raise ValueError(f"{path}: {K.shape[0]} layers, expected NL {nl} - 2 x {zdrop}")
    out = fit(profile(K), dz, z0=zdrop * dz)
    if out is not None:
        out.update(dz=dz, nl=nl, layers=K.shape[0])
    return out
