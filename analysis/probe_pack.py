#!/usr/bin/env python
"""
@file probe_pack.py
@brief Pack the true, start and fitted probes of released-probe runs into small files, for the "ptychography fitted the phase" panel.

Runs on Blythe beside the reconstructions (numpy, scipy, h5py; read-only on the recon dirs). For each recon dir it reads, as
ptychoshelves-clean/analysis/kick_probe.py does (the same frame for all three, checked there 2026-09-23):
  true    <recon>/01/probe_initial_true.mat     the probe the data were simulated with
  start   <recon>/01/probe_initial.mat          the experimentalist's start (the search's errors), when the run wrote one
  fitted  the dominant mode of reconstruction/probes in the newest analysis/**/*_recons.h5
and writes <out>/<recon name>.npz: each probe in the aperture plane (centred FFT), cropped to +-(convergence + 10) mrad, complex64;
|probe|^2 in real space cropped to +-REAL_HALF_A; d_alpha_mrad, the convergence angle, the object pixel; and the shift- and
phase-invariant overlaps with the truth (kick_probe's measure). A few MB per run.

    $SHARE/phucrh/envs/abtem/bin/python analysis/probe_pack.py --out <dir> <recon dir> [<recon dir> ...]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

REAL_HALF_A = 3.0
MARGIN_MRAD = 10.0


def modes_of(P):
    """(n_modes, N, N) from the h5's probes array, whichever axis carries the modes (as kick_probe.modes_of)."""
    P = np.asarray(P).squeeze()
    if P.ndim == 2:
        return P[None]
    return np.moveaxis(P, -1, 0) if P.shape[-1] < P.shape[0] else P


def overlap(a, b):
    """|<a, b>| maximised over a real-space shift and a global phase, normalised (kick_probe.shift_invariant_overlap)."""
    x = np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(b)))
    return float(np.abs(x).max() / (np.linalg.norm(a) * np.linalg.norm(b)))


def aperture_plane(P):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(P)))


def crop(A, half):
    c = A.shape[0] // 2
    return A[c - half:c + half, c - half:c + half]


def pack(recon_dir, out):
    from scipy.io import loadmat
    import h5py
    d01 = os.path.join(recon_dir, "01")
    meta = loadmat(os.path.join(d01, "sim_meta.mat"), simplify_cells=True)["meta"]
    da = float(meta["d_alpha_rad"]) * 1e3
    conv = float(meta["convergence_mrad"])
    T = np.asarray(loadmat(os.path.join(d01, "probe_initial_true.mat"))["probe"]).astype(np.complex128).squeeze()
    S = None
    sp = os.path.join(d01, "probe_initial.mat")
    if os.path.exists(sp) and not os.path.islink(sp):          # a written start; a link means the true probe was used
        S = np.asarray(loadmat(sp)["probe"]).astype(np.complex128)
        if S.ndim > 2:
            Sm = S.reshape(S.shape[0], S.shape[1], -1)
            S = Sm[:, :, int(np.argmax((np.abs(Sm) ** 2).sum(axis=(0, 1))))]
    h5 = sorted(glob.glob(os.path.join(recon_dir, "analysis", "**", "*_recons.h5"), recursive=True), key=os.path.getmtime)
    if not h5:
        return dict(recon=os.path.basename(recon_dir), status="NO H5")
    with h5py.File(h5[-1], "r") as f:
        M = modes_of(f["reconstruction/probes"][...]).astype(np.complex128)
    if not np.isfinite(M).all():
        return dict(recon=os.path.basename(recon_dir), status="PROBE NOT FINITE")
    pw = (np.abs(M) ** 2).sum(axis=(1, 2))
    F = M[int(np.argmax(pw))]
    if F.shape != T.shape:
        return dict(recon=os.path.basename(recon_dir), status=f"SHAPE {F.shape} vs true {T.shape}")
    n = T.shape[0]
    dx = float(meta.get("dx_object_A", 0)) or None
    hk = min(n // 2, int(np.ceil((conv + MARGIN_MRAD) / da)))
    rec = dict(true_k=crop(aperture_plane(T), hk).astype(np.complex64), fitted_k=crop(aperture_plane(F), hk).astype(np.complex64),
               d_alpha_mrad=da, convergence_mrad=conv, n=n, ov_fitted=overlap(F, T), mode_power=(pw / pw.sum()).astype(np.float32))
    if dx:
        hr = min(n // 2, int(np.ceil(REAL_HALF_A / dx)))
        rec.update(dx_A=dx, true_I=crop(np.abs(T) ** 2, hr).astype(np.float32), fitted_I=crop(np.abs(F) ** 2, hr).astype(np.float32))
    if S is not None and S.shape == T.shape:
        rec.update(start_k=crop(aperture_plane(S), hk).astype(np.complex64), ov_start=overlap(S, T))
        if dx:
            rec.update(start_I=crop(np.abs(S) ** 2, hr).astype(np.float32))
    name = os.path.basename(os.path.normpath(recon_dir))
    np.savez_compressed(os.path.join(out, f"{name}.npz"), **rec)
    js = os.path.join(d01, "probe_initial.json")
    info = json.load(open(js)) if os.path.exists(js) else {}
    return dict(recon=name, status="ok", ov_start=rec.get("ov_start"), ov_fitted=rec["ov_fitted"], n_modes=int(len(M)),
                start=info.get("note") or info.get("kind"), h5=os.path.relpath(h5[-1], recon_dir))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("recon_dirs", nargs="+")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for d in a.recon_dirs:
        try:
            r = pack(d, a.out)
        except Exception as e:              # one bad dir must not stop the others
            r = dict(recon=os.path.basename(os.path.normpath(d)), status=f"ERROR {type(e).__name__}: {e}")
        rows.append(r)
        ov = "" if r.get("status") != "ok" else f"  overlap with the truth: start {r['ov_start'] if r['ov_start'] is None else round(r['ov_start'], 3)}, fitted {r['ov_fitted']:.3f}"
        print(f"{r['recon']}: {r['status']}{ov}")
    json.dump(rows, open(os.path.join(a.out, "probe_pack.json"), "w"), indent=1, default=str)
    return 0 if all(r.get("status") == "ok" for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
