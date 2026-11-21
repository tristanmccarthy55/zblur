#!/usr/bin/env python
"""
@file probe_export.py
@brief Write a finished reconstruction's fitted probe as a start probe for another (the second pass of "fit, then reconstruct").

Reads the dominant mode of reconstruction/probes in the newest <recon>/analysis/**/*_recons.h5 and writes it as probe_initial.mat in
the simulator's and sim/make_probe.py's format ({"probe": N x N complex64, "p": {...}}). The h5's probe reads back in the frame of
probe_initial_true.mat untransposed (ptychoshelves-clean/analysis/kick_probe.py, checked 2026-09-23); the grid is checked against
<recon>/01/probe_initial_true.mat. The engine rescales a start probe's intensity to the data at iteration 0, so its scale does not
matter. Prints the shift- and phase-invariant overlap with the true probe. Used by runs/zb_submit.sh PROBE_FROM_TAG (2026-11-21).

    python analysis/probe_export.py --recon <recon dir> --out <dir>/01/probe_initial.mat
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--recon", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    import h5py
    from scipy.io import loadmat, savemat
    h5 = sorted(glob.glob(os.path.join(a.recon, "analysis", "**", "*_recons.h5"), recursive=True), key=os.path.getmtime)
    if not h5:
        raise SystemExit(f"probe_export: no *_recons.h5 under {a.recon}/analysis")
    with h5py.File(h5[-1], "r") as f:
        P = np.asarray(f["reconstruction/probes"][...]).squeeze()
    if P.ndim == 3:                                            # modes on whichever axis is shortest (kick_probe.modes_of)
        P = np.moveaxis(P, -1, 0) if P.shape[-1] < P.shape[0] else P
        P = P[int(np.argmax((np.abs(P) ** 2).sum(axis=(1, 2))))]
    if not np.isfinite(P).all():
        raise SystemExit(f"probe_export: the probe in {h5[-1]} is not finite")
    T = np.asarray(loadmat(os.path.join(a.recon, "01", "probe_initial_true.mat"))["probe"]).squeeze()
    if P.shape != T.shape:
        raise SystemExit(f"probe_export: probe {P.shape} is not on the true probe's grid {T.shape}")
    x = np.fft.ifft2(np.fft.fft2(P) * np.conj(np.fft.fft2(T)))
    ov = float(np.abs(x).max() / (np.linalg.norm(P) * np.linalg.norm(T)))
    if os.path.lexists(a.out):
        os.remove(a.out)                                       # never write through a link (it could be a sim's own probe)
    savemat(a.out, {"probe": P.astype(np.complex64), "p": {"binning": False, "detector": {"binning": False}}})
    print(f"probe_export: {os.path.relpath(h5[-1], a.recon)} -> {a.out} (overlap with the true probe {ov:.3f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
