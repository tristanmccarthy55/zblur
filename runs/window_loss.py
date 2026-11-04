#!/usr/bin/env python
"""
@file window_loss.py
@brief The share of a probe's intensity outside a reconstruction window: the number runs/zb_submit.sh refuses a row on.

Numpy only, so it runs on the login node before anything is submitted. The probe is built from the row's own values on a large grid:
a disc in k-space of semi-angle alpha, the phase of the ptychoshelves-clean tableau (abTEM's polar notation, defocus = -C10), an inverse
FFT, and the share of |psi|^2 outside a square window centred on the beam axis, which is where the solver's window sits (the same
definition as campaign/plan_probe.py window_loss). Checked against abTEM's own Probe on the Mac; tests/test_submit_dry.sh checks the
numbers below (ARM 80 mrad: 32 % at 35 A, 15 % at 52.5 A, 9 % at 70 A, 3.4 % at 105 A; 60 mrad 0.6 %, 90 mrad 43 % at 35 A).

    python runs/window_loss.py --alpha 80 --defocus -958.6 --aber '{"C30":-406600.0,...}' --window 35.004      # prints the percentage
"""
import argparse
import json

import numpy as np

LAM = 0.0196877    # A, 300 kV


def chi(k_x, k_y, defocus, ab):
    """The aberration phase (rad) in abTEM's polar notation: C_nm with its angle phi_nm, m-fold, order n + 1 in the angle."""
    a = LAM * np.hypot(k_x, k_y)
    p = np.arctan2(k_y, k_x)
    def g(n, m):
        c = ab.get(f"C{n}{m}", 0.0)
        return c * np.cos(m * (p - ab.get(f"phi{n}{m}", 0.0))) if c else 0.0
    w = (a ** 2 / 2 * (g(1, 0) - defocus + g(1, 2)) + a ** 3 / 3 * (g(2, 1) + g(2, 3)) + a ** 4 / 4 * (g(3, 0) + g(3, 2) + g(3, 4))
         + a ** 5 / 5 * (g(4, 1) + g(4, 3) + g(4, 5)) + a ** 6 / 6 * (g(5, 0) + g(5, 2) + g(5, 4) + g(5, 6)))
    return 2 * np.pi / LAM * w


def outside(alpha, defocus, ab, window, ext=None):
    """Fraction of the intensity outside the WINDOW-A square. The grid holds the probe's tails (2.5 x its 100 mrad d99 ~ 450 A at most) and
    samples the aperture 1.25 x, with a smooth, even size; a wider box only moves the result in the second decimal."""
    ext = ext or max(420.0, 3.0 * window)
    n = int(np.ceil(ext * 2 * 1.25 * (alpha / 1000.0) / LAM / 64) * 64)
    k = np.fft.fftfreq(n, ext / n).astype(np.float32)
    kx, ky = np.meshgrid(k, k, indexing="ij", sparse=True)
    disc = LAM * np.hypot(kx, ky) <= alpha / 1000.0
    psi = np.zeros((n, n), np.complex64)
    psi[disc] = np.exp(-1j * chi(kx, ky, defocus, ab))[disc].astype(np.complex64)
    psi = np.fft.ifft2(psi).astype(np.complex64)
    I = np.abs(psi) ** 2
    del psi
    h = int(round(window / 2 / (ext / n)))
    inside = I[:h, :h].sum(dtype=np.float64) + I[:h, -h:].sum(dtype=np.float64) + I[-h:, :h].sum(dtype=np.float64) + I[-h:, -h:].sum(dtype=np.float64)
    return float(1.0 - inside / I.sum(dtype=np.float64))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, required=True, help="convergence semi-angle, mrad")
    ap.add_argument("--defocus", type=float, required=True, help="abTEM defocus, A (the tsv's c1)")
    ap.add_argument("--aber", default="{}", help="aberration JSON (ceos_sweep.tsv aber_json)")
    ap.add_argument("--window", type=float, required=True, help="the solver's window side, A")
    a = ap.parse_args()
    print(f"{100 * outside(a.alpha, a.defocus, json.loads(a.aber), a.window):.1f}")
