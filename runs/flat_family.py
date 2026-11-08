#!/usr/bin/env python
"""
@file flat_family.py
@brief The "phase flat to X" probes of the flat experiment (2026-11-08): an 80 mrad aperture over a corrector that is flat only to X.

The user's plot: depth extent of one atom against angle, two lines. Line A: the perfect probe at aperture X. Line B: an 80 mrad
aperture whose phase is flat only to X. If B sits at the 80 mrad perfect probe's extent, aberrations beyond the flat region cost
nothing and opening the aperture is worth it.

Line B's probes are built the way a corrector is used. The HARDWARE terms (C5, A2, S3, A3, B4, D4, A5, A6, ...) are the ARM200F-class
residuals of campaign/ceos_sweep.tsv (identical in every armf row), scaled by s, a better or worse corrector. The TUNED terms are what
an operator sets: A1 and B2 aligned to zero, C1 and C3 chosen by least squares to flatten the phase over the disc of radius X (piston,
C1, C3 against chi on that disc). Flat to X means |chi - piston| <= pi/4 over that disc, the pi/2 peak-to-valley of
campaign/plan_probe.py FLAT_TOL. s is found by bisection for each X. At s = 1 the ARM's own residuals are flat to ~21 mrad.

Writes runs/flat_tableaux.tsv in campaign/ceos_sweep.tsv's columns (label flat_a080_f<XXX>), read by runs/zb_submit.sh for aber
flat<X>, and prints each probe's window need. chi is runs/window_loss.py's, checked against abTEM's Probe (tests/test_submit_dry.sh).

    python runs/flat_family.py                      # X = 30 40 60, writes runs/flat_tableaux.tsv
    python runs/flat_family.py --x 30 40 50 60      # other angles
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import window_loss as wl  # noqa: E402
LAM = wl.LAM
ALPHA = 80.0                                    # the aperture every line-B probe is opened to (mrad)
TUNED = ("C30", "C12", "phi12", "C21", "phi21")  # set by the operator, not the hardware
TOL = np.pi / 4


def hardware(sweep):
    """The ARM200F-class hardware residuals: an armf row's tableau without the tuned terms (they differ row to row, the rest do not)."""
    rows = [l.rstrip("\n").split("\t") for l in open(sweep) if l.startswith("armf_a")]
    tabs = [json.loads(r[8]) for r in rows]
    hw = {k: v for k, v in tabs[0].items() if k not in TUNED}
    for t in tabs[1:]:
        other = {k: v for k, v in t.items() if k not in TUNED}
        if other != hw:
            raise SystemExit("the armf rows differ in a hardware term; pick one by hand")
    return hw


def disc(n=241):
    t = np.linspace(-1, 1, n)
    kx, ky = np.meshgrid(t, t)
    m = np.hypot(kx, ky) <= 1
    return kx[m], ky[m]


KX, KY = disc()


def tune(ab, x_mrad):
    """(abTEM defocus, C30, half peak-to-valley) that flatten chi over the disc of radius x_mrad, least squares on piston, C1, C3."""
    k = x_mrad / 1000 / LAM
    x, y = KX * k, KY * k
    base = wl.chi(x, y, 0.0, {**ab, "C30": 0.0})
    a = LAM * np.hypot(x, y)
    A = np.c_[np.ones_like(a), 2 * np.pi / LAM * a ** 2 / 2, 2 * np.pi / LAM * a ** 4 / 4]
    coef = np.linalg.lstsq(A, -base, rcond=None)[0]           # chi = base + C10 (2 pi / lam) a^2 / 2 + C30 (2 pi / lam) a^4 / 4
    return -coef[1], coef[2], float(np.ptp(base + A @ coef)) / 2  # abTEM defocus = -C10


def flat_angle(ab):
    x = 5.0
    while x < 120 and tune(ab, x + 0.25)[2] <= TOL:
        x += 0.25
    return x


def scaled(hw, s):
    return {k: (v * s if k.startswith("C") else v) for k, v in hw.items()}


def solve(hw, x_target):
    """the largest s whose tuned probe is flat to x_target (flat angle falls as s grows)."""
    lo, hi = -6.0, 1.0                                          # log10 s
    for _ in range(40):
        mid = (lo + hi) / 2
        if flat_angle(scaled(hw, 10 ** mid)) >= x_target:
            lo = mid
        else:
            hi = mid
    s = 10 ** lo
    ab = scaled(hw, s)
    df, c3, half = tune(ab, x_target)
    full = {**ab, "C30": c3, "C12": 0.0, "phi12": 0.0, "C21": 0.0, "phi21": 0.0}
    return s, df, c3, half, full


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x", type=float, nargs="+", default=[30, 40, 60], help="flat angles (mrad)")
    ap.add_argument("--sweep", default=os.path.join(HERE, "..", "..", "ptychoshelves-clean", "campaign", "ceos_sweep.tsv"))
    ap.add_argument("--out", default=os.path.join(HERE, "flat_tableaux.tsv"))
    a = ap.parse_args()
    hw = hardware(a.sweep)
    print(f"hardware residuals as measured (s = 1): flat to {flat_angle(hw):.2f} mrad with C1, C3 tuned, A1 = B2 = 0")
    lines = ["# THE FLAT PROBES (runs/flat_family.py, 2026-11-08): an 80 mrad aperture over ARM200F-class hardware residuals scaled by s,",
             "# C1 and C3 tuned to flatten the phase over the disc of radius X, A1 = B2 = 0. Flat to X = |chi - piston| <= pi/4 there.",
             "# Columns as campaign/ceos_sweep.tsv; zb_submit.sh reads c5, c3, c1 and aber_json for aber flat<X>.",
             "label\talpha\tc5\tc3\tc1\tdf_perf\tbin\tnl\taber_json\tnote"]
    for x in a.x:
        s, df, c3, half, full = solve(hw, x)
        k80 = ALPHA / 1000 / LAM
        pv80 = float(np.ptp(wl.chi(KX * k80, KY * k80, df, full))) / (2 * np.pi)
        out = {w: wl.outside(ALPHA, df, full, w) * 100 for w in (35.004, 52.5, 70.0, 105.0)}
        label = f"flat_a{int(ALPHA):03d}_f{int(round(x)):03d}"
        note = (f"flat to {x:g} mrad: hardware x {s:.4g}, C1 {df:.2f} A, C3 {c3 / 1e4:.3f} um, {pv80:.1f} waves p-v at 80 mrad; "
                f"outside 35/52.5/70/105 A: " + "/".join(f"{v:.1f}" for v in out.values()) + " %")
        print(f"{label}: {note}")
        lines.append("\t".join([label, f"{ALPHA:g}", f"{full['C50']:.6g}", f"{c3:.6g}", f"{df:.4f}", "-", "-", "-",
                                json.dumps({k: round(v, 6) for k, v in full.items()}, separators=(",", ":")), note]))
    with open(a.out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
