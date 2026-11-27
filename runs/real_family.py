#!/usr/bin/env python
"""
@file real_family.py
@brief A Spectra-class probe corrector at 300 kV, flat to 30 and to 40 mrad, and the priors the calibration search starts from (2026-11-27).

The flat family (runs/flat_family.py) scaled every residual of our ARM200F's tableau by one factor, so C5 and A5 shrank with the rest.
In a hexapole corrector they do not: C5 and the six-fold astigmatism A5 are its intrinsic fifth-order residuals; the rest is what
alignment leaves (Mueller 2006; CEOS: C5-free alignment "is considered the preferred mode of operation" for a STEM corrector, and
DCOR/ASCOR keep |A5| < 0.2 mm). The user (2026-11-27): base the instrument on a common modern one, the aberrations from the literature.

THE BASE: one measured residual set of a modern Thermo Fisher probe corrector, Yucelen, Lazic & Bosch, Sci. Rep. 8, 2676 (2018)
(PMC5805791): FEI Themis 80-300 at 300 kV, notation after Haider et al. 2000, usable convergence angle 29.4 mrad, C1 and A1
optimised by hand, no angles given:
    A2 21 nm, B2 9 nm, C3 -200 nm, A3 102 nm, S3 211 nm, A4 10.7 um, D4 8.2 um, B4 5.7 um, C5 614 um, A5 199 um, S5 5 um, R5 21 um
Haider -> abTEM (chi = 2 pi / lambda sum C_nm theta^(n+1) / (n+1) cos m(phi - phi_nm)): the round terms and the highest fold of each
order x 1 (C3, C5, A2, A3, A4, A5); the others x (n + 1): C21 = 3 B2, C32 = 4 S3, C41 = 5 B4, C43 = 5 D4, C52 = 6 S5, C54 = 6 R5.
Angles: our tableau's (campaign/ceos_sweep.tsv) where it has the term, chosen otherwise; stated as such.
As measured, C1 tuned: flat to 25.2 mrad by the pi/4 rule (|chi - piston| <= pi/4 on the disc), against their usable 29.4 mrad.
Its C5 and A5 alone cap the probe at 39.2 mrad, which matches Hartel et al. 2019 (Ultramicroscopy 206, 112821): in a correctly
aligned DCOR/ASCOR the round aberrations up to fifth order are nulled through the objective lens and the first transfer doublet,
and even the largest residual A5 (0.4 mm, large-gap pole piece) allows at least ~39 mrad at 300 kV. An S-CORR on a Themis Z keeps
the phase within pi/2 to ~50 mrad at 30 kV (arXiv 2102.12159), ~36-40 mrad scaled to 300 kV. So a well-aligned Spectra-class
instrument is flat to ~35-40 mrad (the family's 38 mrad row); flat to 30 is the same corrector less well aligned, the harder
case the calibration is run on. Not modelled: the sixth-order three-lobe D6 (abTEM's polar terms stop at fifth order; the S-CORR
corrects it, Uhlemann et al. 2022 offset it with D4).

THE FAMILY: C5 and A5 kept as measured (intrinsic) while that can reach X (to ~39 mrad); beyond, a better corrector, better aligned: every
measured term and the aligned leftovers x s, C5 and A5 included (the old flat family's rule), stated in the row's note. Within the measured corrector: A1 = 1 A and C21 = 30 A (aligned leftovers, inside the pi/4 tolerance at 40 mrad,
replacing the hand-optimised A1 and the measured B2, which the operator re-aligns); the parasitic terms A2, S3, A3, B4, D4, A4, S5, R5
scaled by one factor s (a better or worse alignment of the same corrector); C1 and C3 tuned over the flat disc (flat_family.tune).
s is bisected so the probe is flat to X. The printout sets each parasitic beside its measured value and Sawada 2017's (JEOL ARM300F).

THE PRIORS (--priors): one row per searched term: its true value (the record only) and the band the search starts in and is bounded by.
The user's rules: 4th and 5th order within literature bands; 3rd within our tableau's 95 % intervals; C1, A1, B2, A2 essentially
unknown (up to 100 x the tableau's intervals). 4th / 5th-order band = the larger of the two published measured states, with Sawada's
error: max(|Yucelen|, |Sawada| + its error); A5 capped at Hartel 2019's 0.4 mm for DCOR/ASCOR-class correctors.

    python runs/real_family.py                     # X = 30 35 38 60: writes runs/real_tableaux.tsv and runs/real_priors.tsv
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import flat_family as ff  # noqa: E402
import window_loss as wl  # noqa: E402

LAM = wl.LAM
# Yucelen 2018, abTEM length (A) and our tableau's angles where it has the term (rad)
YUCELEN = {"C23": 210.0, "phi23": 0.67777, "C21": 3 * 90.0, "phi21": 1.1, "C30": -2000.0,
           "C34": 1020.0, "phi34": 0.833831, "C32": 4 * 2110.0, "phi32": 3.038618,
           "C45": 1.07e5, "phi45": 0.939336, "C43": 5 * 8.2e4, "phi43": 0.058178, "C41": 5 * 5.7e4, "phi41": 2.460914,
           "C50": 6.14e6, "C56": 1.99e6, "phi56": 0.436332, "C52": 6 * 5.0e4, "phi52": 0.2, "C54": 6 * 2.1e5, "phi54": 0.7}
INTRINSIC = ("C50", "C56", "phi56")
PARASITIC = ("C23", "C32", "C34", "C41", "C43", "C45", "C52", "C54")
ALIGNED = {"C12": 1.0, "phi12": 0.4, "C21": 30.0, "phi21": 1.1}
# Sawada 2017 Table 3 (ARM300F, 80 kV), value and error, abTEM length (A)
SAWADA = {"C23": (279.3, 84.9), "C32": (4 * 25270.0, 4 * 3000.0), "C34": (3050.0, 1500.0), "C41": (5 * 219300.0, 5 * 186500.0),
          "C43": (5 * 135300.0, 5 * 113500.0), "C45": (98730.0, 56700.0), "C50": (-3.02e6, 7.6e6), "C56": (1.194e7, 1.5e6)}
NM = {"C10": (1, 0), "C12": (1, 2), "C21": (2, 1), "C23": (2, 3), "C30": (3, 0), "C32": (3, 2), "C34": (3, 4),
      "C41": (4, 1), "C43": (4, 3), "C45": (4, 5), "C50": (5, 0), "C52": (5, 2), "C54": (5, 4), "C56": (5, 6)}
NAME = {"C10": "C1", "C12": "A1", "C21": "B2", "C23": "A2", "C30": "C3", "C32": "S3", "C34": "A3", "C41": "B4", "C43": "D4",
        "C45": "A4", "C50": "C5", "C52": "S5", "C54": "R5", "C56": "A5"}


def waves(key, c, alpha_mrad=80.0):
    n, _ = NM[key]
    t = alpha_mrad / 1000.0
    return abs(c) * t ** (n + 1) / ((n + 1) * LAM)


def build(s):
    ab = {k: YUCELEN[k] for k in INTRINSIC}
    for p in PARASITIC:
        ab[p] = YUCELEN[p] * s
        ab["phi" + p[1:]] = YUCELEN["phi" + p[1:]]
    ab.update(ALIGNED)
    return ab


def build_all(s):
    """a better corrector, better aligned: every measured term (C5 and A5 included) and the aligned leftovers x s"""
    ab = {k: (v * s if k.startswith("C") else v) for k, v in YUCELEN.items() if k not in ("C21", "phi21", "C30")}
    ab.update({k: (v * s if k.startswith("C") else v) for k, v in ALIGNED.items()})
    return ab


def bisect(builder, x):
    lo, hi = -6.0, 1.5
    for _ in range(40):
        mid = (lo + hi) / 2
        if ff.flat_angle(builder(10 ** mid)) >= x:
            lo = mid
        else:
            hi = mid
    if ff.flat_angle(builder(10 ** lo)) < x:          # never a silent floor: the row must really be flat to X
        raise SystemExit(f"real_family: no scale makes the probe flat to {x:g} mrad")
    return 10 ** lo


def solve(x):
    """(s, scope, defocus, C3, tableau): the parasitics scaled while that can reach X; beyond what C5 and A5 allow, every term"""
    if ff.flat_angle(build(1e-9)) >= x:
        s, scope, ab = bisect(build, x), "parasitics", None
        ab = build(s)
    else:
        s, scope = bisect(build_all, x), "every term (a better corrector)"
        ab = build_all(s)
    df, c3, half = ff.tune(ab, x)
    return s, scope, df, c3, {**ab, "C30": c3}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x", type=float, nargs="+", default=[30, 35, 38, 60])
    ap.add_argument("--out", default=os.path.join(HERE, "real_tableaux.tsv"))
    ap.add_argument("--priors", default=os.path.join(HERE, "real_priors.tsv"))
    ap.add_argument("--low-width", type=float, default=100.0, help="C1, A1, B2, A2 bound as a multiple of the tableau's 95 %% intervals")
    a = ap.parse_args()
    print(f"C5 {YUCELEN['C50'] / 1e7:g} mm, A5 {YUCELEN['C56'] / 1e7:g} mm + aligned leftovers alone: flat to {ff.flat_angle(build(1e-9)):.1f} mrad;"
          f" as measured (s = 1): flat to {ff.flat_angle(build(1.0)):.1f} mrad")
    lines = ["# SPECTRA-CLASS PROBES (runs/real_family.py, 2026-11-27): the measured residuals of a Themis 80-300 at 300 kV (Yucelen et al.",
             "# 2018), C5 / A5 as measured, parasitics x s, A1 / B2 aligned leftovers, C1 / C3 tuned over the flat disc. ceos_sweep columns.",
             "label\talpha\tc5\tc3\tc1\tdf_perf\tbin\tnl\taber_json\tnote"]
    rows = {}
    for x in a.x:
        s, scope, df, c3, full = solve(x)
        k80 = 0.080 / LAM
        pv80 = float(np.ptp(wl.chi(ff.KX * k80, ff.KY * k80, df, full))) / (2 * np.pi)
        out = " / ".join(f"{100 * wl.outside(80, df, full, w):.1f}" for w in (45.0, 70.0, 105.0))
        label = f"real_a080_f{int(round(x)):03d}"
        note = (f"flat to {x:g} mrad: {scope} x {s:.4g} of Yucelen 2018, C1 {df:.2f} A, C3 {c3 / 1e4:.3f} um, {pv80:.1f} waves p-v at 80 mrad; "
                f"outside 45 / 70 / 105 A: {out} %")
        print(f"{label}: {note}")
        print("   term   here (A)       Yucelen 2018 (A)   Sawada 2017 (A)      waves at 80 mrad")
        for k in ("C12", "C21", "C23", "C30", "C32", "C34", "C41", "C43", "C45", "C50", "C52", "C54", "C56"):
            sw = SAWADA.get(k)
            yv = YUCELEN.get(k)
            print(f"   {NAME[k]:4s} {full[k]:13.5g}   {('%.4g' % yv) if yv is not None else '-':>14s}   "
                  f"{('%.4g +- %.2g' % sw) if sw else '-':>18s}   {waves(k, full[k]):8.2f}")
        rows[x] = (df, full)
        lines.append("\t".join([label, "80", f"{full['C50']:.6g}", f"{c3:.6g}", f"{df:.4f}", "-", "-", "-",
                                json.dumps({k: round(v, 6) for k, v in full.items()}, separators=(",", ":")), note]))
    with open(a.out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {a.out}")
    df, full = rows[min(a.x)]                         # the priors belong to the calibration probe (the flattest-to-least one)
    tab = {"C10": 39.5, "C12": 59.0, "C21": 3 * 188.0, "C23": 334.0, "C30": 8040.0, "C32": 4 * 3880.0, "C34": 9540.0}
    pl = ["# PRIORS for the calibration search (runs/real_family.py, 2026-11-27). true_C_A is the record only; the search draws its start",
          "# inside the bound and is confined to it. bound_A: +- per component (vectors) or +- (scalars). C10 row: C1 as the search sets",
          "# it (abTEM defocus = -C10, the tableaux' c1 column).",
          "key\tname\tn\tm\ttrue_C_A\ttrue_phi_rad\tbound_A\tbound_waves_80\tsource"]
    for k in ("C10", "C12", "C21", "C23", "C30", "C32", "C34", "C41", "C43", "C45", "C50", "C52", "C54", "C56"):
        n, m = NM[k]
        c = df if k == "C10" else full.get(k, 0.0)
        phi = full.get("phi" + k[1:], 0.0)
        if k in tab and k not in ("C30", "C32", "C34"):
            b, src = a.low_width * tab[k], f"essentially unknown: {a.low_width:g} x our tableau's 95 % interval ({tab[k]:g} A)"
        elif k in ("C30", "C32", "C34"):
            b, src = tab[k], "our tableau's 95 % interval (campaign/arm200f_tableau.tsv, run 3)"
        elif k == "C56":
            b, src = 4.0e6, "Hartel 2019: 0.03-0.4 mm for DCOR/ASCOR-class correctors"
        elif k in SAWADA:
            sv, se = SAWADA[k]
            b = max(abs(YUCELEN[k]), abs(sv) + se)
            src = "max(|Yucelen 2018|, |Sawada 2017| + its error)"
        else:
            b, src = 3 * abs(YUCELEN[k]), "3 x Yucelen 2018 (the only measured value found)"
        pl.append(f"{k}\t{NAME[k]}\t{n}\t{m}\t{c:.6g}\t{phi:.6g}\t{b:.6g}\t{waves(k, b):.3g}\t{src}")
    with open(a.priors, "w") as f:
        f.write("\n".join(pl) + "\n")
    print(f"wrote {a.priors}")
    print("\n".join(pl[3:]))


if __name__ == "__main__":
    main()
