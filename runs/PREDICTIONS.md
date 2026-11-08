# Predictions, fixed 2026-11-04, before any of the runs they concern

Written so the runs after this date can test them. Do not edit a threshold or a prediction after a run lands: add the outcome below
it, and if it fails, say so. Everything here was chosen after seeing the runs of `results/predict_groups.tsv`, so those runs cannot
test it.

## The ruler
`zshare`, as in `analysis/predict_groups.py`: on a kernel's core column (disc r <= 0.148 A minus ring 0.74-1.48 A, per layer), the
share of the squared phase further than 1.9525 A (half of 3.905 A) from the layer of the maximum. Kernels from the campaign's recipe:
sparse grid, 80 mrad, 210 A region, 105 A window, 18 slices of 1.529 A. Seven definitions were tried before this one was picked.

**Threshold: zshare <= 0.038 means Pb right >= 50 % (official `atom_scores` score, against the truth file).**
Standing evidence, not a test: 20 of 20 fixed-probe runs. It does not hold for mode legs (2 of 10), so it makes no claim about them.

## P1. The threshold, on new fixed-probe runs
Each new fixed-probe run with its own kernel (`v2T20_fixed`, `v2T100_fixed`, `v2T295_fixed`, vibration alone at true 20 / 100 / 295 K,
noise alone, fine slices, N = 64): Pb >= 50 % exactly when zshare <= 0.038.
Reads as holding if at most one run in every six disagrees. A run within 5 points of 50 % that disagrees counts as a disagreement.

## P2. The spacing rule, on vibration alone at true motion
E = exp(-2 pi^2 sigma^2 / (c alpha)^2), c = 3.905 A, alpha = 80 mrad; numbers soft (2.5-3.2 1/A depending on the treatment).
Wall near E = 0.1-0.2 in the lateral-only runs so far.
| run | Pb u (A) | E | prediction |
|---|---|---|---|
| true 20 K alone | 0.0385 | 0.74 | Pb >= 50 |
| true 100 K alone | 0.064 | 0.44 | Pb >= 50 |
| true 295 K alone | 0.107 | 0.10 | Pb about 50, as the old u = 0.110 run (49) |
With the Mono point's beam (fs 2.5 + src 0.1, fixed): true 20 K E = 0.51, between fs 2.5 + src 0.1 (E 0.69, Pb 91) and fs 2.5 + src 0.2
(E 0.23, Pb 45), so uncommitted; true 100 K E = 0.30, Pb < 50; true 295 K E = 0.07, Pb < 50.

## P3. Near zone axis, fixed probe (the two routes disagree on purpose)
From `results/predict_groups_tilt.tsv` (superposition assumed) against the aperture route (alpha -> alpha + theta along the tilt).
| point | on axis | 35 mrad, kernel route | 35 mrad, aperture route |
|---|---|---|---|
| fs 5 | Pb 39, zshare 0.070 | zshare 0.008: Pb >= 50 | no focal-spread term: no change |
| fs 2.5 + src 0.2 | Pb 45, zshare 0.042 | zshare 0.042: Pb < 50 | E 0.23 to 0.49: Pb >= 50 |

## P4. Slice control: Mono point, fixed, sim slices ~0.3 A against 0.9 A
- Scored against the truth file it approaches today's as-simulated score: Pb about 100, Ti about 62 (today 91 / 60); Pb rms depth about
  0.35 A (today 0.45).
- The kernel stays clean: Pb zshare near 0.023. A jump towards 0.3 would mean the sim's slices, not the window, made zblur's rings.

## P5. The ARM atom at the right window
Campaign recipe (105 A window) at recon slices of ~1 A (NL 28): no side lobes of the +-0.3-0.5 seen in zblur's 35 A window runs.
If they appear, the window is not the cause and `runs/PLAN.md`'s correction is wrong.

## Reading note, added 2026-11-04 before any run of `g1w`
P5 stands as written, but its numbers cannot decide it alone: the cut lone atom at NL 20 itself shows only -0.22 / +0.22 (`analysis/blobs.py`), under P5's +-0.3-0.5. So read P5 on `runs/g1w.tsv`'s two rows side by side (35 and 105 A windows, otherwise the same, `dip` and `side` from `analysis/blobs.py`, at NL 20; NL 28 is 0.71 A, at the slice floor). Row 1 rings and row 2 does not: the window did it. Both ring: P5 fails. Neither rings: the recording or the step did it. This does not separate the 70 A box from the 35 A window.

## P6. The flat experiment (added 2026-11-08, before any run of `gflat`)
One lone Pb atom, probe update on from the search's misses, recon slices of 1 A. Line A: the perfect probe at aperture x. Line B: an 80 mrad
aperture flat only to x (`runs/flat_family.py`). y: depth FWHM from `analysis/blobs.py`.
- Line A follows the angle down to the slice floor: about 5.6 A at 30 mrad (G0c's 3.1 A at 40, scaled as 1 / alpha^2), 3.1 A at 40, and
  the floor at 60 and 80 (ideal 60, 80 and 100 mrad all gave 1.25-1.32 A at 1 A slices).
- Aberrations cost nothing: line B within 0.3 A of line A's 80 mrad point at x = 30, 40 and 60. They cost everything: line B at line A's
  value for the same x. Only x = 30 and 40 tell the two apart; at 60 both lines sit on the floor.
- Basis: g1w (Outcomes): the ARM probe at 80 mrad, tens of waves at the edge, known and fixed, gives 1.28 A against the perfect probe's
  1.25 A. Untested: whether a released probe from a wrong start finds the same probe on a lone atom, which pins the probe far less than a
  crystal does. `runs/gflat_fixed.tsv` (true probe fixed) separates that from the aberrations.

## Outcomes
- **2026-11-08, P5 holds** (`runs/g1w.tsv`, `results/blobs_g1w.tsv`). ARM 80 mrad lone atom, NL 20 / 28. 35 A window: dip -0.51 / -0.79, side
  +0.36 / +0.43, rings. 105 A window: dip -0.11 / -0.08, side +0.06 / +0.06, FWHM 1.28 / 1.25 A, as the perfect probe (1.25 A, -0.05 / +0.05).
  The rings were the window; `runs/PLAN.md`'s correction stands.
- **2026-11-08, P1 on two new fixed-probe runs** (`ptychoshelves-clean/aberration_experiment/results/block1_scores.tsv`): true 20 K + Mono
  point, z-share 0.028, Pb 80.0 (predicted >= 50: holds); fs 7.29 at 64 configurations, z-share 0.147, Pb 30.7 (predicted < 50: holds).
- P4 pending: the slice control's analysis tarball did not arrive.
