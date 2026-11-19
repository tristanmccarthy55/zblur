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

## P7. The experimentalist's flat experiment (added 2026-11-12, before any run of `ggrid`)
The same two lines on a grid of Pb atoms on one plane (`runs/ggrid.tsv`), the probe released from the crystal search's misses, against the
same grid with the true probe fixed. Read on W80 (the depth holding 80 % of the atom's phase on its centre column) and z-share, declared
2026-11-12 after gflat_fixed's images showed the FWHM missing a 30 mrad atom's 10 A tail (see Outcomes); fixed here before ggrid runs.
- ptychography fits the phase: every released kernel within one slice (1 A) of the fixed one in W80 and within 0.02 in z-share, and the
  80 mrad aperture flat to 30 at W80 <= 4 A while the perfect 30 mrad probe sits at >= 8 A.
- If the released kernels break as the lone atom's did, the misses are too large for one plane of atoms as well, and the claim moves to the
  crystal, or to a smaller start error (`gflat_pu0`, released from the truth, separates the two).

## P8. How close the start must be (added 2026-11-12, before any run at half or a quarter of the misses)
The lone atom, `runs/gflat.tsv`, released from 0.5 and 0.25 times the search's misses (RTAG _pu050, _pu025; A1 1.8 and 0.9 A, about 0.3
and 0.15 waves at the 80 mrad edge). Read against the fixed probe on W80 and z-share. The basin is the largest fraction at which every point
matches the fixed one within one slice in W80 and 0.02 in z-share; the experimentalist claim then reads "from a start within that, the fit
finds the phase". No outcome is predicted: 1 failed (0.6 waves), 0 held.

## P9. Which miss, and does a later release help? (added 2026-11-13, before any of these runs)
The lone atom at 80 mrad (perfect) and the 80 mrad aperture flat to 30, released three ways: from the misses with the probe held fixed through
the whole presolve and released at iteration 20 of the full engine (_pu_late); from the C1 and C3 misses with A1 right (_pu_noa1); from the A1
miss alone (_pu_a1). Why A1: for a lone atom a defocus error is close to moving the atom along the beam, which the kernel (read from its own
peak) absorbs, but twofold astigmatism splits the focus into two lines 2 A1 apart along the beam and smears the atom by about that much.
- If A1 is the culprit: _pu_noa1 matches the fixed probe (W80 within one slice, z-share within 0.02) and _pu_a1 fails as _pu did.
- If the late release helps: _pu_late comes closer to the fixed probe than _pu (W80 4 A at 80 mrad, 21 A at flat to 30 at half the misses).

## P10-P13. The three figures with the probe released (added 2026-11-16, before any of these runs)
Every run: probe released from the start, from the search's C1 and C3 misses (+1.3 A, +49 A) with A1 right (P9's `_pu_noa1`), higher
orders given; 1 A recon slices; W80 and z-share as in `analysis/flat_plot.py`. "Matches" = W80 within one slice and z-share within 0.02.
- **P10, the angle line** (`gflat_noa1` + P9's two points): all seven points match the fixed probe (fixed: perfect 9 / 8 / 4 / 3 A at
  30 / 40 / 60 / 80, flat to 30 / 40 / 60 at 80: 4 / 3 / 3 A). Basis: P9 at perfect 80 and flat to 30.
- **P11, how well A1 must be known** (`gflat_a1`, perfect 80 and flat to 30; A1 0.25 / 0.5 / 1.0 A = 0.04 / 0.08 / 0.16 waves at the
  80 mrad edge): 0.25 A matches; 1.0 A does not (z-share above 0.07, as the quarter misses' 0.9 A gave); 0.5 A is the open one.
- **P12, focal spread 2.5 / 5 / 7.29 A** (`gflat_coh`, one released mode): at 2.5 A each point within one slice of the fixed probe's
  (perfect 10 / 7 / 5 / 4 A, flat to 30: 7 A). At 7.29 A the 80 mrad point is no better than the 60 mrad one (a wall below 80 mrad),
  and every flat-to-x point sits with the perfect 80 mrad point at the same spread, not with the perfect x point. My estimate, untested.
- **P13, dose 1e6 / 1e7 / 1e8 e/A^2** (`gflat_dose`, noiseless sims with a Poisson copy each, DOSE_SEED 0; the user, 2026-11-16,
  before any dose result: 1e6-1e8 is the regime that matters with the probe released and modes on, so 1e5 and 1e4 are dropped):
  every point within one slice of its noiseless value at all three doses, so no curve-up in this range (a lone Pb atom gets 1e4 e
  per position even at 1e6). The user's expectation, the flat-to-x points rising toward the perfect-x points as dose falls, would then
  need doses under 1e6 to show. Nothing measured yet.

## P14. Does an earlier release widen the A1 basin? (added 2026-11-18, before any of these runs)
Perfect 80 and flat to 30 at 80, from C1 +1.3 A, C3 +49 A with A1 1.0 A (0.16 waves at the edge, just outside the basin) and 3.6 A (the
crystal search's miss), the probe released at iteration 5 of the presolve instead of 40 (PROBE_START=5; full engine unchanged, 20).
Basis: P9's late release (probe fixed through the presolve) was worse than the default, so the trend points this way; nothing else.
- If it helps: 1.0 A comes back within one slice and 0.02 z-share of the known probe (perfect 80 3 A / 0.012, flat to 30 4 A / 0.049).
- My expectation: at 1.0 A it moves toward the known probe; at 3.6 A it still fails (W80 at flat to 30 above 6 A), because the start,
  not the schedule, is the limit there. A guess, written to be checked.

## P15. A gentler probe, and fitting C1 and A1 only (added 2026-11-18, before any of these runs)
Same two atoms and starts as P14 (C1 +1.3 A, C3 +49 A, with A1 1.0 A or the search's 3.6 A), ptyco_pipeline branch `aberfit` (a3a762e,
run from a separate checkout so queued jobs keep main): `_pb01` the probe's step 0.1 x the object's (PROBE_BETA=0.1, release as
default); `_af` only C1 and A1 move in the presolve (PROBE_ABER_FIT=C1,A1, from iteration 5; applied as an exact phase), then the
full probe from iteration 20 of the full engine; `_afb` only C1 and A1, in both engines (the higher orders given and never touched:
the user's "fit to the fast moving C1 A1, give it the higher orders"). "Matches" as in P10.
- `_afb` matches at 1.0 A and at 3.6 A: with three free numbers the probe cannot take up the atom, and the error surface in A1 alone
  should have one minimum this close. Basis: the numpy check of the projection and the crystal search finding A1 from 40 A away on
  a fixed-probe error; nothing run.
- `_af` matches at 1.0 A; at 3.6 A it matches if the presolve's fit gets A1 inside the basin before the full release.
- `_pb01`: at 1.0 A closer to the known probe than the default; at 3.6 A it still fails (the start, not the step, is the limit). My
  weakest guess.
The smoke test (`_smoke_af`, perfect 80, 30 + 30 iterations) only shows the code runs and A1 moves; it decides nothing.
- **Revised 2026-11-19, before any P15 run.** The smoke test ran cleanly (both settings parsed, only A1a moved, A1b stayed 0, error
  fell 11.2 -> 10.1) but A1 moved 0.015 A in 30 iterations: at beta_LSQ 0.05 x PROBE_BETA 0.5 a 3.6 A miss needs ~7000 iterations, so
  `_af` and `_afb` as written would only measure slowness. Branch `aberfit` 4112cff adds PROBE_ABER_GAIN (and PROBE_ABER_MAX_STEP).
  First a gain bracket on perfect 80 from the full misses, 60 + 60 iterations, the fit in both engines (`_sg10`, `_sg30`, `_sg100`):
  my expectation, gain 30 takes |A1a| past 2.5 A without oscillating and the error below the smoke's; gain 100 may overshoot. Then
  `_af` / `_afb` at the chosen gain. `_pb01` (gentler probe) is dropped: the smoke says the probe moves too little already. In its
  place, larger steps for everything (BETA_LSQ 0.3, the default release), with the known probe at 0.3 as its own reference.

## P17. Is the dose cost the noise model? (added 2026-11-19, before any of these runs)
Dose audit (the user: "id expect 1e8/AA to be much better than this"): no fault found in the chain. add_poisson_noise.py takes the
step from positions in A (0.1), so 1e8 gives 1e6 e/pattern; Poisson-samples clean x scale; compression lossless and read by MATLAB
(the fixed 1e6 control read one); the engine rescales the probe to the data at iteration 0, fixed or released (LSQML); sim_meta gives
geometry only; links as intended. A model mismatch instead: 1e6 e/pattern on 712 px leaves the dark field at ~0.07 e/px (if ~3 % of
the beam leaves the 80 mrad disc, my guess for a lone Pb atom), where the engine's amplitude (L1) likelihood targets E[sqrt n] ~ 25 %
of sqrt(lambda) (8 % at 1e7, 3 % at 1e6): it reads the high-angle signal as nearly absent. A 256 px camera would hold ~0.5 e/px.
Test: `gflat_lk`, perfect 80 and flat to 30 at 1e8 and 1e7, true probe fixed, LIKELIHOOD=poisson (aberfit a25ddde), the same noisy
copies, against P16's L1 runs.
- If the noise model is the cause: at 1e8 perfect 80 within one slice of noiseless (3 A) and clearly better than L1's fixed 1e8;
  at 1e7 clearly better than L1's too.
- If Poisson matches L1, the mechanism is wrong (or the lone atom really is dose-limited): next a realistic camera (712 px summed
  4 x 4, UPSAMPLING=2, sim/rebin_detector.py), which raises the counts per pixel 16 x.

## P18. A real camera's pixels (added 2026-11-19, before any of these runs)
The user: "trying with a 256*256 camera to see if it improves is doable. we need a nice small step tho". zb_submit.sh CAMERA_BIN=4:
the 712 px patterns (914 at side 90) summed 4 x 4 into a 178 px (228 px) camera at the same +-200 mrad, the recon modelling the
fine grid (UPSAMPLING=2), the 0.1 A scan, window and object pixel unchanged. An exact 256 would need a crop to +-144 mrad (+-112 at
side 90), which also cuts the recorded angle, differently for the two atoms, so not this test. With DOSE the camera sums the same
noisy copies, so each pair sees the same noise. `gflat_cam`: perfect 80 and flat to 30, true probe fixed, L1, noiseless and 1e8 and
1e7, against P16 (the same runs at 712 px).
- Noiseless: within one slice and 0.02 z-share of the 712 px runs (3 A / 0.012, 4 A / 0.049): summing pixels loses little for one
  atom at +-200 mrad. Basis: detector256, where a 256 px camera with the pixels modelled kept the crystal's scores.
- At 1e8 and 1e7: if L1 misreading sub-1-count pixels is the cost (P17's mechanism), the camera beats the 712 px run at the same
  dose, most at 1e7; perfect 80 at 1e8 within one slice of noiseless. If it matches the 712 px run, counts per pixel are not it.

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
- **2026-11-08, P4 holds** (`ptychoshelves-clean/aberration_experiment/results/block1_scores.tsv`): the Mono point at 0.3 A sim slices scores
  Pb 96.0, Ti 61.1 against the truth file (0.9 A slices: 90.7, 59.7), Pb rms depth 0.37 A, z-share 0.019 (predicted ~100, ~62, ~0.35 A, near
  0.023). Within the 15-point rule: the pivot stands. The kernel stayed clean, so the sim's slices did not make zblur's rings either.
- **2026-11-12, P6 with the true probe fixed** (`gflat_fixed`; the perfect probe at 40, 60, 80 mrad from G1 / G0d, same recipe;
  `results/flat_plot.tsv`, `figs/flat_plot.png`). By P6's own metric, the FWHM, neither outcome: the perfect probe gives 3.12 A at 30 mrad,
  no wider than 3.14 A at 40, because the FWHM reads only the core; flat to 30 / 40 / 60 gives 1.93 / 1.70 / 1.26 A against 1.25 A at
  80 mrad, so "costs nothing" (within 0.3 A) holds only at 60. The images: the perfect 30 mrad atom is a weak peak (0.46 rad) with a
  ~10 A tail downstream; flat to 30 is a compact core (1.15 rad). Counting tails (W80 and z-share, chosen AFTER seeing these images, so
  not a test here): perfect probe 9 / 8 / 4 / 3 A at 30 / 40 / 60 / 80 mrad (z-share 0.51 / 0.30 / 0.031 / 0.012); flat to 30 / 40 / 60:
  4 / 3 / 3 A (0.049 / 0.025 / 0.012). With the probe known, opening past the flat region recovers most of the depth.
- **2026-11-12, P6 with the probe released from the search's misses** (`gflat`): failed at every point, the perfect 80 mrad probe's too
  (peak 0.66 against 1.85 rad, column stretched along z; `analysis/kick_probe.py`: overlap with the truth 0.40 at the start, 0.44 after).
  It says nothing about the aberrations: a lone atom did not let the probe be fitted from that start. Next: P7.
- **2026-11-12, released from the truth** (`gflat_pu0`): a lone atom holds the probe once released. Perfect 80 mrad probe: W80 3 A,
  z-share 0.012, peak 1.85 rad (fixed: 3 A, 0.012, 1.70; corrected 2026-11-13: 1.86, the G0d atom peaks a pixel off the crop centre and 1.70 is its centre pixel); flat to 30: 4 A, 0.044, 1.12 (fixed: 4 A, 0.049, 1.15). So gflat's failure was
  the start (A1 0.6 waves off at the edge), not the release or the lone atom. Next: P8.
- **2026-11-13, P8 at half the misses** (`gflat` _pu050): no point matches the fixed probe. W80 / z-share, fixed -> released: perfect 80 mrad
  3 -> 4 A / 0.012 -> 0.116; flat to 30 4 -> 21 A / 0.049 -> 0.66 (NOT CLEAN); flat to 40 3 -> 7 A; flat to 60 3 -> 4 A. With the truth as the
  start it holds, so on a lone atom the basin is under ~0.3 waves of A1 at the 80 mrad edge. A quarter pending.
- **2026-11-13, focal spread 2.5 A, true probe fixed** (`gflat_fs2p5`): perfect probe W80 10 / 7 / 5 / 4 A at 30 / 40 / 60 / 80 mrad (z-share
  0.47 / 0.28 / 0.18 / 0.13; without it 0.51 / 0.30 / 0.031 / 0.012); flat to 30 at 80 mrad: 7 A (0.26; without it 4 A, 0.049). Focal spread
  costs most at high angle but leaves 60 -> 80 still improving, so not the hard cap near sqrt(lambda / pi fs) = 50 mrad I estimated; with it,
  opening a flat-to-30 probe gains 10 -> 7 A rather than reaching the 80 mrad probe's 4 A.
- **2026-11-13, P7 VOID by design** (`runs/ggrid.tsv`): the grid recipe scans at 0.5 A, and its perfect and near-perfect probes are in focus,
  0.15-0.4 A wide, so neighbouring positions barely overlap. With the true probe fixed only flat to 30 (the one wide probe) reconstructed
  (`check_outputs.py`: 6 of 7 with no h5); released from the misses, 6 of 7 went to NaN at iteration 40 at grouping "32;16" and "8;4" alike.
  This is G0b's lesson again (`runs/PLAN.md` G0c: "G0b's failure was scan density, not the lone atom"), which my grid design ignored. The
  grid says nothing about fitting a probe. The lone atom, scanned at 0.1 A, stays the instrument.
- **2026-11-13, P8 at a quarter of the misses** (`gflat` _pu025; A1 0.9 A, ~0.15 waves at the 80 mrad edge): W80 within one slice of the
  fixed probe at every point (perfect 30 / 40 / 60 / 80: 10 / 8 / 4 / 4 A against 9 / 8 / 4 / 3; flat to 30 / 40 / 60: 5 / 4 / 4 against
  4 / 3 / 3), but z-share outside 0.02 at five of seven (perfect 80 mrad 0.079 against 0.012; flat to 30 0.076 against 0.049; flat to 40 and
  60 0.091 against 0.025 and 0.012), and peaks lower (perfect 80: 1.29 against 1.70 rad; corrected 2026-11-13 to each atom's own peak, 1.33 against 1.86: the peak column in results/flat_plot.tsv is the centre pixel, up to a pixel off the atom; W80 and z-share average a 0.15 A disc and do not move). So by P8's rule a quarter is not yet inside the
  basin; on W80 alone it is. Opening still wins: flat to 30 opened to 80, released, 5 A against 10 A at a 30 mrad aperture.
- **2026-11-16, P9** (`gflat_two`): A1 is the culprit, as predicted. Released from the C1 and C3 misses with A1 right (`_pu_noa1`): perfect
  80 mrad W80 3 A, z-share 0.013 (fixed 3 A, 0.012); flat to 30: 4 A, 0.059 (fixed 4 A, 0.049): both match. From the A1 miss alone
  (`_pu_a1`, 3.6 A): perfect 80 4 A, 0.21; flat to 30 18 A, 0.89, an in-plane cross with +-0.8 rad swings along z (check images), as bad
  as the full misses (13 A, 0.55). The late release (`_pu_late`, probe fixed through the presolve) is worse, not better: 15 A, 0.95 and
  12 A, 0.72, both NOT CLEAN; the object forms around the wrong probe first. So the probe update fits the search's C1 and C3 errors on a
  lone atom and not its A1. Next: P10-P13 from the A1-right start, and P11 for how well A1 must be known.
- **2026-11-18, P10** (`gflat_noa1` + P9's two points; probe released from C1 +1.3 A, C3 +49 A, A1 right): W80 within one slice of the
  fixed probe at all seven (perfect 30 / 40 / 60 / 80: 10 / 8 / 4 / 3 A against 9 / 8 / 4 / 3; flat to 30 / 40 / 60 at 80: 4 / 3 / 3 against
  4 / 3 / 3). z-share within 0.02 at five of seven; outside at the two weak low-angle atoms, where it is 0.3-0.5 to begin with (perfect 30
  0.385 against 0.507, less tail; perfect 40 0.328 against 0.297). So on W80 the experimentalist's line is the known-probe line: opening a
  probe flat to 30 mrad to 80 gives 4 A against 10 A at a 30 mrad aperture. Figure: `analysis/fig_flat.py --released`.
- **2026-11-18, P11** (`gflat_a1`: the same start with A1 errors at 0 rad): holds. A1 0.25 A and 0.5 A (0.04, 0.08 waves at the 80 mrad
  edge) match: perfect 80 3 A / 0.015 and 3 A / 0.018 (fixed 3 A, 0.012), flat to 30 4 A / 0.050 and 4 A / 0.056 (fixed 4 A, 0.049).
  1.0 A (0.16 waves) does not: perfect 80 4 A / 0.061, flat to 30 6 A / 0.131, the in-plane cross starting (check image). So A1 must be
  known to between 0.5 and 1.0 A at 80 mrad. A1 alone reaches pi/4 at the 80 mrad edge at 0.77 A, inside that bracket: consistent with the
  flatness rule, a derivation tested by two points either side, not a measured threshold. The crystal search missed A1 by 3.6 A.
- **2026-11-18, P13 at 1e6** (`gflat_dose`): undecided, the reconstructions FAILED at all seven. Check images: noise over the whole field,
  half the field dark, the strongest phase (+3.1 rad) in the box's last slice, not at the atom; W80 12-28 A, three NOT CLEAN, the rest
  passed extract_psf's checks while showing no atom. So these numbers measure nothing. 1e6 at 0.1 A steps is 1e4 e per pattern, 0.02 e
  per pixel on the 711 px frame; the crystal noise runs that worked had ~100x more per pixel. Two candidate causes, both untested: the
  engine's amplitude likelihood (L1, the default; Poisson exists) at those counts, and the probe released at low counts. The perfect 80
  mrad point at 1e6 with the probe fixed (one ~4 h run) separates them.
- **2026-11-19, P13 at 1e7 and 1e8, and the fixed-probe control at 1e6** (`gflat_dose` 224826 / 224841, compressed copies;
  `gflat_d1e6fix`): P13 FAILS, the released probe does not hold its noiseless values even at 1e8. W80 released, 1e8 / 1e7 (noiseless):
  perfect 30 / 40 / 60 / 80: 23 / 18 / 8 / 6 and 22 / 17 / 8 / 7 A (10 / 8 / 4 / 3); flat to 30 / 40 / 60 at 80: 7 / 7 / 6 and 12 / 8 / 7
  (4 / 3 / 3). Check images: at 1e8 perfect 80 and flat to 30 are atoms with a second lobe 3-6 A along the beam (flat to 30 also a
  diagonal streak in plane); perfect 30 is no atom, two broad lobes over the box. The fixed probe at 1e6 (perfect 80, a compressed
  copy regenerated with the same seed): one clean atom, W80 7 A, z-share 0.18, peak 0.77 rad, against the released probe's failure at
  the same dose. So these runs measure the pixel-wise probe update failing under noise as much as the dose cost; the dose figure needs
  the fixed-probe references. As a hint only (confounded): at 1e7 flat to 30 rises to 12 A while perfect 80 stays at 7, the user's
  "dose curves the aberrated line up".

## P16. The dose cost with the probe known (added 2026-11-19, before any of these runs)
`gflat_dosefix`: all seven points at 1e8 and 1e7, true probe fixed, the same noisy copies (DOSE_SEED 0). Basis: the fixed 1e6 control
(perfect 80: 7 A, one clean atom) and the noiseless line.
- At 1e8 the strong atoms (perfect 60 / 80, flat to 30 / 40 / 60; peaks above 1 rad) come within one slice of noiseless; the weak
  perfect 30 and 40 (peaks ~0.5 rad) lose more but stay single atoms.
- At 1e7, perfect 80 between its noiseless 3 A and the 1e6 control's 7 A; flat to 30 above perfect 80, because its atom is weaker
  (1.15 against 1.86 rad noiseless). My estimate.
