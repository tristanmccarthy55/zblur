# zblur: the depth blur of a single atom, source by source

## The question

For multislice ptychography at 300 kV, how much does each source of blur widen a single atom along the beam, and
do the sources combine as a convolution? The story it would give the paper: each limiting factor (aberrations,
chromatic focal spread, source size, phonons) maps to an amount of depth blur, hence a depth resolution, and the
campaign's walls fall where that blur reaches the size of a unit cell (3.9 Å).

Hypotheses, each to be tested, not assumed:

- **H1, focal spread is a depth blur of its own size.** A chromatic focal spread of δ Å rms is, to first order, the
  probe's focus wandering along z by δ, so it should blur the atom along z by about δ.
- **H2, the blurs combine as a convolution, whatever the base shape.** The base response (for an alpha, the aberrations
  and the recorded angle) rings rather than falling like a Gaussian (G0), so H2 is tested on transfer functions, not
  widths: T(k) = FT of the 3-D PSF; each source's envelope is E(k) = T_source(k) / T_base(k) (where |T_base| is well
  above noise); convolution means T_combination = T_base x E_1 x E_2 x ... Predicted envelopes, each with its own
  direction:

  | source | predicted envelope | acts along |
  |---|---|---|
  | chromatic focal spread delta (rms) | exp(-2 pi^2 delta^2 kz^2) | depth only |
  | source size, sigma_s = FWHM / 2.355 | exp(-2 pi^2 sigma_s^2 (kx^2 + ky^2)) | sideways only |
  | phonons, u per axis | exp(-2 pi^2 u^2 abs(k)^2), plus thermal diffuse scattering | all directions |
  | aberrations (probe known) | not an envelope: they change T_base itself | - |

  Needs every condition on one grid (same alpha, recorded angle, slice count, box) so transfer functions divide point
  by point, and a lateral crop wide enough for the low lateral frequencies (the hourglass): the lone atom (+-3.9 A),
  which works with a dense scan (G0c), rather than the 3 A grid's +-1.5 A.
- **H3, the wall is at a unit cell.** The focal-spread wall (good to 3.5 Å, off from 5 Å) sits at about one unit
  cell. If the total depth blur reaching ~3.9 Å is what kills it, the other sources' walls should sit at the same
  total blur.

## The end goal: blobs you can compute with (the user, 2026-10-27)

Treat a reconstruction as a superposition of single-atom responses ("blobs"), each a sum of simple, parameterised
blurs, and answer: near a zone axis, **how far must the crystal be tilted for one atom's depth blur to miss the next
atom's lateral blur?** The response belongs to the probe and beam, not the specimen, so for a few degrees of tilt the
blob is the same in beam coordinates and tilting only moves the atoms: neighbours along a column (spacing c ~ 3.9 A)
sit dz ~ c and dx ~ c theta apart, and their separation is the overlap of one blob with a shifted copy of itself.

So the analysis must fit each kernel to a small model, not just quote widths:
- **core**: a Gaussian, sigma_xy and sigma_z (G0d ideal: FWHM_xy ~0.13 A, FWHM_z ~1.25 A at dz 1.0);
- **cone**: the hourglass, its opening angle and its weight against |z| (G0's sections suggest ~0.1 rad);
- **lobes**: for the aberrated probe, side lobes along z (ARM: +-0.22 at ~+3 A, close to one unit cell);
and how each blur source moves those parameters (focal spread -> sigma_z, source -> sigma_xy, phonons -> amplitude
and both). Two criteria fall out: core clears core when c theta > FWHM_xy (~33 mrad, ~2 deg); the next atom leaves the
previous one's cone only when theta exceeds the cone's opening (~0.1 rad, ~5 deg). The cone's weight says which
matters.

**What G1 cannot test is the superposition itself**: at a zone axis the beam channels along the column, which is
where adding single-atom blobs could fail. G2 (below) adds that test: a two-atom column, 3.9 A apart in depth, offset
sideways by 0 / 0.1 / 0.2 / 0.3 / 0.5 A (0-7 deg), against the sum of two single-atom blobs. It needs a two-atom
geometry in ptychoshelves' simulator (not there yet).

## Phase 0, done 2026-10-21: the campaign's own kernels (`analysis/pilot_kernels.py`)

Every campaign run was scored with Pb and Ti kernels made with the same physics as the run, so single-atom PSFs at
18 slices of 1.53 Å (80 mrad, ARM aberrations known) already exist for every condition. What they show, coarsely:

- **The blur does not widen a Gaussian.** The atom's core stays one to two slices wide in every condition; each blur
  source moves weight out of the core into a broad pedestal and into the vacuum.
- **Focal spread is the depth blur:** the second-moment width of the core phase grows 1.8 → 4.4 Å from 0 to 7.29 Å
  (capped by the ±8 Å crop, so the large ones are underestimates); peak contrast falls to 0.25×.
- **Source size is a sideways blur:** the core peak falls to 0.21× at 0.4 Å but the depth width stays 1.75 Å up to
  0.2 Å.
- **Phonons cost contrast** (peak 0.34× at 100 K, 0.21× at 295 K) with little depth broadening for a lone atom.
- **Combinations: core contrasts roughly multiply** (within ~10 %), as H2 predicts.
- **Probe modes** leave a smaller pedestal but a broader, slightly shifted core.

Limits: 1.53 Å slices cannot resolve the clean core (it sits in one layer, below the slicing); the ±8 Å crop clips
the large blurs; the sims used 0.9 Å slices.

### Phase 0 addendum (2026-10-27): the wall is a core-contrast threshold?

Ranking every fixed-probe pilot kernel (Pb, 80 mrad) by the core contrast it keeps relative to the coherent beam, next
to the campaign's own verdict for that condition (PHASE1.md), one number sorts all 21 conditions -- energy spread,
source size, phonons and their combinations:

| core contrast kept | verdict | conditions |
|---|---|---|
| 1.00 - 0.59 | good, all | fs 0.75-3.5; src 0.05-0.1; fs 2.5 + src 0.1 |
| 0.52, 0.48 | okay | fs 3.5 + src 0.1; src 0.2 |
| <= 0.51 | off, all | fs 5, 7.29; src 0.4; fs 2.5 / 3.5 + src 0.2; phonons 100 / 295 K; the cooled runs; the cold FEG |

Hypothesis **H3'**: the wall is where an atom's core keeps about half its contrast, whatever the source. A correlation
from coarse kernels (grid, ARM probe, 1.53 A slices), one fuzzy pair at the boundary (fs 5: 0.51 off; 0.52 and 0.48
okay); probe-mode runs do not follow it (a different reconstruction). The two single-source walls sit where each blur
reaches the blob's own scale in its direction: sideways the source wall (0.1-0.4 A) against the core's lateral width
(~0.13 A); in depth the focal-spread wall (3.5-5 A) against the LATTICE spacing (3.9 A), not the core's depth width
(1.25 A), because focal spread moves weight into a pedestal rather than widening the core. G1 tests the rule on clean
blobs; G2's pairs at several spacings test whether the half-contrast point moves with the spacing.

## What Lei & Wang's SI gives us (arXiv 2506.22352, SI Notes 4-5, Table 3, Fig. 6)

The group's own single-atom depth test, with the same PtychoShelves LSQML engine (probe updated from the first slice):

- **One Pr atom alone in vacuum, no grid**; aberration-free probe at **defocus 0**, 300 kV; detector to 150 mrad;
  **160 x 160 scan positions**; dose 10^12 e/A^2 (noiseless in effect); **4 nm box**; 80 iterations. Slice thickness
  and scan step are not stated; the profiles are smooth over 20 A, so the slicing is fine (~0.5 A or less).
- **Depth FWHM 9.5 / 3.5 / 2.8 A at 20 / 40 / 80 mrad** (full pattern), measured as the normalised phase through the
  atom's centre against z (their Fig. 6m-o); stable up to a 0.95 alpha hollow. The x-z sections show the hourglass
  (missing cone) around a sharp core.
- In the crystal runs (Note 6) they set **defocus -0.7 nm**, a ~1 A illumination, "avoiding defocus-induced depth
  blurring"; thermal motion from refined ADPs scaled by sqrt(T / 300), 4 K taken as 23 % of room temperature.

So a lone atom DOES reconstruct cleanly with an ideal probe and a dense scan; the campaign's "a lone atom is too
sparse" came with the large ARM-aberrated probe and a coarse scan. Their 2.8 A at 80 mrad is the number our
pipeline must reproduce before any blur is measured.

## G0: fine slicing, and reproducing Lei & Wang (`runs/g0.tsv`)

- **G0a, the campaign kernel at fine slices.** One Pb grid sim (80 mrad, ARM aberrations known, coherent, 0.25 A
  simulation slices, standard 27.525 A box), reconstructed at **dz 1.53 / 0.98 / 0.74 / 0.50 A** (NL 18 / 28 / 37 / 55)
  with the S1 recipe. NL 18 against the campaign's own NL 18 kernel shows what the sim's slicing (0.25 vs 0.9 A) does.
- **G0b, Lei & Wang's test.** One Pb atom alone, ideal probe, 80 mrad, 40 A box, 0.25 A simulation slices,
  reconstructed at **dz 1.0 and 0.5 A**. One deliberate change: their crystal-run defocus (-0.7 nm, ~1 A probe) with a
  0.25 A step over a 10 A field (1681 positions, the campaign's per-iteration cost), instead of defocus 0 with a
  160 x 160 scan, which an in-focus 80 mrad probe (~0.15 A) needs for overlap.
- **Pass:** every kernel clean by `extract_psf.py`'s checks; each arm's FWHM_z and sigma_z stable between its two
  finest slicings; G0b's FWHM_z near 2.8 A. Read with `analysis/g0_report.py` (depth profiles normalised with the 2.8 A
  reference, x-z sections as their Fig. 6).
- **Cost:** two sims (G0b's has ~5x the slices of a campaign sim) + reconstructions ~5 / 8 / 11 / 16 GPU-h (G0a) and
  ~12 / 24 GPU-h (G0b): **~80-100 GPU-h**, six GPUs at once, a day or two.
- **Then:** if G0b reproduces 2.8 A, the single-atom recipe becomes G1's base (it also gives the full x-z PSF, which
  the 3 A grid crop cannot); the grid arm stays as the link to the campaign's kernels.

## G0 results (2026-10-22, pack zblur_g0_20261022_134338)

- **G0a reconstructs cleanly at fine slices** (18 / 28 / 37 layers): the x-z sections show Lei & Wang's shape, a sharp
  core inside the hourglass (missing cone), flat vacuum. 55 layers died on NaNs at beta_LSQ 0.05, batches of 32.
- **The grid's depth profile is not a Gaussian.** Through the core it is a narrow central lobe with a dip just after
  it (-0.2 to -0.3 of the peak) and side lobes 0.3 high 2-6 A away. (First read as band-limited ringing intrinsic to
  multislice ptychography; G0c's clean lone atom rules that out as the main cause, see below.)
- **The central lobe narrows with the slicing:** FWHM 2.52 / 1.75 / 1.45 A at dz 1.53 / 0.98 / 0.74, fitted by
  FWHM^2 = W^2 + (1.54 dz)^2 with **W = 0.89 A** (two parameters on three points: a hint). Far under lambda/alpha^2
  (3.1 A) and Lei & Wang's 2.8 A. A possible reason: the core is built from the highest angles recorded, resolved in
  depth roughly as lambda / (alpha theta_max) ~ 1.2 A at 200 mrad, not lambda / alpha^2. G0c tests it.
- **G0b (lone atom, -0.7 nm probe, 1600 positions) never locked on:** error flat ~280 for 90 iterations, then wrapped
  noise; 80 layers died on NaNs. The campaign's "lone atom too sparse" again; Lei & Wang used 16x the positions.
- Runs take 1-2 h (~30 s per iteration at 40 layers), ten times faster than estimated.

## G0c (`runs/g0c_slices.tsv`, `runs/g0c_det.tsv`), ~20 GPU-h

1. **Rescue the fine slices:** G0a's sim, 46 and 55 layers, `BETA_LSQ=0.02 GROUPING=64 RTAG=_b002g64` -> W from five
   slicings.
2. **Is the core width set by the recorded angle?** G0a's grid recorded to 160 and 120 mrad (28 and 37 layers): it
   should widen ~1.25x and ~1.7x if W ~ lambda / (alpha theta_max).
3. **Lei & Wang's lone atom as they ran it:** in-focus ideal probe, 0.1 A steps over 8 A (6561 positions), 40 layers.

## G0c results (2026-10-24, packs zblur_g0c_det / zblur_g0c_slices_20261023_123321)

- **Lei & Wang's lone atom works when run their way** (in-focus ideal probe, 0.1 A steps over 8 A, 6561 positions,
  40 A box, 40 layers): a clean, nearly Gaussian depth profile, **FWHM 1.48 A** at dz 1.0 (about 1.5 slices, so
  slice-limited), core share 0.91, lobes only -0.07 / +0.04; the x-z section a compact blob, faint hourglass. G0b's
  failure was scan density, not the lone atom.
- **So the grid's ringing is not intrinsic.** It appears only with the grid + ARM probe + 0.5 A scan. Prime suspect:
  the ARM-aberrated probe (its balanced C1 / C3 give it structure along z; the side lobe sits 3.5-4.5 A downstream),
  i.e. aberrations change the base shape (H2's table), as expected. G0d separates probe from grid / scan.
- **Recorded angle: inconclusive.** At 160 mrad the grid's core widens ~1.3x (2.29 vs 1.75 A at NL 28; 1.97 vs 1.45
  at NL 37), close to the predicted 1.25x; at 120 mrad it does not (1.77, 1.62 A) but its side lobes grow to 0.4-0.5
  and its core share falls.
- **Fine slices still fail:** 46 and 55 layers at beta_LSQ 0.02, batches of 64 left no kernel (logs pending).

## G0d (`runs/g0d.tsv`): the lone atom in a 20 A box, ideal and ARM probes

Every crash so far was at 46+ layers. A 20 A box gives dz 0.5 A at 40 layers: (1) whether the crash follows the
layer count or the slice thickness, (2) the lone atom's own width at dz 1.0 / 0.74 / 0.5, (3) the ARM probe on the
same lone atom (0.2 A steps over 16 A, as the probe is wide): if the ringing appears, it is the aberrated probe.
~30 GPU-h.

## G0d results (2026-10-26, pack zblur_g0d_20261026_131446)

- **The engine has a slice-thickness floor, not a layer-count limit.** In the 20 A box, 40 layers (dz 0.5 A) did not
  crash but came out as wrapped noise for the ideal probe and heavily degraded for the ARM probe; dz 1.0 and 0.74 A are
  clean. With the 27.5 A box's NaNs at 46 / 55 layers (dz 0.60 / 0.50), the floor is ~0.6-0.7 A: below it the data no
  longer pin neighbouring slices down separately (about half the atom's own depth width). So the slicing-law
  extrapolation (W ~ 0.9 A) cannot be pushed to dz -> 0 and is dropped.
- **The base: an ideal probe on a lone atom is clean and near-Gaussian.** FWHM_z 1.25 A at dz 1.0 (1.53 A at 0.74),
  lobes <= 0.13, a faint hourglass.
- **The ARM probe makes the ringing.** On the same lone atom: the same central width (1.32 A) but dips and side lobes
  of +-0.22 (ideal +-0.05) and a strong hourglass. Aberrations reshape the base (H2's table); the grid adds a little
  (0.3-0.5).
- Timing: a coherent lone-atom sim ~40 min, a 20-layer reconstruction ~4.6 h (6561 positions).

### First tilt answer, from the two G0d bases (`analysis/blobs.py`, superposition assumed)

C(theta) = how much of one atom's blob lands on the next atom's core 3.905 A down a column tilted by theta:
- **ideal probe**: ~3 % at the zone axis and <= 5 % at every tilt -- superposition alone predicts no confusion, so any
  zone-axis failure with a near-perfect probe would have to be channelling (G2);
- **ARM probe**: 30 % at the zone axis (its downstream lobe sits ~one cell away), falling in two stages: the lobe
  slides off by ~50 mrad (3 deg, to ~12-13 %), then the cone holds until ~120 mrad (7 deg) for 10 % and 150 mrad
  (8.6 deg) for 5 %.

## G1 (`runs/g1.tsv`): each source alone, combinations, and alpha, on the base

The base (ideal probe, lone Pb atom, 80 mrad, 0.1 A steps over 8 A, 20 A box, 0.25 A sim slices, dz 1.0 A) is G0d's
own run. Twenty points, each one sim and one 20-layer reconstruction:

| set | points |
|---|---|
| chromatic focal spread alone | 0.75, 1.5, 2.5, 3.5, 5 A rms |
| source size alone | 0.05, 0.1, 0.2, 0.4 A FWHM |
| phonons alone | 295, 100, 20 K |
| combinations | fs 2.5 + src 0.1; fs 2.5 + src 0.2; fs 2.5 + 20 K; src 0.1 + 20 K; fs 2.5 + src 0.1 + 20 K |
| the perfect atom against alpha | 60, 100 mrad (20 A box); 40 mrad (40 A box, 40 layers) |

Cost: 17 sims carry 16 configurations (coherence or phonons), so ~3-16x a coherent sim's 40 min each (the campaign's
own walltime rule assumes ~3x); 20 reconstructions ~4.6 h (the 40 mrad one ~9 h): **~150-290 GPU-h**, a day or two
with the queue free. The ARM base (G0d) gives the aberrations' own row; G1b later puts the key sources on it.

## G1b (`runs/g1b.tsv`): the same sources on the aberrated base

The ARM200F tableau is the probe the paper is about, so G1 is mirrored on G0d's ARM base (lone Pb atom, ARM probe at
80 mrad, 0.2 A steps over 16 A, 20 A box, dz 1.0): the same 17 sources and combinations, and the ARM-aberrated atom
against alpha (40, 60, 90 mrad, each with its own retuned tableau row; 100 mrad left out, its ~11 A probe needs a
larger window).

The test it adds: to first order the blurs ignore the aberrations. A defocus offset is free-space propagation, so a
focal spread slides the (aberrated) focus along z, a pure depth convolution; a source shift slides the probe sideways,
a pure lateral one. So each source's envelope E = T_source / T_base should come out the same on the ARM base as on the
ideal one, while the aberrations only set the base (its ringing). Where the envelopes differ, multiple scattering
through the aberrated probe's structure breaks the first-order picture. Same cost as G1, ~150-290 GPU-h.

## G2: the superposition test (two-atom columns), then the tilt answer

A pair of Pb atoms 3.9 A apart in depth, offset sideways by 0, 0.1, 0.2, 0.3, 0.5 A, ideal and ARM probes: the
reconstruction against the sum of two G0d/G1 single-atom blobs. If the sum holds off axis but fails on axis, the
tilt at which it starts to hold is the answer to the end goal, with the blob model giving it for any spacing and
blur budget. Needs a `--atoms` (explicit positions) option in sim/simulate_4dstem.py.

## (was G2) combinations (the convolution test)

fs 2.5 + src 0.1; fs 2.5 + src 0.2; fs 3.5 + src 0.1; fs 2.5 + 20 K; src 0.1 + 20 K; fs 2.5 + src 0.1 + 20 K; and
the campaign's own point, ARM + fs 2.5 + src 0.1 + 20 K. Measured against the prediction from G1's single sources.

## G3: the method

Probe modes on a few G1/G2 points (does releasing the probe narrow the PSF?), and Ti throughout.

## What is measured (`analysis/zprofile.py`)

Per kernel (`zprofile.psf3d`): the phase through the atom's centre column against depth (as Lei & Wang plot it; never
the laterally integrated phase, the zero-frequency part, which multislice ptychography does not resolve in depth),
its FWHM read off directly, a slice-averaged Gaussian width fitted within ±4 Å of the peak (so the pedestal does not
drag it), the share of the profile in that core, the lateral FWHM in the peak layer, and the x-z section. Phase 0's
core-minus-ring metrics (`pilot_kernels.py`) stay for the coarse kernels.

## The figure it builds to

A depth-blur budget: each source's added blur against its size (and against instruments), the combinations'
measured blur against the convolution prediction, and a line at one unit cell.
