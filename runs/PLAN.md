# zblur: the depth blur of a single atom, source by source

## The question

For multislice ptychography at 300 kV, how much does each source of blur widen a single atom along the beam, and
do the sources combine as a convolution? The story it would give the paper: each limiting factor (aberrations,
chromatic focal spread, source size, phonons) maps to an amount of depth blur, hence a depth resolution, and the
campaign's walls fall where that blur reaches the size of a unit cell (3.9 Å).

Hypotheses, each to be tested, not assumed:

- **H1, focal spread is a depth blur of its own size.** A chromatic focal spread of δ Å rms is, to first order, the
  probe's focus wandering along z by δ, so it should blur the atom along z by about δ.
- **H2, the blurs combine as a convolution.** In the transfer function the partial-coherence envelopes multiply; in
  real space that is a convolution of the PSFs, so widths add in quadrature (variances add) and core contrasts multiply.
- **H3, the wall is at a unit cell.** The focal-spread wall (good to 3.5 Å, off from 5 Å) sits at about one unit
  cell. If the total depth blur reaching ~3.9 Å is what kills it, the other sources' walls should sit at the same
  total blur.

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
- **The depth profile is not a Gaussian.** Through the core it is a narrow central lobe with a dip just after it
  (-0.2 to -0.3 of the peak) and side lobes 0.3 high 2-6 A away: ringing, the depth version of an Airy pattern. Each
  lateral frequency is resolved in depth over a limited band (wide for fine detail, narrow for coarse: the hourglass),
  and sharp-edged bands ring. A single Gaussian sigma is the wrong summary; report the central-lobe FWHM and the lobes.
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

## G1: each source alone (designed after G0 picks the slicing)

Base: the **ideal (aberration-free) probe**, coherent, no phonons, single atom (if G0b passes), at the slicing G0
picks; Pb first, Ti after.

| set | points |
|---|---|
| α, the perfect atom | 40, 60, 80, 100 mrad |
| aberrations | ARM tableau at 80 mrad (known) |
| chromatic focal spread | 0.75, 1.5, 2.5, 3.5, 5, 7.29 Å |
| source size | 0.05, 0.1, 0.2, 0.4 Å |
| phonons | 295, 100, 20 K |

About 18 points per element. At NL ~37-55 that is **~200-300 GPU-h for Pb**. If G0's large-blur tails reach the box,
G1 uses a deeper box (`zvac` 8 Å instead of 4).

## G2: combinations (the convolution test)

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
