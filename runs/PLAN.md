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

## G0: does fine slicing reconstruct cleanly? (`runs/g0.tsv`)

The campaign's kernels went to junk at 23-28 slices before the S1 fix (`PSF_KERNELS.md`: a single plane is a weak
constraint). Finer slicing is the same risk, so test it first, on the cleanest case.

- One sim: Pb grid, 80 mrad, ARM aberrations, coherent, **0.25 Å simulation slices**, the standard 27.525 Å box.
- Four reconstructions of it: **dz 1.53 / 0.98 / 0.74 / 0.50 Å** (NL 18 / 28 / 37 / 55), true probe fixed, the S1
  recipe otherwise unchanged.
- **Pass:** each kernel clean by `extract_psf.py`'s own checks (25 atoms, argmax at the crop centre and on the atom
  plane, peak/background tens, same-layer agreement high), and the core profile CONVERGES as dz shrinks (its fitted
  width stable between 0.74 and 0.50 Å). The NL 18 run against the campaign's NL 18 kernel shows what the sim's own
  slicing (0.25 vs 0.9 Å) changes.
- **Cost:** one sim (24 h wall booked; probably a few GPU-hours) + reconstructions ~5 / 8 / 11 / 16 GPU-h (scaled from
  5.3 h at NL 18): **~40-60 GPU-h**, five GPUs at most at once, about a day.

## G1: each source alone (designed after G0 picks the slicing)

Base: the **ideal (aberration-free) probe**, coherent, no phonons, at the slicing G0 picks; Pb first, Ti after.

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

Per kernel: the atom's core phase per layer (the centre minus a local ring; the laterally integrated phase is the
zero-frequency part, which multislice ptychography does not resolve in depth). From it: peak contrast, the share
within ±1.5 Å of the true depth, a model-free second-moment width, and a slice-averaged Gaussian width where the
profile spans several slices. Lateral width from the in-plane profile.

## The figure it builds to

A depth-blur budget: each source's added blur against its size (and against instruments), the combinations'
measured blur against the convolution prediction, and a line at one unit cell.
