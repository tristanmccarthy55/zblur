# zblur

How much does each source of blur (aberrations, chromatic focal spread, source size, phonons) widen a single atom
along the beam in multislice ptychography at 300 kV, and do they combine as a convolution? A side project of the
PX915 aberration campaign (`../ptychoshelves-clean`), started 2026-10-21. The plan, hypotheses and gates are in
[`runs/PLAN.md`](runs/PLAN.md).

## Layout

| path | what |
|---|---|
| `runs/PLAN.md` | the question, the hypotheses, the gates G0-G3 with costs |
| `runs/zb_submit.sh` | submits sims -> recons -> kernels -> pack on Blythe, through ptychoshelves-clean's own job scripts |
| `runs/g0.tsv` | G0: fine slicing on one Pb grid (G0a) and Lei & Wang's lone-atom test (G0b) |
| `runs/g0c_slices.tsv`, `runs/g0c_det.tsv` | G0c: the fine slices rescued, the recorded-angle test, the lone atom retried |
| `runs/g0d.tsv` | G0d: the lone atom in a 20 Å box (fine slices, few layers), ideal and ARM probes |
| `runs/g1.tsv` | G1: every blur source alone, combinations, and α, on the clean base |
| `tests/test_submit_dry.sh` | dry run of the submitter, checks the job wiring (run before every change) |
| `analysis/zprofile.py` | the depth profile of a kernel and its widths |
| `analysis/pilot_kernels.py` | phase 0: the campaign's existing kernels, every condition |
| `analysis/g0_report.py` | G0: the 3-D PSF of each kernel in a pack, convergence with slicing, the 2.8 Å reference |
| `results/`, `figs/` | tables and figures |

## Needs

- A **ptychoshelves-clean checkout** at `PTYCHO_ROOT` (on Blythe the campaign's
  `/springbrook/share/physics/phucrh/ptyco_baseline/ptyco_pipeline`). The submitter uses its `sim/run_sim.slurm`,
  `run_recon_synthetic_ML.slurm`, `analysis/atomfind/extract_psf.py` and `campaign/ceos_sweep.tsv` unchanged, refuses if
  those are modified, and records the commit in every manifest.
- Locally, `~/hyperspy-bundle/bin/python` (numpy, scipy, matplotlib) for the analysis.

## Running

```bash
# on Blythe, from this repo
bash tests/test_submit_dry.sh runs/g0.tsv "$PTYCHO_ROOT"     # dry run first
bash runs/zb_submit.sh runs/g0.tsv                           # prints the pack and its scp line
# on the Mac
scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/zblur_g0_2*.tgz' ~/Desktop/zblur/
~/hyperspy-bundle/bin/python analysis/pilot_kernels.py      # phase 0, from the campaign's kernels on ~/Desktop
```

**Do not `git pull` the ptychoshelves checkout while zblur jobs are queued**: MATLAB reads the engine when a job runs,
not when it is submitted. Each kernel records the commit it actually ran under (`psf_<name>_pin.txt`) and
`g0_report.py` flags any that differ from the manifest's.

Results land under `$SHARE/phucrh/zblur_runs` (one folder per sim and per reconstruction, every physics value in its
name); the pack holds the kernels, their check images and the manifest.

## Status

- 2026-10-21: phase 0 done from the campaign's kernels (`runs/PLAN.md`, `figs/pilot_profiles.png`,
  `results/pilot_table.tsv`). G0 built and dry-run, not yet submitted: after reading Lei & Wang's SI it has a second arm,
  their lone-atom test (target: 2.8 Å FWHM at 80 mrad).
- 2026-10-22/23: G0 ran (`runs/PLAN.md`, G0 results): the grid reconstructs cleanly to 37 layers; the depth response
  rings (central lobe + dips + side lobes), its central lobe narrowing with the slicing (W ~ 0.9 Å); the lone atom
  failed; 55 / 80 layers hit NaNs. G0c built and dry-run.
- 2026-10-26: G0c: the lone atom works with a dense in-focus scan (clean, near-Gaussian, FWHM 1.48 Å at dz 1.0), so
  the grid's ringing is not intrinsic (suspect: the ARM probe); the recorded-angle test is inconclusive; 46 / 55
  layers still fail. G0d built and dry-run.
- 2026-10-27: G0d: the engine has a slice-thickness floor (~0.6-0.7 Å; 0.5 Å gives noise whatever the layer count);
  the ideal-probe lone atom is a clean near-Gaussian base (FWHM 1.25 Å at dz 1.0); the ARM probe makes the ringing.
  G1 built and dry-run.
