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
| `runs/g0.tsv` | G0: fine slicing on one Pb grid |
| `tests/test_submit_dry.sh` | dry run of the submitter, checks the job wiring (run before every change) |
| `analysis/zprofile.py` | the depth profile of a kernel and its widths |
| `analysis/pilot_kernels.py` | phase 0: the campaign's existing kernels, every condition |
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

Results land under `$SHARE/phucrh/zblur_runs` (one folder per sim and per reconstruction, every physics value in its
name); the pack holds the kernels, their check images and the manifest.

## Status

- 2026-10-21: phase 0 done from the campaign's kernels (`runs/PLAN.md`, `figs/pilot_profiles.png`,
  `results/pilot_table.tsv`); G0 built and dry-run, not yet submitted.
