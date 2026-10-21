# zblur -- notes for agents

Read `README.md` and `runs/PLAN.md` first. This is a side project of the PX915 campaign in `../ptychoshelves-clean`;
it reuses that repo's simulator, reconstruction and kernel extractor unchanged (never edit them from here).

- The user runs every cluster job. Give exact command blocks with their `scp -O` pull lines; never submit yourself.
- Dry-run every change to `runs/zb_submit.sh` with `tests/test_submit_dry.sh` before handing it over.
- Look at the depth profiles and check images before judging a kernel; a clean kernel passes `extract_psf.py`'s own
  checks (25 atoms, argmax at the crop centre, peak/background tens).
- Never regularize the layers (REGLAYER stays 0): it low-passes the depth axis, which is what this measures.
- Plain language for a professor-level reader; honest numbers only.
