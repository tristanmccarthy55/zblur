#!/usr/bin/env bash
# @file test_submit_dry.sh
# @brief Dry-run runs/zb_submit.sh on a points file against a local ptychoshelves checkout and check the wiring:
#        one sim per point, one recon per slice count waiting on its sim, one kernel job per recon waiting on it, one
#        pack waiting on every kernel, unique names, and zdrop = round(zvac / dz). Then, whatever the points file, the guards
#        of 2026-11-04: the window (calibrated against abTEM), old files' lines against HEAD's script, a side, a truncated sim,
#        a failed recon's leftovers, REDO, kernels only, PACK_ALSO, CLEAN_RAW, the disk check.
#   bash tests/test_submit_dry.sh [runs/g0.tsv] [path/to/ptychoshelves-clean]
#   BETA_LSQ=0.02 GROUPING=64 RTAG=_x bash tests/test_submit_dry.sh runs/g0c_slices.tsv   # overrides are checked too
set -euo pipefail
cd "$(dirname "$0")/.."
PTS="${1:-runs/g0.tsv}"; PR="${2:-}"
if [ -z "$PR" ]; then for p in ../ptychoshelves-clean "${SHARE:-/springbrook/share/physics}/phucrh/ptyco_baseline/ptyco_pipeline"; do [ -d "$p" ] && { PR="$p"; break; }; done; fi
[ -d "$PR" ] || { echo "no ptychoshelves checkout (give it as the second argument)"; exit 1; }
PYG=""; for p in "${ZB_GUARD_PY:-}" "$HOME/hyperspy-bundle/bin/python" "${SHARE:-/springbrook/share/physics}/phucrh/envs/abtem/bin/python" python3; do
  [ -n "$p" ] && "$p" -c 'import numpy, h5py' 2>/dev/null && { PYG="$p"; break; }; done
[ -n "$PYG" ] || { echo "no python with numpy and h5py for the guards (set ZB_GUARD_PY)"; exit 1; }
export ZB_GUARD_PY="$PYG"; export PYTHONPYCACHEPREFIX="${PYTHONPYCACHEPREFIX:-$(mktemp -d)}"
# ALLOW_CUT: this first part tests the wiring of any file, ARM rows included; the window guard has its own checks below
out="$(ALLOW_CUT=1 DRY_RUN=1 ALLOW_DIRTY=1 PTYCHO_ROOT="$PR" ZB_ROOT=/zbroot PACK_DIR=/pack SHARE=/share bash runs/zb_submit.sh "$PTS" 2>&1)"
fail=0; ok(){ echo "ok   $1"; }; bad(){ echo "FAIL $1"; fail=1; }
npts=$(grep -vc '^#' "$PTS" || true); nrec=$(grep -v '^#' "$PTS" | awk -F'\t' '{n += split($14, a, ",")} END{print n}')
nsim=$(echo "$out" | grep -c 'job-name=zb_sim_' || true); nr=$(echo "$out" | grep -c 'job-name=zb_rec_' || true)
ne=$(echo "$out" | grep -c 'job-name=zb_psf_' || true); np=$(echo "$out" | grep -c 'job-name=zb_pack_' || true)
[ "$nsim" = "$npts" ] && ok "one sim per point ($nsim)" || bad "sims $nsim, points $npts"
[ "$nr" = "$nrec" ] && ok "one recon per slice count ($nr)" || bad "recons $nr, expected $nrec"
[ "$ne" = "$nrec" ] && ok "one kernel job per recon ($ne)" || bad "kernel jobs $ne"
[ "$np" = "1" ] && ok "one pack" || bad "packs $np"
# every recon waits on the sim submitted before it; every kernel job on the recon just before it
echo "$out" | awk '/^DRY\[/{id=substr($1,5); sub(/\].*/,"",id); if ($0 ~ /zb_sim_/) sim=id; if ($0 ~ /zb_rec_/) {rec=id; if ($0 !~ "afterok:dry" sim " ") print "recon " id " not on sim " sim} if ($0 ~ /zb_psf_/ && $0 !~ "afterok:dry" rec " ") print "kernel " id " not on recon " rec}' >/tmp/zbdep.$$ || true
[ ! -s /tmp/zbdep.$$ ] && ok "dependencies wired" || { bad "dependencies"; cat /tmp/zbdep.$$; }; rm -f /tmp/zbdep.$$
d=$(echo "$out" | grep -o 'job-name=zb_rec_[^ ]*' | sort | uniq -d); [ -z "$d" ] && ok "recon names unique" || bad "duplicate $d"
echo "$out" | grep -q 'cd '"$PR"')' && ok "sbatch from the ptychoshelves checkout" || bad "submit dir"
# zdrop: grid = round(zvac / dz) (the vacuum bands); single = round(2 / dz) (2 A at each face)
echo "$out" | grep -E 'dz [0-9.]+ A, zdrop [0-9]+ \(' | sed -E 's/.*dz ([0-9.]+) A, zdrop ([0-9]+) \((grid|single), zvac ([0-9.]+)\).*/\1 \2 \3 \4/' | \
  awk '{z = ($3 == "grid") ? int($4 / $1 + 0.5) : int(2 / $1 + 0.5); if (z != $2) {print "zdrop " $2 " at dz " $1 " (" $3 "), expected " z; e = 1}} END{exit e}' \
  && ok "zdrop per geometry" || bad "zdrop"
# crops: a single atom is extracted with --single, a grid never; both crops are set in A, so half-xy x dx stays
# 1.476 A (grid) / 3.937 A (single) at any recorded angle, and --dx is always passed
cr=$(echo "$out" | grep 'job-name=zb_psf_' | sed -E 's/.*job-name=zb_psf_([^ ]*).*--half-xy ([0-9]+) --dx ([0-9.]+).*/\1 \2 \3/' | \
     awk '{g = ($1 ~ /_grid_/); w = $2 * $3; want = g ? 1.476 : 3.937; if (w < want - $3 || w > want + $3) print "crop " $1 " " w " A, want " want}')
ns=$(echo "$out" | grep -c 'job-name=zb_psf_.*_single_' || true); ng=$(echo "$out" | grep -c 'job-name=zb_psf_.*_grid_' || true)
s1=$(echo "$out" | grep 'job-name=zb_psf_.*_single_' | grep -c -- '--single' || true)
g1=$(echo "$out" | grep 'job-name=zb_psf_.*_grid_' | grep -c -- '--single' || true)
[ -z "$cr" ] && [ "$s1" = "$ns" ] && [ "$g1" = 0 ] && ok "kernel crops in A: ${ns} single, ${ng} grid" || { bad "kernel crops"; echo "$cr"; }
# a recorded angle other than 200 mrad reaches the sim and the name; 200 is never passed (those sims are reused as they are)
grep -v '^#' "$PTS" | awk -F'\t' '{print $8}' | while read -r d; do
  if [ "$d" = 200 ]; then echo "$out" | grep 'job-name=zb_sim_' | grep -v -E 'job-name=zb_sim_[^ ]*_d[0-9]+ ' | grep -q 'DETECTOR_MAX_ANGLE' && echo "BAD det200"
  else echo "$out" | grep "job-name=zb_sim_.*_d${d} " | grep -q "DETECTOR_MAX_ANGLE=${d}" || echo "BAD det${d}"; fi
done | grep -q BAD && bad "detector angle" || ok "detector angle: 200 default, others passed and named"
# engine overrides reach every reconstruction, and the tag every reconstruction and kernel name
bl="${BETA_LSQ:-0.05}"; gr="${GROUPING:-32}"; rt="${RTAG:-}"
nbad=$(echo "$out" | grep 'job-name=zb_rec_' | grep -v -c "BETA_LSQ=${bl}" || true); nbad2=$(echo "$out" | grep 'job-name=zb_rec_' | grep -v -c "GROUPING=${gr},BETA" || true)
[ "$nbad" = 0 ] && [ "$nbad2" = 0 ] && ok "beta_LSQ ${bl}, grouping ${gr} on every recon" || bad "engine overrides"
if [ -n "$rt" ]; then
  nt=$(echo "$out" | grep -E 'job-name=zb_(rec|psf)_' | grep -v -c -- "NL[0-9]*${rt} " || true)
  [ "$nt" = 0 ] && ok "tag ${rt} on every recon and kernel" || bad "tag missing on ${nt}"
fi
# probe modes (PROBE_MODES > 1): every recon released with the campaign's mode-leg switches and tagged _m<modes> (unless
# RTAG is given); a point with focal spread or source size starts from its mixed state, a coherent one does not; 48 h
pm="${PROBE_MODES:-1}"
if [ "$pm" -gt 1 ]; then
  rec=$(echo "$out" | grep 'job-name=zb_rec_')
  nm=$(echo "$rec" | grep -v -c "PROBE_MODES=${pm},.*PROBE_START=${PROBE_START:-40},PROBE_START2=${PROBE_START2:-20},PROBE_SUPPORT_FFT=1" || true)
  [ "$nm" = 0 ] && ok "${pm} modes, probe released (${PROBE_START:-40} / ${PROBE_START2:-20}), aperture on, every recon" || bad "mode switches missing on ${nm}"
  nmix=$(grep -v '^#' "$PTS" | awk -F'\t' '$11 != "0" || $12 != "0" {n += split($14, a, ",")} END{print n + 0}')
  gmix=$(echo "$rec" | grep -c "PROBE_MIXED=${pm}" || true)
  cmix=$(echo "$rec" | grep 'PROBE_MIXED' | grep -v -c -E '_(fs|src)[0-9p]+' || true)
  [ "$gmix" = "$nmix" ] && [ "$cmix" = 0 ] && ok "mixed-state start on the ${nmix} blurred recons, none on a coherent one" \
    || bad "PROBE_MIXED on ${gmix} (want ${nmix}), on ${cmix} coherent"
  tg="${RTAG:-_m${pm}}"; nt=$(echo "$out" | grep -E 'job-name=zb_(rec|psf)_' | grep -v -c -- "NL[0-9]*${tg} " || true)
  [ "$nt" = 0 ] && ok "tag ${tg} on every recon and kernel" || bad "mode tag missing on ${nt}"
  nw=$(echo "$rec" | grep -v -c -- '--time=48:00:00' || true); [ "$nw" = 0 ] && ok "48 h walltime on every mode recon" || bad "walltime on ${nw}"
else
  echo "$out" | grep 'job-name=zb_rec_' | grep -q -E 'PROBE_MIXED|PROBE_START' && bad "a 1-mode recon got mode switches" \
    || ok "1 mode: probe held fixed (no PROBE_START / PROBE_MIXED)"
fi
# phonon sims (the sigma fix, ptychoshelves a8b9264): named _T<T>v2 so a pre-fix _T<T> folder is never reused, and
# refused, before anything is submitted, when the checkout lacks the fix (here: a commit it does not have)
nT=$(grep -v '^#' "$PTS" | awk -F'\t' '$13 != "-"' | wc -l | tr -d ' ')
if [ "$nT" -gt 0 ]; then
  nv=$(echo "$out" | grep 'job-name=zb_sim_' | grep -c -E '_T[0-9]+v2[_ ]' || true)
  nold=$(echo "$out" | grep -E 'job-name=zb_(sim|rec|psf)_' | grep -c -E '_T[0-9]+(_|NL| )' || true)
  [ "$nv" = "$nT" ] && [ "$nold" = 0 ] && ok "phonon sims named _T<T>v2 (${nv}), no pre-fix name" || bad "phonon names: v2 ${nv} of ${nT}, old ${nold}"
fi
tmp="$(mktemp "${TMPDIR:-/tmp}/zbph.XXXXXX")"; { grep -v '^#' "$PTS" | head -1 | awk -F'\t' 'BEGIN{OFS="\t"} {$13 = "20"; print}'; } >"$tmp"
po="$(PHONON_FIX=0123456789abcdef0123456789abcdef01234567 DRY_RUN=1 ALLOW_DIRTY=1 PTYCHO_ROOT="$PR" ZB_ROOT=/zbroot PACK_DIR=/pack SHARE=/share bash runs/zb_submit.sh "$tmp" 2>&1)" && pe=0 || pe=$?
rm -f "$tmp"
[ "$pe" != 0 ] && ! echo "$po" | grep -q '^DRY\[' && echo "$po" | grep -q 'need the frozen-phonon sigma fix' \
  && ok "phonon point refused without the sigma fix, nothing submitted" || bad "phonon guard (exit ${pe})"
echo "$out" | grep 'job-name=zb_sim_.*_single_' | grep -q 'GRID_SPACING=0' && ok "single atom has GRID_SPACING=0" || { [ "$ns" = 0 ] || bad "single GRID_SPACING"; }
echo "$out" | grep -q 'job-name=zb_sim_.*_ideal_' && { echo "$out" | grep 'job-name=zb_sim_.*_ideal_' | grep -q 'CS=0,C5=0,DEFOCUS=' \
  && ok "ideal probe: CS 0, C5 0, its own defocus" || bad "ideal probe exports"; }
echo "$out" | grep -q 'job-name=zb_sim_' && { echo "$out" | grep -q 'SLICE_THICKNESS=0.25' && ok "sim slice 0.25 A reaches the sim" || bad "sim slice"; }

# ============================================================ the guards of 2026-11-04
PRA="$(cd "$PR" && pwd)"; T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
dry(){ ALLOW_DIRTY=1 DRY_RUN=1 PTYCHO_ROOT="$PR" ZB_ROOT=/zbroot PACK_DIR=/pack SHARE=/share "$@" 2>&1; }
# -- the window, against abTEM's own Probe on the same tableau (square window as plan_probe.window_loss; 2026-11-04, abTEM 1.0.5, grids 140-420 A).
# PIPELINE_BIBLE and runs/PLAN.md said 35.7 / 40.8 % at 35 A until 2026-11-04: neither that square nor a disc gives it; 52.5, 70 and 105 A agree with the planner's records.
wl(){ local row; row=$(awk -F'\t' -v L="armf_a0$1" '$1==L' "$PR/campaign/ceos_sweep.tsv"); "$PYG" runs/window_loss.py --alpha "$1" --defocus="$(echo "$row" | cut -f5)" --aber "$(echo "$row" | cut -f9)" --window "$2"; }
e=""; for c in "80 35.004 32.0" "80 52.5 15.0" "80 70 9.0" "80 105 3.35" "60 35.004 0.6" "90 35.004 43.0" "90 140 5.1"; do
  set -- $c; g=$(wl "$1" "$2"); awk -v a="$g" -v b="$3" 'BEGIN{d = a - b; if (d < 0) d = -d; exit !(d <= 0.5)}' || e="${e} ${1}mrad@${2}A=${g} (want ${3})"; done
[ -z "$e" ] && ok "window_loss within 0.5 points of abTEM at seven settings" || bad "window_loss:${e}"
# -- the guard refuses the old ARM rows, before anything; the ideal rows pass
o="$(dry bash runs/zb_submit.sh runs/g1v2.tsv)" && rc=0 || rc=$?
n=$(echo "$o" | grep -c '_arm_single.*outside its' || true)
[ "$rc" != 0 ] && ! echo "$o" | grep -q '^DRY\[' && echo "$o" | grep -q '^NOT SUBMITTED' && [ "$n" = 6 ] \
  && ok "g1v2.tsv as written refused before anything is submitted; its 6 ARM rows named" || bad "g1v2 guard (exit ${rc}, ${n} rows named)"
for f in runs/g0.tsv runs/g0d.tsv runs/g1b.tsv runs/g1m.tsv; do dry bash runs/zb_submit.sh "$f" >/dev/null && bad "${f}: ARM rows not refused"; done; ok "the ARM rows of g0, g0d, g1b, g1m refused too"
grep -E '^#|ideal' runs/g1v2.tsv >"$T/ideal.tsv"
o="$(dry bash runs/zb_submit.sh "$T/ideal.tsv")" && rc=0 || rc=$?
[ "$rc" = 0 ] && [ "$(echo "$o" | grep -c 'job-name=zb_sim_')" = 6 ] && ok "g1v2's six ideal rows pass the guard" || bad "ideal rows refused: $(echo "$o" | grep -m1 'NOT SUB')"
# -- old files: every line as the script before the guards (88e3ad4) writes it, but for the check in the kernel job and the count in the pack
git show 88e3ad4:runs/zb_submit.sh >"$T/zb_ref.sh"     # the last script before the guards; HEAD would compare the guards with themselves
norm(){ grep '^DRY\[' | sed -E "s#--wrap=[^|]*check_outputs\.py' recon '[^']*' \|\| exit 1; #--wrap=#; s#; rc=[$][?]; miss=0; for f in .*##; s/[0-9]{8}_[0-9]{6}/TS/g"; }
e=""; for f in runs/g0.tsv runs/g0c_det.tsv runs/g0c_slices.tsv runs/g0d.tsv runs/g1.tsv runs/g1b.tsv runs/g1m.tsv runs/g1v2.tsv; do
  a="$(DRY_RUN=1 ALLOW_DIRTY=1 PTYCHO_ROOT="$PR" ZB_ROOT=/zbroot PACK_DIR=/pack SHARE=/share bash "$T/zb_ref.sh" "$f" 2>&1 | norm)"
  b="$(dry env ALLOW_CUT=1 bash runs/zb_submit.sh "$f" | norm)"
  [ -n "$a" ] && [ "$a" = "$b" ] || { e="${e} ${f}"; diff <(echo "$a") <(echo "$b") | head -4 | cut -c1-200 || true; }
done
[ -z "$e" ] && ok "all 8 old files: names, exports, times, dependencies as 88e3ad4's script" || bad "lines differ from 88e3ad4's:${e}"
# -- g1v2_redo.tsv is g1v2's first five rows, verbatim
[ "$(grep -v '^#' runs/g1v2_redo.tsv)" = "$(grep -v '^#' runs/g1v2.tsv | head -5)" ] && ok "g1v2_redo.tsv = g1v2 rows 1-5, verbatim" || bad "g1v2_redo.tsv rows"
# -- a side: named, exported, sized from the row; a row that cannot fit refused
o="$(dry env ALLOW_CUT=1 GROUPING='32;16' bash runs/zb_submit.sh runs/g1w.tsv)" && rc=0 || rc=$?
s1=$(echo "$o" | grep 'job-name=zb_sim_.*_w16_s0p25_bz20p000_d133 '); s2=$(echo "$o" | grep 'job-name=zb_sim_.*_rs210_d133 ')
r2=$(echo "$o" | grep 'job-name=zb_rec_.*_rs210_d133_NL20 ')
[ "$rc" = 0 ] && [ "$(echo "$o" | grep -c 'job-name=zb_rec_')" = 4 ] && [[ "$s1" != *REGION_SIDE* ]] && [[ "$s2" == *"REGION_SIDE=210"* && "$s2" == *"DETECTOR_MAX_ANGLE=133"* && "$s2" == *"SLICE_THICKNESS=0.25"* ]] \
  && [[ "$r2" == *"--mem=175G"* && "$r2" == *"GROUPING=32;16"* ]] && [ "$(echo "$s1" | sed -E 's/.*--export=//; s/JOB_DIR=[^,]*,//')" = "$(echo "$s2" | sed -E 's/.*--export=//; s/JOB_DIR=[^,]*,//; s/,REGION_SIDE=210//')" ] \
  && ok "g1w: _rs210 named and exported; the pair's sims differ by REGION_SIDE alone; 175G for 1420 px ($(echo "$r2" | grep -o -- '--time=[^ ]*'))" || { bad "g1w side"; echo "$o" | tail -3 | cut -c1-200; }
o="$(dry env GROUPING='32;16' bash runs/zb_submit.sh runs/g1w.tsv)" && rc=0 || rc=$?
[ "$rc" != 0 ] && echo "$o" | grep -q 'bz20p000_d133: .* outside its 35.004' && ! echo "$o" | grep -q 'rs210_d133: .* outside its' && ok "g1w's cut row refused without ALLOW_CUT; its 105 A row not named" || bad "g1w guard"
o="$(dry env -u GROUPING ALLOW_CUT=1 bash runs/zb_submit.sh runs/g1w.tsv)" && rc=0 || rc=$?
[ "$rc" != 0 ] && echo "$o" | grep -q 'full-engine grouping' && ok "1420 px at grouping 32 refused" || bad "grouping guard"
printf 'Pb\t80\tarm\t-\tsingle\t0.2\t16\t200\t0.25\t0.2375\t0\t0\t-\t20,28\t210\n' >"$T/literal.tsv"
o="$(dry env GROUPING='32;16' bash runs/zb_submit.sh "$T/literal.tsv")" && rc=0 || rc=$?
[ "$rc" != 0 ] && echo "$o" | grep -q 'NOT SUBMITTED.*GB of patterns' && ok "105 A window at +-200 mrad and 6561 positions refused: $(echo "$o" | grep -o '[0-9]* px x [0-9]* positions = [0-9.]* GB')" || bad "size guard"
# -- real paths, stub sbatch: what a rerun finds
sbatch(){ local n; n=$(( $(cat "$STUB_N" 2>/dev/null || echo 1000) + 1 )); echo "$n" >"$STUB_N"; printf '%s\n' "$*" >>"$STUB_LOG"; echo "$n"; }
squeue(){ printf '%s\n' "${STUB_SQUEUE:-}" | grep . || true; }
mmlsquota(){ printf '%s\n' "${STUB_QUOTA-springbrook share-physics GRP 3.356T 3.91T 3.92T 0 none | 100 0 0 0 none}"; }
export -f sbatch squeue mmlsquota
mkdir -p "$T/env/bin"; ln -s "$(command -v "$PYG")" "$T/env/bin/python"
fake(){ "$PYG" - "$1" "$2" <<'PYEOF'
import os, sys
import h5py, numpy as np
what, d = sys.argv[1:3]
if what.startswith("sim"):
    s = os.path.join(d, "01"); os.makedirs(s, exist_ok=True)
    with h5py.File(os.path.join(s, "data_dp.hdf5"), "w") as f: f.create_dataset("dp", data=np.ones((4, 8, 8), np.float32))
    with h5py.File(os.path.join(s, "data_position.hdf5"), "w") as f: f.create_dataset("probe_positions_0", data=np.zeros((2, 4), np.float32))
    for n in ("probe_initial_true.mat", "sim_meta.mat", "aberrations.json"): open(os.path.join(s, n), "w").write("x")
    if what == "sim_short":
        open(os.path.join(s, "data_dp.hdf5"), "r+b").truncate(1200); os.remove(os.path.join(s, "aberrations.json"))
    sys.exit()
m = os.path.join(d, "01", "MLs_L1_p1_g64_betaLSQ0.05_NpbstInf_synthetic_ML_Ndp356_step01"); os.makedirs(m, exist_ok=True)
if what == "stub":
    open(os.path.join(m, "Niter200.mat"), "wb").write(b"\0" * 3000)
    open(os.path.join(d, "slurm_1.out"), "w").write("Iteration: 200 / 200\nMultislice synthetic reconstruction finished.\n")
else:
    a = os.path.join(d, "analysis", "S00000-00999", "S00001"); os.makedirs(a, exist_ok=True)
    with h5py.File(os.path.join(a, "S00001_712x712_b0_ML_recon_recons.h5"), "w") as f:
        f.create_dataset("reconstruction/object", data=np.ones((20, 1, 32, 32), np.complex64)); f.create_dataset("reconstruction/probes", data=np.ones((32, 32), np.complex64))
    open(os.path.join(d, "slurm_1.out"), "w").write("Reconstructed scan S00001\n")
PYEOF
}
real(){ # $1 root, then env and the script: real paths, stub sbatch -> the sbatch log in $1/sb.log, the screen in $1/out.txt
  local z="$1"; shift; mkdir -p "$z"
  ( unset DRY_RUN; env STUB_LOG="$z/sb.log" STUB_N="$z/n" ALLOW_DIRTY=1 PTYCHO_ROOT="$PRA" ZB_ROOT="$z" PACK_DIR="$z" CONDA_ENV="$T/env" DISK_PROBE_MB=1 "$@" ) >"$z/out.txt" 2>&1; }
N6=Pb_a080_ideal_df0_single_st0p1_w8_s0p25_bz20p000_T100v2; awk -F'\t' '/^#/ || ($3 == "ideal" && $13 == "100" && $11 == "0" && $12 == "0")' runs/g1v2.tsv >"$T/one.tsv"
Z="$T/z1"; fake sim_short "$Z/sim_${N6}"; real "$Z" bash runs/zb_submit.sh "$T/one.tsv" && rc=0 || rc=$?
[ "$rc" != 0 ] && [ ! -s "$Z/sb.log" ] && grep -q "not a whole sim.*rm -rf" "$Z/out.txt" && [ -e "$Z/sim_${N6}/01/data_dp.hdf5" ] && ok "a sim truncated by the quota refused, with the line that removes it" || { bad "truncated sim"; tail -2 "$Z/out.txt" | cut -c1-200; }
Z="$T/z2"; fake sim "$Z/sim_${N6}"; fake stub "$Z/recon_${N6}_NL20"; mkdir -p "$Z/psf"; : >"$Z/psf/psf_${N6}_NL20_vol.npy"; real "$Z" bash runs/zb_submit.sh "$T/one.tsv" && rc=0 || rc=$?
[ "$rc" != 0 ] && [ ! -s "$Z/sb.log" ] && grep -q "earlier attempt.*REDO=1" "$Z/out.txt" && ok "a failed recon's leftovers (a 3 KB Niter200.mat) refused without REDO" || { bad "leftovers refused"; tail -2 "$Z/out.txt" | cut -c1-200; }
real "$Z" REDO=1 bash runs/zb_submit.sh "$T/one.tsv" && rc=0 || rc=$?
[ "$rc" = 0 ] && ! grep -q 'job-name=zb_sim_' "$Z/sb.log" && grep -q "job-name=zb_rec_${N6}_NL20 " "$Z/sb.log" && [ -z "$(find "$Z/recon_${N6}_NL20" -name 'Niter*.mat')" ] && [ ! -e "$Z/psf/psf_${N6}_NL20_vol.npy" ] && [ -f "$Z/recon_${N6}_NL20/slurm_1.out" ] \
  && [ "$(readlink "$Z/recon_${N6}_NL20/01/data_dp.hdf5")" = "$Z/sim_${N6}/01/data_dp.hdf5" ] && ok "REDO=1: stub and stale kernel cleared, logs kept, same name, sim reused, no new sim" || { bad "REDO"; tail -3 "$Z/out.txt" | cut -c1-200; }
grep "job-name=zb_psf_" "$Z/sb.log" | grep -q "check_outputs.py' recon '$Z/recon_${N6}_NL20' || exit 1; git -C" && ok "the kernel job checks the recon first" || bad "kernel job check"
# the pack's own wrap, run: a tar shim that succeeds, then fails; the kernels there, then one gone
mkdir -p "$T/shim"; printf '#!/bin/sh\nexit ${SHIM_TAR_RC:-0}\n' >"$T/shim/tar"; chmod +x "$T/shim/tar"
w="$(grep "job-name=zb_pack_" "$Z/sb.log" | sed 's/.*--wrap=//')"; vf="$Z/psf/psf_${N6}_NL20_vol.npy"
echo x >"$vf"; PATH="$T/shim:$PATH" sh -c "$w" >/dev/null 2>&1 && a=$? || a=$?
PATH="$T/shim:$PATH" SHIM_TAR_RC=2 sh -c "$w" >/dev/null 2>&1 && b=$? || b=$?
rm -f "$vf"; PATH="$T/shim:$PATH" sh -c "$w" >/dev/null 2>&1 && c=$? || c=$?
[ "$a" = 0 ] && [ "$b" != 0 ] && [ "$c" != 0 ] && ok "the pack, run: exit 0 when whole, non-zero on a failed tar or a missing kernel" || bad "pack wrap (whole ${a}, tar ${b}, missing ${c})"
Z="$T/z3"; fake sim "$Z/sim_${N6}"; fake leg "$Z/recon_${N6}_NL20"; real "$Z" bash runs/zb_submit.sh "$T/one.tsv" && rc=0 || rc=$?
[ "$rc" != 0 ] && grep -q "whole reconstruction" "$Z/out.txt" && ok "a whole recon refused without REDO" || bad "whole recon refused"
real "$Z" REDO=1 bash runs/zb_submit.sh "$T/one.tsv" && rc=0 || rc=$?
[ "$rc" = 0 ] && ! grep -q 'job-name=zb_rec_' "$Z/sb.log" && grep "job-name=zb_psf_" "$Z/sb.log" | grep -v -q -- '--dependency' && [ -d "$Z/recon_${N6}_NL20/analysis" ] && ok "REDO=1 on a whole recon: kernel only, nothing cleared" || bad "kernel only"
Z="$T/z4"; fake sim "$Z/sim_${N6}"; real "$Z" STUB_SQUEUE="zb_rec_${N6}_NL20" bash runs/zb_submit.sh "$T/one.tsv" && rc=0 || rc=$?
[ "$rc" != 0 ] && grep -q "still in the queue" "$Z/out.txt" && ok "a recon still in the queue refused" || bad "queued recon"
Z="$T/z5"; fake sim "$Z/sim_${N6}"; mkdir -p "$Z/manifests"; K=Pb_a080_ideal_df0_single_st0p1_w8_s0p25_bz20p000_fs2p5_src0p1_T20v2_NL20
awk -v k="$K" 'BEGIN{OFS="\t"; for (i = 1; i <= 23; i++) $i = "x"; $24 = k; $1 = "x"; print}' >"$Z/manifests/g1v2_old.tsv"
real "$Z" PACK_ALSO="$K" PACK_AFTER=9999 bash runs/zb_submit.sh "$T/one.tsv" && rc=0 || rc=$?
grep "job-name=zb_pack_" "$Z/sb.log" | grep -q "afterany:[0-9:]*:9999" && grep "job-name=zb_pack_" "$Z/sb.log" | grep -q "psf/psf_${K}_vol.npy" && grep -q "$K" "$Z"/manifests/one_*.tsv && ok "PACK_ALSO: waits for it, packs its kernel, copies its manifest row" || bad "PACK_ALSO"
# -- CLEAN_RAW: only where the recon is whole and the kernel exists
Z="$T/z6"; fake sim "$Z/sim_${N6}"; real "$Z" CLEAN_RAW=1 bash runs/zb_submit.sh "$T/one.tsv"; [ -e "$Z/sim_${N6}/01/data_dp.hdf5" ] && [ ! -s "$Z/sb.log" ] && a=1 || a=0
fake leg "$Z/recon_${N6}_NL20"; mkdir -p "$Z/psf"; echo x >"$Z/psf/psf_${N6}_NL20_vol.npy"; real "$Z" CLEAN_RAW=1 bash runs/zb_submit.sh "$T/one.tsv"
[ "$a" = 1 ] && [ ! -e "$Z/sim_${N6}/01/data_dp.hdf5" ] && [ -e "$Z/sim_${N6}/01/sim_meta.mat" ] && [ ! -s "$Z/sb.log" ] && ok "CLEAN_RAW: keeps the sim until its recon is whole and its kernel exists; submits nothing" || bad "CLEAN_RAW"
# -- the disk: the group's row, the margin, a real write
o="$(dry env DISK_CHECK=1 DISK_PROBE_MB=1 DISK_PROBE_DIR="$T" bash runs/zb_submit.sh "$T/ideal.tsv")" && rc=0 || rc=$?
[ "$rc" = 0 ] && echo "$o" | grep -q 'test write of 1 MB ok' && ok "disk: quota row read, test write done: $(echo "$o" | grep -m1 'group physics' | cut -c3-)" || { bad "disk ok case"; echo "$o" | tail -2 | cut -c1-200; }
o="$(dry env DISK_CHECK=1 COMMITTED_GB=560 DISK_PROBE_MB=1 DISK_PROBE_DIR="$T" bash runs/zb_submit.sh "$T/ideal.tsv")" && rc=0 || rc=$?
[ "$rc" != 0 ] && ! echo "$o" | grep -q '^DRY\[' && echo "$o" | grep -q 'NOT SUBMITTED.*margin' && ok "disk: refused when the margin would not be left" || bad "disk margin"
o="$(dry env DISK_CHECK=1 STUB_QUOTA= DISK_PROBE_MB=1 bash runs/zb_submit.sh "$T/ideal.tsv")" && rc=0 || rc=$?
[ "$rc" != 0 ] && echo "$o" | grep -q 'NOT SUBMITTED: could not read' && ok "disk: refused when the quota row cannot be read" || bad "disk row"
o="$(dry env DISK_CHECK=1 DISK_PROBE_MB=1 DISK_PROBE_DIR=/nonexistent/dir bash runs/zb_submit.sh "$T/ideal.tsv")" && rc=0 || rc=$?
[ "$rc" != 0 ] && echo "$o" | grep -q 'NOT SUBMITTED: the test write' && ok "disk: refused when the test write fails" || bad "disk write"
# -- the flat experiment (2026-11-08): flat<X> rows take runs/flat_tableaux.tsv's own C1, C3, C5; PROBE_UPDATE=1 releases one mode
#    from the search's misses; a point still being simulated is never simulated twice
o="$(dry env PROBE_UPDATE=1 GROUPING='32;16' bash runs/zb_submit.sh runs/gflat.tsv)" && rc=0 || rc=$?
f30=$(awk -F'\t' '$1=="flat_a080_f030"' runs/flat_tableaux.tsv)
s30=$(echo "$o" | grep 'job-name=zb_sim_Pb_a080_flat30_single_st0p1_w8_s0p25_bz40p000_rs90 ')
r30=$(echo "$o" | grep 'job-name=zb_rec_Pb_a080_flat30_single_st0p1_w8_s0p25_bz40p000_rs90_NL40_pu ')
pc1=$(awk -v a="$(echo "$f30" | cut -f5)" 'BEGIN{printf "%.4f", a + 1.3}')
[ "$rc" = 0 ] && [[ "$s30" == *"DEFOCUS=$(echo "$f30" | cut -f5),"* && "$s30" == *"CS=$(echo "$f30" | cut -f4),"* && "$s30" == *"C5=$(echo "$f30" | cut -f3),"* && "$s30" == *"REGION_SIDE=90"* ]] \
  && ok "flat30: the tableau's own C1, C3, C5 reach the sim, in a 45 A window" || { bad "flat30 sim"; echo "$o" | grep -m3 -E 'NOT SUB|flat30' | cut -c1-200; }
[[ "$r30" == *"PROBE_START=40,PROBE_START2=20,PROBE_SUPPORT_FFT=1,PROBE_C1=${pc1},"* && "$r30" == *"PROBE_SET=C12=3.6000:phi12=0.000000"* && "$r30" == *"PROBE_MODES=1,"* ]] \
  && ok "PROBE_UPDATE=1: one mode released from the search's misses (C1 ${pc1} A, A1 3.6 A)" || { bad "probe update exports"; echo "$r30" | cut -c1-300; }
[ "$(echo "$o" | grep -c 'job-name=zb_rec_.*_pu ')" = 7 ] && [ "$(echo "$o" | grep -c 'job-name=zb_psf_.*_pu ')" = 7 ] && ok "gflat: seven recons and kernels, all _pu" || bad "gflat counts"
o="$(dry env PROBE_UPDATE=1 START_ERR='C1=0:C3=0:A1=0' START_SCALE='C50=0:C23=0' RTAG=_pu_noho GROUPING='32;16' bash runs/zb_submit.sh runs/gflat.tsv)" && rc=0 || rc=$?
echo "$o" | grep 'job-name=zb_rec_Pb_a080_flat40_.*_pu_noho ' | grep -q 'PROBE_SCALE=C50=0:C23=0' && ok "START_SCALE reaches the start probe" || bad "START_SCALE"
o="$(dry env PROBE_UPDATE=1 PROBE_MODES=8 bash runs/zb_submit.sh runs/gflat.tsv)" && rc=0 || rc=$?
[ "$rc" != 0 ] && ! echo "$o" | grep -q '^DRY\[' && ok "PROBE_UPDATE=1 with modes refused" || bad "update + modes"
o="$(dry env PROBE_UPDATE=1 START_ERR='C5=3' bash runs/zb_submit.sh runs/gflat.tsv)" && rc=0 || rc=$?
[ "$rc" != 0 ] && echo "$o" | grep -q 'START_ERR takes' && ok "START_ERR: an unknown term refused" || bad "START_ERR"
o="$(dry env START_SCALE='C50=0' bash runs/zb_submit.sh runs/gflat.tsv)" && rc=0 || rc=$?
[ "$rc" != 0 ] && ok "START_SCALE without PROBE_UPDATE refused" || bad "START_SCALE alone"
printf 'Pb\t80\tflat35\t-\tsingle\t0.1\t8\t200\t0.25\t0.2375\t0\t0\t-\t20\n' >"$T/flat35.tsv"
o="$(dry bash runs/zb_submit.sh "$T/flat35.tsv")" && rc=0 || rc=$?
[ "$rc" != 0 ] && echo "$o" | grep -q 'no tableau row for flat35' && ok "a flat angle with no tableau row refused" || bad "missing flat row"
grep -v '^#' runs/gflat.tsv | head -1 >"$T/a30.tsv"; Z="$T/z7"; mkdir -p "$Z"
real "$Z" STUB_SQUEUE="zb_sim_Pb_a030_ideal_df0_single_st0p1_w8_s0p25_bz40p000" bash runs/zb_submit.sh "$T/a30.tsv" && rc=0 || rc=$?
[ "$rc" != 0 ] && grep -q 'zb_sim_.* is still in the queue' "$Z/out.txt" && ok "a point whose sim is still in the queue refused (no second sim)" || bad "queued sim"
# -- DOSE (2026-11-12): a noise job between sim and recon, names carry _dose<D>, the recon reads the copy's patterns
o="$(dry env DOSE=1e5 GROUPING='32;16' bash runs/zb_submit.sh runs/gflat_fixed.tsv)" && rc=0 || rc=$?
nz=$(echo "$o" | grep 'job-name=zb_noise_Pb_a080_flat30_single_st0p1_w8_s0p25_bz40p000_rs90_dose1e5 ')
rz=$(echo "$o" | grep 'job-name=zb_rec_Pb_a080_flat30_single_st0p1_w8_s0p25_bz40p000_rs90_dose1e5_NL40 ')
nid=$(echo "$o" | grep -o 'DRY\[[0-9]*\] .*job-name=zb_noise_Pb_a080_flat30[^ ]* ' | sed -E 's/DRY\[([0-9]+)\].*/dry\1/')
[ "$rc" = 0 ] && [[ "$nz" == *"add_poisson_noise.py' --in-dir '/zbroot/sim_Pb_a080_flat30_single_st0p1_w8_s0p25_bz40p000_rs90' --out-dir '/zbroot/sim_Pb_a080_flat30_single_st0p1_w8_s0p25_bz40p000_rs90_dose1e5' --dose 1e5 --seed 0"* ]] \
  && [[ "$rz" == *"afterok:${nid} "* ]] && [ "$(echo "$o" | grep -c 'job-name=zb_noise_')" = 4 ] && [ "$(echo "$o" | grep -c 'job-name=zb_rec_.*_dose1e5_NL')" = 4 ] \
  && ok "DOSE=1e5: four noise jobs, each recon waits on its own, names _dose1e5" || { bad "DOSE wiring"; echo "$o" | grep -m3 -E 'NOT SUB|noise' | cut -c1-250; }
o="$(dry env DOSE=lots bash runs/zb_submit.sh runs/gflat_fixed.tsv)" && rc=0 || rc=$?
[ "$rc" != 0 ] && echo "$o" | grep -q 'DOSE is a fluence' && ok "DOSE: a non-number refused" || bad "DOSE parse"
[ "$fail" = 0 ] && echo "ALL OK" || { echo "SOME CHECKS FAILED"; echo "$out" | tail -5; exit 1; }
