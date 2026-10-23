#!/usr/bin/env bash
# @file test_submit_dry.sh
# @brief Dry-run runs/zb_submit.sh on a points file against a local ptychoshelves checkout and check the wiring:
#        one sim per point, one recon per slice count waiting on its sim, one kernel job per recon waiting on it, one
#        pack waiting on every kernel, unique names, and zdrop = round(zvac / dz).
#   bash tests/test_submit_dry.sh [runs/g0.tsv] [path/to/ptychoshelves-clean]
#   BETA_LSQ=0.02 GROUPING=64 RTAG=_x bash tests/test_submit_dry.sh runs/g0c_slices.tsv   # overrides are checked too
set -euo pipefail
cd "$(dirname "$0")/.."
PTS="${1:-runs/g0.tsv}"; PR="${2:-../ptychoshelves-clean}"
out="$(DRY_RUN=1 ALLOW_DIRTY=1 PTYCHO_ROOT="$PR" ZB_ROOT=/zbroot PACK_DIR=/pack SHARE=/share bash runs/zb_submit.sh "$PTS" 2>&1)"
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
echo "$out" | grep 'job-name=zb_sim_.*_single_' | grep -q 'GRID_SPACING=0' && ok "single atom has GRID_SPACING=0" || { [ "$ns" = 0 ] || bad "single GRID_SPACING"; }
echo "$out" | grep -q 'job-name=zb_sim_.*_ideal_' && { echo "$out" | grep 'job-name=zb_sim_.*_ideal_' | grep -q 'CS=0,C5=0,DEFOCUS=' \
  && ok "ideal probe: CS 0, C5 0, its own defocus" || bad "ideal probe exports"; }
echo "$out" | grep -q 'job-name=zb_sim_' && { echo "$out" | grep -q 'SLICE_THICKNESS=0.25' && ok "sim slice 0.25 A reaches the sim" || bad "sim slice"; }
[ "$fail" = 0 ] && echo "ALL OK" || { echo "SOME CHECKS FAILED"; echo "$out" | tail -5; exit 1; }
