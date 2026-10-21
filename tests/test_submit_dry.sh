#!/usr/bin/env bash
# @file test_submit_dry.sh
# @brief Dry-run runs/zb_submit.sh on a points file against a local ptychoshelves checkout and check the wiring:
#        one sim per point, one recon per slice count waiting on its sim, one kernel job per recon waiting on it, one
#        pack waiting on every kernel, unique names, and zdrop = round(zvac / dz).
#   bash tests/test_submit_dry.sh [runs/g0.tsv] [path/to/ptychoshelves-clean]
set -euo pipefail
cd "$(dirname "$0")/.."
PTS="${1:-runs/g0.tsv}"; PR="${2:-../ptychoshelves-clean}"
out="$(DRY_RUN=1 ALLOW_DIRTY=1 PTYCHO_ROOT="$PR" ZB_ROOT=/zbroot PACK_DIR=/pack SHARE=/share bash runs/zb_submit.sh "$PTS" 2>&1)"
fail=0; ok(){ echo "ok   $1"; }; bad(){ echo "FAIL $1"; fail=1; }
npts=$(grep -vc '^#' "$PTS" || true); nrec=$(grep -v '^#' "$PTS" | awk -F'\t' '{n += split($9, a, ",")} END{print n}')
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
echo "$out" | grep 'dz .* zdrop' | awk '{dz=$3; zd=$6; sub(",","",zd)} END{}' >/dev/null
bz=$(echo "$out" | grep -o 'GRID_BOX_Z=[0-9.]*' | head -1 | cut -d= -f2)
echo "$out" | grep -E 'dz [0-9.]+ A, zdrop [0-9]+' | sed -E 's/.*dz ([0-9.]+) A, zdrop ([0-9]+).*/\1 \2/' | \
  awk '{z = int(4 / $1 + 0.5); if (z != $2) {print "zdrop " $2 " at dz " $1 ", expected " z; e = 1}} END{exit e}' \
  && ok "zdrop = round(4 / dz) (box ${bz} A)" || bad "zdrop"
echo "$out" | grep -q 'SLICE_THICKNESS=0.25' && ok "sim slice 0.25 A reaches the sim" || bad "sim slice"
[ "$fail" = 0 ] && echo "ALL OK" || { echo "SOME CHECKS FAILED"; echo "$out" | tail -5; exit 1; }
