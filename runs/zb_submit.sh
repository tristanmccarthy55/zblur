#!/usr/bin/env bash
# @file zb_submit.sh
# @brief Submit single-atom (one-plane grid) depth-PSF runs on Blythe: sim -> one recon per slice count -> kernel -> pack.
#
#   bash runs/zb_submit.sh runs/g0.tsv              # on Blythe, from the zblur checkout
#   DRY_RUN=1 bash runs/zb_submit.sh runs/g0.tsv    # print every sbatch instead of submitting
#
# Reuses ptychoshelves-clean's own job scripts (sim/run_sim.slurm, run_recon_synthetic_ML.slurm) and kernel extractor
# (analysis/atomfind/extract_psf.py) from a checkout at PTYCHO_ROOT, unchanged, so the operator is the campaign's. The
# recon recipe is the campaign's S1 kernel recipe (aberration_experiment/PSF_KERNELS.md): a Pb or Ti grid on ONE depth
# plane, 3 A spacing, 20 A scan at 0.5 A, BIN 2, the true probe fixed, REGLAYER 0, BETA_LSQ 0.05, 200 iterations. Only
# the physics in the simulation and the slice count of the reconstruction change between runs.
#
# Every directory name encodes every value that changes the data (element, alpha, aberrations, sim slice, vacuum,
# focal spread, source, temperature) so two runs can never share a folder -- the campaign driver's sim folders do not
# encode SLICE or ZVAC, which is why this does not go through it. Results land under ZB_ROOT ($SHARE/phucrh/zblur_runs);
# the pack is $SHARE/phucrh/zblur_<points>_<timestamp>.tgz (kernels, check images, manifest).
#
# Points file (tab-separated, '#' comments): el alpha aber c1 geom step win det slice zvac fs src T nl_list [side]
#   el      Pb | Ti
#   alpha   convergence semi-angle (mrad)
#   aber    ideal (no aberrations) | arm (the ARM200F tableau row armf_a<alpha> of campaign/ceos_sweep.tsv)
#   c1      ideal only: the probe defocus in the simulator's (abTEM) convention, A; - for arm (the tableau's own)
#   geom    grid (the campaign's S1 kernel: one plane, 3 A spacing) | single (one atom alone in vacuum, as in Lei &
#           Wang, arXiv 2506.22352 SI Note 5)
#   step    scan step (A);  win  scan window (A)
#   det     the largest angle recorded (mrad); 200 = the campaign's recording (not passed to the sim, so those sims are
#           reused unchanged), anything else is passed as DETECTOR_MAX_ANGLE and named _d<det>
#   slice   the SIMULATION's slice thickness (A)
#   zvac    vacuum above and below the 5-cell box (A); the box is 5 x 3.905 + 2 zvac, the atom plane at its centre
#   fs      chromatic focal spread (A rms), 0 = none
#   src     source size (A FWHM), 0 = none
#   T       phonons: - = none, 295 = room temperature (per species), else a Debye temperature in K (per species)
#   nl_list reconstruction slice counts, comma-separated
#   side    optional: the simulated box's side (A), named _rs<side>; the solver window is side / 2 (see THE WINDOW)
#
# Per submission (environment): BETA_LSQ (0.05), GROUPING (32), NITER (200), RTAG (appended to every reconstruction and
# kernel name, so a rerun with other engine settings never lands in an earlier run's folder).
#
# PROBE_MODES (1; 2026-10-30, the user: give the lone atom the modes the main reconstructions get): above 1, the probe is
# RELEASED and reconstructed as that many incoherent modes, as the campaign's mode legs are (campaign/run_thin_atomfind.sh,
# kick_job with KICK_START=mixed): a point with focal spread or source size starts from the top PROBE_MODES eigenmodes of
# its own sim's partially coherent state (sim/make_mixed_probe.py, run in the job via PROBE_MIXED); a coherent point
# starts from the true probe with the engine's own extra modes. The probe updates from PROBE_START (40) in the presolve and
# PROBE_START2 (20) in the full engine, inside the aperture (PROBE_SUPPORT_FFT=1). RTAG defaults to _m<modes>, so a mode
# run never lands in a fixed-probe run's folder. (The mixed state holds the focal spread and source size, not vibration,
# as in the campaign's cooled mode legs.)
#   PROBE_MODES=8 NITER=250 bash runs/zb_submit.sh runs/g1m.tsv
#
# THE WINDOW (2026-11-04). A row's solver window is its box / BIN (BIN 2). Column 15, side (optional; absent or - = the 70.008 A box, so
# old files and their names are unchanged), makes the box a side x side A region and is named _rs<side>: side 210 is the campaign's 105 A
# window. Every ARM row before this date ran with no side: a 35 A window with 32 % of the 80 mrad probe outside it, and their rings
# are suspect (runs/PLAN.md Correction). So every row's share of the probe outside its window is computed from its own values
# (runs/window_loss.py) before anything is submitted, and a row above MAX_OUTSIDE_PCT (6) is refused. Why 6, from the crystal runs at
# the same probe [measured: HANDOFF_REALISTIC]: 3.3 % at 105 A clean, 5.1 % at 90 mrad / 140 A clean (92 / 50), 9.0 % at 70 A same scores
# but one corner damaged and columns broken, 15 % at 52.5 A failed. It is the last clean number with room for the helper's ~0.2 %, and a
# guard, not a physics limit: ALLOW_CUT=1 lets a deliberate control through, as runs/g1w.tsv does.
# A row with a side also gets its own resources from its own size, and is refused if it cannot fit: a 105 A window recorded to +-200 mrad
# is 2130 px, and at 6561 positions ~120 GB of patterns against a 192 GB node.
#
# RERUNS AND LOSSES (2026-11-04, the disk-quota failures). A recon can end COMPLETED with nothing saved, and a sim that died writing leaves
# a data_dp.hdf5 that exists: campaign/check_outputs.py opens them instead.
#   a sim that is not whole is refused, with the line that removes it;
#   a recon dir holding an earlier attempt is refused unless REDO=1, which clears its analysis/, checkpoints and stale kernel (a 3 KB
#   Niter200.mat is enough to fool extract_psf's fallback) and runs again on the same sim, same name. A recon found whole is not rerun:
#   with REDO=1 it only gets its kernel job. Other engine settings still take an RTAG. A recon still in the queue is refused.
#   each kernel job runs the check first, and fails if the recon is not whole; the pack exits non-zero if a kernel is missing or its tar failed.
#   PACK_ALSO="<kernel> ..." PACK_AFTER=<job id> adds kernels of an earlier submission to this pack, with their manifest rows.
# DISK. Before submitting (also under DRY_RUN with DISK_CHECK=1): refused if group physics on share-physics (mmlsquota) would have less than
# DISK_MARGIN_GB (200) under its soft limit after this submission and COMMITTED_GB (what submissions before it have still to write), or if a
# test write of DISK_PROBE_MB (2048) fails. CLEAN_RAW=1, optional and only if disk is short (the user keeps every clean sim by default): removes the
# raw patterns of a file's sims instead of submitting, only where every recon of the row is whole and its kernel exists.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PTS="${1:?usage: zb_submit.sh POINTS.tsv}"
[ -f "$PTS" ] || { echo "no points file $PTS" >&2; exit 1; }
PTYCHO_ROOT="${PTYCHO_ROOT:-/springbrook/share/physics/phucrh/ptyco_baseline/ptyco_pipeline}"
ZB_ROOT="${ZB_ROOT:-${SHARE:-/springbrook/share/physics}/phucrh/zblur_runs}"
PACK_DIR="${PACK_DIR:-${SHARE:-/springbrook/share/physics}/phucrh}"
PYBIN="${CONDA_ENV:-${SHARE:-/springbrook/share/physics}/phucrh/envs/abtem}/bin/python"
TSV="${PTYCHO_ROOT}/campaign/ceos_sweep.tsv"
DRY_RUN="${DRY_RUN:-0}"
CHECK="${PTYCHO_ROOT}/campaign/check_outputs.py"
GUARD_PY="$PYBIN"; [ -x "$GUARD_PY" ] || GUARD_PY="${ZB_GUARD_PY:-python3}"     # the login-node checks; the jobs use PYBIN
MAX_OUTSIDE_PCT="${MAX_OUTSIDE_PCT:-6}"; ALLOW_CUT="${ALLOW_CUT:-0}"; REDO="${REDO:-0}"; CLEAN_RAW="${CLEAN_RAW:-0}"
PACK_ALSO="${PACK_ALSO:-}"; PACK_AFTER="${PACK_AFTER:-}"; COMMITTED_GB="${COMMITTED_GB:-0}"

GRIDSP=3; BIN=2; CELL_Z=3.905; THIN=5; PHONON_N=16; LAM=0.0196877
NITER="${NITER:-200}"; BETA_LSQ="${BETA_LSQ:-0.05}"; GROUPING="${GROUPING:-32}"; RTAG="${RTAG:-}"
PROBE_MODES="${PROBE_MODES:-1}"; PROBE_START="${PROBE_START:-40}"; PROBE_START2="${PROBE_START2:-20}"
case "$PROBE_MODES" in ''|*[!0-9]*|0) echo "PROBE_MODES must be a whole number >= 1, not '${PROBE_MODES}'" >&2; exit 1;; esac
[ "$PROBE_MODES" -gt 1 ] && [ -z "$RTAG" ] && RTAG="_m${PROBE_MODES}"
# Phonon points need the frozen-phonon sigma fix (ptychoshelves-clean a8b9264, 2026-10-30): before it, abTEM 1.0.9 moved
# every atom sqrt(3) too far. Their sims are named _T<T>v2 so the pre-fix _T<T> folders are never reused, and they are
# refused unless the checkout contains the fix. Checked before anything is submitted.
PHONON_FIX="${PHONON_FIX:-a8b9264}"
ph=$(awk -F'\t' '$0 !~ /^#/ && NF >= 13 && $13 != "-" {print $1, $2, $3, "T=" $13}' "$PTS")
if [ -n "$ph" ] && ! git -C "$PTYCHO_ROOT" merge-base --is-ancestor "$PHONON_FIX" HEAD 2>/dev/null; then
    printf "phonon points need the frozen-phonon sigma fix %s in PTYCHO_ROOT (%s); pull it there first and run\nsim/test_phonon_sigma.py with the sim env. Points:\n%s\n" "$PHONON_FIX" "$PTYCHO_ROOT" "$ph" >&2
    exit 1
fi

for f in sim/run_sim.slurm run_recon_synthetic_ML.slurm analysis/atomfind/extract_psf.py campaign/ceos_sweep.tsv campaign/check_outputs.py; do
    [ -f "${PTYCHO_ROOT}/${f}" ] || { echo "PTYCHO_ROOT=${PTYCHO_ROOT} has no ${f}" >&2; exit 1; }
done
# the SOURCE this uses must be the committed version, so the pin below means something. Compiled MEX binaries are not
# checked (Blythe rebuilds get_LSQ_step_mex.mexa64 for its own MATLAB/CUDA, which every campaign run also uses); their
# checksum goes into the manifest instead, so the exact engine binary is still on record.
if ! git -C "$PTYCHO_ROOT" diff --quiet HEAD -- sim run_recon_synthetic_ML.slurm ':(glob)ptycho/**/*.m' analysis/atomfind/extract_psf.py; then
    [ "${ALLOW_DIRTY:-0}" = "1" ] || { echo "PTYCHO_ROOT has uncommitted changes in sim/, ptycho/*.m or the recon/extract scripts:" >&2
        git -C "$PTYCHO_ROOT" diff --stat HEAD -- sim run_recon_synthetic_ML.slurm ':(glob)ptycho/**/*.m' analysis/atomfind/extract_psf.py >&2
        echo "commit them, or ALLOW_DIRTY=1 if you know why" >&2; exit 1; }
fi
MEX="${PTYCHO_ROOT}/ptycho/+engines/+GPU/+LSQML/private/get_LSQ_step_mex.mexa64"
MEXSUM="$( [ -f "$MEX" ] && cksum "$MEX" | awk '{print $1}' || echo none)"
PIN="$(git -C "$PTYCHO_ROOT" rev-parse --short HEAD)"
num() { echo "$1" | sed 's/\./p/'; }

# ---- one row's names and sizes, for the checks and for the submission alike
row_setup(){
    a3=$(printf '%03d' "$alpha")
    boxz=$(awk -v t="$THIN" -v c="$CELL_Z" -v v="$zvac" 'BEGIN{printf "%.3f", t*c + 2*v}')
    atomz=$(awk -v b="$boxz" 'BEGIN{printf "%.4f", b/2}')
    name="${el}_a${a3}_${aber}"; [ "$aber" = "ideal" ] && name="${name}_df$(num "$c1")"
    name="${name}_${geom}_st$(num "$step")_w$(num "$win")_s$(num "$slice")_bz$(num "$boxz")"
    [ "$side" != "-" ] && name="${name}_rs$(num "$side")"
    [ "$fs" != "0" ] && name="${name}_fs$(num "$fs")"
    [ "$src" != "0" ] && name="${name}_src$(num "$src")"
    [ "$T" != "-" ] && name="${name}_T${T}v2"          # v2: simulated after the sigma fix (the old _T<T> moved sqrt(3) u)
    [ "$det" != "200" ] && name="${name}_d${det}"
    dx=$(awk -v l="$LAM" -v t="$det" 'BEGIN{printf "%.5f", l / (2 * t / 1000)}')    # object pixel = lambda / (2 theta_max)
    sdir="${ZB_ROOT}/sim_${name}"
    box=70.008; [ "$side" != "-" ] && box="$side"
    winA=$(awk -v b="$box" -v n="$BIN" 'BEGIN{printf "%.3f", b / n}')
    # the size of its patterns (the simulator crops the detector to a multiple of 2 BIN); positions as abTEM 1.0.9 counts them, one more per axis
    Ndp=$(awk -v t="$det" -v b="$box" -v l="$LAM" -v n="$BIN" 'BEGIN{u = int(2 * t / 1000 * b / l); print int(u / (2 * n)) * 2}')
    npos=$(awk -v w="$win" -v s="$step" 'BEGIN{m = w / s; k = int(m + 0.999999) + 1; print k * k}')
    dataGB=$(awk -v p="$npos" -v n="$Ndp" 'BEGIN{printf "%.1f", p * n * n * 4 / 1e9}')
}
recon_state(){ # $1 recon dir -> none | whole | leftover. A dir holding only the links a submission made is none
    local d="$1"
    [ -d "$d/analysis" ] || ls "$d"/slurm_*.out >/dev/null 2>&1 || ls -d "$d"/01/MLs_* >/dev/null 2>&1 || { echo none; return; }
    if "$GUARD_PY" "$CHECK" recon "$d" >/dev/null 2>&1; then echo whole; else echo leftover; fi
}
refuse(){ echo "NOT SUBMITTED: $*" >&2; exit 1; }
"$GUARD_PY" -c 'import numpy' 2>/dev/null || refuse "no python with numpy at ${GUARD_PY} for the window guard (set CONDA_ENV, or ZB_GUARD_PY)"

# ---- CLEAN_RAW=1: the raw patterns of this file's sims go, only where every recon of the row is whole and its kernel exists. Nothing is submitted.
if [ "$CLEAN_RAW" = "1" ]; then
    while IFS=$'\t' read -r el alpha aber c1 geom step win det slice zvac fs src T nls side rest; do
        case "${el:-}" in ''|\#*) continue;; esac
        side="${side:--}"; row_setup; ok=1
        IFS=',' read -r -a NLS <<<"$nls"
        for nl in "${NLS[@]}"; do
            [ "$(recon_state "${ZB_ROOT}/recon_${name}_NL${nl}${RTAG}")" = whole ] || { echo "  ${name}: NL${nl}${RTAG} not whole, its sim stays" >&2; ok=0; }
            [ -s "${ZB_ROOT}/psf/psf_${name}_NL${nl}${RTAG}_vol.npy" ] || { echo "  ${name}: no kernel for NL${nl}${RTAG}, its sim stays" >&2; ok=0; }
        done
        if [ "$ok" = 1 ] && [ -e "${sdir}/01/data_dp.hdf5" ]; then
            echo "  ${name}: removing $(du -sh "${sdir}/01/data_dp.hdf5" | cut -f1) of patterns (probes and sim_meta stay)" >&2
            rm -f "${sdir:?}/01/data_dp.hdf5" "${sdir:?}/01/data_position.hdf5"
        fi
    done <"$PTS"
    exit 0
fi

# ---- PREFLIGHT: every row, before anything is written or submitted
QUEUED=""; if command -v squeue >/dev/null 2>&1; then QUEUED="$(squeue --me -h -o %j 2>/dev/null)" || { [ "$REDO" = "1" ] && { echo "squeue failed: with REDO=1 a running recon cannot be told from a dead one" >&2; exit 1; }; }; fi
SEEN=" "; LOSSES=""; CUT=""; NEED=0
while IFS=$'\t' read -r el alpha aber c1 geom step win det slice zvac fs src T nls side rest; do
    case "${el:-}" in ''|\#*) continue;; esac
    [ -z "${rest:-}" ] || refuse "too many columns: $el $alpha ..."
    side="${side:--}"
    case "$el" in Pb|Ti) ;; *) refuse "el must be Pb or Ti, not $el";; esac
    case "$aber" in ideal|arm) ;; *) refuse "aber must be ideal or arm, not $aber";; esac
    case "$geom" in grid|single) ;; *) refuse "geom must be grid or single, not $geom";; esac
    [ "$aber" = "arm" ] && [ "$c1" != "-" ] && refuse "arm takes the tableau's own C1: put - in the c1 column"
    [ "$aber" = "ideal" ] && [ "$c1" = "-" ] && refuse "ideal needs a c1 (0 = in focus)"
    if [ "$side" != "-" ]; then
        case "$side" in *[!0-9.]*|'') refuse "side must be a length in A or -, not '${side}'";; esac
    fi
    row_setup
    case "$SEEN" in *" ${name} "*) refuse "duplicate point ${name} in ${PTS}";; esac
    SEEN="${SEEN}${name} "
    # the window, from the row's own values
    if [ "$aber" = "arm" ]; then
        trow=$(awk -F'\t' -v L="armf_a${a3}" '$1==L' "$TSV")
        [ -n "$trow" ] || refuse "no armf_a${a3} row in ${TSV}"
        gc1=$(echo "$trow" | cut -f5); gaj=$(echo "$trow" | cut -f9)
    else
        gc1="$c1"; gaj='{"C30":0.0,"C50":0.0}'
    fi
    key="${alpha}|${gc1}|${winA}"
    pct=$(printf '%s\n' "$LOSSES" | awk -F'=' -v k="$key" '$1 == k {print $2; exit}')
    if [ -z "$pct" ]; then
        pct=$("$GUARD_PY" "${HERE}/window_loss.py" --alpha "$alpha" --defocus="$gc1" --aber "$gaj" --window "$winA") || refuse "window_loss.py failed for ${name}"
        LOSSES="${LOSSES}${key}=${pct}
"
    fi
    echo "  ${name}: window ${winA} A, ${pct} % of the probe outside" >&2
    if awk -v p="$pct" -v m="$MAX_OUTSIDE_PCT" 'BEGIN{exit !(p > m)}'; then
        CUT="${CUT}
    ${name}: ${pct} % outside its ${winA} A window"
    fi
    # a row with a side: does it fit, at its own size?
    if [ "$side" != "-" ]; then
        peak=$(awk -v d="$dataGB" 'BEGIN{printf "%.0f", 4.3 * d}')        # a recon peaks at ~4.3 x its data [measured: 1305041]
        awk -v p="$peak" 'BEGIN{exit !(p > 160)}' && refuse "${name}: ${Ndp} px x ${npos} positions = ${dataGB} GB of patterns, a recon peak of ~${peak} GB on a 192 GB node (5776 positions at 1420 px, ~190 GB, did not fit; 4096, ~135 GB, did). Fewer positions or a smaller recorded angle."
        gfull=$(echo "$GROUPING" | sed 's/.*[; ]//')
        if [ "$Ndp" -gt 800 ] && [ "$gfull" -gt 16 ]; then refuse "${name}: ${Ndp} px needs a full-engine grouping <= 16 on an L40 (it is ${gfull}): GROUPING=\"32;16\""; fi
    fi
    # the sim: whole, or refused; never a truncated file taken for a finished sim
    if [ -e "${sdir}/01/data_dp.hdf5" ]; then
        "$GUARD_PY" "$CHECK" sim "$sdir" >/dev/null 2>&1 || refuse "${sdir} is not a whole sim: $("$GUARD_PY" "$CHECK" sim "$sdir" 2>&1 | head -1). Remove it, then submit again: rm -rf '${sdir}'"
    else
        NEED=$(awk -v n="$NEED" -v d="$dataGB" 'BEGIN{printf "%.0f", n + d + 0.5}')
    fi
    # the recons: no rerun lands on a failed attempt's leftovers
    IFS=',' read -r -a NLS <<<"$nls"
    for nl in "${NLS[@]}"; do
        rdir="${ZB_ROOT}/recon_${name}_NL${nl}${RTAG}"
        case "
${QUEUED}
" in *"
zb_rec_${name}_NL${nl}${RTAG}
"*) refuse "zb_rec_${name}_NL${nl}${RTAG} is still in the queue";; esac
        st=$(recon_state "$rdir")
        case "$st:$REDO" in
            leftover:1|whole:1|none:*) ;;
            whole:*) refuse "${rdir} already holds a whole reconstruction. REDO=1 only makes its kernel; other engine settings take an RTAG";;
            *) refuse "${rdir} holds an earlier attempt that is not whole: $("$GUARD_PY" "$CHECK" recon "$rdir" 2>&1 | head -1). REDO=1 clears it and runs again on the same sim";;
        esac
        [ "$st" = whole ] || NEED=$(awk -v n="$NEED" 'BEGIN{print n + 3}')     # a recon dir: 2-4 GB [records], kernels ~1 MB
    done
done <"$PTS"
if [ -n "$CUT" ] && [ "$ALLOW_CUT" != "1" ]; then
    refuse "these rows' windows would cut the probe (above ${MAX_OUTSIDE_PCT} %):${CUT}
Crystal at this probe: 3.3 % clean, 5.1 % clean, 9.0 % damaged, 15 % failed. Give the row a side (column 15: 210 is a 105 A window), or ALLOW_CUT=1 for a deliberate control."
fi
[ -z "$CUT" ] || echo "  ALLOW_CUT=1: these rows run with the probe cut, on purpose:${CUT}" >&2
echo "  this submission writes ~${NEED} GB" >&2
if [ "$DRY_RUN" != "1" ] || [ "${DISK_CHECK:-0}" = "1" ]; then
    row=$(mmlsquota -g physics --block-size auto 2>/dev/null | grep 'share-physics ' | head -1) || row=""
    read -r qu qs <<<"$(echo "$row" | awk 'function gb(x,  n, c) {n = x + 0; c = substr(x, length(x)); return (c == "T") ? n * 1099.5 : (c == "G") ? n * 1.0737 : (c == "M") ? n * 0.00105 : n * 1.024e-6}
      {k = 0; for (i = 1; i <= NF; i++) if ($i ~ /^[0-9.]+[KMGTP]?$/) {k++; if (k == 1) a = $i; if (k == 2) b = $i}} END {if (k >= 2) printf "%.0f %.0f", gb(a), gb(b)}')"
    [ -n "${qu:-}" ] && [ -n "${qs:-}" ] && [ "${qs}" != 0 ] || refuse "could not read used / soft from 'mmlsquota -g physics --block-size auto | grep share-physics': '${row}'"
    left=$(( qs - qu - COMMITTED_GB - NEED ))
    echo "  disk: group physics on share-physics ${qu} GB used of ${qs} GB soft; this writes ~${NEED} GB, earlier parts ${COMMITTED_GB} GB still to come; leaves ${left} GB (margin ${DISK_MARGIN_GB:-200} GB)" >&2
    [ "$left" -ge "${DISK_MARGIN_GB:-200}" ] || refuse "${left} GB would be left under the soft limit, under the ${DISK_MARGIN_GB:-200} GB margin: free space first"
    probe="${DISK_PROBE_DIR:-${ZB_ROOT}}/.zb_write_test_$$"
    "$GUARD_PY" -c 'import os, sys
f, n = sys.argv[1], int(sys.argv[2])
b = os.urandom(1 << 20)
try:
    with open(f, "wb") as h:
        for _ in range(n):
            h.write(b)
        h.flush(); os.fsync(h.fileno())
finally:
    if os.path.exists(f):
        os.remove(f)' "$probe" "${DISK_PROBE_MB:-2048}" >/dev/null 2>&1 || refuse "the test write of ${DISK_PROBE_MB:-2048} MB to ${probe} failed"
    echo "  disk: test write of ${DISK_PROBE_MB:-2048} MB ok" >&2
fi

TS="$(date +%Y%m%d_%H%M%S)"; TAG="$(basename "$PTS" .tsv)"
MAN="${ZB_ROOT}/manifests/${TAG}_${TS}.tsv"

# sbatch FROM the ptychoshelves checkout: its job scripts find their code through SLURM_SUBMIT_DIR, which is where sbatch
# was called from (not --chdir) -- the campaign driver cd's to its repo for the same reason. Under DRY_RUN, print it and
# hand out numbered fake ids (a counter file, since each call runs in a $( ) subshell) so the dependencies can be read.
if [ "$DRY_RUN" = "1" ]; then DRY_CTR="$(mktemp "${TMPDIR:-/tmp}/zbdry.XXXXXX")"; echo 0 >"$DRY_CTR"; trap 'rm -f "$DRY_CTR"' EXIT; fi
submit() {
    if [ "$DRY_RUN" = "1" ]; then
        local n; n=$(( $(cat "$DRY_CTR") + 1 )); echo "$n" >"$DRY_CTR"
        echo "DRY[${n}] (cd ${PTYCHO_ROOT}) sbatch $*" >&2; echo "dry${n}"
    else
        ( cd "$PTYCHO_ROOT" && sbatch --parsable "$@" )
    fi
}
[ "$DRY_RUN" = "1" ] || mkdir -p "${ZB_ROOT}/logs" "${ZB_ROOT}/psf" "${ZB_ROOT}/manifests"
[ "$DRY_RUN" = "1" ] || printf "sim\tnl\tdz_A\tzdrop\tboxz_A\tzvac_A\tel\talpha\taber\tc1_A\tgeom\tstep_A\twin_A\tdet_mrad\tdx_A\tslice_A\tfs_A\tsrc_A\tT_K\tbeta_lsq\tgrouping\tniter\trtag\tkernel\tsim_job\trecon_job\textract_job\tptycho_pin\tmex_cksum\tprobe_modes\tprobe_start\tside_A\twindow_A\toutside_pct\n" >"$MAN"

EXTRACT_JOBS=(); PSF_NAMES=()
while IFS=$'\t' read -r el alpha aber c1 geom step win det slice zvac fs src T nls side rest; do
    case "${el:-}" in ''|\#*) continue;; esac
    side="${side:--}"
    row_setup       # (every row passed the preflight above)

    # ---- aberrations
    if [ "$aber" = "arm" ]; then
        row=$(awk -F'\t' -v L="armf_a${a3}" '$1==L' "$TSV")
        [ -n "$row" ] || { echo "no armf_a${a3} row in ${TSV}" >&2; exit 1; }
        c5=$(echo "$row" | cut -f3); c3=$(echo "$row" | cut -f4); c1=$(echo "$row" | cut -f5); aj=$(echo "$row" | cut -f9)
    else
        c5=0; c3=0; aj='{"C30":0.0,"C50":0.0}'
    fi
    export ABERRATIONS_JSON="$aj"     # has commas: read from the environment by run_sim.slurm (--export=ALL)

    # ---- physics
    exp="ALL,JOB_DIR=${sdir},SLICE_THICKNESS=${slice},SCAN_STEP=${step},CONVERGENCE=${alpha},PHONON_SIGMA=0.08,PHONON_SEED=1"
    if [ "$T" = "-" ]; then
        exp="${exp},PHONONS=0,PER_SPECIES_SIGMA=0"; stime=24:00:00     # 0.25 A slices: ~3.6x the 0.9 A sims
    elif [ "$T" = "295" ]; then
        exp="${exp},PHONONS=${PHONON_N},PER_SPECIES_SIGMA=1"; stime=24:00:00
    else
        exp="${exp},PHONONS=${PHONON_N},PER_SPECIES_SIGMA=1,PHONON_TEMPERATURE=${T}"; stime=24:00:00
    fi
    [ "$fs" != "0" ] && exp="${exp},FOCAL_SPREAD=${fs}"
    [ "$src" != "0" ] && exp="${exp},SOURCE_FWHM=${src}"
    exp="${exp},BIN_FACTOR=${BIN},RECON_FULL_BOX=1,Z_VACUUM=${zvac},ABERRATED=1,PROBE_INITIAL=nominal"
    exp="${exp},CS=${c3},C5=${c5},DEFOCUS=${c1},OVERWRITE=0,SCAN_WINDOW=${win}"
    if [ "$geom" = "grid" ]; then gsp="$GRIDSP"; else gsp=0; fi        # GRID_SPACING 0 = one atom (build_single_atom)
    exp="${exp},SINGLE_ATOM=${el},ATOM_Z=${atomz},GRID_SPACING=${gsp},GRID_BOX_Z=${boxz}"
    [ "$det" != "200" ] && exp="${exp},DETECTOR_MAX_ANGLE=${det}"
    [ "$side" != "-" ] && exp="${exp},REGION_SIDE=${side}"
    if [ "$side" != "-" ]; then
        # a region's sim, from its own size: 5.2e-10 s per position, slice and grid point [measured: 17.5 min a configuration, 1600 positions x
        # 31 slices x 6358^2; 19 min, 6561 x 80 x 2122^2], twice that, at least 6 h and at most the 48 h limit
        stime=$(awk -v p="$npos" -v b="$boxz" -v s="$slice" -v r="$side" -v c="$( [ "$T$fs$src" = "-00" ] && echo 1 || echo "$PHONON_N" )" \
                'BEGIN{sl = int(b / s + 0.999999); h = 2 * 5.2e-10 * p * sl * (r / 0.03303) ^ 2 * c / 3600; if (h < 6) h = 6; if (h > 48) h = 48; printf "%02d:00:00", h + 0.999}')
    fi

    if [ -e "${sdir}/01/data_dp.hdf5" ]; then
        sj=""; echo "  ${name}: sim whole (checked), reused" >&2
    else
        [ "$DRY_RUN" = "1" ] || mkdir -p "$sdir"
        sj=$(submit --chdir="$PTYCHO_ROOT" --job-name="zb_sim_${name}" --time="$stime" \
             --output="${ZB_ROOT}/logs/zb_sim_%j.out" --error="${ZB_ROOT}/logs/zb_sim_%j.err" \
             --export="$exp" "${PTYCHO_ROOT}/sim/run_sim.slurm")
    fi
    dep=(); [ -n "$sj" ] && dep=(--dependency="afterok:${sj}")

    # ---- one reconstruction per slice count, each followed by its kernel extraction
    IFS=',' read -r -a NLS <<<"$nls"
    for nl in "${NLS[@]}"; do
        dz=$(awk -v b="$boxz" -v n="$nl" 'BEGIN{printf "%.4f", b/n}')
        # grid: drop the vacuum bands, crop +-1.48 A (30 px at the campaign's 0.0492 A pixel; the 3 A grid allows no
        # more). single: keep the depth (drop 2 A at each face, where the box edges collect junk) and crop +-3.9 A so
        # the x-z section keeps the cone. The pixel follows the recorded angle, so the crops are set in A.
        if [ "$geom" = "grid" ]; then
            zdrop=$(awk -v v="$zvac" -v d="$dz" 'BEGIN{printf "%d", v/d + 0.5}')
            xargs="--half-xy $(awk -v d="$dx" 'BEGIN{printf "%d", 1.476/d + 0.5}') --dx ${dx}"
        else
            zdrop=$(awk -v d="$dz" 'BEGIN{printf "%d", 2/d + 0.5}')
            xargs="--single --half-xy $(awk -v d="$dx" 'BEGIN{printf "%d", 3.937/d + 0.5}') --dx ${dx}"
        fi
        rdir="${ZB_ROOT}/recon_${name}_NL${nl}${RTAG}"
        pname="${name}_NL${nl}${RTAG}"
        state=$(recon_state "$rdir")
        if [ "$state" = leftover ] && [ "$DRY_RUN" != "1" ]; then    # REDO=1 (the preflight refused it otherwise): nothing of a failed attempt stays to fool a rerun
            echo "  ${pname}: clearing an earlier attempt: analysis/, $(find "${rdir}/01" -name 'Niter*.mat' 2>/dev/null | wc -l | tr -d ' ') checkpoint(s), any kernel; slurm logs kept" >&2
            rm -rf "${rdir:?}/analysis"; find "${rdir}/01" -maxdepth 1 -type d -name 'MLs_*' -exec rm -rf {} + 2>/dev/null || true
            find "${rdir}/01" \( -name 'Niter*.mat' -o -name 'exit_wave*.mat' \) -type f -delete 2>/dev/null || true
            rm -f "${ZB_ROOT}/psf/psf_${pname}_vol.npy" "${ZB_ROOT}/psf/psf_${pname}_check.png" "${ZB_ROOT}/psf/psf_${pname}_pin.txt"
        fi
        # the start probe: fixed = the true one; modes = the sim's mixed state (written by the job) or, coherent, the true
        # one with the engine's extra modes
        pstart=fixed
        if [ "$PROBE_MODES" -gt 1 ]; then
            if [ "$fs" != "0" ] || [ "$src" != "0" ]; then pstart=mixed; else pstart=seeded; fi
        fi
        if [ "$DRY_RUN" != "1" ] && [ "$state" != whole ]; then
            mkdir -p "${rdir}/01"
            for f in data_dp.hdf5 data_position.hdf5 sim_meta.mat; do ln -sf "${sdir}/01/${f}" "${rdir}/01/${f}"; done
            if [ "$pstart" = mixed ]; then
                # the job writes probe_initial.mat HERE: it must not be a link, or the writer would write through it and
                # replace the sim's own true probe. The writer reads the true probe for the record.
                ln -sf "${sdir}/01/probe_initial_true.mat" "${rdir}/01/probe_initial_true.mat"
                rm -f "${rdir:?}/01/probe_initial.mat" "${rdir:?}/01/probe_initial.json"
            else
                ln -sf "${sdir}/01/probe_initial_true.mat" "${rdir}/01/probe_initial.mat"    # the true probe (held fixed when 1 mode)
            fi
        fi
        rexp="ALL,NLAYERS=${nl},SIM_BASE=${rdir}/,REGLAYER=0,PROBE_MODES=${PROBE_MODES},NITER=${NITER},SAVE_EVERY=${NITER},GROUPING=${GROUPING},BETA_LSQ=${BETA_LSQ}"
        if [ "$PROBE_MODES" -gt 1 ]; then
            rexp="${rexp},PROBE_START=${PROBE_START},PROBE_START2=${PROBE_START2},PROBE_SUPPORT_FFT=1"
            [ "$pstart" = mixed ] && rexp="${rexp},PROBE_MIXED=${PROBE_MODES}"
        fi
        if [ "$nl" -le 20 ]; then rtime=10:00:00; elif [ "$nl" -le 40 ]; then rtime=24:00:00; else rtime=48:00:00; fi
        [ "$PROBE_MODES" -gt 1 ] && rtime=48:00:00         # modes multiply the work: the crystal's 8-mode legs ran ~2.4x
        rmem=96G
        if [ "$side" != "-" ]; then
            # a recon's peak, 1.2 x, within the node [as the driver's mem_for]; time from the crystal's full engine, 8.56e8 position x px^2 x layer per s
            # at grouping 16, + 32 % for its presolve [measured: 1323756: 4318 + 13628 s], 1.5 x, within the 48 h limit
            rmem=$(awk -v d="$dataGB" 'BEGIN{g = int(1.2 * 4.3 * d + 0.999); if (g < 96) g = 96; if (g > 175) g = 175; printf "%dG", g}')
            rtime=$(awk -v p="$npos" -v n="$Ndp" -v l="$nl" -v i="$NITER" 'BEGIN{h = 1.5 * 1.32 * i * p * n * n * l / 8.56e8 / 3600; if (h < 4) h = 4; if (h > 48) h = 48; printf "%02d:00:00", h + 0.999}')
        fi
        ejdep=()
        if [ "$state" = whole ]; then
            rj=""; echo "  ${pname}: recon whole (checked), kernel only" >&2
        else
            rj=$(submit --chdir="$PTYCHO_ROOT" --job-name="zb_rec_${name}_NL${nl}${RTAG}" --time="$rtime" --mem="$rmem" \
                 ${dep[@]+"${dep[@]}"} --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
                 --export="$rexp" "${PTYCHO_ROOT}/run_recon_synthetic_ML.slurm")
            ejdep=(--dependency="afterok:${rj}")
        fi
        # the check first: COMPLETED on the recon says nothing, and extract_psf would fall back to a checkpoint
        ej=$(submit --job-name="zb_psf_${pname}" --time=00:30:00 --mem=16G --cpus-per-task=2 ${ejdep[@]+"${ejdep[@]}"} \
             --output="${ZB_ROOT}/logs/zb_psf_%j.out" --error="${ZB_ROOT}/logs/zb_psf_%j.err" \
             --wrap="'${PYBIN}' '${CHECK}' recon '${rdir}' || exit 1; git -C '${PTYCHO_ROOT}' rev-parse --short HEAD > '${ZB_ROOT}/psf/psf_${pname}_pin.txt'; '${PYBIN}' '${PTYCHO_ROOT}/analysis/atomfind/extract_psf.py' '${rdir}' '${pname}' --zdrop ${zdrop} ${xargs} --out '${ZB_ROOT}/psf'")
        EXTRACT_JOBS+=("$ej"); PSF_NAMES+=("$pname")
        [ "$DRY_RUN" = "1" ] || printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$name" "$nl" "$dz" "$zdrop" "$boxz" "$zvac" "$el" "$alpha" "$aber" "$c1" "$geom" "$step" "$win" "$det" "$dx" "$slice" \
            "$fs" "$src" "$T" "$BETA_LSQ" "$GROUPING" "$NITER" "${RTAG:--}" "$pname" "${sj:-reused}" "${rj:-reused}" "$ej" "$PIN" "$MEXSUM" \
            "$PROBE_MODES" "$pstart" "$side" "$winA" "$(printf '%s\n' "$LOSSES" | awk -F'=' -v k="${alpha}|${c1}|${winA}" '$1 == k {print $2; exit}')" >>"$MAN"
        echo "  ${pname}: dz ${dz} A, zdrop ${zdrop} (${geom}, zvac ${zvac}), dx ${dx} A; probe ${pstart}, ${PROBE_MODES} mode(s); sim ${sj:-reused}, recon ${rj:-reused}, kernel ${ej}" >&2
    done
done <"$PTS"

[ ${#EXTRACT_JOBS[@]} -gt 0 ] || { echo "no points in $PTS" >&2; exit 1; }
deps=$(IFS=:; echo "${EXTRACT_JOBS[*]}")
[ -z "$PACK_AFTER" ] || deps="${deps}:${PACK_AFTER}"
# _pin.txt: the ptychoshelves commit when the kernel was extracted, just after its reconstruction ran. MATLAB reads the
# engine when the job runs, not when it is submitted, so a pull in between would otherwise go unrecorded (it happened on
# 2026-10-23: 8edc045 at submission, 138deb4 by the time G0c ran; that change was inert without UPSAMPLING).
for p in $PACK_ALSO; do          # a kernel of an earlier submission: into this pack, with its manifest row, so blobs.py can read it
    PSF_NAMES+=("$p")
    if [ "$DRY_RUN" != "1" ]; then
        mrow=$(cat "${ZB_ROOT}"/manifests/*.tsv 2>/dev/null | awk -F'\t' -v k="$p" '$24 == k' | tail -1)
        [ -n "$mrow" ] && printf '%s\n' "$mrow" >>"$MAN" || echo "  PACK_ALSO ${p}: no manifest row found" >&2
    fi
done
files=""; vols=""; for p in "${PSF_NAMES[@]}"; do files="${files} psf/psf_${p}_vol.npy psf/psf_${p}_check.png psf/psf_${p}_pin.txt"; vols="${vols} psf/psf_${p}_vol.npy"; done
pack="${PACK_DIR}/zblur_${TAG}_${TS}.tgz"
# a missing kernel or a failed tar fails the job, after the tar, so what did land still comes home and sacct says so
pj=$(submit --job-name="zb_pack_${TAG}" --time=00:20:00 --mem=4G --cpus-per-task=1 --dependency="afterany:${deps}" \
     --output="${ZB_ROOT}/logs/zb_pack_%j.out" --error="${ZB_ROOT}/logs/zb_pack_%j.err" \
     --wrap="cd '${ZB_ROOT}' && cp '${MAN}' manifest_${TAG}_${TS}.tsv && tar czf '${pack}' --ignore-failed-read manifest_${TAG}_${TS}.tsv ${files}; rc=\$?; miss=0; for f in ${vols}; do [ -s \"\$f\" ] || { echo \"MISSING \$f\"; miss=1; }; done; [ \$rc = 0 ] || echo \"tar exit \$rc\"; [ \$miss = 0 ] && [ \$rc = 0 ]")
echo "ptychoshelves pin ${PIN} (LSQ MEX cksum ${MEXSUM}); beta_LSQ ${BETA_LSQ}, grouping ${GROUPING}, ${NITER} it, tag '${RTAG}'; manifest ${MAN}" >&2
echo "pack job ${pj} -> ${pack}" >&2
echo "PULL: scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:${PACK_DIR}/zblur_${TAG}_2*.tgz' ~/Desktop/zblur/" >&2
