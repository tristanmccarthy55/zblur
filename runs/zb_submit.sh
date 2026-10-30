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
# Points file (tab-separated, '#' comments): el alpha aber c1 geom step win det slice zvac fs src T nl_list
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

for f in sim/run_sim.slurm run_recon_synthetic_ML.slurm analysis/atomfind/extract_psf.py campaign/ceos_sweep.tsv; do
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
num() { echo "$1" | sed 's/\./p/'; }

[ "$DRY_RUN" = "1" ] || mkdir -p "${ZB_ROOT}/logs" "${ZB_ROOT}/psf" "${ZB_ROOT}/manifests"
[ "$DRY_RUN" = "1" ] || printf "sim\tnl\tdz_A\tzdrop\tboxz_A\tzvac_A\tel\talpha\taber\tc1_A\tgeom\tstep_A\twin_A\tdet_mrad\tdx_A\tslice_A\tfs_A\tsrc_A\tT_K\tbeta_lsq\tgrouping\tniter\trtag\tkernel\tsim_job\trecon_job\textract_job\tptycho_pin\tmex_cksum\tprobe_modes\tprobe_start\n" >"$MAN"

EXTRACT_JOBS=(); PSF_NAMES=(); SEEN=" "
while IFS=$'\t' read -r el alpha aber c1 geom step win det slice zvac fs src T nls rest; do
    case "${el:-}" in ''|\#*) continue;; esac
    [ -z "${rest:-}" ] || { echo "too many columns: $el $alpha ..." >&2; exit 1; }
    case "$el" in Pb|Ti) ;; *) echo "el must be Pb or Ti, not $el" >&2; exit 1;; esac
    case "$aber" in ideal|arm) ;; *) echo "aber must be ideal or arm, not $aber" >&2; exit 1;; esac
    case "$geom" in grid|single) ;; *) echo "geom must be grid or single, not $geom" >&2; exit 1;; esac
    [ "$aber" = "arm" ] && [ "$c1" != "-" ] && { echo "arm takes the tableau's own C1: put - in the c1 column" >&2; exit 1; }
    [ "$aber" = "ideal" ] && [ "$c1" = "-" ] && { echo "ideal needs a c1 (0 = in focus)" >&2; exit 1; }
    a3=$(printf '%03d' "$alpha")
    boxz=$(awk -v t="$THIN" -v c="$CELL_Z" -v v="$zvac" 'BEGIN{printf "%.3f", t*c + 2*v}')
    atomz=$(awk -v b="$boxz" 'BEGIN{printf "%.4f", b/2}')
    name="${el}_a${a3}_${aber}"; [ "$aber" = "ideal" ] && name="${name}_df$(num "$c1")"
    name="${name}_${geom}_st$(num "$step")_w$(num "$win")_s$(num "$slice")_bz$(num "$boxz")"
    [ "$fs" != "0" ] && name="${name}_fs$(num "$fs")"
    [ "$src" != "0" ] && name="${name}_src$(num "$src")"
    [ "$T" != "-" ] && name="${name}_T${T}v2"          # v2: simulated after the sigma fix (the old _T<T> moved sqrt(3) u)
    [ "$det" != "200" ] && name="${name}_d${det}"
    dx=$(awk -v l="$LAM" -v t="$det" 'BEGIN{printf "%.5f", l / (2 * t / 1000)}')    # object pixel = lambda / (2 theta_max)
    case "$SEEN" in *" ${name} "*) echo "duplicate point ${name} in ${PTS}" >&2; exit 1;; esac
    SEEN="${SEEN}${name} "
    sdir="${ZB_ROOT}/sim_${name}"

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

    if [ -e "${sdir}/01/data_dp.hdf5" ]; then
        sj=""; echo "  ${name}: sim exists, reused" >&2
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
        # the start probe: fixed = the true one; modes = the sim's mixed state (written by the job) or, coherent, the true
        # one with the engine's extra modes
        pstart=fixed
        if [ "$PROBE_MODES" -gt 1 ]; then
            if [ "$fs" != "0" ] || [ "$src" != "0" ]; then pstart=mixed; else pstart=seeded; fi
        fi
        if [ "$DRY_RUN" != "1" ]; then
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
        rj=$(submit --chdir="$PTYCHO_ROOT" --job-name="zb_rec_${name}_NL${nl}${RTAG}" --time="$rtime" --mem=96G \
             ${dep[@]+"${dep[@]}"} --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
             --export="$rexp" "${PTYCHO_ROOT}/run_recon_synthetic_ML.slurm")
        pname="${name}_NL${nl}${RTAG}"
        ej=$(submit --job-name="zb_psf_${pname}" --time=00:30:00 --mem=16G --cpus-per-task=2 --dependency="afterok:${rj}" \
             --output="${ZB_ROOT}/logs/zb_psf_%j.out" --error="${ZB_ROOT}/logs/zb_psf_%j.err" \
             --wrap="git -C '${PTYCHO_ROOT}' rev-parse --short HEAD > '${ZB_ROOT}/psf/psf_${pname}_pin.txt'; '${PYBIN}' '${PTYCHO_ROOT}/analysis/atomfind/extract_psf.py' '${rdir}' '${pname}' --zdrop ${zdrop} ${xargs} --out '${ZB_ROOT}/psf'")
        EXTRACT_JOBS+=("$ej"); PSF_NAMES+=("$pname")
        [ "$DRY_RUN" = "1" ] || printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$name" "$nl" "$dz" "$zdrop" "$boxz" "$zvac" "$el" "$alpha" "$aber" "$c1" "$geom" "$step" "$win" "$det" "$dx" "$slice" \
            "$fs" "$src" "$T" "$BETA_LSQ" "$GROUPING" "$NITER" "${RTAG:--}" "$pname" "${sj:-reused}" "$rj" "$ej" "$PIN" "$MEXSUM" \
            "$PROBE_MODES" "$pstart" >>"$MAN"
        echo "  ${pname}: dz ${dz} A, zdrop ${zdrop} (${geom}, zvac ${zvac}), dx ${dx} A; probe ${pstart}, ${PROBE_MODES} mode(s); sim ${sj:-reused}, recon ${rj}, kernel ${ej}" >&2
    done
done <"$PTS"

[ ${#EXTRACT_JOBS[@]} -gt 0 ] || { echo "no points in $PTS" >&2; exit 1; }
deps=$(IFS=:; echo "${EXTRACT_JOBS[*]}")
# _pin.txt: the ptychoshelves commit when the kernel was extracted, just after its reconstruction ran. MATLAB reads the
# engine when the job runs, not when it is submitted, so a pull in between would otherwise go unrecorded (it happened on
# 2026-10-23: 8edc045 at submission, 138deb4 by the time G0c ran; that change was inert without UPSAMPLING).
files=""; for p in "${PSF_NAMES[@]}"; do files="${files} psf/psf_${p}_vol.npy psf/psf_${p}_check.png psf/psf_${p}_pin.txt"; done
pack="${PACK_DIR}/zblur_${TAG}_${TS}.tgz"
pj=$(submit --job-name="zb_pack_${TAG}" --time=00:20:00 --mem=4G --cpus-per-task=1 --dependency="afterany:${deps}" \
     --output="${ZB_ROOT}/logs/zb_pack_%j.out" --error="${ZB_ROOT}/logs/zb_pack_%j.err" \
     --wrap="cd '${ZB_ROOT}' && cp '${MAN}' manifest_${TAG}_${TS}.tsv && tar czf '${pack}' --ignore-failed-read manifest_${TAG}_${TS}.tsv ${files}")
echo "ptychoshelves pin ${PIN} (LSQ MEX cksum ${MEXSUM}); beta_LSQ ${BETA_LSQ}, grouping ${GROUPING}, ${NITER} it, tag '${RTAG}'; manifest ${MAN}" >&2
echo "pack job ${pj} -> ${pack}" >&2
echo "PULL: scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:${PACK_DIR}/zblur_${TAG}_2*.tgz' ~/Desktop/zblur/" >&2
