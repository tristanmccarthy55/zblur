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
# Points file (tab-separated, '#' comments): el alpha aber c1 geom step win slice zvac fs src T nl_list
#   el      Pb | Ti
#   alpha   convergence semi-angle (mrad)
#   aber    ideal (no aberrations) | arm (the ARM200F tableau row armf_a<alpha> of campaign/ceos_sweep.tsv)
#   c1      ideal only: the probe defocus in the simulator's (abTEM) convention, A; - for arm (the tableau's own)
#   geom    grid (the campaign's S1 kernel: one plane, 3 A spacing) | single (one atom alone in vacuum, as in Lei &
#           Wang, arXiv 2506.22352 SI Note 5)
#   step    scan step (A);  win  scan window (A)
#   slice   the SIMULATION's slice thickness (A)
#   zvac    vacuum above and below the 5-cell box (A); the box is 5 x 3.905 + 2 zvac, the atom plane at its centre
#   fs      chromatic focal spread (A rms), 0 = none
#   src     source size (A FWHM), 0 = none
#   T       phonons: - = none, 295 = room temperature (per species), else a Debye temperature in K (per species)
#   nl_list reconstruction slice counts, comma-separated
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

GRIDSP=3; BIN=2; CELL_Z=3.905; THIN=5; NITER=200; BETA_LSQ=0.05; GROUPING=32; PHONON_N=16

for f in sim/run_sim.slurm run_recon_synthetic_ML.slurm analysis/atomfind/extract_psf.py campaign/ceos_sweep.tsv; do
    [ -f "${PTYCHO_ROOT}/${f}" ] || { echo "PTYCHO_ROOT=${PTYCHO_ROOT} has no ${f}" >&2; exit 1; }
done
# the parts of ptychoshelves this uses must be the committed version, so the pin below means something
if ! git -C "$PTYCHO_ROOT" diff --quiet HEAD -- sim run_recon_synthetic_ML.slurm ptycho analysis/atomfind/extract_psf.py; then
    [ "${ALLOW_DIRTY:-0}" = "1" ] || { echo "PTYCHO_ROOT has uncommitted changes in sim/, ptycho/ or the recon/extract scripts; commit or ALLOW_DIRTY=1" >&2; exit 1; }
fi
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
[ "$DRY_RUN" = "1" ] || printf "sim\tnl\tdz_A\tzdrop\tboxz_A\tzvac_A\tel\talpha\taber\tc1_A\tgeom\tstep_A\twin_A\tslice_A\tfs_A\tsrc_A\tT_K\tsim_job\trecon_job\textract_job\tptycho_pin\n" >"$MAN"

EXTRACT_JOBS=(); PSF_NAMES=(); SEEN=" "
while IFS=$'\t' read -r el alpha aber c1 geom step win slice zvac fs src T nls rest; do
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
    [ "$T" != "-" ] && name="${name}_T${T}"
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
        # grid: drop the vacuum bands, crop +-30 px (the 3 A grid allows no more). single: keep the depth (drop 2 A
        # at each face, where the box edges collect junk) and crop +-80 px (3.9 A) so the x-z section keeps the cone.
        if [ "$geom" = "grid" ]; then
            zdrop=$(awk -v v="$zvac" -v d="$dz" 'BEGIN{printf "%d", v/d + 0.5}'); xargs="--half-xy 30"
        else
            zdrop=$(awk -v d="$dz" 'BEGIN{printf "%d", 2/d + 0.5}'); xargs="--single --half-xy 80"
        fi
        rdir="${ZB_ROOT}/recon_${name}_NL${nl}"
        if [ "$DRY_RUN" != "1" ]; then
            mkdir -p "${rdir}/01"
            for f in data_dp.hdf5 data_position.hdf5 sim_meta.mat; do ln -sf "${sdir}/01/${f}" "${rdir}/01/${f}"; done
            ln -sf "${sdir}/01/probe_initial_true.mat" "${rdir}/01/probe_initial.mat"    # the true probe, held fixed
        fi
        if [ "$nl" -le 20 ]; then rtime=10:00:00; elif [ "$nl" -le 40 ]; then rtime=24:00:00; else rtime=48:00:00; fi
        rj=$(submit --chdir="$PTYCHO_ROOT" --job-name="zb_rec_${name}_NL${nl}" --time="$rtime" --mem=96G \
             ${dep[@]+"${dep[@]}"} --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
             --export="ALL,NLAYERS=${nl},SIM_BASE=${rdir}/,REGLAYER=0,PROBE_MODES=1,NITER=${NITER},SAVE_EVERY=${NITER},GROUPING=${GROUPING},BETA_LSQ=${BETA_LSQ}" \
             "${PTYCHO_ROOT}/run_recon_synthetic_ML.slurm")
        pname="${name}_NL${nl}"
        ej=$(submit --job-name="zb_psf_${pname}" --time=00:30:00 --mem=16G --cpus-per-task=2 --dependency="afterok:${rj}" \
             --output="${ZB_ROOT}/logs/zb_psf_%j.out" --error="${ZB_ROOT}/logs/zb_psf_%j.err" \
             --wrap="'${PYBIN}' '${PTYCHO_ROOT}/analysis/atomfind/extract_psf.py' '${rdir}' '${pname}' --zdrop ${zdrop} ${xargs} --out '${ZB_ROOT}/psf'")
        EXTRACT_JOBS+=("$ej"); PSF_NAMES+=("$pname")
        [ "$DRY_RUN" = "1" ] || printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$name" "$nl" "$dz" "$zdrop" "$boxz" "$zvac" "$el" "$alpha" "$aber" "$c1" "$geom" "$step" "$win" "$slice" "$fs" "$src" "$T" "${sj:-reused}" "$rj" "$ej" "$PIN" >>"$MAN"
        echo "  ${pname}: dz ${dz} A, zdrop ${zdrop} (${geom}, zvac ${zvac}); sim ${sj:-reused}, recon ${rj}, kernel ${ej}" >&2
    done
done <"$PTS"

[ ${#EXTRACT_JOBS[@]} -gt 0 ] || { echo "no points in $PTS" >&2; exit 1; }
deps=$(IFS=:; echo "${EXTRACT_JOBS[*]}")
files=""; for p in "${PSF_NAMES[@]}"; do files="${files} psf/psf_${p}_vol.npy psf/psf_${p}_check.png"; done
pack="${PACK_DIR}/zblur_${TAG}_${TS}.tgz"
pj=$(submit --job-name="zb_pack_${TAG}" --time=00:20:00 --mem=4G --cpus-per-task=1 --dependency="afterany:${deps}" \
     --output="${ZB_ROOT}/logs/zb_pack_%j.out" --error="${ZB_ROOT}/logs/zb_pack_%j.err" \
     --wrap="cd '${ZB_ROOT}' && cp '${MAN}' manifest_${TAG}_${TS}.tsv && tar czf '${pack}' --ignore-failed-read manifest_${TAG}_${TS}.tsv ${files}")
echo "ptychoshelves pin ${PIN}; manifest ${MAN}" >&2
echo "pack job ${pj} -> ${pack}" >&2
echo "PULL: scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:${PACK_DIR}/zblur_${TAG}_2*.tgz' ~/Desktop/zblur/" >&2
