#!/bin/bash
set -euo pipefail

# TRIGEMINAL SYSTEM TRACTOGRAPHY - Samir Akeb (2022-2023)
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Arnaud Bore (2023-2024)
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Nasrin Rafiei (2025-2026)
#
# SECOND-ORDER ENSEMBLE VERSION
# Compatible with the newer original second-order ROI/BDO logic
# while keeping the organized ensemble structure.
#
# Organized pipeline:
#   0) pool first-order final spinal / remaining_cp bundles
#   1) prepare second-order seed masks and thalamus ROIs
#   2) run second-order tracking for each (step, theta) combo in ORIG
#   3) merge combo outputs per tracking role in ORIG
#   4) register merged tractograms once to MNI
#   5) prepare second-order MNI ROIs and cut masks
#   6) filter merged second-order tractograms into pathway bundles
#   7) cut filtered bundles with label masks
#   8) reject outliers and save final second-order bundles

usage() {
    cat <<EOF 1>&2
Usage:
  $(basename "$0") -s <subjects_parent_or_single_subject_dir> -m <ROIs_clean_dir> -o <out_dir>
                   [-f fa_threshold] [-t threads] [-g true|false]
                   [-p step_size] [-e theta_deg]
                   [--npv_spinal_long N] [--npv_spinal_short N] [--npv_thalamus N]

Tracking params:
  -p step_size    If set with -e, runs SINGLE combo (no ensemble)
  -e theta_deg    If set with -p, runs SINGLE combo (no ensemble)
  If -p/-e not provided, ensemble is used:
    step  = 0.1 0.5 1.0
    theta = 20 30 40

Second-order seed budgets (TOTAL per role; auto-split across combos):
  --npv_spinal_long N   default: 1000
  --npv_spinal_short N  default: 100
  --npv_thalamus N      default: 500
EOF
    exit 1
}

# -------------------------
# Parse short options first
# -------------------------
s=""
m=""
o=""
f=""
t=""
g=""
p=""
e=""

while getopts ":s:m:o:f:t:g:p:e:" args; do
    case "${args}" in
        s) s=${OPTARG} ;;
        m) m=${OPTARG} ;;
        o) o=${OPTARG} ;;
        f) f=${OPTARG} ;;
        t) t=${OPTARG} ;;
        g) g=${OPTARG} ;;
        p) p=${OPTARG} ;;
        e) e=${OPTARG} ;;
        *) usage ;;
    esac
done
shift $((OPTIND-1))

# -------------------------
# Parse optional long args
# -------------------------
npv_spinal_long_total=1000
npv_spinal_short_total=100
npv_thalamus_total=500

while [[ $# -gt 0 ]]; do
    case "$1" in
        --npv_spinal_long)
            npv_spinal_long_total="$2"
            shift 2
            ;;
        --npv_spinal_short)
            npv_spinal_short_total="$2"
            shift 2
            ;;
        --npv_thalamus)
            npv_thalamus_total="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" 1>&2
            usage
            ;;
    esac
done

if [[ -z "${s}" || -z "${m}" || -z "${o}" ]]; then
    usage
fi

subject_dir="${s}"
mni_dir="${m}"
out_dir="${o}"
fa_threshold="${f:-0.15}"
nb_threads="${t:-1}"

gpu=""
if [[ -n "${g}" ]]; then
    gpu="--use_gpu"
fi



trk_is_empty() {
    local f="$1"

    if [[ ! -f "${f}" ]]; then
        return 0
    fi

    local n_str
    n_str=$(scil_tractogram_count_streamlines "${f}" 2>/dev/null | grep -Eo '[0-9]+' | tail -n 1 || true)

    if [[ -z "${n_str}" || "${n_str}" -eq 0 ]]; then
        return 0
    else
        return 1
    fi
}





# -------------------------
# Ensemble grid
# -------------------------
if [[ -n "${p}" && -n "${e}" ]]; then
    step_list=("${p}")
    theta_list=("${e}")
else
    step_list=(0.1 0.5 1.0)
    theta_list=(20 30 40)
fi

n_combos=$(( ${#step_list[@]} * ${#theta_list[@]} ))

# Split total budgets across combos (ceiling division)
npv_spinal_long_per_combo=$(( (npv_spinal_long_total + n_combos - 1) / n_combos ))
npv_spinal_short_per_combo=$(( (npv_spinal_short_total + n_combos - 1) / n_combos ))
npv_thalamus_per_combo=$(( (npv_thalamus_total + n_combos - 1) / n_combos ))

echo "Folder subjects: ${subject_dir}"
echo "Folder MNI: ${mni_dir}"
echo "Output folder: ${out_dir}"
echo "FA threshold: ${fa_threshold}"
echo "GPU: ${gpu}"
echo "Threads: ${nb_threads}"
echo "Tracking grid: steps=${step_list[*]}  thetas=${theta_list[*]}"
echo "Second-order budgets per combo:"
echo "  spinal long:  ${npv_spinal_long_per_combo}  (from total ${npv_spinal_long_total})"
echo "  spinal short: ${npv_spinal_short_per_combo} (from total ${npv_spinal_short_total})"
echo "  thalamus:     ${npv_thalamus_per_combo}     (from total ${npv_thalamus_total})"

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="${nb_threads}"

# Detect if -s is one subject folder or a parent folder of many subjects
if [[ -d "${subject_dir}/tractoflow" ]]; then
    subject_list=("${subject_dir}")
else
    shopt -s nullglob
    subject_list=("${subject_dir}"/*/)
    shopt -u nullglob
fi

if [[ ${#subject_list[@]} -eq 0 ]]; then
    echo "ERROR: No subject folders found in ${subject_dir}"
    exit 1
fi

# ============================================================
# 0) Pool first-order final bundles used by second-order
# ============================================================
mkdir -p "${out_dir}/mni_space/tracking_first_order/final"

for nside in left right; do
    spinal_inputs=()
    remaining_inputs=()

    for nsub_path in "${subject_list[@]}"; do
        nsub=$(basename "${nsub_path}")
        spinal_file="${out_dir}/${nsub}/mni_space/tracking_first_order/final_merged/final/${nsub}_${nside}_spinal.trk"
        remaining_file="${out_dir}/${nsub}/mni_space/tracking_first_order/final_merged/final/${nsub}_${nside}_remaining_cp.trk"

        [[ -f "${spinal_file}" ]] && spinal_inputs+=("${spinal_file}")
        [[ -f "${remaining_file}" ]] && remaining_inputs+=("${remaining_file}")
    done

    if [[ ${#spinal_inputs[@]} -eq 0 ]]; then
        echo "ERROR: No first-order spinal bundles found for side ${nside}."
        exit 1
    fi

    if [[ ${#remaining_inputs[@]} -eq 0 ]]; then
        echo "ERROR: No first-order remaining_cp bundles found for side ${nside}."
        exit 1
    fi

    scil_tractogram_math union "${spinal_inputs[@]}" \
        "${out_dir}/mni_space/tracking_first_order/final/all_${nside}_spinal.trk" -f

    scil_tractogram_compute_density_map \
        "${out_dir}/mni_space/tracking_first_order/final/all_${nside}_spinal.trk" \
        "${out_dir}/mni_space/tracking_first_order/final/all_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
        --binary -f

    scil_tractogram_math union "${remaining_inputs[@]}" \
        "${out_dir}/mni_space/tracking_first_order/final/all_${nside}_remaining_cp.trk" -f

    scil_tractogram_compute_density_map \
        "${out_dir}/mni_space/tracking_first_order/final/all_${nside}_remaining_cp.trk" \
        "${out_dir}/mni_space/tracking_first_order/final/all_${nside}_remaining_cp_density_mni.nii.gz" \
        --binary -f
done

# ============================================================
# 1) Process each subject
# ============================================================
for nsub_path in "${subject_list[@]}"; do
    nsub=$(basename "${nsub_path}")

    mkdir -p "${out_dir}/${nsub}/orig_space/rois"
    mkdir -p "${out_dir}/${nsub}/orig_space/tracking_second_order"/{trials,merged,final}
    mkdir -p "${out_dir}/${nsub}/mni_space/rois"
    mkdir -p "${out_dir}/${nsub}/mni_space/tracking_second_order"/{orig,filtered,final,cut}

    orig_rois_dir="${out_dir}/${nsub}/orig_space/rois"
    mni_rois_dir="${out_dir}/${nsub}/mni_space/rois"
    orig_tracking_dir="${out_dir}/${nsub}/orig_space/tracking_second_order"
    mni_tracking_dir_second_order="${out_dir}/${nsub}/mni_space/tracking_second_order"
    first_order_final_dir="${out_dir}/${nsub}/mni_space/tracking_first_order/final_merged/final"

    orig_trials_root="${orig_tracking_dir}/trials"
    orig_merged_root="${orig_tracking_dir}/merged"

    echo ""
    echo "|------------- PROCESSING SECOND-ORDER ENSEMBLE FOR ${nsub} -------------|"
    echo ""

    # -------------------------
    # Safety checks
    # -------------------------
    for nside in left right; do
        [[ -f "${first_order_final_dir}/${nsub}_${nside}_spinal.trk" ]] || {
            echo "ERROR: Missing first-order file ${first_order_final_dir}/${nsub}_${nside}_spinal.trk"
            exit 1
        }
        [[ -f "${first_order_final_dir}/${nsub}_${nside}_remaining_cp.trk" ]] || {
            echo "ERROR: Missing first-order file ${first_order_final_dir}/${nsub}_${nside}_remaining_cp.trk"
            exit 1
        }
    done

    [[ -f "${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz" ]] || {
        echo "ERROR: Missing ${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz"
        exit 1
    }

    [[ -f "${mni_rois_dir}/${nsub}_aparc.DKTatlas+aseg_mni.nii.gz" ]] || {
        echo "ERROR: Missing ${mni_rois_dir}/${nsub}_aparc.DKTatlas+aseg_mni.nii.gz"
        exit 1
    }

    [[ -f "${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz" ]] || {
        echo "ERROR: Missing ${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz"
        echo "Check that -f matches the first-order FA threshold."
        exit 1
    }

    [[ -f "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" ]] || {
        echo "ERROR: Missing affine transform"
        exit 1
    }

    [[ -f "${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz" ]] || {
        echo "ERROR: Missing forward warp"
        exit 1
    }

    [[ -f "${out_dir}/${nsub}/orig_space/transfo/2orig_1InverseWarp.nii.gz" ]] || {
        echo "ERROR: Missing inverse warp"
        exit 1
    }

    # -------------------------
    # 1) Prepare second-order seed masks and thalamus ROIs
    # -------------------------
    echo "|------------- 1) Prepare second-order seed masks and thalamus ROIs -------------|"

    echo "|------------- 1.1) Copy pooled spinal seed masks and transform them to orig space -------------|"
    for nside in left right; do
        cp "${out_dir}/mni_space/tracking_first_order/final/all_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
           "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz"

        antsApplyTransforms \
            -d 3 \
            -i "${out_dir}/mni_space/tracking_first_order/final/all_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
            -r "${nsub_path}/tractoflow/${nsub}__t1_warped.nii.gz" \
            -t "${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz" \
            -t "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" \
            -o "${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz"
    done

    echo "|------------- 1.2) Extract thalamus masks from the anatomical segmentation -------------|"
    Right_Thalamus=(49)
    Left_Thalamus=(10)

    scil_labels_combine "${orig_rois_dir}/${nsub}_right_thalamus_orig.nii.gz" \
        --volume_ids "${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz" ${Right_Thalamus[*]} \
        --merge_groups -f

    scil_labels_combine "${mni_rois_dir}/${nsub}_right_thalamus_mni.nii.gz" \
        --volume_ids "${mni_rois_dir}/${nsub}_aparc.DKTatlas+aseg_mni.nii.gz" ${Right_Thalamus[*]} \
        --merge_groups -f

    scil_labels_combine "${orig_rois_dir}/${nsub}_left_thalamus_orig.nii.gz" \
        --volume_ids "${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz" ${Left_Thalamus[*]} \
        --merge_groups -f

    scil_labels_combine "${mni_rois_dir}/${nsub}_left_thalamus_mni.nii.gz" \
        --volume_ids "${mni_rois_dir}/${nsub}_aparc.DKTatlas+aseg_mni.nii.gz" ${Left_Thalamus[*]} \
        --merge_groups -f

    # -------------------------
    # 2) Run second-order ensemble tracking in orig space
    # -------------------------
    echo "|------------- 2) Run second-order ensemble tracking in orig space -------------|"
    for step_size in "${step_list[@]}"; do
        for theta in "${theta_list[@]}"; do
            combo_tag="step_${step_size}_theta_${theta}"
            mkdir -p "${orig_trials_root}/${combo_tag}"

            echo "|=== Second-order combo: ${combo_tag} ===|"

            for nside in left right; do
                # Tracking from Spinal bundle - long
                scil_tracking_local \
                    "${nsub_path}/tractoflow/${nsub}__fodf.nii.gz" \
                    "${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz" \
                    "${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz" \
                    "${orig_trials_root}/${combo_tag}/${nsub}_${nside}_from_spinal_track_npv1000_${combo_tag}.trk" \
                    --npv "${npv_spinal_long_per_combo}" \
                    --step "${step_size}" \
                    --theta "${theta}" \
                    ${gpu} -v -f

                # Tracking from Spinal bundle - short
                scil_tracking_local \
                    "${nsub_path}/tractoflow/${nsub}__fodf.nii.gz" \
                    "${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz" \
                    "${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz" \
                    "${orig_trials_root}/${combo_tag}/${nsub}_${nside}_from_spinal_track_npv100_${combo_tag}.trk" \
                    --npv "${npv_spinal_short_per_combo}" \
                    --step "${step_size}" \
                    --theta "${theta}" \
                    ${gpu} -v -f

                # Tracking from Thalamus
                scil_tracking_local \
                    "${nsub_path}/tractoflow/${nsub}__fodf.nii.gz" \
                    "${orig_rois_dir}/${nsub}_${nside}_thalamus_orig.nii.gz" \
                    "${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz" \
                    "${orig_trials_root}/${combo_tag}/${nsub}_${nside}_from_thalamus_npv500_${combo_tag}.trk" \
                    --npv "${npv_thalamus_per_combo}" \
                    --step "${step_size}" \
                    --theta "${theta}" \
                    ${gpu} -v -f
            done
        done
    done

    # -------------------------
    # 3) Merge second-order combo tractograms in orig space
    # -------------------------
    echo "|------------- 3) Merge second-order combo tractograms in orig space -------------|"
    for nside in left right; do
        for role in from_thalamus_npv500 from_spinal_track_npv100 from_spinal_track_npv1000; do
            files=()
            for step_size in "${step_list[@]}"; do
                for theta in "${theta_list[@]}"; do
                    combo_tag="step_${step_size}_theta_${theta}"
                    f="${orig_trials_root}/${combo_tag}/${nsub}_${nside}_${role}_${combo_tag}.trk"
                    [[ -f "${f}" ]] && files+=("${f}")
                done
            done

            if (( ${#files[@]} )); then
                scil_tractogram_math concatenate "${files[@]}" \
                    "${orig_merged_root}/${nsub}_${nside}_${role}.trk" -f
            else
                echo "WARN: no files found for ${nside} ${role}"
            fi
        done
    done

    # -------------------------
    # 4) Register merged second-order tractograms to MNI space
    # -------------------------
    echo "|------------- 4) Register merged second-order tractograms to MNI space -------------|"
    for nside in left right; do
        for role in from_thalamus_npv500 from_spinal_track_npv100 from_spinal_track_npv1000; do
            in_trk="${orig_merged_root}/${nsub}_${nside}_${role}.trk"
            out_trk="${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_${role}.trk"

            if [[ -f "${in_trk}" ]]; then
                scil_tractogram_apply_transform \
                    "${in_trk}" \
                    "${mni_dir}/MNI/mni_masked.nii.gz" \
                    "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" \
                    "${out_trk}" \
                    --in_deformation "${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz" \
                    --remove_invalid \
                    --reverse_operation -f
            else
                echo "WARN: merged orig tractogram not found for ${nside} ${role}"
            fi
        done
    done

    # -------------------------
    # 5) Prepare second-order MNI ROIs and cut masks
    # -------------------------
    echo "|------------- 5) Prepare second-order MNI ROIs and cut masks -------------|"
    echo "|------------- 5.1) Copy VPM and pathway-specific MNI ROIs -------------|"

    for nside in left right; do
        cp "${mni_dir}/MNI/Distal/${nside}/VPM.nii.gz" \
           "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz"

        for nroi in "${mni_dir}/MNI/from_${nside}"/*.nii.gz; do
            [[ -f "${nroi}" ]] || continue
            ROI_basename=$(basename "${nroi}")
            cp "${nroi}" "${mni_rois_dir}/${nsub}_second_order_${ROI_basename/nii/_mni.nii}"
        done
    done

    echo "|------------- 5.2) Build cut masks for second-order pathway extraction -------------|"
    for nside in left right; do
        if [[ "${nside}" == "left" ]]; then
            contra_nside="right"
        else
            contra_nside="left"
        fi

        # DTTT ipsilateral dPSN : remaining_cp to VPM
        scil_volume_math union \
            "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz" \
            "${out_dir}/mni_space/tracking_first_order/final/all_${nside}_remaining_cp_density_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_mni.nii.gz" \
            --data_type uint8 -f
        scil_labels_from_mask \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_labels_mni.nii.gz" \
            -f

        # DTTT ipsilateral CS : spinal to VPM
        scil_volume_math union \
            "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_mni.nii.gz" \
            --data_type uint8 -f
        scil_labels_from_mask \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_labels_mni.nii.gz" \
            -f

        # DTTT contralateral CS : spinal to thalamus
        scil_volume_math union \
            "${mni_rois_dir}/${nsub}_${contra_nside}_thalamus_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_mni.nii.gz" \
            --data_type uint8 -f
        scil_labels_from_mask \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_labels_mni.nii.gz" \
            -f

        # VTTT contralateral OS and IS : spinal to thalamus
        scil_volume_math union \
            "${mni_rois_dir}/${nsub}_${contra_nside}_thalamus_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_mni.nii.gz" \
            --data_type uint8 -f
        scil_labels_from_mask \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_labels_mni.nii.gz" \
            -f

        # VTTT contralateral vPSN : remaining_cp to VPM
        scil_volume_math union \
            "${mni_rois_dir}/${nsub}_${contra_nside}_VPM_mni.nii.gz" \
            "${out_dir}/mni_space/tracking_first_order/final/all_${nside}_remaining_cp_density_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_mni.nii.gz" \
            --data_type uint8 -f
        scil_labels_from_mask \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_labels_mni.nii.gz" \
            -f
    done

    # -------------------------
    # 6) Filter merged second-order tractograms into pathway bundles
    # -------------------------
    echo "|------------- 6) Filter merged second-order tractograms into pathway bundles -------------|"
    for nside in left right; do
        if [[ "${nside}" == "left" ]]; then
            contra_nside="right"
        else
            contra_nside="left"
        fi

        echo "|------------- 6.1) From ${nside} - VTTT (contralateral) - OS/IS and vPSN -------------|"

        # OS and IS
        scil_tractogram_filter_by_roi \
            "${mni_tracking_dir_second_order}/orig/${nsub}_${contra_nside}_from_thalamus_npv500.trk" \
            "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_VTTT_Controlat_OSandIS.trk" \
            --drawn_roi "${mni_rois_dir}/${nsub}_left_cerebellum_wm_mni.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_rois_dir}/${nsub}_right_cerebellum_wm_mni.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" 'any' 'include' \
            --drawn_roi "${mni_rois_dir}/${nsub}_${contra_nside}_thalamus_mni.nii.gz" 'any' 'include' \
            --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_INC_Pons_Controlat.nii.gz" 'any' 'include' \
            --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_CaudalMedulla_Controlat.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_INC_VTT_Area.nii.gz" 'any' 'include' \
            --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Pons_Ipsilat.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_dir}/MNI/cs_plaque.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_INC_VTT_Area.nii.gz" 'any' 'include' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/VTTT_Controlat_OSandIS_1.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/VTTT_Controlat_OSandIS_2.bdo" 'any' 'exclude' \
            -f

        # vPSN
        scil_tractogram_filter_by_roi \
            "${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_spinal_track_npv1000.trk" \
            "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_VTTT_Controlat_vPSN.trk" \
            --drawn_roi "${mni_rois_dir}/${nsub}_left_cerebellum_wm_mni.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_rois_dir}/${nsub}_right_cerebellum_wm_mni.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" 'either_end' 'include' \
            --drawn_roi "${mni_rois_dir}/${nsub}_${contra_nside}_VPM_mni.nii.gz" 'any' 'include' \
            --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_CaudalMedulla_Controlat.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_INC_VTT_Area.nii.gz" 'any' 'include' \
            --bdo "${mni_dir}/MNI/from_${nside}/VTTT_Controlat/VTTT_Controlat_vPSN_1.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/VTTT_Controlat/VTTT_Controlat_vPSN_2.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/VTTT_Controlat/VTTT_Controlat_vPSN_3.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/VTTT_Controlat/VTTT_Controlat_vPSN_4.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/VTTT_Controlat/VTTT_Controlat_vPSN_5.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/VTTT_Controlat/VTTT_Controlat_vPSN_6.bdo" 'any' 'exclude' \
            -f

        echo "|------------- 6.2) ${nside} - DTTT (contralateral) - CS -------------|"
        scil_tractogram_filter_by_roi \
            "${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_thalamus_npv500.trk" \
            "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${contra_nside}_DTTT_Controlat_CS.trk" \
            --drawn_roi "${mni_rois_dir}/${nsub}_${contra_nside}_spinal_density_second_order_seed_mni.nii.gz" 'any' 'include' \
            --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_thalamus_mni.nii.gz" 'any' 'include' \
            --drawn_roi "${mni_dir}/MNI/from_${contra_nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_dir}/MNI/from_${contra_nside}/DTTT_Controlat_INC_CaudalMedulla_Ipsilat.nii.gz" 'any' 'include' \
            --drawn_roi "${mni_dir}/MNI/from_${contra_nside}/DTTT_Controlat_INC_Medulla_Controlat.nii.gz" 'any' 'include' \
            --drawn_roi "${mni_dir}/MNI/from_${contra_nside}/DTTT_Controlat_EXC_Midbrain_Ipsilat.nii.gz" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Controlat_1.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Controlat_2.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Controlat_3.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Controlat_4.bdo" 'any' 'exclude' \
            -f

        echo "|------------- 6.3) ${nside} - DTTT (ipsilateral) - dPSN and CS -------------|"

        # CS
        scil_tractogram_filter_by_roi \
            "${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_spinal_track_npv100.trk" \
            "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk" \
            --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" 'either_end' 'include' \
            --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz" 'any' 'include' \
            --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_dir}/MNI/midsagittal_plane.nii.gz" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_CS_1.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_CS_2.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_CS_3.bdo" 'any' 'exclude' \
            -f

        # dPSN
        scil_tractogram_filter_by_roi \
            "${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_spinal_track_npv1000.trk" \
            "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_dPSN.trk" \
            --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz" 'any' 'include' \
            --drawn_roi "${out_dir}/mni_space/tracking_first_order/final/all_${nside}_remaining_cp_density_mni.nii.gz" 'either_end' 'include' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_1.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_2.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_3.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_4.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_5.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_6.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_7.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_8.bdo" 'any' 'exclude' \
            -f
    done



# -------------------------
# 7) Cut filtered second-order bundles with label masks
# -------------------------
echo "|------------- 7) Cut filtered second-order bundles with label masks -------------|"
for nside in left right; do

    in_trk="${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_VTTT_Controlat_OSandIS.trk"
    out_trk="${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_VTTT_Controlat_OSandIS.trk"
    if trk_is_empty "${in_trk}"; then
        echo "WARN: ${in_trk} is missing or empty, skipping cut."
    else
        scil_tractogram_cut_streamlines \
            "${in_trk}" \
            --labels "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_labels_mni.nii.gz" \
            "${out_trk}" -f

        if trk_is_empty "${out_trk}"; then
            echo "WARN: ${out_trk} is empty after cut, removing."
            rm -f "${out_trk}"
        fi
    fi

    in_trk="${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_VTTT_Controlat_vPSN.trk"
    out_trk="${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_VTTT_Controlat_vPSN.trk"
    if trk_is_empty "${in_trk}"; then
        echo "WARN: ${in_trk} is missing or empty, skipping cut."
    else
        scil_tractogram_cut_streamlines \
            "${in_trk}" \
            --labels "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_labels_mni.nii.gz" \
            "${out_trk}" -f

        if trk_is_empty "${out_trk}"; then
            echo "WARN: ${out_trk} is empty after cut, removing."
            rm -f "${out_trk}"
        fi
    fi

    in_trk="${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Controlat_CS.trk"
    out_trk="${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Controlat_CS.trk"
    if trk_is_empty "${in_trk}"; then
        echo "WARN: ${in_trk} is missing or empty, skipping cut."
    else
        scil_tractogram_cut_streamlines \
            "${in_trk}" \
            --labels "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_labels_mni.nii.gz" \
            "${out_trk}" -f

        if trk_is_empty "${out_trk}"; then
            echo "WARN: ${out_trk} is empty after cut, removing."
            rm -f "${out_trk}"
        fi
    fi

    in_trk="${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk"
    out_trk="${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk"
    if trk_is_empty "${in_trk}"; then
        echo "WARN: ${in_trk} is missing or empty, skipping cut."
    else
        scil_tractogram_cut_streamlines \
            "${in_trk}" \
            --labels "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_labels_mni.nii.gz" \
            "${out_trk}" -f

        if trk_is_empty "${out_trk}"; then
            echo "WARN: ${out_trk} is empty after cut, removing."
            rm -f "${out_trk}"
        fi
    fi

    in_trk="${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_dPSN.trk"
    out_trk="${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Ipsilat_dPSN.trk"
    if trk_is_empty "${in_trk}"; then
        echo "WARN: ${in_trk} is missing or empty, skipping cut."
    else
        scil_tractogram_cut_streamlines \
            "${in_trk}" \
            --labels "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_labels_mni.nii.gz" \
            "${out_trk}" -f

        if trk_is_empty "${out_trk}"; then
            echo "WARN: ${out_trk} is empty after cut, removing."
            rm -f "${out_trk}"
        fi
    fi
done












    # -------------------------
    # 8) Reject outliers and save final second-order bundles
    # -------------------------
    echo "|------------- 8) Reject outliers and save final second-order bundles -------------|"
    for nside in left right; do
        scil_bundle_reject_outliers \
            "${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk" \
            "${mni_tracking_dir_second_order}/final/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk" \
            --alpha 0.30 -f

        for nbundle in DTTT_Ipsilat_dPSN DTTT_Controlat_CS VTTT_Controlat_OSandIS VTTT_Controlat_vPSN; do
            scil_bundle_reject_outliers \
                "${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_${nbundle}.trk" \
                "${mni_tracking_dir_second_order}/final/${nsub}_from_${nside}_${nbundle}.trk" \
                --alpha 0.50 -f
        done
    done

    echo "|------------- SECOND-ORDER ENSEMBLE FOR ${nsub} IS COMPLETED -------------|"
    echo ""
done


echo "|------------- 9) [BACK-TO-ORIG] Register final MNI bundles to orig space -------------|"

for nside in left right; do
    for nbundle in DTTT_Ipsilat_CS DTTT_Ipsilat_dPSN DTTT_Controlat_CS VTTT_Controlat_OSandIS VTTT_Controlat_vPSN; do
        in_trk="${mni_tracking_dir_second_order}/final/${nsub}_from_${nside}_${nbundle}.trk"
        out_trk="${orig_tracking_dir}/final/${nsub}_from_${nside}_${nbundle}_orig.trk"

        if [[ -f "${in_trk}" ]]; then
            scil_tractogram_apply_transform \
                "${in_trk}" \
                "${nsub_path}/tractoflow/${nsub}__t1_warped.nii.gz" \
                "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" \
                "${out_trk}" \
                --inverse \
                --in_deformation "${out_dir}/${nsub}/orig_space/transfo/2orig_1InverseWarp.nii.gz" \
                --remove_invalid -f
        else
            echo "WARN: Final MNI bundle not found for ${nside} ${nbundle}, skipping back-to-orig."
        fi
    done
done



#test for second_order_path