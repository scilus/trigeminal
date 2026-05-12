#!/bin/bash
set -euo pipefail

# TRIGEMINAL SYSTEM TRACTOGRAPHY - Samir Akeb (2022-2023)
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Arnaud Bore (2023-2024)
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Nasrin Rafiei (2025-2026)
#
# SECOND-ORDER ENSEMBLE VERSION
# SUBJECT-BY-SUBJECT VERSION:
#   - no pooled all_left/all_right files
#   - each subject uses its own first-order spinal / remaining_cp bundles
#   - nothing is created outside the subject folder
#   - subject-specific density maps are stored directly in mni_space/rois
#
# Organized pipeline:
#   0) prepare subject-specific first-order density maps
#   1) prepare second-order seed masks and thalamus ROIs
#   2) run second-order tracking for each (step, theta) combo in ORIG
#   3) merge combo outputs per tracking role in ORIG
#   4) register merged tractograms once to MNI
#   5) prepare second-order MNI ROIs and cut masks
#   6) filter merged second-order tractograms into pathway bundles
#   7) cut filtered bundles with label masks
#   8) reject outliers and save final second-order bundles
#   9) bring final MNI bundles back to ORIG

usage() {
    cat <<EOF 1>&2
Usage:
  $(basename "$0") -s <subjects_parent_or_single_subject_dir> -m <ROIs_clean_dir> -o <out_dir>
                   [-f fa_threshold] [-t threads] [-g true|false]
                   [-p step_size] [-e theta_deg]
                   [--npv_spinal_long N] [--npv_spinal_short N] [--npv_thalamus N]
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
if [[ -n "${g}" && "${g}" == "true" ]]; then
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

mask_has_nonzero_voxels() {
    local f="$1"

    if [[ ! -f "${f}" ]]; then
        return 1
    fi

    local nz
    nz=$(python - <<PY
import nibabel as nib
import numpy as np
img = nib.load(r"""${f}""")
print(int((img.get_fdata() > 0).sum()))
PY
)

    [[ "${nz}" -gt 0 ]]
}

get_label_ids_for_cut() {
    local labels_file="$1"

    python - <<PY
import nibabel as nib
import numpy as np

f = r"""${labels_file}"""
vals = np.unique(nib.load(f).get_fdata())
vals = [int(v) for v in vals if v != 0]

if len(vals) < 2:
    print("")
elif len(vals) == 2:
    print(f"{vals[0]} {vals[1]}")
else:
    print(f"{vals[0]} {vals[-1]}")
PY
}

safe_compute_density_map() {
    local in_trk="$1"
    local out_mask="$2"

    if trk_is_empty "${in_trk}"; then
        echo "WARN: ${in_trk} is missing or empty. Skipping density map."
        rm -f "${out_mask}"
        return 0
    fi

    scil_tractogram_compute_density_map \
        "${in_trk}" \
        "${out_mask}" \
        --binary -f

    if ! mask_has_nonzero_voxels "${out_mask}"; then
        echo "WARN: ${out_mask} has no non-zero voxels after density computation. Removing."
        rm -f "${out_mask}"
    fi
}

safe_apply_mask_transform_to_orig() {
    local in_mask="$1"
    local ref_img="$2"
    local warp="$3"
    local affine="$4"
    local out_mask="$5"

    if ! mask_has_nonzero_voxels "${in_mask}"; then
        echo "WARN: ${in_mask} is missing or empty. Skipping mask transform."
        rm -f "${out_mask}"
        return 0
    fi

    antsApplyTransforms \
        -d 3 \
        -i "${in_mask}" \
        -r "${ref_img}" \
        -t "${warp}" \
        -t "${affine}" \
        -o "${out_mask}"

    if ! mask_has_nonzero_voxels "${out_mask}"; then
        echo "WARN: ${out_mask} has no non-zero voxels after transform. Removing."
        rm -f "${out_mask}"
    fi
}

safe_build_cut_labels() {
    local mask_a="$1"
    local mask_b="$2"
    local out_mask="$3"
    local out_labels="$4"

    if ! mask_has_nonzero_voxels "${mask_a}"; then
        echo "WARN: ${mask_a} is missing or empty. Skipping cut-mask construction."
        rm -f "${out_mask}" "${out_labels}"
        return 0
    fi

    if ! mask_has_nonzero_voxels "${mask_b}"; then
        echo "WARN: ${mask_b} is missing or empty. Skipping cut-mask construction."
        rm -f "${out_mask}" "${out_labels}"
        return 0
    fi

    scil_volume_math union \
        "${mask_a}" \
        "${mask_b}" \
        "${out_mask}" \
        --data_type uint8 -f

    if ! mask_has_nonzero_voxels "${out_mask}"; then
        echo "WARN: ${out_mask} has no non-zero voxels after union. Removing."
        rm -f "${out_mask}" "${out_labels}"
        return 0
    fi

    scil_labels_from_mask \
        "${out_mask}" \
        "${out_labels}" \
        -f
}

safe_filter_by_roi() {
    local in_trk="$1"
    local out_trk="$2"
    shift 2

    if trk_is_empty "${in_trk}"; then
        echo "WARN: ${in_trk} is missing or empty. Skipping filter."
        rm -f "${out_trk}"
        return 0
    fi

    local args=("$@")
    local i=0
    local n=${#args[@]}

    while (( i < n )); do
        local key="${args[$i]}"

        if [[ "${key}" == "--drawn_roi" ]]; then
            if (( i + 3 >= n )); then
                echo "WARN: malformed --drawn_roi arguments for ${out_trk}. Skipping filter."
                rm -f "${out_trk}"
                return 0
            fi

            local roi_file="${args[$((i + 1))]}"
            local roi_action="${args[$((i + 3))]}"

            if [[ ! -f "${roi_file}" ]]; then
                echo "WARN: ROI ${roi_file} is missing. Skipping filter for ${out_trk}."
                rm -f "${out_trk}"
                return 0
            fi

            if [[ "${roi_action}" == "include" ]] && ! mask_has_nonzero_voxels "${roi_file}"; then
                echo "WARN: include ROI ${roi_file} is empty. Skipping filter for ${out_trk}."
                rm -f "${out_trk}"
                return 0
            fi

            i=$((i + 4))
            continue
        fi

        if [[ "${key}" == "--bdo" ]]; then
            if (( i + 3 >= n )); then
                echo "WARN: malformed --bdo arguments for ${out_trk}. Skipping filter."
                rm -f "${out_trk}"
                return 0
            fi

            local bdo_file="${args[$((i + 1))]}"

            if [[ ! -f "${bdo_file}" ]]; then
                echo "WARN: BDO ${bdo_file} is missing. Skipping filter for ${out_trk}."
                rm -f "${out_trk}"
                return 0
            fi

            i=$((i + 4))
            continue
        fi

        i=$((i + 1))
    done

    scil_tractogram_filter_by_roi \
        "${in_trk}" \
        "${out_trk}" \
        "${args[@]}" \
        -f

    if trk_is_empty "${out_trk}"; then
        echo "WARN: ${out_trk} is empty after filtering. Removing."
        rm -f "${out_trk}"
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

    for nside in left right; do
        if trk_is_empty "${first_order_final_dir}/${nsub}_${nside}_spinal.trk"; then
            echo "WARN: Missing or empty first-order spinal file: ${first_order_final_dir}/${nsub}_${nside}_spinal.trk"
            echo "WARN: ${nside} spinal-based second-order tracking will be skipped."
        fi

        if trk_is_empty "${first_order_final_dir}/${nsub}_${nside}_remaining_cp.trk"; then
            echo "WARN: Missing or empty first-order remaining_cp file: ${first_order_final_dir}/${nsub}_${nside}_remaining_cp.trk"
            echo "WARN: ${nside} remaining_cp-dependent second-order bundles may be skipped."
        fi
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
    # 0) Prepare subject-specific first-order density maps
    # -------------------------
    echo "|------------- 0) Prepare subject-specific first-order density maps -------------|"
    for nside in left right; do
        safe_compute_density_map \
            "${first_order_final_dir}/${nsub}_${nside}_spinal.trk" \
            "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz"

        safe_compute_density_map \
            "${first_order_final_dir}/${nsub}_${nside}_remaining_cp.trk" \
            "${mni_rois_dir}/${nsub}_${nside}_remaining_cp_density_mni.nii.gz"
    done

    # -------------------------
    # 1) Prepare second-order seed masks and thalamus ROIs
    # -------------------------
    echo "|------------- 1) Prepare second-order seed masks and thalamus ROIs -------------|"

    echo "|------------- 1.1) Transform subject spinal seed masks to orig space -------------|"
    for nside in left right; do
        safe_apply_mask_transform_to_orig \
            "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
            "${nsub_path}/tractoflow/${nsub}__t1_warped.nii.gz" \
            "${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz" \
            "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" \
            "${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz"
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
                spinal_seed="${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz"
                thalamus_seed="${orig_rois_dir}/${nsub}_${nside}_thalamus_orig.nii.gz"
                wm_mask="${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz"

                if mask_has_nonzero_voxels "${spinal_seed}"; then
                    scil_tracking_local \
                        "${nsub_path}/tractoflow/${nsub}__fodf.nii.gz" \
                        "${spinal_seed}" \
                        "${wm_mask}" \
                        "${orig_trials_root}/${combo_tag}/${nsub}_${nside}_from_spinal_track_npv1000_${combo_tag}.trk" \
                        --npv "${npv_spinal_long_per_combo}" \
                        --step "${step_size}" \
                        --theta "${theta}" \
                        ${gpu} -v -f

                    scil_tracking_local \
                        "${nsub_path}/tractoflow/${nsub}__fodf.nii.gz" \
                        "${spinal_seed}" \
                        "${wm_mask}" \
                        "${orig_trials_root}/${combo_tag}/${nsub}_${nside}_from_spinal_track_npv100_${combo_tag}.trk" \
                        --npv "${npv_spinal_short_per_combo}" \
                        --step "${step_size}" \
                        --theta "${theta}" \
                        ${gpu} -v -f
                else
                    echo "WARN: spinal seed is missing or empty for ${nside} at ${combo_tag}. Skipping spinal tracking."
                fi

                if mask_has_nonzero_voxels "${thalamus_seed}"; then
                    scil_tracking_local \
                        "${nsub_path}/tractoflow/${nsub}__fodf.nii.gz" \
                        "${thalamus_seed}" \
                        "${wm_mask}" \
                        "${orig_trials_root}/${combo_tag}/${nsub}_${nside}_from_thalamus_npv500_${combo_tag}.trk" \
                        --npv "${npv_thalamus_per_combo}" \
                        --step "${step_size}" \
                        --theta "${theta}" \
                        ${gpu} -v -f
                else
                    echo "WARN: thalamus seed is missing or empty for ${nside} at ${combo_tag}. Skipping thalamus tracking."
                fi
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
            cleaned_files=()
            clean_dir="${orig_tracking_dir}/trials_cleaned/${nside}_${role}"
            mkdir -p "${clean_dir}"

            for step_size in "${step_list[@]}"; do
                for theta in "${theta_list[@]}"; do
                    combo_tag="step_${step_size}_theta_${theta}"
                    f="${orig_trials_root}/${combo_tag}/${nsub}_${nside}_${role}_${combo_tag}.trk"
                    if ! trk_is_empty "${f}"; then
                        files+=("${f}")
                    fi
                done
            done

            if (( ${#files[@]} == 0 )); then
                echo "WARN: no non-empty files found for ${nside} ${role}"
                rm -f "${orig_merged_root}/${nsub}_${nside}_${role}.trk"
                continue
            fi

            echo "Cleaning and merging ${#files[@]} files for ${nside} ${role}"

            for f in "${files[@]}"; do
                cleaned_f="${clean_dir}/$(basename "${f}")"

python <<PY
from dipy.io.streamline import load_tractogram, save_tractogram
sft = load_tractogram(r"""${f}""", 'same', bbox_valid_check=False)
sft.remove_invalid_streamlines()
save_tractogram(sft, r"""${cleaned_f}""", bbox_valid_check=False)
PY

                cleaned_files+=("${cleaned_f}")
            done

            scil_tractogram_math concatenate "${cleaned_files[@]}" \
                "${orig_merged_root}/${nsub}_${nside}_${role}.trk" -f
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

            if trk_is_empty "${in_trk}"; then
                echo "WARN: merged orig tractogram missing or empty for ${nside} ${role}"
                rm -f "${out_trk}"
            else
                scil_tractogram_apply_transform \
                    "${in_trk}" \
                    "${mni_dir}/MNI/mni_masked.nii.gz" \
                    "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" \
                    "${out_trk}" \
                    --in_deformation "${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz" \
                    --remove_invalid \
                    --reverse_operation -f

                if trk_is_empty "${out_trk}"; then
                    echo "WARN: ${out_trk} is empty after transform to MNI. Removing."
                    rm -f "${out_trk}"
                fi
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

        safe_build_cut_labels \
            "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_remaining_cp_density_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_labels_mni.nii.gz"

        safe_build_cut_labels \
            "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_labels_mni.nii.gz"

        safe_build_cut_labels \
            "${mni_rois_dir}/${nsub}_${contra_nside}_thalamus_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_labels_mni.nii.gz"

        safe_build_cut_labels \
            "${mni_rois_dir}/${nsub}_${contra_nside}_thalamus_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_labels_mni.nii.gz"

        safe_build_cut_labels \
            "${mni_rois_dir}/${nsub}_${contra_nside}_VPM_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_remaining_cp_density_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_mni.nii.gz" \
            "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_labels_mni.nii.gz"
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

        safe_filter_by_roi \
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
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/VTTT_Controlat_OSandIS_2.bdo" 'any' 'exclude'

        safe_filter_by_roi \
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
            --bdo "${mni_dir}/MNI/from_${nside}/VTTT_Controlat/VTTT_Controlat_vPSN_6.bdo" 'any' 'exclude'

        echo "|------------- 6.2) ${nside} - DTTT (contralateral) - CS -------------|"
        safe_filter_by_roi \
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
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Controlat_4.bdo" 'any' 'exclude'

        echo "|------------- 6.3) ${nside} - DTTT (ipsilateral) - dPSN and CS -------------|"

        safe_filter_by_roi \
            "${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_spinal_track_npv100.trk" \
            "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk" \
            --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" 'either_end' 'include' \
            --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz" 'any' 'include' \
            --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz" 'any' 'exclude' \
            --drawn_roi "${mni_dir}/MNI/midsagittal_plane.nii.gz" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_CS_1.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_CS_2.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_CS_3.bdo" 'any' 'exclude'

        safe_filter_by_roi \
            "${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_spinal_track_npv1000.trk" \
            "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_dPSN.trk" \
            --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz" 'any' 'include' \
            --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_remaining_cp_density_mni.nii.gz" 'either_end' 'include' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_1.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_2.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_3.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_4.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_5.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_6.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_7.bdo" 'any' 'exclude' \
            --bdo "${mni_dir}/MNI/from_${nside}/new_ROIs/DTTT_Ipsilat_dPSN_8.bdo" 'any' 'exclude'
    done

    # -------------------------
    # 7) Cut filtered second-order bundles with label masks
    # -------------------------
    echo "|------------- 7) Cut filtered second-order bundles with label masks -------------|"
    for nside in left right; do
        for nbundle in \
            VTTT_Controlat_OSandIS \
            VTTT_Controlat_vPSN \
            DTTT_Controlat_CS \
            DTTT_Ipsilat_CS \
            DTTT_Ipsilat_dPSN; do

            in_trk="${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_${nbundle}.trk"
            out_trk="${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_${nbundle}.trk"
            labels_file="${mni_rois_dir}/${nsub}_${nside}_second_order_${nbundle}_Cuts_labels_mni.nii.gz"

            if trk_is_empty "${in_trk}"; then
                echo "WARN: ${in_trk} is missing or empty, skipping cut."
                continue
            fi

            if [[ ! -f "${labels_file}" ]]; then
                echo "WARN: Missing labels file ${labels_file}, skipping."
                continue
            fi

            label_ids=$(get_label_ids_for_cut "${labels_file}")

            if [[ -z "${label_ids}" ]]; then
                echo "WARN: Could not determine 2 valid label ids for ${labels_file}, skipping."
                continue
            fi

            echo "Cutting ${nside} ${nbundle} with label ids: ${label_ids}"

            scil_tractogram_cut_streamlines \
                "${in_trk}" \
                --labels "${labels_file}" \
                --label_ids ${label_ids} \
                "${out_trk}" -f

            if trk_is_empty "${out_trk}"; then
                echo "WARN: ${out_trk} is empty after cut, removing."
                rm -f "${out_trk}"
            fi
        done
    done

    # -------------------------
    # 8) Reject outliers and save final second-order bundles
    # -------------------------
    echo "|------------- 8) Reject outliers and save final second-order bundles -------------|"
    for nside in left right; do
        in_trk="${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk"
        out_trk="${mni_tracking_dir_second_order}/final/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk"

        if trk_is_empty "${in_trk}"; then
            echo "WARN: ${in_trk} missing or empty, skipping outlier rejection."
        else
            scil_bundle_reject_outliers \
                "${in_trk}" \
                "${out_trk}" \
                --alpha 0.30 -f

            if trk_is_empty "${out_trk}"; then
                echo "WARN: ${out_trk} is empty after outlier rejection, removing."
                rm -f "${out_trk}"
            fi
        fi

        for nbundle in DTTT_Ipsilat_dPSN DTTT_Controlat_CS VTTT_Controlat_OSandIS VTTT_Controlat_vPSN; do
            in_trk="${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_${nbundle}.trk"
            out_trk="${mni_tracking_dir_second_order}/final/${nsub}_from_${nside}_${nbundle}.trk"

            if trk_is_empty "${in_trk}"; then
                echo "WARN: ${in_trk} missing or empty, skipping outlier rejection."
            else
                scil_bundle_reject_outliers \
                    "${in_trk}" \
                    "${out_trk}" \
                    --alpha 0.50 -f

                if trk_is_empty "${out_trk}"; then
                    echo "WARN: ${out_trk} is empty after outlier rejection, removing."
                    rm -f "${out_trk}"
                fi
            fi
        done
    done

    # -------------------------
    # 9) [BACK-TO-ORIG] Register final MNI bundles to orig space
    # -------------------------
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

    echo "|------------- SECOND-ORDER ENSEMBLE FOR ${nsub} IS COMPLETED -------------|"
    echo ""
done

# git checkout -b second_order_incremental_concatenate
