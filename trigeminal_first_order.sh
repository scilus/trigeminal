#!/bin/bash
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Samir Akeb (2022-2023)
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Arnaud Bore (2023-2024)



# This script performs ensemble tractography of the trigeminal system using fODFs and FA maps
# from TractoFlow. Anatomical ROIs from both FreeSurfer and the ROIs_mean MNI folder are used to
# segment and virtually dissect the trigeminal nerve. Tracking is performed in original space over a
# grid of (step size, theta) combinations, merged, transformed to MNI space, filtered and segmented 
# into three components (mesencephalic, spinal, remaining nucleus). Final cleaned bundles are saved
# in both original (orig_space) and MNI (mni_space) coordinate spaces.


# EXAMPLE COMMAND
#
#   bash trigeminal_first_order.sh \
#       -s /path/to/subject/folder/sub-01 \
#       -m /path/to/ROIs_mean \
#       -o /path/to/output_dir \
#       -t 10 \
#       -p 0.5 \
#       -e 30 \
#       -f 0.15 \
#       -n 20000 \
#       -g true
#
# Notes:
#   - If -p (step) and -e (theta) are not provided, the script performs ensemble tracking using:
#       step = {0.1, 0.5, 1.0}
#       theta = {20, 30, 40}
#   - GPU option (-g) accelerates local tracking when supported; omit it to run on CPU.



# Input structure
#
#    [input]
#    ├── sub-01
#    │   ├── freesurfer
#    │   │   └─── aparc.DKTatlas+aseg.mgz
#    │   ├── sub-01__fa.nii.gz
#    │   ├── sub-01__fodf.nii.gz
#    │   └── sub-01__t1_warped.nii.gz
#    │
#    ├── S2
#    .
#    .


# -m : Path to the MNI-space reference folder.
#       In this project, it should point to the path/to/trigeminal/ROIs_clean/ folder
#       within the trigeminal GitHub project.

usage() { 
    echo "$(basename $0) [-s path/to/subject] [-m path/to/mni] [-o output_dir] [-t nb_threads] [-p step_size] [-e theta_deg] [-f fa_threshold] [-n npv_first_order] -g true" 1>&2
    exit 1
}


while getopts "s:m:o:t:p:e:f:n:g:" args; do
    case "${args}" in
        s) subject_dir=${OPTARG} ;;
        m) mni_dir=${OPTARG} ;;
        o) output_dir=${OPTARG} ;;
        t) nb_threads=${OPTARG} ;;
        p) step_size=${OPTARG} ;;
        e) theta=${OPTARG} ;;
        f) fa_threshold=${OPTARG} ;;   
        n) npv_first_order=${OPTARG} ;;    
        g) gpu=${OPTARG} ;;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${subject_dir}" ] || [ -z "${mni_dir}" ] || [ -z "${output_dir}" ]; then
    usage
fi




# If the user provides both step size (-p) and theta (-e) values, the script will use those specific parameters and perform a single run
# of local_tracking (i.e., no ensemble tractography will be performed).

if [ -n "${step_size}" ] && [ -n "${theta}" ]; then
    step_list=(${step_size})
    theta_list=(${theta})
else
    step_list=(0.1 0.5 1.0)
    theta_list=(20 30 40)
fi

# Default values if not set
fa_threshold=${fa_threshold:-0.15}
npv_first_order=${npv_first_order:-20000}


# npv_first_order is the total number of seeds per voxel for the whole first-order tracking.
# It is divided by the number of step/theta combinations.
# The resulting npv_per_run is the number of seeds actually used in each run.

# number of step/theta combos
n_combos=$(( ${#step_list[@]} * ${#theta_list[@]} ))
# seeds per combo (rounded)
npv_per_run=$(( (npv_first_order + n_combos - 1) / n_combos ))  # ceiling division
echo "Using $npv_per_run seeds per voxel per run (based on $npv_first_order total)"



#opposite_side=leftright

# Enable GPU if available (highly recommended for X times faster processing)

if [ ! -z "${gpu}" ]; then
    gpu="--use_gpu"
else
    gpu=""
fi

echo "Folder subject: " ${subject_dir}
echo "Folder MNI: " ${mni_dir}
echo "Output folder: " ${output_dir}
echo "Use GPU: " ${gpu}
echo "Number of threads" ${nb_threads}
echo "Tracking grid: steps=${step_list[*]}  thetas=${theta_list[*]}"


nsub=$(basename "${subject_dir}")
subject_parent=$(dirname "${subject_dir}")

    
mkdir -p ${output_dir}/${nsub}/orig_space/{rois,tracking_first_order,transfo}
mkdir -p ${output_dir}/${nsub}/mni_space/{rois,tracking_first_order}

# Per-combo dirs will be created later, inside the loops.
orig_rois_dir=${output_dir}/${nsub}/orig_space/rois
mni_rois_dir=${output_dir}/${nsub}/mni_space/rois
orig_tracking_dir=${output_dir}/${nsub}/orig_space/tracking_first_order
mni_tracking_root=${output_dir}/${nsub}/mni_space/tracking_first_order
trials_dir=${orig_tracking_dir}/trials
mkdir -p "${trials_dir}"

echo ""
echo "|------------- PROCESSING FIRST ORDER TGN TRACTOGRAPHY FOR ${nsub} -------------|"
echo ""

echo "|------------- 1) Registrations -------------|"
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="${nb_threads}"

echo "Running ANTs registration... (logs saved to ${output_dir}/${nsub}/log.txt)"
antsRegistrationSyN.sh \
    -d 3 \
    -f "${subject_dir}/tractoflow/${nsub}__t1_warped.nii.gz" \
    -m "${mni_dir}/MNI/mni_masked.nii.gz" \
    -t s \
    -o "${output_dir}/${nsub}/orig_space/transfo/2orig_"
    > "${output_dir}/${nsub}/log.txt" 2>&1

## [ORIG-SPACE] Register all ROIs
for nroi in cp_left cp_right; do
    antsApplyTransforms \
    -d 3 \
    -i ${mni_dir}/MNI/${nroi}.nii.gz \
    -r ${subject_dir}/tractoflow/${nsub}__t1_warped.nii.gz \
    -t ${output_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz \
    -t ${output_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat \
    -o ${orig_rois_dir}/${nsub}_${nroi}_orig.nii.gz

    scil_volume_math lower_threshold_eq ${orig_rois_dir}/${nsub}_${nroi}_orig.nii.gz 0.5 ${orig_rois_dir}/${nsub}_${nroi}_orig.nii.gz --data_type int16 -f
done

## [ORIG-SPACE] Reshape aparc.DKTatlas+aseg.mgz
scil_volume_reslice_to_reference ${subject_dir}/freesurfer/aparc.DKTatlas+aseg.mgz \
    ${subject_dir}/tractoflow/${nsub}__t1_warped.nii.gz \
    ${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz \
    --interpolation nearest --keep_dtype -f

echo "|------------- 1) Done -------------|"
echo ""

echo "|------------- 2) Generate exclusions and inclusions ROI -------------|"
## [ORIG-SPACE] Exclusions ROI labels
Left_Cerebral_Cortex=($(seq 1000 1035))
Right_Cerebral_Cortex=($(seq 2000 2035))
Left_Cerebral_WM=(2)
Right_Cerebral_WM=(41)
Cerebellum_Cortex=(8 47)
Right_Cerebellum_WM=(46)
Left_Cerebellum_WM=(7)
Any_Exclusion_ROI=(${Left_Cerebral_Cortex[*]} ${Right_Cerebral_Cortex[*]} ${Left_Cerebral_WM[*]} ${Right_Cerebral_WM[*]} ${Cerebellum_Cortex[*]})

echo "|------------- 2.1) any_exclusion_roi_orig -------------|"
scil_labels_combine ${orig_rois_dir}/${nsub}_any_exclusion_roi_orig.nii.gz \
    --volume_ids ${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz ${Any_Exclusion_ROI[*]} \
    --merge_groups -f

echo "|------------- 2.2) cerebellum_wm_right_orig -------------|"
scil_labels_combine ${orig_rois_dir}/${nsub}_right_cerebellum_wm_orig.nii.gz \
    --volume_ids ${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz ${Right_Cerebellum_WM[*]} \
    --merge_groups -f

echo "|------------- 2.3) cerebellum_wm_left_orig -------------|"
scil_labels_combine ${orig_rois_dir}/${nsub}_left_cerebellum_wm_orig.nii.gz \
    --volume_ids ${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz ${Left_Cerebellum_WM[*]} \
    --merge_groups -f

# WM mask
scil_volume_math lower_threshold ${subject_dir}/tractoflow/${nsub}__fa.nii.gz \
    ${fa_threshold} \
    ${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz \
    --data_type uint8 -f

echo "|------------- 2) Done -------------|"
echo ""

echo "|------------- 3) Register ROI in MNI space -------------|"
for nroi in any_exclusion_roi_orig right_cerebellum_wm_orig left_cerebellum_wm_orig wm_mask_${fa_threshold}_orig aparc.DKTatlas+aseg_orig; do
    antsApplyTransforms \
    -d 3 \
    -i ${orig_rois_dir}/${nsub}_${nroi}.nii.gz \
    -r ${mni_dir}/MNI/mni_masked.nii.gz \
    -t ${output_dir}/${nsub}/orig_space/transfo/2orig_1InverseWarp.nii.gz \
    -t [${output_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat, 1] \
    -o ${mni_rois_dir}/${nsub}_${nroi/orig/mni}.nii.gz \
    -n NearestNeighbor

    scil_volume_math convert ${mni_rois_dir}/${nsub}_${nroi/orig/mni}.nii.gz ${mni_rois_dir}/${nsub}_${nroi/orig/mni}.nii.gz --data_type int16 -f
done

echo "|------------- 3) Done -------------|"
echo ""

# =========================
# Grid over (step, theta)
# =========================
for step_size in "${step_list[@]}"; do
  for theta in "${theta_list[@]}"; do
    combo_tag=step_${step_size}_theta_${theta}
    echo "|=== Running combo: ${combo_tag} ===|"

       

mkdir -p ${trials_dir}/${combo_tag}





    echo "|------------- 4) [ORIG-SPACE] Generate local tractography with inclusions ROI (${combo_tag}) -------------|"
    for nside in left right; do
        # Single run per side per combo 
        scil_tracking_local \
            ${subject_dir}/tractoflow/${nsub}__fodf.nii.gz \
            ${orig_rois_dir}/${nsub}_cp_${nside}_orig.nii.gz \
            ${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz \
            #${orig_tracking_dir}/${combo_tag}/orig/${nsub}_${nside}_from_cp_${combo_tag}.trk \
            ${trials_dir}/${combo_tag}/${nsub}_${nside}_from_cp_${combo_tag}.trk \
            --npv $npv_per_run \
            --step ${step_size} \
            --theta ${theta} \
            ${gpu} -v -f
    done
    echo "|------------- 4) Done -------------|"
    echo ""
  done
done

# =========================
#  MERGE ALL 9 TRACKINGS
# =========================

echo "|============= Concatenating All orig-space tractograms  =============|"
merged_orig_dir=${orig_tracking_dir}/final_merged
mkdir -p "${merged_orig_dir}"

for nside in left right; do
    files=()
    for step_size in "${step_list[@]}"; do
      for theta in "${theta_list[@]}"; do
        combo_tag=step_${step_size}_theta_${theta}
        f="${trials_dir}/${combo_tag}/${nsub}_${nside}_from_cp_${combo_tag}.trk"

        if [[ -f "$f" ]]; then
            files+=("$f")
        fi 
      done
    done

    if (( ${#files[@]} )); then
          scil_tractogram_math -f concatenate "${files[@]}" \
            "${merged_orig_dir}/${nsub}_${nside}_from_cp_merged.trk"
    else
        echo "WARN: No files found for ${nside}"
    fi
done
    echo "|============= Merging done =============|"
    echo ""

   
# =========================
#  REGISTER MERGED TRACTOGRAM
# =========================
    

echo "|------------- 5) Register merged Tracking in MNI space -------------|"
merged_mni_dir="${mni_tracking_root}/final_merged"



mkdir -p \
    "${merged_mni_dir}/orig" \
    "${merged_mni_dir}/filtered" \
    "${merged_mni_dir}/final" \
    "${merged_mni_dir}/segmented/ROIs" \
    "${merged_mni_dir}/segmented/outliers" \
    "${merged_mni_dir}/segmented/length"


for nside in left right; do
  in_trk="${merged_orig_dir}/${nsub}_${nside}_from_cp_merged.trk"
  out_trk="${merged_mni_dir}/orig/${nsub}_${nside}_from_cp_merged.trk"

  if [[ -f "${in_trk}" ]]; then
    scil_tractogram_apply_transform \
      "${in_trk}" \
      "${mni_dir}/MNI/mni_masked.nii.gz" \
      "${output_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" \
      "${out_trk}" \
      --in_deformation "${output_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz" \
      --remove_invalid \
      --reverse_operation -f
  fi
done
echo "|------------- 5) Done -------------|"
echo ""


# =========================
#  FILTER FINAL MERGED ONLY ONCE
# =========================

echo "|------------- 6) [MNI-SPACE] Filtering merged tractograms -------------|"


for nside in left right; do
   in_trk="${merged_mni_dir}/orig/${nsub}_${nside}_from_cp_merged.trk"
   out_trk="${merged_mni_dir}/filtered/${nsub}_${nside}_from_cp_filtered.trk"

 # decide the opposite side
   if [[ "${nside}" == "left" ]]; then
       opp_side="right"
   else
       opp_side="left"
   fi

   if [[ -f "${in_trk}" ]]; then      
   scil_tractogram_filter_by_roi "${in_trk}" "${out_trk}" \
       --drawn_roi "${mni_rois_dir}/${nsub}_any_exclusion_roi_mni.nii.gz" 'any' 'exclude' \
       --drawn_roi "${mni_rois_dir}/${nsub}_${opp_side}_cerebellum_wm_mni.nii.gz" 'any' 'exclude' \
       --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_cerebellum_wm_mni.nii.gz" 'either_end' 'exclude' \
       --drawn_roi "${mni_dir}/MNI/midsagittal_plane.nii.gz" 'any' 'exclude' -f
   else
     echo "WARN: ${in_trk} not found, skipping filtering for ${nside}"
   fi 

done
    echo "|------------- 6) Done -------------|"
    echo ""

# Remove intermediate MNI-space "orig" tractograms to avoid duplication on disk
rm -rf "${merged_mni_dir}/orig"

# =========================
#  SEGMENT ONLY FINAL FILTERED FILE
# =========================

echo "|------------- 7) [MNI-SPACE] Segmentation -------------|"
for nside in left right; do 
    fin="${merged_mni_dir}/filtered/${nsub}_${nside}_from_cp_filtered.trk"

    ## Mesencephalic Tract (Top)
     scil_tractogram_filter_by_roi "${fin}" \
       "${merged_mni_dir}/segmented/ROIs/${nsub}_${nside}_mesencephalic.trk"  \
       --drawn_roi ${mni_dir}/MNI/upper_cut_brainstem.nii.gz 'any' 'include' \
       --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
       --drawn_roi ${mni_dir}/MNI/coronal_plane_for_mesencephalic.nii.gz 'any' 'include' -f


    ## Spinal Tract (bottom)
    scil_tractogram_filter_by_roi "${fin}" \
      "${merged_mni_dir}/segmented/ROIs/${nsub}_${nside}_spinal.trk" \
      --drawn_roi ${mni_dir}/MNI/lowest_cut_brainstem.nii.gz 'any' 'include' \
      --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' -f

    ## Two remaining nucleus/tract
    scil_tractogram_filter_by_roi "${fin}" \
      "${merged_mni_dir}/segmented/ROIs/${nsub}_${nside}_remaining_cp.trk" \
      --drawn_roi ${mni_dir}/MNI/lower_cut_brainstem.nii.gz 'any' 'exclude' \
      --drawn_roi ${mni_dir}/MNI/upper_cut_brainstem.nii.gz 'any' 'exclude' \
      --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
      --bdo ${mni_dir}/MNI/sphere_exclusion_for_remaining_cp.bdo 'any' 'exclude' -f
done

echo "|------------- 7) Done -------------|"
echo ""

# =========================
#  CLEANING FINAL BUNDLES
# =========================

echo "|------------- 8) Cleaning ---------|"

for nside in left right; do
    scil_bundle_reject_outliers \
     
      

      "${merged_mni_dir}/segmented/ROIs/${nsub}_${nside}_mesencephalic.trk" \
      "${merged_mni_dir}/segmented/outliers/${nsub}_${nside}_mesencephalic.trk" -f

    cp "${merged_mni_dir}/segmented/outliers/${nsub}_${nside}_mesencephalic.trk" \
       "${merged_mni_dir}/final/${nsub}_${nside}_mesencephalic.trk"



    scil_bundle_reject_outliers \
      "${merged_mni_dir}/segmented/ROIs/${nsub}_${nside}_spinal.trk" \
      "${merged_mni_dir}/segmented/outliers/${nsub}_${nside}_spinal.trk" -f

    cp "${merged_mni_dir}/segmented/outliers/${nsub}_${nside}_spinal.trk" \
       "${merged_mni_dir}/final/${nsub}_${nside}_spinal.trk"



    # Length-filter ONLY the left spinal bundle: 61 < length < 66 mm
    if [ "${nside}" = "left" ]; then
       scil_tractogram_filter_by_length \
         
         #"${merged_mni_dir}/final/${nsub}_${nside}_spinal_len61_66.trk" \
         

         "${merged_mni_dir}/final/${nsub}_${nside}_spinal.trk" \
         "${merged_mni_dir}/segmented/length/${nsub}_${nside}_spinal_len61_66.trk" \
         --minL 61 --maxL 66 --display_counts


       # Overwrite the original filename so the rest of the pipeline (incl. merge) stays unchanged
        #mv -f \
          #"${merged_mni_dir}/final/${nsub}_${nside}_spinal_len61_66.trk" \
          #"${merged_mni_dir}/final/${nsub}_${nside}_spinal.trk"
       cp -f \
        "${merged_mni_dir}/segmented/length/${nsub}_${nside}_spinal_len61_66.trk" \
        "${merged_mni_dir}/final/${nsub}_${nside}_spinal.trk"

    fi

    

    scil_bundle_reject_outliers \
       "${merged_mni_dir}/segmented/ROIs/${nsub}_${nside}_remaining_cp.trk" \
       "${merged_mni_dir}/segmented/outliers/${nsub}_${nside}_remaining_cp.trk" \
       --alpha 0.97

    cp_final="${merged_mni_dir}/final/${nsub}_${nside}_remaining_cp.trk"
    cp "${merged_mni_dir}/segmented/outliers/${nsub}_${nside}_remaining_cp.trk" "${cp_final}"  

    

    if [[ -f "${cp_final}" ]]; then
        scil_tractogram_filter_by_orientation \
          "${cp_final}" \
          "${cp_final}" \
          --max_z 7 --use_abs -f
    else
        echo "WARN: ${cp_final} does not exist, skipping orientation filtering for ${nside}"
    fi
done

    echo "|------------- 8) Done  -------------|"
    echo ""
    


# =========================
#  9) BRING FINAL MNI BUNDLES BACK TO ORIG SPACE
# =========================

echo "|------------- 9) [BACK-TO-ORIG] Register final MNI bundles to orig space -------------|"

orig_final_dir="${orig_tracking_dir}/final_merged/final"
mkdir -p "${orig_final_dir}"

# We know:
# - MNI final bundles: ${merged_mni_dir}/final/${nsub}_${nside}_<bundle>.trk
# - We want orig-space final bundles in: ${orig_final_dir}/${nsub}_${nside}_<bundle>_orig.trk

for nside in left right; do
    for bundle in mesencephalic spinal remaining_cp; do
        in_trk="${merged_mni_dir}/final/${nsub}_${nside}_${bundle}.trk"
        out_trk="${orig_final_dir}/${nsub}_${nside}_${bundle}_orig.trk"

        if [[ -f "${in_trk}" ]]; then
            scil_tractogram_apply_transform \
              "${in_trk}" \
              "${subject_dir}/tractoflow/${nsub}__t1_warped.nii.gz" \
              "${output_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" \
              "${out_trk}" \
              --in_deformation "${output_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz" \
              --remove_invalid -f
        else
            echo "WARN: Final MNI bundle not found for ${nside} ${bundle}, skipping back-to-orig."
        fi
    done
done

echo "|------------- 9) Done -------------|"
echo ""
