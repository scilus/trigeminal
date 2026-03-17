#!/bin/bash

#
# This version of the script performs a Bundle-Specific Tractography (BST)
# based on the a priori of the trigeminal_nerve atlas. That is, the atlas
# is used as prior to seed aggressively and define masks to help tractography.
# The FA threshold is lowered to 0.1 to track everywhere.

# Input structure
#
#    [input]
#    ├── sub-01
#    │   ├── freesurfer
#    │   │   └─── aparc.DKTatlas+aseg.mgz
#    │   ├── tractoflow
#    │   │   └─── sub-01__fa.nii.gz
#    │   │   └─── sub-01__fodf.nii.gz
#    │   │   └─── sub-01__t1_warped.nii.gz

#
# Example of running the script
# trigeminal_apply_atlas.sh -s input/S1/ -m ~/Research/Source/trigeminal/ROIs_clean/ -a ~/Research/Source/trigeminal/atlas/ -o output_atlas/S1/ -t 8 -g true


usage() { echo "$(basename $0) [-s path/to/subject] [-m path/to/trigeminal/ROIs_clean] [-a path/to/trigeminal/atlas] [-o output_dir] [-i sides] [-t nb_threads] [-p step_size] [-e theta] [-f fa_threshold] [-n npv_first_order] [-g] (if you have a gpu)" 1>&2; exit 1; }

while getopts "s:m:a:o:t:i:p:e:f:n:g:" args; do
    case "${args}" in
        s) s=${OPTARG};;
        m) m=${OPTARG};;
        a) a=${OPTARG};;
        o) o=${OPTARG};;
        t) t=${OPTARG};;
        i) choose_sides=${OPTARG};;
        p) step_size=${OPTARG} ;;
        e) theta=${OPTARG} ;;
        f) fa_threshold=${OPTARG} ;;
        n) npv_first_order=${OPTARG} ;;
        g) g=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${s}" ] || [ -z "${m}" ] || [ -z "${o}" ] || [ -z "${a}" ]; then
    usage
fi

gpu=""
if [ -n "${g}" ]; then
    gpu="--use_gpu"
fi

if [ -n "${step_size}" ] && [ -n "${theta}" ]; then
    step_list=(${step_size})
    theta_list=(${theta})
else
    step_list=(0.1 0.5 1.0)
    theta_list=(20 30 40)
fi

# TODO: Maybe go with 0.15 as Nasrine. To test.
# TODO: "hard" we need a real map-include/map-exclude map to run PFT. The nerve is CSF.
#        PFT would allow the tracking to bounce of the CSF to continue tracking
fa_threshold=${fa_threshold:-0.1}
npv_first_order=${npv_first_order:-450}

# npv_first_order is the total number of seeds per voxel for the whole first-order tracking.
# It is divided by the number of step/theta combinations.
# The resulting npv_per_run is the number of seeds actually used in each run.

# number of step/theta combos
n_combos=$(( ${#step_list[@]} * ${#theta_list[@]} ))
# seeds per combo (rounded)
npv_per_run=$(( (npv_first_order + n_combos - 1) / n_combos ))  # ceiling division
echo "Using $npv_per_run seeds per voxel per run (based on $npv_first_order total)"

subject_dir=${s}
atlas_dir=${a}
mni_dir=${m}
out_dir=${o}
nb_thread=${t}

echo "Folder subject: " ${subject_dir}
echo "Folder MNI: " ${mni_dir}
echo "Atlas MNI: " ${atlas_dir}
echo "Output folder: " ${out_dir}
echo "GPU: " ${gpu}

#      Lets wait to have ensemble tractography before augmenting npv more.
# Observations: at 2mm iso, on Melodie's data, I have a feeling mesencephalic could need NPV > 400. Other parts are ok.
# On 1.7 iso, npv 400 seems like a trade-off.
# HCP, the mesencephalic part is the tiniest. Could take more npv. The turn is hard to make.
# The npv could be adjusted depending on the track of interest. The easiest and thickess is clearly the spinal.
# TODO: once Nasrine's PR is in, we mimick her ensemble tractography with 9 local tracking and npv/9

sides=${choose_sides:-"left right"}
opposite_side=leftright

mkdir -p ${out_dir}/orig_space/{bundles_mask,transfo}
mkdir -p ${out_dir}/{orig_space,mni_space}/rois
mkdir -p ${out_dir}/mni_space/tractograms/{orig,filtered,length,outliers,segmented,final}
mkdir -p ${out_dir}/orig_space/tractograms/{orig,final,template}

orig_rois_dir=${out_dir}/orig_space/rois
mni_rois_dir=${out_dir}/mni_space/rois
orig_tracking_dir=${out_dir}/orig_space/tractograms/
mni_tracking_dir=${out_dir}/mni_space/tractograms/

echo "|------------- 1) Register mni_masked in orig space -------------|"
antsRegistrationSyN.sh -d 3 \
    -f ${subject_dir}/tractoflow/*__t1_warped.nii.gz \
    -m ${mni_dir}/MNI/mni_masked.nii.gz \
    -t s \
    -o ${out_dir}/orig_space/transfo/2orig_ \
    -y 1 \
    -n ${nb_thread} >> ${out_dir}/orig_space/transfo/2orig_log.txt

echo "|------------- 1.1) Reshape aparc.DKTatlas+aseg.mgz orig space -------------|"
## [ORIG-SPACE] Reshape aparc.DKTatlas+aseg.mgz
scil_volume_reslice_to_reference ${subject_dir}/freesurfer/aparc.DKTatlas+aseg.mgz \
    ${subject_dir}/tractoflow/*__t1_warped.nii.gz \
    ${orig_rois_dir}/aparc.DKTatlas+aseg_orig.nii.gz \
    --interpolation nearest --keep_dtype -f
echo "|------------- 1) Done -------------|"
echo ""

echo "|------------- 2) Generate exclusions and inclusions ROI -------------|"
## [ORIG-SPACE] Exclusions ROI labels
Left_Cerebral_Cortex=($(seq 1000 1035)) # warning for missing labels is normal
Right_Cerebral_Cortex=($(seq 2000 2035))
Left_Cerebral_WM=(2)
Right_Cerebral_WM=(41)
Cerebellum_Cortex=(8 47)
Right_Cerebellum_WM=(46)
Left_Cerebellum_WM=(7)
Any_Exclusion_ROI=(${Left_Cerebral_Cortex[*]} ${Right_Cerebral_Cortex[*]} ${Left_Cerebral_WM[*]} ${Right_Cerebral_WM[*]} ${Cerebellum_Cortex[*]}) #Generate bilateral exclusion array

## Exclusions ROIs masks
## Generate bilateral exclusion ROI
echo "|------------- 2.1) any_exclusion_roi_orig -------------|"
scil_labels_combine ${orig_rois_dir}/any_exclusion_roi_orig.nii.gz \
    --volume_ids ${orig_rois_dir}/aparc.DKTatlas+aseg_orig.nii.gz ${Any_Exclusion_ROI[*]} \
    --merge_groups -f

echo "|------------- 2.2) cerebellum_wm_right_orig -------------|"
## WM Cerebellum Right
scil_labels_combine ${orig_rois_dir}/right_cerebellum_wm_orig.nii.gz \
    --volume_ids ${orig_rois_dir}/aparc.DKTatlas+aseg_orig.nii.gz ${Right_Cerebellum_WM[*]} \
    --merge_groups -f

## WM Cerebellum Left
echo "|------------- 2.3) cerebellum_wm_left_orig -------------|"
scil_labels_combine ${orig_rois_dir}/left_cerebellum_wm_orig.nii.gz \
    --volume_ids ${orig_rois_dir}/aparc.DKTatlas+aseg_orig.nii.gz ${Left_Cerebellum_WM[*]} \
    --merge_groups -f

# WM mask
scil_volume_math lower_threshold ${subject_dir}/tractoflow/*__fa.nii.gz \
	${fa_threshold} \
	${orig_rois_dir}/wm_mask_${fa_threshold}_orig.nii.gz \
    --data_type uint8 -f

# Register these ROI in MNI space
echo "|------------- 2.4) Registration in MNI space -------------|"
for roi in right_cerebellum_wm left_cerebellum_wm any_exclusion_roi
do
    antsApplyTransforms -d 3 \
		-i ${orig_rois_dir}/${roi}_orig.nii.gz \
		-r ${mni_dir}/MNI/mni_masked.nii.gz \
		-t [${out_dir}/orig_space/transfo/2orig_0GenericAffine.mat, 1] \
		-t ${out_dir}/orig_space/transfo/2orig_1InverseWarp.nii.gz \
		-o ${mni_rois_dir}/${roi}_mni.nii.gz \
		-n NearestNeighbor;

    scil_volume_math convert \
		${mni_rois_dir}/${roi}_mni.nii.gz \
		${mni_rois_dir}/${roi}_mni.nii.gz --data_type int16 -f
done
echo "|------------- 2) Done -------------|"
echo ""


echo "|------------- 3) Register atlas components  -------------|"
for mask_type in bundles_mask
do
    echo "|------------- 3.1) Component : Register atlas components from folder ${mask_type} in orig space -------------|"
    for component in ${atlas_dir}/${mask_type}/*;
    do
        atlas_component=`basename "$component"`
        echo "|------------- 3.2) Component : ${atlas_component} -------------|"
        antsApplyTransforms -d 3 \
            -i ${component} \
            -r ${subject_dir}/tractoflow/*__t1_warped.nii.gz \
            -t ${out_dir}/orig_space/transfo/2orig_1Warp.nii.gz \
            -t ${out_dir}/orig_space/transfo/2orig_0GenericAffine.mat \
            -o ${out_dir}/orig_space/${mask_type}/${atlas_component} \
            -n NearestNeighbor;

        scil_volume_math convert \
            ${out_dir}/orig_space/${mask_type}/${atlas_component} \
            ${out_dir}/orig_space/${mask_type}/${atlas_component} \
            --data_type int16 -f;
    done
done
echo "|------------- 3) Done -------------|"
echo ""

# TODO: real BST
#     i) generate dilated endpoints mask
#     ii) TODI maps from tempalte
#     iii) e-fodf. I think this step will help tracking take the turn, for e.g. on the mesencephalic part
#     iv) perform local tracking seeding from dilated union masks with endpoints, using e-fodf and a low FA threshold
echo "|------------- 4) Tracking from atlas component  -------------|"
for component in spinal.nii.gz remaining_cp.nii.gz mesencephalic.nii.gz #${atlas_dir}/bundles_mask/*;
do
    for nside in ${sides}
    do
        atlas_component=${nside}_${component}
        for step_size in "${step_list[@]}"; do
            for theta in "${theta_list[@]}"; do
                combo_tag=step_${step_size}_theta_${theta}

                if [[ ${component} == "mesencephalic.nii.gz" ]] || [[ ${component} == "spinal.nii.gz" ]]; then
                    npv_per_run=$(( npv_per_run + 50 ))
                    echo "Note: mesencephalic could benefit from a higher npv now npv_per_run is ${npv_per_run}."
                fi

                echo "|------------- 4.1) Tracking from atlas component ${atlas_component} with npv=${npv_per_run}, step=${step_size}, theta=${theta} -------------|"
                scil_tracking_local ${subject_dir}/tractoflow/*__fodf.nii.gz \
                    ${out_dir}/orig_space/bundles_mask/${atlas_component} \
                    ${orig_rois_dir}/wm_mask_${fa_threshold}_orig.nii.gz \
                    ${out_dir}/orig_space/tractograms/orig__${combo_tag}_${atlas_component/.nii.gz/.trk} \
                    --npv ${npv_per_run} \
                    --step ${step_size} \
                    --theta ${theta} \
                    --min_length 8 --max_length 100 \
                    ${gpu} -f -v
            done
        done
    done
done

for component in mesencephalic.nii.gz spinal.nii.gz remaining_cp.nii.gz # ${atlas_dir}/bundles_mask/*;
do
    for nside in ${sides}
    do
        atlas_component=${nside}_${component}
        echo "|------------- 4.2) Concatenate tracking results for atlas component ${atlas_component} -------------|"
        scil_tractogram_math lazy_concatenate\
            ${out_dir}/orig_space/tractograms/orig__*_${atlas_component/.nii.gz/.trk} \
            ${out_dir}/orig_space/tractograms/orig__${atlas_component/.nii.gz/.trk} -f

        echo "|------------- 4.3) Register tracking to MNI space -------------|"
        scil_tractogram_apply_transform \
            ${out_dir}/orig_space/tractograms/orig__${atlas_component/.nii.gz/.trk} \
            ${mni_dir}/MNI/mni_masked.nii.gz \
            ${out_dir}/orig_space/transfo/2orig_0GenericAffine.mat \
            ${mni_tracking_dir}/orig__${atlas_component/.nii.gz/.trk} \
            --in_deformation ${out_dir}/orig_space/transfo/2orig_1Warp.nii.gz \
            --remove_invalid \
            --reverse_operation -f
    done
done

# moving tractogram to the orig folder
mv ${out_dir}/orig_space/tractograms/orig_*trk ${out_dir}/orig_space/tractograms/orig/

echo "|------------- 4) Done -------------|"
echo ""

echo "|------------- 5) Filtering & Segmenting -------------|"
echo "|------------- 5.1) Major filtering for mesencephalic, remaining cp, spinal -------------|"
# Main filter for mesencephalic, remaining cp, spinal
for atlas_component in mesencephalic spinal remaining_cp
do
    for nside in ${sides}
    do
        echo "|------------- 5.1b) Major filtering for ${atlas_component} nside: ${nside} npv: ${npv} -------------|"
        scil_tractogram_filter_by_roi ${mni_tracking_dir}/orig__${nside}_${atlas_component}.trk \
            ${mni_tracking_dir}/filtered_${nside}_${atlas_component}.trk \
            --drawn_roi ${mni_rois_dir}/any_exclusion_roi_mni.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_rois_dir}/${nside/${nside}/${opposite_side/${nside}/}}_cerebellum_wm_mni.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_rois_dir}/${nside}_cerebellum_wm_mni.nii.gz 'either_end' 'exclude' \
            --drawn_roi ${mni_dir}/MNI/midsagittal_plane.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_dir}/MNI/cp_${nside}_bin.nii.gz 'any' 'include' -f -v # It was not needed when creating atlas since everything was tracked from cp
    done
done

# Filter mesencephalic
echo "|------------- 5.2) Segmentation for mesencephalic -------------|"
for nside in ${sides}
do

    # choose thalamus ROI to exclude based on side
    if [[ "${nside}" == "left" ]]; then
        thal_roi="${mni_dir}/MNI/left_thalamus_231_245_brainnetome__to_mni_masked.nii.gz"
    else
        thal_roi="${mni_dir}/MNI/right_thalamus_232_246_brainnetome__to_mni_masked.nii.gz"
    fi

    scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered_${nside}_mesencephalic.trk \
        ${mni_tracking_dir}/segmented_${nside}_mesencephalic.trk  \
        --drawn_roi ${mni_dir}/MNI/upper_cut_brainstem.nii.gz 'any' 'include' \
        --drawn_roi ${mni_dir}/MNI/new_coronal.nii.gz 'any' 'include' \
        --drawn_roi ${mni_dir}/MNI/coronal_plane_for_mesencephalic.nii.gz 'any' 'include' \
        --drawn_roi ${mni_dir}/MNI/csf_mask.nii.gz 'any' 'exclude' \
        --drawn_roi ${thal_roi} 'any' 'exclude' \
        -f

    # Filter by length (35-70)
    scil_tractogram_filter_by_length \
        ${mni_tracking_dir}/segmented_${nside}_mesencephalic.trk \
        ${mni_tracking_dir}/final_${nside}_mesencephalic.trk \
        --minL 35 --maxL 70 --display_counts -f

#    scil_bundle_reject_outliers \
#        ${mni_tracking_dir}/segmented_${nside}_mesencephalic.trk \
#        ${mni_tracking_dir}/final_${nside}_mesencephalic.trk -f
done

# Filter remaining_cp
echo "|------------- 5.3) Segmentation remaining cp -------------|"
for nside in ${sides}
do
    scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered_${nside}_remaining_cp.trk \
        ${mni_tracking_dir}/segmented_${nside}_remaining_cp.trk  \
        --drawn_roi ${mni_dir}/MNI/new_lower.nii.gz 'any' 'exclude' \
        --drawn_roi ${mni_dir}/MNI/upper_cut_brainstem.nii.gz 'any' 'exclude' \
        --drawn_roi ${mni_dir}/MNI/new_coronal.nii.gz 'any' 'include' \
        --bdo ${mni_dir}/MNI/sphere_exclusion_for_remaining_cp.bdo 'any' 'exclude' \
        --drawn_roi ${mni_dir}/MNI/csf_mask.nii.gz 'any' 'exclude' \
        --drawn_roi ${mni_dir}/MNI/cp_${nside}_bin.nii.gz 'either_end' 'include' \
        -f

    # Filter by length (8-32)
    scil_tractogram_filter_by_length \
        ${mni_tracking_dir}/segmented_${nside}_remaining_cp.trk \
        ${mni_tracking_dir}/final_${nside}_remaining_cp.trk \
        --minL 8 --maxL 32 --display_counts -f

    # TODO: we apply the length and ROI rules of Nasrine in the future.
#   scil_bundle_reject_outliers \
#        ${mni_tracking_dir}/segmented_${nside}_remaining_cp.trk \
#        ${mni_tracking_dir}/final_${nside}_remaining_cp.trk

    scil_tractogram_filter_by_orientation \
        ${mni_tracking_dir}/final_${nside}_remaining_cp.trk  \
        ${mni_tracking_dir}/final_${nside}_remaining_cp.trk  \
        --max_z 7 --use_abs -f
done

# Filter spinal
echo "|------------- 5.4) Segmentation spinal -------------|"
for nside in ${sides}
do
    scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered_${nside}_spinal.trk \
        ${mni_tracking_dir}/segmented_${nside}_spinal.trk \
        --drawn_roi ${mni_dir}/MNI/new_lower.nii.gz 'any' 'include' \
        --drawn_roi ${mni_dir}/MNI/new_coronal.nii.gz 'any' 'include' \
        --drawn_roi ${mni_dir}/MNI/csf_mask.nii.gz 'any' 'exclude' \
        --drawn_roi ${mni_dir}/MNI/cp_${nside}_bin.nii.gz 'either_end' 'include' \
        -f

    # TODO: we need better filetring with ROI of the spinal part. There is a branch that we pick that we should not have.


    # Filter by length (10-55)
    scil_tractogram_filter_by_length \
        ${mni_tracking_dir}/segmented_${nside}_spinal.trk \
        ${mni_tracking_dir}/final_${nside}_spinal.trk \
        --minL 10 --maxL 55 --display_counts -f

#    scil_bundle_reject_outliers \
#        ${mni_tracking_dir}/segmented_${nside}_spinal.trk \
#        ${mni_tracking_dir}/final_${nside}_spinal.trk \
#        -f
done


echo "|------------- 6 Transform back final MNI trk to orig space -------------|"
for atlas_component in mesencephalic spinal remaining_cp
do
    for nside in ${sides}
    do
        scil_tractogram_apply_transform ${mni_tracking_dir}/final_${nside}_${atlas_component}.trk \
            ${subject_dir}/tractoflow/*__t1_warped.nii.gz \
            ${out_dir}/orig_space/transfo/2orig_0GenericAffine.mat \
            ${orig_tracking_dir}/final/final_${nside}_${atlas_component}.trk \
            --inverse \
            --in_deformation ${out_dir}/orig_space/transfo/2orig_1InverseWarp.nii.gz \
            --remove_invalid -f
    done
done

#echo "|------------- 7) Transform back final MNI trk from ATLAS to orig space -------------|"
for atlas_component in mesencephalic spinal remaining_cp
do
    for nside in ${sides}
    do
        scil_tractogram_apply_transform ${atlas_dir}/trks/first_order/${nside}_${atlas_component}.trk \
            ${subject_dir}/tractoflow/*__t1_warped.nii.gz \
            ${out_dir}/orig_space/transfo/2orig_0GenericAffine.mat \
            ${orig_tracking_dir}/template/template_${nside}_${atlas_component}.trk \
            --inverse \
            --in_deformation ${out_dir}/orig_space/transfo/2orig_1InverseWarp.nii.gz \
            --remove_invalid -f
    done
done

echo "|------------- 8) Register bundle proba template  -------------|"
for component in mesencephalic spinal remaining_cp
do
    for nside in ${sides}
    do
        atlas_component=${nside}_${component}
        echo "|------------- 8.1) Component : Register prob masks from template ${atlas_component} in orig space -------------|"
        antsApplyTransforms -d 3 \
            -i ${atlas_dir}/proba/first_order/${atlas_component}_prob.nii.gz \
            -r ${subject_dir}/tractoflow/*__t1_warped.nii.gz \
            -t ${out_dir}/orig_space/transfo/2orig_1Warp.nii.gz \
            -t ${out_dir}/orig_space/transfo/2orig_0GenericAffine.mat \
            -o ${out_dir}/orig_space/bundles_mask/${atlas_component}_prob.nii.gz \
            -n NearestNeighbor \
            -u float;
    done
done

# TODO: need better cleaning up
#  - better ROIs first
#  - cut. too aggressive? Could be optional. If we do, we start with that.
#  - length
#  - recobundle will not work because of field of view. I tried...

# Last organization move of files in the proper output directories in mni_space
mv  ${mni_tracking_dir}/orig_*  ${mni_tracking_dir}/orig/
mv  ${mni_tracking_dir}/filtered_*  ${mni_tracking_dir}/filtered/
mv  ${mni_tracking_dir}/segmented_*  ${mni_tracking_dir}/segmented/
mv  ${mni_tracking_dir}/final_*  ${mni_tracking_dir}/final/


# TODO: option to clean up temporary trk? orig_ trks are big files...
# rm -rf ${mni_tracking_dir}/orig/*trk ${orig_tracking_dir}/orig/*trk

echo "|------------- Done -------------|"
echo ""
