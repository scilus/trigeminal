#!/bin/bash
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Samir Akeb (2022-2023)
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Arnaud Bore (2023-2024)

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

usage() { echo "$(basename $0) [-s path/to/subjects] [-m path/to/mni] [-o output_dir] [-t nb_threads] -g true" 1>&2; exit 1; }

while getopts "s:m:o:t:g:" args; do
    case "${args}" in
        s) s=${OPTARG};;
        m) m=${OPTARG};;
        o) o=${OPTARG};;
        t) t=${OPTARG};;
        g) g=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${s}" ] || [ -z "${m}" ] || [ -z "${o}" ]; then
    usage
fi

subject_dir=${s}
mni_dir=${m}
out_dir=${o}
nb_thread=${t}

fa_threshold=0.20
npv_first_order=10000

opposite_side=leftright

gpu=""
if [ ! -z "${g}" ]; then
    gpu="--use_gpu"
fi


echo "Folder subjects: " ${subject_dir}
echo "Folder MNI: " ${mni_dir}
echo "Output folder: " ${out_dir}
echo "Use GPU: " ${gpu}
echo "Number of threads" ${nb_thread}

for nsub in ${subject_dir}/*/
do
    nsub=`basename "$nsub"`
    rm -rf ${out_dir}/${nsub}
    mkdir -p ${out_dir}/${nsub}/orig_space/{rois,tracking_first_order,transfo}
    mkdir -p ${out_dir}/${nsub}/orig_space/tracking_first_order/orig
    mkdir -p ${out_dir}/${nsub}/mni_space/{rois,tracking_first_order}
    mkdir -p ${out_dir}/${nsub}/mni_space/tracking_first_order/{orig,filtered,segmented,final}

    orig_rois_dir=${out_dir}/${nsub}/orig_space/rois
    mni_rois_dir=${out_dir}/${nsub}/mni_space/rois
    orig_tracking_dir=${out_dir}/${nsub}/orig_space/tracking_first_order
    mni_tracking_dir=${out_dir}/${nsub}/mni_space/tracking_first_order

    echo ""
    echo "|------------- PROCESSING FIRST ORDER TGN TRACTOGRAPHY FOR ${nsub} -------------|"
    echo ""


    echo "|------------- 1) Register mni_masked in orig space -------------|"
    antsRegistrationSyNQuick.sh \
        -d 3 -f ${subject_dir}/${nsub}/tractoflow/*__t1_warped.nii.gz \
        -m ${mni_dir}/MNI/mni_masked.nii.gz \
        -t s -o ${out_dir}/${nsub}/orig_space/transfo/2orig_ \
        -y 1 \
        -n ${nb_thread} >> ${out_dir}/${nsub}/orig_space/transfo/2orig_log.txt

    ## [ORIG-SPACE] Register all ROIs -> ORIG
    for nroi in cp_left cp_right
    do
        antsApplyTransforms \
        -d 3 \
        -i ${mni_dir}/MNI/${nroi}.nii.gz \
        -r ${subject_dir}/${nsub}/tractoflow/*__t1_warped.nii.gz \
        -t ${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz \
        -t ${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat \
        -o ${orig_rois_dir}/${nsub}_${nroi}_orig.nii.gz

        scil_volume_math lower_threshold_eq ${orig_rois_dir}/${nsub}_${nroi}_orig.nii.gz 0.5 ${orig_rois_dir}/${nsub}_${nroi}_orig.nii.gz --data_type int16 -f
    done

    ## [ORIG-SPACE] Reshape aparc.DKTatlas+aseg.mgz
    scil_volume_reslice_to_reference ${subject_dir}/${nsub}/freesurfer/aparc.DKTatlas+aseg.mgz \
        ${subject_dir}/${nsub}/tractoflow/${nsub}__t1_warped.nii.gz \
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
    Any_Exclusion_ROI=(${Left_Cerebral_Cortex[*]} ${Right_Cerebral_Cortex[*]} ${Left_Cerebral_WM[*]} ${Right_Cerebral_WM[*]} ${Cerebellum_Cortex[*]}) #Generate bilateral exclusion array

    ## Exclusions ROIs masks
    ## Generate bilateral exclusion ROI

    echo "|------------- 2.1) any_exclusion_roi_orig -------------|"
    scil_labels_combine ${orig_rois_dir}/${nsub}_any_exclusion_roi_orig.nii.gz \
        --volume_ids ${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz ${Any_Exclusion_ROI[*]} \
        --merge_groups -f

    echo "|------------- 2.2) cerebellum_wm_right_orig -------------|"
    ## WM Cerebellum Right
    scil_labels_combine ${orig_rois_dir}/${nsub}_right_cerebellum_wm_orig.nii.gz \
        --volume_ids ${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz ${Right_Cerebellum_WM[*]} \
        --merge_groups -f

    ## WM Cerebellum Left
    echo "|------------- 2.3) cerebellum_wm_left_orig -------------|"
    scil_labels_combine ${orig_rois_dir}/${nsub}_left_cerebellum_wm_orig.nii.gz \
        --volume_ids ${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz ${Left_Cerebellum_WM[*]} \
        --merge_groups -f

    # WM mask
    scil_volume_math lower_threshold ${subject_dir}/${nsub}/tractoflow/${nsub}__fa.nii.gz \
        ${fa_threshold} \
        ${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz \
        --data_type uint8 -f

    echo "|------------- 2) Done -------------|"
    echo ""

    echo "|------------- 3) Register ROI in MNI space -------------|"
    for nroi in any_exclusion_roi_orig right_cerebellum_wm_orig left_cerebellum_wm_orig wm_mask_${fa_threshold}_orig aparc.DKTatlas+aseg_orig
    do
        antsApplyTransforms \
        -d 3 \
        -i ${orig_rois_dir}/${nsub}_${nroi}.nii.gz \
        -r ${mni_dir}/MNI/mni_masked.nii.gz \
        -t [${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat, 1] \
        -t ${out_dir}/${nsub}/orig_space/transfo/2orig_1InverseWarp.nii.gz \
        -o ${mni_rois_dir}/${nsub}_${nroi/orig/mni}.nii.gz \
        -n NearestNeighbor

        scil_volume_math convert ${mni_rois_dir}/${nsub}_${nroi/orig/mni}.nii.gz ${mni_rois_dir}/${nsub}_${nroi/orig/mni}.nii.gz --data_type int16 -f
    done

    echo "|------------- 3) Done -------------|"
    echo ""

    echo "|------------- 4) [ORIG-SPACE] Generate local tractography with inclusions ROI  -------------|"
    ## Tracking
    for nside in left right
    do
        scil_tracking_local ${subject_dir}/${nsub}/tractoflow/${nsub}__fodf.nii.gz \
            ${orig_rois_dir}/${nsub}_cp_${nside}_orig.nii.gz \
            ${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz \
            ${orig_tracking_dir}/orig/${nsub}_${nside}_from_cp.trk \
            --npv $npv_first_order \
            ${gpu} -v -f
    done
    echo "|------------- 4) Done -------------|"
    echo ""

    echo "|------------- 5) Register Tracking in MNI space -------------|"
    for nside in left right
    do
        scil_tractogram_apply_transform \
            ${orig_tracking_dir}/orig/${nsub}_${nside}_from_cp.trk \
            ${mni_dir}/MNI/mni_masked.nii.gz \
            ${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat \
            ${mni_tracking_dir}/orig/${nsub}_${nside}_from_cp.trk \
            --in_deformation ${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz \
            --remove_invalid \
            --reverse_operation -f
    done
    echo "|------------- 5) Done -------------|"
    echo ""


    echo "|------------- 6) [MNI-SPACE] Filter tractography -------------|"
    ## Filtration for left
    for nside in left right
    do
        scil_tractogram_filter_by_roi ${mni_tracking_dir}/orig/${nsub}_${nside}_from_cp.trk \
            ${mni_tracking_dir}/filtered/${nsub}_${nside}_from_cp_filtered.trk \
            --drawn_roi ${mni_rois_dir}/${nsub}_any_exclusion_roi_mni.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_rois_dir}/${nsub}_${nside/${nside}/${opposite_side/${nside}/}}_cerebellum_wm_mni.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_rois_dir}/${nsub}_${nside}_cerebellum_wm_mni.nii.gz 'either_end' 'exclude' \
            --drawn_roi ${mni_dir}/MNI/midsagittal_plane.nii.gz 'any' 'exclude' -f
        ## --drawn_roi ${mni_dir}/MNI/cp_${nside}.nii.gz any include - DONT NEED THIS SINCE IT'S TRACKED FROM CP


    done

    echo "|------------- 6) Done -------------|"
    echo ""

    echo "|------------- 7) [MNI-SPACE] Segmentation and cleaning tractography -------------|"

    for nside in left right
    do
	echo "|------------- 7) [MNI-SPACE] Running Segmentation - SIDE: ${nside} -------------|"
        ## Mesencephalic Tract (Top)
        scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered/${nsub}_${nside}_from_cp_filtered.trk \
            ${mni_tracking_dir}/segmented/${nsub}_${nside}_mesencephalic.trk  \
            --drawn_roi ${mni_dir}/MNI/upper_cut_brainstem.nii.gz 'any' 'include' \
            --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
            --drawn_roi ${mni_dir}/MNI/coronal_plane_for_mesencephalic.nii.gz 'any' 'include' \
            -f

        ## Spinal Tract (bottom)
        scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered/${nsub}_${nside}_from_cp_filtered.trk \
            ${mni_tracking_dir}/segmented/${nsub}_${nside}_spinal.trk  \
            --drawn_roi ${mni_dir}/MNI/lowest_cut_brainstem.nii.gz 'any' 'include' \
            --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
            -f

        ## Two remaining nucleus/tract from the cisternal portion (Main sensory nucleus & Trigeminal motor nucleus)
        scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered/${nsub}_${nside}_from_cp_filtered.trk \
            ${mni_tracking_dir}/segmented/${nsub}_${nside}_remaining_cp.trk  \
            --drawn_roi ${mni_dir}/MNI/lower_cut_brainstem.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_dir}/MNI/upper_cut_brainstem.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
            --bdo ${mni_dir}/MNI/sphere_exclusion_for_remaining_cp.bdo 'any' 'exclude' \
            -f

        echo "|------------- 7) Done - SIDE: ${nside} -------------|"
        echo ""


        echo "|------------- 8) [MNI-SPACE] Running Cleaning - SIDE: ${nside} -------------|"
        ## Mesencephalic Tract (Top)
        scil_bundle_reject_outliers \
            ${mni_tracking_dir}/segmented/${nsub}_${nside}_mesencephalic.trk \
            ${mni_tracking_dir}/final/${nsub}_${nside}_mesencephalic.trk -f

        ## Spinal Tract
        scil_bundle_reject_outliers \
            ${mni_tracking_dir}/segmented/${nsub}_${nside}_spinal.trk \
            ${mni_tracking_dir}/final/${nsub}_${nside}_spinal.trk -f

        ## Two remaining nucleus
        ## BDO
        ## Longueur

        scil_bundle_reject_outliers \
            ${mni_tracking_dir}/segmented/${nsub}_${nside}_remaining_cp.trk \
            ${mni_tracking_dir}/final/${nsub}_${nside}_remaining_cp.trk  \
            --alpha 0.97

        scil_tractogram_filter_by_orientation \
            ${mni_tracking_dir}/final/${nsub}_${nside}_remaining_cp.trk \
            ${mni_tracking_dir}/final/${nsub}_${nside}_remaining_cp.trk \
            --max_z 7 --use_abs -f

        echo "|------------- 8) Done - SIDE: ${nside} -------------|"
        echo ""
    done
done
