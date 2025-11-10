usage() { echo "$(basename $0) [-s path/to/subjects] [-m path/to/mni] [-a atlas_dir] [-o output_dir] [-g] (if you have a gpu)" 1>&2; exit 1; }

while getopts "s:m:a:o:g:" args; do
    case "${args}" in
        s) s=${OPTARG};;
        m) m=${OPTARG};;
        a) a=${OPTARG};;
        o) o=${OPTARG};;
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

subject_dir=${s}
atlas_dir=${a}
mni_dir=${m}
out_dir=${o}

echo "Folder subjects: " ${subject_dir}
echo "Folder MNI: " ${mni_dir}
echo "Atlas MNI: " ${atlas_dir}
echo "Output folder: " ${out_dir}
echo "GPU: " ${gpu}

npv_list="100 400"
opposite_side=leftright
fa_threshold=0.01
mni_roi_dir

for nsub in ${subject_dir}/*/
do
    nsub=`basename "$nsub"`
    echo "|------------- RUNNING: ${nsub} -------------|"

    mkdir -p ${out_dir}/${nsub}/orig_space/{bundles_mask,transfo}
    mkdir -p ${out_dir}/${nsub}/{orig_space,mni_space}/rois
    mkdir -p ${out_dir}/${nsub}/{orig_space,mni_space}/tractograms/{orig,filtered,segmented,final}

    orig_rois_dir=${out_dir}/${nsub}/orig_space/rois
    mni_rois_dir=${out_dir}/${nsub}/mni_space/rois
    orig_tracking_dir=${out_dir}/${nsub}/orig_space/tractograms/
    mni_tracking_dir=${out_dir}/${nsub}/mni_space/tractograms/

    echo "|------------- 1) Register mni_masked in orig space -------------|"
    antsRegistrationSyNQuick.sh -d 3 \
        -f ${subject_dir}/${nsub}/tractoflow/${nsub}__t1_warped.nii.gz \
        -m ${mni_dir}/MNI/mni_masked.nii.gz \
        -t s \
        -o ${out_dir}/${nsub}/orig_space/transfo/mni2orig_ \
        -y 1 \
        -n 10

    echo "|------------- 1.1) Reshape aparc.DKTatlas+aseg.mgz orig space -------------|"
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

    # Register these ROI in MNI space
    echo "|------------- 2.4) Registration in MNI space -------------|"
    for roi in right_cerebellum_wm left_cerebellum_wm any_exclusion_roi
    do
        antsApplyTransforms -d 3 \
            -i ${orig_rois_dir}/${nsub}_${roi}_orig.nii.gz \
            -r ${mni_dir}/MNI/mni_masked.nii.gz \
            -t [${out_dir}/${nsub}/orig_space/transfo/mni2orig_0GenericAffine.mat, 1] \
            -t ${out_dir}/${nsub}/orig_space/transfo/mni2orig_1InverseWarp.nii.gz \
            -o ${mni_rois_dir}/${nsub}_${roi}_mni.nii.gz \
            -n NearestNeighbor;

        scil_volume_math convert \
        ${mni_rois_dir}/${nsub}_${roi}_mni.nii.gz \
        ${mni_rois_dir}/${nsub}_${roi}_mni.nii.gz --data_type int16 -f
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
                -r ${subject_dir}/${nsub}/tractoflow/${nsub}__t1_warped.nii.gz \
                -t ${out_dir}/${nsub}/orig_space/transfo/mni2orig_1Warp.nii.gz \
                -t ${out_dir}/${nsub}/orig_space/transfo/mni2orig_0GenericAffine.mat \
                -o ${out_dir}/${nsub}/orig_space/${mask_type}/${nsub}__${atlas_component} \
                -n NearestNeighbor;

            scil_volume_math convert \
                ${out_dir}/${nsub}/orig_space/${mask_type}/${nsub}__${atlas_component} \
                ${out_dir}/${nsub}/orig_space/${mask_type}/${nsub}__${atlas_component} \
                --data_type int16 -f;
        done
    done
    echo "|------------- 3) Done -------------|"
    echo ""

    echo "|------------- 4) Tracking from atlas component  -------------|"
    for component in left_mesencephalic.nii.gz left_spinal.nii.gz left_remaining_cp.nii.gz right_mesencephalic.nii.gz right_spinal.nii.gz right_remaining_cp.nii.gz #${atlas_dir}/bundles_mask/*;
    do
        atlas_component=`basename "$component"`
        atlas_component=$component
        for npv in ${npv_list};
        do
            echo "|------------- 4.1) Tracking from atlas component ${atlas_component} with npv=${npv} -------------|"
            scil_tracking_local ${subject_dir}/${nsub}/tractoflow/${nsub}__fodf.nii.gz \
                ${out_dir}/${nsub}/orig_space/bundles_mask/${nsub}__${atlas_component} \
                ${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz \
                ${out_dir}/${nsub}/orig_space/tractograms/orig/${nsub}__${atlas_component/.nii.gz/_npv_${npv}.trk} \
                --npv ${npv} \
                ${gpu} -f;

            echo "|------------- 4.2) Register tracking to MNI space -------------|"
            scil_tractogram_apply_transform \
                ${out_dir}/${nsub}/orig_space/tractograms/orig/${nsub}__${atlas_component/.nii.gz/_npv_${npv}.trk} \
                ${mni_dir}/MNI/mni_masked.nii.gz \
                ${out_dir}/${nsub}/orig_space/transfo/mni2orig_0GenericAffine.mat \
                ${mni_tracking_dir}/orig/${nsub}__${atlas_component/.nii.gz/_npv_${npv}.trk} \
                --in_deformation ${out_dir}/${nsub}/orig_space/transfo/mni2orig_1Warp.nii.gz \
                --remove_invalid \
                --reverse_operation -f
        done
    done
    echo "|------------- 4) Done -------------|"
    echo ""

    echo "|------------- 5) Filtering -------------|"
    echo "|------------- 5.1) Major filtering for mesencephalic, remaining cp, spinal -------------|"
    # Main filter for mesencephalic, remaining cp, spinal
    for atlas_component in mesencephalic spinal remaining_cp
    do
        for nside in left right
        do
            for npv in ${npv_list}
            do
                echo "|------------- 5.1b) Major filtering for ${atlas_component} nside: ${nside} npv: ${npv} -------------|"
                scil_tractogram_filter_by_roi ${mni_tracking_dir}/orig/${nsub}__${nside}_${atlas_component}_npv_${npv}.trk \
                    ${mni_tracking_dir}/filtered/${nsub}_${nside}_${atlas_component}_npv_${npv}.trk \
                    --drawn_roi ${mni_rois_dir}/${nsub}_any_exclusion_roi_mni.nii.gz 'any' 'exclude' \
                    --drawn_roi ${mni_rois_dir}/${nsub}_${nside/${nside}/${opposite_side/${nside}/}}_cerebellum_wm_mni.nii.gz 'any' 'exclude' \
                    --drawn_roi ${mni_rois_dir}/${nsub}_${nside}_cerebellum_wm_mni.nii.gz 'either_end' 'exclude' \
                    --drawn_roi ${mni_dir}/MNI/midsagittal_plane.nii.gz 'any' 'exclude' \
                    --drawn_roi ${mni_dir}/MNI/cp_${nside}_bin.nii.gz any include -f -v # It was not needed when creating atlas since everything was tracked from cp
            done
        done
    done

    # Filter mesencephalic
    echo "|------------- 5.2) Filtering for mesencephalic -------------|"
    for nside in left right
    do
        for npv in ${npv_list}
        do
            scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered/${nsub}_${nside}_mesencephalic_npv_${npv}.trk \
                ${mni_tracking_dir}/segmented/${nsub}_${nside}_mesencephalic_npv_${npv}.trk  \
                --drawn_roi ${mni_dir}/MNI/upper_cut_brainstem.nii.gz 'any' 'include' \
                --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
                --drawn_roi ${mni_dir}/MNI/coronal_plane_for_mesencephalic.nii.gz 'any' 'include' \
                -f

            scil_bundle_reject_outliers \
                ${mni_tracking_dir}/segmented/${nsub}_${nside}_mesencephalic_npv_${npv}.trk \
                ${mni_tracking_dir}/final/${nsub}_${nside}_mesencephalic_npv_${npv}.trk -f
        done
    done

    # Filter remaining_cp
    echo "|------------- 5.3) Filtering remaining cp -------------|"
    for nside in left right
    do
        for npv in ${npv_list}
        do
            scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered/${nsub}_${nside}_remaining_cp_npv_${npv}.trk \
                ${mni_tracking_dir}/segmented/${nsub}_${nside}_remaining_cp_npv_${npv}.trk  \
                --drawn_roi ${mni_dir}/MNI/lower_cut_brainstem.nii.gz 'any' 'exclude' \
                --drawn_roi ${mni_dir}/MNI/upper_cut_brainstem.nii.gz 'any' 'exclude' \
                --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
                --bdo ${mni_dir}/MNI/sphere_exclusion_for_remaining_cp.bdo 'any' 'exclude' \
                -f

            scil_bundle_reject_outliers \
                ${mni_tracking_dir}/segmented/${nsub}_${nside}_remaining_cp_npv_${npv}.trk \
                ${mni_tracking_dir}/final/${nsub}_${nside}_remaining_cp_npv_${npv}.trk  \
                --alpha 0.97

            scil_tractogram_filter_by_orientation \
                ${mni_tracking_dir}/final/${nsub}_${nside}_remaining_cp_npv_${npv}.trk  \
                ${mni_tracking_dir}/final/${nsub}_${nside}_remaining_cp_npv_${npv}.trk  \
                --max_z 7 --use_abs -f
        done
    done

    # Filter spinal
    echo "|------------- 5.4) Filtering spinal -------------|"
    for nside in left right
    do
        for npv in ${npv_list}
        do
            scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered/${nsub}_${nside}_spinal_npv_${npv}.trk \
                ${mni_tracking_dir}/segmented/${nsub}_${nside}_spinal_npv_${npv}.trk \
                --drawn_roi ${mni_dir}/MNI/lowest_cut_brainstem.nii.gz 'any' 'include' \
                --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
                -f
            scil_bundle_reject_outliers \
                ${mni_tracking_dir}/segmented/${nsub}_${nside}_spinal_npv_${npv}.trk \
                ${mni_tracking_dir}/final/${nsub}_${nside}_spinal_npv_${npv}.trk \
                -f
        done
    done
    echo "|------------- 5) Done -------------|"
    echo ""

    echo "|------------- TRACKING FROM TRIGEMINAL ATLAS FOR ${nsub} IS COMPLETED -------------|"
    echo ""
    echo ""
done
