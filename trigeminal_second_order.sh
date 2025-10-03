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

usage() { echo "$(basename $0) [-s path/to/subjects] [-m path/to/mni] [-o output_dir] [-g] (if you have a gpu)" 1>&2; exit 1; }

while getopts "s:m:o:g:" args; do
    case "${args}" in
        s) s=${OPTARG};;
        m) m=${OPTARG};;
        o) o=${OPTARG};;
        g) g=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${s}" ] || [ -z "${m}" ] || [ -z "${o}" ]; then
    usage
fi

gpu=""
if [ -n "${g}" ]; then
    gpu="--use_gpu"
fi

subject_dir=${s}
mni_dir=${m}
out_dir=${o}

fa_threshold=0.20

npv_from_spinal_track_long=1000
npv_from_spinal_track_short=100
npv_from_thalamus_track=500

echo "Folder subjects: " ${subject_dir}
echo "Folder MNI: " ${mni_dir}
echo "Output folder: " ${out_dir}
echo "GPU: " ${gpu}

mkdir -p ${out_dir}/mni_space/tracking_first_order/final/

for nside in left right
do
    scil_tractogram_math union ${out_dir}/*/mni_space/tracking_first_order/final/*_${nside}_spinal.trk \
        ${out_dir}/mni_space/tracking_first_order/final/all_${nside}_spinal.trk -f

    # Building seeding mask by converting spinal bundle TRK to NII
    scil_tractogram_compute_density_map \
        ${out_dir}/mni_space/tracking_first_order/final/all_${nside}_spinal.trk \
        ${out_dir}/mni_space/tracking_first_order/final/all_${nside}_spinal_density_second_order_seed_mni.nii.gz \
        --binary -f

    scil_tractogram_math union ${out_dir}/*/mni_space/tracking_first_order/final/*_${nside}_remaining_cp.trk \
        ${out_dir}/mni_space/tracking_first_order/final/all_${nside}_remaining_cp.trk -f

    # Building mask by converting remaining_cp bundle TRK to NII
    scil_tractogram_compute_density_map \
        ${out_dir}/mni_space/tracking_first_order/final/all_${nside}_remaining_cp.trk \
        ${out_dir}/mni_space/tracking_first_order/final/all_${nside}_remaining_cp_density_mni.nii.gz \
        --binary -f
done

for nsub in ${subject_dir}/*/
do
    nsub=`basename "$nsub"`

    mkdir -p ${out_dir}/${nsub}/orig_space/{rois,tracking_second_order,transfo}
    mkdir -p ${out_dir}/${nsub}/orig_space/tracking_second_order/orig
    mkdir -p ${out_dir}/${nsub}/mni_space/{rois,tracking_second_order}
    mkdir -p ${out_dir}/${nsub}/mni_space/tracking_second_order/{orig,filtered,segmented,final,cut}

    orig_rois_dir=${out_dir}/${nsub}/orig_space/rois
    mni_rois_dir=${out_dir}/${nsub}/mni_space/rois
    orig_tracking_dir=${out_dir}/${nsub}/orig_space/tracking_second_order
    mni_tracking_dir_first_order=${out_dir}/${nsub}/mni_space/tracking_first_order
    mni_tracking_dir_second_order=${out_dir}/${nsub}/mni_space/tracking_second_order

    echo ""
    echo "|------------- PROCESSING SECOND-ORDER FIBERS TRACTOGRAPHY FOR TRIGEMINAL SYSTEM -------------|"
    echo "|------------- FOR DATASET ${nsub} -------------|"
    echo ""

    echo "|------------- 1) Generate local tractography with spinal bundle  -------------|"
    echo "|------------- 1.2) Seeding Mask + Registration in orig space  -------------|"
    for nside in left right
    do
    	cp ${out_dir}/mni_space/tracking_first_order/final/all_${nside}_spinal_density_second_order_seed_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz

        # Register density into orig space
        antsApplyTransforms \
            -d 3 \
            -i ${out_dir}/mni_space/tracking_first_order/final/all_${nside}_spinal_density_second_order_seed_mni.nii.gz \
            -r ${subject_dir}/${nsub}/tractoflow/${nsub}__t1_warped.nii.gz \
            -t ${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz \
            -t ${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat \
            -o ${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz
    done

    echo "|------------- 1.3) Extract thalamus -------------|"

    Right_Thalamus=(49)
    Left_Thalamus=(10)

    scil_labels_combine ${orig_rois_dir}/${nsub}_right_thalamus_orig.nii.gz \
        --volume_ids ${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz ${Right_Thalamus[*]} --merge_groups -f

    scil_labels_combine ${mni_rois_dir}/${nsub}_right_thalamus_mni.nii.gz \
        --volume_ids ${mni_rois_dir}/${nsub}_aparc.DKTatlas+aseg_mni.nii.gz ${Right_Thalamus[*]} --merge_groups -f

    scil_labels_combine ${orig_rois_dir}/${nsub}_left_thalamus_orig.nii.gz \
        --volume_ids ${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz ${Left_Thalamus[*]} --merge_groups -f

    scil_labels_combine ${mni_rois_dir}/${nsub}_left_thalamus_mni.nii.gz \
        --volume_ids ${mni_rois_dir}/${nsub}_aparc.DKTatlas+aseg_mni.nii.gz ${Left_Thalamus[*]} --merge_groups -f

    echo "|------------- 1.4) Tracking  -------------|"
    for nside in left right
    do
        # Tracking from Spinal bundle - npv 1000
        echo "|------------- 1.4a) Tracking from Spinal bundle - npv ${npv_from_spinal_track_long}  -------------|"
        scil_tracking_local \
            ${subject_dir}/${nsub}/tractoflow/${nsub}__fodf.nii.gz \
            ${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz \
            ${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz \
            ${orig_tracking_dir}/orig/${nsub}_${nside}_from_spinal_track_npv1000.trk \
            --npv ${npv_from_spinal_track_long} \
            -v -f ${gpu}

        # Tracking from Spinal bundle - npv 100
        echo "|------------- 1.4b) Tracking from Spinal bundle - npv ${npv_from_spinal_track_short}  -------------|"
        scil_tracking_local \
            ${subject_dir}/${nsub}/tractoflow/${nsub}__fodf.nii.gz \
            ${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz \
            ${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz \
            ${orig_tracking_dir}/orig/${nsub}_${nside}_from_spinal_track_npv100.trk \
            --npv ${npv_from_spinal_track_short} \
            -v -f ${gpu}

        # Tracking from Thalamus
        echo "|------------- 1.4c) Tracking from Thalamus - npv ${npv_from_thalamus_track}  -------------|"
        scil_tracking_local \
            ${subject_dir}/${nsub}/tractoflow/${nsub}__fodf.nii.gz \
            ${orig_rois_dir}/${nsub}_${nside}_thalamus_orig.nii.gz \
            ${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz\
            ${orig_tracking_dir}/orig/${nsub}_${nside}_from_thalamus_npv500.trk \
            --npv ${npv_from_thalamus_track} \
            -v -f ${gpu}
    done

    echo "|------------- 1.5) Register Tracking in MNI space -------------|"
    for ntracking in from_thalamus_npv500.trk from_spinal_track_npv100.trk from_spinal_track_npv1000.trk
    do
        for nside in left right
        do
            scil_tractogram_apply_transform \
                ${orig_tracking_dir}/orig/${nsub}_${nside}_${ntracking} \
                ${mni_dir}/MNI/mni_masked.nii.gz \
                ${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat \
                ${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_${ntracking} \
                --in_deformation ${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz \
                --remove_invalid \
                --reverse_operation -f
        done
    done
    echo "|------------- 1) Done -------------|"
    echo ""

    echo "|------------- 2) Generate ROIs -------------|"
    echo "|------------- 2.1) Registration VPM -------------|"

    for nside in left right
    do
        cp ${mni_dir}/MNI/Distal/${nside}/VPM.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz
        for nroi in ${mni_dir}/MNI/from_${nside}/*;
        do
            ROI_basename=$(basename $nroi)
            cp $nroi ${mni_rois_dir}/${nsub}_second_order_${ROI_basename/nii/_mni.nii}
        done
    done

    echo "|------------- 2.2) Generate masks -------------|"
    for nside in left right
    do
    	if [ "$nside" == "left" ]; then
    		contra_nside="right";
    	else
    		contra_nside="left";
    	fi

       # Generate Cutting masks for DTTT (ipsilateral)
        ## For dPSN : Remaining_CP to VPM
        scil_volume_math union ${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz \
            ${out_dir}/mni_space/tracking_first_order/final/all_${nside}_remaining_cp_density_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_mni.nii.gz \
            --data_type uint8 -f
        scil_labels_from_mask ${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_labels_mni.nii.gz \
            -f

        ## For CS : Spinal bundle to VPM
        scil_volume_math union ${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_mni.nii.gz \
            --data_type uint8 -f
        scil_labels_from_mask ${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_labels_mni.nii.gz \
            -f

        # Generate Cutting masks for DTTT (controlateral)
        ## For CS : Spinal bundle to Thalamus
        scil_volume_math union ${mni_rois_dir}/${nsub}_${contra_nside}_thalamus_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_mni.nii.gz \
            --data_type uint8 -f
        scil_labels_from_mask ${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_labels_mni.nii.gz \
            -f

        ## Generate Cutting masks for VTTT (controlateral)
        ## For OS and IS : Spinal bundle to Thalamus
        scil_volume_math union ${mni_rois_dir}/${nsub}_${contra_nside}_thalamus_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_mni.nii.gz \
            --data_type uint8 -f
        scil_labels_from_mask ${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_labels_mni.nii.gz \
            -f

        #For vPSN : Remaining_CP to VPM
        scil_volume_math union ${mni_rois_dir}/${nsub}_${contra_nside}_VPM_mni.nii.gz \
            ${out_dir}/mni_space/tracking_first_order/final/all_${nside}_remaining_cp_density_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_mni.nii.gz \
            --data_type uint8 -f
        scil_labels_from_mask ${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_mni.nii.gz \
            ${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_labels_mni.nii.gz \
            -f
    done
    echo "|------------- 2) Done -------------|"
    echo ""

    echo "|------------- 3) Generate Second-order bundles for trigeminal system -------------|"
    for nside in left right
    do
	if [ "$nside" == "left" ]; then
    		contra_nside="right";
    	else
    		contra_nside="left";
    	fi

        echo "|------------- 3.1) From ${nside} - VTTT (only controlateral) - OS/IS and vPSN -------------|"
        # OS and IS
        scil_tractogram_filter_by_roi \
            ${mni_tracking_dir_second_order}/orig/${nsub}_${contra_nside}_from_thalamus_npv500.trk \
            ${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_VTTT_Controlat_OSandIS.trk \
            --drawn_roi ${mni_rois_dir}/${nsub}_left_cerebellum_wm_mni.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_rois_dir}/${nsub}_right_cerebellum_wm_mni.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz 'any' 'include' \
            --drawn_roi ${mni_rois_dir}/${nsub}_${contra_nside}_thalamus_mni.nii.gz 'any' 'include' \
            --drawn_roi ${mni_dir}/MNI/from_${nside}/VTTT_Controlat_INC_Pons_Controlat.nii.gz 'any' 'include' \
            --drawn_roi ${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_CaudalMedulla_Controlat.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_dir}/MNI/from_${nside}/VTTT_Controlat_INC_VTT_Area.nii.gz 'any' 'include' \
            --drawn_roi ${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Pons_Ipsilat.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_dir}/MNI/cs_plaque.nii.gz 'any' 'exclude' -f

        ## vPSN
        ## VPM
        scil_tractogram_filter_by_roi \
            ${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_spinal_track_npv1000.trk \
            ${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_VTTT_Controlat_vPSN.trk \
            --drawn_roi ${mni_rois_dir}/${nsub}_left_cerebellum_wm_mni.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_rois_dir}/${nsub}_right_cerebellum_wm_mni.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz 'either_end' 'include' \
            --drawn_roi ${mni_rois_dir}/${nsub}_${contra_nside}_VPM_mni.nii.gz 'any' 'include' \
            --drawn_roi ${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_CaudalMedulla_Controlat.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_dir}/MNI/from_${nside}/VTTT_Controlat_INC_VTT_Area.nii.gz 'any' 'include' -f

        echo "|------------- 3.2) ${nside} - DTTT (controlateral) - CS -------------|"
        scil_tractogram_filter_by_roi \
            ${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_thalamus_npv500.trk \
            ${mni_tracking_dir_second_order}/filtered/${nsub}_from_${contra_nside}_DTTT_Controlat_CS.trk \
            --drawn_roi ${mni_rois_dir}/${nsub}_${contra_nside}_spinal_density_second_order_seed_mni.nii.gz 'any' 'include' \
            --drawn_roi ${mni_rois_dir}/${nsub}_${nside}_thalamus_mni.nii.gz 'any' 'include' \
            --drawn_roi ${mni_dir}/MNI/from_${contra_nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_dir}/MNI/from_${contra_nside}/DTTT_Controlat_INC_CaudalMedulla_Ipsilat.nii.gz 'any' 'include' \
            --drawn_roi ${mni_dir}/MNI/from_${contra_nside}/DTTT_Controlat_INC_Medulla_Controlat.nii.gz 'any' 'include' \
            --drawn_roi ${mni_dir}/MNI/from_${contra_nside}/DTTT_Controlat_EXC_Midbrain_Ipsilat.nii.gz 'any' 'exclude' -f

        echo "|------------- 3.3) ${nside} - DTTT (ipsilateral) - dPSN and CS -------------|"
        ## CS
        ## VPM/Thalamus
        scil_tractogram_filter_by_roi \
            ${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_spinal_track_npv100.trk \
            ${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk \
            --drawn_roi ${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz 'either_end' 'include' \
            --drawn_roi ${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz 'any' 'include' \
            --drawn_roi ${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz 'any' 'exclude' \
            --drawn_roi ${mni_dir}/MNI/midsagittal_plane.nii.gz 'any' 'exclude' -f

        ## dPSN
        ## VPM/Thalamus
        scil_tractogram_filter_by_roi \
            ${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_spinal_track_npv1000.trk \
            ${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_dPSN.trk \
            --drawn_roi ${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz 'any' 'include' \
            --drawn_roi ${out_dir}/mni_space/tracking_first_order/final/all_${nside}_remaining_cp_density_mni.nii.gz 'either_end' 'include' -f
    done
    echo "|------------- 3) Done -------------|"
    echo ""

    echo "|------------- 4) Cutting the Second-order bundles -------------|"
    for nside in left right
    do
    	if [ "$nside" == "left" ]; then
    		contra_nside="right";
    	else
    		contra_nside="left";
    	fi
        echo "|------------- 4.1) ${nside} - VTTT (only controlateral) - vPSN and OS/IS -------------|"
        #OS and IS - STEP 1
        scil_tractogram_cut_streamlines \
            ${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_VTTT_Controlat_OSandIS.trk \
            --labels ${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_labels_mni.nii.gz \
            ${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_VTTT_Controlat_OSandIS.trk \
            -f

        #vPSN
        scil_tractogram_cut_streamlines \
            ${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_VTTT_Controlat_vPSN.trk \
            --labels ${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_labels_mni.nii.gz \
            ${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_VTTT_Controlat_vPSN.trk \
            -f

        echo "|------------- 4.2) ${nside} - DTTT (controlateral) - CS -------------|"
        #CS
        scil_tractogram_cut_streamlines \
            ${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Controlat_CS.trk \
            --labels ${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_labels_mni.nii.gz \
            ${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Controlat_CS.trk \
            -f

        echo "|------------- 4.3) ${nside} - DTTT (ipsilateral) - dPSN and CS -------------|"

        #CS
        scil_tractogram_cut_streamlines \
            ${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk \
            --labels ${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_labels_mni.nii.gz \
            ${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk \
            -f

        #dPSN
        scil_tractogram_cut_streamlines \
            ${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_dPSN.trk \
            --labels ${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_labels_mni.nii.gz \
            ${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Ipsilat_dPSN.trk \
            -f
    done
    echo "|------------- 4) Done -------------|"
    echo ""

    echo "|------------- 5) Cleaning the Second-order bundles -------------|"
    for nside in left right
    do
        #CS
        scil_bundle_reject_outliers \
            ${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk \
            ${mni_tracking_dir_second_order}/final/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk \
            --alpha 0.30 \
            -f

        for nbundle in DTTT_Ipsilat_dPSN DTTT_Controlat_CS VTTT_Controlat_OSandIS VTTT_Controlat_vPSN
        do
            scil_bundle_reject_outliers \
                ${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_${nbundle}.trk \
                ${mni_tracking_dir_second_order}/final/${nsub}_from_${nside}_${nbundle}.trk \
                --alpha 0.50 \
                -f
        done
    done
    echo "|------------- 5) Done -------------|"
    echo ""

    echo "|------------- SECOND-ORDER FIBERS TRACTOGRAPHY FOR ${nsub} IS COMPLETED -------------|"
    echo ""
    echo ""
done
