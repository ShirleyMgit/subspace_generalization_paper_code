#!/bin/sh
# create roi mask
name_mask=LOC #harvardoxford_cortical_prob_Temporal_Fusiform_Cortex_anterior_division
roi=LOC #TemporalFusiform #
name_results=tStat_visual_sameStructSameStimMinusSameStructDiffStim # tStat_hexSameMinusSameStim #tStat_hexStrMinusNothing # tStat_projSameStr_allOthers #tStat_hexStrMinusNothing #tStat_hexSameMinusSameStim #tStat_projSameStr_allOthers #tStat_hexOnHexMinusHexOnCluster

threshold=2

ROOT=/mnt/c/Users/user/Documents/fMRI_EXP
roi_file=$ROOT/roi_masks_fsl/${name_mask}_mask.nii.gz 
brain_side_mask=$ROOT/roi_masks_fsl/left_hemisphere_mask
analysis_dir_output=$ROOT/Alon/subSpaceGener/fsl_normalization/groupStats
blob_mask_dir=$analysis_dir_output/masks
echo blob_mask_dir: $blob_mask_dir
# Check if the directory exists
if [  ! -d "$blob_mask_dir"  ]; then
	# If not, create the directory
	mkdir -p $blob_mask_dir
	echo "${blob_mask_dir} Directory created successfully."
else
	echo "${blob_mask_dir}  Directory already exists."
fi

result_file=$analysis_dir_output/$name_results.nii
threshold_results=$blob_mask_dir/${name_results}_T$threshold
blob_roi=$blob_mask_dir/${name_results}_$roi.nii
new_mask=$blob_mask_dir/${name_results}T${threshold}_${roi}_mask
#new_mask=$blob_mask_dir/${name_results}T2p5_${roi}_mask

echo result file: $result_file
echo roi file: $roi_file


fslmaths $result_file -mul $roi_file $blob_roi
fslmaths $blob_roi -nan -thr $threshold $threshold_results
fslmaths $threshold_results -bin $new_mask
fslmaths $new_mask -mul $brain_side_mask ${new_mask}_L

