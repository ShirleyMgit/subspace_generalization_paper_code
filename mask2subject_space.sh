#!/bin/sh
threshold=2

root=/mnt/c/Users/user/Documents/fMRI_EXP
analysis_dir_output=$root/Alon/subSpaceGener/fsl_normalization/groupStats
masks_dir=$analysis_dir_output/masks
mask_name=tStat_visual_sameStructSameStimMinusSameStructDiffStimT2_LOC_mask_L #tStat_projSameStr_allOthersT${threshold}_Subcallosal_mask
#mask_name=tStat_projSameStr_allOthersT2p5_Subcallosal_mask
mask_in_standard=$masks_dir/$mask_name

for subject_num in 1 2 3 4 5 6 7 8 9;do #10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
	reg_folder=$root/Data4Feat/all_reg_new_names/sub-0$subject_num/reg
	echo mask in standard: $mask_in_standard

	sub_mask_dir=$masks_dir/masks_in_sub_numbers/sub-0$subject_num
	# create subject directory, if it does not exists:
	if [ ! -d $sub_mask_dir ]; then
		mkdir -p $sub_mask_dir
		echo $sub_mask_dir was created
	else
		echo $sub_mask_dir allready exists
	fi

	mask_sub=$sub_mask_dir/$mask_name
	echo mask_sub: $mask_sub
	applywarp -i $mask_in_standard -r $reg_folder/example_func -o $mask_sub -w $reg_folder/standard2highres_warp --postmat=$reg_folder/highres2example_func.mat
	gunzip $mask_sub.nii
done
