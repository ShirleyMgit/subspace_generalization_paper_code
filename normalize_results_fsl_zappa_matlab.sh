#!/bin/sh
# warp fresults from subject funtional space into standrad space
output_dir=$1
inital_ana_name=$2
vnapr=(18 19 20 21 22 24 25 26 28 29 30 31 32 34 35 36 37 38 40 42 43 44 45 46 48 49 50 51)

# Define paths
root=/data/holly-host/smark/fmri_sub_preproc_dir
#analysis_dir=$root/PileNregNativecleaned_MotionCSFonly_newMask/searchLightAllregTim
#analysis_dir_output=$root/fsl_normalized_results/searchLightAllregTim
analysis_dir=$root/$output_dir
analysis_dir_output=$analysis_dir/fsl_norm_new_names

# Check if the directory exists
if [  ! -d "$analysis_dir_output"  ]; then
	# If not, create the directory
	mkdir -p $analysis_dir_output
	echo "${analysis_dir_output} Directory created successfully."
else
	echo "${analysis_dir_output}  Directory already exists."
fi

new_reg_dir=$root/Data4Feat/all_reg_new_names

#num_files=16 #number of files to transform

#inital_ana_name=L100_FSLcleaned280920b #the inital name of the files
brain_mask=${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask

for subjectTag in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
do
	subject_num=$(($subjectTag - 1))
	echo Subject num: $subjectTag , old number ${vnapr[$subject_num]}

	if [ $subjectTag -lt 10 ]; then
		new_reg_dir_sub=$new_reg_dir/sub-0$subjectTag/reg
		analysis_output_dir_sub=$analysis_dir_output/sub-0$subjectTag
	else
		new_reg_dir_sub=$new_reg_dir/sub-$subjectTag/reg
		analysis_output_dir_sub=$analysis_dir_output/sub-$subjectTag
	fi
	
	analysis_dir_sub=$analysis_dir/sub${vnapr[$subject_num]}
	echo analysis_output_dir_sub: $analysis_output_dir_sub

	mkdir -p $analysis_output_dir_sub

	nonlinear_transf=$new_reg_dir_sub/example_func2standard_warp

	for file_num in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
	do 
		echo file number: $file_num
		my_functional=$analysis_dir_sub/${inital_ana_name}${file_num}
		echo my_functional: $my_functional
		my_warped_functional=$analysis_output_dir_sub/standard_${inital_ana_name}$file_num
		my_warped_functional_brain=${my_warped_functional}__brain
		applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_2mm --in=$my_functional --out=$my_warped_functional --warp=$nonlinear_transf 
		#mask with the brain image:
		fslmaths $my_warped_functional -mul $brain_mask $my_warped_functional_brain
		echo done normalization
		# gunzip file for SPM:
		gunzip $my_warped_functional_brain.nii
		#echo start smoothing
		#matlab_command="smooth_fsl_norm_file_name('${my_warped_functional_brain}')"
		#matlab -r -nojvm -nodesktop -nosplash -nodisplay -singleCompThread "smooth_fsl_norm_file_name(${my_warped_functional_brain})"
		#matlab -r -nojvm -nodesktop -nosplash -nodisplay -singleCompThread "smooth_fsl_norm_file_name('hi')"
	done

done