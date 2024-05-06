#!/bin/bash
# Define paths
ROOT=/home/user/Documents/fMRI_Exp
FSLDIR="/home/user/fsl"

for subject in {11..12}; do

	cur_subject=$subject
	if [  "$subject" -lt 10 ]; then
		cur_subject="0${subject}"
		echo "add zero"
	fi
	echo "cur_subject is ${cur_subject}"
	
	FUNC_IMAGE=/home/user/Documents/fMRI_Exp/example_func_all_old/sub-$cur_subject/example_func.nii
	OUTPUT_DIR=$ROOT/MPRAGE_struct_s_only/sub$cur_subject
	STRUCT_IMAGE=$OUTPUT_DIR/s_sub$cur_subject.nii
	STRUCT_IMAGE_BET=$OUTPUT_DIR/s_sub${cur_subject}_brain.nii


	# Check if the brain-extracted file already exists
	if [  -f $STRUCT_IMAGE_BET ]; then
		echo "Brain-extracted file already exists: ${STRUCT_IMAGE_BET}"
	else
		# Perform brain extraction using BET
		bet $STRUCT_IMAGE  $STRUCT_IMAGE_BET
		echo "done bet" 
	fi

	# Linear Registration (FLIRT):

	flirt -in $FUNC_IMAGE -ref $STRUCT_IMAGE_BET -omat $OUTPUT_DIR/example_func2highres.mat -out $OUTPUT_DIR/example_func2highres
	echo "done flirt"
	convert_xfm -inverse -omat $OUTPUT_DIR/highres2example_func.mat $OUTPUT_DIR/example_func2highres.mat
	echo "done convert_xfm"

	#making png images: (for checking the registration)
	slicer $OUTPUT_DIR/example_func2highres $STRUCT_IMAGE_BET  -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; 
	pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $OUTPUT_DIR/example_func2highres1.png ; 
	slicer $STRUCT_IMAGE_BET  $OUTPUT_DIR/example_func2highres -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; 
	pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $OUTPUT_DIR/example_func2highres2.png ; 
	pngappend $OUTPUT_DIR/example_func2highres $OUTOUT_DIR/example_func2highres1.png - $OUTPUT_DIR/example_func2highres2.png $OUTPUT_DIR/example_func2highres.png; 
	rm -f sl?.png example_func2highres2.png
	rm example_func2highres1.png
	echo "saved png"

	flirt -in $STRUCT_IMAGE_BET -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain -out $OUTPUT_DIR/highres2standard -omat $OUTPUT_DIR/highres2standard.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear 
	echo  "done flirt structural"

	fnirt --iout=$OUTPUT_DIR/highres2standard_head --in=$STRUCT_IMAGE --aff=$OUTPUT_DIR/highres2standard.mat --cout=$OUTPUT_DIR/highres2standard_warp --iout=$OUTPUT_DIR/highres2standard --jout=$OUTPUT_DIR/highres2highres_jac --config=T1_2_MNI152_2mm --ref=$FSLDIR/data/standard/MNI152_T1_2mm  --refmask=$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask_dil --warpres=10,10,10
	echo "done fnirt"

	#Apply the Transformations:
	applywarp -i $STRUCT_IMAGE_BET -r $FSLDIR/data/standard/MNI152_T1_2mm_brain -o $OUTPUT_DIR/highres2standard -w $OUTPUT_DIR/highres2standard_warp

	echo "done applywarp"

	convert_xfm -inverse -omat $OUTPUT_DIR/standard2highres.mat $OUTPUT_DIR/highres2standard.mat

	#making png images: (for checking the registration)
	slicer $OUTPUT_DIR/highres2standard $FSLDIR/data/standard/MNI152_T1_2mm_brain -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; 
	pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $OUTPUT_DIR/highres2standard1.png ; 
	slicer $FSLDIR/data/standard/MNI152_T1_2mm_brain $OUTPUT_DIR/highres2standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; 
	pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $OUTPUT_DIR/highres2standard2.png ; 
	pngappend $OUTPUT_DIR/highres2standard1.png - $OUTPUT_DIR/highres2standard2.png $OUTPUT_DIR/highres2standard.png; 
	rm -f sl?.png $OUTPUT_DIR/highres2standard2.png
	echo "saved png"

	convert_xfm -omat $OUTPUT_DIR/example_func2standard.mat -concat $OUTPUT_DIR/highres2standard.mat $OUTPUT_DIR/example_func2highres.mat

	echo "done_convret_xfm"
	convertwarp --ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain --premat=$OUTPUT_DIR/example_func2highres.mat --warp1=$OUTPUT_DIR/highres2standard_warp --out=$OUTPUT_DIR/example_func2standard_warp

	applywarp --ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain  --in=$FUNC_IMAGE  --out=$OUTPUT_DIR/example_func2standard --warp=$OUTPUT_DIR/example_func2standard_warp
	echo "done applywarp"

	convert_xfm -inverse -omat $OUTPUT_DIR/standard2example_func.mat $OUTPUT_DIR/example_func2standard.mat

	#making png images: (for checking the registration)
	slicer $OUTPUT_DIR/example_func2standard $FSLDIR/data/standard/MNI152_T1_2mm_brain -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; 
	pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $OUTPUT_DIR/example_func2standard1.png ; 
	slicer $FSLDIR/data/standard/MNI152_T1_2mm_brain  $OUTPUT_DIR/example_func2standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; 
	pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $OUTPUT_DIR/example_func2standard2.png ; 
	pngappend $OUTPUT_DIR/example_func2standard1.png - $OUTPUT_DIR/example_func2standard2.png example_func2standard.png; 
	rm -f sl?.png %OUTPUT_DIR/example_func2standard2.png
done