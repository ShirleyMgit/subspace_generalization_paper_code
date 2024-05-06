#!/bin/bash
# Define paths
ROOT=/home/user/Documents/fMRI_Exp
SUBJECT=19
if [  "$SUBJECT" -lt 10 ]; then
	SUBJECT="0${SUBJECT}"
	echo "add zero"
fi
	echo "cur_subject is ${cur_subject}"
FUNC_IMAGE=/home/user/Documents/fMRI_Exp/example_func_all_old/sub-$SUBJECT/example_func.nii
OUTPUT_DIR=$ROOT/MPRAGE_struct_s_only/sub$SUBJECT
STRUCT_IMAGE=$OUTPUT_DIR/s_sub$SUBJECT.nii
STRUCT_IMAGE_BET=$OUTPUT_DIR/s_sub${SUBJECT}_brain.nii
FSLDIR="/home/user/fsl"


echo "subject: ${SUBJECT}"

# Check if the brain-extracted file already exists
if [ -f $STRUCT_IMAGE_BET]; then
    echo "Brain-extracted file already exists: ${STRUCT_IMAGE_BET}"
else
	# Perform brain extraction using BET
	bet $STRUCT_IMAGE  $STRUCT_IMAGE_BET
	echo "done bet"
fi

# Linear Registration (FLIRT):
# a. Register the EPI image to the structural image using FLIRT.
# b. Calculate the transformation matrix (affine matrix) from EPI to structural space.

flirt -in $FUNC_IMAGE -ref $STRUCT_IMAGE_BET -omat $OUTPUT_DIR/example_func2highres.mat -out $OUTPUT_DIR/example_func2highres
echo "done flirt"

flirt -in $STRUCT_IMAGE_BET -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain -out $OUTPUT_DIR/highres2standard -omat $OUTPUT_DIR/highres2standard.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear 
echo  "done flirt structural"

fnirt --iout=$OUTPUT_DIR/highres2standard_head --in=$STRUCT_IMAGE --aff=$OUTPUT_DIR/highres2standard.mat --cout=$OUTPUT_DIR/highres2standard_warp --iout=$OUTPUT_DIR/highres2standard --jout=$OUTPUT_DIR/highres2highres_jac --config=T1_2_MNI152_2mm --ref=$FSLDIR/data/standard/MNI152_T1_2mm  --refmask=$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask_dil --warpres=10,10,10
echo "done fnirt"

#Apply the Transformations:
applywarp -i $STRUCT_IMAGE_BET -r $FSLDIR/data/standard/MNI152_T1_2mm_brain -o $OUTPUT_DIR/highres2standard -w $OUTPUT_DIR/highres2standard_warp

echo "done applywarp"

convert_xfm -inverse -omat $OUTPUT_DIR/standard2highres.mat $OUTPUT_DIR/highres2standard.mat

convert_xfm -omat $OUTPUT_DIR/example_func2standard.mat -concat $OUTPUT_DIR/highres2standard.mat $OUTPUT_DIR/example_func2highres.mat
echo "done convert_xfm"

convertwarp --ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain --premat=$OUTPUT_DIR/example_func2highres.mat --warp1=$OUTPUT_DIR/highres2standard_warp --out=$OUTPUT_DIR/example_func2standard_warp

applywarp --ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain  --in=$FUNC_IMAGE  --out=$OUTPUT_DIR/example_func2standard --warp=$OUTPUT_DIR/example_func2standard_warp
echo "done applywarp"

convert_xfm -inverse -omat $OUTPUT_DIR/standard2example_func.mat $OUTPUT_DIR/example_func2standard.mat
