#!/bin/sh

#  get_anat2std_warps.sh <subject_ID>
# This script will take T1.mgz from a freesurfer recon_all output directory, register it to MNI, and output warp files for anat2mni and mni2anat

# Set dirs on lines 25 and 27
#
# This script uses fsl_anat, which performs more steps than is necessary to obtain the warp files (e.g. tissue segmentation). This is only necessary if you are in the peculiar situation of only having T1 images in freesurfer <mgz> space/orientation as fsl_anat does well with handling the reorientation and registration in that case. If this does not apply to you, feel free to just use flirt directly and save time. Just check that your warps work using the images created at the end of this script.
#
#  Created by Parker Singleton && Keith Jamison on 3/9/22.
#

set -e #exit on error
#set -x #print each command as it runs

export FREESURFER_HOME=/Applications/freesurfer/7.1.1
#export SUBJECTS_DIR=$FREESURFER_HOME/subjects
source $FREESURFER_HOME/SetUpFreeSurfer.sh

i=$1 #subject ID

echo ${i}


fsDIR= #folder containing each subject's freesurfer outputs

outDIR= #where to store outputs - ultimately a lot of the outputs from fsl_anat can be deleted if you do not need FSL tissue segmentation, etc. This could be skipped however fsl_anat does not appear currently able to do so.

    
if [[ -d ${outDIR}/${i} ]]
then echo "output dir for ${i} already exists"
else
    mkdir ${outDIR}/${i} #make subject dir in the output dir
fi

#convert t1.mgz to nifti
mri_convert ${fsDIR}/${i}/mri/T1.mgz ${outDIR}/${i}/T1_fs_orient.nii.gz


#brain registration and extraction etc (T1 processing) note: skipping segmentation appears to not be working so this takes longer than ideal
fsl_anat -i ${outDIR}/${i}/T1_fs_orient.nii.gz -o ${outDIR}/${i}/T1 --clobber --nocrop --noseg --nosubcortseg --nocleanup

#convert to warp files
convertwarp -r $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz -w ${outDIR}/${i}/T1.anat/T1_to_MNI_nonlin_coeff.nii.gz -o ${outDIR}/${i}/T1_to_MNI_nonlin_warp.nii.gz

invwarp --ref=${outDIR}/${i}/T1.anat/T1.nii.gz -w ${outDIR}/${i}/T1.anat/T1_to_MNI_nonlin_coeff.nii.gz -o ${outDIR}/${i}/MNI_to_T1_nonlin_warp.nii.gz

#use warp files to create images for sanity check
applywarp -i ${outDIR}/${i}/T1.anat/T1.nii.gz -r $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz -w ${outDIR}/${i}/T1_to_MNI_nonlin_warp.nii.gz -o ${outDIR}/${i}/T1_reg2mni.nii.gz

applywarp -i $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz -r ${outDIR}/${i}/T1.anat/T1.nii.gz -w ${outDIR}/${i}/MNI_to_T1_nonlin_warp.nii.gz -o ${outDIR}/${i}/mni_reg2T1.nii.gz

#do need to hand check these images since the T1 input to fsl_anat is wonky orientation, if there is a problem with reg it will be obvious



##### parallel this script
#shopt -s nullglob
#list_of_participants=(S*)
#shopt -u nullglob
#echo "${list_of_participants[@]}"
#
#for i in "${list_of_participants[@]}"; do echo bash get_anat2std_warps.sh $i; done | parallel -j 3
