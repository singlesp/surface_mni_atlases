#!/bin/bash

#  create_subject_atlases.sh <subject ID> <atlas>
# Must input subject name for 1st argument and atlas for 2nd argument. choices are sch116, sch232, sch454, ls83 ls129, ls234, or ls463
#
# This script takes mni2anat, anat2mni warp files, subject freesurfer outputs, and fsaverage surface labels for schaefer or lausanne atlases and creates those atlases for a given subject in MNI space. This is useful if you have cleaned voxelwise data (eg preprocessed functional scan) in MNI space, and want to parcellate it into different surface-based atlases. Could be modified to output the atlases in native space.
#
#  Created by Parker Singleton && Keith Jamison on 3/9/22.
#

#run these lines in terminal to loop this script.
#list_of_participants=(S*)
#echo "${list_of_participants[@]}"
#
#for i in "${list_of_participants[@]}"; do echo bash create_subject_atlases.sh $i sch116; done | parallel -j 3
#for i in "${list_of_participants[@]}"; do; for j in sch116 sch232 sch454; do echo bash create_subject_atlases.sh $i $j; done;  done | parallel -j 3

set -e #exit on error
#set -x #print each command as it runs

export FREESURFER_HOME=/Applications/freesurfer/7.1.1
#export SUBJECTS_DIR=$FREESURFER_HOME/subjects
source $FREESURFER_HOME/SetUpFreeSurfer.sh

if [[ $1 == "" ]]
then
    echo "must input subject name for 1st argument and atlas for 2nd argument. choices are sch116, sch232, sch454, ls83 ls129, ls234, or ls463"
    exit 0
fi

i=$1

if [[ $2 == "" ]]
then
    echo "must input atlas for second argument. choices are sch116, sch232, sch454, ls83 ls129, ls234, or ls463"
    exit 0
fi

cortname=$2 #input options sch116, sch232, sch454 for the augmented Schaefer/Tian cort/subcort atlases


baseDIR= #depends how your dirs are organized

funcDIR= #points to a dir containing a functional brain-extracted volume registered to MNI space for each subject (used to create a reference volume)


freesurferdir=${baseDIR}/FS_archive/${i} #folder containing /mri /surf, etc...
Subject=${i} #for output naming
isMNI=1 #1=MNI or 0=anatomical space

WD=${baseDIR}/atlases/${i} #working dir (where intermediate surface files will be placed)
OD=${WD}/${cortname} #final output dir

if [[ -f  $OD/final/"$Subject"_${cortname}.nii.gz ]]; then
    while true; do
    read -p "$2 final output for $1 already exist. Proceed anyways? (Y or N)" yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
fi

anat2mni=${baseDIR}/atlases/${Subject}/T1_to_MNI_nonlin_warp.nii.gz #FNIRT .nii.gz warp file for anat->MNI
mni2anat=${baseDIR}/atlases/${Subject}/MNI_to_T1_nonlin_warp.nii.gz #FNIRT .nii.gz warp file for mni->anat
SurfaceAtlasDIR=${baseDIR}/atlases/templates/hcp/standard_mesh_atlases #contains HCP resolution surfaces (downloaded from https://github.com/Washington-University/HCPpipelines/tree/master/global/templates/standard_mesh_atlases)

if [[ -f ${WD}/ref_vol.nii.gz ]]
then echo "ref vol for ${i} already exists"
else
    fslmaths ${funcDIR}/${i}* -Tmean ${WD}/ref_vol.nii.gz
fi

refvol=${WD}/ref_vol.nii.gz #reference volume for output space (eg: an MNI resting state volume: can use fslmaths to avg across time)

#fs average surface templates in the following section are contained in the templates/ dir of this repo
if [ ${cortname} = sch116 ]; then
    atlas_fsavg_lh=${baseDIR}/atlases/templates/sch/lh.Schaefer2018_100Parcels_7Networks_order.annot #fs average surface labels
    atlas_fsavg_rh=${baseDIR}/atlases/templates/sch/rh.Schaefer2018_100Parcels_7Networks_order.annot
    nosubname=sch100 #name of atlas file without sub-cortex
    shift_r=50 #how many values to shift the r/l hemisphere before combining
    sub_cort=${baseDIR}/atlases/templates/Tian/Tian16_sh100.nii.gz #subcortical atlas to be added in (these are MNI templates for Tian/sch // subject specific freesurfer subcortex for lausanne)
elif [ ${cortname} = sch232 ]; then    atlas_fsavg_lh=${baseDIR}/atlases/templates/sch/lh.Schaefer2018_200Parcels_7Networks_order.annot
    atlas_fsavg_rh=${baseDIR}/atlases/templates/sch/rh.Schaefer2018_200Parcels_7Networks_order.annot
    nosubname=sch200
    shift_r=100
    sub_cort=${baseDIR}/atlases/templates/Tian/Tian32_sh200.nii.gz
elif [ ${cortname} = sch454 ]; then
    atlas_fsavg_lh=${baseDIR}/atlases/templates/sch/lh.Schaefer2018_400Parcels_7Networks_order.annot
    atlas_fsavg_rh=${baseDIR}/atlases/templates/sch/rh.Schaefer2018_400Parcels_7Networks_order.annot
    nosubname=sch400
    shift_r=200
    sub_cort=${baseDIR}/atlases/templates/Tian/Tian54_sh400.nii.gz
elif [ ${cortname} = ls83 ]; then
    atlas_fsavg_lh=${baseDIR}/atlases/templates/ls/atl-Cammoun2012_space-fsaverage_res-033_hemi-L_deterministic.annot
    atlas_fsavg_rh=${baseDIR}/atlases/templates/ls/atl-Cammoun2012_space-fsaverage_res-033_hemi-R_deterministic.annot
    nosubname=${cortname}.nosub
    shift_l=41
    sub_cort=${OD}/${cortname}_subcortex.nii.gz
    rh_val=35 #max roi in the rh cortex
    lh_val=35 #max roi in the lh cortex
elif [ ${cortname} = ls129 ]; then
    atlas_fsavg_lh=${baseDIR}/atlases/templates/ls/atl-Cammoun2012_space-fsaverage_res-060_hemi-L_deterministic.annot
    atlas_fsavg_rh=${baseDIR}/atlases/templates/ls/atl-Cammoun2012_space-fsaverage_res-060_hemi-R_deterministic.annot
    nosubname=${cortname}.nosub
    shift_l=64
    sub_cort=${OD}/${cortname}_subcortex.nii.gz
    rh_val=58 #max roi in the rh cortex
    lh_val=58 #max roi in the lh cortex
elif [ ${cortname} = ls234 ]; then
    atlas_fsavg_lh=${baseDIR}/atlases/templates/ls/atl-Cammoun2012_space-fsaverage_res-125_hemi-L_deterministic.annot
    atlas_fsavg_rh=${baseDIR}/atlases/templates/ls/atl-Cammoun2012_space-fsaverage_res-125_hemi-R_deterministic.annot
    nosubname=${cortname}.nosub
    shift_l=115
    sub_cort=${OD}/${cortname}_subcortex.nii.gz
    rh_val=109 #max roi in the rh cortex
    lh_val=112 #max roi in the lh cortex
elif [ ${cortname} = ls463 ]; then
    atlas_fsavg_lh=${baseDIR}/atlases/templates/ls/atl-Cammoun2012_space-fsaverage_res-250_hemi-L_deterministic.annot
    atlas_fsavg_rh=${baseDIR}/atlases/templates/ls/atl-Cammoun2012_space-fsaverage_res-250_hemi-R_deterministic.annot
    nosubname=${cortname}.nosub
    shift_l=230
    sub_cort=${OD}/${cortname}_subcortex.nii.gz
    rh_val=224 #max roi in the rh cortex
    lh_val=226 #max roi in the lh cortex
elif [ ${cortname} = ls1015 ]; then
    atlas_fsavg_lh=${baseDIR}/atlases/templates/ls/atl-Cammoun2012_space-fsaverage_res-500_hemi-L_deterministic.annot
    atlas_fsavg_rh=${baseDIR}/atlases/templates/ls/atl-Cammoun2012_space-fsaverage_res-500_hemi-R_deterministic.annot
    nosubname=${cortname}.nosub
    shift_l=508
    sub_cort=${OD}/${cortname}_subcortex.nii.gz
    rh_val=502 #max roi in the rh cortex
    lh_val=500 #max roi in the lh cortex
fi

#executable dirs:
CARET7DIR=/Applications/workbench/bin_macosx64 #$(dirname $(which wb_command))
FSBIN=${FREESURFER_HOME}/bin
FSLDIR=${FSLDIR}

LowResMesh=32
HighResMesh=164
###########

mkdir -p $OD/final

MatrixXYZ=`mri_info --cras ${freesurferdir}/mri/wmparc.mgz`
MatrixX=`echo ${MatrixXYZ} | awk '{print $1;}'`
MatrixY=`echo ${MatrixXYZ} | awk '{print $2;}'`
MatrixZ=`echo ${MatrixXYZ} | awk '{print $3;}'`
echo "1 0 0 ${MatrixX}" >  $WD/c_ras.mat
echo "0 1 0 ${MatrixY}" >> $WD/c_ras.mat
echo "0 0 1 ${MatrixZ}" >> $WD/c_ras.mat
echo "0 0 0 1"          >> $WD/c_ras.mat

InverseAtlasTransform=${mni2anat/.nii.gz/""}
AtlasTransform=${anat2mni/.nii.gz/""}

AtlasSpaceFolder=$WD/MNINonLinear
T1wFolder=$WD/T1w
NativeFolder=Native
mkdir -p $T1wFolder
mkdir -p "$T1wFolder"/"$NativeFolder"
mkdir -p "$T1wFolder"/fsaverage_LR${LowResMesh}k
mkdir -p $AtlasSpaceFolder
mkdir -p "$AtlasSpaceFolder"/fsaverage
mkdir -p "$AtlasSpaceFolder"/"$NativeFolder"
mkdir -p "$AtlasSpaceFolder"/fsaverage_LR${LowResMesh}k

################################################################

#Map fsaverage surface atlas to subject surfaces
#for each hemisphere:
#1. convert "native" freesurfer sphere and sphere.reg to gii
#2. create 32k (~2mm) version of fsavg-aligned sphere
#3. warp gm+wm surfaces to MNI
#4. downsample anat-space and MNI-space gm+wm surfaces to 32k resolution
#5. map fsaverage atlas label annot to subject-specific 32k surface as label.gii



for Hemisphere in L R; do
    if [ "$Hemisphere" = L ]; then
        hemisphere=l
        StructureName=CORTEX_LEFT
    else
        hemisphere=r
        StructureName=CORTEX_RIGHT
    fi

    if [[ -f "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii ]]; then
        echo "using previously made surface files"
    else
        for Surface in sphere sphere.reg; do
            ${FSBIN}/mris_convert ${freesurferdir}/surf/"$hemisphere"h."$Surface" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
        done

        cp "$SurfaceAtlasDIR"/fs_"$Hemisphere"/fsaverage."$Hemisphere".sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii

        cp "$SurfaceAtlasDIR"/fs_"$Hemisphere"/fs_"$Hemisphere"-to-fs_LR_fsaverage."$Hemisphere"_LR.spherical_std."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".def_sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii

        cp "$SurfaceAtlasDIR"/fsaverage."$Hemisphere"_LR.spherical_std."$HighResMesh"k_fs_LR.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii

        ${CARET7DIR}/wb_command -surface-sphere-project-unproject "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".sphere.reg.native.surf.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".def_sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.surf.gii


        RegSphere="${AtlasSpaceFolder}/${NativeFolder}/${Subject}.${Hemisphere}.sphere.reg.reg_LR.native.surf.gii"


        #mris_convert to create MNINonLinear/Native/"$Subject"."$Hemisphere"."$Surface".native.surf.gii for sphere and sphere.reg
        #copied from pipeline template dir:
        cp "$SurfaceAtlasDIR"/"$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii

        #Warp full resolution wm+gm surfaces from anat to MNI, then downsample both anat and MNI versions to 32k resolution
        for Surface in white pial; do
            ${FSBIN}/mris_convert ${freesurferdir}/surf/"$hemisphere"h."$Surface" "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
            ${CARET7DIR}/wb_command -surface-apply-affine "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii $WD/c_ras.mat "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
            ${CARET7DIR}/wb_command -surface-apply-warpfield "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii "$InverseAtlasTransform".nii.gz "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii -fnirt "$AtlasTransform".nii.gz
        
            #need to downsample the MNI and T1w surfaces
            ${CARET7DIR}/wb_command -surface-resample "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii BARYCENTRIC "$T1wFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere"."$Surface"."$LowResMesh"k_fs_LR.surf.gii

            ${CARET7DIR}/wb_command -surface-resample "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii BARYCENTRIC "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere"."$Surface"."$LowResMesh"k_fs_LR.surf.gii


        done

    fi

    ### get label files for atlas
    if [ ${Hemisphere} = L ]; then
        labelfile_templatespace=${atlas_fsavg_lh}
    elif [ ${Hemisphere} = R ]; then
        labelfile_templatespace=${atlas_fsavg_rh}
    fi
    labelfile_gii=$WD/${hemisphere}h.${cortname}.label.gii
    labelfile_out=$WD/${hemisphere}h.${cortname}_${Subject}.label.gii
    labelfile_out32=$WD/${Hemisphere}.${cortname}_${Subject}_32k.label.gii
    ${FSBIN}/mris_convert --annot ${labelfile_templatespace} ${FREESURFER_HOME}/subjects/fsaverage/surf/${hemisphere}h.white ${labelfile_gii}

    ${CARET7DIR}/wb_command -label-resample ${labelfile_gii} "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".sphere.reg.native.surf.gii BARYCENTRIC ${labelfile_out} -largest
    ${CARET7DIR}/wb_command -set-structure ${labelfile_out} $StructureName

    ${CARET7DIR}/wb_command -label-resample ${labelfile_out} "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.surf.gii  "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii BARYCENTRIC ${labelfile_out32} -largest
    ${CARET7DIR}/wb_command -set-structure ${labelfile_out32} $StructureName
done


############################################
# Map surface labels to voxels
outres=$($FSLDIR/bin/fslinfo $refvol |  grep -Em1 '^pixdim1' | awk '{printf "%g\n",$NF}')
if [ "$isMNI" = 1 ]; then
	res_suffix=".mni${outres}mm"
	anatdir="${AtlasSpaceFolder}"
else
	res_suffix=".${outres}mm"
	anatdir="${T1wFolder}"
fi

for Hemisphere in L R; do

	labelfile_out32=$WD/${Hemisphere}.${cortname}_${Subject}_32k.label.gii
	volfile_out=$OD/"$Subject"."$Hemisphere".${cortname}${res_suffix}.nii.gz
    
	wsurf="$anatdir"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".white."$LowResMesh"k_fs_LR.surf.gii
	psurf="$anatdir"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".pial."$LowResMesh"k_fs_LR.surf.gii
	
	${CARET7DIR}/wb_command -label-to-volume-mapping ${labelfile_out32} $wsurf $refvol ${volfile_out} -ribbon-constrained $wsurf $psurf
	cornerval=$(python -c 'import nibabel as nib; import numpy as np; V=nib.load("'${volfile_out}'"); Vdata=V.get_fdata(); print("%.2f" % (Vdata[0,0,0]))')
	if [ "${cornerval}" != "0.00" ]; then
		#why is this happening?
		python -c 'import nibabel as nib; import numpy as np; f="'${volfile_out}'"; V=nib.load(f); Vdata=V.get_fdata(); Vdata[Vdata==Vdata[0,0,0]]=0; nib.save(nib.Nifti1Image(Vdata,affine=V.affine,header=V.header),f);'
	fi
 
done

#if Lausanne, make subcortical ROIs from freesurfer
if [[ ${cortname} == "ls"* ]]; then

    echo "creating sub-cortical ROI's from aparc+aseg"
    mri_convert --reslice_like ${WD}/T1.anat/T1.nii.gz --resample_type nearest "${freesurferdir}/mri/aparc+aseg.mgz" ${WD}/aparc+aseg.nii.gz -odt float
    applywarp -i ${WD}/aparc+aseg.nii.gz -o ${WD}/aparc+aseg_MNI.nii.gz --ref=${refvol} --warp=${AtlasTransform} --interp=nn
    python fslselect.py ${WD}/aparc+aseg_MNI.nii.gz ${OD}/${cortname}_subcortex.nii.gz $(awk '{ printf "%s=%s\n",$1,$2 }' ${baseDIR}/atlases/templates/ls/${cortname}_subcort_labels.tsv | tr -d \")

fi

hemi_vols=($OD/"$Subject".L.${cortname}${res_suffix}.nii.gz $OD/"$Subject".R.${cortname}${res_suffix}.nii.gz)

if [[ ${cortname} == "sch"* ]]; then
    echo "combining R/L hemisphere atlases"
    fslmaths ${hemi_vols[1]} -add $shift_r -mas ${hemi_vols[1]} $OD/right_shifted.nii.gz
    fslmaths ${hemi_vols[0]} -max $OD/right_shifted.nii.gz $OD/final/"$Subject"_${nosubname}.nii.gz
    echo "adding in Tian subcortex"
    fslmaths $OD/final/"$Subject"_${nosubname}.nii.gz -max ${sub_cort} $OD/final/"$Subject"_${cortname}.nii.gz
fi

if [[ ${cortname} == "ls"* ]]; then
    echo "combining R/L hemisphere atlases and removing corpus callosum"
    
    #remove corpuscallosum roi 4 from atlas
    #make a mask of ROIs 1,2,3
    #subtract 1 from each hemi mask with hemi
    #add 1 to each hemi mask with ROIs1,2,3
    seq 1 $rh_val >> boop.txt
    python fslselect.py ${hemi_vols[1]} ${hemi_vols[1]} $(sed '4d' boop.txt | awk '{ printf "%s=%s\n",$1,NR }')
    rm boop.txt
    seq 1 $lh_val >> boop.txt
    python fslselect.py ${hemi_vols[0]} ${hemi_vols[0]} $(sed '4d' boop.txt | awk '{ printf "%s=%s\n",$1,NR }')
    rm boop.txt
    
    fslmaths ${hemi_vols[0]} -add $shift_l -mas ${hemi_vols[0]} $OD/left_shifted.nii.gz
    fslmaths ${hemi_vols[1]} -max $OD/left_shifted.nii.gz $OD/final/"$Subject"_${nosubname}.nii.gz
    
    echo "adding in subcortex"
    fslmaths $OD/final/"$Subject"_${nosubname}.nii.gz -max ${sub_cort} $OD/final/"$Subject"_${cortname}.nii.gz
fi

echo "complete."


## check multiple outputs in command line:
#from atlases dir:
#subs=(S*)
#at= #sch or ls
#for i in "${subs[@]}"; do; for j in 116 232 454; do fsleyes ${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz ${i}/${at}${j}/final/${i}_${at}${j}.nii.gz -ot label -l random_big templates/mni/*${j}* -ot label -l random_big; done; done
#quit fsleyes and loop continues to next output
