# surface_mni_atlases
Create subject level surface parcellations. Default output is in MNI space, modifcation for other spaces (eg native) is possible.

Curent capabilities are the "augmented Schaefer" (Schaefer cortical + Tian subcortical) at three scales (116, 232, and 454 parcels). 
See: https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal
and: https://github.com/yetianmed/subcortex/tree/master/Group-Parcellation/3T/Subcortex-Only

As well as 5 scales of Lausanne. (Thanks Andrea Luppi of  for providing the fsaverage labels for these.)

Requirements:
python
fsl
workbench
freesurfer



Steps: 
obtain warp files if you do not already have using get_anat2std_warps.sh subject_ID
This script may be modified to suit your particular situation. As currently written it will take a T1.mgz from freesurfer outputs and create warp files however the amount of computation time needed could be reduced if you already have 1 warp file, and only need to create the inverse (near the end of the script), OR if you have T1 in a orientation and space closer to MNI - you can use fnirt directly rather than fsl_anat.

Once you have anat2mni and mni2anat warp files, create_subject_atlases.sh subject_ID atlas 
  available inputs for atlas (second arg) are: sch116 sch232 sch454 ls83 ls129 ls234 ls463 ls1015 (only handles one at a time)

Will need to set dirs at the begining of each script.

Kudos to Keith Jamison of the CoCo Lab for much of the heavy lifting.

https://www.cocolaboratory.com
