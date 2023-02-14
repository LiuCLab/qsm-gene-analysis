#!/bin/bash

export ANTSPATH=/opt/ANTs/bin/
export PATH=${ANTSPATH}:$PATH
for file in LPI_BB_QSM_ORG_*padded.nii
do
   name1=${file%_padded.nii}
   name2=${name1#LPI_BB_QSM_}
   antsApplyTransforms -d 3 -i 40_50_atlas_label.nii.gz -o ${name2}_label.nii.gz -r $file -t [${name2}_40500GenericAffine.mat, 1] -t ${name2}_40501InverseWarp.nii.gz -n GenericLabel
   echo $name2
done
