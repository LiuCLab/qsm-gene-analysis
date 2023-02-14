#!/bin/bash

export ANTSPATH=/opt/ANTs/bin/
export PATH=${ANTSPATH}:$PATH
for file in LPI_BB_QSM_ORG_*padded.nii
do
   name=${file%padded.nii}
   antsRegistrationSyNQuick.sh -d 3 -f 40_50_BB_LPI.nii -m $file -o ${name#LPI_BB_QSM_}4050
done
