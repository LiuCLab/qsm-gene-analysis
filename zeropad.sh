#!/bin/bash

export PATH=$PATH:/usr/local/c3d-1.1.0-Linux-x86_64/bin
for file in QSM_ORG*.img
do
   name=${file%.img}
   c3d $file -pad 0x0x29 0x0x29 0 -o ${name}_padded.nii
done

