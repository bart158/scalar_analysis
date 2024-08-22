#!/bin/bash

subdir="SM_bg_eRpL/"
for file in ${subdir}*.root
do
    extracted=${file/$subdir/}
    root -q -x -l "get_preselec_eff.cc(\""$file"\",\"trees_for_training/"${extracted/".root"/"_new.root"}"\")"
    
done