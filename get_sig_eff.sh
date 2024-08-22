#!/bin/bash

for file in bbll*.root
do
    root -q -x -l "get_preselec_eff.cc(\""$file"\",\"trees_for_training/"${file/".root"/"_new.root"}"\")"
    
done