#!/bin/bash

allFileNames=("training_outcome_20_unpol" "training_outcome_40_unpol" "training_outcome_60_unpol" "training_outcome_70_unpol" "training_outcome_90_unpol" "training_outcome_100_unpol" "training_outcome_120_unpol" )

for i in "${!allFileNames[@]}"; do
    cd 
    root -q -l -x ${allFileNames[$i]}"/train_bdt_qqll.C"
done