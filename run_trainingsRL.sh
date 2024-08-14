#!/bin/bash

#allFileNames=("training_outcome_20_unpol" "training_outcome_40_unpol" "training_outcome_60_unpol" "training_outcome_70_unpol" "training_outcome_90_unpol" "training_outcome_100_unpol" "training_outcome_120_unpol" )

allFileNames=("20" "30" "40" "50" "60" "70" "80" "90" "95" "100" "110" "120")

for i in "${!allFileNames[@]}"; do
    cd "training_outcome_"${allFileNames[$i]}"_eRpL/"
    root -q -l -x "../train_bdt_qqll_v1.C(\"${allFileNames[$i]}\",1)"
    cd ..
done