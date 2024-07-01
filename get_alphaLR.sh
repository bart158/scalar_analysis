#!/bin/bash

#allFileNames=("training_outcome_20_unpol" "training_outcome_40_unpol" "training_outcome_60_unpol" "training_outcome_70_unpol" "training_outcome_90_unpol" "training_outcome_100_unpol" "training_outcome_120_unpol" )

allFileNames=("20" "30" "40" "50" "60" "70" "80" "90" "95" "100" "110" "120")

for i in "${!allFileNames[@]}"; do
    cd "training_outcome_"${allFileNames[$i]}"_eLpR/"
    cp ../training_outcome_20_eLpR/find_limit_with_loop2.C ./find_limit_with_loop2.C
    root -q -l -x "find_limit_with_loop2.C"
    cd ..
done