#!/bin/bash

#allFileNames=("training_outcome_20_unpol" "training_outcome_40_unpol" "training_outcome_60_unpol" "training_outcome_70_unpol" "training_outcome_90_unpol" "training_outcome_100_unpol" "training_outcome_120_unpol" )

allFileNames=("20" "30" "40" "50" "60" "70" "80" "90" "95" "100" "110" "120")

for i in "${!allFileNames[@]}"; do
    root -q -l -x "../calc_lim_v2_with_bmst.C(\"${allFileNames[$i]}\", 0)"
    root -q -l -x "../calc_lim_v2_with_bmst.C(\"${allFileNames[$i]}\", 1)"
    root -q -l -x "../calc_lim_v2_with_bmst.C(\"${allFileNames[$i]}\", 2)"
    root -q -l -x "../calc_lim_v2_with_bmst.C(\"${allFileNames[$i]}\", 3)"
    root -q -l -x "../calc_lim_v2_with_bmst.C(\"${allFileNames[$i]}\", 4)"
    root -q -l -x "../calc_lim_v2_with_bmst.C(\"${allFileNames[$i]}\", 5)"
done