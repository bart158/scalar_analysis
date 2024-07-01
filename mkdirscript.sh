#!/bin/bash

allFileNames=("training_outcome_20_comb" "training_outcome_30_comb" "training_outcome_40_comb" "training_outcome_50_comb" "training_outcome_60_comb" "training_outcome_70_comb" "training_outcome_80_comb" "training_outcome_90_comb" "training_outcome_95_comb" "training_outcome_100_comb" "training_outcome_110_comb" "training_outcome_120_comb")

for i in "${!allFileNames[@]}"; do
    mkdir ${allFileNames[$i]}
done