#!/bin/bash

allInFiles=("bbll_sig_20_eRpL.root" "bbll_sig_40_eRpL.root" "bbll_sig_60_eRpL.root" "bbll_sig_70_eRpL.root" "bbll_sig_90_eRpL.root" "bbll_sig_100_eRpL.root" "bbll_sig_120_eRpL.root")
allOutFiles=("trees_for_training/bbll_sig_20_eRpL_new.root" "trees_for_training/bbll_sig_40_eRpL_new.root" "trees_for_training/bbll_sig_60_eRpL_new.root" "trees_for_training/bbll_sig_70_eRpL_new.root" "trees_for_training/bbll_sig_90_eRpL_new.root" "trees_for_training/bbll_sig_100_eRpL_new.root" "trees_for_training/bbll_sig_120_eRpL_new.root")
allFileNames=("trees_for_training/bbll_sig_20_eRpL_new" "trees_for_training/bbll_sig_40_eRpL_new" "trees_for_training/bbll_sig_60_eRpL_new" "trees_for_training/bbll_sig_70_eRpL_new" "trees_for_training/bbll_sig_90_eRpL_new" "trees_for_training/bbll_sig_100_eRpL_new" "trees_for_training/bbll_sig_120_eRpL_new")
iProc=(-1 8 2 4 0 1 5 3 6 7)
iPol=1
Ms=(5 6 7 8 9 10 11)
for i in "${!allInFiles[@]}"; do
    root -q -l -x "make_new_ttree.C(\""${allInFiles[$i]}"\",\""${allOutFiles[$i]}"\",\""${allFileNames[$i]}"\", -1, $iPol, "${Ms[$i]}")"
done