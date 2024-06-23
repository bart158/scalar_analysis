#!/bin/bash

allInFiles=("bbll_sig_20_eLpR.root" "bbll_sig_40_eLpR.root" "bbll_sig_60_eLpR.root" "bbll_sig_70_eLpR.root" "bbll_sig_90_eLpR.root" "bbll_sig_100_eLpR.root" "bbll_sig_120_eLpR.root")
allOutFiles=("trees_for_training/bbll_sig_20_eLpR_new.root" "trees_for_training/bbll_sig_40_eLpR_new.root" "trees_for_training/bbll_sig_60_eLpR_new.root" "trees_for_training/bbll_sig_70_eLpR_new.root" "trees_for_training/bbll_sig_90_eLpR_new.root" "trees_for_training/bbll_sig_100_eLpR_new.root" "trees_for_training/bbll_sig_120_eLpR_new.root")
allFileNames=("trees_for_training/bbll_sig_20_eLpR_new" "trees_for_training/bbll_sig_40_eLpR_new" "trees_for_training/bbll_sig_60_eLpR_new" "trees_for_training/bbll_sig_70_eLpR_new" "trees_for_training/bbll_sig_90_eLpR_new" "trees_for_training/bbll_sig_100_eLpR_new" "trees_for_training/bbll_sig_120_eLpR_new")
iProc=(-1 8 2 4 0 1 5 3 6 7)
iPol=0
Ms=(5 6 7 8 9 10 11)
for i in "${!allInFiles[@]}"; do
    root -q -l -x "make_new_ttree.C(\""${allInFiles[$i]}"\",\""${allOutFiles[$i]}"\",\""${allFileNames[$i]}"\", -1, $iPol, "${Ms[$i]}")"
done