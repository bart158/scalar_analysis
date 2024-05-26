#!/bin/bash

allInFiles=("bbll_sig_80_eRpL.root" "SM_bg_eRpL/qq_bg_eRpL.root" "SM_bg_eRpL/qqll_bg_eRpL.root" "SM_bg_eRpL/qqlv_bg_eRpL.root" "SM_bg_eRpL/qqqq_bg_eRpL.root" "SM_bg_eRpL/qqtt_bg_eRpL.root" "SM_bg_eRpL/qqtv_bg_eRpL.root" "SM_bg_eRpL/qqvv_bg_eRpL.root" "SM_bg_eRpL/ttll_bg_eRpL.root" "SM_bg_eRpL/tttt_bg_eRpL.root")
allOutFiles=("trees_for_training/bbll_sig_80_eRpL_new.root" "trees_for_training/qq_bg_eRpL_new.root" "trees_for_training/qqll_bg_eRpL_new.root" "trees_for_training/qqlv_bg_eRpL_new.root" "trees_for_training/qqqq_bg_eRpL_new.root" "trees_for_training/qqtt_bg_eRpL_new.root" "trees_for_training/qqtv_bg_eRpL_new.root" "trees_for_training/qqvv_bg_eRpL_new.root" "trees_for_training/ttll_bg_eRpL_new.root" "trees_for_training/tttt_bg_eRpL_new.root")
allFileNames=("trees_for_training/bbll_sig_80_eRpL_new" "trees_for_training/qq_bg_eRpL_new" "trees_for_training/qqll_bg_eRpL_new" "trees_for_training/qqlv_bg_eRpL_new" "trees_for_training/qqqq_bg_eRpL_new" "trees_for_training/qqtt_bg_eRpL_new" "trees_for_training/qqtv_bg_eRpL_new" "trees_for_training/qqvv_bg_eRpL_new" "trees_for_training/ttll_bg_eRpL_new" "trees_for_training/tttt_bg_eRpL_new")
iProc=(-1 8 2 4 0 1 5 3 6 7)
iPol=1
for i in "${!allInFiles[@]}"; do
    root -q -l -x "make_new_ttree.C(\""${allInFiles[$i]}"\",\""${allOutFiles[$i]}"\",\""${allFileNames[$i]}"\", "${iProc[$i]}", $iPol)"
done