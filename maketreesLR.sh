#!/bin/bash

allInFiles=("bbll_sig_80_eLpR.root" "SM_bg_eLpR/qq_bg_eLpR.root" "SM_bg_eLpR/qqll_bg_eLpR.root" "SM_bg_eLpR/qqlv_bg_eLpR.root" "SM_bg_eLpR/qqqq_bg_eLpR.root" "SM_bg_eLpR/qqtt_bg_eLpR.root" "SM_bg_eLpR/qqtv_bg_eLpR.root" "SM_bg_eLpR/qqvv_bg_eLpR.root" "SM_bg_eLpR/ttll_bg_eLpR.root" "SM_bg_eLpR/tttt_bg_eLpR.root")
allOutFiles=("trees_for_training/bbll_sig_80_eLpR_new.root" "trees_for_training/qq_bg_eLpR_new.root" "trees_for_training/qqll_bg_eLpR_new.root" "trees_for_training/qqlv_bg_eLpR_new.root" "trees_for_training/qqqq_bg_eLpR_new.root" "trees_for_training/qqtt_bg_eLpR_new.root" "trees_for_training/qqtv_bg_eLpR_new.root" "trees_for_training/qqvv_bg_eLpR_new.root" "trees_for_training/ttll_bg_eLpR_new.root" "trees_for_training/tttt_bg_eLpR_new.root")
allFileNames=("trees_for_training/bbll_sig_80_eLpR_new" "trees_for_training/qq_bg_eLpR_new" "trees_for_training/qqll_bg_eLpR_new" "trees_for_training/qqlv_bg_eLpR_new" "trees_for_training/qqqq_bg_eLpR_new" "trees_for_training/qqtt_bg_eLpR_new" "trees_for_training/qqtv_bg_eLpR_new" "trees_for_training/qqvv_bg_eLpR_new" "trees_for_training/ttll_bg_eLpR_new" "trees_for_training/tttt_bg_eLpR_new")
iProc=(-1 8 2 4 0 1 5 3 6 7)
iPol=0
for i in "${!allInFiles[@]}"; do
    root -q -l -x "make_new_ttree.C(\""${allInFiles[$i]}"\",\""${allOutFiles[$i]}"\",\""${allFileNames[$i]}"\", "${iProc[$i]}", $iPol)"
done