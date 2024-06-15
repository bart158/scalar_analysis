#!/bin/bash

allInFiles=("SM_bg_eRpR/qqll_bg_eRpR.root" "SM_bg_eRpR/qqlv_bg_eRpR.root" "SM_bg_eRpR/ttll_bg_eRpR.root")
allOutFiles=("trees_for_training/qqll_bg_eRpR_new.root" "trees_for_training/qqlv_bg_eRpR_new.root" "trees_for_training/ttll_bg_eRpR_new.root")
allFileNames=("trees_for_training/qqll_bg_eRpR_new" "trees_for_training/qqlv_bg_eRpR_new" "trees_for_training/ttll_bg_eRpR_new")
iProc=(2 4 6)
iPol=3
for i in "${!allInFiles[@]}"; do
    root -q -l -x "make_new_ttree.C(\""${allInFiles[$i]}"\",\""${allOutFiles[$i]}"\",\""${allFileNames[$i]}"\", "${iProc[$i]}", $iPol)"
done