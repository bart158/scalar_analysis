#!/bin/bash

allInFiles=("SM_bg_eLpL/qqll_bg_eLpL.root" "SM_bg_eLpL/qqlv_bg_eLpL.root" "SM_bg_eLpL/ttll_bg_eLpL.root")
allOutFiles=("trees_for_training/qqll_bg_eLpL_new.root" "trees_for_training/qqlv_bg_eLpL_new.root" "trees_for_training/ttll_bg_eLpL_new.root")
allFileNames=("trees_for_training/qqll_bg_eLpL_new" "trees_for_training/qqlv_bg_eLpL_new" "trees_for_training/ttll_bg_eLpL_new")
iProc=(2 4 6)
iPol=2
for i in "${!allInFiles[@]}"; do
    root -q -l -x "make_new_ttree.C(\""${allInFiles[$i]}"\",\""${allOutFiles[$i]}"\",\""${allFileNames[$i]}"\", "${iProc[$i]}", $iPol)"
done