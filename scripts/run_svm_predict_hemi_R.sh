#!/bin/sh
# run_svm_predict_hemi_R.sh
# run svm_predict_hemi.R in parallel
# Bonhwang Koo
# 3/2/2018
# Command line: cat run_svm_predict_R.sh | parallel -j 7

Rscript --vanilla svm_predict_hemi.R "WASI_FSIQ"
Rscript --vanilla svm_predict_hemi.R "WASI_PRI_Comp"
Rscript --vanilla svm_predict_hemi.R "WASI_VCI_Comp"
Rscript --vanilla svm_predict_hemi.R "WIAT_Comp"
Rscript --vanilla svm_predict_hemi.R "WIAT_Num"
Rscript --vanilla svm_predict_hemi.R "WIAT_Spelling"
Rscript --vanilla svm_predict_hemi.R "WIAT_Word_Reading"
