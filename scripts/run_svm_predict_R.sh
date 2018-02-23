#!/bin/sh
# run_svm_predict_R.sh
# run svm_predict.R in parallel
# Bonhwang Koo
# 2/6/2018
# Command line: cat run_svm_predict_R.sh | parallel -j 7

Rscript --vanilla svm_predict.R "WASI_FSIQ"
Rscript --vanilla svm_predict.R "WASI_PRI_Comp"
Rscript --vanilla svm_predict.R "WASI_VCI_Comp"
Rscript --vanilla svm_predict.R "WIAT_Comp"
Rscript --vanilla svm_predict.R "WIAT_Num"
Rscript --vanilla svm_predict.R "WIAT_Spelling"
Rscript --vanilla svm_predict.R "WIAT_Word_Reading"
