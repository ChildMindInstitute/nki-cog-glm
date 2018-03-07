#!/bin/sh
# run_lm_svm_predict_R.sh
# run lm_svm_predict.R in parallel
# Bonhwang Koo
# 3/5/2018
# Command line: cat run_lm_svm_predict_R.sh | parallel -j 7

Rscript --vanilla lm_svm_predict.R "WASI_FSIQ"
Rscript --vanilla lm_svm_predict.R "WASI_PRI_Comp"
Rscript --vanilla lm_svm_predict.R "WASI_VCI_Comp"
Rscript --vanilla lm_svm_predict.R "WIAT_Comp"
Rscript --vanilla lm_svm_predict.R "WIAT_Num"
Rscript --vanilla lm_svm_predict.R "WIAT_Spelling"
Rscript --vanilla lm_svm_predict.R "WIAT_Word_Reading"
