#!/bin/sh
# run_glm_predict_R.sh
# run glm_predict.R in parallel
# Bonhwang Koo
# 2/6/2018
# Command line: cat run_glm_predict_R.sh | parallel -j 7

Rscript --vanilla glm_predict.R "WASI_FSIQ"
Rscript --vanilla glm_predict.R "WASI_PRI_Comp"
Rscript --vanilla glm_predict.R "WASI_VCI_Comp"
Rscript --vanilla glm_predict.R "WASI_WIAT_Comp"
Rscript --vanilla glm_predict.R "WASI_Num"
Rscript --vanilla glm_predict.R "WASI_Spelling"
Rscript --vanilla glm_predict.R "WASI_Word_Reading"
