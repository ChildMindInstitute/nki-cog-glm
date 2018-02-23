args <- commandArgs(TRUE)
cog_name <- args[1]

setwd("/projects/txu/NKI_lifespan/scripts")
source('/home2/txu/lfcd/R/xt/R_00_source_all_function.R')
library("Rniftilib", lib.loc="/home2/milham/R/x86_64-pc-linux-gnu-library/3.1")
library("R.matlab", lib.loc="/home2/milham/R/x86_64-pc-linux-gnu-library/3.1")
library('dplyr')
library('ggplot2')
library('tidyr')
library('e1071')

##
#groupDir <- '/data2/txu/projects/NKI_lifespan/dsc_2/group'
groupDir <- '/projects/txu/NKI_lifespan/group'
svmdir <- sprintf("%s/boundary/SVM/gradient_age_age2_gm_cov", groupDir)

# read in subject 313
subjects_file = paste0("/data2/txu/projects/NKI_lifespan/scripts/subjects313.list")
subjects_list = as.character(read.table(subjects_file, header=FALSE, sep = "\n") [,])
nsubj <- length(subjects_list)

fname = paste0("/data2/txu/projects/NKI_lifespan/dsc_2/group/info/subjects313_info.mat")
mat <- R.matlab::readMat(fname)
sex <- factor(mat$sex)
age <- mat$age
basic <- data.frame(subject=subjects_list, age, sex)

#
hemis = c("lh","rh")
Hemispheres = c("L", "R")
Nvertex_32k <- 32492 # Total number of vertices in 32k mask
nvertex_on <- c(25406, 25092) # Number of active vertices in each hemisphere
names(nvertex_on) <- c("lh", "rh")

# read in NKI 323 surface
fname <- paste0(groupDir, "/masks/lh.brain.NKI323.wb.32k_fs_LR.nii.gz")
img <- nifti.image.read(fname)[, 1]
surf_idx_lh <- which(img != 0)
fname <- paste0(groupDir, "/masks/rh.brain.NKI323.wb.32k_fs_LR.nii.gz")
img <- nifti.image.read(fname)[, 1]
surf_idx_rh <- which(img != 0)

# read in Glasser parcellation
fname <- "/home/data/Projects/txu/templates/Glasser_et_al_2016_HCP_MMP1.0_RVVG/xt/32k_fs_LR/HCP_MMP_P210.lh.CorticalAreas_dil_Colors.32k_fs_LR.nii.gz"
glasser_lh <- nifti.image.read(fname)[surf_idx_lh, 1]
fname <- "/home/data/Projects/txu/templates/Glasser_et_al_2016_HCP_MMP1.0_RVVG/xt/32k_fs_LR/HCP_MMP_P210.rh.CorticalAreas_dil_Colors.32k_fs_LR.nii.gz"
glasser_rh <- nifti.image.read(fname)[surf_idx_rh, 1]
glasser_parcels_lh <- unique(glasser_lh)[order(unique(glasser_lh))]
glasser_parcels_rh <- unique(glasser_rh)[order(unique(glasser_rh))]
  
# read in the cognitive task
cog_data <- read.csv("/projects/txu/NKI_lifespan/info/cog_asmt_tx.csv")[, -1]
names(cog_data)[1] <- "subject"
cog_data_test <- cog_data[cog_data$subject %in% subjects_list, c("subject", "age", "WASI_VCI_Comp", "WASI_PRI_Comp", "WASI_FSIQ", 
                                                                 "WIAT_Word_Reading", "WIAT_Num", "WIAT_Spelling", "WIAT_Comp")]
cog_df_list <- list(basic, cog_data_test)
cog_data_all <- Reduce(function(x, y) merge(x, y, all=TRUE, 
                                            by=c("subject", "age")), cog_df_list, accumulate=FALSE)


# Remove outliers from a column
# Parameters
# x: Numerical vector
# na.rm: If TRUE, ignore NA values
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# A function to extract gradient map based on list of subjects
# Parameters
# sublist: vector of subjects
# hemi: hemisphere (lh or rh)
# raw: if TRUE, use raw gradient map. if FALSE (default), use age-regressed gradient map
niitogmap <- function(sublist, hemi, raw = FALSE, demeaned = FALSE) {
  # load mask
  print('read in brainmask 32k surface')
  fname <- paste0(groupDir, "/masks/", hemi, ".brain.NKI323.wb.32k_fs_LR.nii.gz")
  img = nifti.image.read(fname)
  Nvertex = img$dim[1]
  mask = img[1:Nvertex];
  idx = which(mask != 0); nvertex = length(idx)
  
  # load gradient residual data
  print('reading data')
  if (raw) {
    fname <- paste0(groupDir, '/boundary/GLM/gradient/y_4D_gradient.', hemi, '.32k_fs_LR.nii.gz')
  } else {
    fname <- paste0(glmdir, '/y_cov_residual.', hemi, '.32k_fs_LR.nii.gz')
  }
  img <- nifti.image.read(fname)
  Nvertex <- img$dim[1]
  gmap <- img[idx,1,1,sublist]
  rownames(gmap) <- idx
  if (demeaned) {
    gmap <- apply(gmap, 2, function(x) {x - mean(x)})
  }
  gmap
}

# A function to convert gradient values into a support vector model
# gmap_df_lh: left hemisphere gradient map
# gmap_df_rh: right hemisphere gradient map
gmaptoGlasserAvg <- function(gmap_df_lh, gmap_df_rh, df_cog) {
  glasser_parcel_avg_lh <- array(data = NA, dim = c(nrow(df_cog), 0))
  glasser_parcel_avg_rh <- array(data = NA, dim = c(nrow(df_cog), 0))
  
  for (i in glasser_parcels_lh) {
    glasser_parcel_avg_lh <- cbind(glasser_parcel_avg_lh, colMeans(subset(gmap_df_lh, rownames(gmap_df_lh) == as.character(i))))
  }
  colnames(glasser_parcel_avg_lh) <- paste0("lh_", glasser_parcels_lh)
  for (i in glasser_parcels_rh) {
    glasser_parcel_avg_rh <- cbind(glasser_parcel_avg_rh, colMeans(subset(gmap_df_rh, rownames(gmap_df_rh) == as.character(i))))
  }
  colnames(glasser_parcel_avg_rh) <- paste0("rh_", glasser_parcels_rh)
  glasser_parcel_avg <- data.frame(cbind(glasser_parcel_avg_lh, glasser_parcel_avg_rh), df_cog[, cog_name])
  names(glasser_parcel_avg) <- c(paste0("lh_", glasser_parcels_lh), paste0("rh_", glasser_parcels_rh), cog_name)
  
  glasser_parcel_avg
}

# output dir
outdir <- sprintf("%s/boundary/GLM/gradient_age_age2_gm_cov/residual_cognitive/%s", groupDir, cog_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Remove outliers from cognitive data
cog_data_rm_outliers <- unlist(lapply(cog_data_all[, cog_name], remove_outliers))
cog_data_all_rm_outliers <- data.frame(index = 1:length(subjects_list), subjects_list, age, sex, cog_data_rm_outliers)
colnames(cog_data_all_rm_outliers) <- c('index', 'subject', 'age', 'sex', cog_name)

idx_sub <- which(is.na(cog_data_all_rm_outliers[,cog_name])==FALSE) # indices of subjects with cognitive score within range
df_cog <- cog_data_all_rm_outliers[idx_sub, ]
df_cog$age <- df_cog$age - mean(df_cog$age) # demean age


set.seed(123)
df_cog$group <- sample(rep(seq_len(10), length.out = nrow(df_cog))) # Randomly split subjects into 10 groups
#gmap_lh <- niitogmap(df_cog$index, 'lh', raw = TRUE)
gmap_lh <- niitogmap(df_cog$index, 'lh', raw = TRUE, demeaned = FALSE)
colnames(gmap_lh) <- df_cog$index
rownames(gmap_lh) <- glasser_lh
#gmap_rh <- niitogmap(df_cog$index, 'rh', raw = TRUE)
gmap_rh <- niitogmap(df_cog$index, 'rh', raw = TRUE, demeaned = FALSE)
colnames(gmap_rh) <- df_cog$index
rownames(gmap_rh) <- glasser_rh

# if (cog_name %in% cog_age_normize_list){
#   df_cog[, cog_name] <- resid(lm(y ~ 0 + df_tmp$age, df_tmp, na.action = na.exclude))
# }
fname <- paste0(outdir, '/data_cog.csv')
write.csv(df_cog, fname)

# Create a dataframe for test prediction data
test_data_all <- data.frame(age = numeric(), actual = numeric(), predicted = numeric(), 
                         group = numeric(), stringsAsFactors = FALSE)
names(test_data_all) <- c("age", cog_name, paste0(cog_name, "_predicted"), "group")


for (i in 1:10) {
  outdir <- sprintf("%s/boundary/SVM/%s/%s", groupDir, cog_name, i)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  df_cog_test <- df_cog[which(df_cog$group == i), ] # Subset cognitive data into test group
  df_cog_train <- df_cog[which(df_cog$group != i),] # Subset cognitive data into training group
  
  # Get gradient data for training group
  gmap_train_lh <- gmap_lh[, as.character(df_cog_train$index)]
  gmap_train_rh <- gmap_rh[, as.character(df_cog_train$index)]
  # Get gradient data for test group
  gmap_test_lh <- gmap_lh[, as.character(df_cog_test$index)]
  gmap_test_rh <- gmap_rh[, as.character(df_cog_test$index)]
  
  glasser_avg_train <- gmaptoGlasserAvg(gmap_train_lh, gmap_train_rh, df_cog_train)
  glasser_avg_test <- gmaptoGlasserAvg(gmap_test_lh, gmap_test_rh, df_cog_test)
  write.csv(glasser_avg_train, paste0(outdir, "/glasser_avg_train.csv"))
  write.csv(glasser_avg_test, paste0(outdir, "/glasser_avg_test.csv"))
  tuneResult <- tune(svm, as.formula(paste0(cog_name, " ~ .")), data = glasser_avg_train,
                     ranges = list(epsilon = seq(0,1,0.1), cost = 2^(0:9)))
  # svm_model <- svm(as.formula(paste0(cog_name, " ~ .")), data = glasser_avg_train, ep)
  svm_model_tuned <- tuneResult$best.model
  #svm_predict <- predict(svm_model, glasser_avg_test[, colnames(glasser_avg_test) != cog_name])
  svm_predict_tuned <- predict(svm_model_tuned, glasser_avg_test[, colnames(glasser_avg_test) != cog_name])
  svm_results <- data.frame(df_cog_test[, c("age", cog_name)], svm_predict_tuned, i)
  names(svm_results) <- c("age", cog_name, paste0(cog_name, "_predicted"), "group")
  
  test_data_all <- rbind(test_data_all, svm_results)
  
  # Create data frame of linear model results
  error <- svm_results[, 2] - svm_results[, 3]
  rmse <- round(sqrt(mean((error)^2)), digits = 3)
  epsilon <- svm_model_tuned$epsilon
  cost <- svm_model_tuned$cost
  
  min <- floor(min(svm_results[, 2:3])/10)*10
  max <- ceiling(max(svm_results[, 2:3])/10)*10
  
  ggplot(svm_results, aes(x = svm_results[, 3], y = svm_results[, 2])) + 
    geom_point() + 
    geom_smooth(method = 'lm', color = "#000000", fill = "#000000", alpha = 0.1) +
    labs(x = cog_name, y = paste(cog_name, " Predicted"), title = paste("SVM Test Dataset ", i, " of 10")) +
    lims(x = c(min, max),
         y = c(min, max)) +
    annotate("text", x = min, y = max - 5,  label = paste0("RMSE = ", rmse, "\nepsilon = ", epsilon, "\ncost = ", cost), hjust = 0)
  ggsave(paste0(outdir, "/svm_", cog_name, "_", i, ".png"))
}

outdir <- sprintf("%s/boundary/SVM/%s", groupDir, cog_name)
# Plot prediction data for all 10 iterations combined
predict <- test_data_all

if (nrow(predict) > 0) {
  predict$age_group <- cut(predict$age, c(-40, -20, 5, 25, 45), include.lowest=TRUE) # Define age groups for subjects
  model <- summary(lm(predict[, 2] ~ predict[, 3]))
  p <- formatC(model$coefficients[2, 4], format = "e", digits = 2)
  r <- round(cor(predict[, 3], predict[, 2], method = "pearson"), digits = 3)
  rmse <- round(sqrt(mean((model$residuals)^2)), digits = 3)
  
  min <- floor(min(predict[, 2:3])/10)*10
  max <- ceiling(max(predict[, 2:3])/10)*10
  
  ggplot(predict, aes(x = predict[, 2], y = predict[, 3])) + 
    geom_point(aes(colour = age_group)) + 
    geom_smooth(method = 'lm', color = "#000000", fill = "#000000", alpha = 0.1) +
    labs(x = cog_name, y = paste(cog_name, " Predicted"), color = "Age Group (Demeaned)", title = "SVM - All Test Datasets") +
    lims(x = c(min, max),
         y = c(min, max)) +
    annotate("text", x = min, y = max - 5,  label = paste0("p = ", p, "\nr = ", r, "\nRMSE = ", rmse), hjust = 0)
  ggsave(paste0(outdir, "/svm_gradient_predict_all.png"))
  write.csv(predict, paste0(outdir, "/svm_gradient_predict_all.csv"))
}

