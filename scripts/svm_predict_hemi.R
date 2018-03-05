args <- commandArgs(TRUE)
cog_name <- args[1]

setwd("/projects/txu/NKI_lifespan/scripts")
source('/home2/txu/lfcd/R/xt/R_00_source_all_function.R')
library("Rniftilib", lib.loc="/home2/milham/R/x86_64-pc-linux-gnu-library/3.1")
library("R.matlab", lib.loc="/home2/milham/R/x86_64-pc-linux-gnu-library/3.1")
library('dplyr')
library('ggplot2')
library('tidyr')
library('e1071') # for svm
library('reshape') # for melt

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

# read in ROI parcellation
fname <- "/projects/txu/NKI_lifespan/hcp/_SampleSphere_32k_400ROIs/resample_400_roi_cluster.lh.32k_fs_LR.nii.gz"
roi_lh <- nifti.image.read(fname)[surf_idx_lh, 1]
fname <- "/projects/txu/NKI_lifespan/hcp/_SampleSphere_32k_400ROIs/resample_400_roi_cluster.rh.32k_fs_LR.nii.gz"
roi_rh <- nifti.image.read(fname)[surf_idx_rh, 1]
roi_parcels_lh <- unique(glasser_lh)[order(unique(glasser_lh))]
roi_parcels_rh <- unique(glasser_rh)[order(unique(glasser_rh))]
  
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

# A function to apply GLM based on list of subjects
# Parameters
# asmt: name of cognitive assessment
# sublist: vector of subjects
# hemi: hemisphere (lh or rh)
# idx: index of vertices
# df_cog: dataframe of cognitive scores
niitoGLM <- function(asmt, sublist, hemi, idx, df_cog) {
  print(asmt)
  
  gmap <- niitogmap(sublist, hemi)
  
  # run the model
  ncolumn <- nvertex_on[[hemi]]
  df <- data.frame(idx,
                   gradient_t = numeric(ncolumn),
                   gradient_p = numeric(ncolumn),
                   sex_t = numeric(ncolumn),
                   sex_p = numeric(ncolumn),
                   age_t = numeric(ncolumn),
                   age_p = numeric(ncolumn),
                   rsquare = numeric(ncolumn),
                   rsquare_adj = numeric(ncolumn))
  df$gradient_p <- 1
  df$sex_p <- 1
  y <- df_cog[,asmt]
  # lm model
  for (i in 1:ncolumn){
    if (i %% 500 == 0) {print(sprintf('->%d ', i))}
    gradient <- gmap[i,]
    dataframe <- data.frame(y, gradient, sex=df_cog$sex, age=df_cog$age)
    try({fm <- lm(y ~ gradient + sex + age, data = dataframe)
    output <- summary(fm)
    df$gradient_t[i] <- output$coefficients['gradient', 't value']
    df$gradient_p[i] <- output$coefficients['gradient', 'Pr(>|t|)']
    df$sex_t[i] <- output$coefficients['sex1', 't value']
    df$sex_p[i] <- output$coefficients['sex1', 'Pr(>|t|)']
    df$age_t[i] <- output$coefficients['age', 't value']
    df$age_p[i] <- output$coefficients['age', 'Pr(>|t|)']
    df$rsquare[i] <- output$r.squared
    df$rsquare_adj[i] <- output$adj.r.squared
    df$lm_failed[i] <- 0
    }, silent = TRUE)
  }
  df
}

# A function to apply GLM based on list of subjects
# Run third order polynomial on gradient
# Parameters
# asmt: name of cognitive assessment
# sublist: vector of subjects
# hemi: hemisphere (lh or rh)
# idx: index of vertices
# df_cog: dataframe of cognitive scores
niitoGLM_poly <- function(asmt, sublist, hemi, idx, df_cog) {
  print(asmt)
  
  gmap <- niitogmap(sublist, hemi)
  
  # run the model
  ncolumn <- nvertex_on[[hemi]]
  df <- data.frame(idx,
                   gradient_t = numeric(ncolumn),
                   gradient_p = numeric(ncolumn),
                   gradient2_t = numeric(ncolumn),
                   gradient2_p = numeric(ncolumn),
                   gradient3_t = numeric(ncolumn),
                   gradient3_p = numeric(ncolumn),
                   sex_t = numeric(ncolumn),
                   sex_p = numeric(ncolumn),
                   age_t = numeric(ncolumn),
                   age_p = numeric(ncolumn),
                   rsquare = numeric(ncolumn),
                   rsquare_adj = numeric(ncolumn))
  df$gradient_p <- 1
  df$sex_p <- 1
  y <- df_cog[,asmt]
  # lm model
  for (i in 1:ncolumn){
    if (i %% 500 == 0) {print(sprintf('->%d ', i))}
    gradient <- gmap[i,]
    dataframe <- data.frame(y, gradient, gradient2 = gradient^2, gradient3 = gradient^3, sex=df_cog$sex, age=df_cog$age)
    try({fm <- lm(y ~ gradient + gradient2 + gradient3 + sex + age, data = dataframe)
    output <- summary(fm)
    df$gradient_t[i] <- output$coefficients['gradient', 't value']
    df$gradient_p[i] <- output$coefficients['gradient', 'Pr(>|t|)']
    df$gradient2_t[i] <- output$coefficients['gradient2', 't value']
    df$gradient2_p[i] <- output$coefficients['gradient2', 'Pr(>|t|)']
    df$gradient3_t[i] <- output$coefficients['gradient3', 't value']
    df$gradient3_p[i] <- output$coefficients['gradient3', 'Pr(>|t|)']
    df$sex_t[i] <- output$coefficients['sex1', 't value']
    df$sex_p[i] <- output$coefficients['sex1', 'Pr(>|t|)']
    df$age_t[i] <- output$coefficients['age', 't value']
    df$age_p[i] <- output$coefficients['age', 'Pr(>|t|)']
    df$rsquare[i] <- output$r.squared
    df$rsquare_adj[i] <- output$adj.r.squared
    df$lm_failed[i] <- 0
    }, silent = TRUE)
  }
  df
}

# A function to convert GLM values to a nifti file
# glmdir: GLM directory
# outdir: output directory
# variable: sex or gradient
# vertices: vector of vertex indices
# hemi: lh or rh
# asmt: name of measure
# df: gmap data frame
GLMtonii <- function(glmdir, outdir, variable = "gradient", vertices, hemi, df) {
  fname <- paste0(glmdir, '/y_cov_residual.', hemi, '.32k_fs_LR.nii.gz')
  img <- nifti.image.read(fname)
  Nvertex <- img$dim[1]
  
  var <- paste0(variable, '_p')
  idx <- vertices
  img$dim = c(Nvertex, 1);
  img[1:Nvertex] = matrix(0,Nvertex,1); img[idx] = df[,var]
  fname <- paste0(outdir, '/', variable, '_p_value.', hemi, '.32k_fs_LR.nii.gz')
  if (hemi == "lh") {
    Hemisphere <- "L"
  } else if (hemi == "rh") {
    Hemisphere <- "R"
  }
  gname <- paste0(outdir, '/', variable, '_p_value.', Hemisphere, '.32k_fs_LR.func.gii')
  xt_R_save_nifti(img, fname)
  scommand <- paste0("sh xt_nii2gii_32k_fs_LR.sh ", fname, " ", gname, " ", Hemisphere)
  system(scommand); file.remove(fname)
  
  scommand <- paste0("sh x_glm_clusters.sh ", Hemisphere, " ", variable, "_p_value ", outdir)
  print(scommand)
  system(scommand)
  
  fname <- paste0(outdir, '/clusters/', hemi, '.', variable, '_p_value_cluster.nii.gz')
  fname
}


# A function to convert gradient values into a support vector model
# gmap_df_lh: left hemisphere gradient map
# gmap_df_rh: right hemisphere gradient map
gmaptoParcelAvg <- function(gmap_df, parcels, hemi, df_cog) {
  parcel_avg <- array(data = NA, dim = c(nrow(df_cog), 0))
  
  for (i in parcels) {
    parcel_avg <- cbind(parcel_avg, colMeans(subset(gmap_df, rownames(gmap_df) == as.character(i))))
  }
  parcel_avg <- cbind(parcel_avg, df_cog[, cog_name])
  colnames(parcel_avg) <- c(paste(hemi, parcels, sep = "_"), cog_name)
  
  parcel_avg
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
df_cog$age_demeaned <- df_cog$age - mean(df_cog$age) # demean age
df_cog <- df_cog[, c("index", "subject", "age", "age_demeaned", "sex", cog_name)]


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
test_data_all <- data.frame(index = numeric(), age = numeric(), actual = numeric(), predicted_svm = numeric(), 
                            predicted_svm = numeric(), group = numeric(), hemi = character(), 
                            stringsAsFactors = FALSE)
names(test_data_all) <- c("index", "age", cog_name, paste0(cog_name, "_predicted_svm"), 
                          paste0(cog_name, "_predicted_svm"), "group", "hemi")



for (i in 1:10) {
  outdir <- sprintf("%s/boundary/SVM/roi_parcellation_lm/%s/%s", groupDir, cog_name, i)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  df_cog_plot <- df_cog
  df_cog_plot$train <- ifelse(df_cog_plot$group == i, 1, 0)
  
  # Plot age vs. IQ, and distinguish between test and training groups
  ggplot(df_cog_plot, aes(x = age, y = df_cog_plot[, cog_name], colour = factor(train))) +
    geom_point() +
    labs(x = "Age", 
         y = cog_name, 
         title = paste("Test Dataset ", i, " of 10"),
         colour = "") +
    scale_colour_discrete(labels = c("Train", "Test")) +
    lims(x = c(0, 100),
         y = c(40, 180))
  ggsave(paste0(outdir, "/", cog_name, "_age_dist.png"))

  df_cog_test <- df_cog[which(df_cog$group == i), ] # Subset cognitive data into test group
  df_cog_train <- df_cog[which(df_cog$group != i),] # Subset cognitive data into training group
  for (hemi in hemis) {
    
    # Get gradient data for training group
    gmap <- get(paste0("gmap_", hemi))
    gmap_train <- gmap[, as.character(df_cog_train$index)]
    
    # Get gradient data for test group
    gmap_test <- gmap[, as.character(df_cog_test$index)]
    
    roi_parcels <- get(paste0("roi_parcels_", hemi))
    
    roi_avg_train <- gmaptoParcelAvg(gmap_train, roi_parcels,
                                     hemi, df_cog_train)
    roi_avg_test <- gmaptoParcelAvg(gmap_test, roi_parcels,
                                    hemi, df_cog_test)
    
    roi_avg_train <- cbind(roi_avg_train, df_cog_train[, c("age_demeaned", "sex")])
    roi_avg_test <- cbind(roi_avg_test, df_cog_test[, c("age_demeaned", "sex")])
    write.csv(roi_avg_train, paste0(outdir, "/roi_avg_train_", hemi, ".csv"))
    write.csv(roi_avg_test, paste0(outdir, "/roi_avg_test_", hemi, ".csv"))
    
    results <- data.frame(df_cog_test[, c("index", "age", cog_name)], i, hemi)
    #results <- data.frame(df_cog_train[, c("index", "age", cog_name)], i)
    names(results) <- c("index", "age", cog_name, "group", "hemi")
    
    # glasser_avg_train <- gmaptoGlasserAvg(gmap_train_lh, gmap_train_rh, df_cog_train)
    # glasser_avg_test <- gmaptoGlasserAvg(gmap_test_lh, gmap_test_rh, df_cog_test)
    # write.csv(glasser_avg_train, paste0(outdir, "/glasser_avg_train.csv"))
    # write.csv(glasser_avg_test, paste0(outdir, "/glasser_avg_test.csv"))
    # tuneResult <- tune(svm, as.formula(paste0(cog_name, " ~ .")), data = roi_avg_train,
    #                    ranges = list(epsilon = seq(0,2,0.05), cost = 2^(0:9)))
    # png(paste0(outdir, "/svm_grid_search_", hemi, ".png"))
    # plot(tuneResult)
    # dev.off()
    # # svm_model <- svm(as.formula(paste0(cog_name, " ~ .")), data = glasser_avg_train, ep)
    # svm_model_tuned <- tuneResult$best.model
    # #svm_predict <- predict(svm_model, glasser_avg_test[, colnames(glasser_avg_test) != cog_name])
    # svm_predict_tuned <- predict(svm_model_tuned, roi_avg_test[, colnames(roi_avg_test) != cog_name])
    # #svm_predict_tuned <- predict(svm_model_tuned, roi_avg_train[, colnames(roi_avg_train) != cog_name])
    # results[, paste0(cog_name, "_predicted_svm")] <- svm_predict_tuned
    
    # Try running linear model
    lm_model <- lm(as.formula(paste0(cog_name, " ~ .")), data = roi_avg_train)
    lm_stats <- data.frame(summary(lm_model)$coef[, 3:4])
    lm_stats$sign <- as.numeric(lm_stats[, 2] < 0.05)
    write.csv(lm_stats, paste0(outdir, "/lm_stats_", hemi, ".csv"))
    
    lm_predict <- predict(lm_model, roi_avg_test[, colnames(roi_avg_test) != cog_name])
    #lm_predict <- predict(lm_model, roi_avg_train[, colnames(roi_avg_train) != cog_name])
    results[, paste0(cog_name, "_predicted_lm")] <- lm_predict
    
    write.csv(results, paste0(outdir, "/results_", hemi, ".csv"))
    
    #results_long <- melt(results,
    #                     measure.vars = paste0(cog_name, c("_predicted_svm", "_predicted_lm")))
    
    test_data_all <- rbind(test_data_all, results)
    
    # SVM statistics
    # error <- results[, cog_name] - results[, paste0(cog_name, "_predicted_svm")]
    # rmse_svm <- round(sqrt(mean((error)^2)), digits = 3)
    # epsilon <- svm_model_tuned$epsilon
    # cost <- svm_model_tuned$cost
    
    # Create data frame of linear model results
    model <- summary(lm(results[, cog_name] ~ results[, paste0(cog_name, "_predicted_lm")]))
    error <- results[, cog_name] - results[, paste0(cog_name, "_predicted_lm")]
    rmse_lm <- round(sqrt(mean((error)^2)), digits = 3)
    p <- formatC(model$coefficients[2, 4], format = "e", digits = 2)
    r <- round(cor(results[, cog_name],  results[, paste0(cog_name, "_predicted_lm")], method = "pearson"), digits = 3)
    
    min <- floor(min(results[, paste0(cog_name, "_predicted_lm")])/10)*10
    max <- ceiling(max(results[, paste0(cog_name, "_predicted_lm")])/10)*10
    
    ggplot(results, aes(x = results[, cog_name],
                        y = results[, paste0(cog_name, "_predicted_lm")])) +
      geom_point() +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      #geom_smooth(method = 'lm', color = "#000000", fill = "#000000", alpha = 0.1) +
      labs(x = cog_name,
           y = paste(cog_name, " Predicted"),
           title = paste("Test Dataset ", i, " of 10 (", hemi, ")")) +
      #scale_colour_discrete(labels = c("SVM", "LM")) +
      lims(x = c(min, max),
           y = c(min, max)) +
      annotate("text",
               x = min,
               y = max - 25,
               label = paste0("LM:\nRMSE = ", rmse_lm, "\np = ", p, "\nr = ", r), hjust = 0, size = 3)
    ggsave(paste0(outdir, "/", cog_name, "_", i, "_", hemi, ".png"))
  }
  predict <- test_data_all[test_data_all$group == i, ]
  
  model <- summary(lm(predict[, cog_name] ~ predict[, paste0(cog_name, "_predicted_lm")]))
  error <- predict[, cog_name] - predict[, paste0(cog_name, "_predicted_lm")]
  rmse_lm <- round(sqrt(mean((error)^2)), digits = 3)
  p <- formatC(model$coefficients[2, 4], format = "e", digits = 2)
  r <- round(cor(predict[, cog_name],  predict[, paste0(cog_name, "_predicted_lm")], method = "pearson"), digits = 3)
  
  min <- floor(min(predict[, paste0(cog_name, "_predicted_lm")])/10)*10
  max <- ceiling(max(predict[, paste0(cog_name, "_predicted_lm")])/10)*10
  
  ggplot(predict, aes(x = predict[, cog_name],
                      y = predict[, paste0(cog_name, "_predicted_lm")],
                      colour = hemi)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    #geom_smooth(method = 'lm', color = "#000000", fill = "#000000", alpha = 0.1) +
    labs(x = cog_name,
         y = paste(cog_name, " Predicted"),
         title = paste("Test Dataset ", i, " of 10 (both hemispheres)")) +
    #scale_colour_discrete(labels = c("SVM", "LM")) +
    lims(x = c(min, max),
         y = c(min, max)) +
    annotate("text",
             x = min,
             y = max - 25,
             label = paste0("LM:\nRMSE = ", rmse_lm, "\np = ", p, "\nr = ", r), hjust = 0, size = 3)
  ggsave(paste0(outdir, "/", cog_name, "_", i, "_both_hemis.png"))
  
}

outdir <- sprintf("%s/boundary/SVM/roi_parcellation_lm/%s", groupDir, cog_name)
# Plot prediction data for all 10 iterations combined
predict <- test_data_all


if (nrow(predict) < 0) {
  #predict_long <- melt(predict,
  #                     measure.vars = paste0(cog_name, c("_predicted_svm", "_predicted_lm")))
  
  # SVM statistics
  #predict_long_svm <- predict_long[which(predict_long$variable == paste0(cog_name, "_predicted_svm")), ]
  #error <- predict_long_svm[, cog_name] - predict_long_svm[, "value"]
  #rmse_svm <- round(sqrt(mean((error)^2)), digits = 3)
  
  # Create data frame of linear model results
  #predict_long_lm <- predict_long[which(predict_long$variable == paste0(cog_name, "_predicted_lm")), ]
  error <- predict[, cog_name] - predict[, paste0(cog_name, "_predicted_lm")]
  rmse_lm <- round(sqrt(mean((error)^2)), digits = 3)  
  model <- summary(lm(predict[, cog_name] ~ predict[, paste0(cog_name, "_predicted_lm")]))
  p <- formatC(model$coefficients[2, 4], format = "e", digits = 2)
  r <- round(cor(predict[, cog_name],  predict[, paste0(cog_name, "_predicted_lm")], method = "pearson"), digits = 3)
  
  
  min <- floor(min(predict[, paste0(cog_name, "_predicted_lm")])/10)*10
  max <- ceiling(max(predict[, paste0(cog_name, "_predicted_lm")])/10)*10
  
  ggplot(predict, aes(x = predict[, cog_name], 
                           y = predict[, paste0(cog_name, "_predicted_lm")], 
                           colour = hemi)) + 
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    #scale_colour_discrete(labels = c("SVM", "LM")) +
    #geom_smooth(method = 'lm', color = "#000000", fill = "#000000", alpha = 0.1) +
    labs(x = cog_name, y = paste(cog_name, " Predicted"), title = paste("All Test Datasets")) +
    lims(x = c(min, max),
         y = c(min, max)) +
    annotate("text", 
             x = min, 
             y = max - 15,  
             label = paste0("LM:\nRMSE = ", rmse_lm, "\np = ", p, "\nr = ", r), hjust = 0, size = 3)
  ggsave(paste0(outdir, "/", cog_name, "_predict_all.png"))
  write.csv(predict, paste0(outdir, "/predict_all.csv"))
}
