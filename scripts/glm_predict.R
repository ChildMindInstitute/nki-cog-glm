setwd("/projects/txu/NKI_lifespan/scripts")
source('/home2/txu/lfcd/R/xt/R_00_source_all_function.R')
library("Rniftilib", lib.loc="/home2/milham/R/x86_64-pc-linux-gnu-library/3.1")
library("R.matlab", lib.loc="/home2/milham/R/x86_64-pc-linux-gnu-library/3.1")
library('dplyr')
library('ggplot2')
library('tidyr')

##
#groupDir <- '/data2/txu/projects/NKI_lifespan/dsc_2/group'
groupDir <- '/projects/txu/NKI_lifespan/group'

# read in subject 313
subjects_file = paste0("/data2/txu/projects/NKI_lifespan/scripts/subjects313.list")
subjects_list = as.character(read.table(subjects_file, header=FALSE, sep = "\n") [,])
nsubj <- length(subjects_list)

fname = paste0("/data2/txu/projects/NKI_lifespan/dsc_2/group/info/subjects313_info.mat")
mat <- R.matlab::readMat(fname)
sex = factor(mat$sex)
age = mat$age
basic <- data.frame(subject=subjects_list, age, sex)

# read in the cognitive task
cog_data <- read.csv("/projects/txu/NKI_lifespan/info/cog_asmt_tx.csv")[, -1]
names(cog_data)[1] <- "subject"
cog_data_test <- cog_data[cog_data$subject %in% subjects_list, c("subject", "age", "WASI_VCI_Comp", "WASI_PRI_Comp", "WASI_FSIQ", 
                                                                 "WIAT_Word_Reading", "WIAT_Num", "WIAT_Spelling", "WIAT_Comp")]
cog_df_list <- list(basic, cog_data_test)
cog_data_all <- Reduce(function(x, y) merge(x, y, all=TRUE, 
                                            by=c("subject", "age")), cog_df_list, accumulate=FALSE)

cog_names <- c("WASI_VCI_Comp", "WASI_PRI_Comp", "WASI_FSIQ", "WIAT_Word_Reading", "WIAT_Num", "WIAT_Spelling", "WIAT_Comp")
cog_age_normize_list <- c("WIAT_Word_Reading", "WIAT_Num", "WIAT_Spelling", "WIAT_Comp")

# Remove outliers from a column
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# 
hemis = c("lh","rh")
Hemispheres = c("L", "R")
Nvertex_32k <- 34292
nvertex_on <- c(25406, 25092)
names(nvertex_on) <- c("lh", "rh")

# First remove outliers
ncog <- length(cog_names)

# A function to extract gradient map based on list of subjects
# Parameters
# sublist: vector of subjects
# hemi: hemisphere (lh or rh)
niitogmap <- function(sublist, hemi) {
  # load mask
  print('read in brainmask 32k surface')
  fname <- paste0(groupDir, "/masks/", hemi, ".brain.NKI323.wb.32k_fs_LR.nii.gz")
  img = nifti.image.read(fname)
  Nvertex = img$dim[1]
  mask = img[1:Nvertex];
  idx = which(mask != 0); nvertex = length(idx)
  
  # load gradient residual data
  print('reading data')
  fname <- paste0(glmdir, '/y_cov_residual.', hemi, '.32k_fs_LR.nii.gz')
  img <- nifti.image.read(fname)
  Nvertex <- img$dim[1]
  gmap <- img[idx,1,1,sublist]
  rownames(gmap) <- idx
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
                   rsquare = numeric(ncolumn),
                   rsquare_adj = numeric(ncolumn))
  df$gradient_p <- 1
  df$sex_p <- 1
  y <- df_cog[,asmt]
  # lm model
  for (i in 1:ncolumn){
    if (i %% 500 == 0) {print(sprintf('->%d ', i))}
    gradient <- gmap[i,]
    dataframe <- data.frame(y, gradient, sex=df_cog$sex)
    try({fm <- lm(y ~ gradient + sex , data = dataframe)
    output <- summary(fm)
    df$gradient_t[i] <- output$coefficients['gradient', 't value']
    df$gradient_p[i] <- output$coefficients['gradient', 'Pr(>|t|)']
    df$sex_t[i] <- output$coefficients['sex1', 't value']
    df$sex_p[i] <- output$coefficients['sex1', 'Pr(>|t|)']
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
GLMtonii <- function(glmdir, outdir, variable, vertices, hemi, df) {
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

for (cog_name in cog_names) {
  # output dir
  outdir <- sprintf("%s/boundary/GLM/gradient_age_age2_gm_cov/residual_cognitive/%s", groupDir, cog_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  glmdir <- sprintf("%s/boundary/GLM/gradient_age_age2_gm_cov", groupDir)
  
  cog_data_rm_outliers <- unlist(lapply(cog_data_all[, cog_name], remove_outliers))
  cog_data_all_rm_outliers <- data.frame(index = 1:length(subjects_list), subjects_list, age, sex, cog_data_rm_outliers)
  colnames(cog_data_all_rm_outliers) <- c('index', 'subject', 'age', 'sex', cog_name)
  
  idx_sub <- which(is.na(cog_data_all_rm_outliers[,cog_name])==FALSE)
  df_cog <- cog_data_all_rm_outliers[idx_sub, ]
  
  set.seed(123)
  #df_cog$group <- sample(1:10, size = nrow(df_cog), replace = TRUE)
  df_cog$group <- sample(rep(seq_len(10), length.out = nrow(df_cog)))
  df_tmp <- data.frame(age = df_cog$age - mean(df_cog$age), y = df_cog[,cog_name])
  
  if (cog_name %in% cog_age_normize_list){
    df_cog[, cog_name] <- resid(lm(y ~ 0 + df_tmp$age, df_tmp, na.action = na.exclude))
  }
  fname <- paste0(outdir, '/data_cog.csv')
  write.csv(df_cog, fname)
  
  for (i in 1:10) {
    for (hemi in hemis) {
      outdir <- sprintf("%s/boundary/GLM/gradient_age_age2_gm_cov/residual_cognitive/%s/%s", groupDir, cog_name, i)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      
      df_cog_test <- df_cog[which(df_cog$group == i), ]
      df_cog_train <- df_cog[which(df_cog$group != i),]
      
      gmap_train <- niitogmap(df_cog_train$index, hemi)
      gmap_test <- niitogmap(df_cog_test$index, hemi)
      train <- niitoGLM(cog_name, df_cog_train$index, hemi, as.numeric(rownames(gmap_train)), df_cog_train)
      
      fname <- paste0(outdir, '/train_', hemi ,'.csv')
      write.csv(train, fname)
      
      xnames <- c('sex', 'gradient')
      stats_names <- c('stats')
      for (xname in xnames) {
        fname <- paste0(glmdir, '/y_cov_residual.', hemi, '.32k_fs_LR.nii.gz')
        img <- nifti.image.read(fname)
        Nvertex <- img$dim[1]
        
        var <- paste0(xname, '_p')
        idx <- train$idx
        
        fname <- GLMtonii(glmdir = glmdir, outdir = outdir, variable = xname, vertices = idx, hemi = hemi, df = train)
        
        img <- nifti.image.read(fname)
        train$cluster <- img[idx]
        n_cluster <- max(img[idx])
        
        if (n_cluster == 0) {
          print(paste0('No significant cluster for ', cog_name, ' ', xname, '_pvalue iteration ', i))
          next
        }
        
        gmap_clusters_train <- gmap_train[which(train$cluster != 0), ]
        clusters <- train$cluster[which(train$cluster != 0)]
        gmap_clusters_train <- gmap_clusters_train[order(clusters), ]
        cluster_order <- clusters[order(clusters)]
        
        gmap_clusters_test <- gmap_test[which(train$cluster != 0), ]
        gmap_clusters_test <- gmap_clusters_test[order(clusters), ]
        # Create a data frame called cluster avg (ncol = number of clusters, nrow = number of training subjects)
        cluster_avg_train <- matrix(data = NA, nrow = ncol(gmap_train), ncol = n_cluster)
        cluster_avg_test <- matrix(data = NA, nrow = ncol(gmap_test), ncol = n_cluster)
        
        for (j in 1:n_cluster) {
          gmap_cluster_train_n <- gmap_clusters_train[which(cluster_order == j), ]
          gmap_cluster_test_n <- gmap_clusters_test[which(cluster_order == j), ]
          cluster_avg_train[, j] <- colMeans(gmap_cluster_train_n)
          cluster_avg_test[, j] <- colMeans(gmap_cluster_test_n)
        }
        cluster_avg_train <- data.frame(cluster_avg_train, df_cog_train[, asmt])
        names(cluster_avg_train) <- c(paste0("cluster", 1:n_cluster), asmt)
        cluster_avg_test <- data.frame(cluster_avg_test, df_cog_test[, asmt])
        names(cluster_avg_test) <- c(paste0("cluster", 1:n_cluster), asmt)
        
        run_glm <- paste0("glm(", asmt, " ~ ", paste(c(names(cluster_avg_train)[1:n_cluster]), collapse = " + "), ", data = cluster_avg_train)")
        glm_result <- eval(parse(text = run_glm))
        coefficients <- glm_result$coefficients
        
        cluster_avg_test[, paste0(asmt, "_predicted")] <- predict.glm(glm_result, newdata = cluster_avg_test[, 1:n_cluster],
                                                                      type = "response")
        
        png(paste0(outdir, "/glm_", xname, "_pvalue_", hemi, ".png"))
        plot(cluster_avg_test[, paste0(asmt, "_predicted")], cluster_avg_test[, asmt], xlab = paste0("Predicted ", asmt), 
             ylab = paste0("Actual ", asmt))
        abline(a=0, b=1, col="red")
        dev.off()
        
        write.csv(cluster_avg_train, paste0(outdir, "/", xname, "_pvalue_", "cluster_avg_train_", hemi, ".csv"))
        write.csv(cluster_avg_test, paste0(outdir, "/", xname, "_pvalue_", "cluster_avg_test_", hemi, ".csv"))
      }
    }
  }
}


# prepare the cog data

# Convert NKI ROI masks
for (h in 1:2){
  hemi <- hemis[h]
  Hemisphere <- Hemispheres[h]
  
  fname = paste0(groupDir, "/masks/", hemi, ".brain.NKI323.wb.32k_fs_LR.nii.gz")
  gname = paste0(groupDir, "/masks/", Hemisphere, ".brain.NKI323.wb.32k_fs_LR.shape.gii") 
  scommand <- paste0("sh xt_nii2gii_32k_fs_LR.sh ", fname, " ", gname, " ", Hemisphere)
  system(scommand)
}

