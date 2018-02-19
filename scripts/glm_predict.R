args <- commandArgs(TRUE)
cog_name <- args[1]

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
glmdir <- sprintf("%s/boundary/GLM/gradient_age_age2_gm_cov", groupDir)

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

ncog <- length(cog_names)

# Convert NKI ROI masks
for (h in 1:2){
  hemi <- hemis[h]
  Hemisphere <- Hemispheres[h]
  
  fname = paste0(groupDir, "/masks/", hemi, ".brain.NKI323.wb.32k_fs_LR.nii.gz")
  gname = paste0(groupDir, "/masks/", Hemisphere, ".brain.NKI323.wb.32k_fs_LR.shape.gii") 
  scommand <- paste0("sh xt_nii2gii_32k_fs_LR.sh ", fname, " ", gname, " ", Hemisphere)
  system(scommand)
}

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
                   hemi,
                   gradient_t = numeric(ncolumn),
                   gradient_p = numeric(ncolumn),
                   sex_t = numeric(ncolumn),
                   sex_p = numeric(ncolumn),
                   age_t = numeric(ncolumn),
                   age_p = numeric(ncolumn),
                   rsquare = numeric(ncolumn),
                   rsquare_adj = numeric(ncolumn))
  df$hemi <- hemi
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
# Run second order polynomial on gradient
# Parameters
# asmt: name of cognitive assessment
# sublist: vector of subjects
# hemi: hemisphere (lh or rh)
# idx: index of vertices
# df_cog: dataframe of cognitive scores
niitoGLM_poly2 <- function(asmt, sublist, hemi, idx, df_cog) {
  print(asmt)
  
  gmap <- niitogmap(sublist, hemi)
  
  # run the model
  ncolumn <- nvertex_on[[hemi]]
  df <- data.frame(idx,
                   hemi,
                   gradient_t = numeric(ncolumn),
                   gradient_p = numeric(ncolumn),
                   gradient2_t = numeric(ncolumn),
                   gradient2_p = numeric(ncolumn),
                   sex_t = numeric(ncolumn),
                   sex_p = numeric(ncolumn),
                   age_t = numeric(ncolumn),
                   age_p = numeric(ncolumn),
                   rsquare = numeric(ncolumn),
                   rsquare_adj = numeric(ncolumn))
  df$hemi <- hemi
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
niitoGLM_poly3 <- function(asmt, sublist, hemi, idx, df_cog) {
  print(asmt)
  
  gmap <- niitogmap(sublist, hemi)
  
  # run the model
  ncolumn <- nvertex_on[[hemi]]
  df <- data.frame(idx,
                   hemi,
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
  df$hemi <- hemi
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
df_tmp <- data.frame(age = df_cog$age - mean(df_cog$age), y = df_cog[,cog_name])

# if (cog_name %in% cog_age_normize_list){
#   df_cog[, cog_name] <- resid(lm(y ~ 0 + df_tmp$age, df_tmp, na.action = na.exclude))
# }
fname <- paste0(outdir, '/data_cog.csv')
write.csv(df_cog, fname)

# Create a dataframe for test predication data
test_data_all <- data.frame(index = numeric(), actual = numeric(), predicted = numeric(), 
                         group = numeric(), stringsAsFactors = FALSE)
names(test_data_all) <- c("subject", cog_name, paste0(cog_name, "_predicted"), "group")

sig_clusters_lh <- matrix(data = NA, nrow = Nvertex_32k, ncol = 10)
sig_clusters2_lh <- matrix(data = NA, nrow = Nvertex_32k, ncol = 10)
# sig_clusters3_lh <- matrix(data = NA, nrow = Nvertex_32k, ncol = 10)
sig_clusters_rh <- matrix(data = NA, nrow = Nvertex_32k, ncol = 10)
sig_clusters2_rh <- matrix(data = NA, nrow = Nvertex_32k, ncol = 10)
# sig_clusters3_rh <- matrix(data = NA, nrow = Nvertex_32k, ncol = 10)

for (i in 1:10) {
  outdir <- sprintf("%s/boundary/GLM/gradient_age_age2_gm_cov/residual_cognitive/%s/%s", groupDir, cog_name, i)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  df_cog_test <- df_cog[which(df_cog$group == i), ] # Subset cognitive data into test group
  df_cog_train <- df_cog[which(df_cog$group != i),] # Subset cognitive data into training group
  
  # Get gradient data for training group
  gmap_train_lh <- niitogmap(df_cog_train$index, 'lh')
  gmap_train_rh <- niitogmap(df_cog_train$index, 'rh')
  # Get gradient data for test group
  gmap_test_lh <- niitogmap(df_cog_test$index, 'lh')
  gmap_test_rh <- niitogmap(df_cog_test$index, 'rh')
  
  # Create data frame of linear model results
  train_lh <- niitoGLM_poly2(cog_name, df_cog_train$index, 'lh', as.numeric(rownames(gmap_train_lh)), df_cog_train)
  train_rh <- niitoGLM_poly2(cog_name, df_cog_train$index, 'rh', as.numeric(rownames(gmap_train_rh)), df_cog_train)
  outdir <- sprintf("%s/boundary/GLM/gradient_age_age2_gm_cov/residual_cognitive/%s/%s", groupDir, cog_name, i)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  n_cluster_total <- 0
  
  
  for (hemi in hemis) {
    gmap_train <- get(paste0('gmap_train_', hemi))
    gmap_test <- get(paste0('gmap_test_', hemi))
    
    train <- get(paste0('train_', hemi))
    
    glasser <- get(paste0('glasser_', hemi))
    train$glasser <- glasser
  
    fname <- paste0(glmdir, '/y_cov_residual.', hemi, '.32k_fs_LR.nii.gz')
    img <- nifti.image.read(fname)
    Nvertex <- img$dim[1]
    
    var <- 'gradient_p'
    idx <- train$idx
    
    # Save linear model dataframe to nifti file, then add clusters to dataframe and save into excel file
    fname <- GLMtonii(glmdir = glmdir, outdir = outdir, variable = "gradient", vertices = idx, hemi = hemi, df = train)
    img <- nifti.image.read(fname)
    sig_clusters <- get(paste0("sig_clusters_", hemi))
    sig_clusters[, i] <- img[,1]
    assign(paste0("sig_clusters_", hemi), sig_clusters)
    
    train$cluster <- img[idx]
    
    fname <- GLMtonii(glmdir = glmdir, outdir = outdir, variable = "gradient2", vertices = idx, hemi = hemi, df = train)
    img <- nifti.image.read(fname)
    sig_clusters <- get(paste0("sig_clusters2_", hemi))
    sig_clusters[, i] <- img[,1]
    assign(paste0("sig_clusters2_", hemi), sig_clusters)
    train$cluster2 <- img[idx]

    # fname <- GLMtonii(glmdir = glmdir, outdir = outdir, variable = "gradient3", vertices = idx, hemi = hemi, df = train)
    # img <- nifti.image.read(fname)
    # sig_clusters <- get(paste0("sig_clusters3_", hemi))
    # sig_clusters[, i] <- img[,1]
    # assign(paste0("sig_clusters3_", hemi), sig_clusters)
    # train$cluster3 <- img[idx]
    
    n_cluster <- max(train$cluster)
    n_cluster2 <- max(train$cluster2)
    # n_cluster3 <- max(train$cluster3)
    
    fname <- paste0(outdir, '/train_', hemi ,'.csv')
    write.csv(train, fname)
    assign(paste0("train_", hemi), train)

    # Skip to next iteration of loop if there are no significant clusters
    #if (n_cluster == 0 & n_cluster2 == 0 & n_cluster3 == 0) {
    if (n_cluster == 0 & n_cluster2 == 0) {
      print(paste0('No significant cluster for ', cog_name, ' gradient_pvalue iteration ', i))
      next
    }
    
    # relabel clusters by glasser parcellation, i.e. if a cluster contains multiple 
    # glasser parcellations, then separate clusters by glasser parcellation
    cluster_idx <- 1
    train$cluster_by_glasser <- 0
    if (n_cluster > 0) {
      for (j in 1:n_cluster) {
        glasser_clusters <- unique(train$glasser[which(train$cluster == j)])
        for (k in glasser_clusters[which(glasser_clusters != 0)]) {
          if (sum(train$cluster == j & train$glasser == k) >= 7) {
            train$cluster_by_glasser[which(train$cluster == j & train$glasser == k)] <- cluster_idx
            cluster_idx <- cluster_idx + 1
          }
        }
      }
    }
    
    cluster_idx <- 1
    train$cluster_by_glasser2 <- 0
    if (n_cluster2 > 0) {
      for (j in 1:n_cluster2) {
        glasser_clusters <- unique(train$glasser[which(train$cluster2 == j)])
        for (k in glasser_clusters[which(glasser_clusters != 0)]) {
          if (sum(train$cluster2 == j & train$glasser == k) >= 7) {
            train$cluster_by_glasser2[which(train$cluster2 == j & train$glasser == k)] <- cluster_idx
            cluster_idx <- cluster_idx + 1
          }
        }
      }
    }
    
    
    # cluster_idx <- 1
    # train$cluster_by_glasser3 <- 0
    # if (n_cluster3 > 0) {
    #   for (j in 1:n_cluster3) {
    #     glasser_clusters <- unique(train$glasser[which(train$cluster3 == j)])
    #     for (k in glasser_clusters[which(glasser_clusters != 0)]) {
    #       if (sum(train$cluster3 == j & train$glasser == k) >= 7) {
    #         train$cluster_by_glasser3[which(train$cluster3 == j & train$glasser == k)] <- cluster_idx
    #         cluster_idx <- cluster_idx + 1
    #       }
    #     }
    #   }
    # }

    
    n_cluster <- max(train$cluster_by_glasser)
    n_cluster2 <- max(train$cluster_by_glasser2)
    # n_cluster3 <- max(train$cluster_by_glasser3)
    
    # For training gradient data frame, subset vertices in significant cluster
    gmap_clusters_train <- gmap_train[which(train$cluster_by_glasser != 0), ]
    clusters <- train$cluster_by_glasser[which(train$cluster_by_glasser != 0)]
    gmap_clusters2_train <- gmap_train[which(train$cluster_by_glasser2 != 0), ]^2
    clusters2 <- train$cluster_by_glasser2[which(train$cluster_by_glasser2 != 0)]
    # gmap_clusters3_train <- gmap_train[which(train$cluster_by_glasser3 != 0), ]^3
    # clusters3 <- train$cluster_by_glasser3[which(train$cluster_by_glasser3 != 0)]
    gmap_clusters_train <- gmap_clusters_train[order(clusters), ] # order rows in dataframe by cluster
    gmap_clusters2_train <- gmap_clusters2_train[order(clusters2), ]
    # gmap_clusters3_train <- gmap_clusters3_train[order(clusters3), ]
    cluster_order <- clusters[order(clusters)]
    cluster_order2 <- clusters2[order(clusters2)]
    # cluster_order3 <- clusters3[order(clusters3)]
    
    # Subset test gradient data frame by significant clusters in test data frame
    gmap_clusters_test <- gmap_test[which(train$cluster_by_glasser != 0), ]
    gmap_clusters_test <- gmap_clusters_test[order(clusters), ]
    gmap_clusters2_test <- gmap_test[which(train$cluster_by_glasser2 != 0), ]^2
    gmap_clusters2_test <- gmap_clusters2_test[order(clusters2), ]
    # gmap_clusters3_test <- gmap_test[which(train$cluster_by_glasser3 != 0), ]^3
    # gmap_clusters3_test <- gmap_clusters3_test[order(clusters3), ]
    
    # Create a data frame called cluster avg (ncol = number of clusters, nrow = number of training subjects)
    cluster_total <- n_cluster + n_cluster2# + n_cluster3
    cluster_avg_train <- matrix(data = NA, nrow = ncol(gmap_train), ncol = cluster_total)
    cluster_avg_test <- matrix(data = NA, nrow = ncol(gmap_test), ncol = cluster_total)
    
    # Calculate avg gradient values in each cluster
    j <- 1
    if (n_cluster > 0) {
      for (k in 1:n_cluster) {
        gmap_cluster_train_n <- gmap_clusters_train[which(cluster_order == k), ]
        gmap_cluster_test_n <- gmap_clusters_test[which(cluster_order == k), ]
        cluster_avg_train[, j] <- colMeans(gmap_cluster_train_n)
        cluster_avg_test[, j] <- colMeans(gmap_cluster_test_n)
        j <- j + 1
      }
    }
    if (n_cluster2 > 0) {
      for (k in 1:n_cluster2) {
        gmap_cluster_train_n <- gmap_clusters2_train[which(cluster_order2 == k), ]
        gmap_cluster_test_n <- gmap_clusters2_test[which(cluster_order2 == k), ]
        cluster_avg_train[, j] <- colMeans(gmap_cluster_train_n)
        cluster_avg_test[, j] <- colMeans(gmap_cluster_test_n)
        j <- j + 1
      }
    }
    # if (n_cluster3 >  0) {
    #   for (k in 1:n_cluster3) {
    #     gmap_cluster_train_n <- gmap_clusters3_train[which(cluster_order3 == k), ]
    #     gmap_cluster_test_n <- gmap_clusters3_test[which(cluster_order3 == k), ]
    #     cluster_avg_train[, j] <- colMeans(gmap_cluster_train_n)
    #     cluster_avg_test[, j] <- colMeans(gmap_cluster_test_n)
    #     j <- j + 1
    #   }
    # }
    # create training and test data frames with cluster averages, subject index, age, sex, and cognitive score
    col_names <- c()
    if (n_cluster > 0) {
      col_names <- c(col_names, paste0("cluster_", 1:n_cluster, "_", hemi))
    }
    if (n_cluster2 > 0) {
      col_names <- c(col_names, paste0("cluster2_", 1:n_cluster2, "_", hemi))
    }
    # if (n_cluster3 > 0) {
    #   col_names <- c(col_names, paste0("cluster3_", 1:n_cluster3, "_", hemi))
    # }
    col_names <- c(col_names, "index", "age", "sex", cog_name)
    cluster_avg_train <- data.frame(cluster_avg_train, df_cog_train[, c("index", "age", "sex", cog_name)])
    names(cluster_avg_train) <- col_names
    cluster_avg_test <- data.frame(cluster_avg_test, df_cog_test[, c("index", "age", "sex", cog_name)])
    names(cluster_avg_test) <- col_names
    assign(paste0("cluster_avg_train_", hemi), cluster_avg_train)
    assign(paste0("cluster_avg_test_", hemi), cluster_avg_test)
    
    remove(train)
    remove(cluster_avg_train)
    remove(cluster_avg_test)
    
    n_cluster_total <- n_cluster_total + cluster_total
  }
  if (n_cluster_total > 0) {
    train <- rbind(train_lh, train_rh)
    cluster_avg_train <- merge(cluster_avg_train_lh, cluster_avg_train_rh, by = c("index", "age", "sex", cog_name))
    cluster_avg_test <- merge(cluster_avg_test_lh, cluster_avg_test_rh, by = c("index", "age", "sex", cog_name))
  
    # Demean cluster averages
    cluster_avg_all <- rbind(cluster_avg_train, cluster_avg_test)
    cluster_avg_all[, 5:ncol(cluster_avg_all)] <- lapply(cluster_avg_all[, 5:ncol(cluster_avg_all)], function(x) {x - mean(x)})
    
    cluster_avg_train <- cluster_avg_all[cluster_avg_all$index %in% cluster_avg_train$index, ]
    cluster_avg_test <- cluster_avg_all[cluster_avg_all$index %in% cluster_avg_test$index, ]
      
    # Run generalized linear model for clusters, age, and sex
    
    run_glm <- paste0("glm(", cog_name, " ~ ", 
                      paste(c(names(cluster_avg_train)[5:ncol(cluster_avg_train)]), collapse = " + "),
                      " + age + sex, data = cluster_avg_train)")
    glm_result <- eval(parse(text = run_glm))
    run_lm <- paste0("lm(", cog_name, " ~ ", 
                      paste(c(names(cluster_avg_train)[5:ncol(cluster_avg_train)]), collapse = " + "),
                      " + age + sex, data = cluster_avg_train)")
    lm_result <- eval(parse(text = run_lm))
    coefficients <- glm_result$coefficients
    lm_coefficients <- lm_result$coefficients
    
    # Use glm to predict cognitive score for test dataset
    #cluster_avg_test[, paste0(cog_name, "_predicted")] <- predict.glm(glm_result, newdata = cluster_avg_test,
    #                                                              type = "response")
    cluster_avg_test[, paste0(cog_name, "_predicted")] <- predict(lm_result, newdata = cluster_avg_test,
                                                                  type = "response")
    
    # png(paste0(outdir, "/glm_", xname, "_pvalue_", hemi, ".png"))
    # plot(cluster_avg_test[, paste0(cog_name, "_predicted")], cluster_avg_test[, cog_name], xlab = paste0("Predicted ", cog_name), 
    #      ylab = paste0("Actual ", cog_name))
    # abline(a=0, b=1, col="red")
    # dev.off()
    
    cluster_avg_train$group <- i
    cluster_avg_test$group <- i
    
    # Combine prediction data with prediction data from previous iterations
    
    predict <- cluster_avg_test[, c("index", cog_name, paste0(cog_name, "_predicted"),
                                    "group")]
    test_data_all <- rbind(test_data_all, predict)
    
    # Plot prediction data=
    model <- summary(lm(predict[, 3] ~ predict[, 2]))
    p <- formatC(model$coefficients[2, 4], format = "e", digits = 2)
    r <- round(cor(predict[, 3], predict[, 2], method = "pearson"), digits = 3)
    
    min <- floor(min(predict[, 2:3])/10)*10
    max <- ceiling(max(predict[, 2:3])/10)*10
    
    ggplot(predict, aes(x = predict[, 2], y = predict[, 3])) + 
      geom_point() + 
      geom_smooth(method = 'lm', color = "#000000", fill = "#000000", alpha = 0.1) +
      labs(x = cog_name, y = paste(cog_name, " Predicted"), title = "gradient p_value") +
      lims(x = c(min, max),
           y = c(min, max)) +
      annotate("text", x = min, y = max - 5,  label = paste0("p = ", p, "\nr = ", r), hjust = 0)
    ggsave(paste0(outdir, "/glm_gradient_pvalue.png"))
  
    write.csv(cluster_avg_train, paste0(outdir, "/gradient_pvalue_", "cluster_avg_train.csv"), row.names = FALSE)
    write.csv(cluster_avg_test, paste0(outdir, "/gradient_pvalue_", "cluster_avg_test.csv"), row.names = FALSE)
  }
}


for (j in c("sig_clusters_lh", "sig_clusters_rh", "sig_clusters2_lh", "sig_clusters2_rh", "sig_clusters3_lh", "sig_clusters3_rh")) {
  hemi <- substring(j, nchar(j) - 1)
  fname <- paste0(groupDir, "/masks/", hemi, ".brain.NKI323.wb.32k_fs_LR.nii.gz")
  img <- nifti.image.read(fname)
  dim(img) <- c(Nvertex_32k, 1)
  sig_clusters <- get(j)
  sig_clusters_combined <- rowSums(sig_clusters > 0)
  img[1:Nvertex_32k] <- sig_clusters_combined
  outdir <- sprintf("%s/boundary/GLM/gradient_age_age2_gm_cov/residual_cognitive/%s", groupDir, cog_name)
  fname <- paste0(outdir, '/', j, '.32k_fs_LR.nii.gz')
  if (hemi == "lh") {
    Hemisphere <- "L"
  } else if (hemi == "rh") {
    Hemisphere <- "R"
  }
  gname <- paste0(outdir, '/', j, '.32k_fs_LR.func.gii')
  xt_R_save_nifti(img, fname)
  scommand <- paste0("sh xt_nii2gii_32k_fs_LR.sh ", fname, " ", gname, " ", Hemisphere)
  system(scommand); file.remove(fname)
}

outdir <- sprintf("%s/boundary/GLM/gradient_age_age2_gm_cov/residual_cognitive/%s", groupDir, cog_name)
# Plot prediction data for all 10 iterations combined
predict <- test_data_all
predict$age <- df_cog$age[match(predict$index, df_cog$index)]
predict$sex <- df_cog$sex[match(predict$index, df_cog$index)]

if (nrow(predict) > 0) {
  predict$age_group <- cut(predict$age, c(-40, -20, 5, 25, 45), include.lowest=TRUE) # Define age groups for subjects
  model <- summary(lm(predict[, 2] ~ predict[, 3]))
  p <- formatC(model$coefficients[2, 4], format = "e", digits = 2)
  r <- round(cor(predict[, 3], predict[, 2], method = "pearson"), digits = 3)
  
  min <- floor(min(predict[, 2:3])/10)*10
  max <- ceiling(max(predict[, 2:3])/10)*10
  
  ggplot(predict, aes(x = predict[, 2], y = predict[, 3])) + 
    geom_point(aes(colour = age_group)) + 
    geom_smooth(method = 'lm', color = "#000000", fill = "#000000", alpha = 0.1) +
    labs(x = cog_name, y = paste(cog_name, " Predicted"), color = "Age Group (Demeaned)", title = "gradient p_value all") +
    lims(x = c(min, max),
         y = c(min, max)) +
    annotate("text", x = min, y = max - 5,  label = paste0("p = ", p, "\nr = ", r), hjust = 0)
  ggsave(paste0(outdir, "/glm_gradient_pvalue_all.png"))
  write.csv(predict, paste0(outdir, "/glm_gradient_predict_all.csv"))
}

