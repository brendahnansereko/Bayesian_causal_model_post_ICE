
# Capture the Slurm array task ID (default to 1 if not provided)
args <- commandArgs(trailingOnly = TRUE)
batch <- ifelse(length(args) > 0, as.numeric(args[1]), 1)

# Set up random seed for reproducibility
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(69012365)  # Base seed
s <- .Random.seed
for (i in 1:batch) {
  s <- nextRNGStream(s)
}
.GlobalEnv$.Random.seed <- s

# Print batch number (for debugging)
print(paste("Running batch:", batch))

# Load necessary library
library(MASS) 
library(tidyr)
library(dplyr)
library(rbmi)
library(rstan)

rstan_options(auto_write = TRUE)
rstan_options(cxx = "/opt/ohpc/pub/compiler/gcc/9.4.0/bin/g++")

if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2", repos = "https://cran.r-project.org")
}
library(reshape2)

setwd("/home/lsh1901704/")

# Load the pre-generated datasets for this batch
# Assuming the file is named: batch_X_datasets.rds
data_list <- readRDS(paste0("/home/lsh1901704/pioneer_paper2/datasets/low/batch_", batch, "_datasets.rds"))

# Parameters
MAR <- MAR_se <- causal <- causal_se <- k_0 <- k_0_se <- mean_impute_se  <- jackknife_se <- boot_se <- mean_impute <- jackknife_se <- numeric(length(data_list))

iter <- length(data_list)

for (j in 1:iter) {
  
  # Use the data from the list
  data_discont <- data_list[[j]]
  
  # _______________________________________
  # 2.introduce monotone MCAR
  # _______________________________________
  # Introduce monotone MCAR in 10% of discontinuers in the Treatment group
  discont_treatment_ids <- unique(data_discont$id[data_discont$groupx == "Treatment" & data_discont$Discontinue == 1])
  n_mcar <- ceiling(length(discont_treatment_ids) * 0.9)
  #set.seed(123)  # Reproducible random selection
  mcar_ids <- sample(discont_treatment_ids, n_mcar)
  
  for (i in mcar_ids) {
    person_data <- data_discont[data_discont$id == i, ]
    r_time <- max(person_data$r)
    if (r_time < 5) {  # If r_time == 5, no future visits exist
      mcar_visitn <- (r_time + 1):5
      data_discont$y_mar[data_discont$id == i & data_discont$visitn %in% mcar_visitn] <- NA
    }
  }
  
  
  # Introduce monotone MCAR in 10% of discontinuers in the Control group
  discont_ctrl_ids <- unique(data_discont$id[data_discont$groupx == "Control" & data_discont$Discontinue == 1])
  n_mcar <- ceiling(length(discont_ctrl_ids) * 0.9)
  #set.seed(123)  # Reproducible random selection
  mcar_ctrl_ids <- sample(discont_ctrl_ids, n_mcar)
  
  for (i in mcar_ctrl_ids) {
    person_data <- data_discont[data_discont$id == i, ]
    r_time <- max(person_data$r)
    if (r_time < 5) {  # If r_time == 5, no future visits exist
      mcar_visitn <- (r_time + 1):5
      data_discont$y_mar[data_discont$id == i & data_discont$visitn %in% mcar_visitn] <- NA
    }
  }
  
  
  data_discont <- data_discont[data_discont$visitn > 0, ]
  data_wide <- data_discont %>%
    dplyr::select(id, groupx, visitn, y_mar, baseline_var, Discontinue,r) %>%
    pivot_wider(names_from = visitn, values_from = y_mar, names_prefix = "Y") %>%
    rename(ID = id, Arm = groupx, Y0 = baseline_var)
  
  
  # missingness indicator
  data_wide <- data_wide %>%
    mutate(
      MCAR_Indicator = if_else(is.na(Y1) | is.na(Y2) | is.na(Y3) | is.na(Y4) | is.na(Y5), 1, 0)
    )
  
  data <- data_wide
  
  mod <- rstan::stan_model(file = "/home/lsh1901704/stan_beta_5visits_miss2.stan")

   N <- nrow(data)
  theta_jack <- numeric(N)  # Store treatment effects per leave-one-out sample
  
  for (p in 1:N) {
    data_loo <- data[-p, ]  # Leave one subject out
    
    # Recalculate all data-dependent inputs using data_loo
    treatVec <- as.numeric(data_loo$Arm == "Treatment")
    discontVec <- as.numeric(data_loo$Discontinue)
    
    discontInd <- which(data_loo$Discontinue == 1 & data_loo$Arm == "Treatment")
    discontonlyInd <- which(data_loo$Discontinue == 1 & data_loo$MCAR_Indicator != 1 & data_loo$Arm == "Treatment")
    missInd <- which(data_loo$MCAR_Indicator == 1 & data_loo$Arm == "Treatment")
    missInd_c <- which(data_loo$MCAR_Indicator == 1 & data_loo$Arm == "Control")
    nondiscontmissInd <- which(data_loo$MCAR_Indicator != 1 & data_loo$Discontinue != 1)
    
    yVec <- as.matrix(data_loo[, c("Y1", "Y2", "Y3", "Y4", "Y5")])
    yVec[is.na(yVec)] <- 1000
    
    t_indiv <- ifelse((data_loo$Discontinue == 1 & data_loo$Arm == "Treatment") |
                        (data_loo$MCAR_Indicator == 1 & data_loo$Arm == "Control"), 
                      data_loo$r, 5)
    
    numMisPerT <- table(factor(data_loo$r[data_loo$Arm == "Treatment"], levels = 0:4))
    numMisPerC <- table(factor(data_loo$r[data_loo$Arm == "Control"], levels = 0:4))
    numMisPerT <- numMisPerT[c("1", "2", "3", "4", "0")]
    numMisPerC <- numMisPerC[c("1", "2", "3", "4", "0")]
    
    datForStan <- list(
      K = 5,
      J = 1,
      N = nrow(data_loo),
      t = 1,
      x = as.matrix(data_loo$Y0),
      y = yVec,
      W = treatVec,
      discont_index = as.array(discontInd),
      miss_index = as.array(missInd),
      miss_index_c = as.array(missInd_c),
      discontonly_index = as.array(discontonlyInd),
      nondiscontmiss_index = nondiscontmissInd,
      num_discont = length(discontInd),
      num_miss = length(missInd),
      num_miss_c = length(missInd_c),
      num_discontonly = length(discontonlyInd),
      num_nondiscontmiss = length(nondiscontmissInd),
      Z = discontVec,
      priorMiss = rep(1, 5),
      numMissPer0 = as.vector(numMisPerC),
      numMissPer1 = as.vector(numMisPerT),
      t_indiv = t_indiv
    )
    
  
  fit <- rstan::optimizing(mod,
    data = datForStan,
    as_vector = FALSE,
    init = "0"
  )
    
    params <- fit$par
    y_pred <- yVec
    Sigma <- params$Sigma
    beta <- params$beta
    alpha0 <- params$alpha0
    alpha1 <- params$alpha1
    k_0 <- params$k_0
    x <- datForStan$x
    W <- datForStan$W
    t_indiv <- datForStan$t_indiv
    
    for (n in 1:nrow(data_loo)) {
      t_n <- t_indiv[n]
      post_n <- 5 - t_n
      x_n <- x[n, , drop = FALSE]
      y_obs <- yVec[n, 1:t_n, drop = FALSE]
      
      y_pred[n, 1:t_n] <- as.vector(y_obs)
      
      if (n %in% missInd && post_n > 0) {
        Sigma_11 <- Sigma[1:t_n, 1:t_n, drop = FALSE]
        Sigma_21 <- Sigma[(t_n + 1):5, 1:t_n, drop = FALSE]
        Sigma_12 <- t(Sigma_21)
        Sigma_22 <- Sigma[(t_n + 1):5, (t_n + 1):5, drop = FALSE]
        Sigma_11_inv <- solve(Sigma_11)
        
        mu_obs <- if (W[n] == 1) alpha1[1:t_n] else alpha0[1:t_n]
        mu_obs <- mu_obs + as.vector(beta[1:t_n, , drop = FALSE] * as.numeric(x_n))
        mu_miss <- if (W[n] == 1) alpha0[(t_n + 1):5] else alpha0[(t_n + 1):5]
        mu_miss <- mu_miss + as.vector(beta[(t_n + 1):5, , drop = FALSE] * as.numeric(x_n))
        
        mu_post <- mu_miss + as.vector(Sigma_21 %*% Sigma_11_inv %*% (as.vector(y_obs) - mu_obs)) +
          k_0 * (alpha1[t_n] - alpha0[t_n])
        
        y_pred[n, (t_n + 1):5] <- mu_post
      }
      
      if (n %in% missInd_c && post_n > 0) {
        Sigma_11 <- Sigma[1:t_n, 1:t_n, drop = FALSE]
        Sigma_21 <- Sigma[(t_n + 1):5, 1:t_n, drop = FALSE]
        Sigma_12 <- t(Sigma_21)
        Sigma_22 <- Sigma[(t_n + 1):5, (t_n + 1):5, drop = FALSE]
        Sigma_11_inv <- solve(Sigma_11)
        
        mu_obs <- if (W[n] == 1) alpha1[1:t_n] else alpha0[1:t_n]
        mu_obs <- mu_obs + as.vector(beta[1:t_n, , drop = FALSE] * as.numeric(x_n))
        mu_miss <- if (W[n] == 1) alpha1[(t_n + 1):5] else alpha0[(t_n + 1):5]
        mu_miss <- mu_miss + as.vector(beta[(t_n + 1):5, , drop = FALSE] * as.numeric(x_n))
        
        mu_post <- mu_miss + as.vector(Sigma_21 %*% Sigma_11_inv %*% (as.vector(y_obs) - mu_obs))
        
        y_pred[n, (t_n + 1):5] <- mu_post
      }
    }
    
    y_pred_df <- as.data.frame(y_pred)
    colnames(y_pred_df) <- paste0("y_pred_", 1:5)
    y_pred_df$ID <- seq_len(nrow(data_loo))
    data_loo$ID <- seq_len(nrow(data_loo))
    data_comb <- full_join(data_loo, y_pred_df, by = "ID")
    
    lm_fit <- lm(y_pred_5 ~ as.factor(Arm) + Y0, data = data_comb)
    theta_jack[p] <- coef(lm_fit)[2]  # Save treatment effect
  }
  

  theta_hat <- mean(theta_jack)
  mean_impute[j] <- theta_hat
  mean_impute_se[j] <- sqrt((N - 1) / N * sum((theta_jack - theta_hat)^2))
  
  print(j)
}

res <- data.frame(mean_impute,mean_impute_se)

# Save results uniquely for each job instance
output_csv <- paste0("results_low_beta_jk90_miss_batch_", batch, ".csv")


write.csv(res, output_csv)
#save(res, file = output_rdata)

#print(paste("Results saved to:", output_csv, "and", output_rdata))