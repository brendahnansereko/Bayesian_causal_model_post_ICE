
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

#Print batch number (for debugging)
print(paste("Running batch:", batch))


# Load necessary library
library(MASS) # For multivariate normal simulation
library(tidyr)
library(dplyr)
library(rbmi)
library(rstan)
library(cmdstanr)

rstan_options(auto_write = TRUE)
rstan_options(cxx = "/opt/ohpc/pub/compiler/gcc/9.4.0/bin/g++")


if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2", repos = "https://cran.r-project.org")
}

library(reshape2)

setwd("/home/lsh1901704/")
# Required libraries
library(rbmi)
library(rstan)
library(dplyr)
library(purrr)
library(mmrm)
library(tidyr)
library(MASS)  # For mvrnorm()

#set.seed(423456)

# Result storage
MAR <- MAR_se <- causal <- causal_se <- k_0 <- k_0_se <- mean_impute_se  <- jackknife_se <- boot_se <- mean_impute <- jackknife_se <- 0

iter <- 50

k <- 1

for (j in 1:iter) {
  n <- 500
  
  # Visits
  visits  <-  c(0, 4, 8, 14, 20, 26)
  n_visit <-  length(visits)
  
  # Mean trajectory control
  muC <- c(7.92, 7.82, 7.8, 7.8, 7.78, 7.78)
  
  # Mean trajectory intervention
  muT <-  c(7.92, 7.55, 7.20, 7.1, 7.05, 7.05)
  
  # create Sigma
  # Spatial power correlation matrix - corresponds to AR(1) for equidistant visits
  rho      <- 0.8
  exponent <- abs(matrix(visits, nrow = n_visit, ncol = n_visit, byrow = TRUE) - visits) / (visits[2] - visits[1])
  corr_mat <- rho^exponent
  
  # Variance at baseline is pooled variance from each group
  sigma2_bl = 0.48 
  
  # Variance for treatment and control group at each visit
  sigma2_trt = c(sigma2_bl, 0.8, 1.1, 1.4, 1.23, 1.48)
  sigma2_ctl = c(sigma2_bl, 0.8, 1.1, 1.4, 1.23, 1.48)
  
  #sigma2_ctl <- sigma2_trt <-  rep(3, 6)
  
  # Percentage of pts. that discontinued randomized trt
  rate_trtdiscont_trt <- (175 - c(175, 169, 164, 156, 151, 149)) / 175 
  rate_trtdiscont_ctl <- (178 - c(178, 172, 165, 158, 140, 133)) / 178 
  
  # covariance matrix for treatment group
  mat1 <- matrix(sigma2_ctl, nrow = n_visit, ncol = n_visit, byrow = TRUE)
  mat2 <-  matrix(sigma2_ctl, nrow = n_visit, ncol = n_visit, byrow = FALSE)
  Sigma <- sqrt(mat1 * mat2) * corr_mat
  
  # Set probability of discontinuation
  probDisc_C <- 0.02
  probDisc_T <- 0.03
  or_outcome <- 1.10 # +1 point increase => +10% odds of discontinuation
  
  # Set drop-out rate following discontinuation
  prob_dropout <- 0.5
  
  # Set simulation parameters of the control group
  parsC <- set_simul_pars(
    mu = muC,
    sigma = Sigma,
    n = n)
  
  # Set simulation parameters of the intervention group
  parsT <- set_simul_pars(
    mu = muT,
    sigma = Sigma,
    n = n)
  
  # Simulate data
  data <- simulate_data(
    pars_c = parsC,
    pars_t = parsT,
    post_ice1_traj = "MAR")
  
  data <- subset(data, select = c("id","visit","group","outcome_bl","outcome_noICE"))
  colnames(data) <- c("id","visit","group","y_bl","y")
  data$trt <- ifelse(data$group=="Control",0,1)
  data$groupx <- ifelse(data$group=="Control","Control","Treatment")
  
  # Beta parameter for logistic regression modeling of DAR trt discont
  beta_discont_dar <- data.frame(groupx = rep(c("Control", "Treatment"), each = 5),
                                 visitn = rep(1:5, times = 2),
                                 beta_main = rep(-13, times = 10),
                                 beta_prev_resp = c(1.42, 1.14, 1.33, 1.51, 1.46, 
                                                    1.42, 1.14, 1.47, 1.48, 1.40), 
                                 beta_current_resp = rep(0, times = 10),
                                 beta_baseline = c(0, 0.3, 0.1, 0.05, 0, 
                                                   0, 0.3, 0.1, 0.05, 0), 
                                 stringsAsFactors = FALSE)
  
  #_____________________________________________
  # simulate Discontinuation
  #______________________________________________
  is_dnar <- any(names(beta_discont_dar) == "beta_current_resp")
  
  data$visitn <- as.numeric(data$visit) -1
  data$baseline_var <- data$y_bl
  data$response_ontrt <- data$y
  
  data_discont <- full_join(data, beta_discont_dar, 
                            by = c("groupx", "visitn")) 
  
  data_discont <- data_discont %>% 
    group_by(id) %>% 
    mutate(prev_response = lag(response_ontrt, default = 0),
           logit_prob = beta_main + beta_prev_resp * prev_response + beta_baseline * baseline_var,
           logit_prob = if(is_dnar){logit_prob + beta_current_resp * response_ontrt} else {logit_prob},
           prob = 1/(1+exp(-logit_prob)),
           prob = if_else(is.na(prob) | visitn==1 , 0, prob),
           ontrt = 1-rbinom(n = n(), size = 1, prob = prob)) %>%  
    mutate(first_offtrt_visit = min(visitn[ontrt == 0], max(visitn)),
           ontrt = if_else(visitn > first_offtrt_visit, 0, ontrt),
           r = ifelse(max(first_offtrt_visit)==5 & min(ontrt)==1,0,first_offtrt_visit-1)) %>% 
    ungroup()
  
  data_discont <- data_discont %>% 
    group_by(id) %>% 
    mutate(Discontinue = if_else(any(ontrt == 0), 1, 0)) %>%
    ungroup()
  
  data_discont <- data_discont %>% 
    mutate(
      y_mar = if_else(ontrt == 1 | (ontrt == 0 & groupx == "Control"), response_ontrt, as.numeric(NA))
    ) 
  
  # _______________________________________
  # 1. Conditional imputation
  # _______________________________________
  
  for (i in unique(data_discont$id[data_discont$groupx == "Treatment" & data_discont$Discontinue == 1])) {
    person_data <- data_discont[data_discont$id == i, ]
    r_time <- max(person_data$r)
    
    if (r_time < 6 & r_time != 0) {
      obs_idx <- 2:(r_time+1)
      fut_idx <- (r_time+2):6
      
      Sigma_11 <- as.matrix(Sigma[obs_idx , obs_idx , drop = FALSE])
      Sigma_22 <- as.matrix(Sigma[fut_idx , fut_idx , drop = FALSE])
      Sigma_21 <- as.matrix(Sigma[fut_idx , obs_idx , drop = FALSE]) 
      Sigma_12 <- t(Sigma_21)
      
      
      beta <- Sigma_21 %*% solve(Sigma_11)
      Sigma_cond <- Sigma_22 - beta %*% (Sigma_12)
      
      
      Y_obs <- as.numeric(person_data$y_mar[person_data$visitn %in% (obs_idx-1)])
      mu_obs <- muT[obs_idx]
      mu_obs_control <- muC[obs_idx]
      mu_future <- muT[fut_idx ]
      mu_control_future <- muC[fut_idx ]
      
      mu_cond <- beta %*% (Y_obs - mu_obs)
      mu_cond <- mu_cond + k * (muT[r_time+1] - muC[r_time+1]) + mu_control_future
      
      imputed_vals <- MASS::mvrnorm(1, mu_cond, Sigma_cond)
      # print(paste("Imputed values for patient", i, ":", imputed_vals))
      
      for (p in seq_along(fut_idx)) {
        data_discont$y_mar[data_discont$id == i & data_discont$visitn == (fut_idx[p]-1)] <- imputed_vals[p]
      }
    }
  }
  
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
  
  mod <- rstan::stan_model(file = "/home/lsh1901704/stan_beta_5visits_miss2_impute.stan")
  
  N <- nrow(data)

   B <- 200  # number of bootstrap samples
  theta_boot <- numeric(B)
  
  for (b in 1:B) {
  #  set.seed(1000 + b)  # optional, for reproducibility
    
    boot_indices <- sample(1:N, size = N, replace = TRUE)
    data_boot <- data[boot_indices, ]
    
    # Recalculate all data-dependent inputs using data_boot
    treatVec <- as.numeric(data_boot$Arm == "Treatment")
    discontVec <- as.numeric(data_boot$Discontinue)
    
    discontInd <- which(data_boot$Discontinue == 1 & data_boot$Arm == "Treatment")
    discontonlyInd <- which(data_boot$Discontinue == 1 & data_boot$MCAR_Indicator != 1 & data_boot$Arm == "Treatment")
    missInd <- which(data_boot$MCAR_Indicator == 1 & data_boot$Arm == "Treatment")
    missInd_c <- which(data_boot$MCAR_Indicator == 1 & data_boot$Arm == "Control")
    nondiscontmissInd <- which(data_boot$MCAR_Indicator != 1 & data_boot$Discontinue != 1)
    
    yVec <- as.matrix(data_boot[, c("Y1", "Y2", "Y3", "Y4", "Y5")])
    yVec[is.na(yVec)] <- 1000
    
    t_indiv <- ifelse((data_boot$Discontinue == 1 & data_boot$Arm == "Treatment") |
                        (data_boot$MCAR_Indicator == 1 & data_boot$Arm == "Control"), 
                      data_boot$r, 5)
    
    numMisPerT <- table(factor(data_boot$r[data_boot$Arm == "Treatment"], levels = 0:4))
    numMisPerC <- table(factor(data_boot$r[data_boot$Arm == "Control"], levels = 0:4))
    numMisPerT <- numMisPerT[c("1", "2", "3", "4", "0")]
    numMisPerC <- numMisPerC[c("1", "2", "3", "4", "0")]
    
    datForStan <- list(
      K = 5,
      J = 1,
      N = nrow(data_boot),
      t = 1,
      x = as.matrix(data_boot$Y0),
      y = yVec,
      W = treatVec,
      discont_index = discontInd,
      miss_index = missInd,
      miss_index_c = missInd_c,
      discontonly_index = discontonlyInd,
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
    
    for (n in 1:nrow(data_boot)) {
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
    y_pred_df$ID <- seq_len(nrow(data_boot))
    data_boot$ID <- seq_len(nrow(data_boot))
    data_comb <- full_join(data_boot, y_pred_df, by = "ID")
    
    lm_fit <- lm(y_pred_5 ~ as.factor(Arm) + Y0, data = data_comb)
    theta_boot[b] <- coef(lm_fit)[2]
    
  }
  
  theta_hat <- mean(theta_boot)
  mean_impute_se[j] <- sd(theta_boot)
  mean_impute[j] <- theta_hat
  
  print(j)
}

res <- data.frame(mean_impute,mean_impute_se)

# Save results uniquely for each job instance
output_csv <- paste0("results_high_beta_miss_batch_", batch, ".csv")


write.csv(res, output_csv)
#save(res, file = output_rdata)

#print(paste("Results saved to:", output_csv, "and", output_rdata))