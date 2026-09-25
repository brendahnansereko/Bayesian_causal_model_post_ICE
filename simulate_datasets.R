# --- sim_generator.R ---

# Capture the Slurm array task ID
args <- commandArgs(trailingOnly = TRUE)
batch <- ifelse(length(args) > 0, as.numeric(args[1]), 1)

# Set up random seed for reproducibility (L'Ecuyer-CMRG for parallel streams)
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(69012365)  # Base seed
s <- .Random.seed
for (i in 1:batch) {
  s <- nextRNGStream(s)
}
.GlobalEnv$.Random.seed <- s


# Output directory
out_dir <- "/home/lsh1901704/pioneer_paper2/datasets/k_1/low"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# Load necessary library
library(MASS) # For multivariate normal simulation
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

# Parameters
#set.seed(123) # For reproducibility
MAR <- MAR_se <- causal <- causal_se <- k_0 <- k_0_se <- 0

iter <- 50

k <- 1

batch_data_list <- list()

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
                                 beta_main = rep(-15, times = 10),
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
  # 1. Conditional imputation for post-ICE data
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

   
  print(j)

# Store the processed dataset for this iteration
  batch_data_list[[j]] <- data_discont
  if(j %% 10 == 0) print(paste("Iteration", j, "complete"))
}

# Save all 50 datasets as one batch file
save_path <- file.path(out_dir, paste0("batch_", batch, "_datasets.rds"))
saveRDS(batch_data_list, file = save_path)
print(paste("Batch", batch, "saved to", save_path))
