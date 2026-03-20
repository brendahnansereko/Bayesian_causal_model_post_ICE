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
data_list <- readRDS(paste0("/home/lsh1901704/pioneer_paper2/datasets/high/batch_", batch, "_datasets.rds"))

# Parameters
MAR <- MAR_se <- causal <- causal_se <- k_0 <- k_0_se <- numeric(length(data_list))

iter <- length(data_list)

for (j in 1:iter) {
  
  # Use the data from the list
  data_discont <- data_list[[j]]

  # _______________________________________
  # 2.introduce monotone MCAR
  # _______________________________________
  # Introduce monotone MCAR in 10% of discontinuers in the Treatment group
  discont_treatment_ids <- unique(data_discont$id[data_discont$groupx == "Treatment" & data_discont$Discontinue == 1])
  n_mcar <- ceiling(length(discont_treatment_ids) * 0.20)
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
  n_mcar <- ceiling(length(discont_ctrl_ids) * 0.20)
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
  
  
  treatVec <- as.numeric(data$Arm == "Treatment")
  discontVec <- as.numeric(data$Discontinue)
  
  
  #discontInd <- which(data$Discontinue== 1 & data$Arm == "Treatment")
  #nondiscontInd <- which((data$Discontinue != 1 & data$Arm == "Treatment") |
                          # data$Arm == "Control")
  
  #identifying indices of that are discontinuing
  discontInd = which(data$Discontinue==1 & data$Arm == "Treatment")
  discontonlyInd = which(data$Discontinue==1 & data$MCAR_Indicator!=1 & data$Arm == "Treatment")
  missInd = which(data$MCAR_Indicator==1)
  nondiscontmissInd = which((data$MCAR_Indicator!=1 & data$Discontinue!=1 & data$Arm == "Treatment") |
                              data$Arm == "Control")
  
  yVec <- as.matrix(data[, c("Y1", "Y2", "Y3", "Y4", "Y5")])
  yVec[is.na(yVec)] <- 1000
  
  levels_r <- as.character(0:4)
  
  numMisPerT <- table(factor(data$r[data$Arm == "Treatment"], levels = 0:4))
  numMisPerC <- table(factor(data$r[data$Arm == "Control"], levels = 0:4))
  
  numMisPerC <- numMisPerC[c("1", "2", "3", "4", "0")]
  numMisPerT <- numMisPerT[c("1", "2", "3", "4", "0")]
  
  K <- 5
  
  t_indiv <- ifelse((data$Discontinue == 1 & data$Arm == "Treatment") |
                    (data$MCAR_Indicator == 1 & data$Arm == "Control") , data$r , 5)
  
  
  #preparing the data that would be passed to STAN
  datForStan = list(
    K = K, #number of post randomization time points
    J = 1, #number of covariates
    N = dim(data)[1], #number of units
    t = 1, # time of discontinuation
    # allSize = length(yVec), #number of timepoints*number of units
    x = as.matrix(data$Y0), #covariate matrix
    y = yVec, #outcomes matrix
    W = treatVec,
    discont_index = as.array(discontInd), #discontinuation indicator
    miss_index = as.array(missInd),  #missing data indicator
    discontonly_index = as.array(discontonlyInd), #discontunation only indicator
    nondiscontmiss_index = as.array(nondiscontmissInd), #non miss non discontinuation indicator
    num_discont = length(discontInd), #number of discontinued patients
    num_miss = length(missInd), #number of discontinued patients with missing data
    num_discontonly = length(discontonlyInd), #number of discontinued patients without missing data
    num_nondiscontmiss = length(nondiscontmissInd), #number of patients who didnot discontinue
    Z = discontVec, #discont indicator,
    priorMiss = rep(1, 5),
    #priorMiss = c(1.5725, 2.9450, 2.0450, 0.6250, 42.8125),
    #priorMiss = c(3.3052, 3.5478, 1.8818, 0.6418, 10.6234),
    numMissPer0 = as.vector(numMisPerC), #number of missing values in the control group per period
    numMissPer1 = as.vector(numMisPerT), #number of missing values in the treatment group per period
    t_indiv = t_indiv
  )
  
  
 
fit <- stan(file = "/home/lsh1901704/stan_beta_5visits_miss2.stan", data = datForStan, 
            iter = 1000, chains = 1,
            #control = list(max_treedepth = 15, adapt_delta = 0.99),
            warmup = 300
           )
  
  matrix_of_draws <- data.frame(as.matrix(fit))
  lm_fit <- summary(lm(Y5 ~ as.factor(Arm) + Y0, data = data))
  
  MAR[j] <- lm_fit$coefficients[2, 1]
  MAR_se[j] <- lm_fit$coefficients[2, 2]
  causal[j] <- mean(matrix_of_draws$causal_J2R)
  causal_se[j] <- sd(matrix_of_draws$causal_J2R)
  k_0[j] <- mean(matrix_of_draws$k_0)
  k_0_se[j] <- sd(matrix_of_draws$k_0)
  
  print(j)
}

res <- as.data.frame(cbind(MAR,MAR_se,causal,causal_se,k_0,k_0_se))

# Save results uniquely for each job instance
output_csv <- paste0("results_high_beta20_miss_batch_", batch, ".csv")


write.csv(res, output_csv)
#save(res, file = output_rdata)

#print(paste("Results saved to:", output_csv, "and", output_rdata))