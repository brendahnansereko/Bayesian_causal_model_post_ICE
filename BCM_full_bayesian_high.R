
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
    rename(id = id, group = groupx, y_bl = baseline_var)
  
  data <- data_discont
  #____________________________________________________
  #        MAR imputation analysis
  #____________________________________________________
  
  data$visit <- as.factor(data$visit)
  
  data <- data %>%
    mutate(visit = fct_drop(visit))  # Drop unused levels from visitn
  
  #expand the data by adding rows with NAs for the visits that patients were missing:
  expandedData <- expand_locf(
    data,
    id = levels(data$id),
    visit = levels(data$visit),
    vars = c("y_bl","group"),
    group = c("id"),
    order = c("id", "visit")
  )
  
  expandedData <- expandedData[expandedData$visit != "0", ]
  # Create another data frame that records, for each patient,
  # the first visit after the  (ICE) occurred and the imputation approach for missing data after
  # the ICE.
  # Then proceed with the conversion and filtering
  ice_info_jr <- expandedData %>%
    group_by(id) %>%
    mutate(
      # Convert visit and first_offtrt_visit to numeric if they are factors
      visit = as.numeric(as.character(visit)),  
      first_offtrt_visit = as.numeric(as.character(first_offtrt_visit))
    ) %>%
    filter(
      (visit > first_offtrt_visit & first_offtrt_visit != 5) |  # Ensure the comparison works without factor issue
        (visit == first_offtrt_visit & first_offtrt_visit != 5)  # First visit (ignoring last visit)
    ) %>%
    mutate(strategy = as.character(ifelse(group == "Intervention", "JR", "MAR"))) %>%
    summarise(
      visit = as.factor(visit[1]),  # First affected visit
      strategy = strategy[1]  # Strategy for the first affected visit
    )
  
  # Impute and analyse - The following code is for using Stan to fit the Bayesian model:
  draw_obj <- draws(
    data = expandedData,
    data_ice = ice_info_jr,
    vars = set_vars(
      outcome = "y_mar",
      visit = "visit",
      subjid = "id",
      group = "group",
      covariates = c("y_bl*visit", "group*visit")),
    method = method_bayes(control=control_bayes(warmup=200,thin=1),
                          n_samples = 5000,
                          same_cov = TRUE
                          #seed = 484381136
    ))
  
  # now thin the draws for the imputation methods
  thinnedDraws <- draw_obj
  # take every 50th draw
  thinnedDraws$samples <- thinnedDraws$samples[seq(50, length(thinnedDraws$samples), by = 50)]
  # now modify method accordingly as if we had thinned in the first place
  thinnedDraws$method$n_samples <- length(thinnedDraws$samples)
  thinnedDraws$method$control$thin <- 50
  class(thinnedDraws$samples) <- c("sample_list", "list")
  
  #____________________________________________________
  #        RBI imputation analysis - JR
  #____________________________________________________
  
  jrImps <- impute(thinnedDraws,
                   references = c("Intervention" = "Control", "Control" = "Control"),
                   update_strategy = ice_info_jr)
  ancovaVars <- set_vars(
    subjid = "id",
    outcome = "y_mar",
    visit = "visit",
    group = "group",
    covariates = c("y_bl"))
  jrMIFits <- analyse(
    imputations = jrImps,
    vars = ancovaVars
  )
  jrMIResults <- pool(jrMIFits)
  jrMIResults
  j2r_mean[j] <- jrMIResults$pars$trt_5$est
  j2r_se[j] <- jrMIResults$pars$trt_5$se
  
  
  
  #____________________________________________________
  #        RBI imputation analysis - CIR
  #____________________________________________________
  
  ice_info_cir <- ice_info_jr %>%
    mutate(strategy = ifelse(strategy == "JR", "CIR", strategy))
  
  #set.seed(837263)
  cirImps <- impute(thinnedDraws,
                    references = c("Intervention" = "Control", "Control" = "Control"),
                    update_strategy = ice_info_cir)
  
  cirMIFits <- analyse(
    imputations = cirImps,
    vars = ancovaVars
  )
  cirMIResults <- pool(cirMIFits)
  cirMIResults
  cir_mean[j] <- cirMIResults$pars$trt_5$est
  cir_se[j] <- cirMIResults$pars$trt_5$se
  
}

#Data sets
res <- as.data.frame( cbind(cir_mean,cir_se,j2r_mean,j2r_se))

write.csv(res, paste0("res_high_", batch, ".csv"))


