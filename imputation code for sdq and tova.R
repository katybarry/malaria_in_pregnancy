# =============================================================================
# Script:  02_multiple_imputation.R
# Author:  Katharine Barry
# Purpose: Multiple imputation of missing covariates for the SDQ and TOVA
#          analysis populations created by 01_data_preparation_MiP_SDQ_TOVA.R
#

# =============================================================================
library(dplyr)
library(mice)

sdq_pop  <- readRDS("malaria_sdq_pop.rds")
tova_pop <- readRDS("malaria_tova_pop.rds")

# ---- Variables ----
analysis_covs <- c("age_at_delivery", "v01bas_gravidity_cat", "pre_pregnancy_BMI",
                   "v01med_iptpgroup_1", "ID_visit1", "helminth_inf1",
                   "mother_edu_1year", "TOV_stot",
                   "child_sex")   # precision variable (not a confounder)


auxiliary     <- c("preterm_birth", "MSE_v01bir_bweight_5", "CRP_first_visit")

exposure_pred <- c("MiP")     # complete: used as predictors

carry_along   <- c("mide",
                   "ANV1_malaria", "ANV2_malaria", "ANV3_malaria", "PCR_malaria",
                   "birth_malaria", "PCR_placenta", "malaria_emerg", "malaria_pre_birth",
                   "malaria_episodes", "malaria_episodes_3cat")

stopifnot(!anyNA(sdq_pop$child_sex), !anyNA(tova_pop$child_sex))

prep <- function(df, outcomes) {
  df %>%
    mutate(malaria_episodes_3cat = factor(case_when(MiP == 0 ~ "0",
                                                    malaria_episodes == 1 ~ "1",
                                                    malaria_episodes >= 2 ~ "2+"),
                                          levels = c("0", "1", "2+"))) %>%
    select(all_of(c(carry_along, exposure_pred, outcomes, analysis_covs, auxiliary))) %>%
    mutate(
      mide = factor(mide),
      across(c(v01bas_gravidity_cat, v01med_iptpgroup_1, ID_visit1, helminth_inf1,
               mother_edu_1year, child_sex, preterm_birth, CRP_first_visit, MiP,
               ANV1_malaria, ANV2_malaria, ANV3_malaria, PCR_malaria, birth_malaria,
               PCR_placenta, malaria_pre_birth), as.factor),
      across(c(age_at_delivery, pre_pregnancy_BMI, TOV_stot, MSE_v01bir_bweight_5), as.numeric)
    )
}

# no_pred: variables kept in the imputed data for analysis only (not imputed,
# not used as predictors), e.g. the TOVA subscales, which have their own
# missing values and are near-collinear with the ADHD index.
run_mi <- function(dat, outcomes, no_pred = character(0), m = 20, maxit = 20, seed = 123) {
  cat("\nMissing values before imputation:\n")
  print(colSums(is.na(dat))[colSums(is.na(dat)) > 0])
  
  # Default methods: pmm (continuous), logreg (binary), polyreg (>2 categories).
  # Complete variables automatically get "" (not imputed).
  meth <- make.method(dat)
  meth[c(carry_along, exposure_pred, outcomes)] <- ""   # never impute exposure/outcome/carried
  
  # Predictor matrix: everything predicts everything, except carried variables
  pred <- make.predictorMatrix(dat)
  pred[, carry_along] <- 0
  pred[carry_along, ] <- 0
  pred[, no_pred] <- 0
  pred[no_pred, ] <- 0
  if ("sdq6ans_totalsdq" %in% colnames(pred)) pred[, "sdq6ans_totalsdq"] <- 0
  cat("\nImputation methods:\n"); print(meth[meth != ""])
  
  imp <- mice(dat, m = m, maxit = maxit, method = meth, predictorMatrix = pred,
              seed = seed, printFlag = FALSE)
  
  cat("\nLogged events (should be NULL or only harmless):\n"); print(imp$loggedEvents)
  cat("\nMissing in analysis covariates after imputation (should all be 0):\n")
  print(colSums(is.na(complete(imp, 1)[, analysis_covs])))
  imp
}

# ---- SDQ population ----
sdq_outcomes <- c("sdq6ans_internalizing", "sdq6ans_externalizing", "sdq6ans_totalsdq")
dat_sdq <- prep(sdq_pop, sdq_outcomes)
imp_sdq <- run_mi(dat_sdq, sdq_outcomes)

# ---- TOVA population ----
tova_outcomes  <- c("EXP_scoreADHD")   # complete; used as a predictor
tova_subscales <- c("EXP_dprimestdscore", "EXP_responsetimemsec",
                    "EXP_responsetimevariabilitymsec", "EXP_commission", "EXP_omission")
dat_tova <- prep(tova_pop, c(tova_outcomes, tova_subscales))
imp_tova <- run_mi(dat_tova, c(tova_outcomes, tova_subscales), no_pred = tova_subscales)


# ---- Save the mids objects (not the stacked long data) ----
saveRDS(imp_sdq,  "imp_sdq.rds")
saveRDS(imp_tova, "imp_tova.rds")





