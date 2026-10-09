# =============================================================================
# Script:  04_sensitivity_complete_malaria.R
# Author:  Katharine Barry
# Purpose: Sensitivity analysis restricted to mothers with complete malaria
#          testing at ANV1, ANV2 and delivery (smear + peripheral PCR +
#          placental PCR). 
# =============================================================================
library(dplyr)
library(mice)


# =============================================================================
# 1. SENSITIVITY POPULATIONS
# =============================================================================
# Complete malaria testing = a result for all five tests below.
# Unscheduled consultations are not required: they only happened when a woman
# was unwell, so no consultation is not a missing test.
malaria_tests <- c(ANV1_malaria = "ANV1 smear",
                   ANV2_malaria = "ANV2 smear",
                   ANV3_malaria = "Delivery smear",
                   PCR_malaria  = "Peripheral PCR at delivery",
                   PCR_placenta = "Placental PCR")

complete_malaria <- function(df) rowSums(is.na(df[, names(malaria_tests)])) == 0

cohort   <- readRDS("cohort_alive_eligible_6y.rds")
sdq_pop  <- readRDS("malaria_sdq_pop.rds")
tova_pop <- readRDS("malaria_tova_pop.rds")

# ---- Flowchart: sequential exclusions from the eligible population ----
cat("\n========== SENSITIVITY POPULATION (complete malaria testing) ==========\n")
cat("Alive and eligible at 6 years:", nrow(cohort), "\n")
remaining <- cohort
for (v in names(malaria_tests)) {
  n_missing <- sum(is.na(remaining[[v]]))
  remaining <- remaining %>% filter(!is.na(.data[[v]]))
  cat(sprintf("  Missing %-28s -%3d  -> %d\n", malaria_tests[[v]], n_missing, nrow(remaining)))
}
cat("Complete malaria testing:      ", nrow(remaining), "\n")
cat("  of whom assessed at 6 years: ", sum(remaining$i_EXPL %in% 1), "\n")

cat("\nMissing each test (not sequential; a woman can miss several):\n")
print(colSums(is.na(cohort[, names(malaria_tests)])))

# ---- Analysis populations: main populations restricted to complete testing ----
sdq_sens  <- sdq_pop  %>% filter(complete_malaria(.))
tova_sens <- tova_pop %>% filter(complete_malaria(.))
cat("\nSDQ sensitivity population: ", nrow(sdq_sens),  "of", nrow(sdq_pop),  "\n")
cat("TOVA sensitivity population:", nrow(tova_sens), "of", nrow(tova_pop), "\n")

# ---- Who is excluded? (SDQ population) ----
sdq_pop %>%
  mutate(group = ifelse(complete_malaria(.), "Complete testing", "Incomplete testing")) %>%
  group_by(group) %>%
  summarise(n           = n(),
            MiP_pct     = round(100 * mean(MiP == 1), 1),
            preterm_pct = round(100 * mean(preterm_birth == 1, na.rm = TRUE), 1),
            int_mean    = round(mean(sdq6ans_internalizing), 2),
            ext_mean    = round(mean(sdq6ans_externalizing), 2)) %>%
  print()

saveRDS(sdq_sens,  "malaria_sdq_sens_pop.rds")
saveRDS(tova_sens, "malaria_tova_sens_pop.rds")


# =============================================================================
# 2. MULTIPLE IMPUTATION WITHIN THE SENSITIVITY POPULATIONS
# =============================================================================
# Same imputation model as 02_multiple_imputation.R.
analysis_covs <- c("age_at_delivery", "v01bas_gravidity_cat", "pre_pregnancy_BMI",
                   "v01med_iptpgroup_1", "ID_visit1", "helminth_inf1",
                   "mother_edu_1year", "TOV_stot",
                   "child_sex")   # precision variable
auxiliary     <- c("preterm_birth", "MSE_v01bir_bweight_5", "CRP_first_visit")
exposure_pred <- c("MiP")
carry_along   <- c("mide",
                   "ANV1_malaria", "ANV2_malaria", "ANV3_malaria", "PCR_malaria",
                   "birth_malaria", "PCR_placenta", "malaria_emerg", "malaria_pre_birth",
                   "malaria_episodes", "malaria_episodes_3cat")

sdq_outcomes   <- c("sdq6ans_internalizing", "sdq6ans_externalizing", "sdq6ans_totalsdq")
tova_outcomes  <- c("EXP_scoreADHD")
tova_subscales <- c("EXP_dprimestdscore", "EXP_responsetimemsec",
                    "EXP_responsetimevariabilitymsec", "EXP_commission", "EXP_omission")

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

run_mi <- function(dat, outcomes, no_pred = character(0), m = 20, maxit = 20, seed = 123) {
  cat("\nMissing values before imputation:\n")
  print(colSums(is.na(dat))[colSums(is.na(dat)) > 0])
  
  meth <- make.method(dat)
  meth[c(carry_along, exposure_pred, outcomes)] <- ""   # never impute exposure/outcome/carried
  
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

imp_sdq_sens  <- run_mi(prep(sdq_sens, sdq_outcomes), sdq_outcomes)
imp_tova_sens <- run_mi(prep(tova_sens, c(tova_outcomes, tova_subscales)),
                        c(tova_outcomes, tova_subscales), no_pred = tova_subscales)

# Diagnostics
plot(imp_sdq_sens)
plot(imp_tova_sens)

saveRDS(imp_sdq_sens,  "imp_sdq_sens.rds")
saveRDS(imp_tova_sens, "imp_tova_sens.rds")


# =============================================================================
# 3. REGRESSION MODELS (same specification as 03_regression_models.R)
# =============================================================================
covs <- analysis_covs

# Exposure = extra adjustment for that exposure only
exposures <- list(
  MiP                   = character(0),
  ANV1_malaria          = character(0),
  ANV2_malaria          = character(0),
  birth_malaria         = "malaria_pre_birth",
  ANV3_malaria          = "malaria_pre_birth",
  PCR_malaria           = "malaria_pre_birth",
  PCR_placenta          = "malaria_pre_birth",
  malaria_episodes_3cat = character(0)
)

# One outcome x one exposure: fit in each imputed dataset, pool with Rubin's
# rules, keep only the exposure row(s)
run_model <- function(imp, outcome, exposure, extra = character(0), adjusted = TRUE) {
  rhs  <- c(exposure, if (adjusted) c(covs, extra))
  fstr <- paste(outcome, "~", paste(rhs, collapse = " + "))
  fit  <- with(imp, lm(as.formula(fstr)))
  summary(pool(fit), conf.int = TRUE) %>%
    filter(startsWith(as.character(term), exposure)) %>%
    transmute(outcome = !!outcome, exposure = !!exposure, term = as.character(term),
              model = ifelse(adjusted, "Adjusted", "Unadjusted"),
              N = nobs(fit$analyses[[1]]),
              beta = estimate, lower = `2.5 %`, upper = `97.5 %`, p = p.value)
}

# Every outcome x every exposure, unadjusted and adjusted
run_all <- function(imp, outcomes) {
  bind_rows(lapply(outcomes, function(o)
    bind_rows(lapply(names(exposures), function(e)
      bind_rows(run_model(imp, o, e, adjusted = FALSE),
                run_model(imp, o, e, exposures[[e]], adjusted = TRUE))))))
}

sens_results <- bind_rows(
  run_all(imp_sdq_sens,  sdq_outcomes),
  run_all(imp_tova_sens, c(tova_outcomes, tova_subscales))
) %>%
  mutate(beta_ci = sprintf("%.2f (%.2f, %.2f)", beta, lower, upper),
         p = round(p, 3))

print(as_tibble(sens_results), n = Inf)


# =============================================================================
# 4. SAVE
# =============================================================================
saveRDS(sens_results, "sensitivity_results.rds")
writexl::write_xlsx(sens_results, "sensitivity_results.xlsx")
