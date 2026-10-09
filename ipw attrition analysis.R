# =============================================================================
# Script:  07_ipw_attrition.R
# Author:  Katharine Barry
# Purpose: Inverse probability weighting (IPW) for attrition .

# =============================================================================
library(dplyr)
library(mice)
library(mitools)    # MIcombine(): Rubin's rules with any variance matrix
library(sandwich)   # robust (sandwich) variance for weighted models

cohort   <- readRDS("cohort_alive_eligible_6y.rds")
sdq_pop  <- readRDS("malaria_sdq_pop.rds")
tova_pop <- readRDS("malaria_tova_pop.rds")
imp_sdq  <- readRDS("imp_sdq.rds")
imp_tova <- readRDS("imp_tova.rds")


# =============================================================================
# 1-2. WEIGHTS
# =============================================================================
# Predictors of being included, measured at inclusion and birth (available for
# almost everyone, including children lost to follow-up). Education and wealth
# (1 year) are NOT used: missing for most children lost to follow-up.
# Missing predictor values: median (continuous) or "missing" category
# (categorical) + a missing indicator, so every eligible child gets a weight.
make_ipw <- function(cohort, included_ids) {
  fill_num <- function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x)
  fill_cat <- function(x) factor(ifelse(is.na(x), "missing", as.character(x)))
  
  d <- cohort %>%
    mutate(included   = as.integer(mide %in% included_ids),
           age_miss   = as.integer(is.na(age_at_delivery)),
           bmi_miss   = as.integer(is.na(pre_pregnancy_BMI)),
           bw_miss    = as.integer(is.na(MSE_v01bir_bweight_5)),
           age_f      = fill_num(age_at_delivery),
           bmi_f      = fill_num(pre_pregnancy_BMI),
           bw_f       = fill_num(MSE_v01bir_bweight_5),
           grav_f     = fill_cat(v01bas_gravidity_cat),
           iptp_f     = fill_cat(v01med_iptpgroup_1),
           iron_f     = fill_cat(ID_visit1),
           helm_f     = fill_cat(helminth_inf1),
           mip_f      = fill_cat(MiP),
           preterm_f  = fill_cat(preterm_birth),
           sex_f      = fill_cat(child_sex))
  
  # Missing indicators only when there is something missing (else glm gives NA)
  miss_terms <- c("age_miss", "bmi_miss", "bw_miss")
  miss_terms <- miss_terms[sapply(miss_terms, function(v) sum(d[[v]]) > 0)]
  
  f <- as.formula(paste("included ~ age_f + grav_f + bmi_f + iptp_f + iron_f + helm_f +",
                        "mip_f + preterm_f + bw_f + sex_f",
                        if (length(miss_terms)) paste("+", paste(miss_terms, collapse = " + "))))
  ps <- glm(f, family = binomial, data = d)
  
  d$p_incl  <- predict(ps, type = "response")      # P(included | characteristics)
  d$ipw_raw <- mean(d$included) / d$p_incl          # stabilised weight
  q <- quantile(d$ipw_raw[d$included == 1], c(0.01, 0.99))
  d$ipw     <- pmin(pmax(d$ipw_raw, q[1]), q[2])     # truncated at 1st/99th pct
  
  cat("Children in weight model:", nrow(d), "| included:", sum(d$included), "\n")
  cat("Weights among included (should have mean close to 1):\n")
  print(summary(d$ipw[d$included == 1]))
  
  list(model = ps, data = d)
}

cat("\n================ SDQ population ================\n")
w_sdq  <- make_ipw(cohort, sdq_pop$mide)
cat("\n================ TOVA population ===============\n")
w_tova <- make_ipw(cohort, tova_pop$mide)

# Which characteristics predict being included? (odds ratios)
round(exp(cbind(OR = coef(w_sdq$model), confint.default(w_sdq$model))), 2)


# =============================================================================
# 3. BALANCE CHECK
# =============================================================================
# After weighting, the included children should look like ALL eligible children.
balance <- function(w) {
  d   <- w$data
  inc <- d$included == 1
  vals <- list(
    "Mother's age (mean)"         = d$age_f,
    "Pre-pregnancy BMI (mean)"    = d$bmi_f,
    "Birthweight, g (mean)"       = d$bw_f,
    "Primigravida (%)"            = 100 * (d$grav_f == "Primigravidae"),
    "4+ previous pregnancies (%)" = 100 * (d$grav_f == "4 or more"),
    "MiP (%)"                     = 100 * (d$mip_f == "1"),
    "Preterm birth (%)"           = 100 * (d$preterm_f == "1")
  )
  data.frame(
    characteristic      = names(vals),
    all_eligible        = sapply(vals, mean),
    included_unweighted = sapply(vals, function(x) mean(x[inc])),
    included_weighted   = sapply(vals, function(x) weighted.mean(x[inc], d$ipw[inc])),
    row.names = NULL
  ) %>% mutate(across(where(is.numeric), ~ round(.x, 2)))
}
cat("\nBalance, SDQ population:\n");  print(balance(w_sdq))
cat("\nBalance, TOVA population:\n"); print(balance(w_tova))


# =============================================================================
# 4. WEIGHTED MODELS ON THE IMPUTED DATA
# =============================================================================
# Add the weights to every imputed dataset
add_weights <- function(imp, w) {
  long <- complete(imp, action = "long", include = TRUE) %>%
    mutate(mide = as.character(mide)) %>%
    left_join(w$data %>% filter(included == 1) %>% select(mide, ipw), by = "mide")
  stopifnot(!anyNA(long$ipw))
  as.mids(long)
}
imp_sdq_w  <- add_weights(imp_sdq,  w_sdq)
imp_tova_w <- add_weights(imp_tova, w_tova)

# Same adjustment set and exposures as 03_regression_models.R
covs <- c("age_at_delivery", "v01bas_gravidity_cat", "pre_pregnancy_BMI",
          "v01med_iptpgroup_1", "ID_visit1", "helminth_inf1",
          "mother_edu_1year", "TOV_stot", "child_sex")
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

# Weighted lm in each imputed dataset, robust (HC0) variance because the weights
# are estimated, then Rubin's rules with MIcombine()
run_weighted <- function(imp, outcome, exposure, extra = character(0)) {
  f <- as.formula(paste(outcome, "~", paste(c(exposure, covs, extra), collapse = " + ")))
  fits <- lapply(seq_len(imp$m), function(i) lm(f, data = complete(imp, i), weights = ipw))
  comb <- MIcombine(lapply(fits, coef), lapply(fits, vcovHC, type = "HC0"))
  est <- coef(comb); se <- sqrt(diag(vcov(comb))); df <- comb$df
  keep <- startsWith(names(est), exposure)
  tibble(outcome = outcome, exposure = exposure, term = names(est)[keep],
         model = "Adjusted + IPW", N = nobs(fits[[1]]),
         beta  = est[keep],
         lower = est[keep] - qt(0.975, df[keep]) * se[keep],
         upper = est[keep] + qt(0.975, df[keep]) * se[keep],
         p     = 2 * pt(-abs(est[keep] / se[keep]), df[keep]))
}

run_all_w <- function(imp, outcomes) {
  bind_rows(lapply(outcomes, function(o)
    bind_rows(lapply(names(exposures), function(e) run_weighted(imp, o, e, exposures[[e]])))))
}

ipw_results <- bind_rows(
  run_all_w(imp_sdq_w,  c("sdq6ans_internalizing", "sdq6ans_externalizing", "sdq6ans_totalsdq")),
  run_all_w(imp_tova_w, "EXP_scoreADHD")
) %>%
  mutate(beta_ci = sprintf("%.2f (%.2f, %.2f)", beta, lower, upper), p = round(p, 3))

# Side by side with the main adjusted results (for the supplementary table)
main <- readRDS("regression_results.rds") %>%
  filter(model == "Adjusted") %>%
  select(outcome, exposure, term, main = beta_ci, p_main = p)

ipw_comparison <- ipw_results %>%
  select(outcome, exposure, term, N, ipw = beta_ci, p_ipw = p) %>%
  left_join(main, by = c("outcome", "exposure", "term")) %>%
  relocate(main, p_main, .after = N)

print(ipw_comparison, n = Inf)

saveRDS(ipw_results, "ipw_results.rds")
writexl::write_xlsx(list(ipw_vs_main = ipw_comparison, ipw_results = ipw_results),
                    "ipw_results.xlsx")