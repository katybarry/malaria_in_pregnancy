# =============================================================================
# Author:  Katharine Barry
# Purpose: Data cleaning, population creation,and variable creation
# =============================================================================


# =============================================================================
# 0. SETUP AND PATHS
# =============================================================================
library(dplyr)
library(stringr)
library(tidyr)

path_tabmev7         <- "tabmev7.rds"
path_live_births_all <- "M:/child mortality paper and data/live_births_all.rds"

# Helper: print a flowchart step (mothers and children)
n_step <- function(df, label) {
  cat(sprintf("%-55s mothers: %4d | children: %4d\n", label,
              n_distinct(df$midm), n_distinct(df$mide)))
  invisible(df)
}


# =============================================================================
# 1. LOAD AND CLEAN tabmev7
# =============================================================================

tabmev7 <- readRDS(path_tabmev7)

tabmev7 <- tabmev7 %>%
  mutate(midm = str_replace(studysubjectid, "^MiPB-", "B"),
         mide = str_replace(subjectide,     "^MiCB-", "B"))

# Manual correction found in DEVINE

# took out id for github but correct id here

# Birth outcome (child level)

tabmev7 <- tabmev7 %>%
  mutate(birth_outcome = case_when(
    v05del_delivoutc_7 == "Live birt"               ~ "Live Birth",
    v05del_delivoutc_7 == "Stillbirt"               ~ "Still Birth",
    v05del_delivoutc_7 %in% c("Spontaneo", "Other") ~ "Abortion",
    TRUE                                            ~ "No delivery recorded"
  ))


# =============================================================================
# 2. FLOWCHART: MiPPAD -> APEC -> LIVE SINGLETONS -> ALIVE AT 6 YEARS
# =============================================================================
cat("\n---------------- FLOWCHART ----------------\n")

tabmev7 %>% n_step("MiPPAD Benin")                                   # 1183 mothers

# MiPPAD outcomes before APEC (for the top box of the flowchart)
tabmev7 %>% distinct(midm, .keep_all = TRUE) %>% count(birth_outcome)
tabmev7 %>% distinct(midm, .keep_all = TRUE) %>%
  count(vxxstu_notcompletersn_5, birth_outcome)                      # migrations, withdrawals, ...

table(tabmev7$vxxstu_other_5)


# ---- 2.1 APEC ----
APEC <- tabmev7 %>% filter(!is.na(i_GEAP)) %>% n_step("APEC")       # 1005 mothers

cat("\nAPEC mothers by birth outcome:\n")
print(APEC %>% distinct(midm, .keep_all = TRUE) %>% count(birth_outcome))

APEC %>% distinct(midm, .keep_all = TRUE) %>%
  count(vxxstu_notcompletersn_5, birth_outcome)     


APEC %>%
  distinct(midm, .keep_all = TRUE) %>%
  filter(vxxstu_notcompletersn_5 == "Other", birth_outcome == "No delivery recorded") %>%
  select(midm, mide, vxxstu_notcompletersn_5, vxxstu_other_5, v05del_delivoutc_7, v05del_deloutc_7oth)

# ---- 2.2 Live births ----
# CHANGED: the old APEC dataset was never restricted to live births, so 44
# stillbirths/abortions and 60 women with no delivery recorded stayed in
# APEC_noprob (967). They are now excluded explicitly.
live <- APEC %>% filter(birth_outcome == "Live Birth") %>% n_step("APEC live births")

# ---- 2.3 Singletons ----
twins_mothers <- live %>% filter(v05del_multbirth_10 %in% "Yes") %>% distinct(midm)
cat("Mothers with multiple births:", nrow(twins_mothers), "\n")       # 28?

# CHANGED: `!= "Yes"` silently dropped rows where v05del_multbirth_10 is NA.
singletons <- live %>%
  filter(!midm %in% twins_mothers$midm) %>%
  distinct(midm, .keep_all = TRUE) %>%
  n_step("Singleton live births")

# ---- 2.4 Neurological / genetic disorders (Amanda's exclusions) ----
cat("\nexcl_Amanda:\n"); print(table(singletons$excl_Amanda, useNA = "ifany"))
no_disorder <- singletons %>%
  filter(is.na(excl_Amanda) | excl_Amanda != 1) %>%
  n_step("No neurological/genetic disorder")

# ---- 2.5 Deaths before the 6-year assessment ----
# Uses the cleaned death classification from the child mortality paper
# (date-based and visit-coded deaths, with the 7 known false deaths overridden).
live_births_all <- readRDS(path_live_births_all)

death_info <- live_births_all %>%
  mutate(dead_before_6y_visit = (died_by_6y | coalesce(coded_death_6y, FALSE)) &
           !(override_bad_deathdate | override_false_stil)) %>%
  select(mide, dead_before_6y_visit)

no_disorder <- no_disorder %>%
  left_join(death_info, by = "mide") %>%
  # A child who attended EXPLORE was alive at the assessment, whatever the death record says
  mutate(dead = coalesce(dead_before_6y_visit, FALSE) & !(i_EXPL %in% 1))

cat("Children not found in live_births_all:", sum(is.na(no_disorder$dead_before_6y_visit)), "\n")
cat("Deaths before the 6-year assessment:", sum(no_disorder$dead), "\n")

cohort <- no_disorder %>%
  filter(!dead) %>%
  n_step("Alive and eligible at 6 years")


nrow(no_disorder)                         # should be 796 + 67 = 863
all(no_disorder$mide %in% APEC$mide)      # TRUE: everyone is from APEC

# Compare with all MiPPAD (twins included), from the mortality paper
sum(live_births_all$died_by_6y | coalesce(live_births_all$coded_death_6y, FALSE), na.rm = TRUE)


table(cohort$DEC_genet_dc_n)

# =============================================================================
# 3. EXPOSURE: MALARIA IN PREGNANCY
# =============================================================================
# Coding in the data:
#   AME_GEResultat1-3 (ANV1, ANV2, delivery smear): 1 = negative, 2 = positive, 3 = unknown
#   AME_GE1-4 (unscheduled consultations):          1 = positive, 2 = negative, 3 = unknown
#   (AME_GE5 only contains "3", AME_GE6 is empty)
#   AME_PCR_Acc1/2, AME_PCR_Plac1/2: "pos" / "neg"

cohort <- cohort %>%
  mutate(
    # ---- Scheduled smears ----
    ANV1_malaria = case_when(AME_GEResultat1 == 2 ~ 1, AME_GEResultat1 == 1 ~ 0, TRUE ~ NA_real_),
    ANV2_malaria = case_when(AME_GEResultat2 == 2 ~ 1, AME_GEResultat2 == 1 ~ 0, TRUE ~ NA_real_),
    ANV3_malaria = case_when(AME_GEResultat3 == 2 ~ 1, AME_GEResultat3 == 1 ~ 0, TRUE ~ NA_real_),  # delivery smear
    
    # ---- Peripheral PCR at delivery ----
    PCR_malaria1 = case_when(AME_PCR_Acc1 == "pos" ~ 1, AME_PCR_Acc1 == "neg" ~ 0, TRUE ~ NA_real_),
    PCR_malaria2 = case_when(AME_PCR_Acc2 == "pos" ~ 1, AME_PCR_Acc2 == "neg" ~ 0, TRUE ~ NA_real_),
    PCR_malaria  = case_when(PCR_malaria1 %in% 1 | PCR_malaria2 %in% 1 ~ 1,
                             PCR_malaria1 %in% 0 | PCR_malaria2 %in% 0 ~ 0,
                             TRUE ~ NA_real_),
    
    # ---- Placental PCR ----
    PCR_plac1    = case_when(AME_PCR_Plac1 == "pos" ~ 1, AME_PCR_Plac1 == "neg" ~ 0, TRUE ~ NA_real_),
    PCR_plac2    = case_when(AME_PCR_Plac2 == "pos" ~ 1, AME_PCR_Plac2 == "neg" ~ 0, TRUE ~ NA_real_),
    PCR_placenta = case_when(PCR_plac1 %in% 1 | PCR_plac2 %in% 1 ~ 1,
                             PCR_plac1 %in% 0 | PCR_plac2 %in% 0 ~ 0,
                             TRUE ~ NA_real_),
    
    # ---- Unscheduled (emergency) consultations ----

    # and malaria_emerg = 1 if any positive, 0 otherwise.
    n_emerg_pos   = rowSums(across(all_of(paste0("AME_GE", 1:4)), ~ .x %in% 1)),
    n_emerg_neg   = rowSums(across(all_of(paste0("AME_GE", 1:4)), ~ .x %in% 2)),
    malaria_emerg = as.numeric(n_emerg_pos > 0),
    
    # ---- Malaria at delivery (smear or peripheral PCR; placenta kept separate) ----
    # CHANGED: the old code coded birth_malaria = 0 when BOTH tests were missing.
    birth_malaria = case_when(ANV3_malaria %in% 1 | PCR_malaria %in% 1 ~ 1,
                              ANV3_malaria %in% 0 | PCR_malaria %in% 0 ~ 0,
                              TRUE ~ NA_real_),
    
    # ---- Malaria before delivery (ANV1 or ANV2) ----
    malaria_pre_birth = case_when(ANV1_malaria %in% 1 | ANV2_malaria %in% 1 ~ 1,
                                  ANV1_malaria %in% 0 | ANV2_malaria %in% 0 ~ 0,
                                  TRUE ~ NA_real_),
    
    # ---- MiP at any time (main exposure) ----
    # Positive: any positive test (ANV1, ANV2, unscheduled, delivery smear,
    #           peripheral PCR, placental PCR).
    # Negative: no positive test and at least one negative test (as in the Methods).
    any_pos  = ANV1_malaria %in% 1 | ANV2_malaria %in% 1 | n_emerg_pos > 0 |
      ANV3_malaria %in% 1 | PCR_malaria %in% 1 | PCR_placenta %in% 1,
    any_test = !is.na(ANV1_malaria) | !is.na(ANV2_malaria) | (n_emerg_pos + n_emerg_neg) > 0 |
      !is.na(ANV3_malaria) | !is.na(PCR_malaria) | !is.na(PCR_placenta),
    MiP = case_when(any_pos ~ 1, any_test ~ 0, TRUE ~ NA_real_),
    
    # ---- Number of malaria episodes ----

    malaria_dose_old = (ANV1_malaria %in% 1) + (ANV2_malaria %in% 1) + (birth_malaria %in% 1) +
      malaria_emerg + (PCR_malaria %in% 1) + (PCR_placenta %in% 1),
    # PROPOSED definition (see question 1 in the message): one episode per time point.
    delivery_pos     = ANV3_malaria %in% 1 | PCR_malaria %in% 1 | PCR_placenta %in% 1,
    malaria_episodes = (ANV1_malaria %in% 1) + (ANV2_malaria %in% 1) + n_emerg_pos + delivery_pos
  ) %>%
  mutate(
    malaria_episodes_4cat = factor(case_when(MiP == 0 ~ "0",
                                             malaria_episodes == 1 ~ "1",
                                             malaria_episodes == 2 ~ "2",
                                             malaria_episodes >= 3 ~ "3+"),
                                   levels = c("0", "1", "2", "3+")),
    malaria_dose_old_4cat = factor(case_when(malaria_dose_old == 0 ~ "0",
                                             malaria_dose_old == 1 ~ "1",
                                             malaria_dose_old == 2 ~ "2",
                                             malaria_dose_old >= 3 ~ "3+"),
                                   levels = c("0", "1", "2", "3+"))
  )


table(cohort$malaria_dose_old, cohort$malaria_episodes)

# Checks
cat("\n---------------- MALARIA ----------------\n")
for (v in c("ANV1_malaria", "ANV2_malaria", "ANV3_malaria", "PCR_malaria",
            "PCR_placenta", "malaria_emerg", "birth_malaria", "MiP")) {
  cat(v, ":\n"); print(table(cohort[[v]], useNA = "ifany"))
}
cat("Old vs proposed number of episodes:\n")
print(table(old = cohort$malaria_dose_old_4cat, new = cohort$malaria_episodes_4cat, useNA = "ifany"))


# =============================================================================
# 4. COVARIATES
# =============================================================================
cohort <- cohort %>%
  mutate(
    # ---- Mother's age at delivery ----
    age_at_delivery = as.numeric(difftime(as.Date(v05del_delivdate_4m),
                                          as.Date(AME_DateNaiss), units = "days")) / 365.25,
    age_at_delivery = ifelse(age_at_delivery < 10, NA, age_at_delivery),  # missing birthdates
    
    # ---- Gravidity ----
    v01bas_gravidity_cat = factor(v01bas_gravidity_cat,
                                  levels = c("Primigravidae", "1-3 previous", "4 or more")),
    
    # ---- Pre-pregnancy BMI (Ouedraogo 2012: +1 kg per month after 12 weeks) ----
    AME_Taille_m          = AME_Taille / 100,
    estimated_weight_gain = pmax(AME_AgeGesta - 12, 0) / 4.345,
    pre_pregnancy_BMI     = (AME_Poids - estimated_weight_gain) / AME_Taille_m^2,
    
    # ---- IPTp group ----
    v01med_iptpgroup_1 = factor(v01med_iptpgroup_1),
    
    # ---- CRP inflammation at inclusion (> 5 mg/L, Mireku 2016) ----
    CRP_first_visit = if_else(AME_CRP1 > 5, 1, 0),
    
    # ---- Iron deficiency at inclusion (ferritin was stored x10) ----
    iron_divided_1 = AME_Ferritine1 / 10,
    ID_visit1 = if_else(iron_divided_1 < 12 |
                          (CRP_first_visit == 1 & iron_divided_1 >= 12 & iron_divided_1 <= 70), 1, 0),
    
    # ---- Helminth infection at inclusion (Kato-Katz: 1 = neg, 2 = pos, 3 = missing) ----
    # CHANGED: was coded twice in the old script; once is enough.
    helminth_inf1 = case_when(AME_Kato1 == 1 ~ 0, AME_Kato1 == 2 ~ 1, TRUE ~ NA_real_),
    
    # ---- Mother's education at 1 year (0 = none, 1 = any) ----
    mother_edu_1year = case_when(TOV_scol == 0 ~ 0, TOV_scol %in% 1:3 ~ 1, TRUE ~ NA_real_),
    
    # ---- Household wealth index at 1 year ----
    TOV_stot = as.numeric(TOV_stot),
    
    # ---- Preterm birth ----
    preterm_birth = case_when(v05del_delivgestage_1 < 37 ~ 1,
                              v05del_delivgestage_1 >= 37 ~ 0, TRUE ~ NA_real_),
    
    # ---- Birthweight ----
    MSE_v01bir_bweight_5 = as.numeric(MSE_v01bir_bweight_5),
    
    # ---- Child sex ----
    child_sex = case_when(
      MSE_v01bir_sex_4 %in% c("Female", "Male") ~ MSE_v01bir_sex_4,   # birth record first
      EXP_v114 == 1 | TOV_sexe == 0 | DEV_qpa_esexe == 1 ~ "Male",      # then later waves
      EXP_v114 == 2 | TOV_sexe == 1 | DEV_qpa_esexe == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    # ---- Gestational age at ANV1 and ANV2 (for Table 1 / reviewer 1.12) ----
    GA_ANV1 = as.numeric(AME_AgeGesta),
    GA_ANV2 = as.numeric(AME_AgeGestamv1)   # visits table: mv1 = visit 2 (IPTp2)
  )

cat("\n---------------- COVARIATES ----------------\n")
summary(cohort[, c("age_at_delivery", "pre_pregnancy_BMI", "TOV_stot",
                   "MSE_v01bir_bweight_5", "GA_ANV1", "GA_ANV2")])


# Youngest mothers and extreme BMI
cohort %>% filter(age_at_delivery < 15) %>%
  select(mide, age_at_delivery, AME_DateNaiss, v05del_delivdate_4m)

cohort %>% filter(pre_pregnancy_BMI > 40 | pre_pregnancy_BMI < 15) %>%
  select(mide, pre_pregnancy_BMI, AME_Poids, AME_Taille, AME_AgeGesta)

cohort %>%
  filter(mide %in% c(# took out id for github)) %>%
  select(mide, AME_Taille, AME_Taillemv1, AME_Taillemv2,
         AME_Poids, AME_Poidsmv1, AME_Poidsmv2, GA_ANV1, GA_ANV2)


# Inclusion before 14 weeks
sum(cohort$GA_ANV1 < 14, na.rm = TRUE)

# ANV2 implausible: after 40 weeks, or before ANV1
cohort %>% filter(GA_ANV2 > 40 | GA_ANV2 <= GA_ANV1) %>%
  select(mide, GA_ANV1, GA_ANV2, AME_DateVisitmv1, v05del_delivdate_4m)

cohort <- cohort %>%
  mutate(
    gap_weeks = as.numeric(as.Date(AME_DateVisitmv1) - as.Date(AME_DateVisit)) / 7,
    GA_ANV2 = case_when(
      !AME_NumVisitmv1 %in% 2                         ~ NA_real_,                         # record is delivery, not ANV2
      gap_weeks > 0 & gap_weeks <= 20                 ~ GA_ANV1 + gap_weeks,              # dates plausible
      as.numeric(AME_AgeGestamv1) > GA_ANV1           ~ as.numeric(AME_AgeGestamv1),      # dates wrong: use recorded GA
      TRUE                                            ~ NA_real_                          # no usable information
    )
  )

summary(cohort$GA_ANV2)
sum(cohort$GA_ANV2 <= cohort$GA_ANV1, na.rm = TRUE)

# =============================================================================
# 5. SDQ SCORING (complete case)
# =============================================================================

orig_vars <- paste0("EXP_v", 843:867)
items     <- paste0("sdq6ans_", 1:25)

for (i in 1:25) {
  cohort[[items[i]]] <- case_when(cohort[[orig_vars[i]]] == "P" ~ 0,
                                  cohort[[orig_vars[i]]] == "U" ~ 1,
                                  cohort[[orig_vars[i]]] == "T" ~ 2,
                                  TRUE ~ NA_real_)
}

rev_items <- paste0("sdq6ans_", c(7, 11, 14, 21, 25))
cohort <- cohort %>% mutate(across(all_of(rev_items), ~ 2 - .x))

score_scale <- function(df, nums) rowSums(df[paste0("sdq6ans_", nums)], na.rm = FALSE)

cohort <- cohort %>%
  mutate(
    sdq6ans_emotional     = score_scale(., c(3, 8, 13, 16, 24)),
    sdq6ans_conduct       = score_scale(., c(5, 7, 12, 18, 22)),
    sdq6ans_hyperactive   = score_scale(., c(2, 10, 15, 21, 25)),
    sdq6ans_peerproblems  = score_scale(., c(6, 11, 14, 19, 23)),
    sdq6ans_prosocial     = score_scale(., c(1, 4, 9, 17, 20)),
    sdq6ans_internalizing = sdq6ans_emotional + sdq6ans_peerproblems,
    sdq6ans_externalizing = sdq6ans_conduct + sdq6ans_hyperactive,
    sdq6ans_totalsdq      = sdq6ans_internalizing + sdq6ans_externalizing,
    sdq_n_answered        = rowSums(!is.na(across(all_of(items))))
  )

cat("\n---------------- SDQ ----------------\n")
cat("Items answered (assessed children):\n")
print(table(cohort$sdq_n_answered[cohort$i_EXPL %in% 1]))
table(cohort$EXP_v849, cohort$sdq6ans_7, useNA = "ifany")   # P -> 2, T -> 0


# =============================================================================
# 6. TOVA
# =============================================================================
# Two tests with ADHD scores below the possible range (< -10) are invalid
# (on the advice of M. Boivin). Raw D prime > 10 are data-entry errors -> NA.
tova_invalid_ids <- c("B4501", "B4615")

cohort <- cohort %>%
  mutate(
    EXP_dprime  = ifelse(EXP_dprime > 10, NA, EXP_dprime),
    tova_done   = !is.na(EXP_scoreADHD),
    tova_valid  = tova_done & !mide %in% tova_invalid_ids,
    tova_missing_reason = case_when(
      tova_done ~ NA_character_,
      EXP_pasTOVAraison %in% c(14, 18) |
        grepl("DEFICIEN|Attention", EXP_TOVAautreraison, ignore.case = TRUE) |
        grepl("Attention", EXP_dautreraisonpasTOVA, ignore.case = TRUE) ~ "Child-related",
      EXP_pasTOVAraison %in% c(22, 99) |
        grepl("COURANT|MATERIEL", EXP_dautreraisonpasTOVA, ignore.case = TRUE) |
        EXP_TOVAraison %in% 6 ~ "Logistical",
      TRUE ~ "Not recorded"
    )
  )

cat("\n---------------- TOVA ----------------\n")
assessed <- cohort %>% filter(i_EXPL %in% 1)
cat("TOVA done:", sum(assessed$tova_done), "| invalid:", sum(assessed$tova_done & !assessed$tova_valid), "\n")
print(table(assessed$tova_missing_reason, useNA = "ifany"))

# =============================================================================
# 7. ANALYSIS POPULATIONS + FLOWCHART SUMMARY
# =============================================================================
# CHANGED: the TOVA population is now taken from all children assessed at
# EXPLORE (not only those with a valid SDQ), so children with a TOVA but no
# SDQ are included.
Explore  <- cohort  %>% filter(i_EXPL %in% 1)
sdq_pop  <- Explore %>% filter(!is.na(sdq6ans_totalsdq), !is.na(MiP))
tova_pop <- Explore %>% filter(tova_valid, !is.na(MiP))

cat("\n================ FLOWCHART SUMMARY ================\n")
cat("Alive and eligible at 6 years:      ", nrow(cohort), "\n")
cat("  Lost to follow-up (not assessed): ", sum(!cohort$i_EXPL %in% 1), "\n")
cat("Assessed at 6 years (EXPLORE):      ", nrow(Explore), "\n")
cat("\nSDQ branch\n")
cat("  SDQ not completed (0 items):      ", sum(Explore$sdq_n_answered == 0), "\n")
cat("  SDQ incomplete:                   ", sum(Explore$sdq_n_answered > 0 & is.na(Explore$sdq6ans_totalsdq)), "\n")
cat("  Missing MiP:                      ", sum(!is.na(Explore$sdq6ans_totalsdq) & is.na(Explore$MiP)), "\n")
cat("  SDQ ANALYSIS POPULATION:          ", nrow(sdq_pop), "\n")
cat("\nTOVA branch\n")
cat("  TOVA not done:                    ", sum(!Explore$tova_done), "\n")
print(table(Explore$tova_missing_reason))
cat("  TOVA invalid:                     ", sum(Explore$tova_done & !Explore$tova_valid), "\n")
cat("  Missing MiP:                      ", sum(Explore$tova_valid & is.na(Explore$MiP)), "\n")
cat("  TOVA ANALYSIS POPULATION:         ", nrow(tova_pop), "\n")
cat("  ...of whom without a valid SDQ:   ", sum(is.na(tova_pop$sdq6ans_totalsdq)), "\n")


# =============================================================================
# 8. SAVE
# =============================================================================
saveRDS(cohort,   "cohort_alive_eligible_6y.rds")   # attrition table (Supp. Table 8)
saveRDS(sdq_pop,  "malaria_sdq_pop.rds")            # SDQ analyses
saveRDS(tova_pop, "malaria_tova_pop.rds")           # TOVA analyses


cohort=cohort_alive_eligible_6y
sdq_pop=malaria_sdq_pop
tova_pop=malaria_tova_pop



# gt summary

library(dplyr)
library(gtsummary)
library(gt)

# ---- Episodes in 3 categories (0, 1, 2+) ----
add_ep3 <- function(df) df %>%
  mutate(malaria_episodes_3cat = factor(case_when(MiP == 0 ~ "0",
                                                  malaria_episodes == 1 ~ "1",
                                                  malaria_episodes >= 2 ~ "2+"),
                                        levels = c("0", "1", "2+")))
cohort   <- add_ep3(cohort)
sdq_pop  <- add_ep3(sdq_pop)
tova_pop <- add_ep3(tova_pop)

# ---- Variables and labels ----
vars <- c("age_at_delivery", "v01bas_gravidity_cat", "pre_pregnancy_BMI", "mother_edu_1year",
          "v01med_iptpgroup_1", "ID_visit1", "helminth_inf1", "TOV_stot", "CRP_first_visit",
          "child_sex", "preterm_birth", "MSE_v01bir_bweight_5",
          "GA_ANV1", "ANV1_malaria", "GA_ANV2", "ANV2_malaria", "ANV3_malaria", "PCR_malaria",
          "birth_malaria", "malaria_emerg", "PCR_placenta", "MiP", "malaria_episodes_3cat",
          "sdq6ans_internalizing", "sdq6ans_externalizing", "sdq6ans_totalsdq", "EXP_scoreADHD")

labels <- list(
  age_at_delivery       ~ "Mother's age at delivery (years)",
  v01bas_gravidity_cat  ~ "Mother's gravidity",
  pre_pregnancy_BMI     ~ "Mother's pre-pregnancy BMI (kg/m²)",
  mother_edu_1year      ~ "Mother had any formal education",
  v01med_iptpgroup_1    ~ "IPTp treatment group",
  ID_visit1             ~ "Iron deficiency at inclusion",
  helminth_inf1         ~ "Helminth infection at inclusion",
  TOV_stot              ~ "Household wealth index at 1 year",
  CRP_first_visit       ~ "Inflammation (CRP > 5 mg/L) at inclusion",
  child_sex             ~ "Child's sex",
  preterm_birth         ~ "Preterm birth (< 37 weeks)",
  MSE_v01bir_bweight_5  ~ "Birthweight (g)",
  GA_ANV1               ~ "Gestational age at first antenatal visit (weeks)",
  ANV1_malaria          ~ "Malaria at first antenatal visit (smear)",
  GA_ANV2               ~ "Gestational age at second antenatal visit (weeks)",
  ANV2_malaria          ~ "Malaria at second antenatal visit (smear)",
  ANV3_malaria          ~ "Malaria at delivery (smear)",
  PCR_malaria           ~ "Malaria at delivery (peripheral qPCR)",
  birth_malaria         ~ "Malaria at delivery (smear or qPCR)",
  malaria_emerg         ~ "Malaria at an unscheduled visit",
  PCR_placenta          ~ "Placental malaria (qPCR)",
  MiP                   ~ "Malaria in pregnancy (at least once)",
  malaria_episodes_3cat ~ "Number of malaria episodes",
  sdq6ans_internalizing ~ "SDQ internalizing score",
  sdq6ans_externalizing ~ "SDQ externalizing score",
  sdq6ans_totalsdq      ~ "SDQ total difficulties score",
  EXP_scoreADHD         ~ "TOVA ADHD score"
)

outcomes <- c("sdq6ans_internalizing", "sdq6ans_externalizing", "sdq6ans_totalsdq", "EXP_scoreADHD")

make_tbl <- function(df) {
  df %>%
    select(all_of(vars)) %>%
    tbl_summary(
      label = labels,
      type = list(c(mother_edu_1year, ID_visit1, helminth_inf1, CRP_first_visit, preterm_birth,
                    ANV1_malaria, ANV2_malaria, ANV3_malaria, PCR_malaria, birth_malaria,
                    malaria_emerg, PCR_placenta, MiP) ~ "dichotomous",
                  c(TOV_stot, GA_ANV1, GA_ANV2) ~ "continuous"),
      value = list(child_sex ~ "Female"),
      statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                       all_of(outcomes) ~ "{mean} ({sd})",
                       all_categorical() ~ "{n} ({p}%)"),
      digits = list(all_of(outcomes) ~ 2),
      missing = "ifany",
      missing_text = "Missing"
    )
}

# ---- Table 1: eligible cohort | SDQ population | TOVA population ----
table1 <- tbl_merge(
  tbls = list(make_tbl(cohort), make_tbl(sdq_pop), make_tbl(tova_pop)),
  tab_spanner = c(paste0("**Eligible at 6 years**, N = ", nrow(cohort)),
                  paste0("**SDQ population**, N = ", nrow(sdq_pop)),
                  paste0("**TOVA population**, N = ", nrow(tova_pop)))
) %>%
  modify_footnote(everything() ~ NA) %>%
  modify_caption("Median (Q1, Q3) for covariates; mean (SD) for outcome scores; n (%)")

table1
table1 %>% as_gt() %>% gtsave("Table1_populations.docx")



supp8_sdq  <- attrition_tbl(sdq_pop$mide,  "Children included vs not included: SDQ population")
supp8_tova <- attrition_tbl(tova_pop$mide, "Children included vs not included: TOVA population")

supp8_sdq
supp8_tova


library(flextable)

table1 %>%
  as_flex_table() %>%
  save_as_docx(path = "Table1_populations.docx")

supp8_sdq %>%
  as_flex_table() %>%
  save_as_docx(path = "SuppTable8_attrition_SDQ.docx")





