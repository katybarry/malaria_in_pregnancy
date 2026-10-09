# =============================================================================
# Script:  04_sensitivity_complete_malaria.R
# Author:  Katharine Barry
# Purpose: looking for invalid ADHD scores
# =============================================================================
library(dplyr)
library(mice)

filter(!is.na(EXP_scoreADHD))
table(malaria_sdq_pop$EXP_TOVAraison, useNA = "always")

#1-Poor compréhension of questions, 2-Attention limitée, 3-Manque de motivation, 4-Impatience et manque de repos, 5-Attitude négative envers l’enquêteur, 6-Interruption pendant le test, 9-Autre, spécifier

table(malaria_sdq_pop$EXP_TOVAautreraison, useNA = "always")


library(dplyr)

malaria_sdq_pop <- malaria_sdq_pop %>%
  mutate(tova_status = case_when(
    !is.na(EXP_scoreADHD) ~ "Score available",
    TRUE ~ "No score"
  ))

# 1. How many have / don't have a score (should be 450 + 55, or 452 + 53?)
table(malaria_sdq_pop$tova_status, useNA = "always")

# 2. Reasons by score status
table(malaria_sdq_pop$EXP_TOVAraison, malaria_sdq_pop$tova_status, useNA = "always")

# 3. Free-text reasons by score status
malaria_sdq_pop %>%
  filter(!is.na(EXP_TOVAraison) | !is.na(EXP_TOVAautreraison)) %>%
  select(EXP_TOVAraison, EXP_TOVAautreraison, tova_status, EXP_scoreADHD) %>%
  print(n = Inf)

# 4. How were the 2 invalid tests identified? List any TOVA validity variables
grep("TOVA|tova|ADHD", names(malaria_sdq_pop), value = TRUE)


# Validity flag and the 2 excluded tests
table(malaria_sdq_pop$EXP_TOVAvalide, malaria_sdq_pop$tova_status, useNA = "always")
malaria_sdq_pop %>%
  filter(tova_status == "Score available", EXP_TOVAvalide != 1 | is.na(EXP_TOVAvalide)) %>%
  select(EXP_TOVAvalide, EXP_scoreADHD, EXP_TOVAraison, EXP_TOVAautreraison)

# Reasons for NOT doing the TOVA, among the 53
table(malaria_sdq_pop$EXP_pasTOVAraison, malaria_sdq_pop$tova_status, useNA = "always")
table(malaria_sdq_pop$EXP_dautreraisonpasTOVA, malaria_sdq_pop$tova_status, useNA = "always")
table(malaria_sdq_pop$EXP_TOVAmiss, malaria_sdq_pop$tova_status, useNA = "always")


malaria_sdq_pop %>%
  filter(tova_status == "No score") %>%
  select(EXP_pasTOVAraison, EXP_dautreraisonpasTOVA, EXP_TOVAraison,
         EXP_TOVAautreraison, EXP_TOVAvalide, EXP_TOVAmiss) %>%
  arrange(EXP_pasTOVAraison) %>%
  print(n = Inf)

# Find the D prime variable name
grep("prime|dprime|Dprime", names(malaria_sdq_pop), value = TRUE, ignore.case = TRUE)

# Replace EXP_dprime below with the right name
malaria_sdq_pop %>%
  filter(!is.na(EXP_scoreADHD)) %>%
  mutate(
    excluded_by_you   = EXP_dprime > 10,
    flagged_invalid   = EXP_TOVAvalide == 0 | is.na(EXP_TOVAvalide)
  ) %>%
  filter(excluded_by_you | flagged_invalid) %>%
  select(EXP_dprime, EXP_scoreADHD, EXP_TOVAvalide, excluded_by_you, flagged_invalid)


malaria_sdq_pop %>%
  filter(mide %in% c("B4501", "B4615") |
           (!is.na(EXP_scoreADHD) & (EXP_TOVAvalide != 1 | is.na(EXP_TOVAvalide)))) %>%
  select(mide, EXP_scoreADHD, EXP_dprime, EXP_dprimestdscore,
         EXP_commission, EXP_omission, EXP_TOVAvalide)

summary(TOVA_data$EXP_dprimestdscore)
