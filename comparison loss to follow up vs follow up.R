# =============================================================================
# Script:  05_attrition_table.R
# Author:  Katharine Barry
# Purpose: Supplementary table: compare children included in each analysis
#          population with all other eligible children (lost to follow-up or
#          not analysed). 
# =============================================================================
library(dplyr)
library(gtsummary)
library(flextable)

cohort   <- readRDS("cohort_alive_eligible_6y.rds")   # all eligible (796)
sdq_pop  <- readRDS("malaria_sdq_pop.rds")
tova_pop <- readRDS("malaria_tova_pop.rds")

# ---- Continuous variables: medians + Wilcoxon (as in the other tables) ----
# For means instead: cont_stat <- "{mean} ({sd})"; cont_test <- "t.test"
cont_stat <- "{median} ({p25}, {p75})"
cont_test <- "wilcox.test"

vars <- c("age_at_delivery", "v01bas_gravidity_cat", "pre_pregnancy_BMI",
          "mother_edu_1year", "v01med_iptpgroup_1", "ID_visit1", "helminth_inf1",
          "TOV_stot", "CRP_first_visit", "child_sex", "preterm_birth",
          "MSE_v01bir_bweight_5",
          "ANV1_malaria", "ANV2_malaria", "ANV3_malaria", "PCR_malaria",
          "malaria_emerg", "PCR_placenta", "MiP", "malaria_episodes_3cat")

# Labels matching Table 1
labels <- list(
  age_at_delivery       ~ "Mother's age at delivery (years)",
  v01bas_gravidity_cat  ~ "Mother's gravidity",
  pre_pregnancy_BMI     ~ "Mother's estimated pre-pregnancy BMI (kg/m²)",
  mother_edu_1year      ~ "Mother had any formal education",
  v01med_iptpgroup_1    ~ "IPTp treatment group",
  ID_visit1             ~ "Iron deficiency at inclusion",
  helminth_inf1         ~ "Helminth infection at inclusion",
  TOV_stot              ~ "Household wealth index at 1 year",
  CRP_first_visit       ~ "Inflammation (CRP > 5 mg/L) at inclusion",
  child_sex             ~ "Child's sex",
  preterm_birth         ~ "Preterm birth (< 37 weeks)",
  MSE_v01bir_bweight_5  ~ "Birthweight (g)",
  ANV1_malaria          ~ "Malaria at first antenatal visit (smear)",
  ANV2_malaria          ~ "Malaria at second antenatal visit (smear)",
  ANV3_malaria          ~ "Malaria at delivery (smear)",
  PCR_malaria           ~ "Malaria at delivery (peripheral PCR)",
  malaria_emerg         ~ "Malaria at an unscheduled visit",
  PCR_placenta          ~ "Placental malaria (PCR)",
  MiP                   ~ "Malaria in pregnancy (at least once)",
  malaria_episodes_3cat ~ "Number of malaria episodes"
)

dat <- cohort %>%
  mutate(across(where(haven::is.labelled), haven::as_factor),
         malaria_episodes_3cat = factor(case_when(MiP == 0 ~ "0",
                                                  malaria_episodes == 1 ~ "1",
                                                  malaria_episodes >= 2 ~ "2+"),
                                        levels = c("0", "1", "2+")))

summ <- function(d, by = NULL) {
  d %>%
    tbl_summary(by = {{ by }},
                include = all_of(vars),
                label = labels,
                statistic = list(all_continuous()  ~ cont_stat,
                                 all_categorical() ~ "{n} ({p}%)"),
                digits = all_continuous() ~ 1,
                missing = "ifany", missing_text = "Missing")
}

# ---- Included vs not included, per population ----
compare <- function(analysis_pop) {
  dat %>%
    mutate(group = factor(ifelse(mide %in% analysis_pop$mide, "Included", "Not included"),
                          levels = c("Included", "Not included"))) %>%
    summ(by = group) %>%
    # Pearson's chi-square (Fisher's exact when expected counts < 5) and
    # Wilcoxon rank-sum; missing values excluded from the tests
    add_p(test = all_continuous() ~ cont_test,
          pvalue_fun = label_style_pvalue(digits = 2))
}

tab_attrition <- tbl_merge(
  tbls        = list(compare(sdq_pop), compare(tova_pop)),
  tab_spanner = c("**SDQ population**", "**TOVA population**")
) %>%
  bold_labels()

tab_attrition

save_as_docx(
  "Supplementary Table 5: Characteristics of children included in and excluded from the analysis populations" =
    as_flex_table(tab_attrition) %>% autofit(),
  path = "Supp_Table_attrition.docx",
  pr_section = officer::prop_section(page_size = officer::page_size(orient = "landscape"))
)
