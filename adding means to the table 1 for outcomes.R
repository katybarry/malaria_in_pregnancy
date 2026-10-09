library(dplyr)
library(tidyr)

describe <- function(x) {
  c(n      = sum(!is.na(x)),
    mean   = mean(x, na.rm = TRUE),
    sd     = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q1     = unname(quantile(x, 0.25, na.rm = TRUE)),
    q3     = unname(quantile(x, 0.75, na.rm = TRUE)),
    min    = min(x, na.rm = TRUE),
    max    = max(x, na.rm = TRUE))
}

# SDQ scores (N = 505)
sdq_vars <- c("sdq6ans_internalizing", "sdq6ans_externalizing", "sdq6ans_totalsdq")
round(sapply(malaria_sdq_pop[sdq_vars], function(x) describe(as.numeric(x))), 2)

# TOVA ADHD score (N = 450)
round(describe(as.numeric(TOVA_data$EXP_scoreADHD)), 2)

# Check: total SDQ should equal internalizing + externalizing
malaria_sdq_pop %>%
  mutate(diff = as.numeric(sdq6ans_totalsdq) -
           (as.numeric(sdq6ans_internalizing) + as.numeric(sdq6ans_externalizing))) %>%
  count(diff)


mean(TOVA_data$EXP_scoreADHD < -1.80) * 100
