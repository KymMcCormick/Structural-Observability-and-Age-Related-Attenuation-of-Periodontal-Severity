# =============================================================================
# NHANES 2009-2014 PERIODONTAL SURVEY ANALYSIS — QUADRATIC MODELS AND
# MONOTONICITY TESTS
# Structural Observability and Age-Related Attenuation of Periodontal Severity
# McCormick KM
# =============================================================================
# Description:
#   Fits survey-weighted quadratic regression models for each periodontal
#   severity outcome and tests for non-monotonicity in the CAL >=6mm
#   trajectory using the delta method and a pairwise adjacent-group contrast.
#
# Input:  nhanes_fig_df — person-level analytic dataset with person-level
#         summary measures. One row per participant. age_group is an ordered
#         factor with 11 levels (30-34 through 80+); age_index is its numeric
#         equivalent (1-11).
#
# Key variables:
#   perc_sites_ge3     — percentage of examined sites with CAL >=3 mm
#   perc_sites_ge5     — percentage of examined sites with CAL >=5 mm
#   perc_sites_ge6     — percentage of examined sites with CAL >=6 mm
#   mean_CAL           — mean clinical attachment loss (mm) across examined sites
#   pct_sites_observed — percentage of eligible sites with recorded CAL
#   age_index          — continuous age group index (1=30-34, ..., 11=80+)
#
# Outputs:
#   pred_all   — data frame of quadratic model predictions for Figure 2
#   F-test results for linear vs quadratic model comparison (printed)
#   Delta method test: predicted peak vs age 80+ for CAL >=6mm (printed)
#   Pairwise contrast: age groups 55-59 vs 60-64 for CAL >=6mm (printed)
#
# Packages: dplyr, survey
# =============================================================================

library(dplyr)
library(survey)

options(survey.lonely.psu = "adjust")

nhanes_fig_df <- nhanes_fig_df %>%
  mutate(age_index = as.numeric(age_group))

svy_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  data = nhanes_fig_df,
  nest = TRUE
)

m3 <- svyglm(perc_sites_ge3 ~ poly(age_index, 2, raw = TRUE), design = svy_design)
m5   <- svyglm(perc_sites_ge5 ~ poly(age_index, 2, raw = TRUE), design = svy_design)
m6   <- svyglm(perc_sites_ge6 ~ poly(age_index, 2, raw = TRUE), design = svy_design)
mCAL <- svyglm(mean_CAL ~ poly(age_index, 2, raw = TRUE), design = svy_design)
mOBS <- svyglm(pct_sites_observed ~ poly(age_index, 2, raw = TRUE), design = svy_design)

newdat <- data.frame(
  age_index = seq_along(levels(nhanes_fig_df$age_group))
)

newdat$age_group <- factor(
  levels(nhanes_fig_df$age_group),
  levels = levels(nhanes_fig_df$age_group),
  ordered = TRUE
)

predict_model <- function(model, label) {
  
  pr <- predict(model, newdat, se.fit = TRUE)
  
  fit <- as.numeric(pr)
  se  <- sqrt(as.numeric(attr(pr, "var")))
  
  data.frame(
    age_group = newdat$age_group,
    measure = label,
    fit = fit,
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se
  )
}

pred_all <- bind_rows(
  predict_model(m3,   "% sites ≥3 mm"),
  predict_model(m5,   "% sites ≥5 mm"),
  predict_model(m6,   "% sites ≥6 mm"),
  predict_model(mCAL, "Mean CAL (mm)"),
  predict_model(mOBS, "% sites observed")
)

anova(svyglm(perc_sites_ge3 ~ age_index, design = svy_design), m3)
anova(svyglm(perc_sites_ge5 ~ age_index, design = svy_design), m5)
anova(svyglm(perc_sites_ge6 ~ age_index, design = svy_design), m6)
anova(svyglm(mean_CAL ~ age_index, design = svy_design), mCAL)

# ------------------------
# MONOTONICITY TEST FOR >=6mm CAL
# Delta method test: predicted peak vs predicted value at age index 11 (80+)
# ------------------------

# Extract coefficients from survey-weighted quadratic model
b <- coef(m6)
b0 <- b["(Intercept)"]
b1 <- b["poly(age_index, 2, raw = TRUE)1"]
b2 <- b["poly(age_index, 2, raw = TRUE)2"]

# Age index at predicted peak (vertex of quadratic)
age_peak <- -b1 / (2 * b2)
cat("Predicted peak at age index:", round(age_peak, 2), "\n")
cat("Approximate age group:", levels(nhanes_fig_df$age_group)[round(age_peak)], "\n")

# Predicted values at peak and at age index 11 (80+)
pred_peak <- b0 + b1 * age_peak + b2 * age_peak^2
pred_80   <- b0 + b1 * 11       + b2 * 11^2

cat("Predicted % ge6 at peak:", round(pred_peak, 3), "\n")
cat("Predicted % ge6 at age 80+:", round(pred_80, 3), "\n")
cat("Difference (peak - age 80+):", round(pred_peak - pred_80, 3), "\n")

# Delta method SE
# gradient of (pred_peak - pred_80) with respect to (b0, b1, b2)
# intercept terms cancel: d/db0 = 0
# d/db1 = age_peak - 11
# d/db2 = age_peak^2 - 121
grad <- c(0, age_peak - 11, age_peak^2 - 121)

V        <- vcov(m6)
se_diff  <- as.numeric(sqrt(t(grad) %*% V %*% grad))
diff_est <- as.numeric(pred_peak - pred_80)
z_stat   <- diff_est / se_diff
p_value  <- 2 * pnorm(-abs(z_stat))

cat("\nDelta method test: predicted peak vs age 80+\n")
cat("Difference:", round(diff_est, 3), "percentage points\n")
cat("SE:", round(se_diff, 3), "\n")
cat("Z statistic:", round(z_stat, 3), "\n")
cat("P value:", round(p_value, 4), "\n")
cat("95% CI:", round(diff_est - 1.96 * se_diff, 3),
    "to", round(diff_est + 1.96 * se_diff, 3), "\n")

# ------------------------
# PAIRWISE CONTRAST: 55-59 (age index 6) vs 60-64 (age index 7)
# ------------------------

pred_55 <- b0 + b1 * 6 + b2 * 36
pred_60 <- b0 + b1 * 7 + b2 * 49

diff_55_60 <- as.numeric(pred_55 - pred_60)

grad_55_60 <- c(0, 6 - 7, 36 - 49)
se_55_60   <- as.numeric(sqrt(t(grad_55_60) %*% V %*% grad_55_60))
z_55_60    <- diff_55_60 / se_55_60
p_55_60    <- 2 * pnorm(-abs(z_55_60))

cat("Contrast: 55-59 vs 60-64\n")
cat("Difference:", round(diff_55_60, 3), "percentage points\n")
cat("SE:", round(se_55_60, 3), "\n")
cat("Z statistic:", round(z_55_60, 3), "\n")
cat("P value:", round(p_55_60, 4), "\n")
cat("95% CI:", round(diff_55_60 - 1.96 * se_55_60, 3),
    "to", round(diff_55_60 + 1.96 * se_55_60, 3), "\n")
# Flip for intuitive reporting (60-64 minus 55-59)
cat("Decline from 55-59 to 60-64:", round(-diff_55_60, 3), "percentage points\n")
