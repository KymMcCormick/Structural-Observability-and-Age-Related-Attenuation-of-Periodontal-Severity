# =============================================================================
# NHANES 2009-2014 PERIODONTAL SURVEY ANALYSIS
# Structural Observability and Age-Related Attenuation of Periodontal Severity
# McCormick KM
# =============================================================================
# Description:
#   Survey-weighted descriptive analysis of NHANES 2009-2014 periodontal data.
#   Constructs age groups, sets up complex survey design, and estimates
#   weighted population characteristics for Table 1.
#
# Input:  nhanes_perio_site — person-level analytic dataset (derived from
#         NHANES periodontal examination files, cycles 2009-2010, 2011-2012,
#         2013-2014). One row per participant. Restricted to dentate adults
#         aged >=30 years who completed the periodontal examination.
#
# Key variables:
#   RIDAGEYR   — age in years at time of examination
#   SDMVPSU    — masked variance pseudo-PSU (primary sampling unit)
#   SDMVSTRA   — masked variance pseudo-stratum
#   WTMEC6YR   — 6-year combined examination weight (constructed by dividing
#                the 2-year weight by 3, per NCHS guidelines for combining
#                3 consecutive 2-year cycles)
#   sex        — character variable ("Male" / "Female")
#
# Packages: survey, dplyr
# =============================================================================

library(survey)
library(dplyr)

options(survey.lonely.psu = "adjust")


unweighted_N <- nrow(nhanes_perio_site)
unweighted_N

nhanes_perio_site <- nhanes_perio_site %>%
  mutate(
    age_group = cut(
      RIDAGEYR,
      breaks = c(30,35,40,45,50,55,60,65,70,75,80, Inf),
      right = FALSE,
      labels = c(
        "30–34","35–39","40–44","45–49","50–54","55–59",
        "60–64","65–69","70–74","75–79","80+"
      )
    )
  )

svy_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  data = nhanes_perio_site,
  nest = TRUE
)

age_dist <- svymean(~age_group, svy_design, na.rm = TRUE)

age_df <- data.frame(
  age_group = gsub("^age_group", "", names(coef(age_dist))),
  percent   = 100 * as.numeric(coef(age_dist)),
  SE        = 100 * as.numeric(SE(age_dist))
)

age_df

mean_age <- svymean(~RIDAGEYR, svy_design, na.rm = TRUE)

mean_age_value <- as.numeric(coef(mean_age))
mean_age_se    <- as.numeric(SE(mean_age))

mean_age_value
mean_age_se

female_est <- svymean(~I(sex == "Female"), svy_design, na.rm = TRUE)

female_pct <- 100 * as.numeric(coef(female_est))
female_se  <- 100 * as.numeric(SE(female_est))

female_pct
female_se

