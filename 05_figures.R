# =============================================================================
# NHANES 2009-2014 PERIODONTAL SURVEY ANALYSIS — WEIGHTED MEANS AND FIGURE 2
# Structural Observability and Age-Related Attenuation of Periodontal Severity
# McCormick KM, Amarasena N, Jamieson L
# =============================================================================
# Description:
#   Computes survey-weighted age-specific means and 95% confidence intervals
#   for each periodontal severity outcome, combines with quadratic model
#   predictions from the modelling script, and produces Figure 2.
#
# Input:  nhanes_fig_df — person-level analytic dataset saved as an RDS file.
#         pred_all      — data frame of quadratic model predictions produced
#                         by the modelling script (must be present in the
#                         environment before running this script).
#
# Key variables:
#   perc_sites_ge3     — percentage of examined sites with CAL >=3 mm
#   perc_sites_ge5     — percentage of examined sites with CAL >=5 mm
#   perc_sites_ge6     — percentage of examined sites with CAL >=6 mm
#   mean_CAL           — mean clinical attachment loss (mm) across examined sites
#   pct_sites_observed — percentage of eligible sites with recorded CAL
#   age_group          — ordered factor with 11 levels (30-34 through 80+)
#
# Output:
#   Figure 2 — observed weighted age-specific means (solid line with
#   pointranges) overlaid with survey-weighted quadratic model predictions
#   (dashed line), faceted by severity measure.
#   Saved to: outputs/figures/Figure_observed_NHANES.png
#
# Dependencies:
#   Run modelling script first to generate pred_all in the environment.
#
# Packages: dplyr, survey, ggplot2, tidyr
# =============================================================================

library(dplyr)
library(survey)
library(ggplot2)
library(tidyr)

options(survey.lonely.psu = "adjust")

nhanes_fig_df <- readRDS("data/derived/nhanes_fig_df.rds")

svy_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  data = nhanes_fig_df,
  nest = TRUE
)

get_weighted <- function(varname) {
  out <- svyby(
    as.formula(paste0("~", varname)),
    ~age_group,
    svy_design,
    svymean,
    na.rm = TRUE,
    vartype = "se"
  )
  
  out <- as.data.frame(out)
  
  names(out)[names(out) == varname] <- "mean"
  names(out)[names(out) == "se"] <- "se"
  
  out$measure <- varname
  
  out %>%
    mutate(
      lwr = mean - 1.96 * se,
      upr = mean + 1.96 * se
    )
}

wt_ge3  <- get_weighted("perc_sites_ge3")
wt_ge5  <- get_weighted("perc_sites_ge5")
wt_ge6  <- get_weighted("perc_sites_ge6")
wt_obs  <- get_weighted("pct_sites_observed")
wt_mean <- get_weighted("mean_CAL")

wt_all <- bind_rows(wt_ge3, wt_ge5, wt_ge6, wt_obs, wt_mean)

wt_all <- wt_all %>%
  mutate(
    measure = recode(measure,
                     perc_sites_ge3     = "% sites ≥3 mm",
                     perc_sites_ge5     = "% sites ≥5 mm",
                     perc_sites_ge6     = "% sites ≥6 mm",
                     pct_sites_observed = "% sites observed",
                     mean_CAL           = "Mean CAL (mm)"
    ),
    age_group = factor(age_group,
                       levels = levels(nhanes_fig_df$age_group),
                       ordered = TRUE)
  )

p <- ggplot() +
  
  # Quadratic fitted curve (dashed)
  geom_line(
    data = pred_all,
    aes(x = age_group, y = fit, group = 1),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  
  # Weighted observed means (solid)
  geom_line(
    data = wt_all,
    aes(x = age_group, y = mean, group = 1),
    linewidth = 0.9
  ) +
  
  # Observed 95% CI
  geom_pointrange(
    data = wt_all,
    aes(x = age_group, y = mean, ymin = lwr, ymax = upr),
    linewidth = 0.4
  ) +
  
  facet_wrap(~measure, scales = "free_y", ncol = 2) +
  labs(x = "Age group (years)", y = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(face = "bold")
  )

print(p)

ggsave(
  "outputs/figures/Figure_observed_NHANES.png",
  p,
  width = 7,
  height = 10,
  dpi = 300
)

