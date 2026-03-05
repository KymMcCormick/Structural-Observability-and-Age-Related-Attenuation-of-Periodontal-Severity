# =============================================================================
# SIMULATION OF LATENT PERIODONTAL DISEASE SEVERITY AND STRUCTURAL TRUNCATION
# Structural Observability and Age-Related Attenuation of Periodontal Severity
# McCormick KM, Amarasena N, Jamieson L
# =============================================================================
# Description:
#   Monte Carlo simulation of latent periodontal disease (clinical attachment
#   loss, CAL) accumulating monotonically with age across 28 anatomical
#   subunits (teeth), under three observability mechanisms:
#
#     Model 1 (full)           — no tooth loss; all teeth observed
#     Model 2 (noninformative) — tooth loss depends on age only
#     Model 3 (informative)    — tooth loss depends on age and disease
#                                severity, with severity-selectivity
#                                increasing with age
#
#   Compares latent burden with observed burden among retained teeth across
#   age groups to quantify distortion induced by each mechanism. Includes
#   a simulation-based monotonicity test for the >=6mm trajectory and
#   produces Figure 1 (panels A and B) for the manuscript.
#
# Design:
#   10,000 individuals per age group x 11 age groups (30-34 through 80+)
#   x 28 teeth = 3,080,000 tooth-level observations per simulation run.
#   Three runs (full, noninformative, informative). set.seed(1) ensures
#   reproducibility.
#
# Key parameters (primary simulation):
#   base_mu   = 1.5    — baseline CAL (mm)
#   age_slope = 0.4    — age progression coefficient
#   sd_id     = 1.5    — person-level random effect SD
#   sd_tooth  = 1.0    — tooth-level random effect SD
#   sd_eps    = 1.0    — residual SD
#   cal_cap   = 12     — upper bound for CAL (mm)
#   a3        = 0.25   — severity penalty multiplier for CAL >=3mm
#   a6        = 1.00   — severity penalty multiplier for CAL >=6mm
#   w0/w1/w2  = 0.05/0.04/0.004 — age-strengthening weight coefficients
#   p_floor   = 0.06   — minimum retention probability
#   r0        = 2.84   — baseline retention intercept (log-odds)
#   r_age     = -0.12  — age coefficient for baseline retention (log-odds)
#   Baseline retention calibrated to NHANES: 93% at age 30-34, 68% at 80+
#
# Sensitivity analyses (separate runs, Supplementary Material):
#   Flat age interaction: w0=0.10, w1=0, w2=0
#   Weak penalty:         a3=0.125, a6=0.50
#   Strong penalty:       a3=0.50,  a6=2.00
#
# Outputs:
#   Figure 1A — latent vs observed trajectories under informative truncation
#   Figure 1B — observed trajectories under all three truncation mechanisms
#               vs latent reference
#   Saved to: outputs/figures/Figure1_AB.png
#   Monotonicity test results printed to console
#
# Packages: dplyr, tidyr, ggplot2, patchwork
# =============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

set.seed(1) # reproducible simulation
logit <- function(x) 1/(1 + exp(-x)) # the inverse-logit mapping real numbers to probabilities in (0,1)

# ---------------------------
# PRIMARY SIMULATION FUNCTION
#
# Inputs latent disease model, incluing baseline CAL level (base_mu), age effect before susceptibility scaling (age_slope),
# person-level heterogeneity (sd_id), tooth-level heterogeneity (sd_tooth), residual noise (sd_eps), cap CAL at 10mm after truncating below 0 (cal_cap)
#
# ---------------------------

simulate_latent_linear_observed_concave_tooth <- function(
    n_per_age = 20000,
    n_teeth = 28,
    age_groups = 1:11, # five year bins from 30 to 80+
    
    # -------------------------
    # LATENT disease (stable slope)
    # -------------------------
    base_mu   = 1.5,
    age_slope = 0.4,   # smaller slope to avoid global saturation
    sd_id     = 1.5,
    sd_tooth  = 1,
    sd_eps    = 1,
    cal_cap   = 12,
    
    # -------------------------
    # Tooth susceptibility (key idea)
    # -------------------------
    # weights multiply the age-driven progression term; 0.2 = "rarely affected", 1.0 = "highly affected"
    susc_low  = 0.20,
    susc_mid  = 0.55,
    susc_high = 1.00, 
    
    # pattern: 28-tooth mouth excluding 3rd molars
    # We’ll encode a plausible pattern:
    # - molars: high
    # - lower incisors: mid-high (more susceptible)
    # - upper incisors/canines/premolars: low
    
    susc_vec = NULL,
    
    trunc_mode = c("full", "noninformative", "informative"),
    
    # -------------------------
    # Baseline retention (age-only)
    # -------------------------
    r0   = 2.84,
    r_age = -0.12,
    
    # -------------------------
    # Informative retention: selectivity against disease increases with age
    # -------------------------
    # Two-band penalty so you can get a gap at >=3 and bigger gap at >=6
    a3 = 0.25,      # moderate penalty multiplier (>=3)
    a6 = 1.00,      # severe penalty multiplier (>=6)
    w0 = 0.05,      # The penalty strengthens with age through w_age = w0 + w1*age + w2*age^2
    w1 = 0.04,
    w2 = 0.004,
    p_floor = 0.06  # prevents retention probability from going essentially to 0 in extreme cases
) {
  trunc_mode <- match.arg(trunc_mode)
  
  # default susceptibility pattern 
  # Here: 8 molars (high), 4 lower incisors (mid), everything else low.
  if (is.null(susc_vec)) {
    susc_vec <- rep(susc_low, n_teeth)
    molar_idx <- c(1, 2, 13, 14, 15, 16, 27, 28)
    susc_vec[molar_idx] <- susc_high
    li_idx <- c(20, 21, 22, 23)
    susc_vec[li_idx] <- susc_mid
  } else {
    stopifnot(length(susc_vec) == n_teeth)
  }
  
  # persons balanced by age
  df_id <- expand_grid(age_index = age_groups, k = 1:n_per_age) %>%
    mutate(
      id  = row_number(),
      u_i = rnorm(n(), 0, sd_id)
    ) %>%
    select(id, age_index, u_i)
  
  # person × tooth table
  df_tooth <- df_id %>%
    uncount(n_teeth, .id = "tooth") %>%
    mutate(
      mode = trunc_mode,
      
      # tooth susceptibility weight
      susc = susc_vec[tooth],
      
      # random tooth + residual
      v_it = rnorm(n(), 0, sd_tooth),
      eps  = rnorm(n(), 0, sd_eps),
      
      # -------------------------
      # LATENT CAL: stable slope, but progression depends on susceptibility
      # -------------------------
      # baseline component: everyone can have some CAL
      # progression component: scaled by tooth susceptibility (so only some teeth really progress)
      CAL_lat_raw = base_mu + (age_slope * age_index * susc) + u_i + v_it + eps,  # multilevel linear model with heterogeneity
      CAL_lat     = pmin(pmax(0, CAL_lat_raw), cal_cap), # apply truncation to prevent negative CAL or unrealistic huge CAL
      
      # -------------------------
      # baseline retention (other-cause)
      # -------------------------
      p_base = logit(r0 + r_age * age_index), # the probability a tooth is retained if loss depends only on age
      
      # -------------------------
      # informative retention: penalise disease, more strongly with age
      # -------------------------
      ex3 = pmax(0, CAL_lat - 3), # moderate penalty
      ex6 = pmax(0, CAL_lat - 6), # stronger penalty
      w_age = (w0 + w1 * age_index + w2 * (age_index^2)), # age-strengthened penalty for extraction of high CAL tooth
      sev_pen = -(w_age * (a3 * ex3 + a6 * ex6)), # log-odds penalty for tooth extraction
      
      p_ret = case_when(
        mode == "full" ~ 1,
        mode == "noninformative" ~ p_base,
        mode == "informative" ~ logit(qlogis(p_base) + sev_pen) # converts baseline probability to log-odds, add severity penalty, 
                                                                # then convert back to probability
                                                                # Severity modifies the log-odds of retention, not the probability directly
                                                                # keeps everything mathematically coherent and bounded in (0,1)
      ),
      p_ret = if_else(mode == "informative", pmax(p_ret, p_floor), p_ret), # Apply floor to prevent very severe teeth, too strong penalty etc
      R = rbinom(n(), 1, p_ret),            # realised observation: draw 1 with probability p_ret
      CAL_obs = ifelse(R == 1, CAL_lat, NA_real_), # if retained, then observed CAL, otherwise missing
  # --- decompose loss mechanisms on the probability scale ---
  p_inf = logit(qlogis(p_base) + sev_pen),  # informative retention probability (before floor)
  p_inf = pmax(p_inf, p_floor),             # apply same floor you use for informative
  
  loss_other = 1 - p_base,                  # "other reasons" loss
  loss_inf   = p_base - p_inf,              # additional loss attributable to severity selection
  loss_total = 1 - p_inf,                   # total loss under informative mechanism
  prop_inf_among_lost = loss_inf / loss_total,
  
  prop_loss_inf = loss_inf / pmax(loss_total, 1e-12)  # fraction of total loss that is informative
    )
  
  # person summaries
  person <- df_tooth %>%
    group_by(id, age_index, mode) %>%
    summarise(
      pct_obs = 100 * mean(R == 1),
      
      perc_ge3_lat = mean(CAL_lat >= 3) * 100,
      perc_ge5_lat = mean(CAL_lat >= 5) * 100,
      perc_ge6_lat = mean(CAL_lat >= 6) * 100,
      
      perc_ge3_obs = mean(CAL_obs >= 3, na.rm = TRUE) * 100,
      perc_ge5_obs = mean(CAL_obs >= 5, na.rm = TRUE) * 100,
      perc_ge6_obs = mean(CAL_obs >= 6, na.rm = TRUE) * 100,
      
      mean_CAL_lat = mean(CAL_lat),
      mean_CAL_obs = mean(CAL_obs, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  list(tooth = df_tooth, person = person)
}


# ------------------------
# Run 3 modes
# ------------------------
sim_full   <- simulate_latent_linear_observed_concave_tooth(trunc_mode = "full",           n_per_age = 10000)
sim_noninf <- simulate_latent_linear_observed_concave_tooth(trunc_mode = "noninformative", n_per_age = 10000)
sim_inf    <- simulate_latent_linear_observed_concave_tooth(trunc_mode = "informative",    n_per_age = 10000)

person_all <- bind_rows(sim_full$person, sim_noninf$person, sim_inf$person) %>%
  mutate(
    mode = factor(mode,
                  levels = c("informative","noninformative","full"),
                  labels = c("Informative","Age-only","None"))
  )

age_summary <- person_all %>%
  group_by(mode, age_index) %>%
  summarise(
    pct_obs = mean(pct_obs),
    mean_CAL_lat = mean(mean_CAL_lat),
    mean_CAL_obs = mean(mean_CAL_obs, na.rm = TRUE),
    perc_ge3_lat = mean(perc_ge3_lat),
    perc_ge3_obs = mean(perc_ge3_obs, na.rm = TRUE),
    perc_ge5_lat = mean(perc_ge5_lat),
    perc_ge5_obs = mean(perc_ge5_obs, na.rm = TRUE),
    perc_ge6_lat = mean(perc_ge6_lat),
    perc_ge6_obs = mean(perc_ge6_obs, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------
# MONOTONICITY TEST
# Fit quadratic to observed >=6mm trajectory
# Test whether predicted value at peak is significantly higher than at age 11
# Uses delta method for SE of the difference
# ------------------------

# Aggregate to age-group means for observed >=6mm under informative truncation
ge6_obs <- age_summary %>%
  filter(mode == "Informative") %>%
  select(age_index, perc_ge6_obs)

# Fit quadratic model
mod_ge6 <- lm(perc_ge6_obs ~ age_index + I(age_index^2), data = ge6_obs)

# Extract coefficients
b <- coef(mod_ge6)
b0 <- b["(Intercept)"]
b1 <- b["age_index"]
b2 <- b["I(age_index^2)"]

# Age index at predicted peak (vertex of quadratic)
age_peak <- -b1 / (2 * b2)
cat("Predicted peak at age index:", round(age_peak, 2), "\n")
cat("Approximate age group:", age_labs[round(age_peak)], "\n")

# Predicted values at peak and at age index 11 (80+)
pred_peak <- b0 + b1 * age_peak + b2 * age_peak^2
pred_80   <- b0 + b1 * 11       + b2 * 11^2

cat("Predicted % ge6 at peak:", round(pred_peak, 3), "\n")
cat("Predicted % ge6 at age 80+:", round(pred_80, 3), "\n")
cat("Difference (peak - age 80+):", round(pred_peak - pred_80, 3), "\n")

# Delta method SE for (predicted at peak - predicted at age 80)
# gradient of difference with respect to (b0, b1, b2)
# d/db0 = 1 - 1 = 0
# d/db1 = age_peak - 11
# d/db2 = age_peak^2 - 11^2
grad <- c(0, age_peak - 11, age_peak^2 - 121)

V        <- vcov(mod_ge6)
se_diff  <- sqrt(t(grad) %*% V %*% grad)
diff_est <- pred_peak - pred_80
t_stat   <- diff_est / se_diff
p_value  <- 2 * pt(-abs(t_stat), df = nrow(ge6_obs) - 3)

cat("\nDelta method test: predicted peak vs age 80+\n")
cat("Difference:", round(diff_est, 3), "\n")
cat("SE:", round(se_diff, 3), "\n")
cat("t statistic:", round(t_stat, 3), "\n")
cat("p value:", round(p_value, 4), "\n")
cat("95% CI:", round(diff_est - 1.96 * se_diff, 3),
    "to", round(diff_est + 1.96 * se_diff, 3), "\n")

# ------------------------
# Plot Informative only (easier to debug)
# ------------------------
age_labs <- c("30–34","35–39","40–44","45–49","50–54","55–59",
              "60–64","65–69","70–74","75–79","80+")

plot_df <- age_summary %>%
  pivot_longer(cols = -c(mode, age_index), names_to = "name", values_to = "value") %>%
  mutate(
    type = case_when(
      grepl("_lat$", name) ~ "latent",
      grepl("_obs$", name) ~ "observed",
      name == "pct_obs"    ~ "observed",
      TRUE ~ NA_character_
    ),
    metric = case_when(
      name == "pct_obs"              ~ "% teeth observed",
      grepl("^perc_ge3_", name)      ~ "% teeth ≥3mm",
      grepl("^perc_ge5_", name)      ~ "% teeth ≥5mm",
      grepl("^perc_ge6_", name)      ~ "% teeth ≥6mm",
      grepl("^mean_CAL_", name)      ~ "Mean CAL (mm)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(type), !is.na(metric)) %>%
  mutate(
    metric = factor(metric, levels = c("% teeth ≥3mm","% teeth ≥5mm","% teeth ≥6mm",
                                       "Mean CAL (mm)","% teeth observed")),
    type = factor(type, levels = c("latent","observed"))
  )

p1 <- ggplot(plot_df %>% filter(mode == "Informative"),
       aes(age_index, value, linetype = type, group = type)) +
  geom_line(linewidth = 0.9) +
  facet_grid(metric ~ ., scales = "free_y") +
  scale_x_continuous(breaks = 1:11, labels = age_labs) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = "Age group (years)", y = NULL, linetype = NULL)

loss_age <- sim_inf$tooth %>%
  group_by(age_index) %>%
  summarise(
    loss_other = mean(1 - p_base),
    loss_inf   = mean(p_base - p_inf),
    loss_total = mean(1 - p_inf),
    .groups = "drop"
  ) %>%
  mutate(
    prop_inf_among_lost = loss_inf / pmax(loss_total, 1e-12)
  )

ggplot(loss_age, aes(age_index, 100 * prop_inf_among_lost)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = 1:11, labels = age_labs) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = "Age group (years)",
       y = "Informative component of extractions (%)")

loss_age <- sim_inf$tooth %>%   # informative run only
  group_by(age_index) %>%
  summarise(
    loss_other = mean(loss_other),
    loss_inf   = mean(loss_inf),
    loss_total = mean(loss_total),
    prop_loss_inf = mean(prop_loss_inf),
    .groups = "drop"
  ) %>%
  mutate(
    pct_other = 100 * loss_other,
    pct_inf   = 100 * loss_inf,
    pct_total = 100 * loss_total,
    pct_prop_inf = 100 * prop_loss_inf
  )

loss_long <- loss_age %>%
  select(age_index, pct_other, pct_inf) %>%
  pivot_longer(-age_index, names_to = "component", values_to = "pct") %>%
  mutate(
    component = recode(component,
                       pct_other = "Other reasons (age-only)",
                       pct_inf   = "Informative (severity-related)")
  )

ggplot(loss_long, aes(age_index, pct, fill = component)) +
  geom_col() +
  scale_x_continuous(breaks = 1:11, labels = age_labs) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = "Age group (years)", y = "Percentage of teeth lost", fill = NULL)

ggplot(loss_age, aes(age_index, pct_prop_inf)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = 1:11, labels = age_labs) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = "Age group (years)",
       y = "Informative component of total loss (%)")

trunc_df <- age_summary %>%
  filter(mode != "None") %>%
  select(mode, age_index, pct_obs) %>%
  mutate(
    mode = factor(mode,
                  levels = c("Informative","Age-only"))
  )

ggplot(trunc_df,
       aes(age_index, pct_obs, linetype = mode)) +
  geom_line(linewidth = 1.1) +
  scale_x_continuous(breaks = 1:11, labels = age_labs) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = "Age group (years)",
       y = "Percentage of teeth observed",
       color = "Truncation mechanism")


plot_long <- age_summary %>%
  pivot_longer(cols = -c(mode, age_index), names_to = "name", values_to = "value") %>%
  mutate(
    type = case_when(
      grepl("_lat$", name) ~ "Latent (pre-truncation)",
      grepl("_obs$", name) ~ "Observed (among retained)",
      name == "pct_obs"    ~ "Observed (among retained)",
      TRUE ~ NA_character_
    ),
    metric = case_when(
      name == "pct_obs"              ~ "% teeth observed",
      grepl("^perc_ge3_", name)      ~ "% teeth with CAL ≥3 mm",
      grepl("^perc_ge5_", name)      ~ "% teeth with CAL ≥5 mm",
      grepl("^perc_ge6_", name)      ~ "% teeth with CAL ≥6 mm",
      grepl("^mean_CAL_", name)      ~ "Mean CAL (mm)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(type), !is.na(metric)) %>%
  mutate(
    metric = factor(metric, levels = c("% teeth with CAL ≥3 mm",
                                       "% teeth with CAL ≥5 mm",
                                       "% teeth with CAL ≥6 mm",
                                       "Mean CAL (mm)",
                                       "% teeth observed")),
    mode = factor(mode, levels = c("None","Age-only","Informative"))
  )

fig1_df <- plot_long %>%
  filter(mode == "Informative") %>%
  mutate(type = factor(type, levels = c("Latent (pre-truncation)",
                                        "Observed (among retained)")))

p2 <- ggplot(fig1_df, aes(age_index, value, linetype = type, group = type)) +
  geom_line(linewidth = 0.9) +
  facet_grid(metric ~ ., scales = "free_y") +
  scale_x_continuous(breaks = 1:11, labels = age_labs) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = "Age group (years)", y = NULL, linetype = NULL)


latent_ref <- plot_long %>%
  filter(type == "Latent (pre-truncation)", mode == "None") %>%
  select(age_index, metric, latent_value = value)

efig_df <- plot_long %>%
  filter(type == "Observed (among retained)") %>%
  left_join(latent_ref, by = c("age_index","metric")) %>%
  mutate(mode = factor(mode, levels = c("None","Age-only","Informative")))

p_efig <- ggplot(efig_df, aes(x = age_index, y = value, linetype = mode, group = mode)) +
  geom_line(linewidth = 0.9) +
  geom_line(
    data = distinct(efig_df, age_index, metric, latent_value),
    aes(x = age_index, y = latent_value, group = 1),
    linewidth = 0.9,
    linetype = "solid",
    inherit.aes = FALSE
  ) +
  facet_grid(metric ~ ., scales = "free_y") +
  scale_x_continuous(breaks = 1:11, labels = age_labs) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = "Age group (years)", y = NULL, linetype = "Truncation mechanism")

p_efig


library(patchwork)

# Panel A (latent vs observed)
pA <- p1 +
  scale_linetype_manual(values = c("latent" = "solid", "observed" = "dashed")) +
  guides(linetype = guide_legend(title = NULL, ncol = 1)) +
  theme(legend.position = "bottom")

# Panel B (truncation mechanism)
pB <- p_efig +
  scale_linetype_manual(values = c("None" = "solid",
                                   "Age-only" = "dotted",
                                   "Informative" = "dashed")) +
  guides(linetype = guide_legend(title = NULL, nrow = 1)) +
  theme(legend.position = "bottom")

pB <- p_efig +
  scale_linetype_manual(
    values = c("None" = "solid",
               "Age-only" = "dotted",
               "Informative" = "dashed"),
    labels = c("None" = "No truncation",
               "Age-only" = "Background (non-perio) loss",
               "Informative" = "Informative (severity-dependent)")
  ) +
  guides(linetype = guide_legend(title = NULL, ncol = 1)) +
  theme(legend.position = "bottom")

combined <- (pA + pB) +
  plot_layout(ncol = 2) +                 # <-- no guides="collect"
  plot_annotation(tag_levels = "A")

combined

ggsave(
  "outputs/figures/Figure1_AB strong multipliers.png",
  combined,
  width = 7,
  height = 10,
  dpi = 300
)

