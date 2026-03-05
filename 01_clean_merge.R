# ------------------------------------------------------------
# Build NHANES 2009–2014 analytic dataset:
# - DEMO (age, sex, survey design vars, MEC weights)
# - OHXPER (site-level CAL = attachment loss variables)
# - Exclude edentulous participants
# Output: nhanes_perio_site (one row per person, wide site-level CAL columns)
# ------------------------------------------------------------

library(haven)
library(dplyr)
library(purrr)

# ---- helpers ----
read_xpt_cycle <- function(path, cycle_letter) {
  read_xpt(path) %>% mutate(source_cycle = cycle_letter)
}

# Update these paths if needed
demo_paths <- c(F = "data/raw/DEMO_F.XPT",
                G = "data/raw/DEMO_G.XPT",
                H = "data/raw/DEMO_H.XPT")

perio_paths <- c(F = "data/raw/OHXPER_F.XPT",
                 G = "data/raw/OHXPER_G.XPT",
                 H = "data/raw/OHXPER_H.XPT")

# ------------------------------------------------------------
# 1) DEMO: keep only required variables
# ------------------------------------------------------------
demo <- imap_dfr(demo_paths, ~ read_xpt_cycle(.x, .y)) %>%
  transmute(
    SEQN,
    source_cycle,
    RIDAGEYR,   # age in years
    RIAGENDR,   # sex (1=Male, 2=Female)
    SDMVPSU,    # PSU
    SDMVSTRA,   # strata
    WTMEC2YR    # 2-year MEC exam weight
  ) %>%
  filter(!is.na(RIDAGEYR), RIDAGEYR >= 30) %>%
  mutate(
    WTMEC6YR = WTMEC2YR / 3,                      # combined 2009–2014 weight
    sex = factor(RIAGENDR, levels = c(1,2), labels = c("Male","Female"))
  )

# ------------------------------------------------------------
# 2) OHXPER: keep site-level CAL variables
#    In NHANES periodontal files, CAL/AL site variables typically look like:
#    OHX##LA(D/M/S/P/L/A) for attachment loss (CAL) at 6 sites per tooth
# ------------------------------------------------------------
ohxper <- imap_dfr(perio_paths, ~ read_xpt_cycle(.x, .y))

# Identify site-level CAL columns programmatically
# This matches: OHX02LAD ... OHX31LAA (6 sites per tooth)
cal_site_vars <- grep("^OHX\\d{2}LA(D|M|S|P|L|A)$", names(ohxper), value = TRUE)

if (length(cal_site_vars) == 0) {
  stop("No CAL site variables detected. Check OHXPER variable names / pattern.")
}

perio_site <- ohxper %>%
  transmute(
    SEQN,
    source_cycle,
    across(all_of(cal_site_vars), ~ as.numeric(.x))
  ) %>%
  # NHANES commonly uses 99 (and sometimes 9/99 etc) as missing codes.
  # For CAL in OHXPER, 99 is the key one.
  mutate(across(all_of(cal_site_vars), ~ na_if(.x, 99)))

# ------------------------------------------------------------
# 3) Exclude edentulous participants
#    For your purposes: edentulous means "no periodontal sites observed"
#    i.e., all CAL sites are NA.
# ------------------------------------------------------------
perio_site <- perio_site %>%
  mutate(
    n_sites_observed = rowSums(!is.na(across(all_of(cal_site_vars)))),
    edentulous = (n_sites_observed == 0)
  ) %>%
  filter(!edentulous) %>%
  select(-edentulous)

# ------------------------------------------------------------
# 4) Merge DEMO + OHXPER into one analytic dataset
# ------------------------------------------------------------
nhanes_perio_site <- demo %>%
  inner_join(perio_site, by = c("SEQN", "source_cycle"))

rm(demo, ohxper, perio_site)

# Optional: sanity checks
message("Rows in final dataset: ", nrow(nhanes_perio_site))
message("Unique participants: ", dplyr::n_distinct(nhanes_perio_site$SEQN))
message("CAL site variables: ", length(cal_site_vars))

# Save (optional)
 saveRDS(nhanes_perio_site, "data/derived/nhanes_perio_site_2009_2014.rds")
 readr::write_csv(nhanes_perio_site, "data/derived/nhanes_perio_site_2009_2014.csv")
 