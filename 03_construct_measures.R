# code/02_construct_measures.R
# Construct person-level periodontal measures FROM site-level CAL
# using the unified dataset nhanes_perio_site.

source("00_packages.R")
dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)

# ---- Identify site-level CAL columns ----
cal_site_vars <- grep("^OHX\\d{2}LA(D|M|S|P|L|A)$", names(nhanes_perio_site), value = TRUE)
if (length(cal_site_vars) == 0) {
  stop("No site-level CAL variables found (expected OHX##LA[DMSPLA]).")
}

# ---- Construct measures (person-level summaries of site-level CAL) ----
nhanes_measures <- nhanes_perio_site %>%
  dplyr::mutate(
    # count observed sites
    n_sites_observed = rowSums(!is.na(dplyr::across(dplyr::all_of(cal_site_vars)))),
    
    # denominator for % observed
    pct_sites_observed = 100 * n_sites_observed / length(cal_site_vars),
    
    # mean CAL across observed sites
    mean_CAL = dplyr::if_else(
      n_sites_observed > 0,
      rowMeans(dplyr::across(dplyr::all_of(cal_site_vars)), na.rm = TRUE),
      NA_real_
    ),
    
    # % sites above thresholds (computed among observed sites)
    perc_sites_ge3 = dplyr::if_else(
      n_sites_observed > 0,
      100 * rowMeans(dplyr::across(dplyr::all_of(cal_site_vars)) >= 3, na.rm = TRUE),
      NA_real_
    ),
    perc_sites_ge5 = dplyr::if_else(
      n_sites_observed > 0,
      100 * rowMeans(dplyr::across(dplyr::all_of(cal_site_vars)) >= 5, na.rm = TRUE),
      NA_real_
    ),
    perc_sites_ge6 = dplyr::if_else(
      n_sites_observed > 0,
      100 * rowMeans(dplyr::across(dplyr::all_of(cal_site_vars)) >= 6, na.rm = TRUE),
      NA_real_
    ),
    
    # Age groups for the figure
    age_group = cut(
      RIDAGEYR,
      breaks = c(30,35,40,45,50,55,60,65,70,75,80, Inf),
      right = FALSE,
      labels = c(
        "30–34","35–39","40–44","45–49","50–54","55–59",
        "60–64","65–69","70–74","75–79","80+"
      )
    )
  ) %>%
  dplyr::filter(!is.na(age_group)) %>%
  # Optional safety: exclude edentulous again (should already be gone upstream)
  dplyr::filter(n_sites_observed > 0)

# ---- Keep only what you need for the figure (plus design vars) ----
nhanes_fig_df <- nhanes_measures %>%
  dplyr::select(
    SEQN, source_cycle,
    RIDAGEYR, age_group, sex,
    SDMVPSU, SDMVSTRA, WTMEC6YR,
    n_sites_observed, pct_sites_observed,
    mean_CAL,
    perc_sites_ge3, perc_sites_ge5, perc_sites_ge6
  ) %>%
  dplyr::arrange(source_cycle, SEQN)

# ---- Save ----
saveRDS(nhanes_fig_df, "data/derived/nhanes_fig_df.rds")
readr::write_csv(nhanes_fig_df, "data/derived/nhanes_fig_df.csv")

message("Saved: data/derived/nhanes_fig_df.rds and .csv")
message("Rows: ", nrow(nhanes_fig_df))
message("Unique participants: ", dplyr::n_distinct(nhanes_fig_df$SEQN))