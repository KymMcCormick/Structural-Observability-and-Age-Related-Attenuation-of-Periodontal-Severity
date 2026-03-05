# Structural-Observability-and-Age-Related-Attenuation-of-Periodontal-Severity
Simulation code and analysis scripts for "Structural Observability and Age-Related Attenuation of Periodontal Severity" (McCormick, Guzzo, Amarasena, Jamieson). Contains a Monte Carlo simulation framework demonstrating how severity-dependent tooth loss distorts cross-sectional age–severity trajectories, and survey-weighted analysis code. 

Contents

00_packages.R - Package loading for the reproducibility pipeline
01_clean_merge.R - Build NHANES 2009–2014 analytic dataset: DEMO (age, sex, survey design vars, MEC weights); OHXPER (site-level CAL = attachment loss variables). Excludes edentulous participants. Output: nhanes_perio_site (one row per person, wide site-level CAL columns)
02_demographics.R - Survey-weighted descriptive analysis of NHANES 2009-2014 periodontal data.
03_construct_measures.R - Construct person-level periodontal measures FROM site-level CAL
04__fit model.R - Fits survey-weighted quadratic regression models for each periodontal severity outcome and tests for non-monotonicity in the CAL >=6mm trajectory using the delta method and a pairwise adjacent-group contrast.
05_figures.R - Computes survey-weighted age-specific means and 95% confidence intervals for each periodontal severity outcome, combines with quadratic model predictions from the modelling script, and produces Figure 2.
06_Simulation with critical test.R - Monte Carlo simulation of latent periodontal disease (clinical attachment loss, CAL) accumulating monotonically with age across 28 anatomical subunits (teeth), under three observability mechanisms (sensitivity testing).

Data
NHANES 2009–2014 data are publicly available at https://www.cdc.gov/nchs/nhanes. No raw data are included in this repository.
Requirements
R packages: dplyr, tidyr, ggplot2, patchwork, survey
Contact
kym.mccormick@adelaide.edu.au
