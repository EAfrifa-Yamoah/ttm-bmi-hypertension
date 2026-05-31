# ============================================================================
# 00_run_all.R
#
# Master script: runs the complete analysis pipeline in dependency order.
# Working directory must be the repository root.
#
# Usage:
#   Rscript R/00_run_all.R
# or, from an interactive R session at the repository root:
#   source("R/00_run_all.R")
#
# Expects the raw input file at the path defined by RAW_CSV in R/config.R
# (default: data/ttmcombined_utf8.csv). See data/README.md.
# ============================================================================

source(file.path("R", "config.R"))

# Scripts in strict dependency order. Each writes intermediates to data/
# and/or figures to figures/, consumed by later scripts.
scripts <- c(
  "01_data_prep.R",                 # -> analytic/adults/youth cohorts
  "02_table1.R",                    # -> Table 1
  "03_cox_primary.R",               # -> Cox models, adults_cohort_v2
  "04_cindex_bootstrap.R",          # -> C-index bootstrap CIs
  "05_calibration_dca.R",           # -> calibration + DCA (point estimates)
  "06_figures.R",                   # -> Figures 1-3 (4,5 superseded below)
  "07_sensitivity.R",               # -> sensitivity analyses
  "08_multiple_imputation.R",       # -> MI pooled estimates
  "09_interaction_evalue_cuminc.R", # -> interaction, E-values, cumulative incidence
  "10_schoenfeld_plot.R",           # -> Figure 6
  "11_dca_bootstrap.R",             # -> Figure 4 (final, with bootstrap CIs)
  "12_competing_risks.R",           # -> competing-risks sensitivity
  "13_calibration_deciles.R"        # -> Figure 5 (final, decile calibration)
)

# Check the raw data file exists before starting.
if (!file.exists(RAW_CSV)) {
  stop(sprintf(
    "Raw data not found at '%s'.\nPlace the UTF-8 converted Ten to Men / PBS linkage file there first.\nSee data/README.md for details.",
    RAW_CSV))
}

t_start <- Sys.time()
for (s in scripts) {
  message("\n", strrep("=", 70))
  message("Running: ", s)
  message(strrep("=", 70))
  t0 <- Sys.time()
  source(file.path(R_DIR, s), echo = FALSE)
  message(sprintf(">> %s completed in %.1f sec", s,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

message("\n", strrep("=", 70))
message(sprintf("Pipeline complete in %.1f min. Figures in '%s/', intermediates in '%s/'.",
                as.numeric(difftime(Sys.time(), t_start, units = "mins")),
                FIG_DIR, DATA_DIR))
message(strrep("=", 70))
