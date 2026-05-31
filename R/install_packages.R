# ============================================================================
# install_packages.R
# One-shot installer for all packages used by the analysis pipeline.
# Run once before the first pipeline run:  Rscript R/install_packages.R
# ============================================================================

required <- c(
  "dplyr", "tidyr", "readr",   # data manipulation
  "survival",                   # Cox models, Surv(), cox.zph()
  "rms",                        # calibration, cph()
  "boot",                       # bootstrap CIs
  "ggplot2", "survminer", "cowplot", "scales",  # plotting
  "mice",                       # multiple imputation
  "cmprsk"                      # competing risks (Fine–Gray)
)

missing <- required[!(required %in% installed.packages()[, "Package"])]

if (length(missing) == 0) {
  message("All required packages are already installed.")
} else {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

# Report versions for the record
message("\nInstalled versions:")
for (p in required) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
  message(sprintf("  %-12s %s", p, v))
}
