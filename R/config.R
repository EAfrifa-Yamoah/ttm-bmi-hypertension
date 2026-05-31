# ============================================================================
# config.R
# Central configuration: paths, constants, and shared options.
# Sourced at the top of each analysis script.
# ============================================================================

# ---- Directory layout (relative to repository root) -------------------------
# All scripts assume the working directory is the repository root.
# Set this in your R session with setwd("/path/to/ttm-bmi-hypertension")
# or open the project via the .Rproj file (if using RStudio).

DATA_DIR    <- "data"      # raw + derived data (.rds intermediates land here)
FIG_DIR     <- "figures"   # output figures (.png)
R_DIR       <- "R"         # analysis scripts

# Raw input file (UTF-8 converted; see data/README.md for the conversion step).
# This file is NOT distributed in the repository (see data/README.md and the
# data-availability statement). Place your copy at the path below.
RAW_CSV     <- file.path(DATA_DIR, "ttmcombined_utf8.csv")

# ---- Path helpers -----------------------------------------------------------
# Use these so every script reads/writes intermediates in one place.
data_path <- function(name) file.path(DATA_DIR, name)
fig_path  <- function(name) file.path(FIG_DIR,  name)

# ---- Analysis constants -----------------------------------------------------
BASELINE_DATE <- as.Date("2014-01-01")  # start of follow-up
CENSOR_DATE   <- as.Date("2025-02-28")  # administrative censoring date
HORIZON       <- 5                       # years, for risk prediction / DCA / calibration
SEED          <- 2026                    # global RNG seed for reproducibility

# Antihypertensive ATC classes used to define the treated-hypertension outcome
ANTIHTN_ATC   <- c("C02", "C03", "C07", "C08", "C09")

# Bootstrap replicates
B_CINDEX      <- 500   # C-index bootstrap
B_DCA         <- 150   # decision-curve net-benefit bootstrap
B_CALIB       <- 200   # calibration bootstrap (rms::calibrate)

# BMI categories (Normal is the reference)
BMI_LEVELS    <- c("Normal", "Underweight", "Overweight", "Obese")

# ---- Convenience ------------------------------------------------------------
# Create output dirs if missing
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, showWarnings = FALSE)
if (!dir.exists(FIG_DIR))  dir.create(FIG_DIR,  showWarnings = FALSE)

set.seed(SEED)
