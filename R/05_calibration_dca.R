# ============================================================================
# 05_calibration_dca.R
#
# Bootstrap-corrected calibration (rms) and a first-principles decision-curve
#   analysis (point estimates). See 11/13 for the bootstrap-CI and decile versions.
#
# Part of: BMI category and incident treated hypertension in Australian men
#          (Ten to Men cohort, PBS/MBS linkage)
# Run order: source scripts in numeric order, or use 00_run_all.R
# Working directory must be the repository root.
# ============================================================================

source(file.path("R", "config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(ggplot2)
  library(rms)
})

adults <- readRDS(data_path("adults_cohort_v2.rds"))

cc_vars <- c("bmi_cat","age_w1","diabetes_treat_f","pa_meets",
             "smoking_3","high_chol_3","sleep_3","gpmp_tca",
             "time_years","event")
cc <- adults %>% select(all_of(cc_vars)) %>% na.omit()

# Set up datadist for rms
dd <- datadist(cc); options(datadist = "dd")

# ---- Fit Cox models with rms::cph for calibration ----
# Time horizon: 5 years
HORIZON <- 5

mod_bmi_age <- cph(Surv(time_years, event) ~ bmi_cat + age_w1,
                   data = cc, x = TRUE, y = TRUE, surv = TRUE, time.inc = HORIZON)
mod_full <- cph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f +
                  pa_meets + smoking_3 + high_chol_3 + sleep_3 + gpmp_tca,
                data = cc, x = TRUE, y = TRUE, surv = TRUE, time.inc = HORIZON)

# ---- Calibration at 5 years (using rms::calibrate) ----
cat("\n--- Calibration: BMI + age model ---\n")
cal_bmi_age <- calibrate(mod_bmi_age, cmethod = "KM", method = "boot", u = HORIZON,
                         m = 200, B = 200)

cat("\n--- Calibration: Full model ---\n")
cal_full <- calibrate(mod_full, cmethod = "KM", method = "boot", u = HORIZON,
                      m = 200, B = 200)

# Save calibration outputs
saveRDS(list(cal_bmi_age = cal_bmi_age, cal_full = cal_full), data_path("calibration.rds"))

# ---- Decision Curve Analysis ----
# Computes net benefit across threshold probabilities

# Get predicted 5-year risk from a fitted Cox model
predict_risk_5yr <- function(model, newdata) {
  # survfit gives the survival curve for each subject
  sf <- survfit(model, newdata = newdata)
  # Find the survival at time = HORIZON
  st <- summary(sf, times = HORIZON, extend = TRUE)$surv
  1 - as.numeric(st)
}

# Refit coxph (not cph) for survfit
m_bmi_age <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1, data = cc)
m_clinical<- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f, data = cc)
m_full    <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f +
                     pa_meets + smoking_3 + high_chol_3 + sleep_3 + gpmp_tca,
                   data = cc)

cat("\nPredicting 5-year risks...\n")
p_bmi_age  <- predict_risk_5yr(m_bmi_age,  cc)
p_clinical <- predict_risk_5yr(m_clinical, cc)
p_full     <- predict_risk_5yr(m_full,     cc)

# Observed: KM-based survival probability — for DCA we use the original event indicator
# at horizon (administrative censoring approach)
# For DCA with survival, the simplest "intent" version uses observed events by horizon,
# treating censored-before-horizon as non-events (note this is biased toward null;
# more rigorous would use the Pencina/Vickers KM-based approach, but for this n it's reasonable)

# Use IPCW-free DCA via the standard approach: net benefit at threshold p_t
# NB(p_t) = TP/N - FP/N * (p_t / (1-p_t))

# But for survival, we need to account for censoring. Use the survival-adjusted form:
# NB(p_t) = (event-rate among treated * I[treated]) - (1 - event-rate among treated) * I[treated] * p_t/(1-p_t)
# The Vickers et al. method using KM among classified-treated is the standard.

# Implement survival-DCA manually
dca_survival <- function(time, event, risk, thresholds, horizon = 5) {
  out <- data.frame(threshold = thresholds,
                    nb_model = NA_real_,
                    nb_all   = NA_real_,
                    nb_none  = 0)
  # Overall KM-based event rate at horizon
  km_all <- survfit(Surv(time, event) ~ 1)
  s_all <- summary(km_all, times = horizon, extend = TRUE)$surv
  p_all <- 1 - s_all

  for (i in seq_along(thresholds)) {
    p_t <- thresholds[i]
    treated <- risk >= p_t
    n_treat <- sum(treated)
    if (n_treat == 0) {
      out$nb_model[i] <- 0
    } else {
      km_treat <- survfit(Surv(time[treated], event[treated]) ~ 1)
      s_t <- tryCatch(summary(km_treat, times = horizon, extend = TRUE)$surv,
                      error = function(e) NA)
      p_treat <- 1 - s_t
      n <- length(time)
      # Net benefit = (TP/N) - (FP/N) * (p_t/(1-p_t))
      # TP/N = (n_treat/N) * p_treat
      # FP/N = (n_treat/N) * (1 - p_treat)
      tp_n <- (n_treat / n) * p_treat
      fp_n <- (n_treat / n) * (1 - p_treat)
      out$nb_model[i] <- tp_n - fp_n * (p_t / (1 - p_t))
    }
    # Treat-all
    out$nb_all[i] <- p_all - (1 - p_all) * (p_t / (1 - p_t))
  }
  out
}

thr <- seq(0.02, 0.40, by = 0.01)

cat("Computing DCA curves...\n")
dca_bmi_age  <- dca_survival(cc$time_years, cc$event, p_bmi_age,  thr, HORIZON)
dca_clinical <- dca_survival(cc$time_years, cc$event, p_clinical, thr, HORIZON)
dca_full     <- dca_survival(cc$time_years, cc$event, p_full,     thr, HORIZON)

dca_long <- bind_rows(
  data.frame(model = "BMI + age",         dca_bmi_age),
  data.frame(model = "BMI + age + diab.", dca_clinical),
  data.frame(model = "Full model",        dca_full)
)

cat("\nDCA results (selected thresholds):\n")
print(dca_long %>% filter(threshold %in% c(0.05, 0.10, 0.15, 0.20, 0.30)) %>%
        mutate(across(where(is.numeric), ~round(., 4))))

saveRDS(list(thr = thr, dca_long = dca_long,
             dca_bmi_age = dca_bmi_age, dca_clinical = dca_clinical, dca_full = dca_full,
             pred_5yr = list(bmi_age = p_bmi_age, clinical = p_clinical, full = p_full)),
        data_path("dca_results.rds"))

cat("\nDCA complete.\n")

