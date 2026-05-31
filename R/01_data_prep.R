# ============================================================================
# 01_data_prep.R
#
# Build the analytic cohort: apply exclusions, derive the treated-hypertension
#   survival outcome, and construct adult / all-ages / adolescent cohorts.
#
# Part of: BMI category and incident treated hypertension in Australian men
#          (Ten to Men cohort, PBS/MBS linkage)
# Run order: source scripts in numeric order, or use 00_run_all.R
# Working directory must be the repository root.
# ============================================================================

source(file.path("R", "config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# Load raw data (latin-1 encoding to handle non-ASCII)
raw <- read.csv(RAW_CSV, stringsAsFactors = FALSE)

cat("Raw data shape: ", nrow(raw), "rows x", ncol(raw), "columns\n")
cat("Unique IDs:", length(unique(raw$zdcid0001md)), "\n")
cat("Waves available:\n")
print(table(raw$wave))

# ---- Baseline (Wave 1) characteristics ----
# Most covariates are only collected at Wave 1
w1 <- raw %>%
  filter(wave == 1) %>%
  mutate(
    # Outcome flags
    htn_event_ever = as.logical(antihypertensive_any),
    htn_first_date = as.Date(first_htn_rx_date),
    # Wave 1 BMI category
    bmi_cat = factor(BMI_category,
                     levels = c("Underweight","Normal","Overweight","Obese")),
    # Demographic / clinical baseline
    age_w1 = age,
    # PA, smoking, etc - use exact factor levels from data
    pa_meets = (physical_activity == "Met threshold"),
    smoking_current = (smoking == "Current smoker"),
    smoking_former  = (smoking == "Former smoker"),
    smoking_never   = (smoking == "Non-smoker"),
    high_chol = (high_cholesterol == "Yes"),
    sleep_diff = (sleep_issues == "Yes"),  # NA in data means missing; only Yes coded
    diabetes_treat_f = factor(diabetes_treat,
                              levels = c("no diabetes","non-insulin treated","insulin-treated")),
    diabetes_any = diabetes_treat %in% c("non-insulin treated","insulin-treated"),
    gpmp_tca = (gpmp_tca_any == "TRUE"),
    married_defacto = (marital_status == "Married/De facto"),
    never_married   = (marital_status == "Never married"),
    rurality_f      = factor(rurality, levels = c("Major Cities","Inner Regional","Outer Regional/Remote")),
    heart_disease   = (heart_condition == "Yes"),
    stroke_yn       = (stroke == "Yes"),
    self_highBP_w1  = (highBP == "Yes")
  ) %>%
  select(zdcid0001md, age_w1, bmi_cat, htn_event_ever, htn_first_date,
         pa_meets, smoking_current, smoking_former, smoking_never,
         high_chol, sleep_diff,
         diabetes_treat_f, diabetes_any,
         gpmp_tca, married_defacto, never_married,
         rurality_f, heart_disease, stroke_yn, self_highBP_w1,
         employment, household_income, highest_edu,
         antihypertensive_treatment_initiation_date, highbpmedpbs)

cat("\nWave 1 records:", nrow(w1), "\n")

# ---- Define hypertension-free at baseline ----
# Baseline = 1 Jan 2014. Anyone with antihypertensive dispensing before this is excluded.
BASELINE_DATE  <- as.Date("2014-01-01")
CENSOR_DATE    <- as.Date("2025-02-28")

w1 <- w1 %>%
  mutate(
    rx_before_baseline = !is.na(htn_first_date) & htn_first_date < BASELINE_DATE,
    incident_event     = !is.na(htn_first_date) & htn_first_date >= BASELINE_DATE,
    event_date         = if_else(incident_event, htn_first_date, CENSOR_DATE),
    time_years         = as.numeric(difftime(event_date, BASELINE_DATE, units = "days")) / 365.25,
    event              = as.integer(incident_event)
  )

cat("\nDispensings before baseline (excluded):",
    sum(w1$rx_before_baseline, na.rm = TRUE), "\n")
cat("Incident events after baseline:", sum(w1$incident_event, na.rm = TRUE), "\n")

# ---- Analytic cohort: HTN-free at baseline + valid BMI ----
analytic <- w1 %>%
  filter(!rx_before_baseline) %>%        # HTN-free at baseline
  filter(!is.na(bmi_cat)) %>%             # valid BMI
  filter(time_years > 0)                  # positive follow-up

cat("\n=== Analytic cohort ===\n")
cat("N:", nrow(analytic), "\n")
cat("Events:", sum(analytic$event), "\n")
cat("Total person-years:", round(sum(analytic$time_years), 0), "\n")
cat("Overall IR per 1000 PY:",
    round(1000 * sum(analytic$event) / sum(analytic$time_years), 2), "\n")
cat("\nBaseline BMI distribution:\n")
print(table(analytic$bmi_cat))

# Set Normal as reference (clinically meaningful; addresses R2.4)
analytic$bmi_cat <- relevel(analytic$bmi_cat, ref = "Normal")

# Save
saveRDS(analytic, data_path("analytic_cohort.rds"))

# ---- Adults-only cohort (R1.2, R2.5: pre-specified primary analysis restriction) ----
adults <- analytic %>% filter(age_w1 >= 18)
cat("\n=== Adults-only (≥18) cohort ===\n")
cat("N:", nrow(adults), "\n")
cat("Events:", sum(adults$event), "\n")
cat("Total person-years:", round(sum(adults$time_years), 0), "\n")
cat("\nBMI by adults:\n")
print(table(adults$bmi_cat))

saveRDS(adults, data_path("adults_cohort.rds"))

# ---- Older adolescents 10-17 ----
youth <- analytic %>% filter(age_w1 < 18)
cat("\n=== Adolescents (<18) cohort ===\n")
cat("N:", nrow(youth), "\n")
cat("Events:", sum(youth$event), "\n")
saveRDS(youth, data_path("youth_cohort.rds"))

cat("\n=== Data prep complete ===\n")

