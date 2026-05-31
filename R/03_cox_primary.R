# ============================================================================
# 03_cox_primary.R
#
# Primary Cox proportional-hazards models (M1 unadjusted, M2 +age+diabetes,
#   M3 fully adjusted) and the proportional-hazards (Schoenfeld) test.
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
})

adults <- readRDS(data_path("adults_cohort.rds"))

# Make sure Normal is reference (R2.4)
adults$bmi_cat <- relevel(factor(adults$bmi_cat,
                                 levels = c("Normal","Underweight","Overweight","Obese")),
                          ref = "Normal")

# ---- Helper to format Cox output ----
format_cox <- function(model, label = "") {
  s <- summary(model)
  hr  <- s$conf.int[, "exp(coef)"]
  lo  <- s$conf.int[, "lower .95"]
  hi  <- s$conf.int[, "upper .95"]
  p   <- s$coefficients[, "Pr(>|z|)"]
  out <- data.frame(
    var = rownames(s$coefficients),
    HR  = round(hr, 3),
    lo  = round(lo, 3),
    hi  = round(hi, 3),
    p   = signif(p, 3),
    text = sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)
  )
  cat("\n--- ", label, " ---\n", sep="")
  print(out, row.names = FALSE)
  cat(sprintf("n events / n: %d / %d\n", s$nevent, s$n))
  cat(sprintf("Concordance: %.4f (SE %.4f)\n", s$concordance[1], s$concordance[2]))
  cat(sprintf("Log-likelihood: %.2f\n", s$loglik[2]))
  out
}

# ---- Model 1: Unadjusted ----
m1 <- coxph(Surv(time_years, event) ~ bmi_cat, data = adults)
r1 <- format_cox(m1, "Model 1: BMI category, unadjusted")

# ---- Model 2: + age + diabetes ----
m2 <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f, data = adults)
r2 <- format_cox(m2, "Model 2: + age + diabetes treatment")

# ---- Model 3: Fully adjusted ----
# pa_meets, smoking (with NAs - we'll keep as 3-level), high_chol, sleep, gpmp/HbA1c proxy
adults <- adults %>%
  mutate(
    smoking_3 = factor(case_when(
      smoking_current ~ "Current",
      smoking_former  ~ "Former",
      smoking_never   ~ "Never",
      TRUE ~ NA_character_
    ), levels = c("Never","Former","Current")),
    high_chol_3 = factor(case_when(
      high_chol ~ "Yes",
      !high_chol ~ "No",
      TRUE ~ "Missing"
    ), levels = c("No","Yes","Missing")),
    sleep_3 = factor(case_when(
      sleep_diff ~ "Yes",
      TRUE ~ "No/Missing"
    ), levels = c("No/Missing","Yes"))  # only Yes is recorded, others = absent
  )

m3 <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f +
              pa_meets + smoking_3 + high_chol_3 + sleep_3 + gpmp_tca,
            data = adults)
r3 <- format_cox(m3, "Model 3: Fully adjusted")

# ---- Proportional hazards check on M3 ----
ph_test <- cox.zph(m3)
cat("\n--- Proportional hazards test (M3) ---\n")
print(ph_test)

# ---- Save models ----
saveRDS(list(m1 = m1, m2 = m2, m3 = m3, results = list(r1=r1, r2=r2, r3=r3),
             ph_test = ph_test),
        data_path("cox_models_adults.rds"))
saveRDS(adults, data_path("adults_cohort_v2.rds"))  # with smoking_3 etc

cat("\nPrimary Cox models complete.\n")

