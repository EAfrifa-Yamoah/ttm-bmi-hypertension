# ============================================================================
# 04_cindex_bootstrap.R
#
# Harrell's C-index for nested models with 500-replicate bootstrap 95% CIs,
#   including incremental discrimination (delta-C) between models.
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
  library(boot)
})

adults <- readRDS(data_path("adults_cohort_v2.rds"))

# Use complete cases for fair comparison across models
cc_vars <- c("bmi_cat","age_w1","diabetes_treat_f","pa_meets",
             "smoking_3","high_chol_3","sleep_3","gpmp_tca",
             "time_years","event")
cc <- adults %>% select(all_of(cc_vars)) %>% na.omit()
cat("Complete-case N:", nrow(cc), "Events:", sum(cc$event), "\n")

# ---- Build model formulas ----
f_bmi_only <- Surv(time_years, event) ~ bmi_cat
f_age_only <- Surv(time_years, event) ~ age_w1
f_bmi_age  <- Surv(time_years, event) ~ bmi_cat + age_w1
f_clinical <- Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f
f_full     <- Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f +
                pa_meets + smoking_3 + high_chol_3 + sleep_3 + gpmp_tca

# ---- C-index for a model on a data subset ----
c_idx <- function(formula, dat) {
  mod <- coxph(formula, data = dat)
  s <- summary(mod)
  as.numeric(s$concordance[1])
}

# ---- Boot function: returns c-index for each model + differences ----
boot_fun <- function(data, indices) {
  d <- data[indices, ]
  c_bmi      <- c_idx(f_bmi_only,  d)
  c_age      <- c_idx(f_age_only,  d)
  c_bmi_age  <- c_idx(f_bmi_age,   d)
  c_clinical <- c_idx(f_clinical,  d)
  c_full     <- c_idx(f_full,      d)
  c(
    c_bmi      = c_bmi,
    c_age      = c_age,
    c_bmi_age  = c_bmi_age,
    c_clinical = c_clinical,
    c_full     = c_full,
    d_full_vs_bmi      = c_full - c_bmi,
    d_full_vs_clinical = c_full - c_clinical,
    d_clinical_vs_bmi  = c_clinical - c_bmi,
    d_bmi_age_vs_age   = c_bmi_age - c_age   # incremental value of BMI over age
  )
}

cat("\nRunning bootstrap (B=500)...\n")
boot_res <- boot(data = cc, statistic = boot_fun, R = 500)

# ---- Summary ----
nm <- names(boot_res$t0)
cat("\n=== C-index estimates (bootstrap-percentile 95% CI) ===\n\n")
ci_tbl <- data.frame(metric = nm, est = NA, lo = NA, hi = NA)
for (i in seq_along(nm)) {
  est <- boot_res$t0[i]
  ci  <- boot.ci(boot_res, index = i, type = "perc", conf = 0.95)$percent[4:5]
  ci_tbl$est[i] <- est
  ci_tbl$lo[i]  <- ci[1]
  ci_tbl$hi[i]  <- ci[2]
}
ci_tbl$text <- sprintf("%.4f (%.4f, %.4f)", ci_tbl$est, ci_tbl$lo, ci_tbl$hi)
print(ci_tbl[, c("metric","text")], row.names = FALSE)

saveRDS(list(boot = boot_res, ci_tbl = ci_tbl), data_path("cindex_bootstrap.rds"))
cat("\nDone.\n")

