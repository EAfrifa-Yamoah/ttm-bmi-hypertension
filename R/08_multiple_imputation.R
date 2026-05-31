# ============================================================================
# 08_multiple_imputation.R
#
# Multiple imputation by chained equations (mice, m=20) for missing covariates;
#   pooled Cox estimates via Rubin's rules, compared with complete-case.
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
  library(mice)
})

adults <- readRDS(data_path("adults_cohort_v2.rds"))

# ---- Identify missingness pattern ----
miss_vars <- c("bmi_cat","age_w1","diabetes_treat_f","pa_meets",
               "smoking_3","high_chol_3","sleep_3","gpmp_tca")
imp_data <- adults %>% select(all_of(miss_vars), time_years, event)

cat("=== Missingness pattern ===\n")
print(sapply(imp_data, function(x) sum(is.na(x))))
cat(sprintf("\nComplete cases: %d / %d (%.1f%%)\n",
            sum(complete.cases(imp_data)), nrow(imp_data),
            100*sum(complete.cases(imp_data))/nrow(imp_data)))

# The "Missing" level in high_chol_3 needs to become NA for imputation
imp_data$high_chol_3 <- ifelse(imp_data$high_chol_3 == "Missing", NA,
                                as.character(imp_data$high_chol_3))
imp_data$high_chol_3 <- factor(imp_data$high_chol_3, levels = c("No","Yes"))

# Create Nelson-Aalen estimate to include in imputation model (best practice for survival)
imp_data$na_est <- nelsonaalen(imp_data, time_years, event)

cat("\n=== Re-tally missingness ===\n")
print(sapply(imp_data, function(x) sum(is.na(x))))

# ---- Run MICE (m=20 imputations, 10 iterations) ----
cat("\nRunning MICE...\n")
predmat <- make.predictorMatrix(imp_data)
# Use time_years + event + na_est as predictors; don't impute them
predmat[, "time_years"] <- 0  # don't use time_years directly (na_est captures this)
predmat[, "event"] <- 1
predmat[, "na_est"] <- 1
# Don't impute time/event/na_est
methvec <- make.method(imp_data)
methvec["time_years"] <- ""
methvec["event"] <- ""
methvec["na_est"] <- ""

imp <- mice(imp_data, m = 20, maxit = 10, method = methvec,
            predictorMatrix = predmat, printFlag = FALSE, seed = 2026)

cat("Imputation complete.\n")

# ---- Fit Cox model on each imputation and pool ----
cat("\nFitting Cox models on imputed data...\n")
fit_imp <- with(imp, coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f +
                             pa_meets + smoking_3 + high_chol_3 + sleep_3 + gpmp_tca))

pooled <- pool(fit_imp)
pooled_summary <- summary(pooled, conf.int = TRUE, exponentiate = TRUE)
rownames(pooled_summary) <- as.character(pooled_summary$term)

cat("\n=== Pooled Cox results (MI, m=20) ===\n")
print(pooled_summary[, c("estimate","2.5 %","97.5 %","p.value")], digits = 3)

# Format for reporting
bmi_rows <- pooled_summary[grep("bmi_cat", pooled_summary$term), ]
cat("\n=== BMI HRs (MI vs complete-case) ===\n")
cox_complete <- readRDS(data_path("cox_models_adults.rds"))$results$r3
cat("Complete-case (n=9,309):\n")
print(cox_complete[grep("bmi_cat", cox_complete$var), c("var","HR","lo","hi")])
cat("\nMI-pooled (n=12,742, m=20):\n")
mi_out <- data.frame(
  var = bmi_rows$term,
  HR = round(bmi_rows$estimate, 3),
  lo = round(bmi_rows$`2.5 %`, 3),
  hi = round(bmi_rows$`97.5 %`, 3),
  text = sprintf("%.2f (%.2f-%.2f)", bmi_rows$estimate, bmi_rows$`2.5 %`, bmi_rows$`97.5 %`)
)
print(mi_out)

saveRDS(list(imp = imp, pooled = pooled, summary = pooled_summary, mi_out = mi_out),
        data_path("mi_results.rds"))

cat("\nMI analysis complete.\n")

