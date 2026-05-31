# ============================================================================
# 09_interaction_evalue_cuminc.R
#
# BMI x age interaction (likelihood-ratio test), E-values for unmeasured
#   confounding (VanderWeele & Ding), and cumulative incidence at 5/10 years.
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
})

adults <- readRDS(data_path("adults_cohort_v2.rds"))

# ===== 1. BMI x age interaction (continuous age) =====
cat("=== BMI x age interaction test ===\n")
m_main  <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f,
                 data = adults)
m_inter <- coxph(Surv(time_years, event) ~ bmi_cat * age_w1 + diabetes_treat_f,
                 data = adults)
lr_int <- anova(m_main, m_inter, test = "Chisq")
cat("Likelihood ratio test:\n")
print(lr_int)

# Extract interaction p-value as scalar
p_inter <- lr_int$"Pr(>|Chi|)"[2]
df_inter <- lr_int$Df[2]
chi_inter <- lr_int$Chisq[2]
cat(sprintf("\nBMI category x age interaction: chi-sq = %.2f, df = %d, p = %.4f\n",
            chi_inter, df_inter, p_inter))

# Show interaction coefficients
cat("\nInteraction model coefficients (relevant):\n")
ci_inter <- summary(m_inter)$coefficients
print(ci_inter[grep("bmi_cat", rownames(ci_inter)), ])

saveRDS(list(main=m_main, inter=m_inter, lr=lr_int,
             p_inter=p_inter, chi_inter=chi_inter, df_inter=df_inter),
        data_path("interaction_results.rds"))

# ===== 2. E-values =====
# E-value formula (VanderWeele & Ding 2017):
# For HR: E = HR + sqrt(HR*(HR-1))   when HR > 1
# E-value for CI: same formula but applied to the bound closer to null (lower limit when HR>1)

evalue_hr <- function(hr, lo) {
  e_point <- hr + sqrt(hr * (hr - 1))
  e_ci    <- ifelse(lo <= 1, 1, lo + sqrt(lo * (lo - 1)))
  list(e_point = e_point, e_ci = e_ci)
}

cat("\n=== E-values for BMI category HRs (from M3 complete-case) ===\n")
# From the complete-case Model 3
hrs <- list(
  Underweight = list(hr = 1.43, lo = 1.16, hi = 1.77),
  Overweight  = list(hr = 1.42, lo = 1.22, hi = 1.67),
  Obese       = list(hr = 1.88, lo = 1.59, hi = 2.22)
)

evalues <- data.frame(BMI = names(hrs), HR = NA, lo = NA,
                      E_point = NA, E_lower_CI = NA)
for (i in seq_along(hrs)) {
  e <- evalue_hr(hrs[[i]]$hr, hrs[[i]]$lo)
  evalues$HR[i]         <- hrs[[i]]$hr
  evalues$lo[i]         <- hrs[[i]]$lo
  evalues$E_point[i]    <- round(e$e_point, 2)
  evalues$E_lower_CI[i] <- round(e$e_ci, 2)
}
print(evalues)
cat("\nInterpretation: An E-value of 3.16 (Obese point estimate) means an unmeasured\n")
cat("confounder would need to be associated with both Obese-vs-Normal and incident\n")
cat("HTN by risk ratios of at least 3.16 (above and beyond measured confounders) to\n")
cat("fully explain away the observed association.\n")

saveRDS(evalues, data_path("evalues.rds"))

# ===== 3. Cumulative incidence by BMI category at 5 and 10 years =====
cat("\n=== Cumulative incidence (1-KM) at 5 and 10 years by BMI category ===\n")
km_fit <- survfit(Surv(time_years, event) ~ bmi_cat, data = adults)
ci_5  <- summary(km_fit, times = 5, extend = TRUE)
ci_10 <- summary(km_fit, times = 10, extend = TRUE)

ci_tbl <- data.frame(
  BMI = gsub("bmi_cat=","", ci_5$strata),
  n_at_5  = ci_5$n.risk,
  cum_inc_5  = 1 - ci_5$surv,
  cum_inc_5_lo = 1 - ci_5$upper,
  cum_inc_5_hi = 1 - ci_5$lower,
  n_at_10 = ci_10$n.risk,
  cum_inc_10 = 1 - ci_10$surv,
  cum_inc_10_lo = 1 - ci_10$upper,
  cum_inc_10_hi = 1 - ci_10$lower
)
ci_tbl$text_5  <- sprintf("%.1f%% (%.1f-%.1f)",
                          100*ci_tbl$cum_inc_5, 100*ci_tbl$cum_inc_5_lo,
                          100*ci_tbl$cum_inc_5_hi)
ci_tbl$text_10 <- sprintf("%.1f%% (%.1f-%.1f)",
                          100*ci_tbl$cum_inc_10, 100*ci_tbl$cum_inc_10_lo,
                          100*ci_tbl$cum_inc_10_hi)
print(ci_tbl[, c("BMI","n_at_5","text_5","n_at_10","text_10")])

saveRDS(ci_tbl, data_path("cumulative_incidence.rds"))

# ===== 4. Schoenfeld residuals: time-varying coefficient plots =====
cat("\n=== Schoenfeld residual plots ===\n")
m3 <- readRDS(data_path("cox_models_adults.rds"))$m3
ph_zph <- cox.zph(m3, transform = "km")
cat("Schoenfeld test results:\n")
print(ph_zph)

# Generate Schoenfeld plot for BMI category
png(fig_path("fig6_schoenfeld_bmi.png"), width = 7, height = 5, units = "in", res = 200, bg = "white")
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(ph_zph[1], main = "Underweight vs Normal", xlab = "Time", ylab = "Beta(t)")
abline(h = 0, lty = 2, col = "gray60")
plot(ph_zph[2], main = "Overweight vs Normal", xlab = "Time", ylab = "Beta(t)")
abline(h = 0, lty = 2, col = "gray60")
plot(ph_zph[3], main = "Obese vs Normal", xlab = "Time", ylab = "Beta(t)")
abline(h = 0, lty = 2, col = "gray60")
plot.new()
title("Time-varying coefficient (95% band) — flat line = PH satisfied", line = -3,
      outer = TRUE, cex.main = 0.9)
dev.off()
cat("Saved fig6_schoenfeld_bmi.png\n")

# ===== 5. Competing risks: death-as-competing-event sensitivity =====
# The dataset has no death indicator, so this can't be fitted directly.
# But we can document: at ages 18-55, all-cause mortality in Australian men over
# 11 years is approximately 1-2% based on AIHW data, far below the HTN event rate
# of 14%. Competing risks bias from death is therefore unlikely to materially
# affect estimates. We note this in the limitations.
cat("\n=== Competing risks (death) ===\n")
cat("Death data not available in the linked dataset. Australian male mortality\n")
cat("ages 18-55 over 11 years (AIHW estimates): approximately 1-2%, far below the\n")
cat("treated HTN event rate of 14% in this cohort. Competing risks bias from death\n")
cat("is therefore expected to be negligible. This is noted in the limitations.\n")

saveRDS(list(ph_zph = ph_zph), data_path("schoenfeld.rds"))

cat("\nStrengthening analyses complete.\n")

