# ============================================================================
# 13_calibration_deciles.R
#
# Figure 5 (final): decile-based calibration plot with Greenwood 95% CIs.
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
cc_vars <- c("bmi_cat","age_w1","diabetes_treat_f","pa_meets",
             "smoking_3","high_chol_3","sleep_3","gpmp_tca",
             "time_years","event")
cc <- adults %>% select(all_of(cc_vars)) %>% na.omit()

HORIZON <- 5

predict_5yr_fast <- function(model, newdata) {
  lp <- predict(model, newdata = newdata, type = "lp", reference = "sample")
  basehaz_df <- basehaz(model, centered = TRUE)
  idx <- max(which(basehaz_df$time <= HORIZON))
  H0 <- basehaz_df$hazard[idx]
  S0 <- exp(-H0)
  surv <- S0^exp(lp)
  1 - as.numeric(surv)
}

m_bmi_age <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1, data = cc)
m_full <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f +
                  pa_meets + smoking_3 + high_chol_3 + sleep_3 + gpmp_tca, data = cc)

cc$risk_bmi <- predict_5yr_fast(m_bmi_age, cc)
cc$risk_full <- predict_5yr_fast(m_full, cc)

calibrate_deciles <- function(time, event, risk, horizon = 5, n_bins = 10) {
  bins <- cut(risk, quantile(risk, probs = seq(0, 1, length.out = n_bins+1)),
              include.lowest = TRUE, labels = FALSE)
  out <- data.frame(bin = 1:n_bins, predicted = NA, observed = NA, lo = NA, hi = NA, n = NA)
  for (b in 1:n_bins) {
    idx <- which(bins == b)
    out$n[b] <- length(idx)
    out$predicted[b] <- mean(risk[idx])
    km <- survfit(Surv(time[idx], event[idx]) ~ 1)
    ks <- summary(km, times = horizon, extend = TRUE)
    out$observed[b] <- 1 - ks$surv
    # 95% CI on observed risk (Greenwood)
    se <- ks$std.err
    out$lo[b] <- pmax(0, 1 - (ks$surv + 1.96*se))
    out$hi[b] <- pmin(1, 1 - (ks$surv - 1.96*se))
  }
  out
}

cal_bmi  <- calibrate_deciles(cc$time_years, cc$event, cc$risk_bmi, HORIZON)
cal_full <- calibrate_deciles(cc$time_years, cc$event, cc$risk_full, HORIZON)
cal_bmi$model  <- "BMI + age"
cal_full$model <- "Full model"
cal <- bind_rows(cal_bmi, cal_full)

cat("=== Calibration deciles ===\n")
print(cal %>% mutate(across(where(is.numeric), ~round(., 4))))

p_cal <- ggplot(cal, aes(x = predicted, y = observed, color = model, fill = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.005, alpha = 0.6) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = c("BMI + age" = "#3B82F6",
                                "Full model" = "#DC2626")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.32)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.32)) +
  labs(x = "Predicted 5-year hypertension risk (decile mean)",
       y = "Observed 5-year hypertension risk (KM, 95% CI)",
       color = NULL, fill = NULL,
       caption = "Each point represents one decile of predicted risk.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.caption = element_text(size = 9, color = "gray40"))

ggsave(fig_path("fig5_calibration_deciles.png"), p_cal, width = 6.5, height = 5, dpi = 200, bg = "white")
cat("\nDecile calibration plot saved.\n")

saveRDS(cal, data_path("calibration_deciles.rds"))

