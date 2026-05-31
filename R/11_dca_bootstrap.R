# ============================================================================
# 11_dca_bootstrap.R
#
# Figure 4 (final): decision-curve analysis with 150-replicate bootstrap 95%
#   confidence ribbons around each net-benefit curve.
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
  library(boot)
})

adults <- readRDS(data_path("adults_cohort_v2.rds"))
cc_vars <- c("bmi_cat","age_w1","diabetes_treat_f","pa_meets",
             "smoking_3","high_chol_3","sleep_3","gpmp_tca",
             "time_years","event")
cc <- adults %>% select(all_of(cc_vars)) %>% na.omit()

HORIZON <- 5
thr <- seq(0.02, 0.30, by = 0.01)

# Fast predicted 5y risk
predict_5yr_fast <- function(model, newdata) {
  lp <- predict(model, newdata = newdata, type = "lp", reference = "sample")
  basehaz_df <- basehaz(model, centered = TRUE)
  idx <- max(which(basehaz_df$time <= HORIZON))
  H0 <- basehaz_df$hazard[idx]
  S0 <- exp(-H0)
  surv <- S0^exp(lp)
  1 - as.numeric(surv)
}

dca_compute <- function(time, event, risk, thresholds, horizon = 5) {
  out <- data.frame(threshold = thresholds, nb_model = NA, nb_all = NA)
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
      km_t <- survfit(Surv(time[treated], event[treated]) ~ 1)
      s_t <- tryCatch(summary(km_t, times = horizon, extend = TRUE)$surv,
                      error = function(e) NA)
      p_t_obs <- 1 - s_t
      n <- length(time)
      tp_n <- (n_treat / n) * p_t_obs
      fp_n <- (n_treat / n) * (1 - p_t_obs)
      out$nb_model[i] <- tp_n - fp_n * (p_t / (1 - p_t))
    }
    out$nb_all[i] <- p_all - (1 - p_all) * (p_t / (1 - p_t))
  }
  out
}

dca_boot <- function(data, indices, thresholds, horizon) {
  d <- data[indices, ]
  m_full <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f +
                    pa_meets + smoking_3 + high_chol_3 + sleep_3 + gpmp_tca, data = d)
  m_clin <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f, data = d)
  m_bmi  <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1, data = d)
  p_full <- predict_5yr_fast(m_full, d)
  p_clin <- predict_5yr_fast(m_clin, d)
  p_bmi  <- predict_5yr_fast(m_bmi, d)
  dca_full <- dca_compute(d$time_years, d$event, p_full, thresholds, horizon)
  dca_clin <- dca_compute(d$time_years, d$event, p_clin, thresholds, horizon)
  dca_bmi  <- dca_compute(d$time_years, d$event, p_bmi,  thresholds, horizon)
  c(dca_full$nb_model, dca_clin$nb_model, dca_bmi$nb_model, dca_full$nb_all)
}

# Quick timing
cat("Timing single iteration...\n")
t0 <- Sys.time()
test <- dca_boot(cc, 1:nrow(cc), thr, HORIZON)
cat(sprintf("Single iteration: %.2f sec\n", as.numeric(difftime(Sys.time(),t0,units="secs"))))

B <- 150
cat(sprintf("\nRunning DCA bootstrap (B=%d)...\n", B))
flush.console()
t0 <- Sys.time()
boot_dca <- boot(data = cc, statistic = dca_boot, R = B,
                 thresholds = thr, horizon = HORIZON)
cat(sprintf("Bootstrap done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

nthr <- length(thr)
get_ci <- function(idx) {
  ci <- boot.ci(boot_dca, index = idx, type = "perc", conf = 0.95)
  c(ci$percent[4], ci$percent[5])
}

extract_curve <- function(start_idx) {
  pts <- boot_dca$t0[(start_idx-1)*nthr + 1:nthr]
  los <- his <- numeric(nthr)
  for (i in 1:nthr) {
    ci <- get_ci((start_idx-1)*nthr + i)
    los[i] <- ci[1]; his[i] <- ci[2]
  }
  data.frame(threshold = thr, nb = pts, lo = los, hi = his)
}

dca_full_ci <- extract_curve(1) %>% mutate(model = "Full model")
dca_clin_ci <- extract_curve(2) %>% mutate(model = "Clinical (BMI+age+diabetes)")
dca_bmi_ci  <- extract_curve(3) %>% mutate(model = "BMI + age")
dca_all_ci  <- extract_curve(4) %>% mutate(model = "Treat all")
dca_none_ci <- data.frame(threshold = thr, nb = 0, lo = 0, hi = 0, model = "Treat none")

dca_combined <- bind_rows(dca_full_ci, dca_clin_ci, dca_bmi_ci, dca_all_ci, dca_none_ci) %>%
  mutate(model = factor(model, levels = c("Treat all","Treat none","BMI + age",
                                           "Clinical (BMI+age+diabetes)","Full model")))

cat("\nDCA at clinically meaningful thresholds (with 95% bootstrap CI):\n")
selected <- dca_combined %>% filter(threshold %in% c(0.05, 0.10, 0.15, 0.20),
                                     !model %in% c("Treat all","Treat none")) %>%
  arrange(threshold, model) %>%
  mutate(text = sprintf("%.4f (%.4f, %.4f)", nb, lo, hi))
print(selected[, c("threshold","model","text")], row.names = FALSE)

p_dca_ci <- ggplot(dca_combined %>% filter(!model %in% c("Treat all","Treat none")),
                   aes(x = threshold, y = nb, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_line(data = dca_combined %>% filter(model == "Treat all"),
            aes(x = threshold, y = nb), color = "#9CA3AF", linetype = "dotted",
            linewidth = 0.7, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("BMI + age" = "#3B82F6",
                                "Clinical (BMI+age+diabetes)" = "#F59E0B",
                                "Full model" = "#DC2626")) +
  scale_fill_manual(values = c("BMI + age" = "#3B82F6",
                               "Clinical (BMI+age+diabetes)" = "#F59E0B",
                               "Full model" = "#DC2626")) +
  scale_y_continuous(limits = c(-0.05, 0.10), breaks = seq(-0.05, 0.10, 0.025)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0.05, 0.30, 0.05)) +
  labs(x = "Threshold probability",
       y = "Net benefit (95% bootstrap CI shaded)",
       color = NULL, fill = NULL,
       caption = "Treat-all curve shown for reference (dotted). Treat-none = 0 line.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        plot.caption = element_text(size = 9, color = "gray40"))

ggsave(fig_path("fig4_dca_with_ci.png"), p_dca_ci, width = 7.5, height = 5, dpi = 200, bg = "white")
cat("\nDCA with CI plot saved.\n")

saveRDS(list(dca_combined = dca_combined, boot_dca = boot_dca, thr = thr),
        data_path("dca_bootstrap.rds"))

