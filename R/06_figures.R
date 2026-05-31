# ============================================================================
# 06_figures.R
#
# Generate manuscript Figures 1-5 (KM curves, forest plot, C-index, DCA,
#   calibration). Figures 4/5 are superseded by 11/13 in the final manuscript.
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
  library(survminer)
})

adults    <- readRDS(data_path("adults_cohort_v2.rds"))
cox_res   <- readRDS(data_path("cox_models_adults.rds"))
boot_res  <- readRDS(data_path("cindex_bootstrap.rds"))
dca_res   <- readRDS(data_path("dca_results.rds"))
calib_res <- readRDS(data_path("calibration.rds"))

# Theme
theme_pub <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 12),
          strip.text = element_text(face = "bold"))
}

# ============================================================================
# Figure 1: Kaplan-Meier survival by BMI category (Adults >=18)
# ============================================================================

km_fit <- survfit(Surv(time_years, event) ~ bmi_cat, data = adults)

p_km <- ggsurvplot(km_fit,
                   data = adults,
                   conf.int = TRUE,
                   risk.table = TRUE,
                   palette = c("#3B82F6","#9CA3AF","#F59E0B","#DC2626"),
                   legend.title = "Baseline BMI",
                   legend.labs = c("Normal","Underweight","Overweight","Obese"),
                   xlab = "Years since baseline",
                   ylab = "Hypertension-free probability",
                   xlim = c(0, 11),
                   break.x.by = 2,
                   ylim = c(0.65, 1.0),
                   ggtheme = theme_pub(),
                   risk.table.height = 0.28,
                   tables.theme = theme_cleantable())

ggsave(fig_path("fig1_km_adults.png"),
       plot = p_km$plot,
       width = 7, height = 5, dpi = 200, bg = "white")
# Combined version with risk table
png(fig_path("fig1_km_adults_full.png"), width = 7, height = 7, units = "in", res = 200, bg = "white")
print(p_km)
dev.off()

cat("Figure 1 (KM curves) saved.\n")

# ============================================================================
# Figure 2: Forest plot - HRs from M1, M2, M3
# ============================================================================

forest_df <- bind_rows(
  cox_res$results$r1 %>% mutate(Model = "M1: Unadjusted"),
  cox_res$results$r2 %>% mutate(Model = "M2: + age + diabetes"),
  cox_res$results$r3 %>% mutate(Model = "M3: Fully adjusted")
) %>%
  filter(grepl("bmi_cat", var)) %>%
  mutate(BMI = factor(gsub("bmi_cat","", var),
                      levels = c("Underweight","Overweight","Obese"))) %>%
  filter(!is.na(BMI))

p_forest <- ggplot(forest_df, aes(x = HR, y = BMI, color = Model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray60") +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.5),
                  size = 0.4) +
  scale_x_continuous(trans = "log", breaks = c(0.7, 1, 1.5, 2, 3),
                     limits = c(0.7, 3.5)) +
  scale_color_manual(values = c("#9CA3AF","#3B82F6","#1E3A8A")) +
  labs(x = "Hazard ratio (95% CI, log scale)\nReference: Normal BMI",
       y = NULL) +
  theme_pub() +
  theme(legend.position = "bottom")

ggsave(fig_path("fig2_forest.png"), p_forest, width = 7, height = 4, dpi = 200, bg = "white")
cat("Figure 2 (forest plot) saved.\n")

# ============================================================================
# Figure 3: C-index comparison with bootstrap CIs
# ============================================================================

ci_tbl <- boot_res$ci_tbl
ci_models <- ci_tbl[1:5, ]
ci_models$label <- c("BMI alone","Age alone","BMI + age",
                     "BMI + age + diabetes","Full model")
ci_models$label <- factor(ci_models$label, levels = ci_models$label)

p_cindex <- ggplot(ci_models, aes(y = label, x = est)) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.4, color = "#1E3A8A") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray60") +
  scale_x_continuous(limits = c(0.5, 0.8), breaks = seq(0.5, 0.8, 0.05)) +
  labs(x = "Harrell's C-index (bootstrap 95% CI, B=500)",
       y = NULL) +
  theme_pub()

ggsave(fig_path("fig3_cindex.png"), p_cindex, width = 7, height = 3.5, dpi = 200, bg = "white")
cat("Figure 3 (C-index) saved.\n")

# ============================================================================
# Figure 4: Decision Curve Analysis
# ============================================================================

dca_plot_df <- dca_res$dca_long %>%
  bind_rows(data.frame(model = "Treat all",  threshold = dca_res$thr,
                       nb_model = dca_res$dca_full$nb_all,
                       nb_all = NA, nb_none = 0)) %>%
  bind_rows(data.frame(model = "Treat none", threshold = dca_res$thr,
                       nb_model = 0, nb_all = NA, nb_none = 0))

dca_plot_df$model <- factor(dca_plot_df$model,
                            levels = c("Treat all","Treat none","BMI + age",
                                       "BMI + age + diab.","Full model"))

p_dca <- ggplot(dca_plot_df, aes(x = threshold, y = nb_model, color = model, linetype = model)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c("Treat all" = "#9CA3AF",
                                "Treat none" = "#000000",
                                "BMI + age" = "#3B82F6",
                                "BMI + age + diab." = "#F59E0B",
                                "Full model" = "#DC2626")) +
  scale_linetype_manual(values = c("Treat all" = "dotted",
                                   "Treat none" = "dashed",
                                   "BMI + age" = "solid",
                                   "BMI + age + diab." = "solid",
                                   "Full model" = "solid")) +
  scale_y_continuous(limits = c(-0.05, 0.10),
                     breaks = seq(-0.05, 0.10, 0.025)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0.05, 0.4, 0.05)) +
  labs(x = "Threshold probability",
       y = "Net benefit",
       color = NULL, linetype = NULL) +
  theme_pub() +
  theme(legend.position = "bottom")

ggsave(fig_path("fig4_dca.png"), p_dca, width = 7, height = 4.5, dpi = 200, bg = "white")
cat("Figure 4 (DCA) saved.\n")

# ============================================================================
# Figure 5: Calibration plot
# ============================================================================

# rms::calibrate returns a matrix-like object; extract observed-vs-predicted points
cal_bmi <- calib_res$cal_bmi_age
cal_ful <- calib_res$cal_full

# Each calibration object has columns: mean.predicted, KM, KM.corrected, lower, upper
df_cal_bmi <- as.data.frame(unclass(cal_bmi))
df_cal_ful <- as.data.frame(unclass(cal_ful))

cat("Calibration object columns (bmi+age):\n"); print(names(df_cal_bmi))
cat("Calibration object columns (full):\n"); print(names(df_cal_ful))

# Reshape: 1 - KM = predicted risk (since calibrate works on survival)
# rms::calibrate columns: mean.predicted (predicted SURVIVAL), KM (observed survival)
# convert to risk
df_cal_bmi <- df_cal_bmi %>%
  mutate(pred_risk = 1 - mean.predicted,
         obs_risk  = 1 - KM,
         obs_lo    = pmax(0, 1 - (KM + 1.96 * std.err)),
         obs_hi    = pmin(1, 1 - (KM - 1.96 * std.err)),
         model = "BMI + age")

df_cal_ful <- df_cal_ful %>%
  mutate(pred_risk = 1 - mean.predicted,
         obs_risk  = 1 - KM,
         obs_lo    = pmax(0, 1 - (KM + 1.96 * std.err)),
         obs_hi    = pmin(1, 1 - (KM - 1.96 * std.err)),
         model = "Full model")

cal_df <- bind_rows(df_cal_bmi, df_cal_ful)

p_cal <- ggplot(cal_df, aes(x = pred_risk, y = obs_risk, color = model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60") +
  geom_errorbar(aes(ymin = obs_lo, ymax = obs_hi), width = 0.005) +
  geom_point(size = 2) +
  geom_line() +
  scale_color_manual(values = c("BMI + age" = "#3B82F6",
                                "Full model" = "#DC2626")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.35)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.35)) +
  labs(x = "Predicted 5-year hypertension risk",
       y = "Observed 5-year hypertension risk (KM)",
       color = NULL) +
  theme_pub() +
  theme(legend.position = "bottom")

ggsave(fig_path("fig5_calibration.png"), p_cal, width = 7, height = 5, dpi = 200, bg = "white")
cat("Figure 5 (calibration) saved.\n")

cat("\nAll figures complete.\n")

