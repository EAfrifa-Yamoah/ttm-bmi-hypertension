# ============================================================================
# 10_schoenfeld_plot.R
#
# Publication Figure 6: scaled Schoenfeld residual plots with loess smoothers
#   and 95% bands for each BMI contrast.
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
  library(tidyr)
})

# Use individual dummy variables so cox.zph treats them separately
adults <- readRDS(data_path("adults_cohort_v2.rds"))
adults$bmi_Underweight <- as.numeric(adults$bmi_cat == "Underweight")
adults$bmi_Overweight  <- as.numeric(adults$bmi_cat == "Overweight")
adults$bmi_Obese       <- as.numeric(adults$bmi_cat == "Obese")

m3b <- coxph(Surv(time_years, event) ~ bmi_Underweight + bmi_Overweight + bmi_Obese +
               age_w1 + diabetes_treat_f + pa_meets + smoking_3 + high_chol_3 +
               sleep_3 + gpmp_tca, data = adults)
ph_zph <- cox.zph(m3b, transform = "identity")
cat("Schoenfeld test (with individual BMI dummies):\n")
print(ph_zph$table[1:3, ])

sch <- as.data.frame(ph_zph$y)
sch$time <- ph_zph$x

bmi_cols <- c("bmi_Underweight","bmi_Overweight","bmi_Obese")
sch_long <- sch %>%
  select(time, all_of(bmi_cols)) %>%
  pivot_longer(-time, names_to = "var", values_to = "resid") %>%
  mutate(comparison = factor(
    case_when(
      var == "bmi_Underweight" ~ "Underweight vs Normal",
      var == "bmi_Overweight"  ~ "Overweight vs Normal",
      var == "bmi_Obese"       ~ "Obese vs Normal"
    ),
    levels = c("Underweight vs Normal","Overweight vs Normal","Obese vs Normal")
  ))

pvals <- ph_zph$table[bmi_cols, "p"]
panel_labels <- sprintf("%s (p = %.3f)",
                        c("Underweight vs Normal","Overweight vs Normal","Obese vs Normal"),
                        pvals)
names(panel_labels) <- c("Underweight vs Normal","Overweight vs Normal","Obese vs Normal")

p_sch <- ggplot(sch_long, aes(x = time, y = resid)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(alpha = 0.18, size = 0.5, color = "#475569") +
  geom_smooth(method = "loess", se = TRUE, color = "#DC2626", fill = "#FCA5A5",
              span = 0.75) +
  facet_wrap(~ comparison, scales = "free_y",
             labeller = labeller(comparison = panel_labels)) +
  labs(x = "Time (years since baseline)",
       y = expression(beta(t))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(fig_path("fig6_schoenfeld_bmi.png"), p_sch, width = 9, height = 4.5, dpi = 200, bg = "white")
cat("\nPolished Schoenfeld plot saved.\n")

# Save the p-values for reporting
saveRDS(list(ph_zph = ph_zph, pvals_bmi = pvals), data_path("schoenfeld_bmi.rds"))

