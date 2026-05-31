# ============================================================================
# 02_table1.R
#
# Descriptive Table 1: baseline characteristics by BMI category (adult cohort).
#
# Part of: BMI category and incident treated hypertension in Australian men
#          (Ten to Men cohort, PBS/MBS linkage)
# Run order: source scripts in numeric order, or use 00_run_all.R
# Working directory must be the repository root.
# ============================================================================

source(file.path("R", "config.R"))

suppressPackageStartupMessages({
  library(dplyr)
})

adults <- readRDS(data_path("adults_cohort.rds"))

n_pct_logical <- function(x) {
  n <- sum(x, na.rm = TRUE); denom <- sum(!is.na(x))
  if (denom == 0) return("--")
  sprintf("%d (%.1f)", n, 100 * n / denom)
}

n_pct_eq <- function(x, level) {
  n <- sum(x == level, na.rm = TRUE); denom <- sum(!is.na(x))
  if (denom == 0) return("--")
  sprintf("%d (%.1f)", n, 100 * n / denom)
}

mean_sd <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))

bmi_levels <- c("Normal","Underweight","Overweight","Obese")
adults$bmi_cat <- factor(adults$bmi_cat, levels = bmi_levels)

p_continuous <- function(var) {
  m <- aov(as.formula(paste(var, "~ bmi_cat")), data = adults)
  pv <- summary(m)[[1]][["Pr(>F)"]][1]
  if (pv < 0.001) "<0.001" else sprintf("%.3f", pv)
}

p_logical <- function(var) {
  x <- adults[[var]]; if (sum(!is.na(x)) == 0) return("--")
  tab <- table(adults$bmi_cat, x); if (any(rowSums(tab) == 0)) return("--")
  p <- tryCatch(chisq.test(tab)$p.value, error = function(e) NA)
  if (is.na(p)) "--" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
}

p_factor <- function(var) {
  x <- adults[[var]]; if (sum(!is.na(x)) == 0) return("--")
  tab <- table(adults$bmi_cat, x); if (any(rowSums(tab) == 0)) return("--")
  p <- tryCatch(chisq.test(tab)$p.value, error = function(e) NA)
  if (is.na(p)) "--" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
}

cat("\n=== Table 1: Adults (>=18) by BMI category ===\n\n")
cat(sprintf("Total N = %d\n\n", nrow(adults)))

header <- sprintf("%-30s | %12s %12s %12s %12s | %12s | %10s",
                  "Variable", "Normal", "Underwt", "Overweight", "Obese", "Total", "p-value")
cat(header, "\n"); cat(strrep("-", nchar(header)), "\n", sep="")

counts <- table(adults$bmi_cat)
cat(sprintf("%-30s | %12s %12s %12s %12s | %12s |\n",
  "N (%)",
  sprintf("%d (%.1f)", counts["Normal"], 100*counts["Normal"]/nrow(adults)),
  sprintf("%d (%.1f)", counts["Underweight"], 100*counts["Underweight"]/nrow(adults)),
  sprintf("%d (%.1f)", counts["Overweight"], 100*counts["Overweight"]/nrow(adults)),
  sprintf("%d (%.1f)", counts["Obese"], 100*counts["Obese"]/nrow(adults)),
  sprintf("%d", nrow(adults))))

cat(sprintf("%-30s | %12s %12s %12s %12s | %12s | %10s\n",
  "Age, mean (SD)",
  mean_sd(adults$age_w1[adults$bmi_cat=="Normal"]),
  mean_sd(adults$age_w1[adults$bmi_cat=="Underweight"]),
  mean_sd(adults$age_w1[adults$bmi_cat=="Overweight"]),
  mean_sd(adults$age_w1[adults$bmi_cat=="Obese"]),
  mean_sd(adults$age_w1),
  p_continuous("age_w1")))

print_log_row <- function(label, var) {
  cat(sprintf("%-30s | %12s %12s %12s %12s | %12s | %10s\n",
    label,
    n_pct_logical(adults[[var]][adults$bmi_cat=="Normal"]),
    n_pct_logical(adults[[var]][adults$bmi_cat=="Underweight"]),
    n_pct_logical(adults[[var]][adults$bmi_cat=="Overweight"]),
    n_pct_logical(adults[[var]][adults$bmi_cat=="Obese"]),
    n_pct_logical(adults[[var]]),
    p_logical(var)))
}
print_eq_row <- function(label, var, level) {
  cat(sprintf("%-30s | %12s %12s %12s %12s | %12s | %10s\n",
    label,
    n_pct_eq(adults[[var]][adults$bmi_cat=="Normal"], level),
    n_pct_eq(adults[[var]][adults$bmi_cat=="Underweight"], level),
    n_pct_eq(adults[[var]][adults$bmi_cat=="Overweight"], level),
    n_pct_eq(adults[[var]][adults$bmi_cat=="Obese"], level),
    n_pct_eq(adults[[var]], level),
    p_factor(var)))
}

print_log_row("Married/de facto", "married_defacto")
print_log_row("Never married", "never_married")
print_log_row("PA met guidelines", "pa_meets")
print_log_row("Current smoker", "smoking_current")
print_log_row("Sleep difficulties (Yes)", "sleep_diff")
print_log_row("High cholesterol", "high_chol")
print_eq_row("Non-insulin diabetes", "diabetes_treat_f", "non-insulin treated")
print_eq_row("Insulin-treated diabetes", "diabetes_treat_f", "insulin-treated")
print_log_row("Any diabetes treatment", "diabetes_any")
print_log_row("GPMP/TCA recorded", "gpmp_tca")
print_log_row("Self-report high BP at W1", "self_highBP_w1")
print_log_row("Heart disease", "heart_disease")
print_log_row("Stroke", "stroke_yn")

cat(strrep("-", nchar(header)), "\n", sep="")
cat(sprintf("%-30s | %12s %12s %12s %12s | %12s | %10s\n",
  "Incident HTN, n (%)",
  sprintf("%d (%.1f)", sum(adults$event[adults$bmi_cat=="Normal"]), 100*sum(adults$event[adults$bmi_cat=="Normal"])/sum(adults$bmi_cat=="Normal")),
  sprintf("%d (%.1f)", sum(adults$event[adults$bmi_cat=="Underweight"]), 100*sum(adults$event[adults$bmi_cat=="Underweight"])/sum(adults$bmi_cat=="Underweight")),
  sprintf("%d (%.1f)", sum(adults$event[adults$bmi_cat=="Overweight"]), 100*sum(adults$event[adults$bmi_cat=="Overweight"])/sum(adults$bmi_cat=="Overweight")),
  sprintf("%d (%.1f)", sum(adults$event[adults$bmi_cat=="Obese"]), 100*sum(adults$event[adults$bmi_cat=="Obese"])/sum(adults$bmi_cat=="Obese")),
  sprintf("%d (%.1f)", sum(adults$event), 100*sum(adults$event)/nrow(adults)),
  "<0.001"))

py_by <- tapply(adults$time_years, adults$bmi_cat, sum)
events_by <- tapply(adults$event, adults$bmi_cat, sum)
ir_by <- 1000 * events_by / py_by
cat(sprintf("%-30s | %12.0f %12.0f %12.0f %12.0f | %12.0f |\n",
  "Person-years",
  py_by["Normal"], py_by["Underweight"], py_by["Overweight"], py_by["Obese"], sum(adults$time_years)))
cat(sprintf("%-30s | %12.1f %12.1f %12.1f %12.1f | %12.1f |\n",
  "IR per 1000 PY",
  ir_by["Normal"], ir_by["Underweight"], ir_by["Overweight"], ir_by["Obese"],
  1000 * sum(adults$event) / sum(adults$time_years)))

table1_data <- list(N = counts, events = events_by, py = py_by, ir = ir_by,
                    total_n = nrow(adults), total_events = sum(adults$event),
                    total_py = sum(adults$time_years))
saveRDS(table1_data, data_path("table1_data.rds"))
cat("\nTable 1 data saved.\n")

