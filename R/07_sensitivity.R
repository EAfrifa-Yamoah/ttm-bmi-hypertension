# ============================================================================
# 07_sensitivity.R
#
# Sensitivity analyses: all-ages cohort, age strata, 2-year landmark,
#   90-day event exclusion, and combined PBS-or-self-report outcome.
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

adults    <- readRDS(data_path("adults_cohort_v2.rds"))
analytic  <- readRDS(data_path("analytic_cohort.rds"))  # includes adolescents

# helper to extract BMI HRs
extract_bmi <- function(model, label) {
  s <- summary(model)
  ix <- grep("bmi_cat", rownames(s$coefficients))
  data.frame(
    analysis = label,
    var = rownames(s$coefficients)[ix],
    HR  = s$conf.int[ix, "exp(coef)"],
    lo  = s$conf.int[ix, "lower .95"],
    hi  = s$conf.int[ix, "upper .95"],
    p   = s$coefficients[ix, "Pr(>|z|)"]
  )
}

sens_results <- list()

# ---- Sensitivity 1: Primary cohort (Adults, M2-equivalent) ----
m_primary <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f,
                   data = adults)
sens_results$primary <- extract_bmi(m_primary, "Primary (Adults, M2)")
cat("Primary:\n"); print(sens_results$primary)

# ---- Sensitivity 2: All ages (includes 10-17) ----
analytic$bmi_cat <- relevel(factor(analytic$bmi_cat,
                                   levels = c("Normal","Underweight","Overweight","Obese")),
                            ref = "Normal")
analytic$diabetes_treat_f <- factor(analytic$diabetes_treat_f,
                                    levels = c("no diabetes","non-insulin treated","insulin-treated"))
m_all <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f,
               data = analytic)
sens_results$all_ages <- extract_bmi(m_all, "All ages (10-55)")
cat("\nAll ages:\n"); print(sens_results$all_ages)

# ---- Sensitivity 3: Age-stratified ----
adults$age_grp <- cut(adults$age_w1, breaks = c(17, 29, 44, 100),
                       labels = c("18-29","30-44","45+"))

for (g in levels(adults$age_grp)) {
  d <- adults %>% filter(age_grp == g)
  if (sum(d$event) < 20) next
  m <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f, data = d)
  res <- extract_bmi(m, sprintf("Age %s (n=%d, ev=%d)", g, nrow(d), sum(d$event)))
  sens_results[[paste0("age_",g)]] <- res
  cat("\n", g, "(n=", nrow(d), ", events=", sum(d$event), "):\n", sep="")
  print(res)
}

# ---- Sensitivity 4: Excluding first 2 years (latency / reverse causation) ----
late <- adults %>%
  filter(time_years >= 2 | event == 0)  # keep censored, drop early events
# Actually proper landmark: shift baseline to year 2
# Approach: among those still event-free at year 2, refit
landmark2 <- adults %>%
  filter(event == 0 | time_years >= 2) %>%
  mutate(
    event_lmk      = ifelse(time_years >= 2, event, 0),
    time_years_lmk = pmax(time_years - 2, 0)
  ) %>%
  filter(time_years_lmk > 0)

m_late <- coxph(Surv(time_years_lmk, event_lmk) ~ bmi_cat + age_w1 + diabetes_treat_f,
                data = landmark2)
sens_results$landmark2 <- extract_bmi(m_late, "Landmark at 2y")
cat("\nLandmark at 2y:\n"); print(sens_results$landmark2)

# ---- Sensitivity 5: Stricter hypertension definition ----
# Use highbpmedpbs which has timing strata. Keep only those whose treatment was
# initiated in 2014-2015 or 2016-2021 (i.e., during follow-up) as confirmed events
# This is a stricter operationalisation
adults$highbpmedpbs_strict <- adults$event
# Actually antihypertensive_treatment_initiation_date == valid date is the same as event
# We can use a sensitivity restricting to those whose initiation was post-baseline by >= 90d
# (excluding very-early dispensings that might reflect existing condition)
adults_strict <- adults %>%
  filter(event == 0 | time_years > (90/365.25))
m_strict <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f,
                  data = adults_strict)
sens_results$strict_90d <- extract_bmi(m_strict, "Excl events <90d (treated HTN)")
cat("\nExcluding events within 90d of baseline:\n"); print(sens_results$strict_90d)

# ---- Sensitivity 6: Alternative outcome - self-reported high BP ever during FU ----
# Build new outcome: any self-reported high BP at W2 or W3 (or PBS event)
raw <- read.csv(RAW_CSV, stringsAsFactors = FALSE)
raw$wave_int <- as.integer(raw$wave)
selfBP <- raw %>%
  group_by(zdcid0001md) %>%
  summarise(self_BP_followup = any(highBP == "Yes" & wave_int > 1, na.rm = TRUE),
            .groups = "drop")
adults_self <- adults %>% left_join(selfBP, by = "zdcid0001md") %>%
  mutate(combined_event = pmax(event, as.integer(self_BP_followup), na.rm = TRUE))
cat("\nCombined outcome (PBS event OR self-report high BP at W2/W3):\n")
cat("Events:", sum(adults_self$combined_event), "/", nrow(adults_self), "\n")
# Use simple binary outcome (no time at risk for self-report); use logistic-style approach
# Or just refit with combined event as a coarse "ever-HTN" indicator
# For Cox: keep time_years as-is, expand event to combined
m_self <- coxph(Surv(time_years, combined_event) ~ bmi_cat + age_w1 + diabetes_treat_f,
                data = adults_self)
sens_results$self_report <- extract_bmi(m_self, "Self-rep HTN combined outcome")
cat("\nSelf-report combined outcome:\n"); print(sens_results$self_report)

# ---- Save ----
saveRDS(sens_results, data_path("sensitivity_results.rds"))

# ---- Combined summary table ----
sens_df <- bind_rows(sens_results)
sens_df$bmi_level <- gsub("bmi_cat","", sens_df$var)
sens_df$text <- sprintf("%.2f (%.2f-%.2f)", sens_df$HR, sens_df$lo, sens_df$hi)
sens_summary <- sens_df %>%
  select(analysis, bmi_level, text, p) %>%
  mutate(p = signif(p, 3))

cat("\n=== Sensitivity analysis summary ===\n")
print(sens_summary, row.names = FALSE)

saveRDS(sens_summary, data_path("sensitivity_summary.rds"))
cat("\nSensitivity analyses complete.\n")

