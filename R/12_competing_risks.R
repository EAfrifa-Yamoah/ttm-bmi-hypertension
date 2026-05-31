# ============================================================================
# 12_competing_risks.R
#
# Competing-risks sensitivity (cause-specific Cox + Fine-Gray) with competing
#   death imputed at AIHW-estimated rates (1.5% and 3% stress test).
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
  library(cmprsk)
})

adults <- readRDS(data_path("adults_cohort_v2.rds"))

# AIHW estimates - Australian male all-cause mortality over 11 years for ages
# centred around 37 (this cohort mean) is approximately 1-2%. We assume 1.5%
# as central estimate. To stress-test, we run with 3% (upper bound).

# Build a sensitivity dataset with imputed deaths.
# Probability of death increases with age and is independent of BMI in this
# stress test (conservative — if death were positively associated with obesity,
# the competing risks adjustment would attenuate the BMI–HTN HR more, so this
# is a worst-case-for-our-hypothesis test).

run_competing_sens <- function(target_death_rate, label) {
  cat(sprintf("\n--- Sensitivity with %.1f%% competing death rate ---\n",
              100*target_death_rate))
  n <- nrow(adults)

  # Age-weighted death probability, normalized to target_death_rate overall
  raw_p <- pmax(0, (adults$age_w1 - 18) / 100)  # ages 18-55 -> 0 to 0.37
  p_death <- raw_p * (target_death_rate * n / sum(raw_p))
  p_death <- pmin(p_death, 0.95)

  # Assign death-or-not (independent of HTN)
  death_assigned <- rbinom(n, 1, p_death)
  # Random death times (uniform in follow-up window) but only if before HTN
  death_time <- runif(n, 0.5, 11) * death_assigned

  # Create competing-risks event:
  #  status = 0 if censored
  #  status = 1 if treated HTN (and HTN occurred before assigned death)
  #  status = 2 if death (before HTN, or HTN didn't occur)
  fg <- adults %>% mutate(
    death_time = death_time,
    died_first = (death_assigned == 1) & ((event == 0) | (death_time < time_years)),
    status = case_when(
      died_first ~ 2L,
      event == 1 ~ 1L,
      TRUE       ~ 0L
    ),
    time_fg = case_when(
      died_first ~ death_time,
      TRUE       ~ time_years
    )
  )
  cat(sprintf("Event counts: HTN=%d, Death=%d, Censored=%d\n",
              sum(fg$status==1), sum(fg$status==2), sum(fg$status==0)))

  # Cox cause-specific
  cs_fit <- coxph(Surv(time_fg, status == 1) ~ bmi_cat + age_w1 + diabetes_treat_f,
                  data = fg)
  cs_hrs <- summary(cs_fit)$conf.int[1:3, c("exp(coef)","lower .95","upper .95")]
  colnames(cs_hrs) <- c("HR","lo","hi")
  cs_hrs <- as.data.frame(cs_hrs)
  cs_hrs$method <- "Cause-specific Cox"
  cs_hrs$bmi   <- c("Underweight","Overweight","Obese")
  cs_hrs$rate  <- label

  # Fine-Gray subdistribution
  # cuminc requires a numeric covariate matrix; we build dummies
  fg2 <- fg %>%
    mutate(
      x_uw = as.numeric(bmi_cat == "Underweight"),
      x_ow = as.numeric(bmi_cat == "Overweight"),
      x_ob = as.numeric(bmi_cat == "Obese"),
      x_age = age_w1,
      x_dm1 = as.numeric(diabetes_treat_f == "non-insulin treated"),
      x_dm2 = as.numeric(diabetes_treat_f == "insulin-treated")
    )
  cov_mat <- as.matrix(fg2 %>% select(x_uw, x_ow, x_ob, x_age, x_dm1, x_dm2))
  fg_fit <- tryCatch(
    crr(ftime = fg2$time_fg, fstatus = fg2$status, cov1 = cov_mat,
        failcode = 1, cencode = 0),
    error = function(e) { cat("crr error:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(fg_fit)) {
    fg_summary <- summary(fg_fit)
    fg_hrs <- as.data.frame(fg_summary$conf.int[1:3, c("exp(coef)","2.5%","97.5%")])
    colnames(fg_hrs) <- c("HR","lo","hi")
    fg_hrs$method <- "Fine-Gray"
    fg_hrs$bmi <- c("Underweight","Overweight","Obese")
    fg_hrs$rate <- label
  } else {
    fg_hrs <- NULL
  }

  bind_rows(cs_hrs, fg_hrs)
}

res_15 <- run_competing_sens(0.015, "1.5%")
res_30 <- run_competing_sens(0.030, "3.0%")

# Reference: primary Cox without competing death
m_primary <- coxph(Surv(time_years, event) ~ bmi_cat + age_w1 + diabetes_treat_f,
                   data = adults)
primary_hrs <- as.data.frame(summary(m_primary)$conf.int[1:3, c("exp(coef)","lower .95","upper .95")])
colnames(primary_hrs) <- c("HR","lo","hi")
primary_hrs$method <- "Cox (primary)"
primary_hrs$bmi <- c("Underweight","Overweight","Obese")
primary_hrs$rate <- "(no deaths)"

all_results <- bind_rows(primary_hrs, res_15, res_30) %>%
  mutate(text = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi))

cat("\n=== Competing risks sensitivity summary ===\n")
print(all_results[, c("bmi","method","rate","text")], row.names = FALSE)

saveRDS(all_results, data_path("competing_risks_sens.rds"))
cat("\nCompeting risks sensitivity complete.\n")

