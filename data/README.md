# Data

This directory holds the input data and the derived intermediate objects
produced by the pipeline. **No data files are committed to the repository.**

## Required input

The pipeline expects a single UTF-8 encoded CSV at:

```
data/ttmcombined_utf8.csv
```

This is the linked Ten to Men / PBS / MBS analysis file, with one row per
participant-wave. Key columns used by the pipeline include:

| Column | Description |
|---|---|
| `zdcid0001md` | participant identifier |
| `wave` | survey wave (1, 2, 3) |
| `age` | age at wave |
| `BMI_category` | Underweight / Normal / Overweight / Obese (Wave 1 only) |
| `antihypertensive_any` | ever dispensed an antihypertensive |
| `first_htn_rx_date` | date of first antihypertensive dispensing |
| `diabetes_treat` | no diabetes / non-insulin treated / insulin-treated |
| `physical_activity` | meets / does not meet guidelines |
| `smoking` | Current / Former / Non-smoker |
| `high_cholesterol`, `sleep_issues`, `gpmp_tca_any`, `highBP`, ... | covariates |

## Encoding conversion

The source extract is Latin-1 encoded and contains non-ASCII bytes (e.g. in the
income field) that break R's CSV reader. Convert once before running the pipeline:

```bash
iconv -f LATIN1 -t UTF-8 ttmcombined.csv > data/ttmcombined_utf8.csv
```

## Data access

The Ten to Men Australian Longitudinal Study on Male Health is managed by the
Australian Institute of Family Studies (AIFS). Data are available to approved
researchers under a data access agreement: https://ten2men.org.au

PBS and MBS linkage data are provided by Services Australia under separate ethics
and data-custodian approvals. Neither the raw data nor the linked extract can be
redistributed; researchers must obtain their own approvals.

## Derived files (generated, git-ignored)

Running the pipeline writes these intermediates here:

```
analytic_cohort.rds          adults_cohort.rds        adults_cohort_v2.rds
youth_cohort.rds             table1_data.rds          cox_models_adults.rds
cindex_bootstrap.rds         calibration.rds          dca_results.rds
sensitivity_results.rds      sensitivity_summary.rds  mi_results.rds
interaction_results.rds      evalues.rds              cumulative_incidence.rds
schoenfeld.rds               schoenfeld_bmi.rds       dca_bootstrap.rds
competing_risks_sens.rds     calibration_deciles.rds
```

These are reproducible from the raw input and are intentionally excluded from
version control (see `.gitignore`).
