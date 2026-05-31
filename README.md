# Baseline BMI Category and Incident Treated Hypertension in Australian Men

Analysis code for a prospective cohort study of the association between baseline
body mass index (BMI) category and incident treated hypertension in Australian
men, using the **Ten to Men Australian Longitudinal Study on Male Health** linked
with **Pharmaceutical Benefits Scheme (PBS)** and **Medicare Benefits Schedule
(MBS)** records.

The pipeline estimates the BMI–hypertension association with progressive
confounder adjustment, quantifies incremental discrimination and clinical utility,
and stress-tests the findings against detection bias, missing data, latency,
effect modification by age, and competing risks.

---

## Summary of the analysis

- **Design:** prospective cohort, 11-year follow-up (baseline 1 Jan 2014, censor 28 Feb 2025)
- **Primary cohort:** 12,742 men aged ≥18, free of treated hypertension at baseline, with a valid Wave 1 BMI category
- **Exposure:** baseline BMI category (Underweight, Normal [reference], Overweight, Obese)
- **Outcome:** incident treated hypertension (first PBS dispensing of an antihypertensive, ATC C02/C03/C07/C08/C09)
- **Primary model:** Cox proportional hazards with progressive adjustment
- **Key secondary analyses:** bootstrap C-index, decision curve analysis, calibration, multiple imputation, E-values, BMI × age interaction, competing risks

### Headline results

| BMI category | Fully adjusted HR (95% CI) | E-value (CI bound) |
|---|---|---|
| Underweight | 1.43 (1.16–1.77) | 2.21 (1.59) |
| Overweight  | 1.42 (1.22–1.67) | 2.19 (1.74) |
| Obese       | 1.88 (1.59–2.22) | 3.17 (2.56) |

- Incremental discrimination of BMI over age alone: ΔC = 0.028 (95% CI 0.020–0.038)
- Fully adjusted model C-index: 0.754 (95% CI 0.740–0.769)
- BMI × age interaction: χ² = 12.79, df = 3, **p = 0.005**
- 10-year cumulative incidence: 7.5% (Normal) → 20.5% (Obese)

---

## Repository structure

```
ttm-bmi-hypertension/
├── README.md                  # this file
├── LICENSE                    # MIT licence (code)
├── CITATION.cff               # how to cite this repository
├── renv.lock                  # pinned package versions (optional, see Reproducibility)
├── .gitignore
├── R/
│   ├── config.R               # central paths, constants, seed
│   ├── 00_run_all.R           # master runner (sources 01–13 in order)
│   ├── 01_data_prep.R         # cohort construction + outcome derivation
│   ├── 02_table1.R            # Table 1 descriptive statistics
│   ├── 03_cox_primary.R       # Cox models M1–M3 + Schoenfeld test
│   ├── 04_cindex_bootstrap.R  # C-index with bootstrap 95% CIs
│   ├── 05_calibration_dca.R   # calibration + DCA (point estimates)
│   ├── 06_figures.R           # Figures 1–3 (+ superseded 4,5)
│   ├── 07_sensitivity.R       # age strata, landmark, 90-day, combined outcome
│   ├── 08_multiple_imputation.R   # MICE (m=20) pooled Cox
│   ├── 09_interaction_evalue_cuminc.R  # interaction, E-values, cumulative incidence
│   ├── 10_schoenfeld_plot.R   # Figure 6 (Schoenfeld residuals)
│   ├── 11_dca_bootstrap.R     # Figure 4 (DCA with bootstrap ribbons)
│   ├── 12_competing_risks.R   # cause-specific Cox + Fine–Gray
│   └── 13_calibration_deciles.R   # Figure 5 (decile calibration)
├── data/
│   └── README.md              # data access + the UTF-8 conversion step
├── figures/                   # generated figures (git-ignored except .gitkeep)
└── manuscript/
    └── README.md              # pointers to the manuscript + figure mapping
```

---

## Requirements

- **R ≥ 4.3**
- CRAN packages: `dplyr`, `tidyr`, `readr`, `survival`, `rms`, `boot`, `ggplot2`,
  `survminer`, `cowplot`, `mice`, `cmprsk`, `scales`

Install them with:

```r
install.packages(c(
  "dplyr", "tidyr", "readr", "survival", "rms", "boot", "ggplot2",
  "survminer", "cowplot", "mice", "cmprsk", "scales"
))
```

(See `R/install_packages.R` for a one-shot installer.)

---

## How to run

1. **Obtain the data.** This study uses linked Ten to Men / PBS / MBS data, which
   are **not** publicly distributable (see [`data/README.md`](data/README.md) and
   the Data Availability statement). Place your UTF-8 converted input file at:

   ```
   data/ttmcombined_utf8.csv
   ```

2. **Set the working directory to the repository root** (or open the project in
   RStudio at the root).

3. **Run the whole pipeline:**

   ```bash
   Rscript R/00_run_all.R
   ```

   or, interactively:

   ```r
   source("R/00_run_all.R")
   ```

   Intermediates (`.rds`) are written to `data/`; figures (`.png`) to `figures/`.

You can also run any single script (after its prerequisites have been run once),
e.g. `Rscript R/04_cindex_bootstrap.R`.

### Approximate runtime

Most scripts run in seconds. The bootstrap steps dominate:
`04_cindex_bootstrap.R` (~1–2 min, B = 500) and `11_dca_bootstrap.R` (~1 min,
B = 150). The full pipeline completes in well under 10 minutes on a typical laptop.

---

## Reproducibility

- A global seed (`SEED = 2026`) is set in `R/config.R` and applied before every
  stochastic step (bootstrap, multiple imputation, competing-risks imputation).
- All analysis constants (dates, horizon, ATC codes, bootstrap replicates) live
  in `R/config.R` — change them in one place.
- An `renv.lock` is provided for exact package-version pinning; run
  `renv::restore()` to reconstruct the environment. (Optional — the pipeline runs
  on any R ≥ 4.3 with the package versions listed above.)

---

## Data availability

The Ten to Men data are managed by the Australian Institute of Family Studies and
are available to researchers under a data access agreement. PBS/MBS linkage data
are provided under separate approvals. None of these data are included in this
repository. See [`data/README.md`](data/README.md) for the access process and the
one-line encoding conversion required before running the pipeline.

---

## Citation

If you use this code, please cite the associated manuscript and this repository.
See [`CITATION.cff`](CITATION.cff).

## Licence

Code is released under the MIT Licence ([`LICENSE`](LICENSE)). The licence applies
to the analysis code only, not to the underlying data.
