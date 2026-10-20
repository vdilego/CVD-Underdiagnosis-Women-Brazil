# DYD project — R code

Diagnosed years per CVD death (DYD), Brazil 2019, ages 50+, PNS 2019 + SIM 2019.

Open the folder as an RStudio project (or `setwd()` into it). Run the scripts in order:

| Script | What it does | Main output |
|---|---|---|
| `R/00_setup.R` | Packages, paths, **all parameters** (PNS codes, ICD sets, ages) | — |
| `R/functions_dyd.R` | Indicator maths (base R, tested) | — |
| `R/01_pns_prepare.R` | Downloads PNS 2019, recodes, 50+ design, cell counts | `output/pns_cell_counts.csv` |
| `R/02_pns_prevalence.R` | Prevalence of heart disease and stroke by sex × age; age at diagnosis | `output/prevalence.csv` |
| `R/03_sim_prepare.R` | Downloads SIM, deaths by sex × age × cause set | `data/sim_counts_2019.rds` |
| `R/04_lifetable.R` | IBGE 2019 life tables → age groups | `data/lt_groups_2019.rds` |
| `R/05_dyd.R` | dyd(x), DYD(50/60), Kitagawa, stroke control, tipping point | `output/dyd_*.csv` |
| `R/06_sensitivity.R` | Six scenarios | `output/sensitivity.csv` |
| `R/07_riskfactors.R` | Burden, prevalence ratios, decomposition | `output/rf_*.csv` |
| `R/08_figures.R` | Figures for the slides | `output/fig*.png` |

Before anything else:

1. Check every code in `pns_vars` (`00_setup.R`) against the PNS 2019 dictionary.
2. Put `data/ibge_lt_2019.csv` in place (columns `sex, age, mx, Lx, lx`; sex 1 = male, 2 = female).
3. Run `Rscript tests/test_functions.R` and `Rscript tests/test_pipeline.R`: every line must say ok.

Packages: `survey`, `dplyr`, `tidyr`, `readr`, `ggplot2`, `PNSIBGE` (CRAN) and `microdatasus` (`remotes::install_github("rfsaldanha/microdatasus")`).

Notes
- `share_rates` in the Kitagawa table can fall outside 0–1 when the two parts have opposite signs; report both parts in years.
- DYD intervals ignore the covariance of PNS estimates across age groups; fine for GSA, revisit with replicate weights for the paper.
