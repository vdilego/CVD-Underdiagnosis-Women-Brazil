# =============================================================================
# 00_setup.R  —  packages, paths and ALL analysis parameters
# DYD project: diagnosed years per CVD death, Brazil 2019 (PNS 2019 + SIM 2019)
#
# Run this first in every session. Every other script sources it.
# Anything that may need correcting (PNS variable codes, ICD sets, ages)
# is defined HERE, once.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
# install.packages(c("survey", "dplyr", "tidyr", "readr", "ggplot2", "PNSIBGE"))
# remotes::install_github("rfsaldanha/microdatasus")
pkgs <- c("survey", "dplyr", "tidyr", "readr", "ggplot2")
invisible(lapply(pkgs, library, character.only = TRUE))
# PNSIBGE and microdatasus are called with :: in the scripts that need them.

# Lonely PSUs after subsetting to 50+ (standard choice for PNS)
options(survey.lonely.psu = "adjust")

# ---- Paths ------------------------------------------------------------------
dir_data <- "data"      # raw and intermediate data
dir_out  <- "output"    # tables (csv) and figures (png)
dir.create(dir_data, showWarnings = FALSE)
dir.create(dir_out,  showWarnings = FALSE)

source("R/functions_dyd.R")   # all indicator maths (base R, tested)

# ---- Ages -------------------------------------------------------------------
age_breaks <- c(50, 60, 70, 80, Inf)
age_labels <- c("50-59", "60-69", "70-79", "80+")
age_lo     <- c(50, 60, 70, 80)          # lower bound of each group
min_pos_cases <- 30                       # flag cells with fewer positives

# ---- Year(s) ---------------------------------------------------------------
sim_years_main <- 2019
sim_years_alt  <- 2018:2019               # sensitivity: pooled deaths

# ---- PNS 2019 variable codes  (VERIFY EACH ONE IN THE 2019 DICTIONARY) ------
# Left: name used in our code. Right: code in the PNS 2019 microdata.
pns_vars <- c(
  sex          = "C006",    # 1 man, 2 woman
  age          = "C008",    # years
  race         = "C009",    # 1 white 2 black 3 yellow 4 brown 5 indigenous 9 ign.
  educ         = "VDD004A", # 1..7 highest level completed
  uf           = "V0001",   # state code (first digit = macro-region)
  urban        = "V0026",   # 1 urban, 2 rural
  heart        = "Q06306",  # ever told heart disease (infarction, angina, HF, other)
  age_dx_heart = "Q06307",  # age at first heart disease diagnosis   (VERIFY)
  stroke       = "Q068",    # ever told stroke                       (VERIFY)
  age_dx_stroke= "Q069",    # age at first stroke diagnosis          (VERIFY)
  hyp          = "Q00201",  # hypertension: 1 yes, 2 only in pregnancy, 3 no
  diab         = "Q03001",  # diabetes:     1 yes, 2 only in pregnancy, 3 no
  chol         = "Q06001",  # high cholesterol: 1 yes                (VERIFY)
  smoke        = "P050",    # currently smokes: 1 daily, 2 less than daily, 3 no
  exercise     = "P034",    # any exercise in last 3 months: 1 yes, 2 no
  weight       = "P00104",  # self-reported weight, kg               (VERIFY)
  height       = "P00404",  # self-reported height, cm               (VERIFY)
  visit12      = "J01101",  # last doctor visit: 1 = within 12 months (VERIFY)
  plan         = "I00102",  # has private health plan: 1 yes         (VERIFY)
  esf          = "B001"     # household registered in Family Health Strategy (VERIFY)
)

# ---- ICD-10 cause sets (underlying cause) ------------------------------------
# Each set = list of c(from, to) three-character ranges.
cause_sets <- list(
  all          = list(c("A00", "Z99")),
  illdef       = list(c("R00", "R99")),
  cvd_all      = list(c("I00", "I99")),
  ihd          = list(c("I20", "I25")),
  i50          = list(c("I50", "I50")),
  heart_main   = list(c("I20", "I25"), c("I50", "I50")),
  heart_broad  = list(c("I05", "I09"), c("I11", "I11"), c("I13", "I13"),
                      c("I20", "I51")),
  stroke_main  = list(c("I60", "I64")),
  stroke_broad = list(c("I60", "I69"))
)

# Which PNS outcome goes with which death sets
outcome_sets <- list(
  heart  = c(main = "heart_main",  broad = "heart_broad"),
  stroke = c(main = "stroke_main", broad = "stroke_broad")
)
