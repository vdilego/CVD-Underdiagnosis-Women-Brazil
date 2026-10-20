# =============================================================================
# 07_riskfactors.R  —  is there a diagnostic paradox? (PNS 2019, ages 50+)
# 1. Burden: diagnosed risk factors by sex and age
# 2. Diagnosis given risk: prevalence ratios (women vs men) for heart disease
#    and stroke, adjusted step by step (Poisson, log link, design-based SEs)
# 3. Decomposition of the sex gap in diagnosed prevalence
# Outputs (output/): rf_burden.csv, rf_prevalence_ratios.csv,
#                    rf_heart_by_profile.csv, rf_decomposition.csv
# Reminder: risk factors here are self-reported DIAGNOSES (no measured BP in
# 2019). So this compares women and men already known to the health system.
# =============================================================================
source("R/00_setup.R")
des50 <- readRDS(file.path(dir_data, "pns2019_design_50plus.rds"))

rfs <- c("hyp", "diab", "chol", "smoke", "obese", "inactive")

# ---- 1. Burden ---------------------------------------------------------------
burden <- do.call(rbind, lapply(rfs, function(v) {
  r <- svyby(as.formula(paste0("~", v)), ~sex + agegrp, des50, svymean,
             vartype = c("se", "ci"), na.rm = TRUE)
  data.frame(risk_factor = v, sex = r$sex, agegrp = r$agegrp, prev = r[[v]],
             se = r$se, lo = r$ci_l, hi = r$ci_u, row.names = NULL)
}))
# Number of risk factors (0, 1, 2, 3+)
nrf <- svyby(~rf_cat, ~sex + agegrp, subset(des50, !is.na(rf_cat)), svymean)
write_csv(burden, file.path(dir_out, "rf_burden.csv"))
write_csv(as.data.frame(nrf), file.path(dir_out, "rf_count_distribution.csv"))

# Also descriptive: missingness of each risk factor by sex (check BMI!)
miss <- sapply(rfs, function(v) tapply(is.na(des50$variables[[v]]), des50$variables$sex, mean))
print(round(100 * miss, 1))

# ---- 2. Prevalence ratios -----------------------------------------------------
# Models (same sample in all: complete cases on every covariate)
covs <- c("female", "agegrp", rfs, "visit12", "plan", "esf",
          "educ", "race", "region", "urban")
dcc <- subset(des50, complete.cases(des50$variables[, covs]))
cat("Complete-case n:", nrow(dcc$variables), "of", nrow(des50$variables), "\n")

models <- list(
  M1_age      = "female + agegrp",
  M2_risk     = "female + agegrp + hyp + diab + chol + smoke + obese + inactive",
  M3_care     = "female + agegrp + hyp + diab + chol + smoke + obese + inactive + visit12 + plan + esf",
  M4_full     = "female + agegrp + hyp + diab + chol + smoke + obese + inactive + visit12 + plan + esf + educ + race + region + urban",
  M5_by_age   = "female * agegrp + hyp + diab + chol + smoke + obese + inactive + visit12 + plan + esf + educ + race + region + urban"
)

fit_pr <- function(outcome, rhs) {
  m <- svyglm(as.formula(paste(outcome, "~", rhs)), design = dcc,
              family = quasipoisson(link = "log"))
  b <- coef(m); se <- sqrt(diag(vcov(m)))
  data.frame(term = names(b), PR = exp(b), lo = exp(b - 1.96 * se),
             hi = exp(b + 1.96 * se), row.names = NULL)
}

pr <- list()
for (outcome in c("heart", "stroke")) for (mn in names(models)) {
  t <- fit_pr(outcome, models[[mn]])
  t <- t[grepl("^female", t$term), ]       # keep the sex terms only
  pr[[length(pr) + 1]] <- cbind(outcome = outcome, model = mn, t)
}
pr <- do.call(rbind, pr)
print(pr, digits = 3)
write_csv(pr, file.path(dir_out, "rf_prevalence_ratios.csv"))
# In M5: PR for women at age group k = exp(b_female + b_female:agegrp_k).
# Use svycontrast() on the M5 fit if you need intervals for each age group.

# ---- 3. Diagnosis given risk: profiles + decomposition --------------------------
# Profile = combination of diagnosed hypertension, diabetes, high cholesterol
# (8 profiles, e.g. "101" = hyp yes, diab no, chol yes). Sensitivity: rf_cat.
decomp_by <- function(outcome, profile_var, ages = NULL) {
  keep <- !is.na(dcc$variables[[profile_var]])
  if (!is.null(ages)) keep <- keep & dcc$variables$agegrp %in% ages
  dd <- dcc[keep, ]                     # keeps the full design for variances
  f_p  <- as.formula(paste0("~", profile_var))
  share <- svyby(f_p, ~sex, dd, svymean)                 # p_s(r)
  piR   <- svyby(as.formula(paste0("~", outcome)),
                 as.formula(paste0("~sex + ", profile_var)), dd, svymean)
  # reshape shares to long
  sh_cols <- grep(paste0("^", profile_var), names(share), value = TRUE)
  sh <- do.call(rbind, lapply(sh_cols, function(cn) data.frame(
    sex = share$sex, profile = sub(paste0("^", profile_var), "", cn), p = share[[cn]])))
  pi_df <- data.frame(sex = piR$sex, profile = as.character(piR[[profile_var]]),
                      pi = piR[[outcome]])
  x <- merge(sh, pi_df, by = c("sex", "profile"))
  list(table = x, decomp = cbind(outcome = outcome, profile_var = profile_var,
                                 ages = if (is.null(ages)) "50+" else paste(ages, collapse = ","),
                                 rf_decomp(x)))
}

dec <- list(); tabs <- list()
for (outcome in c("heart", "stroke")) for (pv in c("profile", "rf_cat")) {
  r <- decomp_by(outcome, pv)
  dec[[length(dec) + 1]] <- r$decomp
  tabs[[length(tabs) + 1]] <- cbind(outcome = outcome, profile_var = pv, r$table)
  for (ag in age_labels) {                     # by age group
    r2 <- decomp_by(outcome, pv, ages = ag)
    dec[[length(dec) + 1]] <- r2$decomp
  }
}
write_csv(do.call(rbind, tabs), file.path(dir_out, "rf_heart_by_profile.csv"))
dec <- do.call(rbind, dec)
print(dec, digits = 3)
write_csv(dec, file.path(dir_out, "rf_decomposition.csv"))
# Reading: if risk_composition > 0 (women's profiles predict MORE diagnoses)
# while diagnosis_given_risk < 0 for heart disease but not for stroke,
# the diagnostic paradox holds.
