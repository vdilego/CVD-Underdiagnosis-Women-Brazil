# =============================================================================
# 02_pns_prevalence.R  —  diagnosed prevalence by sex x age group (PNS 2019)
# Output: output/prevalence.csv  (outcome, sex, agegrp, pi, se_pi, lo, hi)
#         output/age_at_diagnosis.csv (descriptive only)
# =============================================================================
source("R/00_setup.R")
des50 <- readRDS(file.path(dir_data, "pns2019_design_50plus.rds"))

prev_one <- function(outcome) {
  f <- as.formula(paste0("~", outcome))
  r <- svyby(f, ~sex + agegrp, des50, svymean,
             vartype = c("se", "ci"), na.rm = TRUE)
  data.frame(outcome = outcome, sex = r$sex, agegrp = r$agegrp,
             pi = r[[outcome]], se_pi = r$se, lo = r$ci_l, hi = r$ci_u,
             row.names = NULL)
}

prev <- rbind(prev_one("heart"), prev_one("stroke"))
prev$age_lo <- age_lo[match(prev$agegrp, age_labels)]
print(prev)
write_csv(prev, file.path(dir_out, "prevalence.csv"))

# ---- Age at first diagnosis (descriptive; survivors only) -------------------
# Among people alive in 2019 who report a diagnosis. Favours long durations:
# do NOT read as age at onset or as a delay estimate.
agedx <- rbind(
  data.frame(outcome = "heart",
             svyby(~age_dx_heart, ~sex + agegrp,
                   subset(des50, heart == 1 & !is.na(age_dx_heart)),
                   svymean, vartype = "se"), row.names = NULL) |>
    rename(mean_age_dx = age_dx_heart),
  data.frame(outcome = "stroke",
             svyby(~age_dx_stroke, ~sex + agegrp,
                   subset(des50, stroke == 1 & !is.na(age_dx_stroke)),
                   svymean, vartype = "se"), row.names = NULL) |>
    rename(mean_age_dx = age_dx_stroke)
)
print(agedx)
write_csv(agedx, file.path(dir_out, "age_at_diagnosis.csv"))
