# =============================================================================
# 05_dyd.R  —  main results: dyd(x), DYD(50), DYD(60), Kitagawa, stroke control
# Main scenario: main cause sets, R00-R99 redistributed, I50 kept,
#                life-table x cause-share rates, SIM 2019.
# Outputs (output/): dyd_age.csv, dyd_summary.csv, dyd_ratio_fm.csv,
#                    dyd_kitagawa.csv, stroke_did.csv, tipping_point.csv
# =============================================================================
source("R/00_setup.R")

prev   <- read_csv(file.path(dir_out, "prevalence.csv"), show_col_types = FALSE)
counts <- readRDS(file.path(dir_data, "sim_counts_2019.rds"))
lt_g   <- readRDS(file.path(dir_data, "lt_groups_2019.rds"))

#' Build the input table (prevalence + rates) for one outcome and scenario
build_input <- function(outcome, set, counts, lt_g, ...) {
  r <- cause_rates(counts, lt_g, set, ...)
  p <- prev[prev$outcome == outcome, c("sex", "agegrp", "pi", "se_pi")]
  d <- merge(p, r, by = c("sex", "agegrp"))
  d$sex <- as.character(d$sex); d$agegrp <- as.character(d$agegrp)
  d[order(d$sex, d$age_lo), ]
}

run_main <- function(outcome) {
  d <- build_input(outcome, outcome_sets[[outcome]]["main"], counts, lt_g,
                   illdef = "prop", i50 = "keep", method = "lt_fraction")
  list(
    age     = cbind(outcome = outcome, dyd_age(d)),
    summary = cbind(outcome = outcome, rbind(dyd_summary(d, 50), dyd_summary(d, 60))),
    ratio   = cbind(outcome = outcome, dyd_ratio_fm(d)),
    kit     = cbind(outcome = outcome, rbind(dyd_kitagawa(d, 50), dyd_kitagawa(d, 60)))
  )
}

h <- run_main("heart"); s <- run_main("stroke")

write_csv(rbind(h$age, s$age),         file.path(dir_out, "dyd_age.csv"))
write_csv(rbind(h$summary, s$summary), file.path(dir_out, "dyd_summary.csv"))
write_csv(rbind(h$ratio, s$ratio),     file.path(dir_out, "dyd_ratio_fm.csv"))
write_csv(rbind(h$kit, s$kit),         file.path(dir_out, "dyd_kitagawa.csv"))

did <- stroke_did(h$ratio[, -1], s$ratio[, -1])
write_csv(did, file.path(dir_out, "stroke_did.csv"))

# Tipping point (K gap unknown: 0 = lower bound on the survival gap needed)
tp <- with(h$kit, data.frame(age_from, DYD_F, DYD_M,
                             min_survival_gap_years = tipping_point(DYD_F, DYD_M)))
write_csv(tp, file.path(dir_out, "tipping_point.csv"))

print(h$summary); print(s$summary); print(h$kit); print(did)
