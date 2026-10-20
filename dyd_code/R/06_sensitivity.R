# =============================================================================
# 06_sensitivity.R  —  which assumptions move DYD(60) and the F/M gap?
# Each scenario changes ONE assumption relative to the main scenario.
# Output: output/sensitivity.csv
# =============================================================================
source("R/00_setup.R")
source("R/05_dyd.R")      # reuses prev, counts, lt_g, build_input()

counts_1819 <- readRDS(file.path(dir_data, "sim_counts_2018_2019.rds"))
f_pop <- file.path(dir_data, "pop_groups_2019.rds")
pop_g <- if (file.exists(f_pop)) readRDS(f_pop) else NULL

scenarios <- list(
  main          = list(set = "main",  illdef = "prop", i50 = "keep",    method = "lt_fraction", counts = counts),
  broad_icd     = list(set = "broad", illdef = "prop", i50 = "keep",    method = "lt_fraction", counts = counts),
  no_illdef_red = list(set = "main",  illdef = "none", i50 = "keep",    method = "lt_fraction", counts = counts),
  i50_garbage   = list(set = "main",  illdef = "prop", i50 = "garbage", method = "lt_fraction", counts = counts),
  sim_over_pop  = list(set = "main",  illdef = "prop", i50 = "keep",    method = "sim_pop",     counts = counts),
  pooled_18_19  = list(set = "main",  illdef = "prop", i50 = "keep",    method = "lt_fraction", counts = counts_1819)
)

res <- list()
for (sc in names(scenarios)) {
  p <- scenarios[[sc]]
  if (p$method == "sim_pop" && is.null(pop_g)) { message("Skipping ", sc, ": no population file"); next }
  for (outcome in c("heart", "stroke")) {
    # (the I50 scenario only changes heart disease; it runs for stroke too as a check)
    d <- build_input(outcome, outcome_sets[[outcome]][p$set], p$counts, lt_g,
                     illdef = p$illdef, i50 = p$i50, method = p$method, pop = pop_g)
    k <- dyd_kitagawa(d, 60)
    sm <- dyd_summary(d, 60)
    res[[length(res) + 1]] <- data.frame(
      scenario = sc, outcome = outcome,
      DYD60_F = k$DYD_F, DYD60_M = k$DYD_M, ratio_FM = k$DYD_F / k$DYD_M,
      DYD60_F_common = k$DYD_F_common, DYD60_M_common = k$DYD_M_common,
      share_gap_same_ages = k$share_rates,
      F_lo = sm$lo[sm$sex == "Female"], F_hi = sm$hi[sm$sex == "Female"],
      M_lo = sm$lo[sm$sex == "Male"],   M_hi = sm$hi[sm$sex == "Male"])
  }
}
sens <- do.call(rbind, res)
print(sens, digits = 3)
write_csv(sens, file.path(dir_out, "sensitivity.csv"))
