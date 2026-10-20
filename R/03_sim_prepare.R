# =============================================================================
# 03_sim_prepare.R  —  SIM deaths by sex, age group and cause set
# Output: data/sim_counts_2019.rds and data/sim_counts_2018_2019.rds
#   one row per sex x agegrp; columns = death counts:
#   all, illdef, cvd_all, i50, ihd, heart_main, heart_broad, stroke_main,
#   stroke_broad, and <set>_noI50 for every set (needed for the I50 scenario)
# =============================================================================
source("R/00_setup.R")

sim_counts <- function(years) {
  f_raw <- file.path(dir_data, paste0("sim_raw_", paste(range(years), collapse = "_"), ".rds"))
  if (!file.exists(f_raw)) {
    raw <- microdatasus::fetch_datasus(year_start = min(years), year_end = max(years),
                                       information_system = "SIM-DO",
                                       vars = c("CAUSABAS", "SEXO", "IDADE"))
    saveRDS(raw, f_raw)
  }
  raw <- readRDS(f_raw)

  d <- data.frame(
    sex   = factor(as.character(raw$SEXO), levels = c("1", "2"),
                   labels = c("Male", "Female")),
    age   = sim_age_years(raw$IDADE),
    code3 = icd3(raw$CAUSABAS)
  )
  # Report and drop records we cannot classify
  cat("Records:", nrow(d), "\n",
      "Missing/ignored sex:", sum(is.na(d$sex)), "\n",
      "Missing age:", sum(is.na(d$age)), "\n",
      "Invalid cause code:", sum(is.na(d$code3)), "\n")
  d <- d[!is.na(d$sex) & !is.na(d$age), ]
  d$agegrp <- cut(d$age, age_breaks, right = FALSE, labels = age_labels)
  d <- d[!is.na(d$agegrp), ]

  is_i50 <- in_set(d$code3, cause_sets$i50)
  sets <- setdiff(names(cause_sets), "all")
  for (s in sets) {
    d[[s]] <- as.integer(in_set(d$code3, cause_sets[[s]]))
    d[[paste0(s, "_noI50")]] <- as.integer(d[[s]] == 1 & !is_i50)
  }
  d$all <- 1L   # every death with valid sex and age (invalid cause codes kept)

  count_cols <- c("all", sets, paste0(sets, "_noI50"))
  out <- aggregate(d[count_cols], by = list(sex = d$sex, agegrp = d$agegrp), FUN = sum)
  out$n_years <- length(years)
  out
}

c19 <- sim_counts(sim_years_main)
saveRDS(c19, file.path(dir_data, "sim_counts_2019.rds"))
print(c19[, c("sex", "agegrp", "all", "illdef", "heart_main", "heart_broad",
              "i50", "stroke_main", "stroke_broad")])

# Share of ill-defined and of I50 among heart deaths: report these numbers
c19$pct_illdef  <- round(100 * c19$illdef / c19$all, 1)
c19$pct_i50_in_heart <- round(100 * c19$i50 / c19$heart_main, 1)
write_csv(c19[, c("sex", "agegrp", "all", "pct_illdef", "pct_i50_in_heart")],
          file.path(dir_out, "sim_quality_2019.csv"))

# Sensitivity: pooled 2018-2019 (divide counts by 2 is NOT needed: only
# cause shares are used with the life-table method)
c1819 <- sim_counts(sim_years_alt)
saveRDS(c1819, file.path(dir_data, "sim_counts_2018_2019.rds"))
