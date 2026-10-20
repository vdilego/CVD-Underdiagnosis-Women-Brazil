# =============================================================================
# functions_dyd.R  —  indicator maths for the DYD project (base R only)
#
# Notation (as in the report):
#   pi(x)    diagnosed prevalence in age group x (PNS)
#   mc(x)    CVD death rate in age group x (life table x SIM cause share)
#   L(x)     life-table person-years in age group x
#   dyd(x) = pi(x) / mc(x)                      age-specific DYD, in years
#   DYD(a) = sum pi L / sum mc L, ages >= a     = sum w(x) dyd(x),
#            w(x) = mc L / sum mc L             (life-table CVD death weights)
# =============================================================================

# ---- ICD-10 helpers ---------------------------------------------------------

#' Three-character ICD-10 code ("I219" -> "I21"); invalid codes -> NA
icd3 <- function(x) {
  x <- toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
  out <- substr(x, 1, 3)
  out[!grepl("^[A-Z][0-9]{2}$", out)] <- NA_character_
  out
}

#' TRUE if a three-character code falls in any of the ranges of a cause set
in_set <- function(code3, set) {
  hit <- rep(FALSE, length(code3))
  for (r in set) hit <- hit | (!is.na(code3) & code3 >= r[1] & code3 <= r[2])
  hit
}

#' Age in completed years from the raw SIM IDADE field
#' (first digit = unit: 0 min, 1 h, 2 days, 3 months, 4 years, 5 = 100+ years)
sim_age_years <- function(idade) {
  idade <- trimws(as.character(idade))           # always 3 characters in SIM
  unit  <- substr(idade, 1, 1)
  val   <- suppressWarnings(as.numeric(substr(idade, 2, 3)))
  out   <- rep(NA_real_, length(idade))
  out[unit %in% c("0", "1", "2", "3")] <- 0
  out[unit == "4"] <- val[unit == "4"]
  out[unit == "5"] <- 100 + val[unit == "5"]
  out
}

# ---- Life table -------------------------------------------------------------

#' Collapse a single-age life table to the analysis age groups.
#' lt: data.frame with sex, age, mx, Lx (radix 100,000). Returns sex, agegrp,
#' age_lo, L (person-years) and mx_all (deaths / person-years in the group).
lt_to_groups <- function(lt, breaks, labels) {
  lt$agegrp <- cut(lt$age, breaks, right = FALSE, labels = labels)
  lt <- lt[!is.na(lt$agegrp), ]
  lt$dx_L <- lt$mx * lt$Lx                        # life-table deaths
  g <- aggregate(cbind(Lx, dx_L) ~ sex + agegrp, data = lt, FUN = sum)
  data.frame(sex = g$sex, agegrp = g$agegrp,
             age_lo = breaks[as.integer(g$agegrp)],
             L = g$Lx, mx_all = g$dx_L / g$Lx)
}

# ---- Cause-specific death rates under a scenario ---------------------------

#' Death rates for one cause set under one scenario.
#' counts: data.frame with sex, agegrp and death counts per column:
#'   all, illdef, cvd_all, i50, <set>, and <set>_noI50 (the set without I50).
#' lt_g:   output of lt_to_groups().
#' illdef: "prop" = redistribute R00-R99 proportionally over defined causes;
#'         "none" = leave them out of every cause.
#' i50:    "keep" = I50 counted where the set includes it;
#'         "garbage" = I50 deaths redistributed proportionally over the
#'         defined cardiovascular causes (I00-I99 except I50).
#' method: "lt_fraction" = mx_all (IBGE life table) x cause share in SIM;
#'         "sim_pop"     = SIM deaths / population estimate (needs pop column).
cause_rates <- function(counts, lt_g, set, illdef = c("prop", "none"),
                        i50 = c("keep", "garbage"),
                        method = c("lt_fraction", "sim_pop"), pop = NULL) {
  illdef <- match.arg(illdef); i50 <- match.arg(i50); method <- match.arg(method)
  d <- merge(counts, lt_g, by = c("sex", "agegrp"))
  Dc <- d[[set]]
  if (i50 == "garbage") {
    no50    <- d[[paste0(set, "_noI50")]]
    cvd_no50 <- d$cvd_all - d$i50
    Dc <- no50 + d$i50 * no50 / cvd_no50
  }
  defined <- if (illdef == "prop") d$all - d$illdef else d$all
  Dc_adj  <- Dc * d$all / defined                 # deaths after redistribution
  if (method == "lt_fraction") {
    mc <- d$mx_all * Dc_adj / d$all
  } else {
    if (is.null(pop)) stop("method 'sim_pop' needs a population table")
    d <- merge(d, pop, by = c("sex", "agegrp"))
    mc <- Dc_adj / d$pop
  }
  data.frame(sex = d$sex, agegrp = d$agegrp, age_lo = d$age_lo,
             L = d$L, mc = mc, Dc = Dc_adj)
}

# ---- The indicator ----------------------------------------------------------

#' Age-specific dyd(x) with delta-method 95% intervals.
#' d: data.frame with sex, agegrp, pi, se_pi, mc, Dc (and L, age_lo).
#' Variance of log dyd = (se_pi/pi)^2 + 1/Dc  (survey + Poisson deaths).
dyd_age <- function(d) {
  d$dyd    <- d$pi / d$mc
  d$se_log <- sqrt((d$se_pi / d$pi)^2 + 1 / d$Dc)
  d$lo     <- d$dyd * exp(-1.96 * d$se_log)
  d$hi     <- d$dyd * exp( 1.96 * d$se_log)
  d
}

#' Summary DYD(a) for each sex, with 95% intervals.
#' Ignores the covariance of PNS estimates across age groups (small; check
#' with replicate weights if time allows).
dyd_summary <- function(d, a) {
  d <- d[d$age_lo >= a, ]
  out <- lapply(split(d, d$sex), function(s) {
    TD <- sum(s$pi * s$L); C <- sum(s$mc * s$L)
    var_TD <- sum((s$L * s$se_pi)^2)
    var_C  <- sum((s$mc * s$L)^2 / s$Dc)
    se_log <- sqrt(var_TD / TD^2 + var_C / C^2)
    DYD <- TD / C
    data.frame(sex = s$sex[1], age_from = a, DYD = DYD, se_log = se_log,
               lo = DYD * exp(-1.96 * se_log), hi = DYD * exp(1.96 * se_log))
  })
  do.call(rbind, out)
}

#' Life-table CVD death weights w(x) within each sex, ages >= a
dyd_weights <- function(d, a) {
  d <- d[d$age_lo >= a, ]
  d$w <- d$mc * d$L
  d$w <- d$w / ave(d$w, d$sex, FUN = sum)
  d
}

#' Kitagawa split of DYD_F - DYD_M into
#'   rates: sum wbar (dyd_F - dyd_M)        (differences at the same ages)
#'   comp:  sum dydbar (w_F - w_M)          (age distribution of CVD deaths)
#' Also returns DYD for each sex with COMMON weights (wbar).
#' Sex labels must be "Female" and "Male".
dyd_kitagawa <- function(d, a) {
  d <- dyd_weights(dyd_age(d), a)
  f <- d[d$sex == "Female", ]; m <- d[d$sex == "Male", ]
  f <- f[order(f$age_lo), ]; m <- m[order(m$age_lo), ]
  wbar <- (f$w + m$w) / 2; dbar <- (f$dyd + m$dyd) / 2
  gap   <- sum(f$w * f$dyd) - sum(m$w * m$dyd)
  rates <- sum(wbar * (f$dyd - m$dyd))
  comp  <- sum(dbar * (f$w - m$w))
  data.frame(age_from = a,
             DYD_F = sum(f$w * f$dyd), DYD_M = sum(m$w * m$dyd),
             DYD_F_common = sum(wbar * f$dyd), DYD_M_common = sum(wbar * m$dyd),
             gap = gap, due_to_rates = rates, due_to_age_comp = comp,
             share_rates = rates / gap)
}

#' Female/male ratio of dyd(x) by age, with 95% intervals (log scale)
dyd_ratio_fm <- function(d) {
  d <- dyd_age(d)
  f <- d[d$sex == "Female", ]; m <- d[d$sex == "Male", ]
  mm <- merge(f, m, by = c("agegrp", "age_lo"), suffixes = c("_f", "_m"))
  mm$ratio  <- mm$dyd_f / mm$dyd_m
  mm$se_log <- sqrt(mm$se_log_f^2 + mm$se_log_m^2)
  mm$lo <- mm$ratio * exp(-1.96 * mm$se_log)
  mm$hi <- mm$ratio * exp( 1.96 * mm$se_log)
  mm[order(mm$age_lo), c("agegrp", "age_lo", "ratio", "se_log", "lo", "hi")]
}

#' Stroke negative control: log(F/M ratio, heart) - log(F/M ratio, stroke),
#' by age. Negative = cardiac-specific female disadvantage.
#' Treats the two ratios as independent (conservative-ish; PNS estimates are
#' positively correlated within persons, so the true SE is likely smaller).
stroke_did <- function(ratio_heart, ratio_stroke) {
  mm <- merge(ratio_heart, ratio_stroke, by = c("agegrp", "age_lo"),
              suffixes = c("_heart", "_stroke"))
  mm$did <- log(mm$ratio_heart) - log(mm$ratio_stroke)
  mm$se  <- sqrt(mm$se_log_heart^2 + mm$se_log_stroke^2)
  mm$lo  <- mm$did - 1.96 * mm$se; mm$hi <- mm$did + 1.96 * mm$se
  mm[order(mm$age_lo), c("agegrp", "did", "se", "lo", "hi")]
}

#' Tipping point: minimum |difference in survival with disease| (years) that
#' equal diagnosis would require, given the DYD gap and an assumed gap in the
#' competing term K = Comp x Dur_O (default 0).
tipping_point <- function(DYD_F, DYD_M, K_gap = 0) {
  abs((DYD_F - DYD_M) - K_gap)
}

# ---- Risk-factor decomposition ----------------------------------------------

#' Kitagawa split of the sex gap in diagnosed prevalence into
#'   diagnosis given risk: sum pbar(r) [pi_F(r) - pi_M(r)]
#'   risk composition:     sum pibar(r) [p_F(r) - p_M(r)]
#' x: data.frame with sex ("Female"/"Male"), profile, p (share of the sex in
#'    the profile, summing to 1 within sex) and pi (prevalence in the profile).
rf_decomp <- function(x) {
  f <- x[x$sex == "Female", ]; m <- x[x$sex == "Male", ]
  mm <- merge(f, m, by = "profile", suffixes = c("_f", "_m"))
  pbar <- (mm$p_f + mm$p_m) / 2; pibar <- (mm$pi_f + mm$pi_m) / 2
  gap  <- sum(mm$p_f * mm$pi_f) - sum(mm$p_m * mm$pi_m)
  dx   <- sum(pbar * (mm$pi_f - mm$pi_m))
  comp <- sum(pibar * (mm$p_f - mm$p_m))
  data.frame(prev_F = sum(mm$p_f * mm$pi_f), prev_M = sum(mm$p_m * mm$pi_m),
             gap = gap, diagnosis_given_risk = dx, risk_composition = comp,
             expected_F_with_male_rates = sum(mm$p_f * mm$pi_m))
}
