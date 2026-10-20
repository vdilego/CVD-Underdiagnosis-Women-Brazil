# =============================================================================
# tests/test_functions.R  —  checks the indicator maths on synthetic data.
# Base R only. Run from the project folder:  Rscript tests/test_functions.R
# Every check must print "ok". Run it again after any change to functions_dyd.R
# =============================================================================
source("R/functions_dyd.R")
ok <- function(cond, msg) { if (!isTRUE(cond)) stop("FAILED: ", msg); cat("ok  ", msg, "\n") }
age_breaks <- c(50, 60, 70, 80, Inf); age_labels <- c("50-59", "60-69", "70-79", "80+")

# ---- ICD and age helpers ------------------------------------------------------
codes <- c("I219", "I50.0", "I109", "I110", "I64", "I690", "R99", "C34", "x", NA)
c3 <- icd3(codes)
ok(identical(c3, c("I21", "I50", "I10", "I11", "I64", "I69", "R99", "C34", NA, NA)), "icd3")
heart_main <- list(c("I20", "I25"), c("I50", "I50"))
heart_broad <- list(c("I05", "I09"), c("I11", "I11"), c("I13", "I13"), c("I20", "I51"))
ok(identical(in_set(c3, heart_main), c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)), "in_set main")
ok(identical(in_set(c3, heart_broad), c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)), "in_set broad (I10 out, I11 in)")
ok(identical(sim_age_years(c("465", "501", "310", "220", "400")), c(65, 101, 0, 0, 0)), "SIM age decoding")

# ---- Synthetic life table (Gompertz) -------------------------------------------
mk_lt <- function(sex, a, b) {
  age <- 0:100; mx <- pmin(a * exp(b * age), 1.5)
  qx <- 1 - exp(-mx); lx <- 1e5 * cumprod(c(1, 1 - qx[-length(qx)]))
  Lx <- c(lx[-1] + (lx[-length(lx)] - lx[-1]) / 2, lx[length(lx)] / mx[length(mx)])
  data.frame(sex = sex, age = age, mx = mx, Lx = Lx, lx = lx)
}
lt <- rbind(mk_lt("Male", 2e-4, 0.085), mk_lt("Female", 1e-4, 0.088))
lt_g <- lt_to_groups(lt, age_breaks, age_labels)
ok(nrow(lt_g) == 8 && all(lt_g$L > 0), "life table grouped to 4 age groups x 2 sexes")
ok(all(diff(lt_g$mx_all[lt_g$sex == "Male"][order(lt_g$age_lo[lt_g$sex == "Male"])]) > 0),
   "grouped death rates increase with age")

# ---- Synthetic SIM counts --------------------------------------------------------
counts <- expand.grid(sex = c("Male", "Female"), agegrp = factor(age_labels, levels = age_labels))
counts$all <- c(60e3, 40e3, 90e3, 70e3, 120e3, 110e3, 130e3, 200e3)
counts$illdef <- round(0.06 * counts$all)
counts$cvd_all <- round(0.30 * counts$all)
counts$ihd <- round(0.12 * counts$all)
counts$i50 <- round(0.03 * counts$all) + ifelse(counts$sex == "Female", 2000, 0)
counts$heart_main <- counts$ihd + counts$i50
counts$heart_main_noI50 <- counts$ihd
counts$stroke_main <- round(0.09 * counts$all); counts$stroke_main_noI50 <- counts$stroke_main

r_none <- cause_rates(counts, lt_g, "heart_main", illdef = "none")
chk <- merge(r_none, merge(counts, lt_g), by = c("sex", "agegrp"))
ok(isTRUE(all.equal(chk$mc, chk$mx_all * chk$heart_main / chk$all)), "rate = mx_all x cause share")
r_prop <- cause_rates(counts, lt_g, "heart_main", illdef = "prop")
ok(all(r_prop$mc > r_none$mc[match(paste(r_prop$sex, r_prop$agegrp), paste(r_none$sex, r_none$agegrp))]),
   "ill-defined redistribution raises cause rates")
r_g <- cause_rates(counts, lt_g, "heart_main", illdef = "none", i50 = "garbage")
ok(all(r_g$mc < r_none$mc[match(paste(r_g$sex, r_g$agegrp), paste(r_none$sex, r_none$agegrp))]),
   "I50 garbage scenario lowers heart rates (part of I50 goes to other CVD)")
# total CVD deaths are preserved by the I50 redistribution
cvd_no50 <- counts$cvd_all - counts$i50
ok(isTRUE(all.equal(sum(counts$i50 * counts$ihd / cvd_no50 + counts$i50 * (cvd_no50 - counts$ihd) / cvd_no50),
                    sum(counts$i50))), "I50 redistribution conserves deaths")

# ---- The indicator ----------------------------------------------------------------
build <- function(r, pi_fun) {
  d <- r; d$pi <- pi_fun(d); d$se_pi <- 0.1 * d$pi
  d$sex <- as.character(d$sex); d$agegrp <- as.character(d$agegrp); d
}
# (a) prevalence proportional to the death rate -> dyd constant -> DYD equals it
d_const <- build(r_prop, function(d) 7 * d$mc)
ok(all(abs(dyd_age(d_const)$dyd - 7) < 1e-10), "dyd(x) = pi/mc")
ok(all(abs(dyd_summary(d_const, 60)$DYD - 7) < 1e-10), "DYD(a) of a constant dyd = that constant")

# (b) DYD(a) equals the w-weighted average of dyd(x)
d_var <- build(r_prop, function(d) pmin(0.05 + 0.002 * (d$age_lo - 50) - 0.01 * (d$sex == "Female"), 0.9))
w <- dyd_weights(dyd_age(d_var), 50)
avg <- tapply(w$w * w$dyd, w$sex, sum)
sm <- dyd_summary(d_var, 50)
ok(isTRUE(all.equal(as.numeric(avg[sm$sex]), sm$DYD)), "DYD(a) = sum w(x) dyd(x)")
ok(all(sm$lo < sm$DYD & sm$DYD < sm$hi), "intervals contain the estimate")

# (c) Kitagawa: the two parts add up to the gap
k <- dyd_kitagawa(d_var, 50)
ok(isTRUE(all.equal(k$gap, k$due_to_rates + k$due_to_age_comp)), "Kitagawa parts sum to the gap")
ok(isTRUE(all.equal(k$DYD_F, sm$DYD[sm$sex == "Female"])), "Kitagawa DYD_F matches dyd_summary")

# (d) ratio and stroke difference-in-differences
rt <- dyd_ratio_fm(d_var)
ok(nrow(rt) == 4 && all(rt$lo < rt$ratio & rt$ratio < rt$hi), "F/M ratio by age with CIs")
did0 <- stroke_did(rt, rt)
ok(all(abs(did0$did) < 1e-12), "stroke DiD is zero when heart and stroke ratios are equal")
ok(tipping_point(5, 7) == 2 && tipping_point(5, 7, K_gap = -1) == 1, "tipping point")

# ---- Risk-factor decomposition ------------------------------------------------------
x <- data.frame(sex = rep(c("Female", "Male"), each = 3), profile = rep(c("0", "1", "2"), 2),
                p  = c(0.2, 0.4, 0.4, 0.4, 0.4, 0.2),
                pi = c(0.03, 0.08, 0.15, 0.04, 0.10, 0.20))
dec <- rf_decomp(x)
ok(isTRUE(all.equal(dec$gap, dec$diagnosis_given_risk + dec$risk_composition)), "risk-factor decomposition adds up")
ok(dec$risk_composition > 0 && dec$diagnosis_given_risk < 0,
   "synthetic paradox detected: women's risk predicts more, diagnosis-given-risk lower")

cat("\nAll checks passed.\n")
