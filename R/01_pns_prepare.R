# =============================================================================
# 01_pns_prepare.R  —  PNS 2019: download, survey design, recodes
# Output: data/pns2019_design.rds (survey design, selected residents, all ages,
#         with analysis variables added) and output/pns_cell_counts.csv
# =============================================================================
source("R/00_setup.R")

# ---- 1. Download (first time only; ~ several hundred MB) ---------------------
# selected = TRUE keeps the selected resident (modules P, Q, J ...) with the
# right weight and post-stratification; labels = FALSE keeps numeric codes.
f_raw <- file.path(dir_data, "pns2019_design_raw.rds")
if (!file.exists(f_raw)) {
  des_raw <- PNSIBGE::get_pns(year = 2019, selected = TRUE, labels = FALSE,
                              design = TRUE, savedir = dir_data)
  saveRDS(des_raw, f_raw)
}
des <- readRDS(f_raw)

# ---- 2. Check that every variable exists --------------------------------------
missing <- setdiff(pns_vars, names(des$variables))
if (length(missing) > 0) {
  stop("Not found in PNS 2019 file (fix pns_vars in 00_setup.R): ",
       paste(missing, collapse = ", "))
}

# ---- 3. Recodes ----------------------------------------------------------------
num <- function(code) suppressWarnings(as.numeric(as.character(des$variables[[code]])))
v   <- des$variables
yes <- function(code, yes_values = 1) {           # 1/0, NA kept as NA
  x <- num(code); ifelse(is.na(x), NA_integer_, as.integer(x %in% yes_values))
}

v$age    <- num(pns_vars["age"])
v$sex    <- factor(num(pns_vars["sex"]), levels = c(1, 2), labels = c("Male", "Female"))
v$female <- as.integer(v$sex == "Female")
v$agegrp <- cut(v$age, age_breaks, right = FALSE, labels = age_labels)

# Outcomes
v$heart  <- yes(pns_vars["heart"])
v$stroke <- yes(pns_vars["stroke"])
v$age_dx_heart  <- num(pns_vars["age_dx_heart"])
v$age_dx_stroke <- num(pns_vars["age_dx_stroke"])
v$age_dx_heart[ which(v$age_dx_heart  > v$age | v$heart  != 1)] <- NA  # inconsistent
v$age_dx_stroke[which(v$age_dx_stroke > v$age | v$stroke != 1)] <- NA

# Risk factors (diagnosed; "only during pregnancy" counts as NO)
v$hyp   <- yes(pns_vars["hyp"])
v$diab  <- yes(pns_vars["diab"])
v$chol  <- yes(pns_vars["chol"])
v$smoke <- yes(pns_vars["smoke"], c(1, 2))
v$inactive <- yes(pns_vars["exercise"], 2)        # no exercise in last 3 months
w <- num(pns_vars["weight"]); h <- num(pns_vars["height"])
w[which(w < 25 | w > 250)] <- NA; h[which(h < 120 | h > 220)] <- NA # implausible
v$bmi   <- w / (h / 100)^2
v$obese <- ifelse(is.na(v$bmi), NA_integer_, as.integer(v$bmi >= 30))
v$rf_n  <- rowSums(cbind(v$hyp, v$diab, v$chol, v$smoke, v$obese), na.rm = FALSE)
v$rf_cat <- cut(v$rf_n, c(-Inf, 0, 1, 2, Inf), labels = c("0", "1", "2", "3+"))
v$profile <- interaction(v$hyp, v$diab, v$chol, sep = "")   # 8 profiles, e.g. "101"

# Contact with care
v$visit12 <- yes(pns_vars["visit12"])
v$plan    <- yes(pns_vars["plan"])
v$esf     <- yes(pns_vars["esf"])

# Controls
ed <- num(pns_vars["educ"])
v$educ <- factor(ifelse(ed <= 2, "Low", ifelse(ed <= 5, "Middle", "High")),
                 levels = c("Low", "Middle", "High"))
rc <- num(pns_vars["race"])
v$race <- factor(ifelse(rc == 1, "White", ifelse(rc == 2, "Black",
                 ifelse(rc == 4, "Brown", "Other"))),
                 levels = c("White", "Brown", "Black", "Other"))
v$region <- factor(substr(as.character(num(pns_vars["uf"])), 1, 1),
                   levels = as.character(1:5),
                   labels = c("North", "Northeast", "Southeast", "South", "Center-West"))
v$urban <- factor(num(pns_vars["urban"]), levels = c(1, 2), labels = c("Urban", "Rural"))

des$variables <- v

# ---- 4. Analysis subset: 50+ (subset() keeps the full design for variances) --
des50 <- subset(des, !is.na(agegrp))
saveRDS(des50, file.path(dir_data, "pns2019_design_50plus.rds"))

# ---- 5. Unweighted cell counts — CHECK BEFORE ESTIMATING ---------------------
v50 <- des50$variables
cells <- v50 |>
  group_by(sex, agegrp) |>
  summarise(n = n(),
            heart_pos  = sum(heart  == 1, na.rm = TRUE),
            heart_na   = sum(is.na(heart)),
            stroke_pos = sum(stroke == 1, na.rm = TRUE),
            stroke_na  = sum(is.na(stroke)), .groups = "drop") |>
  mutate(flag_small = heart_pos < min_pos_cases | stroke_pos < min_pos_cases)
print(cells)
write_csv(cells, file.path(dir_out, "pns_cell_counts.csv"))
