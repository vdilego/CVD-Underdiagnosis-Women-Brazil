# End-to-end check of 05_dyd.R logic on synthetic inputs (base R only)
source("R/functions_dyd.R"); ok <- function(...) invisible(TRUE)
age_breaks <- c(50, 60, 70, 80, Inf); age_labels <- c("50-59", "60-69", "70-79", "80+")
outcome_sets <- list(heart = c(main = "heart_main"), stroke = c(main = "stroke_main"))
lines <- readLines("tests/test_functions.R")
eval(parse(text = lines[grep("^mk_lt", lines):grep("^counts\\$stroke_main", lines)]))
lt_g <- lt_to_groups(lt, age_breaks, age_labels)
# prevalence as read from CSV: character keys
prev <- expand.grid(outcome = c("heart", "stroke"), sex = c("Male", "Female"),
                    agegrp = age_labels, stringsAsFactors = FALSE)
prev$pi <- with(prev, ifelse(outcome == "heart", 0.06, 0.03) *
                  (1 + 0.4 * (match(agegrp, age_labels) - 1)) * ifelse(sex == "Female", 0.8, 1))
prev$se_pi <- 0.08 * prev$pi
src <- readLines("R/05_dyd.R")
eval(parse(text = src[grep("^build_input", src):(grep("^h <- run_main", src) - 1)]))
h <- run_main("heart"); s <- run_main("stroke")
stopifnot(nrow(h$age) == 8, nrow(h$ratio) == 4, all(is.finite(h$age$dyd)))
print(h$summary); print(h$kit); print(stroke_did(h$ratio[, -1], s$ratio[, -1]))
cat("Pipeline OK\n")
