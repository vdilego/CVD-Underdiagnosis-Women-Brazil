# =============================================================================
# 04_lifetable.R  —  IBGE 2019 life tables by sex -> analysis age groups
# Input:  data/ibge_lt_2019.csv  (prepared by hand from the IBGE files)
#         columns: sex (1 male, 2 female), age (single years; last row = open
#         group), mx, Lx, lx  (radix 100,000)
# Output: data/lt_groups_2019.rds  (sex, agegrp, age_lo, L, mx_all)
# =============================================================================
source("R/00_setup.R")

lt <- read_csv(file.path(dir_data, "ibge_lt_2019.csv"), show_col_types = FALSE)
lt$sex <- factor(lt$sex, levels = c(1, 2), labels = c("Male", "Female"))

# Checks: radix, monotone lx, open age group
stopifnot(all(abs(lt$lx[lt$age == 0] - 1e5) < 1))
stopifnot(all(tapply(lt$lx, lt$sex, function(x) all(diff(x) <= 0))))
cat("Open age group starts at:", max(lt$age), "\n")

lt_g <- lt_to_groups(as.data.frame(lt), age_breaks, age_labels)
print(lt_g)
saveRDS(lt_g, file.path(dir_data, "lt_groups_2019.rds"))

# ---- Optional: population estimates for the 'sim_pop' scenario ---------------
f_pop <- file.path(dir_data, "pop2019.csv")
if (file.exists(f_pop)) {
  pop <- read_csv(f_pop, show_col_types = FALSE)   # sex, age (single), pop
  pop$sex <- factor(pop$sex, levels = c(1, 2), labels = c("Male", "Female"))
  pop$agegrp <- cut(pop$age, age_breaks, right = FALSE, labels = age_labels)
  pop_g <- aggregate(pop ~ sex + agegrp, data = pop[!is.na(pop$agegrp), ], FUN = sum)
  saveRDS(pop_g, file.path(dir_data, "pop_groups_2019.rds"))
}
