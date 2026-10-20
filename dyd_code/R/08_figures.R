# =============================================================================
# 08_figures.R  —  figures for the GSA slides (16:9, 300 dpi)
# =============================================================================
source("R/00_setup.R")

col_heart <- "#A8432F"; col_stroke <- "#2A7F72"
theme_dyd <- theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.title.position = "plot")

# ---- Fig 1: female/male ratio of dyd(x) by age, heart vs stroke --------------
rat <- read_csv(file.path(dir_out, "dyd_ratio_fm.csv"), show_col_types = FALSE)
g1 <- ggplot(rat, aes(agegrp, ratio, colour = outcome, group = outcome)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.4), size = 0.8) +
  scale_y_log10() +
  scale_colour_manual(values = c(heart = col_heart, stroke = col_stroke),
                      labels = c(heart = "Heart disease", stroke = "Stroke")) +
  labs(x = "Age group", y = "Female / male dyd(x)  (log scale)", colour = NULL,
       title = "Diagnosed years per death, women relative to men",
       subtitle = "Below 1: fewer diagnosed years per death among women. Brazil, 2019") +
  theme_dyd
ggsave(file.path(dir_out, "fig1_ratio_fm.png"), g1, width = 10, height = 5.6, dpi = 300)

# ---- Fig 2: DYD(60) by sex ------------------------------------------------------
sm <- read_csv(file.path(dir_out, "dyd_summary.csv"), show_col_types = FALSE)
g2 <- ggplot(subset(sm, age_from == 60), aes(sex, DYD, colour = sex)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.9) +
  facet_wrap(~outcome, scales = "free_y",
             labeller = as_labeller(c(heart = "Heart disease", stroke = "Stroke"))) +
  scale_colour_manual(values = c(Female = col_heart, Male = "#1B2A3A"), guide = "none") +
  labs(x = NULL, y = "DYD(60), years",
       title = "Diagnosed years per CVD death above age 60") +
  theme_dyd
ggsave(file.path(dir_out, "fig2_dyd60.png"), g2, width = 10, height = 5.6, dpi = 300)

# ---- Fig 3: prevalence ratios, women vs men ---------------------------------------
pr <- read_csv(file.path(dir_out, "rf_prevalence_ratios.csv"), show_col_types = FALSE)
pr <- subset(pr, term == "female" & model != "M5_by_age")
g3 <- ggplot(pr, aes(PR, model, colour = outcome)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(0.5), size = 0.8) +
  scale_colour_manual(values = c(heart = col_heart, stroke = col_stroke),
                      labels = c(heart = "Heart disease", stroke = "Stroke")) +
  labs(x = "Prevalence ratio, women vs men (diagnosed)", y = NULL, colour = NULL,
       title = "Diagnosis given risk: adjusted prevalence ratios, ages 50+") +
  theme_dyd
ggsave(file.path(dir_out, "fig3_prevalence_ratios.png"), g3, width = 10, height = 5.6, dpi = 300)
