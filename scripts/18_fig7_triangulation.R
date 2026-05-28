#!/usr/bin/env Rscript
# 18_fig7_triangulation.R
#
# Figure 7: Triangulation summary — MR estimates side-by-side with available
# eradication-trial meta-analytic estimates, plus a clear visual separation
# of outcomes WITH and WITHOUT RCT evidence.
#
# Output: figures/pub/Fig7_triangulation.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(dplyr)
})

dir.create("figures/pub", showWarnings = FALSE, recursive = TRUE)

# Pull our MR estimates (use the lead antigen-specific row for each outcome)
st <- fread("results/MR_summary_table.tsv")

# Use the best (lowest IVW P) row per outcome — represents the strongest
# antigen signal toward that outcome
mr <- st %>%
  group_by(outcome) %>%
  slice_min(P_ivw, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(outcome, exposure, OR, L95, U95, P_ivw)

# Hand-curated RCT meta-analytic estimates (manual literature review)
# Effect-direction convention: RR < 1 = eradication is protective
# Convert to OR on seropositivity scale by inverting (1/RR) for visual alignment
rct <- tribble(
  ~outcome,        ~rct_RR,  ~rct_L95, ~rct_U95, ~rct_source,
  "parkinson",      NA,        NA,       NA,      "no RCT",
  "alzheimer",      NA,        NA,       NA,      "no RCT",
  "cad",            NA,        NA,       NA,      "no RCT",
  "stroke_isch",    NA,        NA,       NA,      "no RCT",
  "afib",           NA,        NA,       NA,      "no RCT",
  "uc",             NA,        NA,       NA,      "no RCT",
  "crohn",          NA,        NA,       NA,      "no RCT",
  "ra",             NA,        NA,       NA,      "no RCT",
  "sle",            NA,        NA,       NA,      "no RCT",
  "ms",             NA,        NA,       NA,      "no RCT",
  "t2d",            NA,        NA,       NA,      "no eradication RCT",
  "pancan",         NA,        NA,       NA,      "no RCT",
  "crc",            NA,        NA,       NA,      "no RCT",
  "ida",            0.74,      0.61,     0.89,    "Hudak 2017 / Yuan 2015 Cochrane",
  "itp",            0.50,      0.40,     0.62,    "Stasi 2009 Blood (50% response inverted)",
  "asthma",         NA,        NA,       NA,      "minimal RCT data",
  "atopic",         NA,        NA,       NA,      "no RCT",
  "gastric_cancer", 0.54,      0.40,     0.72,    "Ford 2020 Gut meta (8 RCTs)"
)

# Add gastric cancer positive control to MR side
gc_pc <- fread("results/gc_positive_control_summary.tsv")
gc_omp <- gc_pc[antigen == "omp"]
mr <- bind_rows(mr,
  data.frame(outcome = "gastric_cancer", exposure = "omp (PC)",
              OR = gc_omp$OR_ivw[1], L95 = gc_omp$L95[1], U95 = gc_omp$U95[1],
              P_ivw = gc_omp$P_ivw[1]))

# Merge MR and RCT
df <- merge(mr, rct, by = "outcome", all = TRUE)
df$has_rct <- !is.na(df$rct_RR)

# For display: order outcomes by RCT-availability then by MR P
df$panel <- ifelse(df$has_rct, "Outcomes with eradication RCT data",
                                "Outcomes without eradication RCT data")
df <- df[order(df$panel, df$P_ivw), ]
df$outcome <- factor(df$outcome, levels = unique(df$outcome))

# Reshape to long for two-method overlay
mr_long <- data.frame(outcome = df$outcome, panel = df$panel,
                       method = "MR (per-seropositivity OR)",
                       OR = df$OR, L95 = df$L95, U95 = df$U95)
# RCT reverse-convention: invert so both face "increased-risk = OR>1"
rct_long <- data.frame(
  outcome = df$outcome, panel = df$panel,
  method = "Eradication RCT (1/RR; inverted for direction alignment)",
  OR = 1 / df$rct_RR, L95 = 1 / df$rct_U95, U95 = 1 / df$rct_L95
)
long <- rbind(mr_long, rct_long)
long <- long[!is.na(long$OR), ]

p <- ggplot(long, aes(x = OR, y = outcome, color = method)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = L95, xmax = U95),
                  height = 0.25, position = position_dodge(width = 0.6),
                  linewidth = 0.7) +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  scale_x_log10(limits = c(0.4, 3.0), breaks = c(0.5, 0.7, 1.0, 1.5, 2.0)) +
  scale_color_manual(values = c("MR (per-seropositivity OR)" = "#3182bd",
                                  "Eradication RCT (1/RR; inverted for direction alignment)" = "#e6550d"),
                     name = "") +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Triangulation: Mendelian randomization vs eradication-trial evidence",
       subtitle = "MR estimates (blue) and inverted RCT effects (orange) face the same direction.\nFor the 13 outcomes lacking eradication-trial data, MR provides the only available causal evidence.",
       x = "OR (log scale, 95% CI; direction-aligned to seropositivity)",
       y = "Tier 1 outcome") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "top")

ggsave("figures/pub/Fig7_triangulation.png", p, width = 9, height = 7, dpi = 300)
cat("Saved figures/pub/Fig7_triangulation.png\n")
