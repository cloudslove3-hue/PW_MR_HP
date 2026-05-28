#!/usr/bin/env Rscript
# 11_publication_figures.R  Publication-grade figures
#
# Outputs (figures/pub/):
#   Fig1_heatmap.png       7 antigens x 17 outcomes OR heatmap with significance stars
#   Fig2_forest_crohn.png  Crohn focused forest (the lone borderline signal)
#   Fig3_mvmr_compare.png  Univariate vs MVMR direct effect (top candidate pairs)
#   Fig4_coloc.png         OMP loci coloc PP barplot (HLA vs non-HLA)
#   Fig5_volcano.png       improved volcano with annotations
#
# Designed for Nat Commun-style figures: white background, sans-serif, legible

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(scales)
})

dir.create("figures/pub", showWarnings = FALSE, recursive = TRUE)
sink("logs/11_publication_figures.log", split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 11_publication_figures.R  ", format(Sys.time()), " ===\n")

# ---- Common theme -----------------------------------------------------------
theme_pub <- function(base = 11) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      strip.background   = element_rect(fill = "grey95", color = NA),
      strip.text         = element_text(face = "bold"),
      plot.title         = element_text(face = "bold", size = base + 2),
      legend.position    = "top"
    )
}

st <- fread("results/MR_summary_table.tsv")

# ---- Fig 1: Heatmap of OR (7 antigens × 17 outcomes) -----------------------
ant_order <- c("hpylori_general","cagA","catalase","groEL","omp","ureA","vacA")
out_order <- unique(st$outcome)
hm <- st[, .(exposure, outcome, OR, P_ivw, P_fdr)]
hm$exposure <- factor(hm$exposure, levels = ant_order)
hm$outcome  <- factor(hm$outcome,  levels = rev(out_order))
hm$log_or   <- log(hm$OR)
hm$sig <- with(hm, ifelse(!is.na(P_fdr) & P_fdr < 0.05, "***",
                  ifelse(!is.na(P_ivw) & P_ivw < 0.05, "*", "")))

p1 <- ggplot(hm, aes(x = exposure, y = outcome, fill = log_or)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig), size = 4) +
  scale_fill_gradient2(low = "#3182bd", mid = "white", high = "#de2d26",
                        midpoint = 0, name = "log(OR)",
                        limits = c(-0.3, 0.3), oob = squish) +
  labs(title = "H. pylori antigens vs Tier 1 outcomes (IVW)",
       subtitle = "* P<0.05 (raw)  *** FDR<0.05",
       x = "H. pylori antigen", y = "Outcome") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("figures/pub/Fig1_heatmap.png", p1, width = 7.5, height = 6, dpi = 300)
cat("Saved Fig1_heatmap.png\n")

# ---- Fig 2: Crohn forest (the lone borderline signal) ----------------------
sub <- st[outcome == "crohn"][order(OR)]
sub$exposure <- factor(sub$exposure, levels = sub$exposure)
p2 <- ggplot(sub, aes(x = OR, y = exposure)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = L95, xmax = U95), height = 0.2, linewidth = 0.7) +
  geom_point(size = 3.2, color = "#e6550d") +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)\nP=%.3f", OR, L95, U95, P_ivw)),
            hjust = -0.1, size = 3, family = "sans") +
  scale_x_log10(limits = c(0.4, 2.5), breaks = c(0.5, 0.7, 1.0, 1.5, 2.0)) +
  labs(title = "Crohn's disease ~ H. pylori antigens (IVW)",
       subtitle = "Borderline OMP signal subsequently falsified by colocalization (PP.H4<0.05)",
       x = "OR (log scale, 95% CI)", y = "Antigen") +
  theme_pub()
ggsave("figures/pub/Fig2_forest_crohn.png", p2, width = 8, height = 4.5, dpi = 300)
cat("Saved Fig2_forest_crohn.png\n")

# ---- Fig 3: MVMR vs univariate comparison ---------------------------------
if (file.exists("results/mvmr_results.tsv")) {
  mv <- fread("results/mvmr_results.tsv")
  ant_only <- mv[is_antigen == TRUE]
  ant_only$OR_mvmr <- exp(ant_only$b)
  ant_only$L95_mvmr <- exp(ant_only$b - 1.96 * ant_only$se)
  ant_only$U95_mvmr <- exp(ant_only$b + 1.96 * ant_only$se)
  comp <- merge(ant_only[, .(exposure = mvmr_exposure_label,
                              outcome  = mvmr_outcome_label,
                              OR_mvmr, L95_mvmr, U95_mvmr, P_mvmr = pval)],
                st[, .(exposure, outcome, OR_uni = OR,
                        L95_uni = L95, U95_uni = U95, P_uni = P_ivw)],
                by = c("exposure","outcome"))
  comp$pair <- paste(comp$exposure, comp$outcome, sep = " → ")
  comp <- comp[order(P_uni)]
  comp$pair <- factor(comp$pair, levels = rev(comp$pair))

  long <- rbind(
    data.frame(pair = comp$pair, OR = comp$OR_uni,
               L95 = comp$L95_uni, U95 = comp$U95_uni,
               P = comp$P_uni, model = "Univariate IVW"),
    data.frame(pair = comp$pair, OR = comp$OR_mvmr,
               L95 = comp$L95_mvmr, U95 = comp$U95_mvmr,
               P = comp$P_mvmr, model = "MVMR (BMI+smoking+edu)")
  )

  p3 <- ggplot(long, aes(x = OR, y = pair, color = model)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = L95, xmax = U95),
                    height = 0.25, position = position_dodge(width = 0.6),
                    linewidth = 0.7) +
    geom_point(size = 2.8, position = position_dodge(width = 0.6)) +
    scale_x_log10() +
    scale_color_manual(values = c("Univariate IVW" = "#e6550d",
                                    "MVMR (BMI+smoking+edu)" = "#3182bd"),
                       name = "") +
    labs(title = "Univariate MR vs MVMR-adjusted direct effects",
         subtitle = "All antigen direct effects null after adjustment",
         x = "OR (log scale, 95% CI)", y = "") +
    theme_pub()
  ggsave("figures/pub/Fig3_mvmr_compare.png", p3, width = 8.5, height = 5, dpi = 300)
  cat("Saved Fig3_mvmr_compare.png\n")
}

# ---- Fig 4: Coloc PP barplot ----------------------------------------------
if (file.exists("results/coloc_omp_crohn.tsv")) {
  co <- fread("results/coloc_omp_crohn.tsv")
  co$region <- with(co, ifelse(chr == 6, paste0(region_lead, " (HLA)"),
                                paste0(region_lead, " (chr", chr, ")")))
  co$region <- factor(co$region, levels = co$region)
  ppm <- melt(co[, .(region, PP_H0, PP_H1, PP_H2, PP_H3, PP_H4)],
              id.vars = "region", variable.name = "Hypothesis",
              value.name = "PP")
  ppm$Hypothesis <- factor(ppm$Hypothesis,
    levels = c("PP_H0","PP_H1","PP_H2","PP_H3","PP_H4"),
    labels = c("H0: no causal var",
                "H1: OMP only",
                "H2: Crohn only",
                "H3: distinct vars",
                "H4: shared var"))

  p4 <- ggplot(ppm, aes(x = region, y = PP, fill = Hypothesis)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("#d9d9d9","#3182bd","#9ecae1",
                                  "#fd8d3c","#e6550d")) +
    labs(title = "Bayesian colocalization at OMP instrument loci",
         subtitle = "PP.H4 (shared causal variant) <0.05 at every locus",
         x = "Locus (lead SNP)", y = "Posterior probability") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
  ggsave("figures/pub/Fig4_coloc.png", p4, width = 9, height = 5, dpi = 300)
  cat("Saved Fig4_coloc.png\n")
}

# ---- Fig 5: Improved volcano ------------------------------------------------
vt <- st
vt$sig_label <- with(vt,
  ifelse(!is.na(P_fdr) & P_fdr < 0.05, "FDR<0.05",
  ifelse(!is.na(P_ivw) & P_ivw < 0.05, "P<0.05 raw", "ns")))
vt$pair <- with(vt, paste(exposure, outcome, sep = " → "))
top_label <- vt[order(P_ivw)][1:5, ]

p5 <- ggplot(vt, aes(x = log(OR), y = -log10(P_ivw), color = sig_label)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.7, size = 1.7) +
  facet_wrap(~ exposure, nrow = 2) +
  scale_color_manual(values = c("FDR<0.05" = "#e6550d",
                                "P<0.05 raw" = "#fd8d3c",
                                "ns" = "grey60"), name = "") +
  labs(title = "Phenome-wide MR: H. pylori 7 antigens × 17 outcomes",
       subtitle = "No pair survives FDR correction",
       x = "log(OR)", y = expression(-log[10](P))) +
  theme_pub()
ggsave("figures/pub/Fig5_volcano.png", p5, width = 10, height = 6, dpi = 300)
cat("Saved Fig5_volcano.png\n")

cat("\nAll publication figures written to figures/pub/\n")
