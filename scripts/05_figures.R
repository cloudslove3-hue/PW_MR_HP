#!/usr/bin/env Rscript
# 05_figures.R  Forest + volcano + session log
#
# Input:   results/MR_summary_table.tsv
# Output:  figures/forest/<outcome>.png  (one per outcome with >=2 antigens)
#          figures/volcano_overall.png
#          logs/sessionInfo.txt

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

dir.create("figures/forest", showWarnings = FALSE, recursive = TRUE)
dir.create("logs",           showWarnings = FALSE, recursive = TRUE)
sink(file.path("logs","05_figures.log"), split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 05_figures.R  ", format(Sys.time()), " ===\n")

st <- fread("results/MR_summary_table.tsv")

# ---- Forest per outcome -----------------------------------------------------
plot_forest <- function(df, title) {
  ggplot(df, aes(x = OR, y = exposure)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = L95, xmax = U95), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_x_log10() +
    labs(title = title, x = "OR (log scale)", y = "H. pylori antigen") +
    theme_minimal(base_size = 12)
}

per_outcome <- split(st, st$outcome)
n_plots <- 0
for (nm in names(per_outcome)) {
  sub <- per_outcome[[nm]]
  if (nrow(sub) < 2) next
  fn <- sprintf("figures/forest/%s.png", gsub("[^A-Za-z0-9]+", "_", nm))
  ggsave(fn, plot_forest(sub, nm), width = 7, height = 4, dpi = 300)
  n_plots <- n_plots + 1
}
cat("Forest plots written:", n_plots, "\n")

# ---- Volcano (overall) ------------------------------------------------------
v <- ggplot(st, aes(x = log(OR), y = -log10(P_ivw),
                     color = !is.na(P_fdr) & P_fdr < 0.05)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  facet_wrap(~ exposure) +
  scale_color_manual(values = c("grey60","red"),
                     name = "FDR<0.05") +
  labs(title = "H. pylori antigens vs Tier 1 outcomes",
       x = "log(OR)", y = "-log10(P)") +
  theme_minimal(base_size = 12)
ggsave("figures/volcano_overall.png", v, width = 10, height = 7, dpi = 300)
cat("Volcano written: figures/volcano_overall.png\n")

# ---- Session info -----------------------------------------------------------
con <- file("logs/sessionInfo.txt", open = "w")
writeLines(capture.output(sessionInfo()), con); close(con)
writeLines(format(Sys.time()), "logs/last_run_timestamp.txt")

cat("Done.\n")
