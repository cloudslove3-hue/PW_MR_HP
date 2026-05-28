#!/usr/bin/env Rscript
# 04_summary.R  Aggregate results -> MR_summary_table.tsv + signals
#
# Inputs:  results/all_mr_results.rds
# Output:  results/MR_summary_table.tsv
#          results/signals_FDR_significant.tsv
#          results/MR_secondary_methods.tsv   (Egger, median, mode side-by-side)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("logs",    showWarnings = FALSE, recursive = TRUE)
sink(file.path("logs","04_summary.log"), split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 04_summary.R  ", format(Sys.time()), " ===\n")

all_results <- readRDS("results/all_mr_results.rds")

safe_pval_pleio <- function(p) {
  if (is.null(p) || inherits(p, "try-error") || !is.data.frame(p) || nrow(p) == 0) return(NA_real_)
  p$pval[1]
}
safe_q_pval <- function(h) {
  if (is.null(h) || inherits(h, "try-error") || !is.data.frame(h) || nrow(h) == 0) return(NA_real_)
  idx <- which(h$method == "Inverse variance weighted")[1]
  if (is.na(idx)) return(NA_real_)
  h$Q_pval[idx]
}

rows <- lapply(all_results, function(x) {
  if (is.null(x)) return(NULL)
  res <- x$mr_res
  if (is.null(res) || !nrow(res)) return(NULL)
  pick <- function(m) {
    r <- res[res$method == m, ]
    if (!nrow(r)) return(c(b = NA, se = NA, p = NA))
    c(b = r$b[1], se = r$se[1], p = r$pval[1])
  }
  ivw <- pick("Inverse variance weighted")
  egg <- pick("MR Egger")
  wmd <- pick("Weighted median")
  wmo <- pick("Weighted mode")
  data.frame(
    exposure = x$exposure, outcome = x$outcome, n_snp = x$n_snp,
    # IVW primary
    OR    = exp(ivw["b"]),
    L95   = exp(ivw["b"] - 1.96 * ivw["se"]),
    U95   = exp(ivw["b"] + 1.96 * ivw["se"]),
    P_ivw = ivw["p"],
    # secondary
    OR_egger     = exp(egg["b"]),  P_egger     = egg["p"],
    OR_wmedian   = exp(wmd["b"]),  P_wmedian   = wmd["p"],
    OR_wmode     = exp(wmo["b"]),  P_wmode     = wmo["p"],
    # sensitivity
    P_egger_intercept = safe_pval_pleio(x$pleiotropy),
    Q_pval            = safe_q_pval(x$heterogeneity),
    # MR-PRESSO (may return "<0.001" as character — coerce)
    presso_global_p = {
      v <- if (!is.null(x$presso) && !inherits(x$presso, "try-error")) {
        tryCatch(x$presso$`MR-PRESSO results`$`Global Test`$Pvalue,
                 error = function(e) NA)
      } else NA
      if (is.character(v)) {
        v <- suppressWarnings(as.numeric(sub("^<", "", v)))
      }
      as.numeric(v[1])
    },
    row.names = NULL
  )
})
summary_tbl <- bind_rows(rows)

summary_tbl$P_bonf <- p.adjust(summary_tbl$P_ivw, method = "bonferroni")
summary_tbl$P_fdr  <- p.adjust(summary_tbl$P_ivw, method = "fdr")
summary_tbl <- summary_tbl[order(summary_tbl$P_ivw), ]

fwrite(summary_tbl, "results/MR_summary_table.tsv", sep = "\t")
cat("Wrote results/MR_summary_table.tsv  (", nrow(summary_tbl), "rows )\n")

sig <- subset(summary_tbl, !is.na(P_fdr) & P_fdr < 0.05)
fwrite(sig, "results/signals_FDR_significant.tsv", sep = "\t")
cat("Wrote results/signals_FDR_significant.tsv (", nrow(sig), "rows )\n")

cat("\nTop 10 by IVW P:\n")
print(head(summary_tbl[, c("exposure","outcome","n_snp","OR","L95","U95",
                            "P_ivw","P_fdr","P_egger_intercept","Q_pval")], 10))
cat("\nDone. Proceed to scripts/05_figures.R\n")
