#!/usr/bin/env Rscript
# 07_loo_outliers.R  Leave-one-out + MR-PRESSO outlier handling
#
# For every pair flagged by QC (Cochran Q P<0.05 OR MR-PRESSO global P<0.05),
# run leave-one-out IVW and report PRESSO outlier-corrected estimate.
#
# Input:   results/all_mr_results.rds, results/MR_summary_table.tsv
# Output:  results/leave_one_out.tsv
#          results/presso_corrected_estimates.tsv
#          logs/07_loo_outliers.log

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("logs",    showWarnings = FALSE, recursive = TRUE)
sink("logs/07_loo_outliers.log", split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 07_loo_outliers.R  ", format(Sys.time()), " ===\n")

all_results <- readRDS("results/all_mr_results.rds")
st          <- fread("results/MR_summary_table.tsv")

flagged <- st$exposure[!is.na(st$Q_pval) & st$Q_pval < 0.05 |
                      !is.na(st$presso_global_p) & st$presso_global_p < 0.05]
flagged_keys <- paste(st$exposure, st$outcome, sep = "__")[
  (!is.na(st$Q_pval) & st$Q_pval < 0.05) |
  (!is.na(st$presso_global_p) & st$presso_global_p < 0.05)]
flagged_keys <- unique(flagged_keys)
cat("Flagged pairs:", length(flagged_keys), "\n")

# ---- Leave-one-out per flagged pair -----------------------------------------
loo_rows <- list()
for (key in flagged_keys) {
  x <- all_results[[key]]
  if (is.null(x) || is.null(x$dat) || nrow(x$dat) < 4) next
  res <- tryCatch(mr_leaveoneout(x$dat, method = TwoSampleMR::mr_ivw),
                  error = function(e) NULL)
  if (is.null(res) || !nrow(res)) next
  res$exposure_label <- x$exposure
  res$outcome_label  <- x$outcome
  loo_rows[[key]] <- res
}
loo_all <- bind_rows(loo_rows)
if (nrow(loo_all)) {
  loo_all$OR  <- exp(loo_all$b)
  loo_all$L95 <- exp(loo_all$b - 1.96 * loo_all$se)
  loo_all$U95 <- exp(loo_all$b + 1.96 * loo_all$se)
  cols <- c("exposure_label","outcome_label","SNP","samplesize",
            "b","se","p","OR","L95","U95")
  cols <- intersect(cols, names(loo_all))
  fwrite(loo_all[, cols], "results/leave_one_out.tsv", sep = "\t")
  cat("Wrote results/leave_one_out.tsv  (", nrow(loo_all), "rows )\n")
} else {
  cat("No LOO rows produced\n")
}

# ---- Identify dominant SNPs ------------------------------------------------
# Per pair: which SNP, when removed, swings P most? (P_loo - P_all_in)
dominant <- list()
for (key in flagged_keys) {
  x <- all_results[[key]]
  if (is.null(x) || is.null(x$mr_res)) next
  pall <- x$mr_res$pval[x$mr_res$method == "Inverse variance weighted"][1]
  if (is.na(pall)) next
  sub <- loo_all[loo_all$exposure_label == x$exposure &
                  loo_all$outcome_label == x$outcome &
                  loo_all$SNP != "All", ]
  if (!nrow(sub)) next
  sub$delta_p <- sub$p - pall
  worst <- sub[which.max(abs(sub$delta_p)), ]
  dominant[[key]] <- data.frame(
    exposure = x$exposure, outcome = x$outcome,
    p_full = pall,
    most_influential_SNP = worst$SNP,
    p_without = worst$p,
    p_changed_by = round(worst$delta_p, 4),
    flips_significance = (pall < 0.05) != (worst$p < 0.05)
  )
}
dom <- bind_rows(dominant)
if (nrow(dom)) {
  fwrite(dom, "results/loo_dominant_snps.tsv", sep = "\t")
  cat("Wrote results/loo_dominant_snps.tsv\n")
  cat("Pairs whose significance flips on removing a single SNP:",
      sum(dom$flips_significance), "\n")
  if (any(dom$flips_significance)) {
    cat("\n=== Flipping pairs ===\n")
    print(dom[dom$flips_significance, ])
  }
}

# ---- PRESSO outlier-corrected estimates -------------------------------------
presso_rows <- list()
for (key in names(all_results)) {
  x <- all_results[[key]]
  if (is.null(x) || is.null(x$presso) || inherits(x$presso, "try-error")) next
  pr <- tryCatch(x$presso, error = function(e) NULL)
  if (is.null(pr)) next
  main <- pr$`Main MR results`
  if (is.null(main) || !is.data.frame(main)) next
  coerce_num <- function(v) {
    if (is.null(v) || length(v) == 0) return(NA_real_)
    if (is.character(v)) v <- suppressWarnings(as.numeric(sub("^<", "", v)))
    as.numeric(v[1])
  }
  for (i in seq_len(nrow(main))) {
    presso_rows[[length(presso_rows) + 1]] <- data.frame(
      exposure = x$exposure, outcome = x$outcome,
      mr_type  = as.character(main$`MR Analysis`[i]),
      causal_beta = coerce_num(main$`Causal Estimate`[i]),
      sd          = coerce_num(main$Sd[i]),
      t_stat      = coerce_num(main$`T-stat`[i]),
      pval        = coerce_num(main$`P-value`[i]),
      global_p    = coerce_num(pr$`MR-PRESSO results`$`Global Test`$Pvalue),
      n_outliers  = if (!is.null(pr$`MR-PRESSO results`$`Distortion Test`)) {
        length(pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
      } else 0L
    )
  }
}
presso_tbl <- bind_rows(presso_rows)
if (nrow(presso_tbl)) {
  presso_tbl$OR  <- exp(presso_tbl$causal_beta)
  presso_tbl$L95 <- exp(presso_tbl$causal_beta - 1.96 * presso_tbl$sd)
  presso_tbl$U95 <- exp(presso_tbl$causal_beta + 1.96 * presso_tbl$sd)
  # Wide: compare raw vs outlier-corrected
  raw <- presso_tbl[presso_tbl$mr_type == "Raw", ]
  oc  <- presso_tbl[presso_tbl$mr_type == "Outlier-corrected", ]
  oc$pval <- suppressWarnings(as.numeric(sub("^<", "", as.character(oc$pval))))
  wide <- merge(
    raw[, c("exposure","outcome","causal_beta","sd","pval","OR","L95","U95","global_p")],
    oc[,  c("exposure","outcome","causal_beta","sd","pval","OR","L95","U95","n_outliers")],
    by = c("exposure","outcome"), suffixes = c("_raw","_corrected"), all = TRUE
  )
  fwrite(wide, "results/presso_corrected_estimates.tsv", sep = "\t")
  cat("Wrote results/presso_corrected_estimates.tsv  (",
      nrow(wide), "pairs )\n")
  # Pairs where outlier removal changes significance
  flip <- subset(wide, !is.na(pval_raw) & !is.na(pval_corrected) &
                       (pval_raw < 0.05) != (pval_corrected < 0.05))
  cat("Pairs where outlier removal flips P<0.05:", nrow(flip), "\n")
  if (nrow(flip)) print(flip[, c("exposure","outcome","pval_raw","pval_corrected",
                                  "OR_raw","OR_corrected")])
}

cat("\nDone.\n")
