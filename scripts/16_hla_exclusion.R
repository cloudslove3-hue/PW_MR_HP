#!/usr/bin/env Rscript
# 16_hla_exclusion.R
#
# Round 2 / TASK 2 (lite): HLA-exclusion sensitivity analysis.
# Re-run Tier 1 IVW MR after dropping all SNPs in chr6:25,000,000-35,000,000
# from every antigen instrument set. The purpose is to test whether univariate
# signals (especially OMP-driven ones) persist outside the HLA pleiotropy
# region, complementing the coloc.abf finding (PP.H4<0.05 at every OMP locus).
#
# Outputs:
#   results/MR_summary_noHLA.tsv
#   results/qc_hla_exclusion.tsv
#   logs/16_hla_exclusion.log

set.seed(20260521)

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("logs",    showWarnings = FALSE, recursive = TRUE)
sink("logs/16_hla_exclusion.log", split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 16_hla_exclusion.R  ", format(Sys.time()), " ===\n")

# ---- Load primary data -----------------------------------------------------
exposure_list     <- readRDS("data/exposure/exposure_list.rds")
outcome_data_list <- readRDS("data/outcome/outcome_data_list.rds")

HLA_CHR <- "6"
HLA_LO  <- 25e6
HLA_HI  <- 35e6

# ---- Filter instruments to drop HLA SNPs -----------------------------------
exposure_list_noHLA <- list()
qc_rows <- list()
for (lab in names(exposure_list)) {
  d <- exposure_list[[lab]]
  if (is.null(d)) next
  in_hla <- d$chr.exposure == HLA_CHR &
            d$pos.exposure >= HLA_LO &
            d$pos.exposure <= HLA_HI
  d_kept <- d[!in_hla, ]
  qc_rows[[lab]] <- data.frame(
    antigen = lab, n_total = nrow(d),
    n_HLA_dropped = sum(in_hla),
    n_kept = nrow(d_kept),
    HLA_SNPs = paste(d$SNP[in_hla], collapse = ",")
  )
  if (nrow(d_kept) >= 3) {
    exposure_list_noHLA[[lab]] <- d_kept
  } else {
    cat(lab, ": only ", nrow(d_kept), "SNPs after HLA exclusion, skipping\n")
  }
}
qc <- bind_rows(qc_rows)
fwrite(qc, "results/qc_hla_exclusion.tsv", sep = "\t")
cat("\nHLA-exclusion QC:\n"); print(qc[, c("antigen","n_total","n_HLA_dropped","n_kept")])

# ---- Re-run MR for all antigen × outcome pairs -----------------------------
run_pair <- function(exp_df, out_df, exp_lab, out_lab) {
  if (is.null(out_df) || nrow(out_df) == 0) return(NULL)
  dat <- tryCatch(harmonise_data(exp_df, out_df, action = 2),
                  error = function(e) NULL)
  if (is.null(dat)) return(NULL)
  dat <- subset(dat, mr_keep == TRUE)
  if (nrow(dat) < 3) return(NULL)
  res <- tryCatch(mr(dat, method_list = c("mr_ivw","mr_weighted_median",
                                            "mr_egger_regression")),
                  error = function(e) NULL)
  if (is.null(res) || !nrow(res)) return(NULL)
  pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  het   <- tryCatch(mr_heterogeneity(dat),    error = function(e) NULL)
  list(exposure = exp_lab, outcome = out_lab, dat = dat, mr_res = res,
       pleiotropy = pleio, heterogeneity = het, n_snp = nrow(dat))
}

rows <- list()
for (lab in names(exposure_list_noHLA)) {
  exp_df <- exposure_list_noHLA[[lab]]
  for (out_lab in names(outcome_data_list)) {
    out_df <- outcome_data_list[[out_lab]]
    x <- run_pair(exp_df, out_df, lab, out_lab)
    if (is.null(x)) next
    ivw <- subset(x$mr_res, method == "Inverse variance weighted")
    wmd <- subset(x$mr_res, method == "Weighted median")
    egg <- subset(x$mr_res, method == "MR Egger")
    rows[[length(rows) + 1]] <- data.frame(
      exposure = lab, outcome = out_lab, n_snp = x$n_snp,
      OR    = if (nrow(ivw)) exp(ivw$b[1]) else NA,
      L95   = if (nrow(ivw)) exp(ivw$b[1] - 1.96 * ivw$se[1]) else NA,
      U95   = if (nrow(ivw)) exp(ivw$b[1] + 1.96 * ivw$se[1]) else NA,
      P_ivw = if (nrow(ivw)) ivw$pval[1] else NA,
      OR_wmd = if (nrow(wmd)) exp(wmd$b[1]) else NA,
      P_wmd  = if (nrow(wmd)) wmd$pval[1] else NA,
      OR_egger = if (nrow(egg)) exp(egg$b[1]) else NA,
      P_egger  = if (nrow(egg)) egg$pval[1] else NA,
      egger_intercept_P = if (!is.null(x$pleiotropy) && nrow(x$pleiotropy))
                            x$pleiotropy$pval[1] else NA,
      Q_pval = if (!is.null(x$heterogeneity) && nrow(x$heterogeneity)) {
        idx <- which(x$heterogeneity$method == "Inverse variance weighted")[1]
        if (is.na(idx)) NA else x$heterogeneity$Q_pval[idx]
      } else NA,
      row.names = NULL
    )
  }
}
no_hla <- bind_rows(rows)
no_hla$P_fdr  <- p.adjust(no_hla$P_ivw, method = "fdr")
no_hla$P_bonf <- p.adjust(no_hla$P_ivw, method = "bonferroni")
no_hla <- no_hla[order(no_hla$P_ivw), ]
fwrite(no_hla, "results/MR_summary_noHLA.tsv", sep = "\t")
cat("\nWrote results/MR_summary_noHLA.tsv  (", nrow(no_hla), "pairs )\n")

# ---- Compare with HLA-inclusive --------------------------------------------
primary <- fread("results/MR_summary_table.tsv")
joined <- merge(no_hla[, c("exposure","outcome","n_snp","OR","L95","U95","P_ivw","P_fdr")],
                primary[, c("exposure","outcome","OR","L95","U95","P_ivw","P_fdr")],
                by = c("exposure","outcome"), suffixes = c("_noHLA","_primary"))
joined$delta_OR <- joined$OR_noHLA - joined$OR_primary
joined$P_changed_dir <- !is.na(joined$P_ivw_noHLA) & !is.na(joined$P_ivw_primary) &
  ((joined$P_ivw_noHLA < 0.05) != (joined$P_ivw_primary < 0.05))
joined <- joined[order(joined$P_ivw_primary), ]
fwrite(joined, "results/MR_HLA_comparison.tsv", sep = "\t")
cat("\nTop pairs by primary P (with HLA vs without HLA):\n")
print(head(joined[, c("exposure","outcome","n_snp","OR_primary","P_ivw_primary",
                       "OR_noHLA","P_ivw_noHLA","P_changed_dir")], 15))

# Headline: how many primary-P<0.05 signals survive HLA exclusion?
prim_sig <- joined[joined$P_ivw_primary < 0.05, ]
n_survive <- sum(!is.na(prim_sig$P_ivw_noHLA) & prim_sig$P_ivw_noHLA < 0.05)
cat("\nPrimary P<0.05 signals:", nrow(prim_sig), "\n")
cat("Survive HLA exclusion (P<0.05 still):", n_survive, "\n")
cat("Lost on HLA exclusion:", nrow(prim_sig) - n_survive, "\n")

# Save FDR-significant noHLA signals
sig <- subset(no_hla, !is.na(P_fdr) & P_fdr < 0.05)
fwrite(sig, "results/signals_FDR_noHLA.tsv", sep = "\t")
cat("\nFDR<0.05 signals after HLA exclusion:", nrow(sig), "\n")

cat("\nDone.\n")
