#!/usr/bin/env Rscript
# 99_qc.R  Automated QC checklist per md spec section "검증 단계에서 주의할 점"
#
# Inputs:
#   data/exposure/exposure_list.rds
#   results/all_mr_results.rds
#   results/MR_summary_table.tsv
#
# Output:
#   logs/qc_report.txt   human-readable
#   results/qc_flags.tsv machine-readable per-pair flags

suppressPackageStartupMessages({
  library(data.table)
})

dir.create("logs",    showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

report_path <- "logs/qc_report.txt"
sink(report_path, split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== QC Report  ", format(Sys.time()), " ===\n\n")

exposure_list <- readRDS("data/exposure/exposure_list.rds")
all_results   <- readRDS("results/all_mr_results.rds")
st            <- fread("results/MR_summary_table.tsv")

# --- Check 1: F-stat >= 10 ---------------------------------------------------
cat("### 1. F-statistic (>=10) per antigen ###\n")
for (nm in names(exposure_list)) {
  d <- exposure_list[[nm]]
  if (is.null(d)) { cat("  ", nm, ": MISSING\n"); next }
  weak <- sum(d$F_stat < 10)
  cat(sprintf("  %-18s n=%d  meanF=%.1f  minF=%.1f  weak(<10)=%d  %s\n",
              nm, nrow(d), mean(d$F_stat), min(d$F_stat), weak,
              if (weak > 0) "[WARN]" else "OK"))
}

# --- Check 2: SNP count >= 3 per antigen ------------------------------------
cat("\n### 2. SNP count >=3 per antigen ###\n")
for (nm in names(exposure_list)) {
  d <- exposure_list[[nm]]
  n <- if (is.null(d)) 0 else nrow(d)
  cat(sprintf("  %-18s n=%d  %s\n", nm, n, if (n < 3) "[FAIL]" else "OK"))
}

# --- Check 3: harmonization losses ------------------------------------------
cat("\n### 3. Harmonization: kept vs supplied per pair ###\n")
for (k in names(all_results)) {
  x <- all_results[[k]]
  if (is.null(x)) next
  cat(sprintf("  %-50s n_after=%d\n", k, x$n_snp))
}

# --- Check 4-5: Egger intercept + Q ----------------------------------------
cat("\n### 4. MR-Egger intercept P<0.05 (horizontal pleiotropy) ###\n")
flag_egg <- subset(st, !is.na(P_egger_intercept) & P_egger_intercept < 0.05)
if (nrow(flag_egg)) {
  print(flag_egg[, c("exposure","outcome","n_snp","P_egger_intercept")])
} else cat("  none\n")

cat("\n### 5. Cochran's Q P<0.05 (heterogeneity) ###\n")
flag_q <- subset(st, !is.na(Q_pval) & Q_pval < 0.05)
if (nrow(flag_q)) {
  print(flag_q[, c("exposure","outcome","n_snp","Q_pval")])
} else cat("  none\n")

# --- Check 6: IVW vs weighted median direction concordance -----------------
cat("\n### 6. IVW vs weighted median direction mismatch ###\n")
mismatch <- subset(st, !is.na(OR) & !is.na(OR_wmedian) &
                      sign(log(OR)) != sign(log(OR_wmedian)))
if (nrow(mismatch)) {
  print(mismatch[, c("exposure","outcome","OR","OR_wmedian","P_ivw","P_wmedian")])
} else cat("  none\n")

# --- Check 7: MR-PRESSO global test P<0.05 ----------------------------------
cat("\n### 7. MR-PRESSO global test P<0.05 ###\n")
flag_p <- subset(st, !is.na(presso_global_p) & presso_global_p < 0.05)
if (nrow(flag_p)) {
  print(flag_p[, c("exposure","outcome","n_snp","presso_global_p","P_ivw")])
} else cat("  none\n")

# --- Check 8: FDR-significant signals --------------------------------------
cat("\n### 8. FDR-significant (P_fdr<0.05) signals ###\n")
sig <- subset(st, !is.na(P_fdr) & P_fdr < 0.05)
cat("  ", nrow(sig), "signals\n")
if (nrow(sig)) print(sig[, c("exposure","outcome","OR","L95","U95","P_ivw","P_fdr")])

# --- Per-pair flags table ---------------------------------------------------
flags <- data.frame(
  pair = paste(st$exposure, st$outcome, sep = "__"),
  exposure = st$exposure, outcome = st$outcome, n_snp = st$n_snp,
  flag_few_snp     = st$n_snp < 3,
  flag_egger_intc  = !is.na(st$P_egger_intercept) & st$P_egger_intercept < 0.05,
  flag_q           = !is.na(st$Q_pval) & st$Q_pval < 0.05,
  flag_dir_mismatch = !is.na(st$OR) & !is.na(st$OR_wmedian) &
                       sign(log(st$OR)) != sign(log(st$OR_wmedian)),
  flag_presso      = !is.na(st$presso_global_p) & st$presso_global_p < 0.05,
  flag_fdr_sig     = !is.na(st$P_fdr) & st$P_fdr < 0.05
)
fwrite(flags, "results/qc_flags.tsv", sep = "\t")
cat("\nWrote results/qc_flags.tsv\n")
cat("=== QC done ===\n")
