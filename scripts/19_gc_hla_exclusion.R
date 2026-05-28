#!/usr/bin/env Rscript
# 19_gc_hla_exclusion.R
# Test whether OMP -> gastric cancer (positive control) survives HLA exclusion.
# If yes: positive control is HLA-independent -> Tier 1 HLA-exclusion falsification is decisive
# If no: both genuine and pleiotropic OMP signals depend on HLA -> reframe interpretation

set.seed(20260521)
suppressPackageStartupMessages({
  library(TwoSampleMR); library(data.table)
})
sink("logs/19_gc_hla_exclusion.log", split = TRUE); on.exit(sink(NULL), add = TRUE)
cat("=== 19_gc_hla_exclusion.R  ", format(Sys.time()), " ===\n")

omp <- readRDS("data/exposure/exposure_list.rds")$omp
gc  <- readRDS("data/outcome/gc_outcome_eur.rds")

HLA_lo <- 25e6; HLA_hi <- 35e6
omp_hla    <- omp[omp$chr.exposure == "6" &
                   omp$pos.exposure >= HLA_lo & omp$pos.exposure <= HLA_hi, ]
omp_noHLA  <- omp[!(omp$chr.exposure == "6" &
                     omp$pos.exposure >= HLA_lo & omp$pos.exposure <= HLA_hi), ]
cat("OMP SNPs total:", nrow(omp), "  HLA:", nrow(omp_hla),
    "  non-HLA:", nrow(omp_noHLA), "\n")
cat("Dropped HLA SNPs:", paste(omp_hla$SNP, collapse = ", "), "\n\n")

# MR with HLA-only and non-HLA-only
run_subset <- function(exp_df, lab) {
  if (nrow(exp_df) < 3) { cat(lab, ": <3 SNPs, skip\n"); return(NULL) }
  dat <- harmonise_data(exp_df, gc, action = 2)
  dat <- subset(dat, mr_keep == TRUE)
  cat(lab, ": harmonized SNPs =", nrow(dat), "\n")
  if (nrow(dat) < 3) return(NULL)
  res <- mr(dat, method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression"))
  print(res)
  res
}

cat("\n--- OMP (all 10 SNPs) vs Gastric cancer ---\n")
r_all <- run_subset(omp, "all")
cat("\n--- OMP (non-HLA, 8 SNPs) vs Gastric cancer ---\n")
r_noHLA <- run_subset(omp_noHLA, "non-HLA")
cat("\n--- OMP (HLA-only, 2 SNPs) vs Gastric cancer ---\n")
r_HLA <- run_subset(omp_hla, "HLA-only")

# Conclusion
cat("\n=== INTERPRETATION ===\n")
pick_ivw <- function(r) if (!is.null(r) && nrow(subset(r, method=="Inverse variance weighted"))) {
    x <- subset(r, method=="Inverse variance weighted")
    sprintf("OR=%.3f (%.3f-%.3f) P=%.4g",
            exp(x$b), exp(x$b-1.96*x$se), exp(x$b+1.96*x$se), x$pval)
  } else "NA"
cat("All     : ", pick_ivw(r_all), "\n")
cat("non-HLA : ", pick_ivw(r_noHLA), "\n")
cat("HLA-only: ", pick_ivw(r_HLA), "\n")

verdict <- {
  if (!is.null(r_noHLA) &&
      nrow(subset(r_noHLA, method == "Inverse variance weighted")) &&
      subset(r_noHLA, method == "Inverse variance weighted")$pval < 0.05) {
    "POSITIVE CONTROL SURVIVES HLA EXCLUSION -> robust; HLA exclusion is a clean test"
  } else {
    "POSITIVE CONTROL FAILS HLA EXCLUSION -> both genuine (GC) and candidate (Crohn) OMP signals depend on HLA; HLA exclusion cannot distinguish them; triangulation with eradication trial evidence becomes the decisive sensitivity test"
  }
}
cat("VERDICT: ", verdict, "\n")
