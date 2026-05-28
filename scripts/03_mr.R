#!/usr/bin/env Rscript
# 03_mr.R  Harmonization + primary MR + sensitivity
#
# Inputs:
#   data/exposure/exposure_list.rds
#   data/outcome/outcome_data_list.rds
#
# Output:
#   results/all_mr_results.rds
#
# Methods per pair:
#   IVW, MR-Egger, weighted median, weighted mode
#   pleiotropy_test (Egger intercept)
#   heterogeneity   (Cochran Q)
#   Steiger filtering (skipped pairs with NA direction kept conservatively)
#   MR-PRESSO (only when n_SNP >= 4)

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("logs",    showWarnings = FALSE, recursive = TRUE)
log_path <- file.path("logs", "03_mr.log")
sink(log_path, split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 03_mr.R  ", format(Sys.time()), " ===\n")

exposure_list     <- readRDS("data/exposure/exposure_list.rds")
outcome_data_list <- readRDS("data/outcome/outcome_data_list.rds")

has_presso <- requireNamespace("MRPRESSO", quietly = TRUE)
if (!has_presso) cat("!!! MRPRESSO not installed; skipping PRESSO step\n")

run_mr_pair <- function(exp_df, out_df, exp_lab, out_lab) {
  if (is.null(out_df) || nrow(out_df) == 0) return(NULL)

  dat <- tryCatch(harmonise_data(exp_df, out_df, action = 2),
                  error = function(e) { cat("  harmonise err:", e$message, "\n"); NULL })
  if (is.null(dat)) return(NULL)
  n_pre <- nrow(dat)
  dat <- subset(dat, mr_keep == TRUE)
  n_kept <- nrow(dat)
  cat("  harmonised:", n_pre, " -> kept:", n_kept, "\n")
  if (n_kept < 3) { cat("  <3 SNPs, skipping\n"); return(NULL) }

  # Steiger filter (best-effort)
  dat2 <- try(steiger_filtering(dat), silent = TRUE)
  if (!inherits(dat2, "try-error") && is.data.frame(dat2)) dat <- dat2

  # Primary MR
  res <- try(mr(dat, method_list = c("mr_ivw","mr_egger_regression",
                                       "mr_weighted_median","mr_weighted_mode")),
             silent = TRUE)
  if (inherits(res, "try-error")) { cat("  mr() failed\n"); return(NULL) }

  pleio <- try(mr_pleiotropy_test(dat),  silent = TRUE)
  het   <- try(mr_heterogeneity(dat),    silent = TRUE)

  presso <- NULL
  if (has_presso && nrow(dat) >= 4) {
    presso <- try(MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome   = "se.outcome",   SdExposure   = "se.exposure",
      OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
      data = as.data.frame(dat),
      NbDistribution = 1000, SignifThreshold = 0.05
    ), silent = TRUE)
  }

  list(exposure = exp_lab, outcome = out_lab,
       dat = dat, mr_res = res,
       pleiotropy = pleio, heterogeneity = het,
       presso = presso, n_snp = nrow(dat))
}

all_results <- list()
for (exp_lab in names(exposure_list)) {
  exp_df <- exposure_list[[exp_lab]]
  if (is.null(exp_df) || nrow(exp_df) == 0) next
  for (out_lab in names(outcome_data_list)) {
    out_df <- outcome_data_list[[out_lab]]
    key <- paste0(exp_lab, "__", out_lab)
    cat("\n>>", key, "\n")
    all_results[[key]] <- run_mr_pair(exp_df, out_df, exp_lab, out_lab)
  }
}

saveRDS(all_results, "results/all_mr_results.rds")
cat("\nDone. Pairs run:", length(all_results),
    " non-null:", sum(!sapply(all_results, is.null)), "\n")
cat("Proceed to scripts/04_summary.R\n")
