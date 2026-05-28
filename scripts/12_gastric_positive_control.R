#!/usr/bin/env Rscript
# 12_gastric_positive_control.R
#
# TASK 1 (Round 1 revision):
# Replicate gastric cancer signal — H. pylori's established causal target —
# using the same Butler-Laporte instruments and the Sakaue Pan-UKB gastric
# cancer GWAS (GCST90018849), matching Rao W et al. Cureus 2025;17:e89185.
#
# PASS criteria (per task list 1.3):
#   (a) OMP -> gastric cancer IVW P < 0.05
#   (b) OR > 1.0 (risk-increasing direction, matching Rao 2025)
#   (c) Weighted median same direction as IVW
#   (d) Weighted mode same direction (or NA acceptable)
#   (e) MR-Egger intercept P > 0.05 (no directional pleiotropy)
#   (f) Steiger directionality consistent (exposure -> outcome) or NA
#
# FAIL => phenome-wide null conclusion is underpowered;
# the Tier 1 narrative must be downgraded.

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(ggplot2)
})

set.seed(20260521)

dir.create("data/outcome", showWarnings = FALSE, recursive = TRUE)
dir.create("results",      showWarnings = FALSE, recursive = TRUE)
dir.create("figures",      showWarnings = FALSE, recursive = TRUE)
dir.create("logs",         showWarnings = FALSE, recursive = TRUE)

sink("logs/12_gastric_positive_control.log", split = TRUE)
on.exit(sink(NULL), add = TRUE)

cat("=== 12_gastric_positive_control.R  ", format(Sys.time()), " ===\n")
cat("PURPOSE: positive control replication — gastric cancer\n")
cat("REFERENCE: Rao W et al. Cureus 2025;17:e89185 used same Butler-Laporte\n")
cat("           antigens against gastric cancer GWAS GCST90018849.\n\n")

# ---- Load exposures ---------------------------------------------------------
exposure_list <- readRDS("data/exposure/exposure_list.rds")
ant_names <- names(exposure_list)
cat("Antigens loaded:", paste(ant_names, collapse = ", "), "\n")
all_snps <- unique(unlist(lapply(exposure_list, function(d) d$SNP)))
cat("Unique exposure SNPs:", length(all_snps), "\n\n")

# ---- Search OpenGWAS for gastric cancer candidates --------------------------
ao <- available_outcomes()
cat("Searching OpenGWAS for gastric cancer outcomes...\n")
gc_pool <- ao[grepl("gastric|stomach", ao$trait, ignore.case = TRUE), ]
gc_pool <- gc_pool[order(-gc_pool$ncase), ]
cat("Top 10 candidates:\n")
print(head(gc_pool[, c("id","trait","sample_size","ncase","population","year")], 10))

# Target ID per task list: GCST90018849 (Sakaue 2021 Pan-UKB)
target_id <- "ebi-a-GCST90018849"
target_info <- tryCatch(gwasinfo(target_id), error = function(e) NULL)
if (is.null(target_info) || nrow(target_info) == 0) {
  cat("\n!!! ebi-a-GCST90018849 not found on OpenGWAS\n")
  cat("Falling back to largest available EUR gastric cancer outcome\n")
  eur <- gc_pool[!is.na(gc_pool$population) & gc_pool$population == "European", ]
  eur <- eur[!is.na(eur$ncase) & eur$ncase >= 500, ]
  if (nrow(eur) == 0) stop("No usable gastric cancer outcome with ncase>=500")
  target_id <- eur$id[1]
  target_info <- gwasinfo(target_id)
}
cat("\nUSING outcome: ", target_id, "\n")
keep_cols <- intersect(c("id","trait","sample_size","ncase","population","year","unit"),
                       names(target_info))
print(target_info[, keep_cols, drop = FALSE])

# ---- Extract outcome data ---------------------------------------------------
cat("\nExtracting outcome data for", length(all_snps), "exposure SNPs...\n")
gc_outcome <- tryCatch(
  extract_outcome_data(snps = all_snps, outcomes = target_id, proxies = TRUE),
  error = function(e) { cat("ERR:", e$message, "\n"); NULL }
)
if (is.null(gc_outcome) || nrow(gc_outcome) == 0) {
  stop("Failed to extract outcome data")
}
cat("Outcome rows fetched:", nrow(gc_outcome), "\n")
saveRDS(gc_outcome, "data/outcome/gc_outcome_eur.rds")

# ---- Per-antigen MR ---------------------------------------------------------
run_pair <- function(exp_df, out_df, exp_lab) {
  if (is.null(exp_df) || nrow(exp_df) == 0) return(NULL)
  dat <- tryCatch(harmonise_data(exp_df, out_df, action = 2),
                  error = function(e) { cat("  harm err:", e$message, "\n"); NULL })
  if (is.null(dat)) return(NULL)
  dat <- subset(dat, mr_keep == TRUE)
  if (nrow(dat) < 3) { cat("  <3 SNPs, skip\n"); return(NULL) }

  # Steiger
  dat2 <- try(steiger_filtering(dat), silent = TRUE)
  steiger_ok <- !inherits(dat2, "try-error") && is.data.frame(dat2)
  if (steiger_ok) dat <- dat2

  res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression",
                                   "mr_weighted_median","mr_weighted_mode"))
  pleio <- try(mr_pleiotropy_test(dat), silent = TRUE)
  het   <- try(mr_heterogeneity(dat),    silent = TRUE)

  presso <- NULL
  if (nrow(dat) >= 4 && requireNamespace("MRPRESSO", quietly = TRUE)) {
    presso <- try(MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome   = "se.outcome",   SdExposure   = "se.exposure",
      OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
      data = as.data.frame(dat),
      NbDistribution = 1000, SignifThreshold = 0.05
    ), silent = TRUE)
  }

  list(exposure = exp_lab, n_snp = nrow(dat), dat = dat, mr_res = res,
       pleiotropy = pleio, heterogeneity = het, presso = presso)
}

cat("\nRunning per-antigen MR vs gastric cancer...\n")
gc_results <- list()
for (lab in ant_names) {
  cat("\n--- ", lab, " ---\n")
  gc_results[[lab]] <- run_pair(exposure_list[[lab]], gc_outcome, lab)
}
saveRDS(gc_results, "results/gc_positive_control.rds")

# ---- Summary table ----------------------------------------------------------
pick <- function(res, m) {
  r <- res[res$method == m, ]
  if (!nrow(r)) return(c(b = NA, se = NA, p = NA))
  c(b = r$b[1], se = r$se[1], p = r$pval[1])
}
rows <- lapply(gc_results, function(x) {
  if (is.null(x)) return(NULL)
  ivw <- pick(x$mr_res, "Inverse variance weighted")
  egg <- pick(x$mr_res, "MR Egger")
  wmd <- pick(x$mr_res, "Weighted median")
  wmo <- pick(x$mr_res, "Weighted mode")
  pleio_p <- if (!is.null(x$pleiotropy) && !inherits(x$pleiotropy, "try-error") &&
                  is.data.frame(x$pleiotropy) && nrow(x$pleiotropy) > 0) x$pleiotropy$pval[1] else NA
  q_p <- if (!is.null(x$heterogeneity) && !inherits(x$heterogeneity, "try-error") &&
              is.data.frame(x$heterogeneity)) {
    idx <- which(x$heterogeneity$method == "Inverse variance weighted")[1]
    if (is.na(idx)) NA else x$heterogeneity$Q_pval[idx]
  } else NA
  presso_p <- NA
  if (!is.null(x$presso) && !inherits(x$presso, "try-error")) {
    v <- x$presso$`MR-PRESSO results`$`Global Test`$Pvalue
    if (is.character(v)) v <- suppressWarnings(as.numeric(sub("^<","",v)))
    presso_p <- as.numeric(v)
  }
  data.frame(
    antigen = x$exposure, n_snp = x$n_snp,
    OR_ivw   = exp(ivw["b"]), L95 = exp(ivw["b"] - 1.96*ivw["se"]),
    U95      = exp(ivw["b"] + 1.96*ivw["se"]),
    P_ivw    = ivw["p"],
    OR_egger = exp(egg["b"]), P_egger = egg["p"],
    OR_wmd   = exp(wmd["b"]), P_wmd  = wmd["p"],
    OR_wmode = exp(wmo["b"]), P_wmode = wmo["p"],
    egger_intercept_P = pleio_p,
    Q_P = q_p,
    PRESSO_global_P = presso_p,
    row.names = NULL
  )
})
summary_tbl <- do.call(rbind, rows)
summary_tbl <- summary_tbl[order(summary_tbl$P_ivw), ]
fwrite(summary_tbl, "results/gc_positive_control_summary.tsv", sep = "\t")
cat("\n=== Per-antigen vs Gastric cancer ===\n")
print(summary_tbl)

# ---- PASS/FAIL decision (focused on OMP) ------------------------------------
omp <- summary_tbl[summary_tbl$antigen == "omp", ]
safe_get <- function(df, col, idx = 1, default = NA) {
  if (!is.data.frame(df) || !col %in% names(df) || nrow(df) < idx) return(default)
  df[[col]][idx]
}
decision_lines <- c(
  "=== GASTRIC CANCER POSITIVE CONTROL — PASS/FAIL DECISION ===",
  paste0("Timestamp: ", format(Sys.time())),
  paste0("Outcome GWAS: ", target_id, "  (", safe_get(target_info, "trait"), ")"),
  paste0("Outcome N=", safe_get(target_info, "sample_size"),
         "  ncase=", safe_get(target_info, "ncase")),
  ""
)
if (nrow(omp) == 0) {
  decision_lines <- c(decision_lines,
    "RESULT: OMP row missing — pipeline failed before reaching OMP.",
    "VERDICT: FAIL (pipeline error, not a power issue).")
} else {
  c1 <- !is.na(omp$P_ivw)   && omp$P_ivw   < 0.05
  c2 <- !is.na(omp$OR_ivw)  && omp$OR_ivw  > 1.0
  c3 <- !is.na(omp$OR_wmd)  && omp$OR_wmd  > 1.0
  c4 <- is.na(omp$OR_wmode) || omp$OR_wmode > 1.0     # NA acceptable
  c5 <- is.na(omp$egger_intercept_P) || omp$egger_intercept_P > 0.05
  pass <- c1 && c2 && c3 && c4 && c5
  decision_lines <- c(decision_lines,
    sprintf("OMP n_SNP                = %d", omp$n_snp),
    sprintf("OMP -> GC IVW OR (95CI)  = %.3f (%.3f - %.3f)",
            omp$OR_ivw, omp$L95, omp$U95),
    sprintf("OMP -> GC IVW P          = %.4g  -> %s",
            omp$P_ivw, ifelse(c1, "PASS", "FAIL")),
    sprintf("OMP -> GC IVW OR > 1     = %s   -> %s",
            ifelse(c2, "yes", "no"), ifelse(c2, "PASS", "FAIL")),
    sprintf("Weighted median OR > 1   = %.3f -> %s",
            ifelse(is.na(omp$OR_wmd), -Inf, omp$OR_wmd), ifelse(c3, "PASS", "FAIL")),
    sprintf("Weighted mode  OR > 1    = %s   -> %s",
            ifelse(is.na(omp$OR_wmode), "NA", sprintf("%.3f", omp$OR_wmode)),
            ifelse(c4, "PASS (or NA accepted)", "FAIL")),
    sprintf("Egger intercept P > 0.05 = %s   -> %s",
            ifelse(is.na(omp$egger_intercept_P), "NA",
                   sprintf("%.4g", omp$egger_intercept_P)),
            ifelse(c5, "PASS", "FAIL")),
    "",
    paste0("OVERALL: ", ifelse(pass, "PASS — proceed to TASK 2-4",
                                "FAIL — investigate before any further work"))
  )
}
writeLines(decision_lines, "logs/gc_positive_control_decision.txt")
cat("\n", paste(decision_lines, collapse = "\n"), "\n", sep = "")

# ---- Forest plot ------------------------------------------------------------
fp <- summary_tbl[!is.na(summary_tbl$OR_ivw), ]
fp$antigen <- factor(fp$antigen, levels = fp$antigen[order(fp$OR_ivw)])
g <- ggplot(fp, aes(x = OR_ivw, y = antigen)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = L95, xmax = U95), height = 0.2, linewidth = 0.7) +
  geom_point(size = 3.2, color = "#e6550d") +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)\nP=%.3g", OR_ivw, L95, U95, P_ivw)),
            hjust = -0.05, size = 3, family = "sans") +
  scale_x_log10() +
  labs(title = "Positive control: H. pylori antigens -> Gastric cancer (IVW)",
       subtitle = sprintf("Outcome: %s (ncase=%s)",
                           safe_get(target_info, "trait"),
                           safe_get(target_info, "ncase")),
       x = "OR (log scale, 95% CI)", y = "Antigen") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave("figures/gc_positive_control_forest.png", g, width = 9, height = 4.5, dpi = 300)
cat("\nForest plot: figures/gc_positive_control_forest.png\n")
cat("\nDone.\n")
