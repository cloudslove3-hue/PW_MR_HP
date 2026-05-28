#!/usr/bin/env Rscript
# 08_mvmr.R  Multivariable MR with BMI / smoking / education as covariates
#
# Strategy:
#   For every pair with IVW P < 0.1 in the (post-A) primary table, run MVMR
#   conditioning on BMI, smoking initiation, and years of education. Report
#   direct effect of the antigen on the outcome controlling for these.
#
# Covariate GWAS:
#   BMI:                 ieu-b-40         Yengo 2018 UKB+GIANT, n~700k
#   Smoking initiation:  ieu-b-4877       GSCAN Liu 2019, n~1.2M
#   Education years:     ieu-a-1239       Okbay 2016, n~290k
#
# Input:   results/MR_summary_table.tsv (post-A)
#          data/exposure/exposure_list.rds
# Output:  results/mvmr_results.tsv
#          logs/08_mvmr.log

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("logs",    showWarnings = FALSE, recursive = TRUE)
sink("logs/08_mvmr.log", split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 08_mvmr.R  ", format(Sys.time()), " ===\n")

# ---- Covariates -------------------------------------------------------------
cov_ids <- c(bmi = "ieu-b-40",
             smoke = "ieu-b-4877",
             edu   = "ieu-a-1239")

# ---- Antigen -> OpenGWAS id map --------------------------------------------
ant_ids <- c(hpylori_general = "ebi-a-GCST90006910",
             cagA            = "ebi-a-GCST90006911",
             catalase        = "ebi-a-GCST90006912",
             groEL           = "ebi-a-GCST90006913",
             omp             = "ebi-a-GCST90006914",
             ureA            = "ebi-a-GCST90006915",
             vacA            = "ebi-a-GCST90006916")

# ---- Outcome label -> id map (read from config) ----------------------------
read_cfg <- function(p) {
  raw <- readLines(p, warn = FALSE)
  raw <- raw[!grepl("^\\s*#", raw) & nzchar(trimws(raw))]
  fread(text = paste(raw, collapse = "\n"), sep = "\t")
}
cfg <- read_cfg("config/outcomes_tier1.tsv")
out_ids <- setNames(cfg$id, cfg$outcome_label)

# ---- Candidate pairs (IVW P < 0.1) -----------------------------------------
st <- fread("results/MR_summary_table.tsv")
cand <- st[!is.na(P_ivw) & P_ivw < 0.1, ][order(P_ivw)]
cat("Candidate pairs (P_ivw<0.1):", nrow(cand), "\n")
if (nrow(cand) == 0) {
  cat("No candidates; nothing to do.\n")
  quit(status = 0)
}
print(head(cand[, c("exposure","outcome","n_snp","OR","P_ivw","P_fdr")], 20))

# ---- MVMR per candidate -----------------------------------------------------
run_mvmr <- function(exp_lab, out_lab) {
  exp_id <- ant_ids[exp_lab]
  out_id <- out_ids[out_lab]
  if (is.na(exp_id) || is.na(out_id) || is.null(exp_id) || is.null(out_id)) {
    cat("  unknown id\n"); return(NULL)
  }
  ids <- c(exp_id, cov_ids)
  cat("\n>> MVMR:", exp_lab, "+", paste(names(cov_ids), collapse = "+"),
      " ->", out_lab, "\n")
  mv_exp <- tryCatch(
    mv_extract_exposures(id_exposure = ids, pval_threshold = 5e-6,
                          clump_r2 = 0.001, clump_kb = 10000),
    error = function(e) { cat("  mv_extract_exposures err:", e$message, "\n"); NULL }
  )
  if (is.null(mv_exp) || nrow(mv_exp) == 0) return(NULL)
  cat("  mv_exp SNPs:", length(unique(mv_exp$SNP)), "\n")

  mv_out <- tryCatch(
    extract_outcome_data(snps = unique(mv_exp$SNP), outcomes = out_id, proxies = TRUE),
    error = function(e) { cat("  mv extract_outcome err:", e$message, "\n"); NULL }
  )
  if (is.null(mv_out) || nrow(mv_out) == 0) return(NULL)

  dat <- tryCatch(mv_harmonise_data(mv_exp, mv_out),
                  error = function(e) { cat("  mv_harmonise err:", e$message, "\n"); NULL })
  if (is.null(dat)) return(NULL)

  res <- tryCatch(mv_multiple(dat),
                  error = function(e) { cat("  mv_multiple err:", e$message, "\n"); NULL })
  if (is.null(res)) return(NULL)
  r <- res$result
  r$OR  <- exp(r$b)
  r$L95 <- exp(r$b - 1.96 * r$se)
  r$U95 <- exp(r$b + 1.96 * r$se)
  r$mvmr_exposure_label <- exp_lab
  r$mvmr_outcome_label  <- out_lab
  r$is_antigen <- r$id.exposure == exp_id
  r
}

all_mvmr <- list()
for (i in seq_len(nrow(cand))) {
  res <- run_mvmr(cand$exposure[i], cand$outcome[i])
  if (!is.null(res)) {
    all_mvmr[[paste(cand$exposure[i], cand$outcome[i], sep = "__")]] <- res
  }
  Sys.sleep(2)
}

mvmr_tbl <- bind_rows(all_mvmr)
if (nrow(mvmr_tbl)) {
  cols <- c("mvmr_exposure_label","mvmr_outcome_label","exposure","is_antigen",
            "nsnp","b","se","pval","OR","L95","U95")
  cols <- intersect(cols, names(mvmr_tbl))
  fwrite(mvmr_tbl[, cols], "results/mvmr_results.tsv", sep = "\t")
  cat("\nWrote results/mvmr_results.tsv  (", nrow(mvmr_tbl), "rows )\n")
  cat("\n=== Antigen direct effects (after BMI/smoking/edu adjustment) ===\n")
  print(mvmr_tbl[mvmr_tbl$is_antigen,
                  c("mvmr_exposure_label","mvmr_outcome_label","nsnp",
                    "OR","L95","U95","pval")])
} else {
  cat("\n!!! No MVMR results produced\n")
}
cat("\nDone.\n")
