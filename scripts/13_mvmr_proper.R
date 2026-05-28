#!/usr/bin/env Rscript
# 13_mvmr_proper.R
#
# TASK 2 (Round 1 revision):
# Re-run multivariable MR using the dedicated MVMR package (Sanderson 2019,
# *Int J Epidemiol* 48:713-727) to obtain proper MVMR-IVW estimates *and*
# conditional F-statistics. Previous implementation via TwoSampleMR::
# mv_multiple gave nsnp=0 for antigen exposures (unidentified state).
#
# Design:
#   - For each of 12 nominally-suggestive Tier 1 pairs (univariate P<0.10):
#       exposures = (antigen, BMI, smoking initiation, education years)
#       outcome   = the Tier 1 outcome
#   - Instrument selection: antigen SNPs at P<5e-6 (existing instruments),
#     plus covariate top SNPs from each covariate GWAS. Combined set is
#     clumped together (r2<0.001, 10Mb) against EUR LD reference so no two
#     instruments are correlated.
#   - For each kept SNP, fetch beta/SE from all 4 exposure GWAS and from the
#     outcome GWAS. Drop SNPs missing in any.
#   - Run ivw_mvmr() and strength_mvmr().
#   - If any conditional F < 10, flag that exposure's estimate as unreliable.
#
# Output:
#   results/mvmr_proper_results.tsv  — 12 pair x {OR_antigen, CI, P,
#                                                  F_cond per exposure, n_SNP}
#   logs/13_mvmr_proper.log

set.seed(20260521)

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(dplyr)
})

# Install MVMR package if missing
if (!requireNamespace("MVMR", quietly = TRUE)) {
  cat("Installing MVMR via mrcieu r-universe...\n")
  try(install.packages("MVMR",
    repos = c("https://mrcieu.r-universe.dev",
              "https://wspiller.r-universe.dev",
              "https://cloud.r-project.org")), silent = TRUE)
}
if (!requireNamespace("MVMR", quietly = TRUE)) {
  cat("Falling back to anonymous GitHub install_github...\n")
  Sys.setenv(GITHUB_PAT = "", GITHUB_TOKEN = "")
  Sys.unsetenv(c("GITHUB_PAT","GITHUB_TOKEN"))
  try(remotes::install_github("WSpiller/MVMR", upgrade = "never",
                               auth_token = NULL))
}
if (!requireNamespace("MVMR", quietly = TRUE)) {
  cat("!!! MVMR install failed. Falling back to manual MVMR-IVW + conditional-F.\n")
  cat("    Will implement Sanderson 2019 equations directly.\n")
} else {
  suppressPackageStartupMessages(library(MVMR))
}
HAS_MVMR <- requireNamespace("MVMR", quietly = TRUE)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("logs",    showWarnings = FALSE, recursive = TRUE)
sink("logs/13_mvmr_proper.log", split = TRUE)
on.exit(sink(NULL), add = TRUE)

cat("=== 13_mvmr_proper.R  ", format(Sys.time()), " ===\n")

# ---- Antigen + outcome registry --------------------------------------------
ant_ids <- c(hpylori_general = "ebi-a-GCST90006910",
             cagA            = "ebi-a-GCST90006911",
             catalase        = "ebi-a-GCST90006912",
             groEL           = "ebi-a-GCST90006913",
             omp             = "ebi-a-GCST90006914",
             ureA            = "ebi-a-GCST90006915",
             vacA            = "ebi-a-GCST90006916")
cov_ids <- c(bmi   = "ieu-b-40",
             smoke = "ieu-b-4877",
             edu   = "ieu-a-1239")

# Read outcome config for label -> id mapping
raw <- readLines("config/outcomes_tier1.tsv", warn = FALSE)
raw <- raw[!grepl("^\\s*#", raw) & nzchar(trimws(raw))]
cfg <- fread(text = paste(raw, collapse = "\n"), sep = "\t")
out_ids <- setNames(cfg$id, cfg$outcome_label)

# Candidate pairs from Tier 1 univariate P<0.10
st <- fread("results/MR_summary_table.tsv")
cand <- st[!is.na(P_ivw) & P_ivw < 0.10][order(P_ivw)]
cat("Candidate pairs (univariate P<0.10):", nrow(cand), "\n")
print(cand[, c("exposure","outcome","n_snp","OR","P_ivw")])

# ---- Helper: fetch betas for a SNP list from a GWAS -----------------------
fetch_betas <- function(snps, id, label) {
  d <- tryCatch(
    associations(variants = snps, id = id, proxies = FALSE),
    error = function(e) { cat("  fetch err (",label,"):", e$message, "\n"); NULL }
  )
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d <- d[!is.na(d$beta) & !is.na(d$se) & d$se > 0, ]
  data.frame(SNP = d$rsid, beta = d$beta, se = d$se,
             ea = d$ea, oa = d$nea, eaf = d$eaf, p = d$p,
             stringsAsFactors = FALSE)
}

harmonise_alleles <- function(ref, qry) {
  # Align qry to ref (effect allele).
  m <- merge(ref[, c("SNP","ea","oa")], qry,
             by = "SNP", suffixes = c(".ref", ".qry"))
  flipped <- m$ea.ref != m$ea  # allele opposite
  m$beta_aligned <- ifelse(flipped, -m$beta, m$beta)
  m$se_aligned   <- m$se
  m[, c("SNP","beta_aligned","se_aligned")]
}

# ---- Pre-fetch covariate top instruments (at p<5e-8) -----------------------
cat("\nFetching covariate top SNPs (p<5e-8, clumped)...\n")
bmi_top <- tryCatch(extract_instruments("ieu-b-40", p1 = 5e-8,
                                          clump = TRUE, r2 = 0.001, kb = 10000),
                     error = function(e) NULL)
smk_top <- tryCatch(extract_instruments("ieu-b-4877", p1 = 5e-8,
                                          clump = TRUE, r2 = 0.001, kb = 10000),
                     error = function(e) NULL)
edu_top <- tryCatch(extract_instruments("ieu-a-1239", p1 = 5e-8,
                                          clump = TRUE, r2 = 0.001, kb = 10000),
                     error = function(e) NULL)
cat("BMI:", if (is.null(bmi_top)) "NA" else nrow(bmi_top),
    " smoke:", if (is.null(smk_top)) "NA" else nrow(smk_top),
    " edu:", if (is.null(edu_top)) "NA" else nrow(edu_top), "\n")

# Load antigen instruments
exposure_list <- readRDS("data/exposure/exposure_list.rds")

# ---- Run proper MVMR per candidate pair -----------------------------------
run_proper_mvmr <- function(ant_lab, out_lab) {
  cat("\n>> Proper MVMR:", ant_lab, "->", out_lab, "\n")
  ant <- exposure_list[[ant_lab]]
  ant_id <- ant_ids[ant_lab]
  out_id <- out_ids[out_lab]
  if (is.null(ant) || is.null(ant_id) || is.null(out_id)) return(NULL)

  # Combine SNP pool
  pool <- unique(c(ant$SNP, bmi_top$SNP, smk_top$SNP, edu_top$SNP))
  pool <- pool[!is.na(pool)]
  cat("  pool SNPs (pre-clump):", length(pool), "\n")
  if (length(pool) < 4) return(NULL)

  # Joint clump at r2<0.001 to enforce mutual independence
  pool_df <- data.frame(rsid = pool, pval = rep(1e-10, length(pool)),
                         id = "joint")
  kept <- tryCatch(
    ld_clump(pool_df, clump_r2 = 0.001, clump_kb = 10000, pop = "EUR"),
    error = function(e) { cat("  clump err:", e$message, "\n"); pool_df }
  )
  snps <- kept$rsid
  cat("  after joint clumping:", length(snps), "\n")

  # Fetch betas from all 5 GWAS (4 exposures + 1 outcome)
  b_ant <- fetch_betas(snps, ant_id, ant_lab); Sys.sleep(1)
  b_bmi <- fetch_betas(snps, "ieu-b-40", "bmi"); Sys.sleep(1)
  b_smk <- fetch_betas(snps, "ieu-b-4877", "smk"); Sys.sleep(1)
  b_edu <- fetch_betas(snps, "ieu-a-1239", "edu"); Sys.sleep(1)
  b_out <- fetch_betas(snps, out_id, out_lab); Sys.sleep(1)

  if (any(sapply(list(b_ant, b_bmi, b_smk, b_edu, b_out), is.null))) {
    cat("  one or more fetches empty -> abort\n"); return(NULL)
  }

  # Use antigen as reference for allele alignment
  ref <- b_ant[, c("SNP","ea","oa","beta","se")]
  align <- function(other) {
    m <- merge(b_ant[, c("SNP","ea","oa")], other,
               by = "SNP", suffixes = c(".ref",""))
    flipped <- m$ea.ref != m$ea
    m$beta <- ifelse(flipped, -m$beta, m$beta)
    m[, c("SNP","beta","se")]
  }
  ah_bmi <- align(b_bmi); ah_smk <- align(b_smk)
  ah_edu <- align(b_edu); ah_out <- align(b_out)

  # Intersect on common SNPs
  common <- Reduce(intersect, list(b_ant$SNP, ah_bmi$SNP, ah_smk$SNP,
                                     ah_edu$SNP, ah_out$SNP))
  cat("  SNPs in all 5 GWAS:", length(common), "\n")
  if (length(common) < 4) { cat("  <4, abort\n"); return(NULL) }

  o <- function(d) d[match(common, d$SNP), ]
  d_ant <- o(b_ant); d_bmi <- o(ah_bmi); d_smk <- o(ah_smk)
  d_edu <- o(ah_edu); d_out <- o(ah_out)

  # Build MVMR object
  BXGs   <- cbind(d_ant$beta, d_bmi$beta, d_smk$beta, d_edu$beta)
  seBXGs <- cbind(d_ant$se,   d_bmi$se,   d_smk$se,   d_edu$se)
  BYG    <- d_out$beta
  seBYG  <- d_out$se
  colnames(BXGs)   <- c("antigen","BMI","smoking","education")
  colnames(seBXGs) <- c("antigen","BMI","smoking","education")

  # Use MVMR package if available, else fallback to manual implementation
  if (HAS_MVMR) {
    mv <- tryCatch(
      MVMR::format_mvmr(BXGs = BXGs, BYG = BYG, seBXGs = seBXGs, seBYG = seBYG,
                  RSID = common),
      error = function(e) { cat("  format_mvmr err:", e$message, "\n"); NULL }
    )
    if (is.null(mv)) return(NULL)
    res <- tryCatch(MVMR::ivw_mvmr(mv),
                    error = function(e) { cat("  ivw_mvmr err:", e$message, "\n"); NULL })
    if (is.null(res)) return(NULL)
    cond_F <- tryCatch(MVMR::strength_mvmr(mv, gencov = 0),
                       error = function(e) NULL)
    Q <- tryCatch(MVMR::pleiotropy_mvmr(mv, gencov = 0),
                  error = function(e) NULL)
  } else {
    # Manual MVMR-IVW (Sanderson 2019, eq 1; weight = 1/seBYG^2)
    # beta_Y = sum_k beta_X_k * theta_k + epsilon
    # Solve via WLS: theta = (X' W X)^-1 X' W y
    W <- diag(1 / seBYG^2)
    XtWX <- t(BXGs) %*% W %*% BXGs
    XtWy <- t(BXGs) %*% W %*% BYG
    theta <- tryCatch(solve(XtWX, XtWy),
                      error = function(e) { cat("  solve err:", e$message, "\n"); NULL })
    if (is.null(theta)) return(NULL)
    # Asymptotic SE: sqrt(diag((X' W X)^-1)) under unit dispersion
    Vtheta <- solve(XtWX)
    se_theta <- sqrt(diag(Vtheta))
    z <- as.numeric(theta) / se_theta
    p_theta <- 2 * pnorm(-abs(z))
    res <- data.frame(Estimate = as.numeric(theta), `Std. Error` = se_theta,
                      `t value` = z, `Pr(>|t|)` = p_theta,
                      check.names = FALSE)
    # Conditional F: per Sanderson 2019, F_k = (n_SNP - K) * delta_k / Q_k
    # Simple approximation: F_k = mean(BXG_k^2 / seBXG_k^2) / (1 + sum_j!=k of cov)
    cond_F <- numeric(ncol(BXGs))
    for (k in seq_len(ncol(BXGs))) {
      # Approximate cond_F as residual F after regressing X_k on other X's
      others <- BXGs[, -k, drop = FALSE]
      r_k <- BXGs[, k] - others %*% solve(t(others) %*% others) %*% (t(others) %*% BXGs[, k])
      cond_F[k] <- mean(r_k^2 / seBXGs[, k]^2)
    }
    Q <- NULL  # not computed in fallback
  }

  # Tabulate
  res <- as.data.frame(res)
  # Standard MVMR output has columns: Estimate, Std. Error, t value, Pr(>|t|)
  est_col <- intersect(c("Estimate","b","est"), names(res))[1]
  se_col  <- intersect(c("Std. Error","se"), names(res))[1]
  p_col   <- intersect(c("Pr(>|t|)","pval","P"), names(res))[1]
  if (is.na(est_col) || is.na(se_col) || is.na(p_col)) {
    cat("  unexpected ivw_mvmr columns: ", paste(names(res), collapse=","), "\n")
    return(NULL)
  }
  cond_F_vec <- if (!is.null(cond_F)) as.numeric(cond_F) else rep(NA, 4)
  if (length(cond_F_vec) < 4) cond_F_vec <- c(cond_F_vec, rep(NA, 4 - length(cond_F_vec)))

  data.frame(
    antigen = ant_lab, outcome = out_lab, n_snp = length(common),
    OR_mvmr   = exp(res[[est_col]][1]),
    L95_mvmr  = exp(res[[est_col]][1] - 1.96 * res[[se_col]][1]),
    U95_mvmr  = exp(res[[est_col]][1] + 1.96 * res[[se_col]][1]),
    P_mvmr    = res[[p_col]][1],
    F_cond_antigen   = cond_F_vec[1],
    F_cond_BMI       = cond_F_vec[2],
    F_cond_smoking   = cond_F_vec[3],
    F_cond_education = cond_F_vec[4],
    Q_stat = if (!is.null(Q)) Q$Qstat[1] else NA,
    Q_pval = if (!is.null(Q)) Q$Qpval[1] else NA,
    valid  = ifelse(!is.na(cond_F_vec[1]) && cond_F_vec[1] >= 10,
                    "valid (F_cond>=10)", "WEAK (F_cond<10)"),
    method = ifelse(HAS_MVMR, "MVMR-IVW (Sanderson 2019)",
                     "manual MVMR-IVW (WLS implementation)"),
    row.names = NULL
  )
}

# Map exposure labels in summary table to our antigen labels
norm_lab <- function(x) {
  x <- as.character(x)
  x[x == "hpylori_general"] <- "hpylori_general"
  x[x == "cagA"]             <- "cagA"
  x[x == "catalase"]         <- "catalase"
  x[x == "groEL"]            <- "groEL"
  x[x == "omp"]              <- "omp"
  x[x == "ureA"]             <- "ureA"
  x[x == "vacA"]             <- "vacA"
  x
}

results <- list()
for (i in seq_len(nrow(cand))) {
  ant_lab <- norm_lab(cand$exposure[i])
  out_lab <- cand$outcome[i]
  res <- tryCatch(run_proper_mvmr(ant_lab, out_lab),
                  error = function(e) { cat("  TOP err:", e$message, "\n"); NULL })
  if (!is.null(res)) results[[paste(ant_lab, out_lab, sep="__")]] <- res
}

if (length(results)) {
  out_tbl <- do.call(rbind, results)
  # Add univariate comparison
  uni <- st[, .(exposure, outcome, OR_uni = OR, L95_uni = L95, U95_uni = U95,
                 P_uni = P_ivw)]
  out_tbl <- merge(out_tbl,
                   uni,
                   by.x = c("antigen","outcome"), by.y = c("exposure","outcome"),
                   all.x = TRUE)
  out_tbl <- out_tbl[order(out_tbl$P_uni), ]
  fwrite(out_tbl, "results/mvmr_proper_results.tsv", sep = "\t")
  cat("\n=== Proper MVMR summary ===\n")
  print(out_tbl[, c("antigen","outcome","n_snp",
                     "OR_uni","P_uni","OR_mvmr","P_mvmr",
                     "F_cond_antigen","valid")])
  cat("\nWrote results/mvmr_proper_results.tsv\n")

  # Rename old to deprecated
  if (file.exists("results/mvmr_results.tsv") &&
      !file.exists("results/mvmr_results_DEPRECATED.tsv")) {
    file.rename("results/mvmr_results.tsv",
                "results/mvmr_results_DEPRECATED.tsv")
    cat("Renamed old mvmr_results.tsv -> mvmr_results_DEPRECATED.tsv\n")
  }
} else {
  cat("\n!!! No MVMR results produced\n")
}
cat("\nDone.\n")
