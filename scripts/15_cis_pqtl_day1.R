#!/usr/bin/env Rscript
# 15_cis_pqtl_day1.R
#
# Round 2 / TASK 1 / Day 1:
# Locate and extract cis-pQTL instruments for TLR1 (chr4) and FCGR2A (chr1)
# from OpenGWAS-indexed pQTL datasets (UKB-PPP Sun 2023, deCODE Ferkingstad
# 2021). Report per-protein n_SNP and F-statistic to enable scenario A/B/C
# decision before committing to full cis-pQTL MR (Day 2-5).
#
# Decision thresholds (per Round 2 task list 1.4):
#   PASS  = n_SNP >= 5 AND mean F-statistic >= 50  -> scenario A possible
#   MARGINAL = n_SNP >= 3 AND mean F-statistic >= 30 -> scenario B
#   FAIL  = n_SNP < 3 OR mean F-statistic < 30  -> scenario C (abandon cis-pQTL)
#
# Output:
#   data/exposure/cis_pqtl_day1_candidates.tsv  -- OpenGWAS IDs found
#   data/exposure/cis_pqtl_instruments.rds      -- per-protein SNP table (if PASS)
#   logs/cis_pqtl_day1.log
#   logs/cis_pqtl_quality_decision.txt

set.seed(20260521)

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(dplyr)
})

dir.create("data/exposure", showWarnings = FALSE, recursive = TRUE)
dir.create("logs",          showWarnings = FALSE, recursive = TRUE)
sink("logs/cis_pqtl_day1.log", split = TRUE)
on.exit(sink(NULL), add = TRUE)

cat("=== 15_cis_pqtl_day1.R  ", format(Sys.time()), " ===\n")
cat("PURPOSE: Round 2 / Day 1 cis-pQTL instrument quality check\n")
cat("TARGETS: TLR1 (chr4), FCGR2A (chr1)\n\n")

# Gene coordinates (GRCh37/hg19)
GENES <- data.frame(
  protein = c("TLR1", "FCGR2A"),
  chr     = c(4, 1),
  start   = c(38792156, 161475219),
  end     = c(38805749, 161498664),
  stringsAsFactors = FALSE
)
WINDOW_MB <- 1
P_THRESH  <- 5e-8

# ---- Search OpenGWAS for protein GWAS -------------------------------------
cat("Querying available_outcomes() for protein candidates...\n")
ao <- available_outcomes()
cat("Total OpenGWAS records:", nrow(ao), "\n")

find_protein <- function(prot_name, alt_names = character(0)) {
  patterns <- c(paste0("\\b", prot_name, "\\b"),
                paste0("^", prot_name, "$"),
                paste0("^", prot_name, " "),
                alt_names)
  patterns <- patterns[nzchar(patterns)]
  rgx <- paste(patterns, collapse = "|")
  hits <- ao[grepl(rgx, ao$trait, ignore.case = TRUE) |
              grepl(rgx, ao$id,    ignore.case = TRUE), ]
  # Prefer pQTL studies (id starts with prot-, large N, recent)
  hits$is_pqtl <- grepl("^prot-|^pqtl|protein|Olink|SomaScan",
                         paste(hits$id, hits$trait, hits$consortium),
                         ignore.case = TRUE)
  hits <- hits[order(-hits$is_pqtl, -hits$sample_size), ]
  hits
}

tlr1_hits <- find_protein("TLR1", c("TLR-1","toll.like.receptor.1"))
cat("\nTLR1 candidate records (top 10):\n")
print(head(tlr1_hits[, c("id","trait","sample_size","consortium","year")], 10))

fcgr2a_hits <- find_protein("FCGR2A", c("CD32","FCGR2","Fc gamma receptor IIA"))
cat("\nFCGR2A candidate records (top 10):\n")
print(head(fcgr2a_hits[, c("id","trait","sample_size","consortium","year")], 10))

# Save full candidate list
all_hits <- rbind(
  cbind(target = "TLR1",   tlr1_hits[, c("id","trait","sample_size","consortium","year","is_pqtl")]),
  cbind(target = "FCGR2A", fcgr2a_hits[, c("id","trait","sample_size","consortium","year","is_pqtl")])
)
fwrite(all_hits, "data/exposure/cis_pqtl_day1_candidates.tsv", sep = "\t")

# ---- Extract cis-pQTL for top candidate ------------------------------------
extract_cis <- function(prot_name, chr, gene_start, gene_end, candidates_df) {
  cat("\n--- Extracting cis-pQTL for", prot_name, "---\n")
  if (nrow(candidates_df) == 0) {
    cat("  no OpenGWAS hits for", prot_name, "\n"); return(NULL)
  }
  # Try candidates in priority order; stop at first that yields usable instrument
  win <- WINDOW_MB * 1e6
  cis_lo <- max(0, gene_start - win)
  cis_hi <- gene_end + win
  cat("  cis-region: chr", chr, ":", cis_lo, "-", cis_hi, "\n", sep = "")

  for (i in seq_len(min(8, nrow(candidates_df)))) {
    id <- candidates_df$id[i]
    cat("\n  trying id:", id, " (", candidates_df$trait[i], ", N=",
        candidates_df$sample_size[i], ")\n", sep = "")
    d <- tryCatch(
      extract_instruments(outcomes = id, p1 = P_THRESH,
                          clump = TRUE, r2 = 0.001, kb = 10000),
      error = function(e) { cat("    err:", e$message, "\n"); NULL }
    )
    if (is.null(d) || nrow(d) == 0) { cat("    no SNPs\n"); next }
    cat("    total SNPs at p<5e-8:", nrow(d), "\n")
    if (!("chr.exposure" %in% names(d)) || !("pos.exposure" %in% names(d))) {
      cat("    missing chr/pos columns\n"); next
    }
    cis <- d[d$chr.exposure == as.character(chr) &
              d$pos.exposure >= cis_lo & d$pos.exposure <= cis_hi, ]
    cat("    cis SNPs (chr", chr, ":", cis_lo, "-", cis_hi, "): ",
        nrow(cis), "\n", sep = "")
    if (nrow(cis) == 0) next
    cis$F_stat <- (cis$beta.exposure / cis$se.exposure)^2
    cis$protein <- prot_name
    cat("    mean F = ", round(mean(cis$F_stat), 1),
        ", min F = ", round(min(cis$F_stat), 1),
        ", max F = ", round(max(cis$F_stat), 1), "\n", sep = "")
    return(cis)
  }
  cat("  no candidate yielded cis-SNPs at p<5e-8\n")
  NULL
}

tlr1_iv <- extract_cis("TLR1", GENES$chr[1], GENES$start[1], GENES$end[1], tlr1_hits)
Sys.sleep(2)
fcgr2a_iv <- extract_cis("FCGR2A", GENES$chr[2], GENES$start[2], GENES$end[2], fcgr2a_hits)

instruments <- list()
if (!is.null(tlr1_iv))   instruments$TLR1   <- tlr1_iv
if (!is.null(fcgr2a_iv)) instruments$FCGR2A <- fcgr2a_iv

if (length(instruments)) {
  saveRDS(instruments, "data/exposure/cis_pqtl_instruments.rds")
  cat("\nSaved data/exposure/cis_pqtl_instruments.rds\n")
}

# ---- Quality decision ------------------------------------------------------
score <- function(iv) {
  if (is.null(iv) || nrow(iv) == 0) return(list(n=0, mF=NA, status="FAIL"))
  n  <- nrow(iv)
  mF <- mean(iv$F_stat)
  status <- if (n >= 5 && mF >= 50) {
    "PASS (scenario A)"
  } else if (n >= 3 && mF >= 30) {
    "MARGINAL (scenario B)"
  } else {
    "FAIL (scenario C)"
  }
  list(n = n, mF = mF, status = status)
}
s_tlr1   <- score(tlr1_iv)
s_fcgr2a <- score(fcgr2a_iv)

dec <- c(
  "=== Round 2 / Day 1 cis-pQTL instrument quality ===",
  paste0("Timestamp: ", format(Sys.time())),
  "",
  sprintf("TLR1   : n_SNP = %d  meanF = %s  -> %s",
          s_tlr1$n,
          if (is.na(s_tlr1$mF)) "NA" else sprintf("%.1f", s_tlr1$mF),
          s_tlr1$status),
  sprintf("FCGR2A : n_SNP = %d  meanF = %s  -> %s",
          s_fcgr2a$n,
          if (is.na(s_fcgr2a$mF)) "NA" else sprintf("%.1f", s_fcgr2a$mF),
          s_fcgr2a$status),
  "",
  "PASS thresholds: n>=5 AND meanF>=50",
  "MARGINAL:        n>=3 AND meanF>=30",
  "FAIL:            either condition violated",
  "",
  "Combined verdict:"
)
verdict <- if (s_tlr1$status == "PASS (scenario A)" && s_fcgr2a$status == "PASS (scenario A)") {
  "PASS -- proceed to scenario A (Gut framing with three converging approaches)"
} else if (grepl("FAIL", s_tlr1$status) && grepl("FAIL", s_fcgr2a$status)) {
  "FAIL -- both proteins inadequate. Recommend scenario C (abandon cis-pQTL; submit v3 to eBioMedicine)"
} else if (grepl("PASS|MARGINAL", s_tlr1$status) || grepl("PASS|MARGINAL", s_fcgr2a$status)) {
  "MIXED -- at least one protein adequate. Proceed to scenario B (Gut attempt, eBioMedicine backup)"
} else {
  "UNCLEAR"
}
dec <- c(dec, paste("  ", verdict))
writeLines(dec, "logs/cis_pqtl_quality_decision.txt")
cat("\n", paste(dec, collapse = "\n"), "\n", sep = "")

# ---- If no IEU pQTL found, document UKB-PPP direct-download path -----------
if (length(instruments) == 0) {
  cat("\n!!! NO usable IEU pQTL records found for either protein.\n")
  cat("    Direct download required from UKB-PPP:\n")
  cat("    https://metabolomips.org/ukbppv1/\n")
  cat("    Look for TLR1 (chr4) and FCGR2A (chr1) summary statistics files.\n")
  cat("    deCODE pQTL fallback: https://www.decode.com/summarydata/\n")
}
cat("\nDone.\n")
