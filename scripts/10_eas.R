#!/usr/bin/env Rscript
# 10_eas.R  East Asian cross-ancestry replication
#
# CRITICAL LIMITATION (must be in paper Discussion):
#   H. pylori antigen serology GWAS exists only in Europeans (Butler-Laporte 2020).
#   We cannot run native EAS MR. Instead we use EUR instruments and EAS outcome
#   GWAS — this is suboptimal because EUR-derived instruments may have different
#   allele frequencies and LD structure in EAS, weakening power. Effects should
#   be directionally consistent if EUR signals are causal.
#
# Strategy:
#   1. For each Tier 1 outcome with an EAS counterpart on OpenGWAS (mostly Biobank Japan):
#   2. Extract outcome data for our 75 EUR exposure SNPs (proxies allowed)
#   3. Run IVW (and weighted median for robustness)
#   4. Compare effect direction and magnitude with EUR primary table
#
# Inputs:   data/exposure/exposure_list.rds
#           results/MR_summary_table.tsv (EUR primary)
# Outputs:  results/eas_replication.tsv
#           logs/10_eas.log

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(dplyr)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("logs",    showWarnings = FALSE, recursive = TRUE)
sink("logs/10_eas.log", split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 10_eas.R  ", format(Sys.time()), " ===\n")

exposure_list <- readRDS("data/exposure/exposure_list.rds")
all_snps <- unique(unlist(lapply(exposure_list, function(d) d$SNP)))
eur <- fread("results/MR_summary_table.tsv")

# ---- Identify EAS outcome candidates ---------------------------------------
# Pull from OpenGWAS: bbj-* IDs are Biobank Japan
ao <- available_outcomes()
eas_pool <- ao[grepl("^bbj-", ao$id) | (!is.na(ao$population) & ao$population == "East Asian"), ]
cat("EAS pool size:", nrow(eas_pool), "\n")

# Keywords matched to our Tier 1 outcomes
keywords <- list(
  parkinson   = "[Pp]arkinson",
  alzheimer   = "[Aa]lzheimer",
  cad         = "[Cc]oronary artery|[Cc]ardiac|[Mm]yocardial infarction",
  stroke_isch = "[Ii]schemic stroke|[Cc]erebral infarction",
  afib        = "[Aa]trial fibrillation",
  uc          = "[Uu]lcerative colitis",
  crohn       = "[Cc]rohn",
  ra          = "[Rr]heumatoid arthritis",
  sle         = "[Ss]ystemic lupus|[Ll]upus",
  ms          = "[Mm]ultiple sclerosis",
  t2d         = "[Tt]ype 2 diabetes",
  pancan      = "[Pp]ancreatic cancer|[Pp]ancreatic",
  crc         = "[Cc]olorectal cancer|[Cc]olorectal",
  ida         = "[Ii]ron deficiency|[Aa]nemia",
  itp         = "[Tt]hrombocytopenic",
  asthma      = "[Aa]sthma",
  atopic      = "[Aa]topic"
)

pick_eas <- function(label, regex) {
  hits <- eas_pool[grepl(regex, eas_pool$trait), ]
  if (!nrow(hits)) return(NULL)
  hits <- hits[order(-hits$sample_size), ]
  hits$tier1_label <- label
  head(hits, 1)
}

eas_targets <- bind_rows(lapply(names(keywords), function(k) pick_eas(k, keywords[[k]])))
cat("\nEAS targets identified:", nrow(eas_targets), "\n")
print(eas_targets[, c("tier1_label","id","trait","sample_size","ncase","population")])

# ---- Extract EAS outcome data (one ID at a time to avoid type-conflict) ----
cat("\nExtracting EAS outcomes for", length(all_snps), "SNPs...\n")
eas_chunks <- list()
for (oid in eas_targets$id) {
  cat("  ", oid, "\n")
  d <- tryCatch(
    extract_outcome_data(snps = all_snps, outcomes = oid, proxies = TRUE),
    error = function(e) { cat("    err:", e$message, "\n"); NULL }
  )
  if (!is.null(d) && nrow(d) > 0) {
    # Coerce numeric cols that sometimes come as character
    for (cn in c("n","samplesize.outcome","ncase.outcome","ncontrol.outcome")) {
      if (cn %in% names(d)) d[[cn]] <- suppressWarnings(as.numeric(as.character(d[[cn]])))
    }
    eas_chunks[[oid]] <- d
  }
  Sys.sleep(1)
}
eas_out <- if (length(eas_chunks)) rbindlist(eas_chunks, fill = TRUE) else NULL
if (is.null(eas_out) || nrow(eas_out) == 0) {
  cat("No EAS data fetched — aborting\n"); quit(status = 0)
}
cat("EAS rows fetched:", nrow(eas_out), "\n")

id2label <- setNames(eas_targets$tier1_label, eas_targets$id)

# ---- Per-pair MR (IVW + weighted median) ----------------------------------
rows <- list()
for (ant_lab in names(exposure_list)) {
  exp_df <- exposure_list[[ant_lab]]
  if (is.null(exp_df) || nrow(exp_df) == 0) next
  for (oid in eas_targets$id) {
    out_df <- eas_out[eas_out$id.outcome == oid, ]
    if (nrow(out_df) < 3) next
    dat <- tryCatch(harmonise_data(exp_df, out_df, action = 2),
                    error = function(e) NULL)
    if (is.null(dat)) next
    dat <- subset(dat, mr_keep == TRUE)
    if (nrow(dat) < 3) next
    res <- tryCatch(
      mr(dat, method_list = c("mr_ivw","mr_weighted_median")),
      error = function(e) NULL
    )
    if (is.null(res) || !nrow(res)) next
    ivw <- res[res$method == "Inverse variance weighted", ]
    wmd <- res[res$method == "Weighted median", ]
    rows[[length(rows) + 1]] <- data.frame(
      exposure   = ant_lab,
      outcome    = id2label[oid],
      eas_id     = oid,
      n_snp      = nrow(dat),
      OR_eas     = if (nrow(ivw)) exp(ivw$b[1]) else NA,
      L95_eas    = if (nrow(ivw)) exp(ivw$b[1] - 1.96*ivw$se[1]) else NA,
      U95_eas    = if (nrow(ivw)) exp(ivw$b[1] + 1.96*ivw$se[1]) else NA,
      P_eas      = if (nrow(ivw)) ivw$pval[1] else NA,
      OR_wmd_eas = if (nrow(wmd)) exp(wmd$b[1]) else NA,
      P_wmd_eas  = if (nrow(wmd)) wmd$pval[1] else NA
    )
  }
}

eas_tbl <- bind_rows(rows)
# Join EUR primary
joined <- merge(
  eas_tbl,
  eur[, c("exposure","outcome","OR","L95","U95","P_ivw")],
  by = c("exposure","outcome"), all.x = TRUE, suffixes = c("","_eur")
)
setnames(joined, "OR", "OR_eur"); setnames(joined, "L95","L95_eur")
setnames(joined, "U95", "U95_eur"); setnames(joined, "P_ivw","P_eur")
joined$dir_concordant <- with(joined,
  !is.na(OR_eur) & !is.na(OR_eas) & sign(log(OR_eur)) == sign(log(OR_eas)))

setorder(joined, P_eas)
fwrite(joined, "results/eas_replication.tsv", sep = "\t")
cat("\nWrote results/eas_replication.tsv  (", nrow(joined), "rows )\n")

cat("\n=== Direction concordance EUR vs EAS ===\n")
print(table(joined$dir_concordant, useNA = "ifany"))

cat("\n=== Top 10 by EAS P ===\n")
print(head(joined[, c("exposure","outcome","n_snp",
                       "OR_eur","P_eur","OR_eas","P_eas","dir_concordant")], 10))

cat("\nDone.\n")
