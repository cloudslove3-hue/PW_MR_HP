#!/usr/bin/env Rscript
# 09_finngen.R  FinnGen phenome-wide Tier 2 scan
#
# Strategy:
#   1. Enumerate all FinnGen IDs in OpenGWAS (finn-b-*, possibly finn-r10-*)
#   2. Filter by ncase >= 200 and dedupe by trait (keep highest ncase release)
#   3. Extract outcome data for all 75 antigen SNPs against this panel
#   4. Run IVW (only) per antigen × outcome (fast, no PRESSO at this scale)
#   5. Bonferroni + FDR over the FULL phenome-wide set
#
# Inputs:   data/exposure/exposure_list.rds
# Outputs:
#   data/outcome/finngen_outcome_data.rds
#   results/finngen_phenome_summary.tsv
#   results/finngen_phenome_signals_FDR.tsv
#   logs/09_finngen.log

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(dplyr)
})

dir.create("data/outcome", showWarnings = FALSE, recursive = TRUE)
dir.create("results",      showWarnings = FALSE, recursive = TRUE)
dir.create("logs",         showWarnings = FALSE, recursive = TRUE)
sink("logs/09_finngen.log", split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 09_finngen.R  ", format(Sys.time()), " ===\n")

exposure_list <- readRDS("data/exposure/exposure_list.rds")
all_snps <- unique(unlist(lapply(exposure_list, function(d) d$SNP)))
cat("Unique exposure SNPs:", length(all_snps), "\n")

# ---- Enumerate FinnGen IDs --------------------------------------------------
ao <- available_outcomes()
fg <- ao[grepl("^finn-", ao$id), ]
cat("Total FinnGen IDs on OpenGWAS:", nrow(fg), "\n")
if (nrow(fg) == 0) {
  cat("No FinnGen IDs available — skipping\n")
  quit(status = 0)
}

# Filter by ncase
fg <- fg[!is.na(fg$ncase) & fg$ncase >= 200, ]
cat("After ncase>=200 filter:", nrow(fg), "\n")

# Dedupe by trait — keep largest ncase per trait
fg <- fg[order(fg$trait, -fg$ncase), ]
fg <- fg[!duplicated(fg$trait), ]
cat("After dedupe by trait:", nrow(fg), "\n")

# Order by ncase desc, take top
N_OUTCOMES <- min(500, nrow(fg))
fg <- head(fg[order(-fg$ncase), ], N_OUTCOMES)
cat("Will scan top", N_OUTCOMES, "FinnGen outcomes (ncase range:",
    min(fg$ncase), "-", max(fg$ncase), ")\n")

fwrite(fg[, c("id","trait","sample_size","ncase","year","consortium")],
       "config/finngen_outcomes_scanned.tsv", sep = "\t")

# ---- Extract outcomes in chunks --------------------------------------------
chunk_extract <- function(snps, ids, chunk = 50, sleep = 1) {
  out <- list()
  for (i in seq(1, length(ids), by = chunk)) {
    sub <- ids[i:min(i + chunk - 1, length(ids))]
    cat("  chunk", i, "-", i + length(sub) - 1, "/", length(ids), "\n")
    d <- tryCatch(extract_outcome_data(snps = snps, outcomes = sub, proxies = TRUE),
                  error = function(e) { cat("  chunk err:", e$message, "\n"); NULL })
    if (!is.null(d)) out[[length(out) + 1]] <- d
    Sys.sleep(sleep)
  }
  rbindlist(out, fill = TRUE)
}

cat("\nExtracting outcomes...\n")
fg_outcomes <- chunk_extract(all_snps, fg$id, chunk = 30, sleep = 1)
cat("Total outcome rows:", nrow(fg_outcomes), "\n")
saveRDS(fg_outcomes, "data/outcome/finngen_outcome_data.rds")

# Map id -> trait for labelling
id2trait <- setNames(fg$trait, fg$id)

# ---- Loop pairs and run IVW -----------------------------------------------
all_exposures <- do.call(rbind, lapply(exposure_list, function(d) {
  d$exposure_label <- d$exposure_label[1] %||% NA
  d
}))
# Helper to avoid null
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 0)) b else a

rows <- list()
ids <- unique(fg_outcomes$id.outcome)
cat("\nIVW MR per antigen × outcome (", length(ids), "outcomes ) ...\n")

for (ant_lab in names(exposure_list)) {
  exp_df <- exposure_list[[ant_lab]]
  if (is.null(exp_df) || nrow(exp_df) == 0) next
  for (oid in ids) {
    out_df <- fg_outcomes[fg_outcomes$id.outcome == oid, ]
    if (nrow(out_df) < 3) next
    dat <- tryCatch(harmonise_data(exp_df, out_df, action = 2),
                    error = function(e) NULL)
    if (is.null(dat)) next
    dat <- subset(dat, mr_keep == TRUE)
    if (nrow(dat) < 3) next
    res <- tryCatch(mr(dat, method_list = "mr_ivw"),
                    error = function(e) NULL)
    if (is.null(res) || !nrow(res)) next
    rows[[length(rows) + 1]] <- data.frame(
      exposure_label = ant_lab,
      outcome_id     = oid,
      outcome_trait  = id2trait[oid],
      n_snp = nrow(dat),
      b     = res$b[1], se = res$se[1], pval = res$pval[1]
    )
  }
}

summary_tbl <- bind_rows(rows)
summary_tbl$OR  <- exp(summary_tbl$b)
summary_tbl$L95 <- exp(summary_tbl$b - 1.96 * summary_tbl$se)
summary_tbl$U95 <- exp(summary_tbl$b + 1.96 * summary_tbl$se)
summary_tbl$P_bonf <- p.adjust(summary_tbl$pval, method = "bonferroni")
summary_tbl$P_fdr  <- p.adjust(summary_tbl$pval, method = "fdr")
summary_tbl <- summary_tbl[order(summary_tbl$pval), ]

fwrite(summary_tbl, "results/finngen_phenome_summary.tsv", sep = "\t")
cat("\nWrote results/finngen_phenome_summary.tsv  (", nrow(summary_tbl), "rows )\n")

sig <- subset(summary_tbl, !is.na(P_fdr) & P_fdr < 0.05)
fwrite(sig, "results/finngen_phenome_signals_FDR.tsv", sep = "\t")
cat("FDR<0.05 signals:", nrow(sig), "\n")

cat("\n=== Top 20 by raw P ===\n")
print(head(summary_tbl[, c("exposure_label","outcome_trait","n_snp",
                            "OR","L95","U95","pval","P_bonf","P_fdr")], 20))

if (nrow(sig)) {
  cat("\n=== FDR<0.05 signals ===\n")
  print(sig[, c("exposure_label","outcome_trait","OR","L95","U95","pval","P_fdr")])
}
cat("\nDone.\n")
