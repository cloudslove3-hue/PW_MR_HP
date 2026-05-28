#!/usr/bin/env Rscript
# 02_outcome.R  Tier 1 outcome extraction
#
# Inputs:
#   data/exposure/exposure_list.rds
#   config/outcomes_tier1.tsv   (researcher-curated; see template)
#
# Output:
#   data/outcome/outcome_data_list.rds
#   data/outcome/outcome_extraction_log.tsv
#
# config/outcomes_tier1.tsv must have at minimum these columns:
#   outcome_label  id  source  path  consortium  n  ncase  year  notes
#
# 'source' is one of: opengwas | finngen | local_file | gwas_catalog
#   opengwas:      id = OpenGWAS ID (e.g. "ieu-b-7"), path = ""
#   finngen:       id = FinnGen phenocode (e.g. "K11_GASTRITDUOD"), path = local .gz
#   local_file:    id = anything informative, path = local file (any format)
#   gwas_catalog:  id = GCST..., path = local file
#
# Discovery helper at the top: prints candidates for keywords from the md spec.
# Run this script once with NO outcomes_tier1.tsv to view candidates,
# then create the config file and re-run.

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
})

dir.create("data/outcome", showWarnings = FALSE, recursive = TRUE)
dir.create("config", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

log_path <- file.path("logs", "02_outcome.log")
sink(log_path, split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 02_outcome.R  ", format(Sys.time()), " ===\n")

# ---- Load config (skip # comments) -----------------------------------------
cfg_path <- "config/outcomes_tier1.tsv"
read_cfg <- function(p) {
  if (!file.exists(p)) return(NULL)
  raw <- readLines(p, warn = FALSE)
  raw <- raw[!grepl("^\\s*#", raw) & nzchar(trimws(raw))]
  if (length(raw) < 2) return(data.frame())  # header only or empty
  tryCatch(fread(text = paste(raw, collapse = "\n"), sep = "\t"),
           error = function(e) data.frame())
}
cfg <- read_cfg(cfg_path)

# ---- Discovery mode ---------------------------------------------------------
if (is.null(cfg) || nrow(cfg) == 0) {
  cat("\nNo populated config/outcomes_tier1.tsv yet. Printing candidates.\n")
  cat("Copy the rows you want into the template at config/outcomes_tier1.tsv.\n\n")
  ao <- available_outcomes()
  find_outcome <- function(keyword, k = 8) {
    h <- ao[grepl(keyword, ao$trait, ignore.case = TRUE), ]
    if (nrow(h) == 0) return(invisible(NULL))
    h <- h[order(-h$sample_size), ]
    cat("\n### ", keyword, " ###\n", sep = "")
    print(head(h[, c("id","trait","sample_size","ncase","consortium","year","population")], k))
  }
  for (kw in c("Parkinson","Alzheimer","coronary","ischemic stroke",
               "atrial fibrillation","ulcerative colitis","Crohn","rheumatoid",
               "lupus","multiple sclerosis","type 2 diabetes","pancreatic cancer",
               "colorectal cancer","iron deficiency","thrombocytopenic",
               "asthma","atopic")) {
    find_outcome(kw)
  }
  cat("\n!!! Stopping. Populate config/outcomes_tier1.tsv and re-run.\n")
  quit(status = 1)
}

# ---- Load instruments + config ----------------------------------------------
exposure_list <- readRDS("data/exposure/exposure_list.rds")
all_snps <- unique(unlist(lapply(exposure_list, function(d) d$SNP)))
cat("Total unique exposure SNPs:", length(all_snps), "\n")

req <- c("outcome_label","id","source")
stopifnot(all(req %in% names(cfg)))

# ---- Outcome extractors -----------------------------------------------------
get_opengwas <- function(id, snps) {
  tryCatch(
    extract_outcome_data(snps = snps, outcomes = id, proxies = TRUE),
    error = function(e) { cat("  opengwas fail:", id, "-", e$message, "\n"); NULL }
  )
}

# FinnGen column convention (R10+): rsids, beta, sebeta, pval, af_alt, alt, ref,
# n_case (per-pheno), n_control
get_finngen <- function(path, id, snps) {
  if (!file.exists(path)) {
    cat("  finngen file missing:", path, "\n"); return(NULL)
  }
  d <- fread(path)
  cn <- tolower(names(d)); names(d) <- cn
  pick <- function(...) { for (x in c(...)) if (x %in% cn) return(x); NA_character_ }
  c_snp  <- pick("rsids","rsid","snp")
  c_beta <- pick("beta")
  c_se   <- pick("sebeta","se","standard_error")
  c_p    <- pick("pval","p_value","pvalue")
  c_ea   <- pick("alt","effect_allele")
  c_oa   <- pick("ref","other_allele")
  c_eaf  <- pick("af_alt","eaf")
  if (any(is.na(c(c_snp,c_beta,c_se,c_p,c_ea,c_oa)))) {
    cat("  finngen header unexpected:", paste(cn, collapse=","), "\n"); return(NULL)
  }
  d <- d[get(c_snp) %in% snps]
  if (nrow(d) == 0) return(NULL)
  format_data(
    as.data.frame(d), type = "outcome", snps = snps,
    snp_col = c_snp, beta_col = c_beta, se_col = c_se,
    effect_allele_col = c_ea, other_allele_col = c_oa,
    eaf_col = if (!is.na(c_eaf)) c_eaf else NULL,
    pval_col = c_p
  )
}

# Generic local file: needs columns SNP, beta, se, ea, oa, pval — adjust with format_data
get_local <- function(path, id, snps) {
  if (!file.exists(path)) return(NULL)
  d <- fread(path); names(d) <- tolower(names(d))
  pick <- function(...) { for (x in c(...)) if (x %in% names(d)) return(x); NA_character_ }
  c_snp  <- pick("snp","rsid","rsids","variant_id")
  c_beta <- pick("beta","effect")
  c_se   <- pick("se","standard_error","sebeta")
  c_p    <- pick("pval","p_value","pvalue")
  c_ea   <- pick("effect_allele","alt","ea")
  c_oa   <- pick("other_allele","ref","nea")
  c_eaf  <- pick("eaf","effect_allele_frequency","af_alt")
  if (any(is.na(c(c_snp,c_beta,c_se,c_p,c_ea,c_oa)))) return(NULL)
  d <- d[get(c_snp) %in% snps]
  if (nrow(d) == 0) return(NULL)
  format_data(
    as.data.frame(d), type = "outcome", snps = snps,
    snp_col = c_snp, beta_col = c_beta, se_col = c_se,
    effect_allele_col = c_ea, other_allele_col = c_oa,
    eaf_col = if (!is.na(c_eaf)) c_eaf else NULL,
    pval_col = c_p
  )
}

# ---- Main loop --------------------------------------------------------------
outcome_data_list <- list()
extraction_log <- data.frame()

for (i in seq_len(nrow(cfg))) {
  lab <- cfg$outcome_label[i]; id <- cfg$id[i]; src <- cfg$source[i]
  pth <- if ("path" %in% names(cfg)) cfg$path[i] else ""
  cat("\n[", i, "/", nrow(cfg), "] ", lab, " <- ", src, " : ", id, "\n", sep = "")

  od <- switch(src,
    opengwas     = get_opengwas(id, all_snps),
    finngen      = get_finngen(pth, id, all_snps),
    local_file   = get_local(pth, id, all_snps),
    gwas_catalog = get_local(pth, id, all_snps),
    { cat("  unknown source:", src, "\n"); NULL }
  )

  if (!is.null(od)) {
    od$outcome <- lab; od$id.outcome <- id
    outcome_data_list[[lab]] <- od
    cat("  -> ", nrow(od), " SNPs harvested\n")
  } else {
    cat("  -> NO data\n")
  }
  extraction_log <- rbind(extraction_log,
    data.frame(outcome_label = lab, id = id, source = src,
               n_snp = if (is.null(od)) 0 else nrow(od)))
  if (src == "opengwas") Sys.sleep(1)  # rate-limit
}

saveRDS(outcome_data_list, "data/outcome/outcome_data_list.rds")
fwrite(extraction_log, "data/outcome/outcome_extraction_log.tsv", sep = "\t")
cat("\n=== Summary ===\n"); print(extraction_log)
cat("\nDone. Proceed to scripts/03_mr.R\n")
