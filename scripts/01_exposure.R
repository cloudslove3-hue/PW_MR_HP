#!/usr/bin/env Rscript
# 01_exposure.R  H. pylori 7-antigen instrument extraction
#
# Output:
#   data/exposure/<antigen>_instruments.tsv  (per antigen)
#   data/exposure/exposure_list.rds          (named list)
#   data/exposure/source_log.tsv             (which path was used per antigen)
#
# Strategy:
#   1. Try OpenGWAS via TwoSampleMR::extract_instruments()
#   2. If ID missing or returns <3 SNPs at p=5e-8, retry at p=5e-6
#   3. If still empty: fall back to direct download from EBI GWAS Catalog FTP
#      (Butler-Laporte 2020 Nat Commun, GCST90006910-90006916), then locally
#      clump via ieugwasr::ld_clump() against the EUR reference.
#
# Butler-Laporte serology phenotype map (GCST IDs are PUBLISHED, do not change):
#   GCST90006910  H. pylori general (any antigen positive)
#   GCST90006911  CagA seropositivity
#   GCST90006912  Catalase
#   GCST90006913  GroEL
#   GCST90006914  OMP (outer membrane protein)
#   GCST90006915  UreA
#   GCST90006916  VacA

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(R.utils)
})

dir.create("data/exposure", showWarnings = FALSE, recursive = TRUE)
dir.create("data/raw_downloads/gwas_catalog", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

log_path <- file.path("logs", "01_exposure.log")
sink(log_path, split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 01_exposure.R  ", format(Sys.time()), " ===\n")

# ---- Antigen registry -------------------------------------------------------
antigens <- data.frame(
  label = c("hpylori_general","cagA","catalase","groEL","omp","ureA","vacA"),
  gcst  = c("GCST90006910","GCST90006911","GCST90006912","GCST90006913",
            "GCST90006914","GCST90006915","GCST90006916"),
  openg = c("ebi-a-GCST90006910","ebi-a-GCST90006911","ebi-a-GCST90006912",
            "ebi-a-GCST90006913","ebi-a-GCST90006914","ebi-a-GCST90006915",
            "ebi-a-GCST90006916"),
  stringsAsFactors = FALSE
)

# ---- Verify OpenGWAS availability -------------------------------------------
ao <- tryCatch(available_outcomes(), error = function(e) {
  cat("!!! available_outcomes() failed: ", e$message, "\n")
  NULL
})

is_on_opengwas <- function(id) !is.null(ao) && id %in% ao$id

# ---- Fallback: GWAS Catalog FTP --------------------------------------------
# EBI GWAS Catalog FTP layout:
#   https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/<bucket>/<GCST>/
# Bucket = GCST90006001-GCST90007000  for our IDs.
gwas_catalog_url <- function(gcst) {
  bucket <- "GCST90006001-GCST90007000"
  sprintf(
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/%s/%s/",
    bucket, gcst
  )
}

download_from_catalog <- function(gcst, label) {
  base <- gwas_catalog_url(gcst)
  # Catalog file naming is not perfectly standardized. Try common patterns.
  patterns <- c(
    sprintf("%s_buildGRCh37.tsv.gz", gcst),
    sprintf("%s.h.tsv.gz", gcst),
    sprintf("harmonised/%s.h.tsv.gz", gcst),
    sprintf("%s.tsv.gz", gcst)
  )
  dest_dir <- file.path("data/raw_downloads/gwas_catalog", label)
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  for (pat in patterns) {
    url <- paste0(base, pat)
    dest <- file.path(dest_dir, basename(pat))
    cat("  trying:", url, "\n")
    ok <- tryCatch({
      download.file(url, dest, mode = "wb", quiet = TRUE)
      file.exists(dest) && file.size(dest) > 1024
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) {
      cat("  downloaded:", dest, "(", file.size(dest), "bytes)\n")
      return(dest)
    }
  }
  cat("  !!! no Catalog file pattern matched for", gcst, "\n")
  NULL
}

format_catalog_file <- function(path, p_thresh = 5e-8) {
  d <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  cn <- tolower(names(d)); names(d) <- cn
  pick <- function(...) {
    for (x in c(...)) if (x %in% cn) return(x)
    NA_character_
  }
  cols <- list(
    snp    = pick("variant_id","rsid","rsids","snp"),
    beta   = pick("beta","effect"),
    se     = pick("standard_error","se","sebeta"),
    ea     = pick("effect_allele","alt"),
    oa     = pick("other_allele","ref"),
    eaf    = pick("effect_allele_frequency","eaf","af_alt"),
    pval   = pick("p_value","pval","pvalue"),
    n      = pick("n","sample_size","n_total")
  )
  if (any(is.na(unlist(cols[c("snp","beta","se","ea","oa","pval")])))) {
    cat("  !!! required columns missing in", path, " - cols seen:",
        paste(cn, collapse = ","), "\n")
    return(NULL)
  }
  d <- d[get(cols$pval) < p_thresh]
  if (nrow(d) == 0) return(NULL)
  format_data(
    as.data.frame(d), type = "exposure",
    snp_col = cols$snp, beta_col = cols$beta, se_col = cols$se,
    effect_allele_col = cols$ea, other_allele_col = cols$oa,
    eaf_col = if (!is.na(cols$eaf)) cols$eaf else NULL,
    pval_col = cols$pval,
    samplesize_col = if (!is.na(cols$n)) cols$n else NULL
  )
}

local_clump <- function(df, r2 = 0.001, kb = 10000, pop = "EUR") {
  if (is.null(df) || nrow(df) < 2) return(df)
  d <- data.frame(rsid = df$SNP, pval = df$pval.exposure, id = df$id.exposure)
  kept <- tryCatch(
    ieugwasr::ld_clump(d, clump_r2 = r2, clump_kb = kb, pop = pop),
    error = function(e) { cat("  ld_clump failed:", e$message, "\n"); NULL }
  )
  if (is.null(kept)) return(df)
  df[df$SNP %in% kept$rsid, , drop = FALSE]
}

# ---- Main extraction loop ---------------------------------------------------
exposure_list <- list()
source_log <- data.frame()

for (i in seq_len(nrow(antigens))) {
  lab  <- antigens$label[i]
  gcst <- antigens$gcst[i]
  oid  <- antigens$openg[i]
  cat("\n--- Antigen:", lab, " (", gcst, ") ---\n")

  d <- NULL; src <- NA; p_used <- NA

  if (is_on_opengwas(oid)) {
    cat("OpenGWAS hit:", oid, "\n")
    d <- try(extract_instruments(outcomes = oid, p1 = 5e-8,
                                  clump = TRUE, r2 = 0.001, kb = 10000),
             silent = TRUE)
    p_used <- 5e-8
    if (inherits(d, "try-error") || is.null(d) || nrow(d) < 3) {
      cat("  <3 SNPs at 5e-8, retrying at 5e-6\n")
      d <- try(extract_instruments(outcomes = oid, p1 = 5e-6,
                                    clump = TRUE, r2 = 0.001, kb = 10000),
               silent = TRUE)
      p_used <- 5e-6
    }
    if (!inherits(d, "try-error") && !is.null(d) && nrow(d) >= 3) {
      src <- "opengwas"
    } else d <- NULL
  }

  if (is.null(d)) {
    cat("Falling back to GWAS Catalog FTP for", gcst, "\n")
    f <- download_from_catalog(gcst, lab)
    if (!is.null(f)) {
      d <- format_catalog_file(f, p_thresh = 5e-8)
      p_used <- 5e-8
      if (is.null(d) || nrow(d) < 3) {
        d <- format_catalog_file(f, p_thresh = 5e-6)
        p_used <- 5e-6
      }
      if (!is.null(d) && nrow(d) >= 1) {
        d$id.exposure <- gcst
        d$exposure    <- lab
        d <- local_clump(d, r2 = 0.001, kb = 10000, pop = "EUR")
        src <- "gwas_catalog"
      }
    }
  }

  if (is.null(d) || nrow(d) == 0) {
    cat("  !!! NO instruments obtained for", lab, "\n")
    source_log <- rbind(source_log,
      data.frame(label = lab, gcst = gcst, source = NA, p_thresh = NA,
                 n_snp = 0, mean_F = NA, min_F = NA))
    next
  }

  d$exposure_label <- lab
  d$F_stat <- (d$beta.exposure / d$se.exposure)^2
  fwrite(d, sprintf("data/exposure/%s_instruments.tsv", lab), sep = "\t")
  exposure_list[[lab]] <- d
  cat("  ", lab, ":", nrow(d), "SNPs,",
      "mean F =", round(mean(d$F_stat), 1),
      ", min F =", round(min(d$F_stat), 1),
      ", source =", src, ", p =", p_used, "\n")
  source_log <- rbind(source_log,
    data.frame(label = lab, gcst = gcst, source = src, p_thresh = p_used,
               n_snp = nrow(d), mean_F = round(mean(d$F_stat),1),
               min_F = round(min(d$F_stat),1)))
}

saveRDS(exposure_list, "data/exposure/exposure_list.rds")
fwrite(source_log, "data/exposure/source_log.tsv", sep = "\t")

cat("\n=== Summary ===\n")
print(source_log)
cat("\nDone. Proceed to scripts/02_outcome.R\n")
