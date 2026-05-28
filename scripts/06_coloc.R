#!/usr/bin/env Rscript
# 06_coloc.R  Colocalization for OMP -> Crohn signal
#
# Per the borderline IVW signal (P=0.024, OR=0.836), we test whether OMP and
# Crohn share causal variants in each of the 10 OMP instrument loci.
#
# Method: Bayesian coloc.abf (Giambartolomei 2014). PP.H4 > 0.75 = strong
# colocalization, 0.5-0.75 = suggestive, <0.5 = no shared causal variant.
#
# CRITICAL: OMP lead SNP is in the HLA region (chr6:32.65Mb). A high PP.H4 in
# HLA does NOT necessarily mean H. pylori biology drives Crohn risk -- HLA
# pleiotropy is the more parsimonious explanation. Non-HLA loci are more
# informative.
#
# Inputs:
#   data/exposure/omp_instruments.tsv  (10 SNPs with chr/pos)
# Outputs:
#   results/coloc_omp_crohn.tsv
#   logs/06_coloc.log

suppressPackageStartupMessages({
  library(data.table)
  library(ieugwasr)
})

# coloc is already installed via 00_setup.R; gwasglue not strictly needed
# because we fetch via ieugwasr::associations() directly
stopifnot(requireNamespace("coloc", quietly = TRUE))

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("logs",    showWarnings = FALSE, recursive = TRUE)
sink("logs/06_coloc.log", split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 06_coloc.R  ", format(Sys.time()), " ===\n")

omp <- fread("data/exposure/omp_instruments.tsv")
cat("OMP instruments:", nrow(omp), "\n")
print(omp[, .(SNP, chr.exposure, pos.exposure, pval.exposure)])

# ---- Manual regional coloc using ieugwasr::associations() ------------------
# We pull both GWAS data at all variants within ±500kb of each OMP lead SNP.
# `associations()` accepts a chr:pos range via the variants arg as "chr:start-end".

run_region_coloc <- function(snp_id, chr, pos, window = 5e5,
                              exp_id = "ebi-a-GCST90006914",
                              out_id = "ieu-a-30") {
  range_str <- sprintf("%s:%d-%d", chr, max(1, pos - window), pos + window)
  cat("\n--- Region:", snp_id, " (", range_str, ") ---\n")

  d1 <- tryCatch(associations(variants = range_str, id = exp_id),
                 error = function(e) { cat("  exp fetch fail:", e$message, "\n"); NULL })
  d2 <- tryCatch(associations(variants = range_str, id = out_id),
                 error = function(e) { cat("  out fetch fail:", e$message, "\n"); NULL })
  if (is.null(d1) || is.null(d2) || nrow(d1) == 0 || nrow(d2) == 0) {
    cat("  empty fetch, skipping\n")
    return(NULL)
  }
  cat("  exp n=", nrow(d1), " out n=", nrow(d2), "\n")

  # Intersect on rsid
  common <- intersect(d1$rsid, d2$rsid)
  if (length(common) < 50) {
    cat("  too few common variants (", length(common), "), skipping\n")
    return(NULL)
  }
  d1 <- d1[match(common, d1$rsid), ]
  d2 <- d2[match(common, d2$rsid), ]

  # Build coloc datasets
  # exposure (OMP) is quantitative continuous (antibody level)
  ds1 <- list(
    snp = common,
    beta = d1$beta, varbeta = d1$se^2,
    type = "quant", N = d1$n[1], MAF = d1$eaf,
    sdY = 1  # standardized
  )
  # outcome (Crohn) is case-control
  ncase <- 5956; ncont <- 14927  # IIBDGC Liu 2015 EUR
  ds2 <- list(
    snp = common,
    beta = d2$beta, varbeta = d2$se^2,
    type = "cc", N = ncase + ncont, s = ncase / (ncase + ncont), MAF = d2$eaf
  )

  res <- tryCatch(coloc::coloc.abf(ds1, ds2),
                  error = function(e) { cat("  coloc.abf err:", e$message, "\n"); NULL })
  if (is.null(res)) return(NULL)
  s <- res$summary
  cat(sprintf("  n_snp=%d  PP.H0=%.3f  H1=%.3f  H2=%.3f  H3=%.3f  H4=%.3f\n",
              s["nsnps"], s["PP.H0.abf"], s["PP.H1.abf"],
              s["PP.H2.abf"], s["PP.H3.abf"], s["PP.H4.abf"]))
  data.frame(
    region_lead = snp_id, chr = chr, pos = pos,
    n_common = s["nsnps"],
    PP_H0 = s["PP.H0.abf"], PP_H1 = s["PP.H1.abf"],
    PP_H2 = s["PP.H2.abf"], PP_H3 = s["PP.H3.abf"], PP_H4 = s["PP.H4.abf"],
    interpretation = ifelse(s["PP.H4.abf"] > 0.75, "strong_coloc",
                     ifelse(s["PP.H4.abf"] > 0.5, "suggestive_coloc",
                     ifelse(s["PP.H3.abf"] > 0.5, "distinct_causal_variants",
                            "underpowered_or_no_shared_causal")))
  )
}

# Loop over all 10 OMP instruments
results <- list()
for (i in seq_len(nrow(omp))) {
  results[[i]] <- run_region_coloc(
    snp_id = omp$SNP[i],
    chr    = omp$chr.exposure[i],
    pos    = omp$pos.exposure[i]
  )
  Sys.sleep(1)
}

out <- do.call(rbind, results)
if (!is.null(out)) {
  fwrite(out, "results/coloc_omp_crohn.tsv", sep = "\t")
  cat("\n=== Summary ===\n"); print(out)
} else {
  cat("\n!!! No regions successfully colocalized\n")
}
cat("\nDone.\n")
