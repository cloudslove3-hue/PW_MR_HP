#!/usr/bin/env Rscript
# 17_coloc_extension.R
#
# Round 2 expansion: coloc.abf at every instrument locus for the three
# remaining nominally-significant Tier 1 pairs AND the gastric cancer
# positive control. Mirrors 06_coloc.R logic.
#
# Pairs:
#   OMP   -> gastric cancer (positive control; ebi-a-GCST90018849)
#   VacA  -> ischaemic stroke (ebi-a-GCST006908)
#   Catalase -> colorectal cancer (ebi-a-GCST90018808)

set.seed(20260521)

suppressPackageStartupMessages({
  library(data.table)
  library(ieugwasr)
  library(coloc)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("logs",    showWarnings = FALSE, recursive = TRUE)
sink("logs/17_coloc_extension.log", split = TRUE)
on.exit(sink(NULL), add = TRUE)
cat("=== 17_coloc_extension.R  ", format(Sys.time()), " ===\n")

# Antigen instrument files
ANT <- list(
  omp      = list(file = "data/exposure/omp_instruments.tsv",
                  id   = "ebi-a-GCST90006914",
                  type = "quant", sdY = 1),
  vacA     = list(file = "data/exposure/vacA_instruments.tsv",
                  id   = "ebi-a-GCST90006916",
                  type = "quant", sdY = 1),
  catalase = list(file = "data/exposure/catalase_instruments.tsv",
                  id   = "ebi-a-GCST90006912",
                  type = "quant", sdY = 1)
)

# Outcome registry
OUT <- list(
  gastric_cancer = list(id = "ebi-a-GCST90018849", type = "cc",
                         ncase = 1029, ncontrol = 475087),
  stroke_isch    = list(id = "ebi-a-GCST006908",   type = "cc",
                         ncase = 34217, ncontrol = 406111),
  crc            = list(id = "ebi-a-GCST90018808", type = "cc",
                         ncase = 6581, ncontrol = 463421)
)

# Pairs to run
PAIRS <- list(
  c("omp",      "gastric_cancer"),
  c("vacA",     "stroke_isch"),
  c("catalase", "crc")
)

run_region <- function(snp_id, chr, pos, exp_lab, out_lab, window = 5e5) {
  exp_id <- ANT[[exp_lab]]$id
  out_id <- OUT[[out_lab]]$id
  range_str <- sprintf("%s:%d-%d", chr, max(1, pos - window), pos + window)
  cat(sprintf("\n--- %s -> %s  [%s]  region %s ---\n",
              exp_lab, out_lab, snp_id, range_str))

  d1 <- tryCatch(associations(variants = range_str, id = exp_id),
                 error = function(e) NULL)
  d2 <- tryCatch(associations(variants = range_str, id = out_id),
                 error = function(e) NULL)
  if (is.null(d1) || is.null(d2) || nrow(d1) == 0 || nrow(d2) == 0) {
    cat("  empty fetch\n"); return(NULL)
  }
  common <- intersect(d1$rsid, d2$rsid)
  if (length(common) < 50) { cat("  <50 common SNPs\n"); return(NULL) }
  d1 <- d1[match(common, d1$rsid), ]; d2 <- d2[match(common, d2$rsid), ]

  ds1 <- list(snp = common, beta = d1$beta, varbeta = d1$se^2,
               type = ANT[[exp_lab]]$type, N = d1$n[1],
               MAF = d1$eaf, sdY = ANT[[exp_lab]]$sdY)
  ncase   <- OUT[[out_lab]]$ncase
  ncontrol<- OUT[[out_lab]]$ncontrol
  ds2 <- list(snp = common, beta = d2$beta, varbeta = d2$se^2,
               type = "cc", N = ncase + ncontrol, s = ncase / (ncase + ncontrol),
               MAF = d2$eaf)

  r <- tryCatch(coloc.abf(ds1, ds2), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  s <- r$summary
  cat(sprintf("  n=%d  PP.H0=%.3f H1=%.3f H2=%.3f H3=%.3f H4=%.3f\n",
              s["nsnps"], s["PP.H0.abf"], s["PP.H1.abf"],
              s["PP.H2.abf"], s["PP.H3.abf"], s["PP.H4.abf"]))
  data.frame(pair = paste(exp_lab, out_lab, sep = "__"),
             region_lead = snp_id, chr = chr, pos = pos,
             n_common = s["nsnps"],
             PP_H0 = s["PP.H0.abf"], PP_H1 = s["PP.H1.abf"],
             PP_H2 = s["PP.H2.abf"], PP_H3 = s["PP.H3.abf"],
             PP_H4 = s["PP.H4.abf"],
             interpretation = ifelse(s["PP.H4.abf"] > 0.75, "strong_coloc",
                              ifelse(s["PP.H4.abf"] > 0.5,  "suggestive_coloc",
                              ifelse(s["PP.H3.abf"] > 0.5,  "distinct_causal",
                                     "underpowered_or_no_shared"))))
}

all_rows <- list()
for (p in PAIRS) {
  exp_lab <- p[1]; out_lab <- p[2]
  iv <- fread(ANT[[exp_lab]]$file)
  cat("\n=== ", exp_lab, " -> ", out_lab, "  (", nrow(iv), "loci) ===\n", sep = "")
  for (i in seq_len(nrow(iv))) {
    r <- run_region(iv$SNP[i], iv$chr.exposure[i], iv$pos.exposure[i],
                    exp_lab, out_lab)
    if (!is.null(r)) all_rows[[length(all_rows) + 1]] <- r
    Sys.sleep(1)
  }
}

out <- do.call(rbind, all_rows)
if (!is.null(out) && nrow(out)) {
  fwrite(out, "results/coloc_extended.tsv", sep = "\t")
  cat("\nWrote results/coloc_extended.tsv  (", nrow(out), "rows )\n")
  cat("\n=== Summary by pair (max PP.H4) ===\n")
  agg <- aggregate(PP_H4 ~ pair, data = out, FUN = max)
  print(agg)
} else cat("\n!!! No coloc rows produced\n")

cat("\nDone.\n")
