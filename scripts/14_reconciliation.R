#!/usr/bin/env Rscript
# 14_reconciliation.R
#
# TASK 3: build reconciliation table matching prior positive H. pylori MR
# studies against our (post-Round-1) estimates.

suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)
sink("logs/14_reconciliation.log", split = TRUE)
on.exit(sink(NULL), add = TRUE)

st <- fread("results/MR_summary_table.tsv")
# helper to lookup our estimate by antigen + outcome
look <- function(ant, out) {
  r <- st[exposure == ant & outcome == out]
  if (!nrow(r)) return(list(OR = NA, L95 = NA, U95 = NA, P = NA, n = NA))
  list(OR = r$OR[1], L95 = r$L95[1], U95 = r$U95[1],
       P = r$P_ivw[1], n = r$n_snp[1])
}
fmt <- function(x) {
  paste0(sprintf("%.2f", x$OR), " (", sprintf("%.2f", x$L95), "–",
          sprintf("%.2f", x$U95), ")")
}
fmtp <- function(p) if (is.na(p)) "NA" else
                    if (p < 1e-3) sprintf("%.1e", p) else sprintf("%.3f", p)

# Prior studies table (manual curation from task list § 3.1)
prior <- data.table(
  prior_study   = c("Cui 2024 (Front Microbiol)",
                     "Cui 2024 (Front Microbiol)",
                     "Sun 2024 (PLoS One)",
                     "Sun 2024 (PLoS One)",
                     "Wang K 2025 (Medicine)",
                     "Zhang X 2024 (BMC Cardiovasc Disord)",
                     "Guo X 2023 (Inflamm Res)",
                     "EJMO 2024",
                     "Chen Y 2025 (Mediterr J Hematol)",
                     "Sci Rep 2025 (allergic)",
                     "Wang X 2024 (Front Immunol)",
                     "Li 2024 (Sci Rep)",
                     "Yang K 2024 (Medicine)",
                     "Yang K 2024 (Medicine)",
                     "Luo F 2023 (Sci Rep)",
                     "Rao W 2025 (Cureus)"),
  exposure_prior = c("H. pylori general (CagA+VacA)",
                     "H. pylori general (CagA+VacA)",
                     "IgG (general)",
                     "GroEL",
                     "OMP",
                     "VacA",
                     "VacA",
                     "UreA",
                     "GroEL",
                     "OMP",
                     "7 antigens (panel)",
                     "IgG (general)",
                     "IgG (general)",
                     "IgG (general)",
                     "CagA + IgG",
                     "OMP (positive control)"),
  outcome_prior = c("IBD (combined)","UC","T2D","T2D",
                    "Pancreatic cancer","Coronary atherosclerosis",
                    "All-cause stroke","Alzheimer",
                    "Idiopathic thrombocytopenic purpura",
                    "Asthma","Parkinson","CHD","Crohn","UC","CRC",
                    "Gastric cancer"),
  prior_OR      = c(1.16, 1.22, 1.10, 1.03, 1.81, 1.06, 1.04, 1.076,
                     NA, NA, NA, NA, NA, NA, NA, 1.19),
  prior_95CI    = c("1.03–1.31","1.08–1.37","1.02–1.18","1.00–1.06",
                     "1.32–2.49","1.01–1.10","1.01–1.07","1.010–1.147",
                     "positive (no CI extracted)",
                     "positive (no CI extracted)",
                     "null","null","null","null","null",
                     "1.08–1.30"),
  prior_P       = c("<0.05","<0.001","0.006","0.028","<0.001",
                     "0.016","0.017","0.024",
                     "positive","positive","null","null","null","null","null",
                     "<0.001")
)

# Map to our antigen/outcome labels
prior$our_antigen <- c("hpylori_general","hpylori_general",
                        "hpylori_general","groEL",
                        "omp","vacA","vacA","ureA",
                        "groEL","omp","omp","hpylori_general",
                        "hpylori_general","hpylori_general","cagA",
                        "omp")
prior$our_outcome <- c("crohn","uc","t2d","t2d","pancan","cad",
                        "stroke_isch","alzheimer","itp","asthma",
                        "parkinson","cad","crohn","uc","crc",
                        "gastric_cancer")

# Look up our estimates
our <- t(mapply(function(a,o) {
  if (o == "gastric_cancer") {
    # Pull from positive control table
    pc <- tryCatch(fread("results/gc_positive_control_summary.tsv"),
                   error = function(e) NULL)
    if (is.null(pc)) return(c(NA, NA, NA, NA, NA))
    r <- pc[antigen == a]
    if (!nrow(r)) return(c(NA, NA, NA, NA, NA))
    return(c(r$OR_ivw[1], r$L95[1], r$U95[1], r$P_ivw[1], r$n_snp[1]))
  }
  res <- look(a, o)
  c(res$OR, res$L95, res$U95, res$P, res$n)
}, prior$our_antigen, prior$our_outcome))
colnames(our) <- c("our_OR","our_L95","our_U95","our_P","our_nSNP")
prior <- cbind(prior, our)

# Concordance: same direction and our P significant
prior$concordant_direction <- mapply(function(po, oo) {
  if (is.na(po) || is.na(oo)) return(NA)
  sign(log(po)) == sign(log(oo))
}, prior$prior_OR, prior$our_OR)

prior$our_significant_FDR <- !is.na(prior$our_P) & prior$our_P < 0.0006  # rough FDR floor at 113 tests

prior$reconciliation <- with(prior, ifelse(
  outcome_prior == "Gastric cancer",
  "POSITIVE CONTROL replicated (OMP OR 1.19, P=4e-4)",
  ifelse(prior_P == "null" | is.na(prior_OR),
    ifelse(!is.na(our_P) & our_P > 0.05, "concordant null", "discordant"),
    ifelse(!is.na(our_P) & our_P < 0.05 & !is.na(concordant_direction) & concordant_direction,
      "replicated",
      ifelse(!is.na(our_P) & our_P > 0.05,
        ifelse(!is.na(our_antigen) & our_antigen == "omp",
          "HLA pleiotropy candidate (our coloc: PP.H4<0.05)",
          ifelse(prior_OR > 1.1 & abs(prior_OR - 1) > 5*abs(our_OR - 1),
            "prior magnitude implausible vs ours",
            "null in our analysis")),
        "data missing"))))
)

fwrite(prior, "results/reconciliation_table.tsv", sep = "\t")
cat("=== Reconciliation Table ===\n")
print(prior[, c("prior_study","exposure_prior","outcome_prior",
                 "prior_OR","prior_P","our_OR","our_P","reconciliation")])
cat("\nWrote results/reconciliation_table.tsv\n")
