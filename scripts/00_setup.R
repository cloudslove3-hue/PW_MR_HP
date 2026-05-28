#!/usr/bin/env Rscript
# 00_setup.R  H. pylori MR pipeline - environment setup
#
# Run once. Idempotent: skips already-installed packages.
# Usage from project root:   Rscript scripts/00_setup.R
#
# Prerequisites:
#   - R >= 4.4
#   - Internet access
#   - OPENGWAS_JWT token saved in ~/.Renviron  (free, https://api.opengwas.io)

stopifnot(getRversion() >= "4.0")

log_path <- file.path("logs", "00_setup.log")
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
sink(log_path, split = TRUE)
on.exit(sink(NULL), add = TRUE)

cat("=== 00_setup.R  ", format(Sys.time()), " ===\n")
cat("R version: ", R.version.string, "\n")
cat("Platform: ",  R.version$platform,  "\n\n")

# ---- User library (Windows: Program Files\R\library is not writable) -------
user_lib <- Sys.getenv("R_LIBS_USER")
if (nchar(user_lib) == 0 || user_lib == "NULL") {
  user_lib <- file.path(Sys.getenv("USERPROFILE"),
                         sprintf("R/win-library/%s.%s",
                                 R.version$major, substr(R.version$minor,1,1)))
}
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
cat("Using user library:", user_lib, "\n\n")

# ---- CRAN packages ----------------------------------------------------------
cran_pkgs <- c(
  "renv", "remotes", "data.table", "dplyr", "ggplot2", "readr",
  "stringr", "MendelianRandomization", "coloc", "R.utils"
)
to_install <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install)) {
  cat("Installing CRAN packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, lib = user_lib, repos = "https://cloud.r-project.org")
} else {
  cat("All CRAN packages present.\n")
}

# ---- TwoSampleMR + ieugwasr via MRCIEU r-universe (no GitHub auth needed) ---
mrcieu_pkgs <- c("TwoSampleMR", "ieugwasr")
mrcieu_missing <- mrcieu_pkgs[!sapply(mrcieu_pkgs, requireNamespace, quietly = TRUE)]
if (length(mrcieu_missing)) {
  cat("Installing from r-universe:", paste(mrcieu_missing, collapse = ", "), "\n")
  install.packages(mrcieu_missing, lib = user_lib,
                   repos = c("https://mrcieu.r-universe.dev",
                             "https://cloud.r-project.org"))
} else {
  cat("MRCIEU packages already installed.\n")
}

# ---- MR-PRESSO: try r-universe, fall back to anonymous GitHub --------------
if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  cat("Trying MRPRESSO via r-universe...\n")
  try(install.packages("MRPRESSO", lib = user_lib,
                        repos = c("https://rondolab.r-universe.dev",
                                  "https://cloud.r-project.org")),
      silent = TRUE)
}
if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  cat("Falling back to anonymous GitHub for MRPRESSO\n")
  # Strip any bad PAT from env to force anonymous (60 req/hr)
  Sys.setenv(GITHUB_PAT = "", GITHUB_TOKEN = "")
  remotes::install_github("rondolab/MR-PRESSO", upgrade = "never",
                           lib = user_lib, auth_token = NULL)
}
if (requireNamespace("MRPRESSO", quietly = TRUE)) {
  cat("MRPRESSO installed.\n")
} else {
  cat("!!! MRPRESSO install failed. Pipeline will skip PRESSO step.\n")
}

# ---- Token sanity check -----------------------------------------------------
tok <- Sys.getenv("OPENGWAS_JWT")
if (nchar(tok) < 20) {
  cat("\n!!! OPENGWAS_JWT not found in environment.\n")
  cat("    Get a free token at https://api.opengwas.io  (Google sign-in)\n")
  cat("    Then add to ~/.Renviron:\n")
  cat("        OPENGWAS_JWT=eyJ.....\n")
  cat("    Restart R and re-run this script.\n")
} else {
  cat("\nOPENGWAS_JWT detected (", nchar(tok), " chars).\n", sep = "")
  ok <- try(ieugwasr::user(), silent = TRUE)
  if (inherits(ok, "try-error")) {
    cat("!!! Token present but ieugwasr::user() failed. Token may be expired.\n")
  } else {
    cat("Token valid. Logged in as: ", ok$jwt_user_email %||% "(unknown)\n", sep = "")
  }
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a) || identical(a, "")) b else a

cat("\nDone. Proceed to scripts/01_exposure.R\n")
