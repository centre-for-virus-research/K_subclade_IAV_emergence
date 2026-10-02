# ==============================================================================
# Entry point for the R test suite.
#
#   Rscript tests/R/testthat.R            # from the repository root
#   Rscript tests/R/testthat.R <dir>      # or point it at the testthat/ folder
#
# Exits non-zero if anything fails, so it can be wired straight into CI.
# ==============================================================================

suppressMessages(library(testthat))

args <- commandArgs(trailingOnly = TRUE)

find_testthat_dir <- function() {
  if (length(args) >= 1 && dir.exists(args[[1]])) return(normalizePath(args[[1]]))

  # Resolve relative to this file when invoked as `Rscript <path>/testthat.R`.
  full <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", full[grepl("^--file=", full)])
  if (length(file_arg) == 1) {
    candidate <- file.path(dirname(normalizePath(file_arg)), "testthat")
    if (dir.exists(candidate)) return(candidate)
  }

  candidate <- file.path(getwd(), "tests", "R", "testthat")
  if (dir.exists(candidate)) return(normalizePath(candidate))

  stop("could not locate tests/R/testthat/; pass it as an argument")
}

testthat_dir <- find_testthat_dir()
cat("Running R tests in:", testthat_dir, "\n\n")

results <- test_dir(testthat_dir, reporter = "summary", stop_on_failure = FALSE)

frame <- as.data.frame(results)
failed <- sum(frame$failed) + sum(frame$error)

cat("\n")
cat(sprintf("passed:  %d\n", sum(frame$passed)))
cat(sprintf("failed:  %d\n", sum(frame$failed)))
cat(sprintf("errors:  %d\n", sum(frame$error)))
cat(sprintf("skipped: %d\n", sum(frame$skipped)))

if (failed > 0) quit(status = 1)
