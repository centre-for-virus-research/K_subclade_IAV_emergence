# ==============================================================================
# Helpers shared by the R test suite
# ==============================================================================

# The analysis scripts in scripts/ are flat, top-to-bottom scripts: sourcing one
# would refit every model, write figures, and (for 02) block on readline().
#
# extract_object() parses a script and evaluates *only* the top-level assignment
# that defines the requested name. Tests therefore exercise the code that ships,
# rather than a copy of it that can silently drift.

project_root <- function() {
  # tests/R/testthat -> tests/R -> tests -> <root>
  normalizePath(file.path(dirname(dirname(dirname(getwd())))), mustWork = TRUE)
}

script_path <- function(name) {
  file.path(project_root(), "scripts", name)
}

data_path <- function(name) {
  file.path(project_root(), "data", name)
}

#' Evaluate one top-level assignment out of a script.
#'
#' @param script   file name inside scripts/
#' @param name     the object to extract
#' @param env_vars named list of values the definition needs in scope
extract_object <- function(script, name, env_vars = list()) {
  exprs <- parse(script_path(script))
  env <- new.env(parent = globalenv())
  for (nm in names(env_vars)) assign(nm, env_vars[[nm]], envir = env)

  found <- FALSE
  for (expr in exprs) {
    if (!is.call(expr)) next
    operator <- as.character(expr[[1]])
    if (!operator %in% c("<-", "=", "<<-")) next
    target <- expr[[2]]
    if (!is.name(target) || as.character(target) != name) next
    eval(expr, envir = env)
    found <- TRUE
    break
  }

  if (!found) stop(sprintf("no top-level definition of '%s' in %s", name, script))
  get(name, envir = env)
}

#' Read a script as plain text (for contract assertions about its source).
script_source <- function(script) {
  readLines(script_path(script), warn = FALSE)
}

#' The scalar numeric constants assigned at the top level of a script.
#'
#' The right-hand side is evaluated in an empty environment rather than tested
#' with is.numeric(): `-50` parses as a call to unary minus, not a numeric
#' literal, so a naive literal check silently misses every negative constant.
extract_numeric_constants <- function(script) {
  exprs <- parse(script_path(script))
  out <- list()
  for (expr in exprs) {
    if (!is.call(expr)) next
    if (!as.character(expr[[1]]) %in% c("<-", "=")) next
    target <- expr[[2]]
    if (!is.name(target)) next

    value <- tryCatch(eval(expr[[3]], envir = new.env(parent = baseenv())),
                      error = function(e) NULL)
    if (is.numeric(value) && length(value) == 1 && !is.na(value)) {
      out[[as.character(target)]] <- value
    }
  }
  out
}

#' Rebuild the analysis frame that 01_neutralisation_analysis.R models.
#'
#' Mirrors sections 1.1 - 3.3 of the script. Kept here so several test files can
#' share one (cached) copy without refitting anything.
build_data_all <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    suppressMessages({
      library(dplyr)
      library(readr)
    })

    neut <- read_csv(data_path("neutralisation_75.csv"), show_col_types = FALSE, progress = FALSE)
    samples <- read_csv(data_path("sample_metadata.csv"),
                        col_types = cols(SampleID = col_character()), progress = FALSE)

    neut <- filter(neut, SampleID != "837SAnti-H3 1:250")

    biorep <- neut %>%
      group_by(SampleID, Virus, Batch) %>%
      summarise(u_Neut = mean(Neut), sd_Neut = sd(Neut), .groups = "drop")

    samples <- filter(samples, SampleID %in% neut$SampleID)
    samples <- samples %>%
      arrange(PatientID, SampleDate) %>%
      distinct(PatientID, .keep_all = TRUE)

    biorep$Batch <- factor(biorep$Batch)
    biorep$u_Neut <- pmax(pmin(1, biorep$u_Neut / 100), 0)
    n <- nrow(biorep)
    biorep$u_Neut <- (biorep$u_Neut * (n - 1) + 0.5) / n

    all_data <- left_join(biorep, samples, by = "SampleID")
    all_data <- filter(all_data, !is.na(PatientID))
    all_data$Virus <- factor(all_data$Virus, levels = c("THA22", "ENG24", "ENG25"))
    all_data$vacc_years_since_last <- factor(all_data$vacc_years_since_last,
                                             levels = c(">4", "4", "3", "2", "1"))
    all_data$vacc_group <- ifelse(all_data$vacc_years_since_last == ">4", ">4",
                           ifelse(all_data$vacc_years_since_last == "1", "1", "2-4"))
    all_data$vacc_group <- factor(all_data$vacc_group, levels = c(">4", "2-4", "1"))
    all_data$Age_policy <- cut(all_data$Age, breaks = c(18, 50, 65, Inf), right = FALSE,
                               labels = c("18-49", "50-64", "65+"))

    cached <<- all_data
    cached
  }
})
