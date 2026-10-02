# ==============================================================================
# Contracts between the shipped CSVs and the R analysis scripts.
#
# The scripts lean heavily on factor(levels = ...) and cut(breaks = ...). Both
# turn unexpected values into NA *silently*, and betareg then drops those rows
# without comment, so a single stray level could shrink the analysis set without
# anything in the console saying so. These tests assert that every observed
# value is one the scripts actually name.
# ==============================================================================

suppressMessages({
  library(readr)
  library(dplyr)
})

neut <- read_csv(data_path("neutralisation_75.csv"), show_col_types = FALSE, progress = FALSE)
meta <- read_csv(data_path("sample_metadata.csv"),
                 col_types = cols(SampleID = col_character()), progress = FALSE)
ic50s <- read_csv(data_path("IC50s.csv"),
                  col_types = cols(Sample = col_character()), progress = FALSE)
fits <- read_csv(data_path("IC50_fits.csv"),
                 col_types = cols(SampleID = col_character()), progress = FALSE)

# ------------------------------------------------------------------------------
# neutralisation_75.csv
# ------------------------------------------------------------------------------

test_that("neutralisation data has the columns the scripts read", {
  expect_true(all(c("SampleID", "Batch", "Virus", "Replicate", "Plate", "Well", "Neut") %in% names(neut)))
})

test_that("neutralisation data is not empty and has no missing Neut values", {
  expect_gt(nrow(neut), 0)
  expect_equal(sum(is.na(neut$Neut)), 0)
})

test_that("Virus takes only the three levels the scripts order", {
  expect_setequal(unique(neut$Virus), c("THA22", "ENG24", "ENG25"))
})

test_that("factoring Virus with the scripts' levels loses nothing", {
  expect_equal(sum(is.na(factor(neut$Virus, levels = c("THA22", "ENG24", "ENG25")))), 0)
})

test_that("the sheep antiserum control is present so that filtering it is meaningful", {
  expect_gt(sum(neut$SampleID == "837SAnti-H3 1:250"), 0)
})

test_that("removing the sheep antiserum leaves only patient samples", {
  kept <- filter(neut, SampleID != "837SAnti-H3 1:250")
  expect_false(any(grepl("Anti-H3", kept$SampleID)))
})

test_that("neutralisation percentages sit in a plausible assay range", {
  # Values below 0 and at exactly 100 are expected; the script clamps them.
  expect_lte(max(neut$Neut), 100)
  expect_gt(min(neut$Neut), -200)
})

test_that("each virus is measured on every batch", {
  counts <- table(neut$Batch, neut$Virus)
  expect_true(all(counts > 0))
})

# ------------------------------------------------------------------------------
# sample_metadata.csv
# ------------------------------------------------------------------------------

test_that("sample metadata has the columns the scripts read", {
  expect_true(all(c("PatientID", "SampleID", "SampleDate", "Age", "Sex",
                    "vacc_years_since_last", "vacc_number_of_doses") %in% names(meta)))
})

test_that("SampleID is unique in the metadata", {
  expect_equal(sum(duplicated(meta$SampleID)), 0)
})

test_that("repeat patients exist, so the distinct() de-duplication is load-bearing", {
  expect_gt(sum(duplicated(meta$PatientID)), 0)
})

test_that("keeping the earliest sample per patient yields one row per patient", {
  deduplicated <- meta %>%
    arrange(PatientID, SampleDate) %>%
    distinct(PatientID, .keep_all = TRUE)

  expect_equal(sum(duplicated(deduplicated$PatientID)), 0)
  expect_equal(nrow(deduplicated), length(unique(meta$PatientID)))
})

test_that("de-duplication really does keep the earliest date", {
  deduplicated <- meta %>%
    arrange(PatientID, SampleDate) %>%
    distinct(PatientID, .keep_all = TRUE)

  earliest <- meta %>% group_by(PatientID) %>% summarise(first = min(SampleDate), .groups = "drop")
  joined <- left_join(deduplicated, earliest, by = "PatientID")

  expect_true(all(joined$SampleDate == joined$first))
})

test_that("Sex takes only the two levels the figures order", {
  expect_setequal(unique(meta$Sex), c("Female", "Male"))
  expect_equal(sum(is.na(factor(meta$Sex, levels = c("Female", "Male")))), 0)
})

test_that("vacc_years_since_last takes only the levels script 01 names", {
  levels_used <- c(">4", "4", "3", "2", "1")
  expect_true(all(unique(meta$vacc_years_since_last) %in% levels_used))
  expect_equal(sum(is.na(factor(meta$vacc_years_since_last, levels = levels_used))), 0)
})

test_that("the coarse vaccination grouping covers every record", {
  grouped <- ifelse(meta$vacc_years_since_last == ">4", ">4",
             ifelse(meta$vacc_years_since_last == "1", "1", "2-4"))

  expect_equal(sum(is.na(grouped)), 0)
  expect_setequal(unique(grouped), c(">4", "2-4", "1"))
  expect_equal(sum(is.na(factor(grouped, levels = c(">4", "2-4", "1")))), 0)
})

test_that("every age falls inside the cut() ranges the scripts define", {
  expect_equal(sum(is.na(meta$Age)), 0)
  expect_gte(min(meta$Age), 18)

  bins <- list(
    Age_10y = list(breaks = c(18, 30, 40, 50, 60, 70, 80, 90, 100),
                   labels = c("18-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+")),
    Age_20y = list(breaks = c(18, 40, 60, 80, Inf), labels = c("18-39", "40-59", "60-79", "80+")),
    Age_policy = list(breaks = c(18, 50, 65, Inf), labels = c("18-49", "50-64", "65+")),
    Age_3groups = list(breaks = c(18, 40, 65, Inf), labels = c("18-39", "40-64", "65+"))
  )

  for (name in names(bins)) {
    binned <- cut(meta$Age, breaks = bins[[name]]$breaks, right = FALSE, labels = bins[[name]]$labels)
    expect_equal(sum(is.na(binned)), 0, info = sprintf("%s drops %d ages", name, sum(is.na(binned))))
  }
})

test_that("Age_10y is the bin at risk if anyone reaches 100", {
  # Documents why the open-ended bins are the safer choice.
  expect_true(is.na(cut(100, breaks = c(18, 30, 40, 50, 60, 70, 80, 90, 100), right = FALSE)))
  expect_false(is.na(cut(100, breaks = c(18, 50, 65, Inf), right = FALSE)))
})

# ------------------------------------------------------------------------------
# Joins
# ------------------------------------------------------------------------------

test_that("most neutralisation samples carry metadata", {
  orphans <- setdiff(neut$SampleID, meta$SampleID)
  expect_lt(length(orphans), 10)
})

test_that("samples without metadata are dropped by the PatientID filter", {
  orphans <- setdiff(neut$SampleID, meta$SampleID)
  joined <- left_join(tibble(SampleID = orphans), meta, by = "SampleID")
  expect_true(all(is.na(joined$PatientID)))
})

test_that("the assembled analysis frame is non-empty and fully specified", {
  all_data <- build_data_all()

  expect_gt(nrow(all_data), 0)
  for (column in c("Virus", "Batch", "Sex", "vacc_group", "Age_policy", "u_Neut")) {
    expect_equal(sum(is.na(all_data[[column]])), 0,
                 info = sprintf("%s contains NA after assembly", column))
  }
})

test_that("every virus x vaccination x age cell the models fit is populated", {
  all_data <- build_data_all()
  counts <- table(all_data$Virus, all_data$vacc_group, all_data$Age_policy)
  expect_true(all(counts > 0), info = "an interaction cell is empty; emmeans would return NA")
})

# ------------------------------------------------------------------------------
# The Smithson & Verkuilen transform
# ------------------------------------------------------------------------------

test_that("the response betareg receives is strictly inside (0, 1)", {
  all_data <- build_data_all()

  expect_true(all(all_data$u_Neut > 0))
  expect_true(all(all_data$u_Neut < 1))
})

test_that("the transform maps the closed unit interval into the open one", {
  n <- 100
  squeezed <- (c(0, 0.5, 1) * (n - 1) + 0.5) / n

  expect_true(all(squeezed > 0))
  expect_true(all(squeezed < 1))
  expect_equal(squeezed[2], 0.5)
})

test_that("clamping happens before the transform, so out-of-range assay values cannot escape", {
  raw <- c(-54.19, 0, 50, 100, 120)
  clamped <- pmax(pmin(1, raw / 100), 0)

  expect_equal(clamped, c(0, 0, 0.5, 1, 1))
})

# ------------------------------------------------------------------------------
# IC50s.csv and IC50_fits.csv
# ------------------------------------------------------------------------------

test_that("IC50 dilution data has the columns script 02 reads", {
  expect_true(all(c("Batch", "Virus", "Sample", "Dilution", "Neut") %in% names(ic50s)))
})

test_that("dilutions form the expected two-fold series", {
  dilutions <- sort(unique(ic50s$Dilution))
  expect_equal(dilutions, c(50, 100, 200, 400, 800, 1600, 3200))
})

test_that("every sample x virus combination has a full dilution series", {
  per_combination <- ic50s %>%
    group_by(Sample, Virus) %>%
    summarise(n_dilutions = n_distinct(Dilution), .groups = "drop")

  expect_true(all(per_combination$n_dilutions == 7))
})

test_that("IC50 fits have the columns script 03 reads", {
  expect_true(all(c("SampleID", "Virus", "Verdict", "Lower", "Upper", "IC50", "Slope") %in% names(fits)))
})

test_that("the BOM on IC50_fits.csv does not corrupt the first column name", {
  expect_true("SampleID" %in% names(fits))
  expect_false(any(grepl("^﻿", names(fits))))
})

test_that("Verdict takes only the two values script 03 branches on", {
  expect_setequal(unique(fits$Verdict), c("Curved", "Flattened"))
})

test_that("flattened fits have an equal lower and upper asymptote", {
  flattened <- filter(fits, Verdict == "Flattened")
  expect_true(all(flattened$Lower == flattened$Upper))
})

test_that("flattened fits carry no slope and no IC50", {
  flattened <- filter(fits, Verdict == "Flattened")
  expect_true(all(flattened$Slope == 0))
  expect_true(all(is.na(flattened$IC50)))
})

test_that("every fitted sample has metadata to join against", {
  expect_equal(length(setdiff(fits$SampleID, meta$SampleID)), 0)
})

test_that("every fit corresponds to a real dilution series", {
  combinations <- distinct(ic50s, Sample, Virus)
  unmatched <- anti_join(distinct(fits, SampleID, Virus), combinations,
                         by = c("SampleID" = "Sample", "Virus" = "Virus"))
  expect_equal(nrow(unmatched), 0)
})
