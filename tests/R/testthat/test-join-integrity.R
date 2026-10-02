# ==============================================================================
# Join integrity across the three R scripts
#
# Both 01 and 03 assemble their analysis frame with a chain of left_join()s on
# keys that are never asserted to match. A left join cannot fail loudly: an
# unmatched key produces NA columns, the model silently drops those rows, and
# the only symptom is a smaller n that nothing prints. These tests pin the
# number of rows each join is expected to carry, so a change in any input file
# that breaks a key shows up here rather than as a quietly different result.
# ==============================================================================

suppressMessages({
  library(readr)
  library(dplyr)
})

neut <- read_csv(data_path("neutralisation_75.csv"), show_col_types = FALSE, progress = FALSE)
meta <- read_csv(data_path("sample_metadata.csv"),
                 col_types = cols(SampleID = col_character()), progress = FALSE)
fits <- read_csv(data_path("IC50_fits.csv"),
                 col_types = cols(SampleID = col_character()), progress = FALSE)

# ------------------------------------------------------------------------------
# Key hygiene
# ------------------------------------------------------------------------------

test_that("join keys are typed as character everywhere they are read", {
  # SampleIDs look numeric in places; readr guessing numeric in one file and
  # character in another would make every join silently fail.
  expect_type(meta$SampleID, "character")
  expect_type(fits$SampleID, "character")
})

test_that("no join key carries stray whitespace", {
  for (key in list(neut$SampleID, meta$SampleID, fits$SampleID)) {
    expect_equal(sum(key != trimws(key)), 0)
  }
})

test_that("no join key is missing", {
  expect_equal(sum(is.na(neut$SampleID)), 0)
  expect_equal(sum(is.na(meta$SampleID)), 0)
  expect_equal(sum(is.na(fits$SampleID)), 0)
})

test_that("the metadata key is unique, so no join can fan out rows", {
  expect_equal(sum(duplicated(meta$SampleID)), 0)
})

test_that("the IC50 fits key is unique per sample and virus", {
  expect_equal(sum(duplicated(fits[, c("SampleID", "Virus")])), 0)
})

# ------------------------------------------------------------------------------
# Script 01: neutralisation -> metadata
# ------------------------------------------------------------------------------

test_that("the script 01 join neither drops nor duplicates biological replicates", {
  biorep <- neut %>%
    filter(SampleID != "837SAnti-H3 1:250") %>%
    group_by(SampleID, Virus, Batch) %>%
    summarise(u_Neut = mean(Neut), .groups = "drop")

  samples <- meta %>%
    filter(SampleID %in% neut$SampleID) %>%
    arrange(PatientID, SampleDate) %>%
    distinct(PatientID, .keep_all = TRUE)

  joined <- left_join(biorep, samples, by = "SampleID")

  expect_equal(nrow(joined), nrow(biorep))
})

test_that("rows dropped by the PatientID filter are exactly the unmatched ones", {
  all_data <- build_data_all()

  biorep <- neut %>%
    filter(SampleID != "837SAnti-H3 1:250") %>%
    group_by(SampleID, Virus, Batch) %>%
    summarise(u_Neut = mean(Neut), .groups = "drop")

  samples <- meta %>%
    filter(SampleID %in% neut$SampleID) %>%
    arrange(PatientID, SampleDate) %>%
    distinct(PatientID, .keep_all = TRUE)

  expected_kept <- sum(biorep$SampleID %in% samples$SampleID)
  expect_equal(nrow(all_data), expected_kept)
})

test_that("the analysis frame keeps a stated number of observations", {
  # A regression anchor: if an input file changes, this is the first thing to
  # move, and every AIC comparison in section 4 depends on it.
  expect_equal(nrow(build_data_all()), 1458)
})

test_that("the analysis frame holds three viruses per retained sample", {
  all_data <- build_data_all()
  per_sample <- count(all_data, SampleID)

  expect_true(all(per_sample$n == 3),
              info = "a sample is missing one of the three virus measurements")
})

test_that("every retained patient appears exactly once per virus", {
  all_data <- build_data_all()
  expect_equal(sum(duplicated(all_data[, c("PatientID", "Virus")])), 0)
})

# ------------------------------------------------------------------------------
# Script 03: fits -> neutralisation -> metadata
# ------------------------------------------------------------------------------

test_that("the script 03 join chain preserves one row per fit", {
  neut_mean <- neut %>%
    group_by(SampleID, Virus) %>%
    summarise(u_Neut = mean(Neut, na.rm = TRUE), .groups = "drop")

  comp <- fits %>%
    left_join(neut_mean, by = c("SampleID", "Virus")) %>%
    left_join(meta, by = "SampleID")

  expect_equal(nrow(comp), nrow(fits))
})

test_that("every fitted sample picks up its 1:75 neutralisation value", {
  neut_mean <- neut %>%
    group_by(SampleID, Virus) %>%
    summarise(u_Neut = mean(Neut, na.rm = TRUE), .groups = "drop")

  comp <- left_join(fits, neut_mean, by = c("SampleID", "Virus"))

  expect_equal(sum(is.na(comp$u_Neut)), 0,
               info = "an IC50 fit has no matching 1:75 measurement")
})

test_that("every fitted sample picks up its vaccination metadata", {
  comp <- left_join(fits, meta, by = "SampleID")

  expect_equal(sum(is.na(comp$vacc_years_since_last)), 0)
  expect_equal(sum(is.na(comp$Age)), 0)
  expect_equal(sum(is.na(comp$Sex)), 0)
})

test_that("the correlation subset is a known size", {
  # Section 1.4 filters to rows with both an IC50 and a neutralisation value.
  neut_mean <- neut %>%
    group_by(SampleID, Virus) %>%
    summarise(u_Neut = mean(Neut, na.rm = TRUE), .groups = "drop")

  comp <- fits %>%
    left_join(neut_mean, by = c("SampleID", "Virus")) %>%
    filter(!is.na(IC50), !is.na(u_Neut))

  expect_equal(nrow(comp), sum(!is.na(fits$IC50)))
  expect_gt(nrow(comp), 250)
})

test_that("the IC50 model frame loses only the rows with no IC50 at all", {
  comp <- left_join(fits, meta, by = "SampleID")
  comp$IC50 <- ifelse(comp$Verdict == "Flattened",
                      ifelse(comp$Lower < 50, -4, -12), comp$IC50)

  usable <- sum(!is.na(comp$IC50))
  expect_equal(usable, nrow(fits) - sum(fits$Verdict == "Curved" & is.na(fits$IC50)))
  expect_equal(usable, 350)
})

# ------------------------------------------------------------------------------
# Script 01 technical replicate structure
# ------------------------------------------------------------------------------

test_that("technical replicate counts are consistent enough to take an sd", {
  replicates <- neut %>%
    filter(SampleID != "837SAnti-H3 1:250") %>%
    count(SampleID, Virus, Batch)

  # sd() of a single observation is NA, which would propagate into CV_Neut.
  expect_true(all(replicates$n >= 2),
              info = sprintf("%d groups have a single technical replicate",
                             sum(replicates$n < 2)))
})

test_that("the technical-error columns are computable for every group", {
  biorep <- neut %>%
    filter(SampleID != "837SAnti-H3 1:250") %>%
    group_by(SampleID, Virus, Batch) %>%
    summarise(u_Neut = mean(Neut), sd_Neut = sd(Neut), .groups = "drop")

  expect_equal(sum(is.na(biorep$sd_Neut)), 0)

  cv <- biorep$sd_Neut / (biorep$u_Neut + abs(min(biorep$u_Neut)) + 1)
  expect_equal(sum(is.na(cv)), 0)
  expect_true(all(is.finite(cv)))
})

test_that("the shifted CV denominator cannot reach zero", {
  # CV_Neut divides by (u_Neut + abs(min(u_Neut)) + 1); the shift is what keeps
  # a negative mean from producing a division by zero.
  biorep <- neut %>%
    filter(SampleID != "837SAnti-H3 1:250") %>%
    group_by(SampleID, Virus, Batch) %>%
    summarise(u_Neut = mean(Neut), .groups = "drop")

  denominator <- biorep$u_Neut + abs(min(biorep$u_Neut)) + 1
  expect_true(all(denominator >= 1))
})
