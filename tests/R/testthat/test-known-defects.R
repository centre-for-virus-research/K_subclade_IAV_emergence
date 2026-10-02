# ==============================================================================
# Known defects
#
# Each test here asserts that a CONFIRMED defect is still present. They are
# deliberately phrased that way so the suite stays green while the defects are
# visible and counted: when one is fixed, its test fails and should be rewritten
# as a normal regression test (or deleted along with the workaround).
#
# These were left unfixed because correcting them changes published numbers, so
# the call belongs to the authors rather than to the test suite.
# ==============================================================================

suppressMessages({
  library(readr)
  library(dplyr)
})

fits <- read_csv(data_path("IC50_fits.csv"),
                 col_types = cols(SampleID = col_character()), progress = FALSE)
meta <- read_csv(data_path("sample_metadata.csv"),
                 col_types = cols(SampleID = col_character()), progress = FALSE)
ic50s <- read_csv(data_path("IC50s.csv"),
                  col_types = cols(Sample = col_character()), progress = FALSE)

ic50_constants <- extract_numeric_constants("02_IC50_fitting.R")

# ------------------------------------------------------------------------------
# DEFECT 1
# 03_assay_comparisons.R: `ifelse(Lower < 50, -4, -12)`
#
# Script 02 produces a "flatten high" curve by setting Lower = Upper = 100, and
# a "flatten low" curve by setting both to the mean of the observations. The
# threshold test in script 03 therefore only matches the intent when the mean of
# a flat-low curve happens to be below 50. Two samples plateau just above 50 and
# are consequently recorded as the most potent possible titre instead of the
# least potent -- an 8 log2 unit (256-fold) error in the wrong direction.
# ------------------------------------------------------------------------------

test_that("DEFECT: no published flattened fit uses the flatten-high sentinel", {
  flattened <- filter(fits, Verdict == "Flattened")

  # If this ever becomes non-zero the <50 test starts behaving as intended for
  # those rows, and the defect below changes shape.
  expect_equal(sum(flattened$Lower == 100), 0)
})

test_that("DEFECT: flat-low curves plateauing above 50% are censored as maximally potent", {
  flattened <- filter(fits, Verdict == "Flattened")
  miscoded <- filter(flattened, Lower >= 50, Lower != 100)

  expect_equal(nrow(miscoded), 2)
  expect_setequal(miscoded$SampleID, c("837S1233", "837S1392"))
  expect_true(all(miscoded$Lower > 50 & miscoded$Lower < 53))

  # As coded they become -12 (1:4096); by intent they are -4 (1:16).
  as_coded <- ifelse(miscoded$Lower < 50, -4, -12)
  by_intent <- ifelse(miscoded$Lower == 100, -12, -4)

  expect_equal(as_coded, c(-12, -12))
  expect_equal(by_intent, c(-4, -4))
})

test_that("DEFECT: the miscoding moves the virus estimates, without flipping conclusions", {
  skip_if_not_installed("emmeans")
  suppressMessages(library(emmeans))

  build <- function(rule) {
    frame <- left_join(fits, meta, by = "SampleID")
    frame$Virus <- factor(frame$Virus, levels = c("THA22", "ENG24", "ENG25"))
    frame$vacc_group <- ifelse(frame$vacc_years_since_last == ">4", ">4",
                        ifelse(frame$vacc_years_since_last == "1", "1", "2-4"))
    frame$vacc_group <- factor(frame$vacc_group, levels = c(">4", "2-4", "1"))
    frame$IC50 <- if (rule == "as_coded") {
      ifelse(frame$Verdict == "Flattened", ifelse(frame$Lower < 50, -4, -12), frame$IC50)
    } else {
      ifelse(frame$Verdict == "Flattened", ifelse(frame$Lower == 100, -12, -4), frame$IC50)
    }
    as.data.frame(suppressMessages(
      emmeans(lm(IC50 ~ Virus * vacc_group, data = frame), ~ Virus, weights = "equal")
    ))
  }

  as_coded <- build("as_coded")
  by_intent <- build("intended")

  # The estimates genuinely differ ...
  expect_false(isTRUE(all.equal(as_coded$emmean, by_intent$emmean)))
  # ... but by less than a tenth of a log2 unit, so the published direction holds.
  expect_lt(max(abs(as_coded$emmean - by_intent$emmean)), 0.1)
})

# ------------------------------------------------------------------------------
# DEFECT 2
# data/IC50s.csv: the stored DilutionLog2 column is reversed for 306 rows
#
# For 14 sample x virus series the stored column runs backwards: the most
# concentrated well carries the value belonging to the most dilute one. Script
# 02 recomputes the column from Dilution, so the analysis is unaffected, but the
# published CSV is wrong for anyone reading it directly.
# ------------------------------------------------------------------------------

test_that("DEFECT: the stored DilutionLog2 column disagrees with Dilution for some rows", {
  mismatch <- abs(ic50s$DilutionLog2 - log2(1 / ic50s$Dilution)) > 1e-9

  expect_equal(sum(mismatch), 306)
  expect_equal(length(unique(ic50s$Sample[mismatch])), 14)
  expect_setequal(unique(ic50s$Virus[mismatch]), c("ENG24", "ENG25"))
})

test_that("DEFECT: the mismatched rows carry an exactly reversed dilution series", {
  ladder <- sort(unique(ic50s$Dilution))
  reversed_value <- log2(1 / rev(ladder))[match(ic50s$Dilution, ladder)]
  mismatch <- abs(ic50s$DilutionLog2 - log2(1 / ic50s$Dilution)) > 1e-9

  expect_true(all(abs(ic50s$DilutionLog2[mismatch] - reversed_value[mismatch]) < 1e-9))
})

test_that("the analysis is shielded because script 02 recomputes the column", {
  expect_true(any(grepl("DilutionLog2\\s*<-\\s*log2\\(1/", script_source("02_IC50_fitting.R"))))
})

# ------------------------------------------------------------------------------
# DEFECT 3
# data/IC50_fits.csv was not produced by the current 02_IC50_fitting.R
#
# The shipped fits contain values the current script's bounds and acceptance
# window would have rejected, so the file cannot be regenerated from the script
# as it stands.
# ------------------------------------------------------------------------------

test_that("DEFECT: published slopes exceed the SlopeMax the script enforces", {
  above <- sum(fits$Slope > ic50_constants$SlopeMax, na.rm = TRUE)

  expect_equal(ic50_constants$SlopeMax, 2)
  expect_equal(above, 11)
  expect_gt(max(fits$Slope, na.rm = TRUE), ic50_constants$SlopeMax)
})

test_that("DEFECT: a published IC50 lies outside the acceptance window", {
  outside <- sum(fits$IC50 > ic50_constants$IC50Max, na.rm = TRUE)

  expect_equal(ic50_constants$IC50Max, -4)
  expect_equal(outside, 1)
})

test_that("DEFECT: some fits are labelled Curved but carry no IC50", {
  # Script 02 refuses the "curved" verdict when the fit did not converge, so
  # these rows cannot have come from the loop as written. lm() drops them.
  expect_equal(sum(fits$Verdict == "Curved" & is.na(fits$IC50)), 10)
})

# ------------------------------------------------------------------------------
# DEFECT 4
# Slope is fitted on the natural-log scale but reused on the base-2 scale
#
# Script 02 fits plogis((x - IC50) * Slope), i.e. a natural-log logistic in
# log2-dilution units. Script 03 then feeds those slopes into base-2
# expressions, which stretches the curve by a factor of log(2) and inflates the
# fold-change returned by log2_shift() by 1/log(2) ~ 1.44.
# ------------------------------------------------------------------------------

test_that("DEFECT: script 03 uses a base-2 logistic for slopes fitted in base e", {
  expect_true(any(grepl("100 / \\(1 \\+ 2\\^\\(mean_Slope", script_source("03_assay_comparisons.R"))))
  expect_true(any(grepl("plogis\\(\\(DilutionLog2", script_source("02_IC50_fitting.R"))))
})

test_that("DEFECT: the two parameterisations differ by a factor of log(2)", {
  slope <- 1.2
  ic50 <- -8
  x <- seq(-12, -4, length.out = 50)

  as_fitted <- plogis((x - ic50) * slope)
  as_plotted <- 1 / (1 + 2^(slope * (ic50 - x)))

  expect_false(isTRUE(all.equal(as_fitted, as_plotted)))
  expect_equal(as_plotted, plogis((x - ic50) * slope * log(2)))
})

test_that("DEFECT: log2_shift overstates the dilution shift by 1/log(2)", {
  mean_Slope <- 1.185
  log2_shift <- extract_object("03_assay_comparisons.R", "log2_shift",
                               env_vars = list(mean_Slope = mean_Slope))

  p1 <- 0.742
  p2 <- 0.538

  as_written <- log2_shift(p1, p2)
  # Consistent with the fitted model, the shift in log2-dilution units is the
  # natural-logit difference divided by the slope.
  consistent <- abs((qlogis(p1) - qlogis(p2)) / mean_Slope)

  expect_equal(as_written / consistent, 1 / log(2), tolerance = 1e-9)
  expect_gt(as_written, consistent)
  expect_equal(2^as_written / 2^consistent, 2^(consistent * (1 / log(2) - 1)), tolerance = 1e-6)
})

# ------------------------------------------------------------------------------
# DEFECT 5
# 02_IC50_fitting.R never saves its output
#
# The loop builds curve_parameters by hand classification, but nothing writes it
# to disk, so data/IC50_fits.csv cannot be reproduced by running the script.
# ------------------------------------------------------------------------------

test_that("DEFECT: script 02 produces curve_parameters but never writes it out", {
  source_lines <- script_source("02_IC50_fitting.R")

  expect_true(any(grepl("curve_parameters\\s*<-", source_lines)))
  expect_length(grep("write_csv|write\\.csv|saveRDS", source_lines), 0)
})

test_that("the file that script 02 should have written does exist in the repository", {
  expect_true(file.exists(data_path("IC50_fits.csv")))
})
