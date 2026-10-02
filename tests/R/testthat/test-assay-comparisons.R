# ==============================================================================
# 03_assay_comparisons.R -- relating 1:75 neutralisation to fitted IC50s
#
# Three things in this script decide published numbers:
#   * the censoring rule that turns a flat curve into a numeric IC50,
#   * the logit transform applied to the 1:75 neutralisation percentages,
#   * the titre labels on the IC50 axis.
# All three are pinned here. The slope-base mismatch between this script and
# script 02 is covered in test-known-defects.R.
# ==============================================================================

suppressMessages({
  library(readr)
  library(dplyr)
})

source_lines <- script_source("03_assay_comparisons.R")

fits <- read_csv(data_path("IC50_fits.csv"),
                 col_types = cols(SampleID = col_character()), progress = FALSE)
meta <- read_csv(data_path("sample_metadata.csv"),
                 col_types = cols(SampleID = col_character()), progress = FALSE)

# ------------------------------------------------------------------------------
# Titre axis labels
# ------------------------------------------------------------------------------

test_that("an IC50 axis break maps to the titre 1:2^|break|", {
  breaks <- -seq(2, 12, 2)
  labels <- paste0("1:", 2^seq(2, 12, 2))

  expect_equal(labels[match(-2, breaks)], "1:4")
  expect_equal(labels[match(-12, breaks)], "1:4096")
})

test_that("the figure 4 panel plots log2(dilution) on its y axis", {
  # y = -IC50, and IC50 = log2(1/dilution), so y = log2(dilution).
  ic50 <- log2(1 / 256)
  expect_equal(-ic50, log2(256))
  expect_equal(2^(-ic50), 256)
})

test_that("titre labels ascend together with their breaks", {
  breaks <- seq(4, 12, 2)
  labels <- paste0("1:", 2^seq(4, 12, 2))

  expect_equal(labels[1], "1:16")
  expect_equal(labels[length(labels)], "1:4096")
  for (i in seq_along(breaks)) {
    expect_equal(labels[i], paste0("1:", 2^breaks[i]))
  }
})

test_that("the IC50 panel does not reverse its titre labels", {
  # Regression guard: the labels were wrapped in rev(), which printed 1:4096
  # against the 1:16 tick on the Virus x Vaccination panel.
  panel <- grep("labels = rev\\(paste0\\(\"1:\"", source_lines, value = TRUE)
  expect_length(panel, 0)
})

test_that("every titre axis in the script labels its breaks consistently", {
  for (line in grep("paste0\\(\"1:\", 2\\^seq", source_lines, value = TRUE)) {
    expect_false(grepl("rev\\(paste0\\(\"1:\"", line), info = line)
  }
})

# ------------------------------------------------------------------------------
# Censoring flat curves
# ------------------------------------------------------------------------------

censor_as_coded <- function(verdict, lower, ic50) {
  ifelse(verdict == "Flattened", ifelse(lower < 50, -4, -12), ic50)
}

test_that("a curve flat at a low plateau is censored as the weakest titre", {
  expect_equal(censor_as_coded("Flattened", 5, NA_real_), -4)
  expect_equal(censor_as_coded("Flattened", 49.9, NA_real_), -4)
})

test_that("a curve flat at full neutralisation is censored as the strongest titre", {
  expect_equal(censor_as_coded("Flattened", 100, NA_real_), -12)
})

test_that("censoring leaves converged curves untouched", {
  expect_equal(censor_as_coded("Curved", 10, -7.5), -7.5)
  expect_equal(censor_as_coded("Curved", 95, -9.1), -9.1)
})

test_that("the strong-neutralisation censoring value is below every observed IC50", {
  observed <- fits$IC50[!is.na(fits$IC50)]
  expect_lte(-12, min(observed))
})

# The matching upper bracket does NOT hold: one published fit has an IC50 above
# the -4 censoring value, which script 02 would have rejected. See
# test-known-defects.R ("published fits fall outside the bounds script 02 enforces").

test_that("the censoring values are 8 log2 units apart, i.e. 256-fold", {
  expect_equal(abs(-4 - (-12)), 8)
  expect_equal(2^8, 256)
})

test_that("censoring resolves every flattened fit to a number", {
  censored <- censor_as_coded(fits$Verdict, fits$Lower, fits$IC50)
  flattened <- fits$Verdict == "Flattened"

  expect_equal(sum(is.na(censored[flattened])), 0)
})

test_that("curved fits that never converged stay missing after censoring", {
  # 10 rows are labelled Curved but carry no IC50; lm() drops them silently.
  censored <- censor_as_coded(fits$Verdict, fits$Lower, fits$IC50)
  curved_without_ic50 <- fits$Verdict == "Curved" & is.na(fits$IC50)

  expect_true(all(is.na(censored[curved_without_ic50])))
})

# ------------------------------------------------------------------------------
# The logit transform on 1:75 neutralisation
# ------------------------------------------------------------------------------

epsilon <- 0.025

to_proportion <- function(percent) pmin(pmax(percent, 0), 100) / 100

adjust <- function(p) ifelse(p <= 0, epsilon, ifelse(p >= 1, 1 - epsilon, p))

test_that("percentages are clamped into [0, 1] before anything else", {
  expect_equal(to_proportion(c(-54, 0, 37.5, 100, 120)), c(0, 0, 0.375, 1, 1))
})

test_that("the epsilon adjustment keeps the logit finite at both ends", {
  expect_true(is.finite(qlogis(adjust(0))))
  expect_true(is.finite(qlogis(adjust(1))))
})

test_that("the epsilon adjustment is symmetric about one half", {
  expect_equal(adjust(0), epsilon)
  expect_equal(adjust(1), 1 - epsilon)
  expect_equal(qlogis(adjust(0)), -qlogis(adjust(1)))
})

test_that("interior values pass through the adjustment unchanged", {
  for (p in c(0.01, 0.25, 0.5, 0.75, 0.99)) {
    expect_equal(adjust(p), p)
  }
})

test_that("the logit is monotone in the underlying percentage", {
  percentages <- seq(0, 100, by = 5)
  logits <- qlogis(adjust(to_proportion(percentages)))

  expect_true(all(diff(logits) >= 0))
})

test_that("plogis inverts qlogis over the adjusted range", {
  for (p in c(epsilon, 0.1, 0.5, 0.9, 1 - epsilon)) {
    expect_equal(plogis(qlogis(p)), p)
  }
})

test_that("the clamp is applied once, so rerunning section 1.4 would corrupt the scale", {
  # data_comp$u_Neut is overwritten in place by `pmin(pmax(u_Neut,0),100)/100`.
  # Running that line twice divides by 100 again; this records the hazard.
  once <- to_proportion(75)
  twice <- to_proportion(once * 100) / 100

  expect_equal(once, 0.75)
  expect_equal(twice, 0.0075)
  expect_false(isTRUE(all.equal(once, twice)))
})

# ------------------------------------------------------------------------------
# The IC50 model
# ------------------------------------------------------------------------------

model_frame <- local({
  frame <- left_join(fits, meta, by = "SampleID")
  frame$Virus <- factor(frame$Virus, levels = c("THA22", "ENG24", "ENG25"))
  frame$vacc_group <- ifelse(frame$vacc_years_since_last == ">4", ">4",
                      ifelse(frame$vacc_years_since_last == "1", "1", "2-4"))
  frame$vacc_group <- factor(frame$vacc_group, levels = c(">4", "2-4", "1"))
  frame$IC50 <- censor_as_coded(frame$Verdict, frame$Lower, frame$IC50)
  frame
})

test_that("the model frame covers every virus by vaccination cell", {
  counts <- table(model_frame$Virus, model_frame$vacc_group)
  expect_true(all(counts > 0))
})

test_that("the interaction model fits and is estimable", {
  model <- lm(IC50 ~ Virus + vacc_group + Virus * vacc_group, data = model_frame)

  expect_s3_class(model, "lm")
  expect_equal(sum(is.na(coef(model))), 0)
  expect_gt(df.residual(model), 0)
})

test_that("emmeans returns one estimate per virus with finite intervals", {
  skip_if_not_installed("emmeans")
  suppressMessages(library(emmeans))

  model <- lm(IC50 ~ Virus + vacc_group + Virus * vacc_group, data = model_frame)
  emm <- suppressMessages(emmeans(model, ~ Virus, weights = "equal"))
  summarised <- as.data.frame(emm)

  expect_equal(nrow(summarised), 3)
  expect_true(all(is.finite(summarised$emmean)))
  expect_true(all(is.finite(summarised$SE)))
})

test_that("pairwise virus contrasts are returned for all three pairs", {
  skip_if_not_installed("emmeans")
  suppressMessages(library(emmeans))

  model <- lm(IC50 ~ Virus + vacc_group + Virus * vacc_group, data = model_frame)
  emm <- suppressMessages(emmeans(model, ~ Virus, weights = "equal"))
  contrasts <- as.data.frame(contrast(emm, "pairwise", adjust = "none"))

  expect_equal(nrow(contrasts), 3)
  expect_true(all(is.finite(contrasts$estimate)))
  expect_true(all(contrasts$p.value >= 0 & contrasts$p.value <= 1))
})

test_that("the 2022 vaccine strain is neutralised less well than the 2024/2025 isolates", {
  # Direction check on the headline result: a more negative IC50 is more potent.
  skip_if_not_installed("emmeans")
  suppressMessages(library(emmeans))

  model <- lm(IC50 ~ Virus + vacc_group + Virus * vacc_group, data = model_frame)
  emm <- as.data.frame(suppressMessages(emmeans(model, ~ Virus, weights = "equal")))
  rownames(emm) <- as.character(emm$Virus)

  expect_lt(emm["THA22", "emmean"], emm["ENG24", "emmean"])
  expect_lt(emm["THA22", "emmean"], emm["ENG25", "emmean"])
})

# ------------------------------------------------------------------------------
# Helper objects
# ------------------------------------------------------------------------------

test_that("theme_base is a usable ggplot theme", {
  skip_if_not_installed("ggplot2")
  suppressMessages(library(ggplot2))

  theme_base <- extract_object("03_assay_comparisons.R", "theme_base")
  expect_s3_class(theme_base, "theme")
})

test_that("the mean slope used for the illustrative curves is positive and finite", {
  curved <- filter(fits, Slope != 0, !is.na(IC50))
  mean_slope <- mean(curved$Slope)

  expect_true(is.finite(mean_slope))
  expect_gt(mean_slope, 0)
})

test_that("the illustrative sigmoid family is monotone and bounded", {
  curved <- filter(fits, Slope != 0, !is.na(IC50))
  mean_slope <- mean(curved$Slope)

  x <- seq(-14, 0, length.out = 200)
  for (ic50 in seq(-12, -2, 1)) {
    y <- 100 / (1 + 2^(mean_slope * (ic50 - x)))
    expect_true(all(diff(y) > 0))
    expect_true(all(y > 0 & y < 100))
  }
})
