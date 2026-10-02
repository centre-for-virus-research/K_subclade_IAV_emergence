# ==============================================================================
# 01_neutralisation_analysis.R -- beta regression on 1:75 neutralisation
#
# The script builds its model forward by AIC and ends at model_2ae, from which
# every published marginal mean and contrast is taken. These tests refit the
# key models on the real data and check the properties the script depends on:
# that the response is admissible for betareg, that the chosen model is
# estimable, and that update() puts added interactions in the mean submodel
# rather than the precision submodel.
# ==============================================================================

suppressMessages({
  library(dplyr)
  library(betareg)
})

all_data <- build_data_all()

# Fitting is the slow part, so the two models the tests share are fitted once.
# betareg warns that it could not find a starting value for the precision
# parameter and falls back to 1; the fit converges regardless, and the script
# itself emits the same warning.
minimal_model <- suppressWarnings(
  betareg(u_Neut ~ Batch + Virus | Batch + Virus, data = all_data)
)

forward_model <- suppressWarnings(betareg(
  u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group +
    Virus * vacc_group +
    Age_policy * vacc_group |
    Batch + Virus,
  data = all_data
))

# ------------------------------------------------------------------------------
# Admissibility of the response
# ------------------------------------------------------------------------------

test_that("the analysis frame is large enough to fit the forward model", {
  expect_gt(nrow(all_data), 100)
  expect_gt(df.residual(forward_model), 0)
})

test_that("betareg accepts the Smithson & Verkuilen transformed response", {
  expect_s3_class(minimal_model, "betareg")
  expect_true(all(is.finite(coef(minimal_model))))
})

test_that("an untransformed response containing 0 or 1 would be rejected", {
  frame <- data.frame(y = c(0, 0.3, 0.6, 1, 0.5, 0.2), g = factor(rep(c("a", "b"), 3)))
  expect_error(betareg(y ~ g, data = frame))
})

# ------------------------------------------------------------------------------
# The models the script selects
# ------------------------------------------------------------------------------

test_that("every coefficient of the forward model is estimable", {
  expect_equal(sum(is.na(coef(forward_model))), 0)
})

test_that("the forward model has both a mean and a precision submodel", {
  mean_terms <- attr(terms(forward_model, model = "mean"), "term.labels")
  precision_terms <- attr(terms(forward_model, model = "precision"), "term.labels")

  expect_true(all(c("Batch", "Virus", "Age_policy", "Sex", "vacc_group") %in% mean_terms))
  expect_true("Virus:vacc_group" %in% mean_terms)
  expect_true("Age_policy:vacc_group" %in% mean_terms)
  expect_setequal(precision_terms, c("Batch", "Virus"))
})

test_that("the forward model beats the batch-and-virus-only model on AIC", {
  go_compare <- extract_object("01_neutralisation_analysis.R", "go_compare")
  comparison <- go_compare(forward_model, minimal_model)

  expect_equal(comparison$AIC.verdict[1], "Strongly Favour")
})

test_that("fitted values stay inside the open unit interval", {
  fitted_values <- fitted(forward_model)

  expect_true(all(fitted_values > 0))
  expect_true(all(fitted_values < 1))
  expect_equal(length(fitted_values), nrow(all_data))
})

test_that("quantile residuals are finite and roughly centred", {
  residual_values <- residuals(forward_model, type = "quantile")

  expect_true(all(is.finite(residual_values)))
  expect_lt(abs(mean(residual_values)), 0.25)
})

# ------------------------------------------------------------------------------
# update() on a two-part formula
# ------------------------------------------------------------------------------

test_that("update() adds an interaction to the mean submodel, not the precision one", {
  # Figures 4B and 4C depend on this: the script adds Virus:Sex and
  # Virus:Age_policy via update(). Because `|` binds more loosely than `+` in a
  # formula, an added term could plausibly have landed in the precision model
  # instead, which would silently change what those contrasts mean.
  updated <- update(forward_model, . ~ . + Virus:Sex)

  mean_terms <- attr(terms(updated, model = "mean"), "term.labels")
  precision_terms <- attr(terms(updated, model = "precision"), "term.labels")

  expect_true("Virus:Sex" %in% mean_terms)
  expect_false("Virus:Sex" %in% precision_terms)
  expect_setequal(precision_terms, c("Batch", "Virus"))
})

test_that("update() leaves the original model untouched", {
  before <- attr(terms(forward_model, model = "mean"), "term.labels")
  invisible(update(forward_model, . ~ . + Virus:Sex))
  after <- attr(terms(forward_model, model = "mean"), "term.labels")

  expect_equal(before, after)
})

test_that("the age interaction update behaves the same way", {
  updated <- update(forward_model, . ~ . + Virus:Age_policy)

  expect_true("Virus:Age_policy" %in% attr(terms(updated, model = "mean"), "term.labels"))
  expect_setequal(attr(terms(updated, model = "precision"), "term.labels"), c("Batch", "Virus"))
})

test_that("an updated model adds parameters rather than replacing them", {
  updated <- update(forward_model, . ~ . + Virus:Sex)

  expect_gt(length(coef(updated)), length(coef(forward_model)))
  expect_equal(sum(is.na(coef(updated))), 0)
})

# ------------------------------------------------------------------------------
# Marginal means
# ------------------------------------------------------------------------------

test_that("emmeans returns one response-scale estimate per virus", {
  skip_if_not_installed("emmeans")
  suppressMessages(library(emmeans))

  emm <- as.data.frame(suppressMessages(
    emmeans(forward_model, ~ Virus, type = "response", weights = "equal")
  ))

  expect_equal(nrow(emm), 3)
  expect_setequal(as.character(emm$Virus), c("THA22", "ENG24", "ENG25"))
})

test_that("response-scale marginal means are proportions", {
  skip_if_not_installed("emmeans")
  suppressMessages(library(emmeans))

  emm <- as.data.frame(suppressMessages(
    emmeans(forward_model, ~ Virus, type = "response", weights = "equal")
  ))
  estimate <- emm[[2]]

  expect_true(all(estimate > 0 & estimate < 1))
})

test_that("nested marginal means cover every virus by vaccination cell", {
  skip_if_not_installed("emmeans")
  suppressMessages(library(emmeans))

  emm <- as.data.frame(suppressMessages(
    emmeans(forward_model, ~ Virus | vacc_group, type = "response", weights = "equal")
  ))

  expect_equal(nrow(emm), 3 * 3)
  expect_true(all(is.finite(emm$SE)))
})

test_that("pairwise contrasts are produced for every virus pair", {
  skip_if_not_installed("emmeans")
  suppressMessages(library(emmeans))

  emm <- suppressMessages(emmeans(forward_model, ~ Virus, type = "response", weights = "equal"))
  contrasts <- as.data.frame(contrast(emm, method = "pairwise"))

  expect_equal(nrow(contrasts), 3)
  expect_true(all(contrasts$p.value >= 0 & contrasts$p.value <= 1))
})

test_that("the 2022 vaccine strain is neutralised best of the three", {
  # Direction check on the headline result: THA22 is the vaccine-matched strain,
  # the 2024 and 2025 isolates are the drifted ones.
  skip_if_not_installed("emmeans")
  suppressMessages(library(emmeans))

  emm <- as.data.frame(suppressMessages(
    emmeans(forward_model, ~ Virus, type = "response", weights = "equal")
  ))
  rownames(emm) <- as.character(emm$Virus)
  estimate <- setNames(emm[[2]], rownames(emm))

  expect_gt(estimate[["THA22"]], estimate[["ENG24"]])
  expect_gt(estimate[["THA22"]], estimate[["ENG25"]])
})

# ------------------------------------------------------------------------------
# Figure output
# ------------------------------------------------------------------------------

test_that("the script creates its figures directory before saving", {
  source_lines <- script_source("01_neutralisation_analysis.R")

  creates <- grep('dir\\.create\\(here\\("figures"\\)', source_lines)
  first_save <- grep('ggsave\\(here\\("figures"', source_lines)

  expect_gt(length(creates), 0)
  expect_gt(length(first_save), 0)
  expect_lt(min(creates), min(first_save))
})

test_that("ggsave refuses a missing directory, which is why that matters", {
  skip_if_not_installed("ggplot2")
  suppressMessages(library(ggplot2))

  target <- file.path(tempfile("missing-dir-"), "plot.svg")
  plot <- ggplot(data.frame(x = 1:3, y = 1:3)) + geom_point(aes(x, y))

  expect_error(ggsave(target, plot = plot, width = 3, height = 2))
})

test_that("a figure panel builds on the real data", {
  skip_if_not_installed("ggplot2")
  suppressMessages(library(ggplot2))

  counts <- count(ungroup(all_data), Virus)
  plot <- ggplot(all_data) +
    geom_boxplot(aes(x = Virus, y = u_Neut * 100, fill = Virus), alpha = 0.75, outlier.shape = NA) +
    geom_text(data = counts, aes(x = Virus, y = -7, label = n), inherit.aes = FALSE) +
    scale_fill_manual(values = c("#4677AA", "#DECC76", "#CC6476"))

  built <- ggplot_build(plot)
  expect_s3_class(built, "ggplot_built")
  expect_equal(nrow(counts), 3)
})

test_that("the figures written match the panels the paper reports", {
  source_lines <- script_source("01_neutralisation_analysis.R")

  for (panel in c("fig4A.svg", "fig4B.svg", "fig4C.svg", "fig4D.svg")) {
    expect_true(any(grepl(panel, source_lines, fixed = TRUE)), info = panel)
  }
})
