# ==============================================================================
# The forward model-building chain in 01_neutralisation_analysis.R
#
# Section 4 of that script picks model_2ae by comparing AIC at each step, and
# records the outcome of each comparison in a trailing comment. Those comments
# are the audit trail for why the published model is the one it is, so they are
# re-derived here from the real data rather than taken on trust.
#
# The other thing checked here is the precondition that makes any of it valid:
# AIC is only comparable across models fitted to the same observations, and
# betareg drops rows with a missing predictor silently. A single NA in one
# candidate predictor would shrink that model's data and make its AIC
# incomparable, with nothing in the console to say so.
# ==============================================================================

suppressMessages({
  library(dplyr)
  library(betareg)
})

all_data <- build_data_all()
all_data$vacc_number_of_doses <- factor(all_data$vacc_number_of_doses)

go_compare <- extract_object("01_neutralisation_analysis.R", "go_compare")

fit <- function(formula) suppressWarnings(betareg(formula, data = all_data))

# The candidate models the script actually compares, in the order it builds them.
candidates <- list(
  doses   = u_Neut ~ Batch + Virus + Age + Sex + vacc_number_of_doses | Batch + Virus,
  recency = u_Neut ~ Batch + Virus + Age + Sex + vacc_years_since_last | Batch + Virus,
  none    = u_Neut ~ Batch + Virus + Age + Sex | Batch + Virus,
  coarse  = u_Neut ~ Batch + Virus + Age + Sex + vacc_group | Batch + Virus,
  m1      = u_Neut ~ Batch + Virus | Batch + Virus,
  m1a     = u_Neut ~ Batch + Virus + Age_policy | Batch + Virus,
  m1b     = u_Neut ~ Batch + Virus + Sex | Batch + Virus,
  m1c     = u_Neut ~ Batch + Virus + vacc_group | Batch + Virus,
  m1ac    = u_Neut ~ Batch + Virus + Age_policy + vacc_group | Batch + Virus,
  m1abc   = u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group | Batch + Virus,
  m2a     = u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group + Virus * vacc_group | Batch + Virus,
  m2e     = u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group + Age_policy * vacc_group | Batch + Virus,
  m2ae    = u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group +
                     Virus * vacc_group + Age_policy * vacc_group | Batch + Virus,
  m3      = u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group +
                     Virus * vacc_group + Age_policy * vacc_group +
                     Virus * vacc_group * Age_policy | Batch + Virus
)

fits <- lapply(candidates, fit)

# ------------------------------------------------------------------------------
# The precondition for comparing AIC at all
# ------------------------------------------------------------------------------

test_that("every candidate model is fitted to the same observations", {
  counts <- vapply(fits, nobs, numeric(1))

  expect_equal(length(unique(counts)), 1L,
               info = paste("differing nobs:", paste(names(counts), counts, collapse = ", ")))
  expect_equal(unname(counts[["m1"]]), nrow(all_data))
})

test_that("no candidate predictor contains a missing value", {
  # This is what keeps the counts above equal; betareg would drop such rows
  # without a message and quietly break every AIC comparison downstream.
  for (column in c("Age", "Age_policy", "Sex", "vacc_group", "vacc_years_since_last",
                   "vacc_number_of_doses", "Batch", "Virus", "u_Neut")) {
    expect_equal(sum(is.na(all_data[[column]])), 0,
                 info = sprintf("%s contains NA and would shrink the model frame", column))
  }
})

test_that("comparing models on different data warns but still returns a verdict", {
  # Documented hazard rather than a defect in this script. stats::AIC() warns
  # when the models were fitted to different numbers of observations, but
  # go_compare() passes the warning straight through and still prints a
  # confident verdict, so the equal-nobs property above has to be maintained by
  # the caller rather than relied on being caught here.
  set.seed(11)
  big <- data.frame(x = rnorm(100))
  big$y <- big$x + rnorm(100, sd = 0.5)
  small <- big[1:40, ]

  wide <- lm(y ~ x, data = big)
  narrow <- lm(y ~ x, data = small)
  expect_false(nobs(wide) == nobs(narrow))

  expect_warning(
    comparison <- go_compare(wide, narrow),
    "not all fitted to the same number of observations"
  )

  # The verdict is produced regardless, which is the part worth knowing about.
  expect_equal(nrow(comparison), 2)
  expect_true(all(comparison$AIC.verdict %in%
                  c("Strongly Favour", "Weakly Favour", "Equivalent",
                    "Weakly Disfavour", "Strongly Disfavour")))
})

test_that("the real model comparisons produce no such warning", {
  expect_no_warning(go_compare(fits$m2ae, fits$m3))
  expect_no_warning(go_compare(fits$m1, fits$m1abc))
})

# ------------------------------------------------------------------------------
# The verdicts recorded in the script's comments
# ------------------------------------------------------------------------------

expect_verdict <- function(first, second, expected) {
  verdict <- go_compare(fits[[first]], fits[[second]])$AIC.verdict[1]
  expect_match(verdict, expected,
               info = sprintf("go_compare(%s, %s) gave '%s'", first, second, verdict))
}

test_that("vaccination recency is preferred over dose count, as recorded", {
  # Script: go_compare(model_doses, model_recency)
  expect_verdict("doses", "recency", "Disfavour")
})

test_that("each first-order term improves on the batch-and-virus model, as recorded", {
  # Script: "1a strongly favoured", "1b strongly favoured", "1c strongly favoured"
  for (candidate in c("m1a", "m1b", "m1c")) {
    expect_verdict("m1", candidate, "Disfavour")
  }
})

test_that("age and vaccination together beat either alone, as recorded", {
  # Script: "1ac strongly favoured" over both 1a and 1c
  expect_verdict("m1a", "m1ac", "Disfavour")
  expect_verdict("m1c", "m1ac", "Disfavour")
})

test_that("adding sex is a further improvement, as recorded", {
  # Script: "1abc weakly favoured"
  expect_verdict("m1ac", "m1abc", "Disfavour")
})

test_that("the virus-by-vaccination interaction is favoured, as recorded", {
  # Script: go_compare(model_2, model_2a) # weakly favoured
  expect_verdict("m1abc", "m2a", "Disfavour")
})

test_that("the two-interaction model beats each single-interaction model, as recorded", {
  # Script: "2ae weakly favoured" over both 2a and 2e
  expect_verdict("m2a", "m2ae", "Disfavour")
  expect_verdict("m2e", "m2ae", "Disfavour")
})

test_that("the three-way interaction is rejected, as recorded", {
  # Script: go_compare(model_2ae, model_3) # 2ae strongly favoured
  expect_verdict("m2ae", "m3", "Strongly Favour")
})

test_that("the selected model is the one with the lowest AIC of the chain", {
  aic_values <- vapply(fits, AIC, numeric(1))
  chain <- c("m1", "m1a", "m1b", "m1c", "m1ac", "m1abc", "m2a", "m2e", "m2ae", "m3")

  expect_equal(names(which.min(aic_values[chain])), "m2ae")
})

test_that("model_2ae is better than its own starting point by a wide margin", {
  expect_lt(AIC(fits$m2ae), AIC(fits$m1) - 6)
})

# ------------------------------------------------------------------------------
# Age binning boundaries
# ------------------------------------------------------------------------------

test_that("Age_policy bins are closed on the left, as right = FALSE specifies", {
  binned <- function(age) {
    as.character(cut(age, breaks = c(18, 50, 65, Inf), right = FALSE,
                     labels = c("18-49", "50-64", "65+")))
  }

  expect_equal(binned(18), "18-49")
  expect_equal(binned(49), "18-49")
  expect_equal(binned(49.9), "18-49")
  expect_equal(binned(50), "50-64")
  expect_equal(binned(64), "50-64")
  expect_equal(binned(64.9), "50-64")
  expect_equal(binned(65), "65+")
  expect_equal(binned(120), "65+")
})

test_that("an age below the first break is dropped rather than binned", {
  expect_true(is.na(cut(17, breaks = c(18, 50, 65, Inf), right = FALSE)))
})

test_that("the policy bins reproduce the UK vaccination age thresholds", {
  # 65 is the eligibility threshold the grouping is named for.
  expect_equal(nlevels(cut(all_data$Age, breaks = c(18, 50, 65, Inf), right = FALSE)), 3)
  expect_equal(
    sum(all_data$Age >= 65),
    sum(cut(all_data$Age, breaks = c(18, 50, 65, Inf), right = FALSE,
            labels = c("18-49", "50-64", "65+")) == "65+")
  )
})

test_that("every age bin the models use is populated", {
  # An empty level would make betareg drop a coefficient to NA.
  for (breaks in list(c(18, 50, 65, Inf), c(18, 40, 65, Inf), c(18, 40, 60, 80, Inf))) {
    counts <- table(cut(all_data$Age, breaks = breaks, right = FALSE))
    expect_true(all(counts > 0), info = paste(breaks, collapse = ","))
  }
})
