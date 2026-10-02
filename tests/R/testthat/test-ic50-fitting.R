# ==============================================================================
# 02_IC50_fitting.R -- four-parameter logistic fitting
#
# The script is interactive (it calls readline() per curve), so it cannot be
# sourced. What is testable, and what matters, is the model it fits: a logistic
# re-parameterised so that the curve passes exactly through the most concentrated
# observation, with IC50 recovered afterwards from the fitted asymptotes.
#
# These tests pin that re-parameterisation, the bound construction around it,
# and the fact that the curve drawn in section 1.3 is the same curve that was
# fitted in section 1.2.
# ==============================================================================

suppressMessages({
  library(readr)
  library(dplyr)
  library(minpack.lm)
})

constants <- extract_numeric_constants("02_IC50_fitting.R")
source_lines <- script_source("02_IC50_fitting.R")

# ------------------------------------------------------------------------------
# Parameter bounds declared by the script
# ------------------------------------------------------------------------------

test_that("the script declares every bound the fit relies on", {
  for (name in c("BottomMin", "BottomMax", "TopMin", "TopMax", "IC50Min", "IC50Max", "SlopeMin", "SlopeMax")) {
    expect_true(name %in% names(constants), info = sprintf("%s is no longer a top-level constant", name))
  }
})

test_that("each bound pair is ordered low-to-high", {
  expect_lt(constants$BottomMin, constants$BottomMax)
  expect_lt(constants$TopMin, constants$TopMax)
  expect_lt(constants$IC50Min, constants$IC50Max)
  expect_lt(constants$SlopeMin, constants$SlopeMax)
})

test_that("the bottom and top asymptote windows do not overlap", {
  expect_lte(constants$BottomMax, constants$TopMin)
})

test_that("asymptote bounds straddle the neutralisation percentage scale", {
  expect_lte(constants$BottomMin, 0)
  expect_gte(constants$TopMax, 100)
})

test_that("the slope is constrained to be positive", {
  # A negative slope would describe a curve that neutralises better when more
  # dilute, which is not a meaningful fit.
  expect_gt(constants$SlopeMin, 0)
})

test_that("the IC50 acceptance window lies in the measurable dilution range", {
  # DilutionLog2 = log2(1/Dilution), so the assayed range is about -11.6 to -5.6.
  expect_lt(constants$IC50Max, 0)
  expect_gt(constants$IC50Min, -30)
})

# ------------------------------------------------------------------------------
# The DilutionLog2 transform
# ------------------------------------------------------------------------------

test_that("the script recomputes DilutionLog2 rather than trusting the column", {
  expect_true(any(grepl("DilutionLog2\\s*<-\\s*log2\\(1/", source_lines)),
              info = "script 02 must not rely on the stored DilutionLog2 column")
})

test_that("log2(1/dilution) is monotone decreasing in dilution", {
  dilutions <- c(50, 100, 200, 400, 800, 1600, 3200)
  transformed <- log2(1 / dilutions)

  expect_true(all(diff(transformed) < 0))
  expect_equal(transformed[1], -log2(50))
})

test_that("the most concentrated well has the largest DilutionLog2", {
  dilutions <- c(50, 100, 200, 400, 800, 1600, 3200)
  expect_equal(which.max(log2(1 / dilutions)), 1L)
})

# ------------------------------------------------------------------------------
# The anchored re-parameterisation
# ------------------------------------------------------------------------------

anchored_curve <- function(x, Bottom, Top, Slope, upperX, upperY) {
  ic50 <- upperX - qlogis((upperY - Bottom) / (Top - Bottom)) / Slope
  Bottom + (Top - Bottom) * plogis((x - ic50) * Slope)
}

implied_ic50 <- function(Bottom, Top, Slope, upperX, upperY) {
  upperX - qlogis((upperY - Bottom) / (Top - Bottom)) / Slope
}

test_that("the fitted curve passes exactly through the anchor point", {
  # This is the whole point of the re-parameterisation: whatever the asymptotes
  # and slope turn out to be, the curve is pinned to the most concentrated
  # observation.
  for (Bottom in c(-50, -10, 0, 20, 39)) {
    for (Top in c(60, 80, 100)) {
      for (Slope in c(0.6, 1, 2)) {
        upperX <- -5.64
        upperY <- 55
        expect_equal(anchored_curve(upperX, Bottom, Top, Slope, upperX, upperY), upperY,
                     info = sprintf("Bottom=%g Top=%g Slope=%g", Bottom, Top, Slope))
      }
    }
  }
})

test_that("the implied IC50 is where the curve reaches the midpoint of its asymptotes", {
  Bottom <- 5; Top <- 95; Slope <- 1.2; upperX <- -5.64; upperY <- 70
  ic50 <- implied_ic50(Bottom, Top, Slope, upperX, upperY)

  expect_equal(anchored_curve(ic50, Bottom, Top, Slope, upperX, upperY), (Bottom + Top) / 2)
})

test_that("the curve is monotone increasing in DilutionLog2 for a positive slope", {
  x <- seq(-12, -5, length.out = 100)
  y <- anchored_curve(x, Bottom = 0, Top = 100, Slope = 1, upperX = -5.64, upperY = 80)

  expect_true(all(diff(y) > 0))
})

test_that("the curve never leaves its asymptotes", {
  x <- seq(-40, 20, length.out = 500)
  y <- anchored_curve(x, Bottom = -20, Top = 90, Slope = 1.5, upperX = -5.64, upperY = 60)

  # Inclusive: far from the midpoint plogis() underflows to exactly 0 or 1.
  expect_true(all(y >= -20))
  expect_true(all(y <= 90))
})

test_that("the curve stays strictly inside its asymptotes across the assayed range", {
  x <- log2(1 / c(50, 100, 200, 400, 800, 1600, 3200))
  y <- anchored_curve(x, Bottom = -20, Top = 90, Slope = 1.5, upperX = -5.64, upperY = 60)

  expect_true(all(y > -20))
  expect_true(all(y < 90))
})

test_that("section 1.3 redraws the same curve that section 1.2 fitted", {
  # 1.2 fits  Bottom + (Top-Bottom) * plogis((x - IC50) * Slope)
  # 1.3 draws Bottom + (Top-Bottom) / (1 + exp(-(x - IC50) * Slope))
  Bottom <- 10; Top <- 90; Slope <- 1.3; IC50 <- -8
  x <- seq(-12, -5, length.out = 50)

  fitted_form <- Bottom + (Top - Bottom) * plogis((x - IC50) * Slope)
  drawn_form <- Bottom + (Top - Bottom) / (1 + exp(-(x - IC50) * Slope))

  expect_equal(fitted_form, drawn_form)
})

test_that("the logistic is on the natural-log scale, not base 2", {
  # Worth pinning explicitly: script 03 reuses these slopes in a base-2
  # expression, and the two differ by a factor of log(2).
  Slope <- 1; IC50 <- -8; x <- -7

  natural <- plogis((x - IC50) * Slope)
  base_two <- 1 / (1 + 2^(Slope * (IC50 - x)))

  expect_false(isTRUE(all.equal(natural, base_two)))
  expect_equal(base_two, plogis((x - IC50) * Slope * log(2)))
})

# ------------------------------------------------------------------------------
# Bound construction per curve
# ------------------------------------------------------------------------------

feasibility <- function(upperY, epsilon = 0.001) {
  bottom_max <- min(constants$BottomMax, upperY - epsilon)
  top_min <- max(constants$TopMin, upperY + epsilon)
  list(bottom_max = bottom_max,
       top_min = top_min,
       feasible = bottom_max > constants$BottomMin && top_min < constants$TopMax)
}

test_that("the per-curve bounds keep the anchor strictly between the asymptotes", {
  for (upperY in c(-10, 0, 25, 50, 75, 95, 99)) {
    bounds <- feasibility(upperY)
    expect_lt(bounds$bottom_max, upperY)
    expect_gt(bounds$top_min, upperY)
  }
})

test_that("a typical anchor is feasible", {
  expect_true(feasibility(80)$feasible)
  expect_true(feasibility(50)$feasible)
  expect_true(feasibility(20)$feasible)
})

test_that("an anchor at full neutralisation makes the fit infeasible by construction", {
  # TopMin_i becomes 100.001, which exceeds TopMax, so no logistic is attempted
  # and the curve must be classified by hand instead.
  expect_false(feasibility(100)$feasible)
  expect_gte(feasibility(100)$top_min, constants$TopMax)
})

test_that("an anchor at or below the bottom bound makes the fit infeasible", {
  expect_false(feasibility(constants$BottomMin)$feasible)
  expect_false(feasibility(constants$BottomMin - 10)$feasible)
})

test_that("the feasibility rule matches the real dilution data", {
  d <- read_csv(data_path("IC50s.csv"), col_types = cols(Sample = col_character()), progress = FALSE)
  d$DilutionLog2 <- log2(1 / d$Dilution)

  anchors <- d %>%
    group_by(Batch, Virus, Sample, DilutionLog2) %>%
    summarise(u_Neut = mean(Neut), .groups = "drop") %>%
    group_by(Sample, Virus) %>%
    slice_max(DilutionLog2, n = 1, with_ties = FALSE) %>%
    ungroup()

  infeasible <- sum(!vapply(anchors$u_Neut, function(y) feasibility(y)$feasible, logical(1)))

  # Only curves saturating at 100% are excluded; if this grows, the hand
  # classification burden in section 1.3 has grown with it.
  expect_lt(infeasible, 0.05 * nrow(anchors))
})

# ------------------------------------------------------------------------------
# Round trip through nlsLM
# ------------------------------------------------------------------------------

test_that("nlsLM recovers the generating parameters from clean simulated data", {
  true_bottom <- 5
  true_top <- 95
  true_slope <- 1.1
  true_ic50 <- -8.5

  x <- log2(1 / c(50, 100, 200, 400, 800, 1600, 3200))
  y <- true_bottom + (true_top - true_bottom) * plogis((x - true_ic50) * true_slope)
  frame <- data.frame(DilutionLog2 = x, u_Neut = y)

  upper_i <- which.max(frame$DilutionLog2)
  upperX <- frame$DilutionLog2[upper_i]
  upperY <- frame$u_Neut[upper_i]
  bounds <- feasibility(upperY)

  # Start values exactly as the script computes them: the clamping is what keeps
  # the starting point inside the per-curve bounds.
  start_bottom <- pmin(pmax(min(frame$u_Neut), constants$BottomMin), bounds$bottom_max)
  start_top <- pmin(pmax(max(frame$u_Neut), bounds$top_min), constants$TopMax)

  fit <- nlsLM(
    u_Neut ~ Bottom + (Top - Bottom) *
      plogis((DilutionLog2 - (upperX - qlogis((upperY - Bottom) / (Top - Bottom)) / Slope)) * Slope),
    data = frame,
    start = list(Bottom = start_bottom, Top = start_top, Slope = 1),
    lower = c(Bottom = constants$BottomMin, Top = bounds$top_min, Slope = constants$SlopeMin),
    upper = c(Bottom = bounds$bottom_max, Top = constants$TopMax, Slope = constants$SlopeMax),
    control = nls.lm.control(maxiter = 1000)
  )

  coefs <- coef(fit)
  recovered <- implied_ic50(coefs[["Bottom"]], coefs[["Top"]], coefs[["Slope"]], upperX, upperY)

  expect_equal(unname(coefs[["Slope"]]), true_slope, tolerance = 1e-3)
  expect_equal(recovered, true_ic50, tolerance = 1e-3)
})

test_that("the script's start values are always inside the per-curve bounds", {
  # nlsLM produces NaNs and parks on a bound when started outside the feasible
  # box, so this clamping is load-bearing rather than cosmetic.
  for (upperY in c(-5, 10, 45, 70, 91.2, 99)) {
    bounds <- feasibility(upperY)
    if (!bounds$feasible) next

    observed <- c(upperY - 30, upperY, upperY + 5)
    start_bottom <- pmin(pmax(min(observed), constants$BottomMin), bounds$bottom_max)
    start_top <- pmin(pmax(max(observed), bounds$top_min), constants$TopMax)

    expect_gte(start_bottom, constants$BottomMin)
    expect_lte(start_bottom, bounds$bottom_max)
    expect_gte(start_top, bounds$top_min)
    expect_lte(start_top, constants$TopMax)
  }
})

test_that("a recovered IC50 outside the acceptance window is rejected by the script's rule", {
  accepted <- function(ic50) ic50 >= constants$IC50Min && ic50 <= constants$IC50Max

  expect_true(accepted(-8))
  expect_true(accepted(constants$IC50Min))
  expect_true(accepted(constants$IC50Max))
  expect_false(accepted(constants$IC50Max + 0.001))
  expect_false(accepted(constants$IC50Min - 0.001))
})

# ------------------------------------------------------------------------------
# Plot helper
# ------------------------------------------------------------------------------

test_that("theme_pub is a usable ggplot theme", {
  skip_if_not_installed("ggplot2")
  suppressMessages(library(ggplot2))

  theme_pub <- extract_object("02_IC50_fitting.R", "theme_pub")

  expect_true(is.function(theme_pub))
  expect_s3_class(theme_pub(), "theme")
})

test_that("a plot built with theme_pub renders without error", {
  skip_if_not_installed("ggplot2")
  suppressMessages(library(ggplot2))

  theme_pub <- extract_object("02_IC50_fitting.R", "theme_pub")
  plot <- ggplot(data.frame(x = 1:5, y = 1:5)) + geom_point(aes(x, y)) + theme_pub()

  expect_silent(invisible(ggplot_build(plot)))
})
