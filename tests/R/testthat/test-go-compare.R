# ==============================================================================
# go_compare() -- the AIC model-selection helper in 01_neutralisation_analysis.R
#
# Every model-building decision in section 4 of that script (and therefore the
# choice of model_forward, which produces the published contrasts) rests on this
# function reading AIC differences the right way round. A sign error here would
# invert every "favoured"/"disfavoured" comment in the script without any other
# visible symptom, so the verdict boundaries are pinned exactly.
# ==============================================================================

go_compare <- extract_object("01_neutralisation_analysis.R", "go_compare")

# A minimal model whose AIC is exactly -2*loglik + 2*df, so the boundaries
# between verdicts can be hit precisely rather than approximately.
fake_model <- function(aic, df = 1) {
  structure(list(), class = "fakemod", loglik = (2 * df - aic) / 2, df = df)
}
registerS3method(
  "logLik", "fakemod",
  function(object, ...) structure(attr(object, "loglik"), df = attr(object, "df"), class = "logLik")
)

test_that("go_compare is a function of two models", {
  expect_true(is.function(go_compare))
  expect_named(formals(go_compare), c("m1", "m2"))
})

test_that("the fake model fixture produces the AIC it was asked for", {
  expect_equal(AIC(fake_model(100)), 100)
  expect_equal(AIC(fake_model(42, df = 3)), 42)
})

test_that("the result carries one row per model and the expected columns", {
  comp <- go_compare(fake_model(100), fake_model(110))

  expect_s3_class(comp, "data.frame")
  expect_equal(nrow(comp), 2)
  expect_true(all(c("df", "AIC", "AIC.diff", "AIC.verdict") %in% names(comp)))
})

test_that("AIC.diff is each model's AIC minus the other's", {
  comp <- go_compare(fake_model(100), fake_model(110))

  expect_equal(comp$AIC.diff[1], -10)  # m1 is 10 AIC better
  expect_equal(comp$AIC.diff[2], 10)
})

test_that("AIC.diff is antisymmetric", {
  comp <- go_compare(fake_model(100), fake_model(137))
  expect_equal(comp$AIC.diff[1], -comp$AIC.diff[2])
})

test_that("identical fits are called equivalent", {
  comp <- go_compare(fake_model(100), fake_model(100))

  expect_equal(comp$AIC.diff, c(0, 0))
  expect_equal(comp$AIC.verdict, c("Equivalent", "Equivalent"))
})

test_that("a clearly better model is strongly favoured and its rival disfavoured", {
  comp <- go_compare(fake_model(100), fake_model(200))

  expect_equal(comp$AIC.verdict[1], "Strongly Favour")
  expect_equal(comp$AIC.verdict[2], "Strongly Disfavour")
})

test_that("the verdict is always stated from the perspective of its own row", {
  comp <- go_compare(fake_model(200), fake_model(100))

  expect_equal(comp$AIC.verdict[1], "Strongly Disfavour")
  expect_equal(comp$AIC.verdict[2], "Strongly Favour")
})

# ------------------------------------------------------------------------------
# Verdict boundaries. delta is m1's AIC minus m2's.
# ------------------------------------------------------------------------------

# The cut points are all strict `>` comparisons, which makes every band
# half-open in the same direction:
#
#        delta <= -6   Strongly Favour
#   -6 < delta <= -2   Weakly Favour
#   -2 < delta <=  2   Equivalent
#    2 < delta <=  6   Weakly Disfavour
#    6 < delta         Strongly Disfavour
boundaries <- list(
  list(delta = -8,   verdict = "Strongly Favour"),
  list(delta = -6.1, verdict = "Strongly Favour"),
  list(delta = -6,   verdict = "Strongly Favour"),
  list(delta = -5.9, verdict = "Weakly Favour"),
  list(delta = -4,   verdict = "Weakly Favour"),
  list(delta = -2.1, verdict = "Weakly Favour"),
  list(delta = -2,   verdict = "Weakly Favour"),
  list(delta = -1.9, verdict = "Equivalent"),
  list(delta = 0,    verdict = "Equivalent"),
  list(delta = 2,    verdict = "Equivalent"),
  list(delta = 2.1,  verdict = "Weakly Disfavour"),
  list(delta = 4,    verdict = "Weakly Disfavour"),
  list(delta = 6,    verdict = "Weakly Disfavour"),
  list(delta = 6.1,  verdict = "Strongly Disfavour"),
  list(delta = 20,   verdict = "Strongly Disfavour")
)

for (case in boundaries) {
  local({
    delta <- case$delta
    expected <- case$verdict
    test_that(sprintf("delta of %+g is '%s'", delta, expected), {
      comp <- go_compare(fake_model(100 + delta), fake_model(100))
      expect_equal(comp$AIC.diff[1], delta)
      expect_equal(comp$AIC.verdict[1], expected)
    })
  })
}

test_that("verdicts are drawn only from the documented vocabulary", {
  vocabulary <- c("Strongly Favour", "Weakly Favour", "Equivalent",
                  "Weakly Disfavour", "Strongly Disfavour")

  for (delta in seq(-20, 20, by = 0.5)) {
    comp <- go_compare(fake_model(100 + delta), fake_model(100))
    expect_true(all(comp$AIC.verdict %in% vocabulary))
  }
})

test_that("the classification is asymmetric exactly on the cut points", {
  # Because every comparison is a strict `>`, a delta of exactly +/-2 or +/-6
  # lands in a different band depending on its sign. This is cosmetic for the
  # script (no comparison in it sits exactly on a cut point) but it is a real
  # asymmetry and worth knowing about before reusing the helper.
  at_minus_two <- go_compare(fake_model(98), fake_model(100))$AIC.verdict[1]
  at_plus_two <- go_compare(fake_model(102), fake_model(100))$AIC.verdict[1]

  expect_equal(at_minus_two, "Weakly Favour")
  expect_equal(at_plus_two, "Equivalent")
  expect_false(at_minus_two == "Equivalent")

  at_minus_six <- go_compare(fake_model(94), fake_model(100))$AIC.verdict[1]
  at_plus_six <- go_compare(fake_model(106), fake_model(100))$AIC.verdict[1]

  expect_equal(at_minus_six, "Strongly Favour")
  expect_equal(at_plus_six, "Weakly Disfavour")
})

test_that("favour and disfavour are exact mirrors away from the cut points", {
  mirror <- c("Strongly Favour" = "Strongly Disfavour",
              "Weakly Favour" = "Weakly Disfavour",
              "Equivalent" = "Equivalent",
              "Weakly Disfavour" = "Weakly Favour",
              "Strongly Disfavour" = "Strongly Favour")

  for (delta in c(-10, -5, -3, -1, 0, 1, 3, 5, 10)) {
    comp <- go_compare(fake_model(100 + delta), fake_model(100))
    expect_equal(unname(mirror[comp$AIC.verdict[1]]), comp$AIC.verdict[2],
                 info = sprintf("delta = %g", delta))
  }
})

test_that("it works on real fitted models", {
  set.seed(42)
  frame <- data.frame(x = rnorm(60), z = rnorm(60))
  frame$y <- 2 * frame$x + rnorm(60, sd = 0.3)

  informative <- lm(y ~ x, data = frame)
  null_model <- lm(y ~ 1, data = frame)

  comp <- go_compare(informative, null_model)

  expect_equal(nrow(comp), 2)
  expect_lt(comp$AIC[1], comp$AIC[2])
  expect_equal(comp$AIC.verdict[1], "Strongly Favour")
})

test_that("adding a useless predictor is not rewarded", {
  set.seed(7)
  frame <- data.frame(x = rnorm(80), noise = rnorm(80))
  frame$y <- 2 * frame$x + rnorm(80, sd = 0.3)

  comp <- go_compare(lm(y ~ x + noise, data = frame), lm(y ~ x, data = frame))

  # One extra parameter that explains nothing costs ~2 AIC.
  expect_gt(comp$AIC.diff[1], 0)
  expect_true(comp$AIC.verdict[1] %in% c("Equivalent", "Weakly Disfavour"))
})
