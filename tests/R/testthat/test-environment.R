# ==============================================================================
# The R half of the environment
#
# scripts/00_setup.R installs whatever CRAN currently ships, so the versions in
# a fresh environment will not match the ones recorded in the README. These
# tests check what actually matters: that every package the scripts load is
# present, and that it is new enough for the features they use.
# ==============================================================================

required_packages <- c("tidyverse", "here", "betareg", "emmeans",
                       "patchwork", "minpack.lm", "plotrix", "viridisLite")

test_that("R is at least 4.5", {
  expect_gte(getRversion(), package_version("4.5.0"))
})

for (pkg in required_packages) {
  local({
    package <- pkg
    test_that(sprintf("%s is installed", package), {
      expect_true(requireNamespace(package, quietly = TRUE))
    })
  })
}

test_that("00_setup.R installs every package the analysis scripts load", {
  setup_lines <- script_source("00_setup.R")
  declared_line <- grep("^packages\\s*<-", setup_lines, value = TRUE)
  expect_length(declared_line, 1)

  declared <- regmatches(declared_line, gregexpr('"[^"]+"', declared_line))[[1]]
  declared <- gsub('"', "", declared)

  loaded <- character(0)
  for (script in c("01_neutralisation_analysis.R", "02_IC50_fitting.R", "03_assay_comparisons.R")) {
    lines <- script_source(script)
    matches <- regmatches(lines, regexpr("^\\s*library\\(([A-Za-z0-9._]+)\\)", lines))
    matches <- matches[nzchar(matches)]
    loaded <- c(loaded, gsub("^\\s*library\\(|\\)$", "", matches))
  }

  expect_true(all(unique(loaded) %in% declared),
              info = paste("not installed by 00_setup.R:",
                           paste(setdiff(unique(loaded), declared), collapse = ", ")))
})

test_that("the README documents every package 00_setup.R installs", {
  readme <- readLines(file.path(project_root(), "README.md"), warn = FALSE)
  setup_lines <- script_source("00_setup.R")
  declared_line <- grep("^packages\\s*<-", setup_lines, value = TRUE)
  declared <- gsub('"', "", regmatches(declared_line, gregexpr('"[^"]+"', declared_line))[[1]])

  undocumented <- declared[!vapply(declared, function(p) any(grepl(p, readme, fixed = TRUE)), logical(1))]
  expect_equal(length(undocumented), 0,
               info = paste("missing from the README:", paste(undocumented, collapse = ", ")))
})

test_that("ggplot2 is new enough for the theme and guide features the scripts use", {
  skip_if_not_installed("ggplot2")

  # theme_classic(ink = ) arrived in ggplot2 4.0.0; axis.minor.ticks.length and
  # guide_axis(minor.ticks = ) arrived in 3.5.0.
  expect_gte(packageVersion("ggplot2"), package_version("4.0.0"))
})

test_that("theme_classic accepts the ink argument script 02 passes", {
  skip_if_not_installed("ggplot2")
  suppressMessages(library(ggplot2))

  expect_s3_class(theme_classic(ink = "#212121"), "theme")
})

test_that("guide_axis accepts minor.ticks, as every scale in the scripts assumes", {
  skip_if_not_installed("ggplot2")
  suppressMessages(library(ggplot2))

  expect_no_error(guide_axis(minor.ticks = TRUE))
})

test_that("every file format the scripts ggsave() to can actually be written", {
  # This is the gap that `library()` scanning cannot see: ggsave() dispatches to
  # a device package by file extension, so .svg needs svglite even though no
  # script ever calls library(svglite). Without it, 01 ran every model and then
  # died at the first ggsave().
  skip_if_not_installed("ggplot2")
  suppressMessages(library(ggplot2))

  extensions <- character(0)
  for (script in c("01_neutralisation_analysis.R", "02_IC50_fitting.R", "03_assay_comparisons.R")) {
    lines <- script_source(script)
    hits <- regmatches(lines, gregexpr('"[^"]*\\.(svg|png|pdf|jpe?g|tiff|eps)"', lines))
    extensions <- c(extensions, tools::file_ext(gsub('"', "", unlist(hits))))
  }
  extensions <- unique(extensions)
  skip_if(length(extensions) == 0, "no figures are saved by the scripts")

  plot <- ggplot(data.frame(x = 1:3, y = 1:3)) + geom_point(aes(x, y))
  for (ext in extensions) {
    target <- file.path(tempdir(), paste0("device-check.", ext))
    expect_no_error(suppressMessages(ggsave(target, plot = plot, width = 3, height = 2)))
    expect_true(file.exists(target), info = ext)
    unlink(target)
  }
})

test_that("svglite is installed, since the figure panels are saved as SVG", {
  expect_true(requireNamespace("svglite", quietly = TRUE))
})

test_that("00_setup.R installs svglite even though no script library()s it", {
  setup_lines <- script_source("00_setup.R")
  declared_line <- grep("^packages\\s*<-", setup_lines, value = TRUE)

  expect_true(grepl("svglite", declared_line, fixed = TRUE))
})

test_that("here() resolves to the repository root", {
  skip_if_not_installed("here")

  expect_true(dir.exists(file.path(project_root(), "data")))
  expect_true(file.exists(file.path(project_root(), "README.md")))
})

test_that("every data file the scripts read is present", {
  for (file in c("neutralisation_75.csv", "sample_metadata.csv", "IC50s.csv", "IC50_fits.csv")) {
    expect_true(file.exists(data_path(file)), info = file)
  }
})

test_that("emmeans can talk to betareg", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("betareg")
  suppressMessages({
    library(emmeans)
    library(betareg)
  })

  set.seed(1)
  frame <- data.frame(y = runif(60, 0.1, 0.9), g = factor(sample(c("a", "b"), 60, TRUE)))
  model <- suppressWarnings(betareg(y ~ g, data = frame))

  # emmGrid is S4.
  emm <- suppressMessages(emmeans(model, ~ g, type = "response"))
  expect_s4_class(emm, "emmGrid")
  expect_equal(nrow(as.data.frame(emm)), 2)
})
