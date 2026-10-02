# ==============================================================================
# ===== (2026) H3N2 K Variant Neutralisation: Dependencies Setup ===============
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script installs all R package dependencies used in this study.

# Note: This script installs the current CRAN release of any package that is not
#       already installed.

# Package versions used in the original analysis:
#
# tidyverse    2.0.0
# here         1.0.2
# patchwork    1.3.2
# betareg      3.2.4
# emmeans      2.0.1
# minpack.lm   1.2.4
# plotrix      3.8.13
# viridisLite  0.4.2
# svglite      2.2.2

# ------------------------------------------------------------------------------
# ----- 1. Install Packages ----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Specify Packages --------------------------------------------------

# Note: svglite is not loaded with library() by any script, but ggsave() requires
#       it to write the .svg figure panels in 01_neutralisation_analysis.R.

packages <- c("tidyverse", "here", "betareg", "emmeans", "patchwork", "minpack.lm", "plotrix", "viridisLite", "svglite")

# ----- 1.2. Install Packages --------------------------------------------------

for (package in packages) {
  
  if (requireNamespace(package, quietly = TRUE)) {
    message(package, " already installed, skipping.")
    next
  }
  
  message("Installing ", package, ".")
  install.packages(package, dependencies = TRUE)
}
