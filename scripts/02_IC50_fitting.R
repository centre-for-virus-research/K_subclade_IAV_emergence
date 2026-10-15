# ==============================================================================
# ===== (2026) H3N2 K Variant Neutralisation: IC50 Curve Fitting ===============
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script reviews neutralisation dilution-series data and fits
# category-specific four-parameter logistic curves to estimate IC50 values.

# ----- 0.2. Dependencies ------------------------------------------------------

library(tidyverse)
library(here)
library(minpack.lm)

# ----- 0.3. Helper Functions --------------------------------------------------

options(scipen = 999)

theme_pub <- function() {
  theme_classic(ink = "#212121") +
    theme(panel.grid.minor = element_line(color = "#f4f5f8"),
          panel.grid.major = element_line(color = "#eceff4"),
          axis.ticks.length = unit(-3, "pt"),
          axis.minor.ticks.length = unit(-2, "pt"),
          strip.background = element_rect(fill = "#eceff4", linewidth = 0.5))
}

# ------------------------------------------------------------------------------
# ----- 1. Model Fitting -------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Data Wrangling ----------------------------------------------------

data <- read_csv(here("data", "IC50s.csv"), col_types = cols(Sample = col_character()))

data$DilutionLog2 <- log2(1/data$Dilution)

data_bio <- data %>%
  group_by(Batch, Virus, Sample, DilutionLog2) %>%
  summarise(u_Neut = mean(Neut), .groups = "drop")

ggplot(data_bio) +
  geom_point(aes(x = DilutionLog2, y = u_Neut, color = Virus, group = Virus)) +
  scale_y_continuous(limits = c(NA, 100)) +
  facet_wrap(~Sample, axes = "all", axis.labels = "margins") +
  theme_pub() +
  theme(aspect.ratio = 1)

# ----- 1.2. Model Constraints -------------------------------------------------

BottomMin <- -50
BottomMax <- 100

TopMin <- -50
TopMax <- 100

IC50Min <- -20
IC50Max <- -4

SlopeMin <- 0.6
SlopeMax <- 2.5

combinations <- distinct(data_bio, Sample, Virus)

# ----- 1.3. Review and Fit Curves ---------------------------------------------

curve_parameters <- transmute(combinations, SampleID = Sample, Virus, Category = NA, Lower = NA, Upper = NA, IC50 = NA, Slope = NA)

for (i in seq_len(nrow(combinations))) {
  
  sample_i <- combinations$Sample[i]
  virus_i <- combinations$Virus[i]
  
  data_i <- filter(data_bio, Sample == sample_i, Virus == virus_i) %>%
    arrange(DilutionLog2)
  
  minY <- min(data_i$u_Neut, na.rm = TRUE)
  maxY <- max(data_i$u_Neut, na.rm = TRUE)
  
  minY_bound <- max(minY, BottomMin)
  maxY_bound <- min(maxY, TopMax)
  
  xmin <- min(data_i$DilutionLog2, na.rm = TRUE)
  xmax <- max(data_i$DilutionLog2, na.rm = TRUE)
  
  plot_ymin <- min(-50, floor(minY / 10) * 10 - 5)
  
  p <- ggplot(data_i) +
    geom_point(aes(x = DilutionLog2, y = u_Neut)) +
    scale_y_continuous(limits = c(plot_ymin, 100)) +
    labs(title = paste0(sample_i, " ", virus_i, " — ", i, " / ", nrow(combinations))) +
    theme_pub() +
    theme(aspect.ratio = 1)
  
  print(p)
  
  repeat {
    
    choice <- readline("1 = fully observed, 2 = lower observed, 3 = upper observed, 4 = flattened low, 5 = flattened high: ")
    
    if (!choice %in% c("1", "2", "3", "4", "5")) {
      message("Please enter 1, 2, 3, 4 or 5.")
      next
    }
    
    fit_i <- NULL
    
    # Fully observed
    if (choice == "1") {
      
      BottomLower_i <- max(BottomMin, minY_bound - 5)
      BottomUpper_i <- min(BottomMax, minY_bound + 5)
      
      TopLower_i <- max(TopMin, maxY_bound - 5)
      TopUpper_i <- min(TopMax, maxY_bound + 5)
      
      midpoint_i <- ((BottomLower_i + BottomUpper_i) / 2 +
                       (TopLower_i + TopUpper_i) / 2) / 2
      
      IC50Start_i <- data_i$DilutionLog2[
        which.min(abs(data_i$u_Neut - midpoint_i))
      ]
      
      IC50Start_i <- pmin(pmax(IC50Start_i, xmin), xmax)
      
      fit_i <- tryCatch(nlsLM(u_Neut ~ Bottom + (Top - Bottom) * plogis((DilutionLog2 - IC50) * Slope), data = data_i, start = list(Bottom = mean(c(BottomLower_i, BottomUpper_i)), Top = mean(c(TopLower_i, TopUpper_i)), IC50 = IC50Start_i, Slope = 1), lower = c(Bottom = BottomLower_i, Top = TopLower_i, IC50 = xmin, Slope = SlopeMin), upper = c(Bottom = BottomUpper_i, Top = TopUpper_i, IC50 = xmax, Slope = SlopeMax), control = nls.lm.control(maxiter = 1000)), error = function(e) NULL)
      
      if (is.null(fit_i)) {
        message("Fit failed. Choose another category.")
        next
      }
    }
    
    # Lower observed
    if (choice == "2") {
      
      BottomLower_i <- max(BottomMin, minY_bound - 5)
      BottomUpper_i <- min(BottomMax, minY_bound + 5)
      
      TopLower_i <- max(TopMin, maxY_bound - 5)
      TopUpper_i <- TopMax
      
      midpoint_i <- ((BottomLower_i + BottomUpper_i) / 2 +
                       (TopLower_i + TopUpper_i) / 2) / 2
      
      IC50Start_i <- data_i$DilutionLog2[
        which.min(abs(data_i$u_Neut - midpoint_i))
      ]
      
      IC50Start_i <- pmin(pmax(IC50Start_i, IC50Min), IC50Max)
      
      fit_i <- tryCatch(nlsLM(u_Neut ~ Bottom + (Top - Bottom) * plogis((DilutionLog2 - IC50) * Slope), data = data_i, start = list(Bottom = mean(c(BottomLower_i, BottomUpper_i)), Top = mean(c(TopLower_i, TopUpper_i)), IC50 = IC50Start_i, Slope = 1), lower = c(Bottom = BottomLower_i, Top = TopLower_i, IC50 = IC50Min, Slope = SlopeMin), upper = c(Bottom = BottomUpper_i, Top = TopUpper_i, IC50 = IC50Max, Slope = SlopeMax), control = nls.lm.control(maxiter = 1000)), error = function(e) NULL)
      
      if (is.null(fit_i)) {
        message("Fit failed. Choose another category.")
        next
      }
    }
    
    # Upper observed
    if (choice == "3") {
      
      BottomLower_i <- BottomMin
      BottomUpper_i <- min(BottomMax, minY_bound + 5)
      
      TopLower_i <- max(TopMin, maxY_bound - 5)
      TopUpper_i <- min(TopMax, maxY_bound + 5)
      
      midpoint_i <- ((BottomLower_i + BottomUpper_i) / 2 +
                       (TopLower_i + TopUpper_i) / 2) / 2
      
      IC50Start_i <- data_i$DilutionLog2[
        which.min(abs(data_i$u_Neut - midpoint_i))
      ]
      
      IC50Start_i <- pmin(pmax(IC50Start_i, IC50Min), IC50Max)
      
      fit_i <- tryCatch(nlsLM(u_Neut ~ Bottom + (Top - Bottom) * plogis((DilutionLog2 - IC50) * Slope), data = data_i, start = list(Bottom = mean(c(BottomLower_i, BottomUpper_i)), Top = mean(c(TopLower_i, TopUpper_i)), IC50 = IC50Start_i, Slope = 1), lower = c(Bottom = BottomLower_i, Top = TopLower_i, IC50 = IC50Min, Slope = SlopeMin), upper = c(Bottom = BottomUpper_i, Top = TopUpper_i, IC50 = IC50Max, Slope = SlopeMax), control = nls.lm.control(maxiter = 1000)), error = function(e) NULL)
      
      if (is.null(fit_i)) {
        message("Fit failed. Choose another category.")
        next
      }
    }
    
    break
  }
  
  if (choice == "1") {
    
    coefs <- coef(fit_i)
    
    curve_parameters$Category[i] <- "Fully observed"
    curve_parameters$Lower[i] <- unname(coefs["Bottom"])
    curve_parameters$Upper[i] <- unname(coefs["Top"])
    curve_parameters$IC50[i] <- unname(coefs["IC50"])
    curve_parameters$Slope[i] <- unname(coefs["Slope"])
  }
  
  if (choice == "2") {
    
    coefs <- coef(fit_i)
    
    curve_parameters$Category[i] <- "Lower observed"
    curve_parameters$Lower[i] <- unname(coefs["Bottom"])
    curve_parameters$Upper[i] <- unname(coefs["Top"])
    curve_parameters$IC50[i] <- unname(coefs["IC50"])
    curve_parameters$Slope[i] <- unname(coefs["Slope"])
  }
  
  if (choice == "3") {
    
    coefs <- coef(fit_i)
    
    curve_parameters$Category[i] <- "Upper observed"
    curve_parameters$Lower[i] <- unname(coefs["Bottom"])
    curve_parameters$Upper[i] <- unname(coefs["Top"])
    curve_parameters$IC50[i] <- unname(coefs["IC50"])
    curve_parameters$Slope[i] <- unname(coefs["Slope"])
  }
  
  if (choice == "4") {
    
    curve_parameters$Category[i] <- "Flattened low"
    curve_parameters$Lower[i] <- 0
    curve_parameters$Upper[i] <- 0
    curve_parameters$IC50[i] <- -4
    curve_parameters$Slope[i] <- 0
  }
  
  if (choice == "5") {
    
    curve_parameters$Category[i] <- "Flattened high"
    curve_parameters$Lower[i] <- 100
    curve_parameters$Upper[i] <- 100
    curve_parameters$IC50[i] <- -20
    curve_parameters$Slope[i] <- 0
  }
}

# ----- 1.4. Save Fits ----------------------------------------------------------

curve_parameters <- curve_parameters %>%
  mutate(across(c(Lower, Upper, IC50, Slope), ~ round(.x, 4)))

write_csv(curve_parameters, here("data", "IC50_fits.csv"))