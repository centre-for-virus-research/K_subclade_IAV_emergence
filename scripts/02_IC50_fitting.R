# ==============================================================================
# ===== (2026) H3N2 K Variant Neutralisation: IC50 Curve Fitting ===============
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script fits and reviews four-parameter logistic curves to
# neutralisation dilution-series data to estimate IC50 values.

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

# ----- 1.2. Attempt to fit logistic curves ------------------------------------

BottomMin <- -50
BottomMax <- 40

TopMin <- 60
TopMax <- 100

IC50Min <- -20
IC50Max <- -4

SlopeMin <- 0.6
SlopeMax <- 2

combinations <- distinct(data_bio, Sample, Virus)

fits <- list()

fit_summary <- mutate(combinations,
                      Converged = FALSE,
                      Bottom = NA,
                      Top = NA,
                      IC50 = NA,
                      Slope = NA)

for (i in seq_len(nrow(combinations))) {
  
  sample_i <- combinations$Sample[i]
  virus_i <- combinations$Virus[i]
  
  data_i <- filter(data_bio, Sample == sample_i, Virus == virus_i) %>%
    arrange(DilutionLog2)
  
  upper_i <- which.max(data_i$DilutionLog2)
  
  upperX <- data_i$DilutionLog2[upper_i]
  upperY <- data_i$u_Neut[upper_i]
  
  epsilon <- 0.001
  
  BottomMax_i <- min(BottomMax, upperY - epsilon)
  TopMin_i <- max(TopMin, upperY + epsilon)
  
  if (BottomMax_i > BottomMin && TopMin_i < TopMax) {
    
    fit_i <- tryCatch(nlsLM(u_Neut ~ Bottom + (Top - Bottom) * plogis((DilutionLog2 - (upperX - qlogis((upperY - Bottom) / (Top - Bottom)) / Slope)) * Slope), data = data_i, start = list(Bottom = pmin(pmax(min(data_i$u_Neut), BottomMin), BottomMax_i), Top = pmin(pmax(max(data_i$u_Neut), TopMin_i), TopMax), Slope = 1), lower = c(Bottom = BottomMin, Top = TopMin_i, Slope = SlopeMin), upper = c(Bottom = BottomMax_i, Top = TopMax, Slope = SlopeMax), control = nls.lm.control(maxiter = 1000)), error = function(e) NULL)
    
  } else {
    
    fit_i <- NULL
  }
  
  fit_name <- paste(sample_i, virus_i, sep = "_")
  fits[[fit_name]] <- fit_i
  
  if (!is.null(fit_i)) {
    
    coefs <- coef(fit_i)
    
    Bottom_i <- unname(coefs["Bottom"])
    Top_i <- unname(coefs["Top"])
    Slope_i <- unname(coefs["Slope"])
    
    IC50_i <- upperX - qlogis((upperY - Bottom_i) / (Top_i - Bottom_i)) / Slope_i
    
    if (IC50_i >= IC50Min && IC50_i <= IC50Max) {
      
      fit_summary$Converged[i] <- TRUE
      
      fit_summary[i, c("Bottom", "Top", "IC50", "Slope")] <- list(Bottom_i, Top_i, IC50_i, Slope_i)
    }
  }
}


# ----- 1.3. Review Curves ------------------------------------------------------

curve_parameters <- transmute(combinations,
                              SampleID = Sample,
                              Virus,
                              Verdict = NA,
                              Lower = NA,
                              Upper = NA,
                              IC50 = NA,
                              Slope = NA)

for (i in seq_len(nrow(combinations))) {
  
  sample_i <- combinations$Sample[i]
  virus_i <- combinations$Virus[i]
  
  data_i <- filter(data_bio, Sample == sample_i, Virus == virus_i) %>%
    arrange(DilutionLog2)
  
  if (fit_summary$Converged[i]) {
    
    curve_i <- tibble(DilutionLog2 = seq(min(data_i$DilutionLog2), max(data_i$DilutionLog2), length.out = 200)) %>%
      mutate(u_Neut = fit_summary$Bottom[i] +
               (fit_summary$Top[i] - fit_summary$Bottom[i]) /
               (1 + exp(-(DilutionLog2 - fit_summary$IC50[i]) * fit_summary$Slope[i])))
    
  } else {
    curve_i <- NULL
  }
  
  p <- ggplot(data_i) +
    geom_point(aes(x = DilutionLog2, y = u_Neut)) +
    scale_y_continuous(limits = c(-50, 100)) +
    labs(title = paste0(sample_i, " ", virus_i, " — ", i, " / ", nrow(combinations)), subtitle = ifelse(fit_summary$Converged[i], "Logistic fit converged", "Logistic fit did not converge")) +
    theme_pub() +
    theme(aspect.ratio = 1)
  
  if (!is.null(curve_i)) {
    p <- p +
      geom_line(data = curve_i, aes(x = DilutionLog2, y = u_Neut), linewidth = 0.8)
  }
  
  print(p)
  
  repeat {
    
    choice <- readline("1 = curved, 2 = flatten low, 3 = flatten high: ")
    
    if (!choice %in% c("1", "2", "3")) {
      message("Please enter 1, 2 or 3.")
      next
    }
    
    if (choice == "1" && !fit_summary$Converged[i]) {
      message("Logistic fit did not converge. Choose flatten low or flatten high.")
      next
    }
    
    break
  }
  
  if (choice == "1") {
    curve_parameters$Verdict[i] <- "Curved"
    curve_parameters$Lower[i] <- fit_summary$Bottom[i]
    curve_parameters$Upper[i] <- fit_summary$Top[i]
    curve_parameters$IC50[i] <- fit_summary$IC50[i]
    curve_parameters$Slope[i] <- fit_summary$Slope[i]
  }
  
  if (choice == "2") {
    intercept_i <- unname(coef(lm(u_Neut ~ 1, data = data_i))[1])
    
    curve_parameters$Verdict[i] <- "Flattened"
    curve_parameters$Lower[i] <- intercept_i
    curve_parameters$Upper[i] <- intercept_i
    curve_parameters$IC50[i] <- NA_real_
    curve_parameters$Slope[i] <- 0
  }
  
  if (choice == "3") {
    curve_parameters$Verdict[i] <- "Flattened"
    curve_parameters$Lower[i] <- 100
    curve_parameters$Upper[i] <- 100
    curve_parameters$IC50[i] <- NA_real_
    curve_parameters$Slope[i] <- 0
  }
}
