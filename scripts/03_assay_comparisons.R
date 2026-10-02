# ==============================================================================
# ===== (2026) H3N2 K Variant Neutralisation: Assay Comparison =================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script compares neutralisation and IC50 assay results and
# analyses IC50 values across virus and vaccination groups.

# ----- 0.2. Dependencies ------------------------------------------------------

library(tidyverse)
library(here)
library(plotrix)
library(patchwork)
library(viridisLite)
library(emmeans)

# ----- 0.3. Helper Functions --------------------------------------------------

options(scipen = 999)

theme_base <- theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        text = element_text(size = 12),
        axis.ticks = element_line(linewidth = 0.5),
        axis.ticks.length = unit(-3, "pt"),
        axis.minor.ticks.length = unit(-2, "pt"))

# ------------------------------------------------------------------------------
# ----- 1. Comparing Neutralisation % to IC50 ----------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Format Data -------------------------------------------------------

data_neut <- read_csv(here("data", "neutralisation_75.csv"))
data_fits <- read_csv(here("data", "IC50_fits.csv"), col_types = cols(SampleID = col_character()))
data_meta <- read_csv(here("data", "sample_metadata.csv"), col_types = cols(SampleID = col_character()))

data_neut <- data_neut %>% group_by(SampleID, Virus) %>%
  summarise(u_Neut = mean(Neut, na.rm = T))

data_comp <- left_join(data_fits, data_neut)
data_comp <- left_join(data_comp, data_meta)

# ----- 1.2. Mean Slope Sigmoid Plots ------------------------------------------

mean(filter(data_fits, Slope != 0)$Slope)

data_curved <- data_fits %>%
  filter(Slope != 0, !is.na(IC50))

mean_Slope <- mean(data_curved$Slope)

x <- seq(-14, 0, length.out = 500)
IC50s <- seq(-12, -2, 1)

plot_data <- expand_grid(x = x, IC50 = IC50s) %>%
  mutate(Mean = 100 / (1 + 2^(mean_Slope * (IC50 - x))))

plot_data$IC50 <- factor(plot_data$IC50, levels = IC50s, labels = paste0("1:", 2^abs(IC50s)))

shifts <- seq(-3.5, 3.5, 1)

bracket_points <- tibble(x = log2(1), y = 100 / (1 + 2^(mean_Slope * shifts)))

p1 <- ggplot(plot_data, aes(x = x, y = Mean, group = IC50, colour = IC50)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = log2(1/75), linetype = "dotted") +
  geom_point(data = bracket_points, aes(x = x, y = y), inherit.aes = FALSE, size = 2) +
  scale_x_continuous(breaks = -seq(0, 14, 2), labels = paste0("1:", 2^seq(0, 14, 2)), minor_breaks = -seq(1, 13, 2), guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(limits = c(0, 100), minor_breaks = seq(0, 100, 5), guide = guide_axis(minor.ticks = TRUE)) +
  scale_color_manual(values = rev(mako(13)[2:12])) +
  labs(x = "Relative Dilution", y = "Neutralisation (%)", colour = expression(IC[50])) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE, label.position = "top", label.theme = element_text(angle = 45, hjust = 0, vjust = 0.5))) +
  theme_base +
  theme(aspect.ratio = 1,
        legend.position = "top",
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))

# ----- 1.3. Delta to Foldchange Function --------------------------------------

log2_shift <- function(p1, p2) {abs((log2(p1 / (1 - p1)) - log2(p2 / (1 - p2))) / mean_Slope)}

log2_shift(0.742, 0.538)
2^log2_shift(0.742, 0.538)

# ----- 1.4. Correlations ------------------------------------------------------

data_comp$u_Neut <- pmin(pmax(data_comp$u_Neut, 0), 100) / 100

epsilon <- 0.025

data_comp_plot <- data_comp %>%
  filter(!is.na(IC50), !is.na(u_Neut)) %>%
  mutate(u_Neut_adj = case_when(u_Neut <= 0 ~ epsilon, u_Neut >= 1 ~ 1 - epsilon, TRUE ~ u_Neut), NeutLogit = qlogis(u_Neut_adj))

cor.test(data_comp_plot$IC50, data_comp_plot$NeutLogit)

model_logit <- lm(NeutLogit ~ IC50, data = data_comp_plot)
summary(model_logit)

pred_data <- tibble(IC50 = seq(min(data_comp_plot$IC50), max(data_comp_plot$IC50), length.out = 500))

pred_data$NeutLogit <- predict(model_logit, newdata = pred_data)
pred_data$u_Neut <- plogis(pred_data$NeutLogit)

# ----- 1.5. Supplementary Plots -----------------------------------------------

p2 <- ggplot(data_comp_plot) +
  geom_point(aes(x = IC50, y = u_Neut * 100), alpha = 0.5) +
  geom_line(data = pred_data, aes(x = IC50, y = u_Neut * 100), linewidth = 1) +
  scale_x_continuous(breaks = -seq(2, 12, 2), labels = paste0("1:", 2^seq(2, 12, 2)), minor_breaks = -seq(3, 11, 2), guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(limits = c(0, 100), minor_breaks = seq(0, 100, 5), guide = guide_axis(minor.ticks = TRUE)) +
  labs(x = expression(IC[50]), y = "Neutralisation (%)") +
  theme_base +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(data_comp_plot) +
  geom_point(aes(x = IC50, y = NeutLogit), alpha = 0.5) +
  geom_line(data = pred_data, aes(x = IC50, y = NeutLogit), linewidth = 1) +
  scale_x_continuous(breaks = -seq(2, 12, 2), labels = paste0("1:", 2^seq(2, 12, 2)), minor_breaks = -seq(3, 11, 2), guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(breaks = seq(-3, 6, 3), minor_breaks = seq(-5, 10, 1), guide = guide_axis(minor.ticks = TRUE)) +
  labs(x = expression(IC[50]), y = "Logit(neutralisation)") +
  theme_base +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 45, hjust = 1))

plots <- p1 + p2 + p3

plots

# ------------------------------------------------------------------------------
# ----- 2. IC50 Models ---------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. Wrangles and Checks -----------------------------------------------

data_comp$Virus <- factor(data_comp$Virus, levels = c("THA22", "ENG24", "ENG25"))

data_comp$vacc_group <- ifelse(data_comp$vacc_years_since_last == ">4", ">4",
                              ifelse(data_comp$vacc_years_since_last == "1", "1", "2-4"))

data_comp$vacc_group <- factor(data_comp$vacc_group, levels = c(">4", "2-4", "1"))

data_comp <- mutate(data_comp, Age_policy = cut(Age, breaks = c(18, 50, 65, Inf), right = FALSE,
                                    labels = c("18-49", "50-64", "65+")))

data_comp$IC50 <- ifelse(data_comp$Verdict == "Flattened",
                         ifelse(data_comp$Lower < 50, -4, -12), data_comp$IC50)

data_comp <- data_comp %>%
  mutate(VirusVacc = interaction(Virus, vacc_group, sep = " × "))

p4 <- ggplot(data_comp) +
  geom_boxplot(aes(x = VirusVacc, y = -IC50, fill = Virus), alpha = 0.75, outlier.shape = NA) +
  theme_base +
  theme(axis.text.x = element_blank(),
        aspect.ratio = 1/2) +
  # y is -IC50 = log2(dilution), so a break at 4 is a titre of 1:16. The labels
  # must therefore ascend with the breaks (they were previously reversed, which
  # printed 1:4096 against the 1:16 tick and vice versa).
  scale_y_continuous(limits = c(3, 12), breaks = seq(4, 12, 2), minor_breaks = seq(4, 13, 1), guide = guide_axis(minor.ticks = TRUE),
                     labels = paste0("1:", 2^seq(4, 12, 2))) +
  scale_fill_manual(values = c("#4677AA", "#DECC76", "#CC6476"),
                    labels = c("H3N2/2022-Vac",
                               "H3N2/2024-J",
                               "H3N2/2025-K")) +
  labs(x = "Virus × Years Since Last Vaccination", y = "IC50")

p4

# ----- 2.2. Model Fit and Contrasts -------------------------------------------

model_IC50 <- lm(IC50 ~ Virus +  vacc_group + Virus * vacc_group, data = data_comp)

summary(model_IC50)

# Virus effects within vaccine group
emm_virus_vacc <- emmeans(model_IC50, ~ Virus | vacc_group, type = "response", weights = "equal")
emm_virus_vacc

contrast(emm_virus_vacc, method = "pairwise", adjust = "none")

# Vaccine group effects within viruses
emm_vacc_virus <- emmeans(model_IC50, ~ vacc_group | Virus, type = "response", weights = "equal")
emm_vacc_virus

contrast(emm_vacc_virus, method = "pairwise")

