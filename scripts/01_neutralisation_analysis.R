# ==============================================================================
# ===== (2026) H3N2 K Variant Neutralisation: Beta Regression ==================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script plots and analyses neutralisation data collected at a
# dilution of 1:75, using a beta regression model and a forward-model-building
# approach.

# ----- 0.2. Dependencies ------------------------------------------------------

library(tidyverse)
library(patchwork)
library(betareg)
library(here)
library(emmeans)
library(viridisLite)

# ----- 0.3. Load Data ---------------------------------------------------------

data_neutralisation <- read_csv(here("data", "neutralisation_75.csv"))
data_samples <- read_csv(here("data", "sample_metadata.csv"), col_types = cols(SampleID = col_character()))

# ----- 0.4. Convenience Functions ---------------------------------------------

theme_base <- theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        text = element_text(size = 12),
        axis.ticks = element_line(linewidth = 0.5),
        axis.ticks.length = unit(-3, "pt"),
        axis.minor.ticks.length = unit(-2, "pt"))

go_compare <- function(m1, m2) {
  
  comp <- AIC(m1, m2)
  comp$AIC.diff <- c(diff(rev(comp$AIC)), diff(comp$AIC))
  comp$AIC.verdict <- ifelse(comp$AIC.diff > 6, "Strongly Disfavour",
                             ifelse(comp$AIC.diff > 2, "Weakly Disfavour",
                                    ifelse(comp$AIC.diff > -2, "Equivalent",
                                           ifelse(comp$AIC.diff > -6, "Weakly Favour", "Strongly Favour"))))
  
  return(comp)
}

# ------------------------------------------------------------------------------
# ----- 1. Technical Replicates & Noise Checks ---------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Wrangle Technical Replicates --------------------------------------

# Remove sheep antisera
data_neutralisation <- filter(data_neutralisation, SampleID != "837SAnti-H3 1:250")

# Number of technical replicates per Virus * SampleID
table(table(data_neutralisation$Virus, data_neutralisation$SampleID))

data_biorep <- data_neutralisation %>% group_by(SampleID, Virus, Batch) %>%
  summarise(u_Neut = mean(Neut),
            sd_Neut = sd(Neut))

data_biorep$CV_Neut <- data_biorep$sd_Neut / (data_biorep$u_Neut + abs(min(data_biorep$u_Neut)) + 1) 
data_biorep$Continuous <- data_biorep$u_Neut > 20 & data_biorep$u_Neut < 80

# ----- 1.2. Remove repeat patient samples -------------------------------------

data_samples <- filter(data_samples, SampleID %in% data_neutralisation$SampleID)

# Keep earliest sample date for repeat patients
data_samples <- data_samples %>%
  arrange(PatientID, SampleDate) %>%
  distinct(PatientID, .keep_all = TRUE)

# ----- 1.2. Heteroscedasticity ------------------------------------------------

ggplot(data_biorep) +
  geom_density(aes(x = sd_Neut)) +
  theme_base +
  ggtitle("Technical Error (sd) distribution")

ggplot(filter(data_biorep, Continuous)) +
  geom_density(aes(x = sd_Neut)) +
  theme_base +
  ggtitle("Technical Error (sd) distribution, continuous data")

ggplot(data_biorep, aes(x = u_Neut, y = sd_Neut)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_base +
  ggtitle("Heteroscedasticity (sd x mean)")

ggplot(filter(data_biorep, Continuous), aes(x = u_Neut, y = sd_Neut)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_base +
  ggtitle("Heteroscedasticity (sd x mean), continuous data")

ggplot(data_biorep, aes(x = u_Neut, y = CV_Neut)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_base +
  ggtitle("Heteroscedasticity (CV x mean)")

ggplot(filter(data_biorep, Continuous), aes(x = u_Neut, y = CV_Neut)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_base +
  ggtitle("Heteroscedasticity (CV x mean), continuous data")

# ----- 1.3. Batch Effects -----------------------------------------------------

ggplot(data_biorep, aes(x = u_Neut, y = sd_Neut)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~Batch) +
  theme_base +
  ggtitle("Heteroscedasticity by batch (sd x mean)")

ggplot(filter(data_biorep, Continuous), aes(x = u_Neut, y = sd_Neut)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~Batch) +
  theme_base +
  ggtitle("Heteroscedasticity by batch (sd x mean), continuous data")

ggplot(data_biorep, aes(x = u_Neut, y = CV_Neut)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~Batch) +
  theme_base +
  ggtitle("Heteroscedasticity by batch (CV x mean)")

ggplot(filter(data_biorep, Continuous), aes(x = u_Neut, y = sd_Neut)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~Batch) +
  theme_base +
  ggtitle("Heteroscedasticity by batch (CV x mean), continuous data")

# ------------------------------------------------------------------------------
# ----- 2. Biological Replicate Checks -----------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. Batch Effects on Means --------------------------------------------

ggplot(data_biorep, aes(x = Virus, y = u_Neut)) +
  geom_boxplot() +
  facet_wrap(~Batch) +
  theme_base

# ----- 2.2. Batch * Virus Effects (Mean & Dispersion) -------------------------

data_biorep$Batch <- factor(data_biorep$Batch)
data_biorep$u_Neut <- pmax(pmin(1, data_biorep$u_Neut/100), 0)

# Smithson & Verkuilen transformation for == 0 or 1
n <- nrow(data_biorep)
data_biorep$u_Neut <- (data_biorep$u_Neut * (n - 1) + 0.5) / n

# Maximal model, Batch * Virus effects on means and dispersions
model__fixed_full__disp_full <- betareg(u_Neut ~ Batch * Virus | Batch * Virus, data = data_biorep)

# Remove dispersion interaction
model__fixed_full__disp_noInt <- betareg(u_Neut ~ Batch * Virus | Batch + Virus, data = data_biorep)
go_compare(model__fixed_full__disp_noInt, model__fixed_full__disp_full)

# Remove dispersion Virus
model__fixed_full__disp_noVirus <- betareg(u_Neut ~ Batch * Virus | Batch, data = data_biorep)
go_compare(model__fixed_full__disp_noVirus, model__fixed_full__disp_noInt)

# Remove dispersion Batch
model__fixed_full__disp_noBatch <- betareg(u_Neut ~ Batch * Virus | Virus, data = data_biorep)
go_compare(model__fixed_full__disp_noBatch, model__fixed_full__disp_noInt)

# Remove fixed effect interaction
model__fixed_noInt__disp_noInt <- betareg(u_Neut ~ Batch + Virus | Batch + Virus, data = data_biorep)
go_compare(model__fixed_noInt__disp_noInt, model__fixed_full__disp_noInt)

# Remove fixed effect Virus
model__fixed_noVirus__disp_noInt <- betareg(u_Neut ~ Batch | Batch + Virus, data = data_biorep)
go_compare(model__fixed_noVirus__disp_noInt, model__fixed_noInt__disp_noInt)

# Remove fixed effect Batch
model__fixed_noBatch__disp_noInt <- betareg(u_Neut ~ Virus | Batch + Virus, data = data_biorep)
go_compare(model__fixed_noBatch__disp_noInt, model__fixed_noInt__disp_noInt)

# ----- 2.3. Minimal Adequate Model --------------------------------------------

summary(model__fixed_noInt__disp_noInt)
pairs(emmeans(model__fixed_noInt__disp_noInt, ~ Virus, type = "response"), adjust = "tukey")
pairs(emmeans(model__fixed_noInt__disp_noInt, ~ Batch, type = "response"), adjust = "tukey")

res <- residuals(model__fixed_noInt__disp_noInt, type = "quantile")
fit <- fitted(model__fixed_noInt__disp_noInt)

ggplot(data.frame(fit, res), aes(x = fit, y = res)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_base

# ------------------------------------------------------------------------------
# ----- 3. Representing Vaccination and Age ------------------------------------
# ------------------------------------------------------------------------------

# ----- 3.1. Join Sample Metadata ----------------------------------------------

data_all <- left_join(data_biorep, data_samples)

data_all <- filter(data_all, !is.na(PatientID)) # Removes repeat patient samples

data_all$Virus <- factor(data_all$Virus, levels = c("THA22", "ENG24", "ENG25"))

data_all$vacc_number_of_doses <- factor(data_all$vacc_number_of_doses)

data_all$vacc_years_since_last <- factor(data_all$vacc_years_since_last, levels = c(">4", "4", "3", "2", "1"))

# ----- 3.2. Representing Vaccination ------------------------------------------

model_doses <- betareg(u_Neut ~ Batch + Virus + Age + Sex + vacc_number_of_doses |
                         Batch + Virus,
                       data = data_all)

model_recency <- betareg(u_Neut ~ Batch + Virus + Age + Sex + vacc_years_since_last |
                           Batch + Virus,
                         data = data_all)

model_none <- betareg(u_Neut ~ Batch + Virus + Age + Sex |
                        Batch + Virus,
                      data = data_all)

go_compare(model_doses, model_none)
go_compare(model_recency, model_none)
go_compare(model_doses, model_recency)


data_all$vacc_group <- ifelse(data_all$vacc_years_since_last == ">4", ">4",
                              ifelse(data_all$vacc_years_since_last == "1", "1", "2-4"))

data_all$vacc_group <- factor(data_all$vacc_group, levels = c(">4", "2-4", "1"))

model_coarse <- betareg(u_Neut ~ Batch + Virus + Age + Sex + vacc_group |
                          Batch + Virus,
                        data = data_all)

go_compare(model_coarse, model_recency) # Some loss with coarseness

table(data_all$Virus, data_all$Sex, data_all$vacc_years_since_last)
table(data_all$Virus, data_all$Sex, data_all$vacc_group) # Coarseness may still be justified in complex models due to subgroup ns

# ----- 3.3. Representing Age --------------------------------------------------

data_all <- mutate(data_all,
                   Age_10y = cut(Age, breaks = c(18, 30, 40, 50, 60, 70, 80, 90, 100), right = FALSE,
                                 labels = c("18-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+")),
                   Age_20y = cut(Age, breaks = c(18, 40, 60, 80, Inf), right = FALSE,
                                 labels = c("18-39", "40-59", "60-79", "80+")),
                   Age_policy = cut(Age, breaks = c(18, 50, 65, Inf), right = FALSE,
                                      labels = c("18-49", "50-64", "65+")),
                   Age_3groups = cut(Age, breaks = c(18, 40, 65, Inf), right = FALSE,
                                          labels = c("18-39", "40-64", "65+")))

model_none <- betareg(u_Neut ~ Batch + Virus + Sex + vacc_group |
                       Batch + Virus,
                     data = data_all)

model_age <- betareg(u_Neut ~ Batch + Virus + Age + Sex + vacc_group |
                       Batch + Virus,
                     data = data_all)

model_age10 <- betareg(u_Neut ~ Batch + Virus + Age_10y + Sex + vacc_group |
                       Batch + Virus,
                     data = data_all)

model_age20 <- betareg(u_Neut ~ Batch + Virus + Age_20y + Sex + vacc_group |
                       Batch + Virus,
                     data = data_all)

model_agepolicy <- betareg(u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group |
                         Batch + Virus,
                       data = data_all)

model_age3groups <- betareg(u_Neut ~ Batch + Virus + Age_3groups + Sex + vacc_group |
                         Batch + Virus,
                       data = data_all)

go_compare(model_none, model_age)
go_compare(model_none, model_age10)
go_compare(model_none, model_age20)
go_compare(model_none, model_agepolicy)
go_compare(model_none, model_age3groups)

go_compare(model_age, model_age10)
go_compare(model_age, model_age20)
go_compare(model_age, model_agepolicy)
go_compare(model_age, model_age3groups)

go_compare(model_age10, model_age20)
go_compare(model_age10, model_agepolicy)
go_compare(model_age20, model_agepolicy)

# ------------------------------------------------------------------------------
# ----- 4. Beta Regression -----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 4.1. Forward Model Building (first-order effects) ----------------------

model_1 <- betareg(u_Neut ~ Batch + Virus |
                     Batch + Virus,
                   data = data_all)

model_1a <- betareg(u_Neut ~ Batch + Virus + Age_policy |
                      Batch + Virus,
                    data = data_all)

model_1b <- betareg(u_Neut ~ Batch + Virus + Sex |
                      Batch + Virus,
                    data = data_all)

model_1c <- betareg(u_Neut ~ Batch + Virus + vacc_group |
                      Batch + Virus,
                    data = data_all)

go_compare(model_1, model_1a) # 1a strongly favoured
go_compare(model_1, model_1b) # 1b strongly favoured
go_compare(model_1, model_1c) # 1c strongly favoured

go_compare(model_1a, model_1b) # 1a strongly favoured
go_compare(model_1a, model_1c) # 1c strongly favoured
go_compare(model_1b, model_1c) # 1c strongly favoured

model_1ac <- betareg(u_Neut ~ Batch + Virus + Age_policy + vacc_group |
                       Batch + Virus,
                     data = data_all)

go_compare(model_1a, model_1ac) # 1ac strongly favoured
go_compare(model_1c, model_1ac) # 1ac strongly favoured

model_1abc <- betareg(u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group |
                        Batch + Virus,
                      data = data_all)

go_compare(model_1ac, model_1abc) # 1abc weakly favoured

# ----- 4.2. Forward Model Building (second-order effects) ---------------------

model_2 <- model_1abc

model_2a <- betareg(u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group + 
                      Virus * vacc_group |
                      Batch + Virus,
                    data = data_all)

model_2b <- betareg(u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group + 
                      Virus * Sex |
                      Batch + Virus,
                    data = data_all)

model_2c <- betareg(u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group + 
                      Virus * Age_policy |
                      Batch + Virus,
                    data = data_all)

model_2d <- betareg(u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group + 
                      Age_policy * Sex |
                      Batch + Virus,
                    data = data_all)

model_2e <- betareg(u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group + 
                      Age_policy * vacc_group |
                      Batch + Virus,
                    data = data_all)

model_2f <- betareg(u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group + 
                      Sex * vacc_group |
                      Batch + Virus,
                    data = data_all)

go_compare(model_2, model_2a) # weakly favoured
go_compare(model_2, model_2b) # weakly disfavoured
go_compare(model_2, model_2c) # equivalent
go_compare(model_2, model_2d) # equivalent
go_compare(model_2, model_2e) # weakly favoured
go_compare(model_2, model_2f) # weakly disfavoured

go_compare(model_2a, model_2c) # 2a favoured over 2c
go_compare(model_2a, model_2d) # 2a favoured over 2d
go_compare(model_2a, model_2e) # 2a and 2e equivalent
go_compare(model_2c, model_2d) # 2c favoured over 2d
go_compare(model_2c, model_2e) # 2e favoured over 2c
go_compare(model_2d, model_2e) # 2e favoured over 2d

model_2ae <- betareg(
  u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group +
    Virus * vacc_group +
    Age_policy * vacc_group |
    Batch + Virus,
  data = data_all
)

go_compare(model_2a, model_2ae) # 2ae weakly favoured
go_compare(model_2e, model_2ae) # 2ae weakly favoured

model_2ace <- betareg(
  u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group +
    Virus * vacc_group +
    Virus * Age_policy +
    Age_policy * vacc_group |
    Batch + Virus,
  data = data_all
)

go_compare(model_2ae, model_2ace) # 2ae weakly favoured

model_2ade <- betareg(
  u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group +
    Virus * vacc_group +
    Age_policy * Sex +
    Age_policy * vacc_group |
    Batch + Virus,
  data = data_all
)

go_compare(model_2ae, model_2ade) # Equivalent

model_2abe <- betareg(
  u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group +
    Virus * vacc_group +
    Virus * Sex +
    Age_policy * vacc_group |
    Batch + Virus,
  data = data_all
)

go_compare(model_2ae, model_2abe) # Equivalent

# ----- 4.3. Forward Model Building (third-order effects) ----------------------

model_3 <- betareg(
  u_Neut ~ Batch + Virus + Age_policy + Sex + vacc_group +
    Virus * vacc_group +
    Age_policy * vacc_group +
    Virus * vacc_group * Age_policy |
    Batch + Virus,
  data = data_all
)

go_compare(model_2ae, model_3) # 2ae strongly favoured

model_forward <- model_2ae

summary(model_forward)


# ------------------------------------------------------------------------------
# ----- 5. Model Outputs -------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 5.1. Marginal Means & Contrasts ----------------------------------------

# Average virus effects
emm_virus <- emmeans(model_forward, ~ Virus, type = "response", weights = "equal")
emm_virus

contrast(emm_virus, method = "pairwise")

# Average sex effects
emm_sex <- emmeans(model_forward, ~ Sex, type = "response", weights = "equal")
emm_sex

contrast(emm_sex, method = "pairwise")

# Average age effects
emm_age <- emmeans(model_forward, ~ Age_policy, type = "response", weights = "equal")
emm_age

contrast(emm_age, method = "pairwise")

# Average vaccine effects
emm_vacc <- emmeans(model_forward, ~ vacc_group, type = "response", weights = "equal")
emm_vacc

contrast(emm_vacc, method = "pairwise")

# Virus effects within vaccine groups
emm_virus_vacc <- emmeans(model_forward, ~ Virus | vacc_group, type = "response", weights = "equal")
emm_virus_vacc

contrast(emm_virus_vacc, method = "pairwise")

# Vaccine effects within viruses
emm_vacc_virus <- emmeans(model_forward, ~ vacc_group | Virus, type = "response", weights = "equal")
emm_vacc_virus

contrast(emm_vacc_virus, method = "pairwise")

# Age effects within vaccine groups
emm_age_vacc <- emmeans(model_forward, ~ Age_policy | vacc_group, type = "response", weights = "equal")
emm_age_vacc

contrast(emm_age_vacc, method = "pairwise")

# Vaccine effects within age groups
emm_vacc_age <- emmeans(model_forward, ~ vacc_group | Age_policy, type = "response", weights = "equal")
emm_vacc_age

contrast(emm_vacc_age, method = "pairwise")

# ----- 5.2. Figure 4A ---------------------------------------------------------

n_all <- count(ungroup(data_all), Virus)

p1 <- ggplot(data_all) +
  geom_point(aes(x = Virus, y = u_Neut * 100), position = position_jitter(width = 0.2, height = 0), color = "pink", size = 0.8) +
  geom_boxplot(aes(x = Virus, y = u_Neut * 100, fill = Virus), alpha = 0.75, outlier.shape = NA) +
  geom_text(data = n_all, aes(x = Virus, y = -7, label = n), inherit.aes = FALSE, size = 3, color = "#4D4D4D") +
  theme_base +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous(limits = c(-10, NA), minor_breaks = seq(0, 100, 5), guide = guide_axis(minor.ticks = TRUE)) +
  scale_fill_manual(values = c("#4677AA", "#DECC76", "#CC6476")) +
  labs(x = "All Sera", y = "Neutralisation (%)")

ggsave(here("figures", "fig4A.svg"), plot = p1, dpi = 300, width = 4, height = 2)

contrast(emm_virus, method = "pairwise")

# ----- 5.3. Figure 4B ---------------------------------------------------------

data_all <- data_all %>%
  mutate(Sex = factor(Sex, levels = c("Female", "Male")),
         VirusSex = interaction(Virus, Sex, sep = " × "))

n_all <- count(ungroup(data_all), VirusSex)

p2 <- ggplot(data_all) +
  geom_point(aes(x = VirusSex, y = u_Neut * 100), position = position_jitter(width = 0.2, height = 0), color = "pink", size = 0.8) +
  geom_boxplot(aes(x = VirusSex, y = u_Neut * 100, fill = Virus), alpha = 0.75, outlier.shape = NA) +
  geom_text(data = n_all, aes(x = VirusSex, y = -7, label = n), inherit.aes = FALSE, size = 3, color = "#4D4D4D") +
  theme_base +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        aspect.ratio = 1/2) +
  scale_y_continuous(limits = c(-10, NA), minor_breaks = seq(0, 100, 5), guide = guide_axis(minor.ticks = TRUE)) +
  scale_fill_manual(values = c("#4677AA", "#DECC76", "#CC6476")) +
  labs(x = "Virus × Sex", y = "Neutralisation (%)")

ggsave(here("figures", "fig4B.svg"), plot = p2, dpi = 300, width = 6, height = 2)

model_sex_interaction <- update(model_forward, . ~ . + Virus:Sex)

emm_sex <- emmeans(model_sex_interaction, ~ Sex, type = "response", weights = "equal")
emm_sex

contrast(emm_sex, method = "pairwise")

emm_sex_within_virus <- emmeans(model_sex_interaction, ~ Sex | Virus, type = "response", weights = "equal")
emm_sex_within_virus

contrast(emm_sex_within_virus, method = "pairwise")

emm_virus_within_sex <- emmeans(model_sex_interaction, ~ Virus | Sex, type = "response", weights = "equal")
emm_virus_within_sex

contrast(emm_virus_within_sex, method = "pairwise")

# ----- 5.4. Figure 4C ---------------------------------------------------------

data_all <- data_all %>%
  mutate(VirusAge = interaction(Virus, Age_policy, sep = " × "))

n_all <- data_all %>%
  ungroup() %>%
  count(VirusAge)

p3 <- ggplot(data_all) +
  geom_point(aes(x = VirusAge, y = u_Neut * 100), position = position_jitter(width = 0.2, height = 0), color = "pink", size = 0.8) +
  geom_boxplot(aes(x = VirusAge, y = u_Neut * 100, fill = Virus), alpha = 0.75, outlier.shape = NA) +
  geom_text(data = n_all, aes(x = VirusAge, y = -7, label = n), inherit.aes = FALSE, size = 3, color = "#4D4D4D") +
  theme_base +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        aspect.ratio = 1/3) +
  scale_y_continuous(limits = c(-10, NA), minor_breaks = seq(0, 100, 5), guide = guide_axis(minor.ticks = TRUE)) +
  scale_fill_manual(values = c("#4677AA", "#DECC76", "#CC6476")) +
  labs(x = "Virus × Age", y = "Neutralisation (%)")

ggsave(here("figures", "fig4C.svg"), plot = p3, dpi = 300, width = 8, height = 2)

model_age_interaction <- update(model_forward, . ~ . + Virus:Age_policy)

emm_age <- emmeans(model_age_interaction, ~ Age_policy, type = "response", weights = "equal")
emm_age

contrast(emm_age, method = "pairwise")

emm_age_within_virus <- emmeans(model_age_interaction, ~ Age_policy | Virus, type = "response", weights = "equal")
emm_age_within_virus

contrast(emm_age_within_virus, method = "pairwise")

emm_virus_within_age <- emmeans(model_age_interaction, ~ Virus | Age_policy, type = "response", weights = "equal")
emm_virus_within_age

contrast(emm_virus_within_age, method = "pairwise")

# ----- 5.4. Figure 4D ---------------------------------------------------------

data_all <- data_all %>%
  mutate(VirusVacc = interaction(Virus, vacc_group, sep = " × "))

n_all <- data_all %>%
  ungroup() %>%
  count(VirusVacc)

p4 <- ggplot(data_all) +
  geom_point(aes(x = VirusVacc, y = u_Neut * 100), position = position_jitter(width = 0.2, height = 0), color = "pink", size = 0.8) +
  geom_boxplot(aes(x = VirusVacc, y = u_Neut * 100, fill = Virus), alpha = 0.75, outlier.shape = NA) +
  geom_text(data = n_all, aes(x = VirusVacc, y = -7, label = n), inherit.aes = FALSE, size = 3, color = "#4D4D4D") +
  theme_base +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        aspect.ratio = 1/3) +
  scale_y_continuous(limits = c(-10, NA), minor_breaks = seq(0, 100, 5), guide = guide_axis(minor.ticks = TRUE)) +
  scale_fill_manual(values = c("#4677AA", "#DECC76", "#CC6476")) +
  labs(x = "Virus × Years Since Last Vaccination", y = "Neutralisation (%)")

ggsave(here("figures", "fig4D.svg"), plot = p4, dpi = 300, width = 10, height = 2)

emm_vacc <- emmeans(model_forward, ~ vacc_group, type = "response", weights = "equal")
emm_vacc

contrast(emm_vacc, method = "pairwise")

emm_virus_vacc <- emmeans(model_forward, ~ Virus | vacc_group, type = "response", weights = "equal")
emm_virus_vacc

contrast(emm_virus_vacc, method = "pairwise")

emm_vacc_virus <- emmeans(model_forward, ~ vacc_group | Virus, type = "response", weights = "equal")
emm_vacc_virus

contrast(emm_vacc_virus, method = "pairwise")




