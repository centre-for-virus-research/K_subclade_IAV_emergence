
# (2026) Understanding the emergence of the influenza A/H3N2 K subclade in its historical and evolutionary context

This repository contains data and analysis scripts from **(2026) Dee K, Imrie RM, MacLean OA, Mojsiejczuk L, Smith EW, Raveendran S, Lamb K, Chen H, Schultz V, Wang Z, Walsh SK. _Understanding the emergence of the influenza A/H3N2 K subclade in its historical and evolutionary context_**.

 https://www.biorxiv.org/content/biorxiv/early/2026/05/24/2026.05.21.726823.full.pdf

## Using this repository

### Platform support
![Windows](https://img.shields.io/badge/Windows-blue?logo=microsoftwindows)
![macOS](https://img.shields.io/badge/macOS-black?logo=apple)
![Linux](https://img.shields.io/badge/Linux-grey?logo=linux)

### Dependencies

![R Version](https://img.shields.io/badge/R-4.5.2-2980b9)
![tidyverse](https://img.shields.io/badge/tidyverse-2.0.0-bbd4f1)
![here](https://img.shields.io/badge/here-1.0.2-bbd4f1)
![patchwork](https://img.shields.io/badge/patchwork-1.3.2-bbd4f1)
![minpack.lm](https://img.shields.io/badge/minpack.lm-1.2.4-bbd4f1)
![plotrix](https://img.shields.io/badge/plotrix-3.8.13-bbd4f1)

![betareg](https://img.shields.io/badge/betareg-3.2.4-1abc9c)
![emmeans](https://img.shields.io/badge/emmeans-2.0.1-1abc9c)

### Running R Scripts

The included script scripts/00_setup.R can be used to install all R package dependencies at once.

Scripts in this repository use the here library to dynamically set paths. For this to work correctly, Rstudio must be opened by double-clicking on one of the files in scripts/. Path errors will appear if Rstudio was first opened using a shortcut or a script from a different location.

## Contents

| Item | Description |
|---|---|
| `data/` | Data used in the analyses |
| └─ `neutralisation_75.csv` | Raw neutralisation assay data measured at a 1:75 serum dilution |
| └─ `sample_metadata.csv` | Serum sample metadata and vaccination history |
| └─ `IC50s.csv` | Serial-dilution neutralisation data used for IC50 curve fitting |
| └─ `IC50_fits.csv` | Fitted logistic-curve parameters and IC50 estimates |
| └─ `gisaid_260126pn.pdf` | GISAID record and DOI information |
| `scripts/` | R scripts used for analysis |
| └─ `00_setup.R` | Installs the R package dependencies used in this repository |
| └─ `01_neutralisation_analysis.R` | Analyses 1:75 neutralisation data using beta regression |
| └─ `02_IC50_fitting.R` | Fits logistic curves to estimate IC50 values |
| └─ `03_assay_comparisons.R` | Compares neutralisation and IC50 assays |
