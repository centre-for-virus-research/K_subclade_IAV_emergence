
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
![viridisLite](https://img.shields.io/badge/viridisLite-0.4.2-1abc9c)

![Python](https://img.shields.io/badge/Python-3.12-3776ab?logo=python)
![PyTorch](https://img.shields.io/badge/PyTorch-2.5.1-ee4c2c?logo=pytorch)
![Transformers](https://img.shields.io/badge/Transformers-4.41.2-ffcc4d?logo=huggingface)

### Running R Scripts

The included script `scripts/00_setup.R` can be used to install all R package dependencies at once.

Scripts in this repository use the `here` library to dynamically set paths. For this to work correctly, RStudio must be opened by double-clicking on one of the files in `scripts/`. Path errors will appear if RStudio was first opened using a shortcut or a script from a different location.

## Contents

| Item | Description |
|---|---|
| `data/` | Data used in the analyses |
| └─ `neutralisation_75.csv` | Raw neutralisation assay data measured at a 1:75 serum dilution |
| └─ `sample_metadata.csv` | Serum sample metadata and vaccination history |
| └─ `IC50s.csv` | Serial-dilution neutralisation data used for IC50 curve fitting |
| └─ `IC50_fits.csv` | Fitted logistic-curve parameters and IC50 estimates |
| └─ `gisaid_260126pn.pdf` | GISAID record and DOI information |
| └─ `PLANT_input_file.csv` | Curated influenza A/H3N2 HA sequences and subclades used for PLANT antigenic cartography |
| `scripts/` | R scripts used for analysis |
| └─ `00_setup.R` | Installs the R package dependencies used in this repository |
| └─ `01_neutralisation_analysis.R` | Analyses 1:75 neutralisation data using beta regression |
| └─ `02_IC50_fitting.R` | Fits logistic curves to estimate IC50 values |
| └─ `03_assay_comparisons.R` | Compares neutralisation and IC50 assays |
| `Plant.run.py` | Projects A/H3N2 HA sequences into 3D antigenic space, computes distances to reference strains, and creates interactive 3D visualizations |
| `Plant_batch_fastas.py` | Batch CLI utility to align, QC, and compute 3D PLANT antigenic embeddings for directories of FASTA files |
| `Functions_HuggingFace.py` | Sequence alignment, reference-based trimming, and quality control helper utilities |
| `external/PLANT` | Git submodule pointing to [TheSatoLab/PLANT](https://github.com/TheSatoLab/PLANT) (model architecture and historical background dataset) |
| `requirements.txt` | Pinned Python package dependencies for PLANT environment reproduction |

---

## PLANT Antigenic Cartography

[PLANT](https://github.com/TheSatoLab/PLANT) (**P**rotein **L**anguage model for **A**ntigenic car**T**ography; Ito et al.) projects influenza A/H3N2 hemagglutinin (HA) protein sequences directly into a 3D latent antigenic space using an ESM-2 protein language model backbone (`facebook/esm2_t33_650M_UR50D`) coupled with specialized projection heads. In this study, PLANT is used to quantify antigenic drift, compute antigenic distances from reference strains, and contextualize the emergence of the A/H3N2 K subclade relative to contemporaneous lineages and historical influenza strains from 1968–2024.

### 1. Repository Submodule Setup

This repository links the core PLANT source code and historical reference dataset via a Git submodule under `external/PLANT`.

- **When cloning for the first time:**
  ```bash
  git clone --recurse-submodules https://github.com/centre-for-virus-research/K_subclade_IAV_emergence.git
  ```
- **If you already cloned without submodules:**
  ```bash
  git submodule update --init --recursive
  ```

> [!NOTE]
> The scripts in this repository automatically detect and load `external/PLANT/src` (for the `plant` package) and `external/PLANT/examples/backgrounds.csv` (historical strains dataset), eliminating the need to configure manual `PYTHONPATH` exports.

---

### 2. PLANT Model Weights from Hugging Face

Pretrained model weights and encoders are hosted on Hugging Face at:
👉 **[TheSatoLab-UTokyo/PLANT](https://huggingface.co/TheSatoLab-UTokyo/PLANT)**

The scripts use the improved checkpoint located under the **`variants/PLANT_fixed`** subfolder, which consists of:
- `model.safetensors` (fine-tuned ESM-2 650M backbone weights)
- `virus_encoder.joblib`, `ref_encoder.joblib`, `vp_encoder.joblib`, `rp_encoder.joblib` (categorical encoders)
- `config.json`

#### Downloading the Weights Locally

You can download the checkpoint using the Hugging Face CLI or Python into `./models/PLANT_model`:

```bash
# Option A: Using huggingface-cli
pip install huggingface_hub
huggingface-cli download TheSatoLab-UTokyo/PLANT \
  --include "variants/PLANT_fixed/*" \
  --local-dir ./models/PLANT_model
```

Or via Python:
```python
# Option B: Using Python
from huggingface_hub import snapshot_download

snapshot_download(
    repo_id="TheSatoLab-UTokyo/PLANT",
    allow_patterns=["variants/PLANT_fixed/*"],
    local_dir="./models/PLANT_model",
)
```

The scripts look for this checkpoint by default in `./models/PLANT_model` (or custom locations specified by `export PLANT_MODEL_DIR="/path/to/PLANT_model"`).

---

### 3. Python Environment & Dependencies

A tested, pinned [requirements.txt](file:///home3/oml4h/K_subclade_IAV_emergence/requirements.txt) is provided in the repository root. To recreate the environment:

```bash
# 1. Create and activate a dedicated conda environment
conda create -n plant_env python=3.12 -y
conda activate plant_env

# 2. Install PyTorch with CUDA support and pinned dependencies
pip install -r requirements.txt --extra-index-url https://download.pytorch.org/whl/cu124
```

#### Optional Environment Variables

The scripts support the following environment variables if you store models or datasets in non-default directories:
- `PLANT_MODEL_DIR`: Directory containing downloaded Hugging Face PLANT checkpoints (default: `./models/PLANT_model` containing `variants/PLANT_fixed`).
- `BACKGROUND_CSV_PATH`: Path to historical strains dataset (default: `./external/PLANT/examples/backgrounds.csv`).
- `PLANT_OUTDIR`: Output directory for generated CSVs and HTML visualizations (default: `./results/PLANT_results`).

---

### 4. Scripts Description & Usage

#### 1. `Plant.run.py` — Sequence Embedding, Distance Calculation & Interactive Visualization

`Plant.run.py` carries out the end-to-end embedding pipeline for a curated set of HA sequences (by default reading `data/PLANT_input_file.csv`):

1. **Alignment & QC Trimming:**
   - Standardises each input sequence to the 329 amino acid HA1 domain using pairwise local-global alignment against a J.2 reference sequence (`Functions_HuggingFace.py`).
   - Discards ambiguous/invalid amino acid residues (`X`, `B`, `*`).
   - Applies strict alignment quality thresholds: minimum sequence identity of 80% (warns below 90%) and minimum reference coverage of 95% (warns below 98%).
2. **3D Latent Embedding:**
   - Tokenizes trimmed sequences with the ESM-2 tokenizer and passes them through the PLANT model.
   - Computes 3D coordinates ($X, Y, Z$) scaled by a factor of 8 to match standard antigenic cartography scaling.
3. **Reference Distance Quantification:**
   - Computes Euclidean distance in 3D antigenic space from the J.2 reference strain (`EPI2981619|HA|A/Croatia/10136RV/2023|EPI_ISL_18856647|J.2`) to evaluate relative antigenic drift.
4. **Interactive 3D Visualizations (Plotly HTML):**
   - `PLANT_embedding_subclade.html`: 3D scatter plot colored by subclade/lineage.
   - `PLANT_embedding_subclade_emphasized.html`: Lineage-focused view highlighting key emerging clades (`K`, `J.2.2Eng24`, and `J.2.vaccine_column`) against a muted background of other lineages.
   - `PLANT_embedding_year.html`: Antigenic map colored continuously by sample collection year.
   - `PLANT_embedding_background_subclade.html` & `PLANT_embedding_background_year.html`: Overlays the sequences of interest onto ~150,000 historical A/H3N2 background strains (1968–2024) from the PLANT reference dataset (`external/PLANT/examples/backgrounds.csv`).
   - `yearly_centroids.csv`: Exports annual mean coordinate centroids of circulating strains over time.

**Execution:**
```bash
python Plant.run.py
```

**Outputs generated:**
- `PLANT_input_filtered.csv`: Quality-filtered and trimmed sequence dataset.
- `PLANT_embeddings_with_distance.csv`: Metadata and 3D coordinates ($X, Y, Z$) with Euclidean distance to the reference strain.
- `yearly_centroids.csv`: Table of historical yearly antigenic centroids and sequence counts.
- Interactive 3D visualization files in HTML format (`*.html`).

---

#### 2. `Plant_batch_fastas.py` — High-Throughput Batch FASTA Embedding CLI

`Plant_batch_fastas.py` is a command-line utility designed to process multiple FASTA files in batch mode (e.g., partitioned by lineage, year, or geographic region) and output coordinate tables per FASTA:

- **Features:**
  - Scans an input directory for FASTA files matching specified patterns (`*.fa`, `*.fasta`, `*.faa`).
  - Automatically loads model weights and encoders once into GPU memory (with FP16 inference) or CPU.
  - Aligns, trims, and quality-filters all sequences against the 329 aa HA1 reference.
  - Batches sequence tokenization and embedding using PyTorch `DataLoader` (default batch size: 64).
  - Writes a separate CSV file (`<fasta_name>.csv`) for each input FASTA containing `sequence_id`, `X`, `Y`, and `Z` coordinates.

**Usage:**
```bash
python Plant_batch_fastas.py \
  --input-dir /path/to/fasta_directory \
  --output-dir /path/to/output_directory \
  --batch-size 64 \
  --suffixes "*.fa" "*.fasta" "*.faa"
```

**Options:**
- `--input-dir`: Path to folder containing input FASTA files.
- `--output-dir`: Path to directory where output CSV files will be saved.
- `--batch-size`: Batch size for embedding inference (default: `64`).
- `--suffixes`: File extension glob patterns to search for (default: `*.fa *.fasta *.faa`).