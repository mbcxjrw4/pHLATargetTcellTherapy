# pHLATargetTcellTherapy
This project implements a reproducible workflow for screening potential TCR-T cell therapy targets using combined **TCGA** and **GTEx** expression data. The method applies tumor vs. normal tissue expression cutoffs, calculates safety/efficacy scores, and ranks gene targets by a combined utility score.

---

## 📂 Project Structure
project_root/ \
├── Data/ \
│ ├── Input/ # Raw TCGA + GTEx input data (read-only) \
│ ├── Processed/ # Intermediate outputs from R scripts \
│ └── Reference/ # Tissue risk metadata, gene annotation files \
├── Scripts/ \
│ ├── 01_extract_expression_data.R \
│ ├── 02_cal_cutoff_percentage_above_cutoff.R \
│ ├── 03_risk_benefit_pipeline.py \
├── Output/ \
│ ├── Tables/ # Final ranked tables per gene/indication \
│ └── Figures/ # Plots and visualizations \
├── Makefile # Workflow orchestration \
├── config.yaml # Pipeline parameters (cutoffs, risk weights, etc.) \
└── README.md # This file \

---

## 🔄 Workflow Overview

1. **Extract expression data** (R)  
   - Reads TCGA + GTEx data from `Data/Input/`
   - Produces combined expression matrix in `Data/Processed/`

2. **Calculate cutoffs and % above cutoff** (R)  
   - Defines cutoffs per gene based on normal tissue expression  
   - Produces per-tissue `% above cutoff` tables

3. **Risk–benefit scoring** (Python)  
   - Reads cutoffs and tissue risk categories  
   - Calculates:
     - **Benefit score** (per tumor indication)
     - **Risk score** (global per gene)
     - **Utility score** = combined safety + efficacy measure
   - Outputs ranked targets and plots
