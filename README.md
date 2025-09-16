# 🎯 pHLATargetTcellTherapy

This repository contains a simple, reproducible workflow for screening potential **TCR-T cell therapy targets** using processed RNA-Seq expression data from the **UCSC Xena platform** (TCGA and GTEx cohorts).

The pipeline:
1. **Extracts & harmonizes** expression data for tumor and normal tissues from UCSC Xena.
2. **Calculates cutoffs** for tumor vs. normal expression and the percentage of samples above these cutoffs.
3. **Applies risk–benefit scoring** to help prioritize genes for further investigation.

Although minimal in complexity, this workflow:
- Automates the steps from data extraction to ranked target lists.
- Uses clear, modular R and Python scripts.
- Organizes intermediate and final results for reproducibility.
- Can be extended with more advanced statistical modeling or visualization.

---

## 📂 Project Structure
```bash
pHLATargetTcellTherapy/ 
├── Data/ 
│ ├── Input/     # Raw TCGA + GTEx input data, antigen searching space (read-only) 
│ ├── Processed/ # Intermediate outputs from R scripts 
│ └── Reference/ # Tissue risk metadata, gene annotation files 
├── Scripts/ 
│ ├── 01_extract_expression_data.R 
│ ├── 02_cal_cutoff_percentage_above_cutoff.R 
│ ├── 03_risk_benefit_pipeline.py 
├── Output/ 
│ ├── Tables/   # Final ranked tables per gene/indication 
│ └── Figures/  # Plots and visualizations 
├── Makefile    # Workflow orchestration 
├── config.yaml # Pipeline parameters (cutoffs, risk weights, etc.) 
└── README.md   # This file 
```

---

## 🔄 Workflow Overview

1. **Extract expression data** (R)  
   - Reads tumor and normal RNA-Seq expression from UCSC Xena TCGA + GTEx datasets.
   - Outputs harmonized expression matrix into `Data/Processed/`.

2. **Calculate cutoffs & percentages** (R)  
   - Defines cutoffs per gene using both absolute and relative normal tissue expression metrics.
   - Calculates the percentage of samples above the cutoff for each tumor and normal tissue type.

3. **Risk–benefit scoring** (Python)  
   - Reads cutoff results and tissue risk categories.
   - Calculates:
     - **Benefit score** (per tumor indication)
     - **Risk score** (global per gene)
     - **Utility score** = combined measure of safety & efficacy
   - Outputs ranked targets and summary plots.

---

## ⚙️ Running the Pipeline

### Requirements
- **R** (≥ 4.0) with packages:
  - `data.table`
- **Python** (≥ 3.9) with packages:
  - `pandas`
  - `numpy`
  - `matplotlib`
  - `pyyaml`

### Run with Makefile
From the project root:
```bash
make
```
This will run:

01_extract_expression_data.R

02_cal_cutoff_percentage_above_cutoff.R

03_risk_benefit_pipeline.py

Intermediate results go into Data/Processed/
Final outputs are stored in Output/

---

## 🧮 Parameters

Pipeline settings are stored in config.yaml:
```
cutoff_percentile: 0.95
relative_k: 1.5
lambda: 0.5
tissue_risk_file: Data/Reference/tissue_risk_weights.csv
```
---

## 📊 Outputs

- Tables: Output/Tables/final_ranked_targets.csv — ranked list of targets per indication

- Figures: Output/Figures/risk_vs_benefit_scatter.png — safety vs. efficacy overview
