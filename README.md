# 🎯 pHLATargetTcellTherapy

## Overview

This repository contains a reproducible bioinformatics pipeline for **identifying and prioritizing potential targets for TCR-T cell therapy** by integrating large-scale transcriptomic datasets from **The Cancer Genome Atlas (TCGA)** and the **Genotype-Tissue Expression (GTEx)** project.

The workflow automates the full discovery process:
1. **Data integration** — harmonizing TCGA and GTEx RNA-Seq expression profiles into a unified, analysis-ready dataset.
2. **Expression cutoff modeling** — calculating both absolute and relative expression cutoffs from normal tissues to minimize off-target toxicity risk.
3. **Tumor target assessment** — quantifying tumor-specific expression across multiple cancer indications.
4. **Risk–benefit scoring** — integrating safety profiles (tissue risk categories) with efficacy metrics to produce a ranked list of candidate targets.
5. **Decision support** — enabling transparent and data-driven prioritization of targets for downstream validation.

---

## ✨ Key Features

- **Reproducible** — Entire analysis is orchestrated via `Makefile` for easy re-runs and automatic step skipping.
- **Data-driven safety thresholds** — Combines absolute max and IQR-based relative max cutoffs for robust expression filtering.
- **Risk-aware** — Incorporates tissue-specific risk stratification for improved translational relevance.
- **Cross-indication utility scoring** — Balances efficacy and safety to allow fair target comparisons across cancer types.
- **Modular** — Each step (R or Python) can be run independently or as part of the full workflow.

---

## Why This Matters

TCR-T cell therapy shows immense promise for precision oncology, but off-target toxicity remains a critical challenge.  
This pipeline addresses that challenge by systematically integrating **tumor expression patterns** with **normal tissue safety profiles** to highlight targets with **high tumor specificity and minimal normal tissue risk**.  
The approach is designed for translational researchers, computational biologists, and immunotherapy teams looking to accelerate target discovery without compromising safety.

---

## 📂 Project Structure
pHLATargetTcellTherapy/ \
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

     - **Utility score** = combined safety + efficacy measure
   - Outputs ranked targets and plots
