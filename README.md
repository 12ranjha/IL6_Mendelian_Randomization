# IL6_Mendelian_Randomization
Drug-target Mendelian Randomization study analyzing genetically proxied IL-6 inhibition on HCC and CRC risk.

# Genetically Proxied IL-6 Receptor Blockade and Cancer Risk

### A Multi-Ancestry Drug-Target Mendelian Randomization Study of Hepatocellular Carcinoma and Colorectal Cancer

Muhammad Muneeb Ahmad Ranjha¹, Hamna Munir², Muhammad Saleem³*, Muzamil Farooq¹
¹ Department of Medicine, King Edward Medical University, Lahore, Pakistan
 ² Department of Computer Science, COMSATS University Islamabad, Lahore Campus, Pakistan
 ³ Medical Faculty, Nangarhar University, Jalalabad, Afghanistan 

---

## Overview

This repository contains the complete codebase, master results tables, and raw intermediate harmonized datasets for our trans-ancestry drug-target Mendelian randomization (MR) study. The project evaluates the lifelong causal implications of genetically proxied Interleukin-6 (IL-6) signaling inhibition (modeling the therapeutic mechanism of the monoclonal antibody tocilizumab) on the risk of Hepatocellular Carcinoma (HCC) and Colorectal Cancer (CRC) across European and East Asian populations.

By providing full access to our raw harmonized data matrices alongside our production execution scripts, this project achieves 100% compliance with the **STROBE-MR** guidelines for open science and reproducible genetic epidemiology.

---

## Directory & File Structure

```text
IL6_Mendelian_Randomization/ (Root Directory)
├── README.md                           # This documentation file
├── il6_mr_pipeline.R                   # Complete automated R execution script
├── Table1_Instruments_Metrics.csv      # Genetic instrument metrics and F-statistics
├── Table2_Master_MR_Estimates.csv     # Main MR causal estimates across all methods
├── Table3_Master_Heterogeneity.csv    # Cochran's Q heterogeneity test results
├── Table4_Master_Pleiotropy.csv       # MR-Egger horizontal pleiotropy intercept tests
└── myfolder/                           # Folder housing raw intermediate data frames
    ├── Raw_Harmonized_Data_CRC_ASN.csv # Raw data matrix: East Asian CRC (BBJ)
    ├── Raw_Harmonized_Data_CRC_EUR.csv # Raw data matrix: European CRC (Huyghe et al.)
    ├── Raw_Harmonized_Data_HCC_ASN.csv # Raw data matrix: East Asian HCC (BBJ)
    └── Raw_Harmonized_Data_HCC_EUR.csv # Raw data matrix: European HCC (FinnGen R10)

## Data Description
Master Results
 TablesTable1_Instruments_Metrics.csv: Contains coordinates, effect alleles, frequencies, betas, standard errors, and derived F-statistics for the 4 independent cis-pQTL instrument variants (rs2228145, rs4129267, rs7529229, rs1800795) extracted from plasma proteomics summary statistics.
  Table2_Master_MR_Estimates.csv: Contains the full causal output metrics (Odds Ratios, 95% Confidence Intervals, and P-values) for the Primary Inverse Variance Weighted (IVW), Weighted Median, and MR-Egger regression models across all four cohort pathways.
  Table3_Master_Heterogeneity.csv & Table4_Master_Pleiotropy.csv: Detail secondary sensitivity bounds checking, confirming the absolute absence of directional horizontal pleiotropy (p > 0.05).
Raw Harmonized Data Matrices (Stored in myfolder/)These files contain the exact intermediate data frames after allele matching, strand alignment, and palindromic variant handling via the TwoSampleMR protocol (action=2).
Raw_Harmonized_Data_HCC_EUR.csv uniquely documents the automated exclusion of the palindromic variant rs1800795 due to intermediate allele frequencies in the FinnGen dataset, demonstrating our strict data hygiene safeguards before running regressions.
 Getting Started & ExecutionPrerequisitesTo execute the computational pipeline, you will need an active installation of R (v4.6 or later) and the following packages:
install.packages("devtools")
install.packages("ggplot2")
devtools::install_github("MRCIEU/TwoSampleMR")
Running the Analysis
The execution script il6_mr_pipeline.R is built to run both online (querying live servers) and completely offline using the provided local harmonized files inside myfolder/[cite: 2, 4]. To reproduce our exact findings and re-generate all publication-grade diagnostic panels (Scatter, Forest, Funnel, and Leave-One-Out plots at 300 DPI), simply execute:

R
source("il6_mr_pipeline.R")

Open Science and Code Availability Statement
This repository is configured to insulate peer-review tracking from external API server maintenance or downtime[cite: 2, 4]. For questions, data inquiries, or collaboration requests regarding our genetic epidemiology frameworks, please open an issue in this repository or contact the corresponding author[cite: 2, 4].
