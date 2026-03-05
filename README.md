# Humanized_mice_EAE

This repository contains the code and resources for the **"Humanized Mice EAE"** project. It provides a comprehensive workflow for single-cell RNA-seq data analysis, figure generation, and preprocessing.

## Contents

1. **Vignette in RMarkdown**: 
   - A detailed walkthrough of the single-cell RNA-seq data analysis pipeline.
   - Includes code for figure generation used in the manuscript.
   - Reproducible and easy to follow for researchers.

2. **Snakemake Workflow**:
   - A Snakemake script for preprocessing single-cell RNA-seq data using Cell Ranger.
   - Automates the alignment, quantification, and quality control steps.

## Project Structure

- **Preprocessing/**: Contains Snakemake workflows for data preprocessing
  - **config/**: Configuration files for the workflow
  - **workflow/**: Snakemake rules and execution scripts
- **analysis.Rmd**: RMarkdown file for downstream analysis and visualization
- **functions.R**: Utility functions for analysis

## Usage

To run the preprocessing workflow:
1. Update the configuration in `Preprocessing/config/config.yaml`
2. Execute Snakemake: `snakemake -s Preprocessing/workflow/Snakefile_preprocess_RNA.smk --cores N`

## Requirements

- Snakemake
- Cell Ranger
- R with required packages (Seurat, ggplot2, etc.)
