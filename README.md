# rcca-multiomics-pipeline
R pipeline for regularized CCA (rCCA) across multi-omics datasets, including feature coordinate export and variable-circle visualization.

This repository contains an R-based pipeline for regularized Canonical Correlation Analysis (rCCA) across multiple omics datasets. The workflow includes feature coordinate export and variable-circle visualization to support interpretation of cross-omics relationships.

## Repository contents

- `01_export_rcca_coordinates.R`  
  Fits rCCA across dataset pairs, computes variable-circle coordinates, and exports feature-level coordinate tables.

- `02_plot_variable_circles.R`  
  Generates variable-circle plots from the exported coordinates and saves them as PDF files.

## Main dependencies

- `mixOmics`
- `ggplot2`
- `ggrepel`
- `tools`

## Expected input files

The scripts expect CSV input matrices such as:

- `MetaTranscriptomics_CCA.csv`
- `BileMetabolites_CCA.csv`
- `MetaGenomics_CCA.csv`
- `RNAseq_CCA.csv`
- `16SMicrobiom_CCA.csv`

## Outputs

The pipeline generates:

- exported feature coordinate CSV files
- labeled feature subsets based on radius cutoffs
- variable-circle PDF plots
- run information file with session details

## Notes

- This repository contains analysis code only.
- Input data are not included in this repository.
- The workflow was used in a collaborative multi-omics analysis setting and is shared here as part of my computational research portfolio.
