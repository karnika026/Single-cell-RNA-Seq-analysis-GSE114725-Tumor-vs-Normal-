# 3' scRNA-Seq-analysis of GSE114725-Tumor-vs-Normal
This repository contains scripts to reproduce and extend the analysis of GSE114725, a 3′ single-cell RNA-seq dataset profiling tumor, normal, blood, and lymph node samples. The workflow includes loading raw counts, creating a Seurat object, and preparing for QC, clustering, and tumor vs normal analysis.

- **Accession:** GSE114725  
- **Data type:** 3′ single-cell RNA-seq  
- **Biological context:** Tumor, normal, blood, and lymph node samples  
- **Source:** NCBI Gene Expression Omnibus (GEO)

## Files used:
- `GSE114725_rna_raw.csv.gz`

### Workflow Overview

1. Load raw gene expression matrix
2. Create Seurat object
3. Quality control and filtering
4. Subset tumor and normal cells
5. Normalize and scale data
6. Differential expression(FindMarkers): Tumor vs Normal

### Differential Expression Strategy

Tumor vs normal differential expression is performed using Seurat’s `FindMarkers()` function after subsetting to tumor and normal cells only.
Normalized (log-transformed) expression values are used, following standard Seurat practice.

## Repository structure

GSE114725/
├── scripts/
│   ├── 01_Load_Data.R
│   ├── 02_QC_Filtering.R
│   ├── 03_Subset_Normal_Tumor.R
│   └── 04_Tumor_vs_Normal_DE.R
├── plots/
└── README.md

### Data Availability
Raw data are available from GEO under accession **GSE114725**. No raw data files are stored in this repository.

