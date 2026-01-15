# 3' scRNA-Seq-analysis of GSE114725-Tumor-vs-Normal
This repository contains scripts to reproduce and extend the analysis of GSE114725, a 3′ single-cell RNA-seq dataset profiling tumor, normal, blood, and lymph node samples. The workflow includes loading raw counts, creating a Seurat object, and preparing for QC, clustering, and tumor vs normal analysis.

- **Accession:** GSE114725  
- **Data type:** 3′ single-cell RNA-seq  
- **Biological context:** Tumor, normal, blood, and lymph node samples  
- **Source:** NCBI Gene Expression Omnibus (GEO)

## Files used:
- `GSE114725_rna_raw.csv.gz`


## Repository structure

```text
.
├── README.md
└── scripts
    └── 01_load_and_create_seurat.R
    └── 02_QC and Filtering

