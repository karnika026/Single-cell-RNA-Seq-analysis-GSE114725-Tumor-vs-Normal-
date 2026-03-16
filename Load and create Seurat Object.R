library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)l
library(readr)
library(purrr)

#Loading files
file_path <- "........Github project/GSE114725/GSE114725_rna_raw.csv.gz"

# Read counts, keeping rownames
counts_df <- read.csv(
  gzfile(file_path, "rt"),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

#Separate Metadata and Expression
meta_cols <- c("tissue", "replicate", "cluster", "cellid")
stopifnot(all(meta_cols %in% colnames(counts_df)))
cell_metadata <- counts_df[, meta_cols]
rownames(cell_metadata) <- counts_df$cellid
expr_df <- counts_df[, !colnames(counts_df) %in% meta_cols]


#Reformat expression matrix (genes* cells)
expr_mat <- t(as.matrix(expr_df))
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- make.unique(colnames(expr_df))
colnames(expr_mat) <- counts_df$cellid

# Quick sanity check
dim(expr_mat)
head(rownames(expr_mat))
head(colnames(expr_mat))

#Creating Seurat Object
seurat_object <- CreateSeuratObject(
  counts = expr_mat,
  project = "GSE114725",
  min.cells = 3,
  min.features = 200
)

#Add Metadata
seurat_object <- AddMetaData(seurat_object, cell_metadata)
head(seurat_object@meta.data)
table(seurat_object$tissue)




