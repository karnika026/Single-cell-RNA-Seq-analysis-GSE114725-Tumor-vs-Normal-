#Creating a subset of Normal and Tumor samples
seurat_data <- subset(seurat_object_qc, subset = tissue %in% c("TUMOR", "NORMAL"))
seurat_data

#Normalizing and Scaling Data
seurat_data <- NormalizeData(seurat_data)
seurat_data <- FindVariableFeatures(seurat_data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_data), 10)
top10

# plot variable features with and without labels
Variablefeaturesplot <- VariableFeaturePlot(seurat_data)
labelpointplot <- LabelPoints(plot = Variablefeaturesplot, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
Variableplot = Variablefeaturesplot + labelpointplot
print(Variableplot) + theme_github
ggsave(filename = "Variableplot.png",
  plot = Variableplot,
  path = "........./GSE114725/PLOTS",
  width = 20,
  height = 10,
  dpi = 300,
   bg = "white"
)

seurat_data <- ScaleData(seurat_data,features = VariableFeatures(seurat_data))
class(seurat_data)
seurat_data
head(seurat_data@meta.data)
tail(seurat_data@meta.data)
unique(seurat_data$tissue)

#Run PCA
seurat_data <- RunPCA(seurat_data, features = VariableFeatures(object = seurat_data))
Vizdimplot <-VizDimLoadings(seurat_data, dims = 1:2, reduction = "pca")
print(Vizdimplot) + theme_github
ggsave(filename = "Vizdimplot.png",
  plot = Vizdimplot,
  path = ".........../GSE114725/PLOTS",
  width = 20,
  height = 10,
  dpi = 300,
  bg = "white"
)
dimplot <- DimPlot(seurat_data, reduction = "pca", raster = FALSE) + NoLegend()
print(dimplot) + theme_github
ggsave(filename = "dimplot.png",
  plot = dimplot,
  path = "........../GSE114725/PLOTS",
  width = 20,
  height = 10,
  dpi = 300,
  bg = "yellow"
)
heatmap <- DimHeatmap(seurat_data, dims = 1:15, cells = 500, balanced = TRUE)
pdf(
  file = ".........../GSE114725/PLOTS/heatmap.pdf",
  width = 12,
  height = 8
)
dev.off()


#Dimensional heatmaps were saved using base R graphics devices due to ComplexHeatmap-based rendering in Seurat.
elbowplot <- ElbowPlot(seurat_data)
print(elbowplot) + theme_github
ggsave(filename = "elbowplot.png",
  plot = elbowplot,
  path = "........./GSE114725/PLOTS",
  width = 20,
  height = 10,
  dpi = 300,
  bg = "yellow"
)
seurat_data <- FindNeighbors(seurat_data, dims = 1:15)
seurat_data <- FindClusters(seurat_data, resolution = 0.5)
head(Idents(seurat_data), 5)
seurat_data <- RunUMAP(seurat_data, dims = 1:15)
umap_plot <- DimPlot(
  seurat_data,
  reduction = "umap",
  raster = FALSE,
  label = TRUE,
  repel = TRUE
) + theme_github
ggsave(filename = "umap_plot.png",
  plot = umap_plot,
  path = "............./GSE114725/PLOTS",
  width = 20,
  height = 10,
  dpi = 300,
  bg = "white"
)
saveRDS(
  seurat_data,
  file = "........./GSE114725/seurat_data.rds"
)
