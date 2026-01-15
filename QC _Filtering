#QC and Filtering
seurat_object <- readRDS("......../GSE114725")
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
plots <- VlnPlot(
  seurat_object,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.6,
  combine = FALSE
)
plots <- lapply(plots, function(p) {
  p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}) 
patchwork::wrap_plots(plots, ncol = 3)
print(plots[[1]])
print(plots[[2]])
print(plots[[3]])
ggsave(
  filename = "nFeature_RNA.png",
  plot = plots[[1]],
  path = "......./GSE114725/PLOTS",
  width = 15,
  height = 10,
  dpi = 300
)
ggsave(
  filename = "nCount_RNA.png",
  plot = plots[[2]],
  path = "........./GSE114725/PLOTS",
  width = 15,
  height = 10,
  dpi = 300
)
ggsave(
  filename = "percent.mt.png",
  plot = plots[[3]],
  path = "............/GSE114725/PLOTS",
  width = 15,
  height = 10,
  dpi = 300
)

plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", raster=FALSE)
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=FALSE)
plot = plot1 + plot2
print(plot1)
ggsave(filename = "plot1.png",
  plot = plot1,
  path = "........./GSE114725/PLOTS",
  width = 15,
  height = 10,
  dpi = 300
)
ggsave(filename = "plot2.png",
  plot = plot2,
  path = "............../GSE114725/PLOTS",
  width = 15,
  height = 10,
  dpi = 300
)
ggsave(filename = "plot.png",
  plot = plot,
  path = ".........../GSE114725/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)
seurat_object_qc <- subset(
  seurat_object,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    nCount_RNA < 20000 &
    percent.mt < 20
)

dim(seurat_object_qc)















