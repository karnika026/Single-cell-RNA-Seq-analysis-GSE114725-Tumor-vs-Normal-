# Differential Expression of Tumor vs Normal
Idents(seurat_data) <- "tissue"    # TUMOR vs NORMAL
levels(Idents(seurat_data))
markers <- FindMarkers(
  seurat_data,
   ident.1 = "TUMOR",
  ident.2 = "NORMAL",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)
head(markers)
tail(markers)

markers <- markers %>%
  tibble::rownames_to_column("gene")
head(markers)
markers <- markers %>%
  mutate(
    sig = p_val_adj < 0.05 & abs(avg_log2FC) > 0.25
  )

p_volcano <- ggplot(markers,
                    aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1) +
  scale_color_manual(values = c("grey60", "firebrick")) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Tumor vs Normal Differential Expression",
    x = "Average log2 Fold Change",
    y = "-log10 adjusted p-value"
  ) +
  theme_classic(base_size = 14)

ggsave(
  filename = "plots/TUMOR_vs_NORMAL_volcano.png",
  plot = p_volcano,
  path =  "........../GSE114725/PLOTS",
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"   
)



