library(Seurat)
library(SeuratDisk)

reference <- LoadH5Seurat("../../../pearson_large_data_files/input/pbmc/pbmc_multimodal.h5seurat")
p1 <-DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = FALSE, label.size = 3, repel = TRUE, raster = FALSE) +
  NoLegend() + theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1, file = "../plots/bridge_embedding.png", width = 8, height = 8, dpi = 400)
