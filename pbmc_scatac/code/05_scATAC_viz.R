library(Seurat)
pbmc <- readRDS("../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")
rna <- readRDS("../../../pearson_mtscatac_large_data_files/output/5March-PearsonRNAseq-integration.rds")

pbmc$refined <- gsub("Activated B-cells", "Naive B-cells", pbmc$predicted.id)
col_vec <- as.character(jdb_palette("corona"))[1:23]
names(col_vec) <- as.character(unique(rna@meta.data$celltype))

ppng <-  DimPlot(pbmc, group.by = "refined") +
  scale_color_manual(values = col_vec)+
  theme_void() + theme(legend.position = 'none') 
cowplot::ggsave2(ppng, file = "../output/scATAC_plot_rna_colors.png", width = 4, height = 4, dpi = 500)
