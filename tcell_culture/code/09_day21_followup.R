library(Seurat)
library(Signac)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(GenomeInfoDb)
library(harmony)
library(BuenColors)
library(ggbeeswarm)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(viridis)
library(Matrix)

pbmc2 <- readRDS( "../../../pearson_large_data_files/output/invitro_tcell/Tcell_scATAC_culture-day21.rds")
FindMarkers(pbmc2, "6", "4", only.pos = TRUE)


mat_var <- cbind(readRDS("../output/day21_called_variants.rds"))
pbmc2$af12631T.C <- mat_var["12631T>C",][colnames(pbmc2)]; pbmc2$af12631T.C  <- ifelse(is.na(pbmc2$af12631T.C), 0, pbmc2$af12631T.C )
pbmc2$af11838T.A <- mat_var["11838T>A",][colnames(pbmc2)]; pbmc2$af11838T.A  <- ifelse(is.na(pbmc2$af11838T.A), 0, pbmc2$af11838T.A )
pbmc2$af14476G.A <- mat_var["14476G>A",][colnames(pbmc2)]; pbmc2$af14476G.A  <- ifelse(is.na(pbmc2$af14476G.A), 0, pbmc2$af14476G.A )
pbmc2$af4225A.G <- mat_var["4225A>G",][colnames(pbmc2)]; pbmc2$af4225A.G  <- ifelse(is.na(pbmc2$af4225A.G), 0, pbmc2$af4225A.G )

p1af <- FeaturePlot(object = pbmc2, features = c("af4225A.G"), #,"af4225A.G" af12631T.C
                    sort.cell = TRUE, max.cutoff = "q97", reduction = "umap") + scale_color_viridis() +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1af, file = "../plots/af4225A.G.day21.png", width = 3, height = 3, dpi = 400)
p2af <- FeaturePlot(object = pbmc2, features = c("af12631T.C"), #,"af4225A.G" af12631T.C
                    sort.cell = TRUE, max.cutoff = "q97", reduction = "umap") + scale_color_viridis() +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p2af, file = "../plots/af12631T.C.day21.png", width = 3, height = 3, dpi = 400)

p3af <- FeaturePlot(object = pbmc2, features = c("heteroplasmy"), #,"af4225A.G" af12631T.C
                    sort.cell = TRUE, max.cutoff = "q97", reduction = "umap") + scale_color_viridis() +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p3af, file = "../plots/del.day21.png", width = 3, height = 3, dpi = 400)


DimPlot(pbmc2, label = TRUE)
pbmc2$caleb_cluster <- case_when(
  pbmc2$seurat_clusters %in% c(4,6) ~ "CD8 naive -like",
  pbmc2$seurat_clusters %in% c(1) ~ "CD4 naive-like",
  pbmc2$seurat_clusters %in% c(0) ~ "CD8 TEM-like",
  pbmc2$seurat_clusters %in% c(2) ~ "CD4 mem-like" ,
  pbmc2$seurat_clusters %in% c(3) ~ "CD4mem-clone1",
  pbmc2$seurat_clusters %in% c(5) ~ "CD4mem-clone2"
)

p1 <- DimPlot(object = pbmc2, group.by = "caleb_cluster",
              reduction = "umap") +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1, file = "../plots/clabels-d21.png", width = 3, height = 3, dpi = 400)

p2 <- ggplot(pbmc2@meta.data, aes(x = heteroplasmy, color = caleb_cluster)) +
  stat_ecdf() +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))
cowplot::ggsave2(p2, file = paste0("../plots/ecdf_tcell_clusters_day21.pdf"), width = 2.9, height = 1.5)

