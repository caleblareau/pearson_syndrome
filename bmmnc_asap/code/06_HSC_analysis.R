library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(BuenColors)
library(dplyr)
library(viridis)
library(patchwork)

pearson_asap <- readRDS( "../../../pearson_large_data_files/output/asap/pearson_asap_master_object.rds")
DimPlot(pearson_asap, label = TRUE)
pearson_asap$UMAP1 <- pearson_asap@reductions$umap@cell.embeddings[,1]
pearson_asap$UMAP2 <- pearson_asap@reductions$umap@cell.embeddings[,2]

df <- data.frame(
  t(pearson_asap@assays$ADT@data)
)

DefaultAssay(pearson_asap) <- "ADT"
pearson_asap$U1 <- pearson_asap@reductions$umap@cell.embeddings[,1]
pearson_asap$U2 <-  pearson_asap@reductions$umap@cell.embeddings[,2]
pearson_asap_ss <- subset(pearson_asap, seurat_clusters %in% "4" & U2 > -3.5 & U1 < -4)

# Threshold matches full data q90
p1 <-FeaturePlot(pearson_asap_ss, c("CD117(c-kit)", "CD71",  "CD34"), max.cutoff = c(1.5,1.5, 1), min.cutoff = 0, sort = TRUE, pt.size = 0.2, ncol = 3) &
  scale_color_gradientn(colors = jdb_palette("solar_extra")) & theme_void() & theme(legend.position = "none") &
  ggtitle("") 
cowplot::ggsave2(p1, file = "../plots/erythroid_viz_prog.png", width = 6, height = 2, dpi = 300)


# Code chunk from ASAP repository
if(FALSE){
  adtbm$U1 <- umap_df[,1]
  adtbm$U2 <- umap_df[,2]
  
  adtbmss <- subset(adtbm, seurat_clusters %in% c("C1", "C2", "C8", "C20") & U1 < 10 & U2 < -2 & U2 > -4.5)
  DimPlot(adtbmss, label = TRUE)
  p1 <- FeaturePlot(adtbmss, c("CD117(c-kit)", "CD71", "CD34"), max.cutoff = "q90", min.cutoff = "q01", sort = TRUE, pt.size = 0.2, ncol = 3) &
    scale_color_gradientn(colors = jdb_palette("solar_extra")) & theme_void() & theme(legend.position = "none") &
    ggtitle("") 
  cowplot::ggsave2(p1, file = "../../../../../Desktop/ASAP_healthy.png", width = 6, height = 2, dpi = 300)
}
