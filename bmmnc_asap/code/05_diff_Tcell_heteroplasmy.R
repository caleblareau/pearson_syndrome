library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(BuenColors)
library(dplyr)
library(viridis)
library(patchwork)

pearson_asap <- readRDS("../../../pearson_large_data_files/output/bone_marrow/pearson_asap_master_object.rds")
DimPlot(pearson_asap, label = TRUE)
pearson_asap$UMAP1 <- pearson_asap@reductions$umap@cell.embeddings[,1]
pearson_asap$UMAP2 <- pearson_asap@reductions$umap@cell.embeddings[,2]

pearson_asap_ss <- subset(pearson_asap, seurat_clusters %in% c("8","9") & UMAP2 > 4 & UMAP1 > -2)
DimPlot(pearson_asap_ss, label = TRUE)

FeaturePlot(object = pearson_asap_ss, "heteroplasmy") +
  scale_color_viridis() + theme_void() + scale_

DefaultAssay(pearson_asap) <- "ADT"
DefaultAssay(pearson_asap_ss) <- "ADT"
FeaturePlot(pearson_asap_ss, c("CD3-2", "CD8","CD16", "CD195"), max.cutoff = "q90", min.cutoff = "q01", sort = TRUE) & 
  scale_color_gradientn(colors = jdb_palette("solar_extra"))


plot_factor_qs_ss <- function(i){
  p1 <- FeaturePlot(object = pearson_asap_ss, i, pt.size = 0.3, sort = TRUE,
                    max.cutoff = "q90") + theme_void() + theme(legend.position = "none") +
    scale_color_gradientn(colors = jdb_palette("solar_extra"))+ ggtitle("")
  return(p1)
}
cowplot::ggsave2(
  cowplot::plot_grid(
    plot_factor_qs_ss("CD3-2"),
    plot_factor_qs_ss("CD8"),
    plot_factor_qs_ss("CD16"),
    plot_factor_qs_ss("CD195"), ncol = 2, scale = 0.9, byrow = TRUE
  ), file = "../plots/adt_small_panels_Tcells_supp_q90.png", width = 1.5*4, height = 1.5*4, dpi = 300)


# make the volcano plot
pval_vec <- sapply(1:dim(pearson_asap_ss@assays$ADT@data), function(i){
  cor.test(pearson_asap_ss$heteroplasmy, (pearson_asap_ss@assays$ADT@data[i,]))$p.value
})
ddf <- data.frame(
  marker = rownames(pearson_asap_ss@assays$ADT@data), 
  correlation = (cor(pearson_asap_ss$heteroplasmy, t(pearson_asap_ss@assays$ADT@data))[1,]),
  pvalue = pval_vec
)
ddf$adjP <- p.adjust(ddf$pvalue)
p1 <- ggplot(ddf, aes(x = correlation, y = -log10(adjP))) +
  geom_point(size = 0.5) +
  pretty_plot(fontsize = 7) + L_border() +  labs(x = "Heteroplasmy Correlation", y = "-log10 adjusted pvalue")

cowplot::ggsave2(cowplot::plot_grid(p1), filename = "../plots/heteroplasmy_markers.pdf", 
                 width = 1.7, height = 1.7)

pA <- FeaturePlot(object = pearson_asap_ss, "heteroplasmy", pt.size = 0.2) +
  scale_color_viridis(begin = 0, end = max(pearson_asap_ss$heteroplasmy)/100 ) + theme_void() + ggtitle("") +
  theme(legend.position = "none")
cowplot::ggsave2(pA, filename = "../plots/zoom_heteroplasmy_nkt.pdf", 
                 width = 1.7, height = 1.7)
