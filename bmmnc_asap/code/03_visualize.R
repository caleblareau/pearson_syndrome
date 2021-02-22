library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(BuenColors)
library(dplyr)
library(viridis)
library(patchwork)

# import the data

pearson_asap <- readRDS( "../../../pearson_large_data_files/output/asap/pearson_asap_master_object.rds")
DimPlot(pearson_asap, label = TRUE)
DefaultAssay(pearson_asap) <- "ADT"

# First make volcano of markers for MDS or not
if(FALSE){
  
  annotate_markers <- function(markers){
    case_when(
      grepl("NCAM", markers) ~ "CD56",
      grepl("HLA-DR", markers) ~ "HLA-DR",
      grepl("CD64", markers) ~ "CD64",
      markers == "CD15" ~ "CD15",
      TRUE ~ "other"
    )
  }
  
  # Find top markers using Seurat for MDS
  Idents(pearson_asap) <- pearson_asap@meta.data$chr7
  pearson_asap_ss_myeloid <- subset(pearson_asap, seurat_clusters %in% c(1))
  mono <- FindMarkers(pearson_asap_ss_myeloid, ident.1 = "Monosomy7", ident.2 = "Wildtype", logfc.threshold = 0.01) 
  pearson_asap_ss_ery <- subset(pearson_asap, seurat_clusters %in% c(4))
  ery <- FindMarkers(pearson_asap_ss_ery, ident.1 = "Monosomy7", ident.2 = "Wildtype", logfc.threshold = 0.01) 
  
  # Annotate for plotting
  mono$color <- annotate_markers(rownames(mono))
  ery$color <- annotate_markers(rownames(ery))
  
  p1 <- ggplot(mono, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(size = 0.5) +
    scale_y_continuous(limits = c(0, 64)) + scale_x_continuous(limits = c(-0.4, 0.8)) +
    pretty_plot(fontsize = 6) + L_border() +  labs(x = "log2FC Monosomy 7 / WT", y = "-log10 adjusted pvalue", color = "") +
    scale_color_manual(values = c("blue2", "purple2", "red2", "purple4", "lightgrey"))
  
  p2 <- ggplot(ery, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(size = 0.5) +
    scale_y_continuous(limits = c(0, 64)) + scale_x_continuous(limits = c(-0.4, 0.8)) +
    pretty_plot(fontsize = 6) + L_border() +  labs(x = "log2FC Monosomy 7 / WT", y = "-log10 adjusted pvalue", color = "") +
    scale_color_manual(values = c("blue2", "purple2", "red2", "purple4", "lightgrey"))
  
  cowplot::ggsave2(cowplot::plot_grid(p1, p2, nrow = 2), filename = "../plots/volcanos_MDS.pdf", 
                   width = 2.4, height = 2.4)
}

lapply(rownames(pearson_asap@assays$ADT), function(i){
  p1 <- FeaturePlot(object = pearson_asap, i, pt.size = 0.1,
                    max.cutoff = "q90") + theme_void() + theme(legend.position = "none") +
    scale_color_gradientn(colors = jdb_palette("solar_extra")) 
  
  iname <- gsub("/", "_", i)
  cowplot::ggsave2(p1, file = paste0("../plots/viz_adts/adt_small_panels_supplement_q90_",iname,".png"), 
                   width = 4, height = 4.3, dpi = 400)
})


plot_factor_qs <- function(i){
  p1 <- FeaturePlot(object = pearson_asap, i, pt.size = 0.1,
                   max.cutoff = "q90") + theme_void() + theme(legend.position = "none") +
    scale_color_gradientn(colors = jdb_palette("solar_extra"))+ ggtitle("")
  return(p1)
}
cowplot::ggsave2(
  cowplot::plot_grid(
    plot_factor_qs("CD71"),
    plot_factor_qs("CD117(c-kit)"),
    plot_factor_qs("CD33"),
    plot_factor_qs("CD35"),
    plot_factor_qs("CD4-1"),
    plot_factor_qs("CD8a"), ncol = 3, scale = 0.8, byrow = FALSE
  ), file = "../plots/adt_small_panels_main_q90.png", width = 2.25*4, height = 1.5*4, dpi = 300)

cowplot::ggsave2(
  cowplot::plot_grid(
    plot_factor_qs("CD45-2"),
    plot_factor_qs("CD19"),
    plot_factor_qs("CD123"),
    plot_factor_qs("CD16"),
    plot_factor_qs("CD56(NCAM)Recombinant"),
    plot_factor_qs("CD335"), ncol = 3, scale = 0.8, byrow = FALSE
  ), file = "../plots/adt_small_panels_supplement_q90.png", width = 2.25*4, height = 1.5*4, dpi = 300)



FeaturePlot(object = pearson_asap, c("CD38-2", "CD71", "CD172a", "CD34", "CD117(c-kit)", "CD59", "CLEC12A"),
            max.cutoff = "q90", cols =jdb_palette("brewer_spectra"))

new_clust = c("CD4T1", "Mono1", "CD8T1", "CD4T2", "Ery1", "CD4T3", "CD8T2", "Mono2", "NKcell", "CD8T3", "Baso", "Ery2", "Bcell1", "Ery3", "Bcell2", "Mono3", "CD4T4", "plasma")
names(new_clust) = as.character(0:17)

# Set up a color scheme
ery <- c("#FB6A4A", "red2", "#A50F15"); names(ery) <- c("Ery1", "Ery2", "Ery3")
cd4t <- jdb_palette("brewer_marine")[c(8,7,6,5)]; names(cd4t) <- c("CD4T1", "CD4T2", "CD4T3", "CD4T4")
short <- c("#CD8500",  "#094078"); names(short) <- c("Baso",  "NKcell")
cd8t <- c("#9DA8D7", "#6D49AB","#581F99"); names(cd8t) <- c("CD8T1", "CD8T2", "CD8T3")
bcell <- c("#50C878", "forestgreen","#28AB87"); names(bcell) <- c("Bcell1", "Bcell2", "plasma")
mono <- jdb_palette("brewer_fire")[c(4,6,8)]; names(mono) <- c("Mono1", "Mono2", "Mono3")

full_color_vec <- c(mono, bcell, cd8t, short, cd4t, ery)

base_plot_df <- data.frame(
  pearson_asap@reductions$umap@cell.embeddings,
  pearson_asap@meta.data
)
base_plot_df$other_anno <- new_clust[as.character(base_plot_df$seurat_clusters)]

p_base <- ggplot(base_plot_df, aes(x = UMAP_1, y = UMAP_2, color = other_anno)) + 
  geom_point(size = 0.3) + scale_color_manual(values = full_color_vec) +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(p_base, file = "../plots/base_umap_colors.png", width = 8, height = 8, dpi = 400)

p_mds <- ggplot(shuf(base_plot_df), aes(x = UMAP_1, y = UMAP_2, color = chr7)) + 
  geom_point(size = 0.3) + scale_color_manual(values = c("red", "dodgerblue3")) +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(p_mds, file = "../plots/mds_embedding.png", width = 8, height = 8, dpi = 400)

p_het <- ggplot(base_plot_df, aes(x = UMAP_1, y = UMAP_2, color = heteroplasmy)) + 
  geom_point(size = 0.3) + 
  scale_color_viridis(na.value="lightgrey", end = max(base_plot_df$heteroplasmy, na.rm = TRUE)/100) +
  theme_void() + theme(legend.position = "none") 
cowplot::ggsave2(p_het, file = "../plots/het_embedding.png", width = 8, height = 8, dpi = 400)
