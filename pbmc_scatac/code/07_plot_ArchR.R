library(Seurat)
library(viridis)

# Import saved data
so <- readRDS("../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")
dat <- readRDS("../../../pearson_mtscatac_large_data_files/output/24May2020_ArchR_smoothed_activity_scores.rds")
sum(colnames(so) == colnames(dat))

# Make simple plotting data frame
df <- data.frame(
  so@reductions$umap@cell.embeddings,
  cell = colnames(so)
)

df_plot <- data.frame(
  df, 
  cluster = so@meta.data$predicted.id,
  t(dat[, as.character(df$cell)])
)

df_plot %>% dplyr::filter(cluster %in% c("NK-like T-cells", "Effector CD8+ T-cells")) %>% group_by(cluster) %>%
  summarize_all(mean) -> mean_all
t(mean_all[, c(-1, -2, -3, -4)]) %>%  data.frame() %>% mutate(diff = X1 -X2, gene = rownames(dat)) %>%
  arrange(diff) %>% head()

fpcl <- function(gene){
  plot_df_one <- df_plot[,c("UMAP_1", "UMAP_2", gene)]
  colnames(plot_df_one) <- c("UMAP_1", "UMAP_2", "gene")
  ggplot(plot_df_one %>% arrange(gene), aes(x = UMAP_1, y = UMAP_2, color = gene)) + 
    geom_point(size = 0.1) + scale_color_gradientn(colors = jdb_palette("solar_basic")) +
    theme_void() + theme(legend.position = "none")
}

ggsave(cowplot::plot_grid(fpcl("CD3D"),  fpcl("CD4"), fpcl("CD8B"), fpcl("NCAM1"), fpcl("CCL5"), ncol = 5), 
       file = "../output/Tcell1.png", width = 14, height = 3, dpi = 500)

ggsave(cowplot::plot_grid(fpcl("GZMK"),  fpcl("KLRK1"), fpcl("EOMES"), fpcl("ZEB2"), fpcl("IKZF2"), ncol = 5), 
       file = "../output/Tcell2.png", width = 14, height = 3, dpi = 500)

ggsave(cowplot::plot_grid(fpcl("MS4A1"),  fpcl("CD14"), fpcl("FCER1A"), fpcl("FCGR3A"), fpcl("TLR7"),  ncol = 5), 
       file = "../output/other_cellstypes.png", width = 14, height = 3, dpi = 500)
