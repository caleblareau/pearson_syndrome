library(BuenColors)

libs <- c("pearson_bci", "pearson_ccf", "pearson_mds", "healthy_pbmc_8k_v2-remap", "healthy_pbmc_4k_v2-remap", "healthy_pbmc_5k_nextgem", "healthy_pbmc_5k")

lapply(libs, function(dir_base){
  readRDS(paste0("../output/seurat_projected_meta_", dir_base, ".rds"))
})%>% rbindlist() %>% data.frame() -> rdf

rdf$split <- ifelse(rdf$pheno == "Pearson", rdf$name, "Healthy")
p1 <- ggplot(rdf, aes(x = refUMAP_1, y = refUMAP_2, color = predicted.celltype.l2)) +
  facet_wrap(~split, nrow = 1)+
  geom_point(size = 0.2) + scale_color_manual(values = (c(jdb_palette("corona"), "salmon"))) +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/four_umap_clusters.png", width = 12, height = 3.4, dpi = 500)

p1 <- ggplot(rdf, aes(x = refUMAP_1, y = refUMAP_2, color = predicted.celltype.l2)) +
  facet_wrap(~split, nrow = 1)+
  geom_point(size = 0.2) + scale_color_manual(values = (c(jdb_palette("corona"), "salmon"))) +
  theme_void() 
p2 <- g_legend(p1)
cowplot::ggsave2(p2, file = "../plots/legend_azimut.pdf", width = 4, height = 7, dpi = 500)
