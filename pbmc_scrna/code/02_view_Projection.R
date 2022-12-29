library(dplyr)
library(data.table)
library(BuenColors)

libs <- gsub(".rds", "", gsub("seurat_projected_meta_", "", list.files("../output/", pattern = "^seurat_projected")))

lapply(libs, function(dir_base){
  readRDS(paste0("../output/seurat_projected_meta_", dir_base, ".rds")) %>%
    mutate(id = dir_base)
})%>% rbindlist() %>% data.frame() -> rdf

rdf$split <- ifelse(rdf$pheno == "Pearson", rdf$name, ifelse(rdf$name %in% c("H1", "H2"), "yadult", "zpediatric"))
table(rdf$split)
p1 <- ggplot(rdf, aes(x = refUMAP_1, y = refUMAP_2, color = predicted.celltype.l2)) +
  facet_wrap(~split, nrow = 1)+
  geom_point(size = 0.2) + scale_color_manual(values = (c(jdb_palette("corona"), "salmon"))) +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/five_umap_clusters.png", width = 12, height = 3.4, dpi = 500)

###############
p1 <- ggplot(rdf, aes(x = refUMAP_1, y = refUMAP_2, color = predicted.celltype.l2)) +
  facet_wrap(~split, nrow = 1)+
  geom_point(size = 0.2) + scale_color_manual(values = (c(jdb_palette("corona"), "salmon"))) +
  theme_void() 
p2 <- g_legend(p1)
cowplot::ggsave2(p2, file = "../plots/legend_azimut.pdf", width = 4, height = 7, dpi = 500)
