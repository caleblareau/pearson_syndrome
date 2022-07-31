library(Seurat)

pbmc2 <- readRDS("../../../pearson_large_data_files/output/invitro_tcell/Tcell_scATAC_culture.rds")


DefaultAssay(pbmc2) < "peaks"
FeaturePlot(object = pbmc2, features = "heteroplasmy" ) +
  scale_color_viridis()

set.seed(1)
cormat <- cor(pbmc2@meta.data$heteroplasmy,t(data.matrix(pbmc2@assays$ADT@data)),
              use = "pairwise.complete")
perm_cormat <- cor(sample(pbmc2@meta.data$heteroplasmy), t(data.matrix(pbmc2@assays$ADT@data)),
                   use = "pairwise.complete")

obs_df <- data.frame(
  gene = colnames(cormat),
  cor = cormat[1,]
) %>% dplyr::filter(!is.na(cor)) %>% arrange(desc(cor)) %>% mutate(rank = 1:n(), what = "observed") 
obs_df


# Pearson correlation to see factors most associated with activation / not
set.seed(1)
cormat <- cor(pbmc2@meta.data$heteroplasmy,t(data.matrix(pbmc2@assays$RNA@data)),
              use = "pairwise.complete")
perm_cormat <- cor(sample(pbmc2@meta.data$heteroplasmy), t(data.matrix(pbmc2@assays$RNA@data)),
                   use = "pairwise.complete")

obs_df <- data.frame(
  gene = colnames(cormat),
  cor = cormat[1,]
) %>% dplyr::filter(!is.na(cor)) %>% arrange(desc(cor)) %>% mutate(rank = 1:n(), what = "observed") 
obs_df

write.table(obs_df, file = "../output/heteroplasmy_Tcell_association_ranking.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

perm_df <- data.frame(
  gene = colnames(perm_cormat),
  cor = perm_cormat[1,]
) %>% dplyr::filter(!is.na(cor)) %>% arrange(desc(cor)) %>% mutate(rank = 1:n(), what = "permuted") 

p1 <- ggplot(obs_df, aes(x = rank, y = cor)) + 
  geom_point(size = 0.2) + scale_color_manual(values = c("black", "firebrick"))  + 
  geom_point(inherit.aes = FALSE, data = perm_df, aes(x = rank, y = cor), color = "lightgrey", size = 0.2) +
  labs(x = "Rank sorted genes", y = "Correlation") + 
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") + theme_void()
cowplot::ggsave2(p1, file = "../plots/Tcell_heteroplasmy_assoc.png", width = 3, height = 3, dpi = 500)

saveRDS(pbmc2, file = "../../../pearson_large_data_files/output/Tcell_scATAC_culture.rds")
