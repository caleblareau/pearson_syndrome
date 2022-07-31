library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(BuenColors)
library(dplyr)
library(viridis)
library(patchwork)

pearson_asap <- readRDS("../../../pearson_large_data_files/output/bone_marrow/pearson_asap_master_object.rds")
DimPlot(pearson_asap, label = TRUE)

DefaultAssay(pearson_asap) <- "ADT"
cd4adt <- FindMarkers(pearson_asap, "5", "0", logfc.threshold = 0)
cd8adt <- FindMarkers(pearson_asap, "6", "2", logfc.threshold = 0)
mdf <- merge(cd4adt, cd8adt, by = "row.names")

mdf %>% filter(p_val_adj.x < 0.01 & p_val_adj.y < 0.01) %>% 
  filter(avg_log2FC.x > 0 & avg_log2FC.y > 0) %>%
  arrange(p_val_adj.x)

p1 <- ggplot(mdf, aes(x = avg_log2FC.x, y = avg_log2FC.y, label = Row.names)) +
  geom_text(size = 2) + pretty_plot(fontsize = 7) + L_border() +
  labs(x = "log2FC CD4 RTE / Other Naive", y= "log2FC CD8 RTE / Other Naive")
#cowplot::ggsave2(p1, file = "../plots/RTE_ADT_scatter.pdf", width = 1.8, height = 1.8)

DefaultAssay(pearson_asap) <- "ACTIVITY"
cd4ga <- FindMarkers(pearson_asap, "5", "0", logfc.threshold = 0)
cd8ga <- FindMarkers(pearson_asap, "6", "2", logfc.threshold = 0)

data.frame(
  pearson_asap@reductions$umap@cell.embeddings,
  t(pearson_asap@assays$ACTIVITY@data[c("TOX", "IKZF2", "ZNF462", "ADAM23"),])
) %>% dplyr::filter(UMAP_2 < 5 & UMAP_1 > 0) -> tcell_plot_df

pA <- ggplot(shuf(tcell_plot_df), aes(x = UMAP_1, y = UMAP_2, color = ZNF462)) +
  geom_point(size = 0.7) +  scale_color_gradientn(colours = c("lightgrey",  "firebrick"),  limits = c(0,1.5), oob = scales::squish) +
  theme_void() + theme(legend.position = "none")

pB <- ggplot(shuf(tcell_plot_df), aes(x = UMAP_1, y = UMAP_2, color = ADAM23)) +
  geom_point(size = 0.7) +  scale_color_gradientn(colours = c("lightgrey",  "firebrick"),  limits = c(0,1.5), oob = scales::squish) +
  theme_void() + theme(legend.position = "none")

pC <- ggplot(shuf(tcell_plot_df), aes(x = UMAP_1, y = UMAP_2, color = TOX)) +
  geom_point(size = 0.7) +  scale_color_gradientn(colours = c("lightgrey",  "firebrick"),  limits = c(0,1.5), oob = scales::squish) +
  theme_void() + theme(legend.position = "none")

pD <- ggplot(shuf(tcell_plot_df), aes(x = UMAP_1, y = UMAP_2, color = IKZF2)) +
  geom_point(size = 0.7) +  scale_color_gradientn(colours = c("lightgrey",  "firebrick"),  limits = c(0,1.5), oob = scales::squish) +
  theme_void() + theme(legend.position = "none")

#cowplot::ggsave2(pA, file = "../plots/GA_ZNF462.png", width = 4, height = 4, dpi = 300)
#cowplot::ggsave2(pB, file = "../plots/GA_ADAM23.png", width = 4, height = 4, dpi = 300)
#cowplot::ggsave2(pC, file = "../plots/GA_TOX.png", width = 4, height = 4, dpi = 300)
#cowplot::ggsave2(pD, file = "../plots/GA_IKZF2.png", width = 4, height = 4, dpi = 300)


rte_markers <- c("AOAH", "ADA", "CACHD1", "DACH1", "FCGRT","ITGA4", "ITGA6", "TCF4", "TLR1", "LRRN3", "TOX", "CR2")            
key <- c("TOX", "IKZF2", "ZNF462", "ADAM23")
cd8ga[rownames(cd8ga) %in% key,]
cd4ga[rownames(cd4ga) %in% rte_markers,]

cd8ga$p_val_adj <- cd8ga$p_val_adj + 1e-190
cd4ga$p_val_adj <- cd4ga$p_val_adj + 1e-190
cd8ga$gene <- rownames(cd8ga)
cd4ga$gene <- rownames(cd4ga)

#write.table(cd4ga, file = "../output/cd4ga_RTEs.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#write.table(cd8ga, file = "../output/cd8ga_RTEs.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cd4ga$color = cd4ga$gene %in% c("TOX", "IKZF2", "ZNF462", "ADAM23")
cd8ga$color = cd8ga$gene %in% c("TOX", "IKZF2", "ZNF462", "ADAM23")

pcd4 <- ggplot(cd4ga %>% arrange((color)), aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
  geom_point(size = 0.5) + scale_color_manual(values = c("lightgrey", "firebrick")) +
  labs(x = "log2 FC", y = "-log10 adj pvalue")+ pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none")
pcd8 <-ggplot(cd8ga %>% arrange((color)), aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
  geom_point(size = 0.5) + scale_color_manual(values = c("lightgrey", "firebrick")) +
  labs(x = "log2 FC", y = "-log10 adj pvalue") + pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none")

#cowplot::ggsave2(pcd4, file = "../plots/GA_CD4.pdf", width = 1.8, height = 1.8)
#cowplot::ggsave2(pcd8, file = "../plots/GA_CD8.pdf", width = 1.8, height = 1.8)

cd4ga %>% dplyr::filter(color)
