library(Seurat)
library(harmony)
library(magrittr)

so <- readRDS("../../../pearson_large_data_files/output/invitro_erythroid/Invitro_erythroid_scRNA_seurat_object.rds")
so <- subset(so, module_score > 0.25)
so <- NormalizeData(so)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
so <- ScaleData(so)

# Now do cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
so <- CellCycleScoring(so, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Run Harmony
so <- RunPCA(so, features = VariableFeatures(object = so))
so <- RunHarmony(
  object = so,
  group.by.vars = c("assignment", "Day", "Phase"),
  reduction = 'pca',
  assay.use = 'RNA',
  project.dim = FALSE
)
so <- RunUMAP(so, dims = 1:30, reduction = 'harmony')
so <- FindNeighbors(object = so, reduction = 'harmony', dims = 1:30)
so <- FindClusters(object = so, verbose = FALSE, resolution = 0.2)
DimPlot(so, label = TRUE)
FeaturePlot(so, "ALG2")

p1 <- DimPlot(so, group.by = "Day", shuffle = TRUE)
p2 <- DimPlot(so, group.by = "assignment", shuffle = TRUE)
p3 <- DimPlot(so, label = TRUE)
p4 <- DimPlot(so, group.by = "MDS")
p5 <- FeaturePlot(so, "HBA1")
p6 <- FeaturePlot(so, "module_score")

library(patchwork)
(p1 + p2 + p3 ) / (p4 + p5 + p6)
sc <- so@meta.data$seurat_clusters
vec <- case_when(
  sc %in% c("0", "5", "7", "9") ~ "Early/cycling",
  sc %in% c("1", "6") ~ "Late",
  sc %in% c("8") ~ "Tcell",
  sc %in% c("2") ~ "Late/NC",
  sc %in% c("4", "3") ~ "Monocytic",
  TRUE ~ "other"
)

so@meta.data$anno <- vec

de.markers <- FindMarkers(so, ident.1 = "Pearson", ident.2 = "Healthy",
                          group.by = "assignment",  subset.ident = "1")

head(de.markers, 40) %>% data.frame()


ggplot(so@meta.data, aes(x = MDS, y = module_score)) +
  geom_boxplot()

FeaturePlot(so, features = c("module_score"))

so@meta.data %>% group_by(assignment, Day, MDS) %>% 
  summarise(count = n()) %>% mutate(prop = count / sum(count))

up2 <- FindMarkers(so, ident.1 = "0", ident.2 = "1", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
up1 <- FindMarkers(so, ident.1 = "1",ident.2 =  "0", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

FeaturePlot(so, c("HBA1", "HBA2", "HBB"), min.cutoff = "q20", max.cutoff = "q95")



markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
top <- markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
DotPlot(so, features = rev(unique(top$gene)), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "assignment") + RotatedAxis()

FeaturePlot(so, features =c( "PPOX", "CPOX", "UROD", "UROS", "HMBS"), split.by = c("assignment"))


FeaturePlot(so, features = c( "HBA1", "FAM83A", "HBM"), split.by = c("assignment"))
FeaturePlot(so, features = c("PPBP", "CD8A", "CD14", "CD3D", "CD4"), split.by = "assignment")
de.markers <- FindMarkers(so, ident.1 = "Pearson", ident.2 = "Healthy",
                          group.by = "assignment" )

FeaturePlot(so, features = c("nCount_RNA", "nFeature_RNA"))



DimPlot(so,  shuffle = TRUE, label = TRUE) &
  scale_color_manual(values = jdb_palette("corona"))

p1 <- DimPlot(so, group.by = "Day", shuffle = TRUE) +
  theme_void() + scale_color_manual(values = c(jdb_palette("brewer_spectra")[c(1,9)])) +
  ggtitle("") + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/day_erythroid.png", width = 5, height = 5.1, dpi = 500)

p2 <- DimPlot(so, group.by = "assignment", shuffle = TRUE) +
  theme_void() + scale_color_manual(values = c("dodgerblue2", "red3")) +
  ggtitle("") + theme(legend.position = "none")
cowplot::ggsave2(p2, file = "../plots/donor_erythroid.png", width = 5, height = 5.1, dpi = 500)

p3 <- DimPlot(so, group.by = "Phase", shuffle = TRUE) +
  theme_void() + scale_color_manual(values = c("dodgerblue3", "firebrick", "green4")) +
  ggtitle("") + theme(legend.position = "none")
cowplot::ggsave2(p3, file = "../plots/phase_erythroid.png", width = 5, height = 5.1, dpi = 500)


px <- DimPlot(so, group.by = "anno", shuffle = TRUE, label = FALSE) +
  scale_color_manual(values = jdb_palette("corona")[c(1,2,3,4,6)]) + theme_void()+
  ggtitle("") + theme(legend.position = "none") 
cowplot::ggsave2(px, file = "../plots/annotation_erythroid.png", width = 5, height = 5.1, dpi = 500)


pg1 <- FeaturePlot(so, features =c("PHGDH"), split.by = c("assignment"), max.cutoff = "q99", min.cutoff = "q01") &
  theme_void() & theme(legend.position = "none")
cowplot::ggsave2(pg1, file = "../plots/gene_split_phghd.png", width = 6, height = 3, dpi = 500)

pg2 <- FeaturePlot(so, features =c("CPOX"), split.by = c("assignment"), max.cutoff = "q99", min.cutoff = "q01") &
  theme_void() & theme(legend.position = "none")
cowplot::ggsave2(pg2, file = "../plots/gene_split_cpox.png", width = 6, height = 3, dpi = 500)

pg3 <- FeaturePlot(so, features =c("HEBP2"), split.by = c("assignment"), max.cutoff = "q99", min.cutoff = "q01") &
  theme_void() & theme(legend.position = "none")
cowplot::ggsave2(pg3, file = "../plots/gene_split_hebp2.png", width = 6, height = 3, dpi = 500)
saveRDS(so@reductions$umap, "../output/umap_coords_scRNA_ery.rds")
