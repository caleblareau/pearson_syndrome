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

p1 <- DimPlot(so, group.by = "Day", shuffle = TRUE)
p2 <- DimPlot(so, group.by = "assignment", shuffle = TRUE)
p3 <- DimPlot(so, label = TRUE)
p4 <- DimPlot(so, group.by = "MDS")
p5 <- FeaturePlot(so, "HBB")
p6 <- FeaturePlot(so, "module_score")

library(patchwork)
(p1 + p2 + p3 ) / (p4 + p5 + p6)
sc <- so@meta.data$seurat_clusters
vec <- case_when(
  sc %in% c("2", "7") ~ "Early/NC",
  sc %in% c("1") ~ "Late",
  sc %in% c("3") ~ "Trans",
  sc %in% c("8") ~ "Tcell",
  sc %in% c("9", "6", "0") ~ "Early/Cycling",
  sc %in% c("4", "5") ~ "Monocytic",
  TRUE ~ "other"
)

so@meta.data$anno <- vec

de.markers <- FindMarkers(so, ident.1 = "Pearson", ident.2 = "Healthy",
                          group.by = "assignment",  subset.ident = "2")

de.markers_4 <- FindMarkers(so, ident.1 = "4", "0")

head(de.markers, 40) %>% data.frame()


ggplot(so@meta.data, aes(x = MDS, y = module_score)) +
  geom_boxplot()

FeaturePlot(so, features = c("module_score"))

so@meta.data %>% group_by(assignment, Day, MDS) %>% 
  summarise(count = n()) %>% mutate(prop = count / sum(count))

up2 <- FindMarkers(so, ident.1 = "2", ident.2 = "1", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
up1 <- FindMarkers(so, ident.1 = "1",ident.2 =  "2", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
top <- markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
DotPlot(so, features = rev(unique(top$gene)), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "assignment") + RotatedAxis()

FeaturePlot(so, features =c("PSAT1", "CBS", "CTH", "PHGDH", "PSPH", "CTSA", "SCPEP1"), split.by = c("assignment"))

FeaturePlot(so, features = c("ACKR1", "GLRX", "RHCE", "HBA1", "FAM83A", "HBM"), split.by = c("assignment"))
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
  scale_color_manual(values = jdb_palette("corona")[c(1,3,2,4,6,5)]) + theme_void()+
  ggtitle("") + theme(legend.position = "none") 
cowplot::ggsave2(px, file = "../plots/annotation_erythroid.png", width = 5, height = 5.1, dpi = 500)
