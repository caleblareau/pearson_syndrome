library(Seurat)
library(harmony)
library(magrittr)

so <- readRDS("../../../pearson_mtscatac_large_data_files/output/invitro_erythroid/Invitro_erythroid_scRNA_seurat_object.rds")
so <- subset(so, module_score > 0.25)
so <- NormalizeData(so)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
so <- ScaleData(so)
so <- RunPCA(so, features = VariableFeatures(object = so))
so <- RunHarmony(
  object = so,
  group.by.vars = c("assignment", "Day"),
  reduction = 'pca',
  assay.use = 'RNA',
  project.dim = FALSE
)
so <- RunUMAP(so, dims = 1:20, reduction = 'harmony')
so <- FindNeighbors(object = so, reduction = 'harmony', dims = 1:20)
so <- FindClusters(object = so, verbose = FALSE, resolution = 0.2)

p1 <- DimPlot(so, group.by = "Day", shuffle = TRUE)
p2 <- DimPlot(so, group.by = "assignment", shuffle = TRUE)
p3 <- DimPlot(so, label = TRUE)
p4 <- DimPlot(so, group.by = "MDS")
p5 <- FeaturePlot(so, "HBB")
p6 <- FeaturePlot(so, "module_score")
library(patchwork)
(p1 + p2 + p3 ) / (p4 + p5 + p6)

de.markers <- FindMarkers(so, ident.1 = "Pearson", ident.2 = "Healthy",
                          group.by = "assignment",  subset.ident = "2")
head(de.markers, 40) %>% data.frame()
de.markers[c("ALAS2", "ABCB7", "SLC19A2", "GLRX5","PSU1", "GLRX"),]


ggplot(so@meta.data, aes(x = MDS, y = module_score)) +
  geom_boxplot()

FeaturePlot(so, features = c("module_score"))

so@meta.data %>% group_by(assignment, Day, MDS) %>% 
  summarise(count = n()) %>% mutate(prop = count / sum(count))

FindMarkers(so, ident.1 = "5", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
top <- markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
DotPlot(so, features = rev(unique(top$gene)), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "assignment") + RotatedAxis()

FeaturePlot(so, features = c("PHGDH", "HEBP2", "ALAS2", "GATA1", "TUBB2A", "TMCC2"), split.by = c("assignment"))
FeaturePlot(so, features = c("ACKR1", "GLRX", "RHCE", "HBA1", "FAM83A", "HBM"), split.by = c("assignment"))
FeaturePlot(so, features = c("PPBP", "CD8A", "CD14", "CD3D", "CD4"), split.by = "assignment")
de.markers <- FindMarkers(so, ident.1 = "Pearson", ident.2 = "Healthy",
                          group.by = "assignment" )

de.markers[c("PSAT1", "CBS", "CTH", "PHGDH", "PSPH", "CTSA", "SCPEP1"),]
