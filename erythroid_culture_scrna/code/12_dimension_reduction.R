library(Seurat)
library(harmony)
library(magrittr)

so <- readRDS("../big_data/init_seurat_object.rds")
so <- NormalizeData(so)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
so <- ScaleData(so)
so <- RunPCA(so, features = VariableFeatures(object = so))
so <- RunHarmony(
  object = so,
  group.by.vars = c("Day", 'assignment'),
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

library(patchwork)
p1 + p2 + p3

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
                          group.by = "assignment", subset.ident = "4" )
head(de.markers, 40) %>% data.frame()
