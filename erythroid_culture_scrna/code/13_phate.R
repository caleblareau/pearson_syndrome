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

phate_harmony <- phate(so@reductions$harmony@cell.embeddings)

phdf <- data.frame(
  phate_harmony$embedding,
  so@meta.data
)
pa <- ggplot(shuf(phdf), aes(x = PHATE1, y = PHATE2, color = Day)) + 
  geom_point(size = 0.2)

pb <- ggplot(shuf(phdf), aes(x = PHATE1, y = PHATE2, color = assignment)) + 
  geom_point(size = 0.2)

pc <- ggplot(shuf(phdf), aes(x = PHATE1, y = PHATE2, color = MDS)) + 
  geom_point(size = 0.2)

library(patchwork)

pa + pb + pc
