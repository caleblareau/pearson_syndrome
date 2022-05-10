library(Seurat)
library(BuenColors)
library(dplyr)
library(viridis)
library(data.table)

so <- readRDS("../../../pearson_large_data_files/output/Tcell_scATAC_culture.rds")
DefaultAssay(so) <- "ADT"
so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR", margin = 2)
FeatureScatter(so, feature1 = "rna_CD4", feature2 = "rna_CD8A", slot = "counts") 
pbmc2 <- ScaleData(pbmc2, assay="ADT")
FeaturePlot(object = so, features = c("CD3D", "CD4", "CD8A", "CD55", "CD27", "CD45RA", "ITGA2", "ITGA4", "CD63"),
            sort.cell = TRUE, max.cutoff = "q95") & scale_color_viridis()
