library(Seurat)
library(BuenColors)
library(dplyr)
library(viridis)
library(data.table)
library(Signac)

so <- readRDS("../../../pearson_large_data_files/output/invitro_tcell/Tcell_scATAC_culture.rds")


DefaultAssay(so) <- "ADT"
so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR", margin = 2)
FeatureScatter(so, feature1 = "adt_CD4", feature2 = "rna_CD8A", slot = "data") 

FeatureScatter(so, feature1 = "rna_CD4", feature2 = "rna_CD8A", slot = "data") 
pbmc2 <- ScaleData(pbmc2, assay="ADT")
FeaturePlot(object = so, features = c("CD3D", "CD4", "CD8A", "CD55", "CD27", "CD45RA", "ITGA2", "ITGA4", "CD63", "CD45RO"),
            sort.cell = TRUE, max.cutoff = "q95") & scale_color_viridis()
FeaturePlot(object = so, features = c("heteroplasmy"),
            sort.cell = TRUE, max.cutoff = "q95") & scale_color_viridis()
mat_var <- cbind(readRDS("../output/day14v1_called_variants.rds"), 
      readRDS("../output/day14v2_called_variants.rds"))
so$af12631T.C <- mat_var["12631T>C",][colnames(so)]; so$af12631T.C  <- ifelse(is.na(so$af12631T.C), 0, so$af12631T.C )
so$af4225A.G <- mat_var["4225A>G",][colnames(so)]; so$af4225A.G  <- ifelse(is.na(so$af4225A.G), 0, so$af4225A.G )
so$af11199C.T <- mat_var["11199C>T",][colnames(so)]; so$af11199C.T  <- ifelse(is.na(so$af11199C.T), 0, so$af11199C.T )
so$af15244A.G <- mat_var["15244A>G",][colnames(so)]; so$af15244A.G  <- ifelse(is.na(so$af15244A.G), 0, so$af15244A.G )

FeaturePlot(object = so, features = c("af12631T.C","af4225A.G","af11199C.T","af15244A.G"),
            sort.cell = TRUE, max.cutoff = "q97") & scale_color_viridis()
