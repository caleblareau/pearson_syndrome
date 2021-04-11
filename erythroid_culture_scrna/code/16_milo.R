library(Seurat)
library(harmony)
library(magrittr)
library(miloR)
library(SingleCellExperiment)
library(scater)

# Preprocess
so <- readRDS("../../../pearson_large_data_files/output/invitro_erythroid/Invitro_erythroid_scRNA_seurat_object.rds")
so <- subset(so, module_score > 0.25)
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so)
sce <- as.SingleCellExperiment(so)

# Run milo
milo <- Milo(sce)
milo <- buildGraph(milo, k = 10, d = 30)
milo <- makeNhoods(milo, prop = 0.1, k = 10, d=30, refined = TRUE)

# Meta data weirdness
md <-  data.frame(colData(milo))
md$Day_assign <- paste0(md$Day, "_", md$assignment)
md$barcode <- rownames(md)
milo <- countCells(milo, meta.data =md, sample="Day_assign")
milo <- buildNhoodGraph(milo)
milo_design <- data.frame(md[,c("Day_assign", "assignment", "Day")])
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$Day_assign
da_results <- testNhoods(milo, design = ~ assignment, design.df = milo_design)

plotNhoodGraphDA(traj_milo, da_results, alpha=0.05)
