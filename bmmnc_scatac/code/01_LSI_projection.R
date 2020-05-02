library(Seurat)
library(Matrix)
library(SummarizedExperiment)
library(uwot)
source("01a_LSI_BMMNC_helper.R")

#Previous Reference Summarized Experiment
#Contains Peaks for Reference Hematopoiesis only
se <- readRDS("../../../pearson_mtscatac_large_data_files/input/bmmnc/scATAC-Healthy-Hematopoiesis-191120.rds")
seDisease <-  readRDS("../../../pearson_mtscatac_large_data_files/output/Pearson-BMMNC-PT3_SummarizedExperiment.rds")

#Load Saved UMAP Manifold
umapManifold <- uwot::load_uwot("../data/scATAC-Projection-UMAP/scATAC-Hematopoiesis-UMAP-model.190505.uwot.tar")

#LSI Projection Matrix
lsiPeaks <- metadata(se)$variablePeaks
matProjectLSI <- assay(seDisease[lsiPeaks,])

#LSI Project
lsiReference <- metadata(se)$LSI
lsiProjection <- projectLSI(matProjectLSI, lsiReference)

#UMAP Projection
#Set Seed Prior to umap_transform (see uwot github)
set.seed(1)
umapProjection <- uwot::umap_transform(data.matrix(lsiProjection), umapManifold, verbose = TRUE)

#Plot Projection
del_df <- fread("../data/del_BMMNC_PT3.deletion_heteroplasmy.tsv")[rep(c(FALSE, FALSE, TRUE), 6000), c("cell_id", "heteroplasmy")]
refDF <- data.frame(barcode = colnames(se), X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2], Type = "reference", 
                    heteroplasmy = NA)
proDF <- merge(data.frame(barcode = colnames(seDisease), X1 = umapProjection[,1], X2 = umapProjection[,2], Type = "Pearson"), 
               del_df, by.x = "barcode", by.y = "cell_id")
projectionDF <- rbind(refDF, proDF)


ggplot(projectionDF, aes(x = X1, y = X2, color = Type)) +
  geom_point() + pretty_plot() + L_border()+
  scale_color_manual(values = c("lightgrey", "firebrick"))

ggplot(projectionDF, aes(x = X1, y = X2, color = heteroplasmy)) +
  geom_point() + pretty_plot() + L_border()+
  scale_color_gradientn(colors = c("blue", "firebrick"))

ggplot(data.frame(colData(se)), aes(x = UMAP1, y = UMAP2, color= Clusters)) +
  geom_point() + scale_color_manual(values=metadata(se)$colorMap$Clusters) 
  