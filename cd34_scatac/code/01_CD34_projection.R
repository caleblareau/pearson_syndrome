library(Matrix)
library(SummarizedExperiment)
library(tidyverse)
library(umap)
library(edgeR)
library(FNN)
library(matrixStats)
library(igraph)
library(BuenColors)
library(GenomicRanges)
library(Seurat)
set.seed(1)
"%ni%" <- Negate("%in%")
source("01a_LSI_project_helpers.R")

nTop = 25000

# Use Granja data as base
SE <- readRDS("../../../pearson_mtscatac_large_data_files/output/granja_10X_CD34.rds")
SE <- SE[,colData(SE)$Group != "CD34_Progenitors_Rep1"]

#Run LSI 1st Iteration
lsi1 <- calcLSI(assay(SE), nComponents = 25, binarize = TRUE, nFeatures = NULL)
in_mat <- lsi1[[1]][,c(2:25)]
rownames(in_mat) <- SE@colData$Group_Barcode
clust1 <- seuratSNN(in_mat)

#Make Pseudo Bulk Library
message("Making PseudoBulk...")
clusterSums <- groupSums(mat = assay(SE), groups = clust1, sparse = TRUE) #Group Sums
logMat <- edgeR::cpm(clusterSums, log = TRUE, prior.count = 3) #log CPM matrix
varPeaks <- head(order(matrixStats::rowVars(logMat), decreasing = TRUE), nTop) #Top variable peaks

#Run LSI 2nd Iteration
lsi2 <- calcLSI(assay(SE)[varPeaks,,drop=FALSE], nComponents = 25, binarize = TRUE, nFeatures = NULL)
in_mat2 <- lsi2[[1]][,c(2:25)]
rownames(in_mat2) <- SE@colData$Group_Barcode
clust2 <- seuratSNN(in_mat2)
length(unique(clust2))

#UMAP
set.seed(1)
umap <- umap::umap(
  lsi2$matSVD[,2:25], 
  n_neighbors = 55, # original 55
  min_dist = 0.45, # original 0.45
  metric = "cosine", 
  verbose = TRUE    )
set.seed(10)

# Multiply by -1 to make the pseudotime read left to right
plot_df <- data.frame(umap$layout, colData(SE), Clusters = clust2, barcode = colData(SE)$Barcodes)

p0 <- ggplot(plot_df, aes(x= X1, y = X2, color = Clusters)) +
   geom_point(size = 0.5) +
   labs(x = "UMAP1", y= "UMAP2", color = "") +
  pretty_plot() + L_border() + theme(legend.position = "bottom")

sel <- readRDS("../data/granja_cd34/granja_published_C1.rds")
lsiProjection <- projectLSI((assay(sel)[varPeaks,]), lsi2)
umapProjection <- round(predict(umap, data.matrix(lsiProjection[,2:25])), 2)

projection_df_C1 <- data.frame(
  celltype = c(gsub("BM_", "", colData(sel)$Group), rep("none", dim(plot_df)[1])),
  umap1 = c(umapProjection[,1], plot_df$X1),
  umap2 = c(umapProjection[,2], plot_df$X2), 
  barcode = c(colData(sel)$Barcodes, plot_df$barcode)
)

p1 <- ggplot(projection_df_C1[dim(projection_df_C1)[1]:1,], aes(x= umap1, y = umap2, color = celltype, label = celltype)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y= "UMAP2", color = "C1 FACS ") +
  pretty_plot() + L_border() + theme(legend.position = "bottom") +
  scale_color_manual(values = c(ejc_color_maps, "none" = "lightgrey"))


sel <- readRDS("../../../pearson_mtscatac_large_data_files/output/Pearson-CD34-PT3_SummarizedExperiment.rds")
lsiProjection <- projectLSI((assay(sel)[varPeaks,]), lsi2)
umapProjection <- round(predict(umap, data.matrix(lsiProjection[,2:25])), 2)

projection_df_pearson_pt <- data.frame(
  celltype = c(rep("Pearson", dim(umapProjection)[1]), rep("base", dim(plot_df)[1])),
  umap1 = c(umapProjection[,1], plot_df$X1),
  umap2 = c(umapProjection[,2], plot_df$X2),
  barcode = c(colnames(sel), plot_df$barcode)
)

p2 <- ggplot(projection_df_pearson_pt[dim(projection_df_pearson_pt)[1]:1,], aes(x= umap1, y = umap2, color = celltype, label = celltype)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y= "UMAP2", color = "Donor") +
  pretty_plot() + L_border() + theme(legend.position = "bottom") +
  scale_color_manual(values = c("Pearson" = "firebrick", "base" = "lightgrey"))


sel <- readRDS("../../../pearson_mtscatac_large_data_files/output/CD34_G10_SummarizedExperiment.rds")
lsiProjection <- projectLSI((assay(sel)[varPeaks,]), lsi2)
umapProjection <- round(predict(umap, data.matrix(lsiProjection[,2:25])), 2)

projection_df_control <- data.frame(
  celltype = c(rep("Healthy", dim(umapProjection)[1]), rep("base", dim(plot_df)[1])),
  umap1 = c(umapProjection[,1], plot_df$X1),
  umap2 = c(umapProjection[,2], plot_df$X2),
  barcode = c(colnames(sel), plot_df$barcode)
)

p3 <- ggplot(projection_df_control[dim(projection_df_control)[1]:1,], aes(x= umap1, y = umap2, color = celltype, label = celltype)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y= "UMAP2", color = "Donor") +
  pretty_plot() + L_border() + theme(legend.position = "bottom") +
  scale_color_manual(values = c("Healthy" = "dodgerblue3", "base" = "lightgrey"))

cowplot:::ggsave2(cowplot::plot_grid(p0, p1, p2, p3, nrow = 2),
                  filename = "../plots/projection_CD34_viz.png", width = 6, height = 6)

save(projection_df_C1, projection_df_pearson_pt, projection_df_control, file = "../output/CD34_umap_embedding_granja_proj3.rda")

projection_df_pearson_pt %>% filter(celltype == "Pearson" & umap1 < -7.5)
projection_df_control %>% filter(celltype == "Healthy" & umap1 < -7.5)
