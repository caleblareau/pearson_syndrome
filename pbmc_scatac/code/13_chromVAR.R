library(Seurat)
library(Signac)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BuenColors)
library(motifmatchr)
library(chromVARmotifs)
data("human_pwms_v2") #filtered collection of human motifs from cisBP database

source("../functions/import_counts_mat_10x.R")
set.seed(1)

# Import ATAC data
counts <- import_counts_mat_10x("../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_filtered_peak_bc_matrix/")
pbmc_pearson_md <- readRDS("../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")@meta.data
pbmc_healthy_md <- readRDS("../../../pearson_mtscatac_large_data_files/output/Healthy_scATAC_labelTransfer_HQcells.rds")@meta.data

a <- readRDS("../../../pearson_mtscatac_large_data_files/output/PBMC_H9-RIP.rds")
b <- readRDS("../../../pearson_mtscatac_large_data_files/output/PBMC_H10-RIP.rds")

colnames(a) <- gsub("-1", "-4", colnames(a))
colnames(b) <- gsub("-1", "-5", colnames(b))

counts_h <- cbind(assays(a)[["counts"]],assays(b)[["counts"]])
md <- rbind(pbmc_pearson_md[,c("Patient", "predicted.id")], pbmc_healthy_md[,c("Patient", "predicted.id")])

# Format counts data for chromVAR
counts_go <- cbind(counts, counts_h)[,rownames(md)]
gr <- rowRanges(a)

# Make a summarized experiment for CV
SE <- SummarizedExperiment(
  assays = list(counts = counts_go), 
  rowRanges = gr,
  colData = md
)

rm(a); rm(b); rm(counts_h); rm(counts)

# Run chromVAR
SE <- filterPeaks(SE)
mm <- motifmatchr::matchMotifs(human_pwms_v2, SE, genome = BSgenome.Hsapiens.UCSC.hg19)
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
dev <- computeDeviations(object = SE,  annotations = mm)
colData(dev)$predicted.id <- gsub("Activated B-cells", "Naive B-cells", colData(dev)$predicted.id )
saveRDS(dev, file = "../../../pearson_mtscatac_large_data_files/output/PBMC_scATAC_healthy_pearson_deviations.rds")

# Compare healthy / normal
dev <- readRDS("../../../pearson_mtscatac_large_data_files/output/PBMC_scATAC_healthy_pearson_deviations.rds")
df <- data.frame(colData(dev),t(assays(dev)[["deviations"]][rowData(dev)$name %in% c("JUN", "FOS", "JUNB", "FOSB", "KLF4"),]))

df %>% group_by(predicted.id, Patient) %>% 
  summarize(
    KLF4 = median(ENSG00000136826_LINE939_KLF4_D),
    FOSB = median(ENSG00000125740_LINE393_FOSB_D), 
    FOS = median(ENSG00000170345_LINE476_FOS_D_N4), 
    JUN = median(ENSG00000177606_LINE517_JUN_D_N5), 
    JUNB = median(ENSG00000171223_LINE487_JUNB_D_N3)) %>% data.frame()
