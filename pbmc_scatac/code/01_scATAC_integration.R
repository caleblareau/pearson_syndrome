library(Seurat)
library(Signac)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(GenomeInfoDb)
library(harmony)
library(BuenColors)
library(ggbeeswarm)
source("../functions/import_counts_mat_10x.R")
set.seed(1)

# Import ATAC data
counts <- import_counts_mat_10x("../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_filtered_peak_bc_matrix/")
metadata <- fread("../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_singlecell.csv") %>%
  filter(cell_id != "None")
rownames(metadata) <- metadata[[1]]
metadata <- data.frame(metadata)

# Create Seurat object
pbmc2 <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)

head(pbmc2@meta.data)
dim(pbmc2)
rm(counts)
rm(metadata)

# import mtDNA depths to filter cells
bci <- fread("../data/depth/Pearson-PBMC-BCI_v12-mtMask_mgatk.depthTable.txt"); bci[["V1"]] <- gsub("-1", "-1",  bci[["V1"]])
ccf <- fread("../data/depth/Pearson-PBMC-CCF_v12-mtMask_mgatk.depthTable.txt"); ccf[["V1"]] <- gsub("-1", "-2",  ccf[["V1"]])
pt3 <- fread("../data/depth/Pearson-PBMC-PT3_v12-mtMask_mgatk.depthTable.txt"); pt3[["V1"]] <- gsub("-1", "-3",  pt3[["V1"]])

# Append mito depth to Seurat object
depth_vec <- c(bci$V2, ccf$V2, pt3$V2)
names(depth_vec) <- as.character(c(bci$V1, ccf$V1, pt3$V1))
pbmc2@meta.data$mtDNA_depth <- depth_vec[rownames(pbmc2@meta.data)]
summary(pbmc2@meta.data$mtDNA_depth)

# Add patient name too 
ptid <- substr(colnames(pbmc2), 18, 18)
pbmc2@meta.data$Patient_ID <- ptid
pbmc2@meta.data$Patient <- case_when(ptid == "1" ~ "BCI", ptid == "2" ~ "CCF", ptid == "3" ~ "PT3")

# Filter cells due to quality control via Signac functions
library(patchwork)
fragment.path <- '../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_fragments.tsv.gz'

pbmc2 <- SetFragments(
  object = pbmc2,
  file = fragment.path
)
pbmc2 <- NucleosomeSignal(object = pbmc2)
pbmc2$pct_reads_in_peaks <- pbmc2$peak_region_fragments / pbmc2$passed_filters * 100
pbmc2$blacklist_ratio <- pbmc2$blacklist_region_fragments / pbmc2$peak_region_fragments

# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v75)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

tss.ranges <- GRanges(
  seqnames = seqnames(gene.ranges),
  ranges = IRanges(start = start(gene.ranges), width = 2),
  strand = strand(gene.ranges)
)

seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')

# to save time use the first 2000 TSSs
pbmc2 <- TSSEnrichment(object = pbmc2, tss.positions = tss.ranges[1:2000])

# Final filter
pbmc2 <- subset(
  x = pbmc2,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 40 &
    nucleosome_signal < 10 &
    TSS.enrichment > 2 & 
    mtDNA_depth > 20
)
dim(pbmc2)

# Do standard dimension
pbmc2 <- RunTFIDF(pbmc2)
pbmc2 <- FindTopFeatures(pbmc2, min.cutoff = 'q75')

pbmc2 <- RunSVD(
  object = pbmc2,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

# Do the clustering, etc. for the basic run
pbmc2 <- FindNeighbors(object = pbmc2, reduction = 'lsi', dims = 2:30)
pbmc2 <- FindClusters(object = pbmc2, verbose = FALSE)
pbmc2 <- RunUMAP(object = pbmc2, reduction = 'lsi', dims = 2:30)
DimPlot(pbmc2)
patient <- as.character(pbmc2@meta.data$Patient)
umap_df <- data.frame(
  pbmc2@reductions$umap@cell.embeddings,
  patient = patient
)
pbmc2 <- AddMetaData(object = pbmc2, metadata = patient, col.name = "patient")
pbmc_no_integration <- pbmc2

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = patient)) +
  geom_point()

# Use Hamrony to fix batch effects
pbmc2 <- RunHarmony(
  object = pbmc2,
  group.by.vars = 'patient',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)

# Recompute characteristics using the Harmony framework 
pbmc2 <- RunUMAP(pbmc2, dims = 2:30, reduction = 'harmony')
pbmc2 <- FindNeighbors(object = pbmc2, reduction = 'harmony', dims = 2:30)
pbmc2 <- FindClusters(object = pbmc2, verbose = FALSE, resolution = 0.5)
DimPlot(object = pbmc2, label = FALSE, group.by = "seurat_clusters") 
DimPlot(object = pbmc2, label = FALSE, group.by = "Patient") 

df2 <- data.frame(pbmc_no_integration@reductions$umap@cell.embeddings); colnames(df2) <- c("UMAP1_bad", "UMAP2_bap")
pbmc2@meta.data$UMAP1_noIntegration <- df2[,1]
pbmc2@meta.data$UMAP1_noIntegration <- df2[,2]

# Do integration
pbmc2[['RNA']] <- CreateAssayObject(counts = readRDS("../../../pearson_mtscatac_large_data_files/output/Perason_aggr_26April2020_geneActivityScores.rds")[,colnames(pbmc2)])
pbmc2 <- NormalizeData(
  object = pbmc2,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc2$nCount_RNA)
)
rm(pbmc_no_integration)
DefaultAssay(object = pbmc2) <- "RNA"
#saveRDS(pbmc2, file = "../../../pearson_mtscatac_large_data_files/output/6May_Pearson_patients_processedSeuratObject_ATAC.rds")


# Import RNA-seq data to handle the co-embedding
pbmc_rna <- readRDS("../../../pearson_mtscatac_large_data_files/output/5March-PearsonRNAseq-integration.rds")
pbmc_rna$tech <- "rna"

# Transfer anchors beteween ATAC-seq and RNA-seq
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc2,
  features = VariableFeatures(object = pbmc_rna),
  reference.assay = "RNA", query.assay = "RNA",
  reduction = 'cca'
)

# Predict cell types using data transfering
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc_rna$seurat_clusters, 
                                     weight.reduction = pbmc2[["lsi"]])
pbmc2 <- AddMetaData(pbmc2, metadata = celltype.predictions)
prediction.score.max <- pbmc2$prediction.score.max
hist(prediction.score.max)
abline(v = 0.5, col = "red")
pbmc2$predicted.id <- pbmc2$predicted.id

DimPlot(pbmc2, group.by = "predicted.id")

saveRDS(pbmc2, file = "../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")
