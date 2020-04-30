library(Seurat)
library(Signac)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(uwot)
library(harmony)
library(BuenColors)
library(ggbeeswarm)
set.seed(1)

# Import ATAC data
SE1 <- readRDS("../data/BCI_PBMC_QC.rds"); colData(SE1)$Patient <- "BCI"; colnames(SE1) <- paste0("BCI_", colnames(SE1))
SE2 <- readRDS("../data/CCF_PBMC_QC.rds"); colData(SE2)$Patient <- "CCF"; colnames(SE2) <- paste0("CCF_", colnames(SE2))
SE3 <- readRDS("../data/PT3_PBMC_QC.rds"); colData(SE3)$Patient <- "PT3"; colnames(SE3) <- paste0("PT3_", colnames(SE3))
SE <- cbind(SE1, SE2, SE3)
rm(SE1); rm(SE2); rm(SE3)

# Pull out counts and meta data to make a Seurat object
counts <- assays(SE)[["counts"]]
drr <- data.frame(rowRanges(SE))
rn_names <-  paste0(drr[,"seqnames"], ":", drr[,"start"], "-", drr[,"end"])
rownames(counts) <- rn_names
metadata <- colData(SE) %>% data.frame()

# Create Seurat object
pbmc2 <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 0,
  meta.data = metadata
)

rm(counts); rm(SE)

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
pbmc2 <- FindNeighbors(object = pbmc2, reduction = 'lsi', dims = 1:30)
pbmc2 <- FindClusters(object = pbmc2, verbose = FALSE)
pbmc2 <- RunUMAP(object = pbmc2, reduction = 'lsi', dims = 1:30)
patient <- as.character(pbmc2@meta.data$Patient)
umap_df <- data.frame(
  pbmc2@reductions$umap@cell.embeddings,
  patient = patient
)
pbmc2 <- AddMetaData(object = pbmc2, metadata = patient, col.name = "patient")
pbmc_no_integration <- pbmc2

# Use Hamrony to fix batch effects
pbmc2 <- RunHarmony(
  object = pbmc2,
  group.by.vars = 'patient',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)

# Recompute characteristics using the Harmony framework 
pbmc2 <- RunUMAP(pbmc2, dims = 1:30, reduction = 'harmony')
pbmc2 <- FindNeighbors(object = pbmc2, reduction = 'harmony', dims = 1:30)
pbmc2 <- FindClusters(object = pbmc2, verbose = FALSE, resolution = 0.5)
DimPlot(object = pbmc2, label = FALSE, group.by = "seurat_clusters") 

DimPlot(object = pbmc_no_integration, label = FALSE, group.by = "patient") 
df2 <- data.frame(pbmc_no_integration@reductions$umap@cell.embeddings); colnames(df2) <- c("UMAP1_bad", "UMAP2_bap")

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = 'patient')) +
  geom_point()


if(FALSE){
  
  fragment_file_filtered <- "../data/filtered_3pbmcs_pearson.fragments.tsv.gz"
  # extract gene coordinates from Ensembl, and ensure name formatting is consistent with Seurat object 
  gene.coords <- genes(EnsDb.Hsapiens.v75, filter = ~ gene_biotype == "protein_coding")
  seqlevelsStyle(gene.coords) <- 'UCSC'
  genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
  genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)
  
  # create a gene by cell matrix
  gene.activities <- FeatureMatrix(
    fragments = fragment_file_filtered,
    features = genebodyandpromoter.coords,
    cells = colnames(pbmc2),
    chunk = 10
  )
  
  # convert rownames from chromsomal coordinates into gene names
  gene.key <- genebodyandpromoter.coords$gene_name
  names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
  rownames(gene.activities) <- gene.key[rownames(gene.activities)]
  saveRDS(gene.activities, file = "../data/filtered_3pbmcs_pearson.gene.activities.rds")
}

pbmc2[['RNA']] <- CreateAssayObject(counts = readRDS("../data/filtered_3pbmcs_pearson.gene.activities.rds"))
pbmc2 <- NormalizeData(
  object = pbmc2,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc2$nCount_RNA)
)

# Import RNA-seq data to handle the co-embedding
pbmc_rna <- readRDS("../scrna/output/30December-PearsonRNAseq-integration.rds")
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
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc_rna$celltype, 
                                     weight.reduction = pbmc2[["lsi"]])
pbmc2 <- AddMetaData(pbmc2, metadata = celltype.predictions)
prediction.score.max <- pbmc2$prediction.score.max
hist(prediction.score.max)
abline(v = 0.5, col = "red")
pbmc2$predicted.id <- pbmc2$predicted.id

#pbmc.atac.filtered <- subset(pbmc2, subset = prediction.score.max > 0.5)
#pbmc.atac.filtered$predicted.id <- factor(pbmc.atac.filtered$predicted.id, levels = unique(pbmc_rna$celltype))  # to make the colors match

DimPlot(object = pbmc2, label = FALSE, group.by = "Patient") 
DimPlot(object = pbmc2, label = FALSE, group.by = "predicted.id") 

df <- data.frame(
  pbmc2@reductions$umap@cell.embeddings,
  df2,
  pbmc2@meta.data
)
df$cell_id <- rownames(df)
het_df <- fread("../data/minimal_3donor_heteroplasmy.tsv"); colnames(het_df) <- c("cell_id", "heteroplasmy", "coverage")
merged_df <- merge(df, het_df)

saveRDS(merged_df, file = "../data/3Perason_master_integration_df.rds")

ggplot(merged_df, aes(x = UMAP_1, y = UMAP_2, color = prediction.score.max > 0.5)) +
  geom_point()

set.seed(6)
ggplot(merged_df, aes(x = UMAP_1, y = UMAP_2, color = predicted.id)) +
  geom_point() + scale_color_manual(values = sample(jdb_palette("lawhoops")))

ggplot(merged_df, aes(x = UMAP_1, y = UMAP_2, color = heteroplasmy)) +
  geom_point() + scale_color_gradientn(colors = jdb_palette("brewer_spectra"))

ggplot(merged_df %>% dplyr::filter(coverage > 20), aes(x = predicted.id, y = heteroplasmy)) +
  facet_wrap(~patient, nrow = 3) +
  geom_quasirandom() + 
  geom_boxplot(fill = NA, outlier.shape = NA, color = "firebrick") + pretty_plot() + L_border()

merged_df %>% dplyr::filter(coverage > 20 & prediction.score.max > 0.5) %>% group_by(Patient, predicted.id) %>% summarize(count = n()) %>% data.frame()
