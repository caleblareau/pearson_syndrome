library(Seurat)
library(Signac)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(GenomeInfoDb)
library(harmony)
library(BuenColors)
library(ggbeeswarm)
library(patchwork)

set.seed(1)

# Make TSS
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


process_healthy_pbmcs <- function(exp){
  
  # Import ATAC data
  SE <- readRDS(paste0("../../../pearson_mtscatac_large_data_files/output/",exp,"-RIP.rds"))
  counts_mat <- assays(SE)[["counts"]]
  rownames(counts_mat) <- paste0("peak", as.character(1:dim(SE)[1]))
  
  metadata <- fread(paste0("../data/control_singlecell/",exp,"_v12-mtMask_singlecell.csv.gz")) %>%
    dplyr::filter(cell_id != "None")
  rownames(metadata) <- metadata[[1]]
  metadata <- data.frame(metadata)
  
  # Create Seurat object
  pbmc2 <- CreateSeuratObject(
    counts = counts_mat,
    assay = 'peaks',
    project = 'ATAC',
    min.cells = 1,
    meta.data = metadata
  )
  
  # Already filtered for depth in previous
  
  # Add patient name too 
  pbmc2@meta.data$Patient_ID <- exp
  pbmc2@meta.data$Patient <- "Healthy"
  
  # Filter cells due to quality control via Signac functions
  fragment.path <- paste0('../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/',exp,'_v12-mtMask_fragments.tsv.gz')
  
  pbmc2 <- SetFragments(
    object = pbmc2,
    file = fragment.path
  )
  
  pbmc2 <- NucleosomeSignal(object = pbmc2)
  pbmc2$pct_reads_in_peaks <- pbmc2$peak_region_fragments / pbmc2$passed_filters * 100
  pbmc2$blacklist_ratio <- pbmc2$blacklist_region_fragments / pbmc2$peak_region_fragments
  
  # to save time use the first 2000 TSSs
  pbmc2 <- TSSEnrichment(object = pbmc2, tss.positions = tss.ranges[1:2000])
  
  # Final filter
  pbmc2 <- subset(
    x = pbmc2,
    subset = peak_region_fragments > 1000 &
      peak_region_fragments < 20000 &
      pct_reads_in_peaks > 40 &
      nucleosome_signal < 10 &
      TSS.enrichment > 2 )
  pbmc2
}
h10 <- process_healthy_pbmcs("PBMC_H10")
h09 <- process_healthy_pbmcs("PBMC_H9")

h10 <- RenameCells(object = h10, new.names = gsub("-1", "-5", colnames(h10)))
h09 <- RenameCells(object = h09, new.names = gsub("-1", "-4", colnames(h09)))
peaks_keep <- intersect(
  rownames(h10@assays$peaks@counts),
  rownames(h09@assays$peaks@counts))

pbmc2  <- CreateSeuratObject(
  counts = cbind(h09@assays$peaks@counts[peaks_keep,], h10@assays$peaks@counts[peaks_keep,]),
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = rbind(h09@meta.data, h10@meta.data)
)

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



# Do integration
pbmc2[['RNA']] <- CreateAssayObject(counts = readRDS("../../../pearson_mtscatac_large_data_files/output/Healthy_26April2020_geneActivityScores.rds")[,colnames(pbmc2)])
pbmc2 <- NormalizeData(
  object = pbmc2,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc2$nCount_RNA)
)
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
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc_rna$celltype, 
                                     weight.reduction = pbmc2[["lsi"]])
pbmc2 <- AddMetaData(pbmc2, metadata = celltype.predictions)
prediction.score.max <- pbmc2$prediction.score.max
hist(prediction.score.max)
abline(v = 0.5, col = "red")
pbmc2$predicted.id <- pbmc2$predicted.id

DimPlot(pbmc2, group.by = "predicted.id")

saveRDS(pbmc2, file = "../../../pearson_mtscatac_large_data_files/output/Healthy_scATAC_labelTransfer_HQcells.rds")
