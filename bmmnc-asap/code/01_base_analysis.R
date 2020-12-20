library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(BuenColors)
library(dplyr)
library(Matrix)
library(viridis)
library(patchwork)
set.seed(1)

# Import data
import_deletion <- function(library_n){
  tab <- read.table(paste0("../data/mgatk_del/del_Pearson_BM_MNC_",library_n,".deletion_heteroplasmy.tsv"), header = TRUE)
  tab$cell_id <- gsub("-1", paste0("-", library_n),tab$cell_id )
  tab %>% filter(reads_all >= 10 & deletion == "del10381-15407" & version == "improved")
}

del_df <- rbind(import_deletion(1), import_deletion(2), import_deletion(3), import_deletion(4), import_deletion(5))
rownames(del_df) <- del_df$cell_id
counts <- Read10X_h5(filename = "../../../pearson_mtscatac_large_data_files/input/bonemarrow_mnc_asap/Pearson_ASAP_filtered_peak_bc_matrix.h5")

# Slower but convenient for the row parsing
metadata <- read.csv(
  file = "../data/Pearson_ASAP_singlecell.csv.gz",
  header = TRUE,
  row.names = 1
)

protein_mat <- readRDS("../data/Pearson_ASAP_allProteinCounts.rds")

# Do some quality control
# Not as stringent as it could be but reasonable

protein_QC <- data.frame(
  barcode = colnames(protein_mat),
  control_count = colSums(protein_mat[grepl("sotype", rownames(protein_mat)),]),
  total_count = colSums(protein_mat),
  pass_mito_qc = colnames(protein_mat) %in% del_df$cell_id
) %>% filter(pass_mito_qc)

# Set some QC thresholds
ggplot(protein_QC, aes(x = log10(total_count), y =control_count)) +
  geom_point()

ggplot(metadata %>% filter(cell_id != "None"), aes(x = log10(passed_filters), y = (peak_region_fragments)/passed_filters *100)) +
  geom_point()

protein_qc_barcodes <- protein_QC$barcode[protein_QC$total_count > 150 & protein_QC$control_count < 10]
mito_qc_barcodes <- protein_QC$barcode[protein_QC$total_count > 150 & protein_QC$control_count < 10]
atac_qc_barcodes <- metadata %>% filter(cell_id != "None") %>% filter(passed_filters > 10^3 & (peak_region_fragments/passed_filters) > 0.25) %>% rownames()

# Derive consensus barcodes based on QC of all 3 modalities 
consensus_barcodes <- intersect(protein_qc_barcodes, intersect(mito_qc_barcodes, atac_qc_barcodes))
length(consensus_barcodes) # over 20k which is nice  
counts_filt <- counts[,consensus_barcodes]
protein_mat_filt <- protein_mat[,consensus_barcodes]
del_df_filt <- del_df[consensus_barcodes,]
metadata_filt <- metadata[consensus_barcodes,]
full_meta_data <- cbind(metadata_filt,del_df_filt)

# Set up Signac/Seurat object
# Note: make sure at least Seurat v 3.2 and Signac v 1.0 is installed...
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"

CA <- CreateChromatinAssay(
  counts = counts_filt,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = '../../../pearson_mtscatac_large_data_files/input/bonemarrow_mnc_asap/pearson_asap_fragments.tsv.gz',
  min.cells = 0,
  min.features = 0
)


# add the gene information to the object
Annotation(CA) <- annotations

full_meta_data_del7 <- merge(full_meta_data, data.table::fread("../../pt3_chr7_del/output/Pearson-ASAP.chr7DelQC.tsv"), by.x = "row.names", by.y = "V4")
rownames(full_meta_data_del7) <- full_meta_data_del7$Row.names
pearson_asap <- CreateSeuratObject(
  counts = CA,
  assay = "peaks",
  meta.data = full_meta_data_del7
)

# Add the protein data to the seurat object and normalize
pearson_asap[["ADT"]] <-  CreateAssayObject(counts = protein_mat_filt)
pearson_asap <- NormalizeData(pearson_asap, assay = "ADT", normalization.method = "CLR")
pearson_asap <- ScaleData(pearson_asap, assay = "ADT")

# Do ATAC-seq dimension reduction
pearson_asap <- RunTFIDF(pearson_asap)
pearson_asap <- FindTopFeatures(pearson_asap, min.cutoff = 'q50')
pearson_asap <- RunSVD(pearson_asap)
DepthCor(pearson_asap)

# Make embedding and clusters
pearson_asap <- RunUMAP(object = pearson_asap, reduction = 'lsi', dims = 2:30)
pearson_asap <- FindNeighbors(object = pearson_asap, reduction = 'lsi', dims = 2:30)
pearson_asap <- FindClusters(object = pearson_asap, verbose = FALSE, algorithm = 3)
p1 <- DimPlot(object = pearson_asap, label = TRUE) + NoLegend()
p1

# Look at heteroplasmy
p2 <- FeaturePlot(object = pearson_asap, "heteroplasmy") +
  scale_color_viridis()
p2 

p1 + p2

# Viz the ASAP tags
DefaultAssay(pearson_asap) <- "ADT"
FeaturePlot(object = pearson_asap, c("CD71", "CD3-1", "CD19", "CD14", "CD8a",
                                     "CD35", "CD4-1", "CD21"),
            max.cutoff = "q90", cols =jdb_palette("brewer_spectra")) 
FindMarkers(pearson_asap, ident.1 = "0", ident.2 = "10") %>% head(20)
FindMarkers(pearson_asap, ident.1 = "2", ident.2 = "8")%>% head(20)

FeaturePlot(object = pearson_asap, c("CD25", "CD31"),
            max.cutoff = "q90", cols =jdb_palette("brewer_spectra")) 

# add the gene activity matrix to the Seurat object as a new assay and normalize it
gene.activities <- GeneActivity(pearson_asap)
pearson_asap[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
pearson_asap <- NormalizeData(
  object = pearson_asap,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(pearson_asap$nCount_ACTIVITY)
)

DefaultAssay(pearson_asap) <- "ACTIVITY"
FeaturePlot(object = pearson_asap, c("TOX", "ADAM23", "ZNF462", "IKZF2", "CR1", "CR2"),
            max.cutoff = "q95")

# Annotate with MDS stuff
FeaturePlot(object = pearson_asap, "pct_in_del", max.cutoff = "q95", min.cutoff = "q05") +
  scale_color_viridis()

saveRDS(pearson_asap, file = "../../../pearson_mtscatac_large_data_files/output/asap/pearson_asap_master_object.rds")





