library(Seurat)
library(Signac)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(GenomeInfoDb)
library(BuenColors)
library(EnsDb.Hsapiens.v75)
library(ggbeeswarm)
library(viridis)
set.seed(1)

# Import ATAC data
metadata <- fread("../../../pearson_large_data_files/input/pbmcs_scatac/aggrs/P3H2_singlecell.csv") %>%
  filter(cell_id != "None")
metadata <- data.frame(metadata)
rownames(metadata) <- metadata[[1]]
mat <- Read10X_h5("../../../pearson_large_data_files/input/pbmcs_scatac/aggrs/P3H2_filtered_peak_bc_matrix.h5")
fragments_file <- "../../../pearson_large_data_files/input/pbmcs_scatac/aggrs/P3H2_fragments.tsv.gz"

# Process meta data
metadata <- metadata[grepl("-3", rownames(metadata)),]
mat <- mat[,grepl("-3", colnames(mat))]
dim(mat)
dim(metadata)

df_PT3 <- fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id)) %>%
  dplyr::filter(reads_all>=10) %>% 
  mutate(scale_heteroplasmy = scale(heteroplasmy))

df_MDS <- fread("../../pt3_chr7_del_scatac/output/Pearson-PBMC.chr7DelQC.tsv")
metadata_full <- merge(merge(df_PT3, metadata, by= "barcode"), df_MDS, by.x = "barcode", by.y  = "V4")
rownames(metadata_full) <- metadata_full$barcode

# Create Seurat object
mat <- mat[,rownames(metadata_full)]

CA <- CreateChromatinAssay(
  counts = mat,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = fragments_file,
  min.cells = 1,
  min.features = 1
)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"
Annotation(CA) <- annotations

# Create Seurat object
pbmc <- CreateSeuratObject(
  counts = CA,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata_full
)
DefaultAssay(pbmc) <- "peaks"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q25')
pbmc <- RunSVD(
  object = pbmc,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:20)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:20)
pbmc <- FindClusters(object = pbmc, resolution = 0.8)

DimPlot(object = pbmc, label = TRUE) 
FeaturePlot(object = pbmc, "heteroplasmy") +
  scale_color_viridis()

FeaturePlot(object = pbmc, "pct_in_del") +
  scale_color_viridis()

gene.activities <- GeneActivity(pbmc)
pbmc[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_ACTIVITY)
)

DefaultAssay(pbmc) <- "ACTIVITY"
FeaturePlot(object = pbmc, c( "CD4", "CD8A", "MS4A1", "CXCL14",
                              "CD3E", "LEF1","TREM1",
                              "CCL5", "CCR7", 
                              "TOX", "ADAM23", "ZNF462", "IKZF2", "CR1", "CR2"),
            max.cutoff = "q95")

DefaultAssay(pbmc) <- "ACTIVITY"
FindMarkers(pbmc, ident.1 = "15") %>% head(20)

FindMarkers(pbmc, ident.1 = "3", ident.2 = "4") %>% head(20)
FindMarkers(pbmc, ident.1 = "3") %>% head(50)
saveRDS(pbmc, file = "../../../pearson_large_data_files/output/pearson_pbmc_pt3.rds")
