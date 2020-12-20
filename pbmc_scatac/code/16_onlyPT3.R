library(Seurat)
library(Signac)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(GenomeInfoDb)
library(BuenColors)
library(ggbeeswarm)
source("../functions/import_counts_mat_10x.R")
set.seed(1)

# Import ATAC data
counts <- import_counts_mat_10x("../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_filtered_peak_bc_matrix/")
metadata <- fread("../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_singlecell.csv") %>%
  filter(cell_id != "None")
metadata <- data.frame(metadata)
rownames(metadata) <- metadata[[1]]

metadata <- metadata[grepl("-3", rownames(metadata)),]
counts <- counts[,grepl("-3", colnames(counts))]

df_BCI <- fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(reads_all>=10) %>% 
  dplyr::filter(deletion == "del6073-13095") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id)) %>%
  mutate(scale_heteroplasmy = scale(heteroplasmy))

df_CCF <- fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(reads_all>=10) %>% 
  dplyr::filter(deletion == "del8482-13447") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id)) %>%
  mutate(scale_heteroplasmy = scale(heteroplasmy))

df_PT3 <- fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id)) %>%
  dplyr::filter(reads_all>=10) %>% 
  mutate(scale_heteroplasmy = scale(heteroplasmy))

metadata_full <- merge(df_PT3, metadata, by = "barcode")
rownames(metadata_full) <- metadata_full$barcode
counts <- counts[,rownames(metadata_full)]

CA <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = '../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_fragments.tsv.gz',
  min.cells = 0,
  min.features = 0
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

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q25')
pbmc <- RunSVD(
  object = pbmc,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)

DefaultAssay(pbmc) <- "peaks"
pbmc <- FindClusters(object = pbmc, resolution = 0.8)
DimPlot(object = pbmc, label = TRUE) 
FeaturePlot(object = pbmc, "heteroplasmy") +
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
FeaturePlot(object = pbmc, c( "CD4", "CD8A", "MS4A1",
                              "CD3E", "LEF1","TREM1",
                              "CCL5", "CCR7", 
                              "TOX", "ADAM23", "ZNF462", "IKZF2", "CR1", "CR2"),
            max.cutoff = "q95")

FeaturePlot(object = pbmc, "reads_all") +
  scale_color_viridis()


FindMarkers(pbmc, ident.1 = "5", ident.2 = "0") %>% head(20)
FindMarkers(pbmc, ident.1 = "6") %>% head(50)

