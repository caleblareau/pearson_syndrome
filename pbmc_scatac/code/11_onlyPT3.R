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
metadata <- fread("../data/patient_singlecell/PBMC_PT3_v12-mtMask.singlecell.csv.gz") %>%
  dplyr::filter(cell_id != "None")
metadata <- data.frame(metadata)
rownames(metadata) <- metadata[[1]]
mat <- Read10X_h5("../data/peaks_by_cells/PT3_filtered_peak_bc_matrix.h5")
fragments_file <- "../../../pearson_large_data_files/input/pbmcs_scatac/fragments/hg19/PT3_fragments.tsv.gz"

# Process meta data
dim(mat)
dim(metadata)

df_PT3 <- fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407") %>%
  dplyr::filter(reads_all>=10) 

metadata_full <- data.frame(merge(df_PT3, metadata, by.x = "cell_id", by.y = "barcode"))
consensus_bc <- intersect(metadata_full$cell_id, colnames(mat))
metadata_full <- metadata_full %>% dplyr::filter(cell_id %in% consensus_bc)
rownames(metadata_full) <- metadata_full$cell_id

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

DefaultAssay(pbmc) <- "peaks"
pbmc <- FindClusters(object = pbmc, resolution = 0.7 )

DimPlot(object = pbmc, label = TRUE) 

FeaturePlot(object = pbmc, "heteroplasmy") +
  scale_color_viridis()

# Annotate clusters
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
                              "CD3E", "LEF1","TREM1", "FCGR3A", 
                              "CCL5", "CCR7", "CD3D", "EOMES", "NKG7", 
                              "TNFRSF17"),
            max.cutoff = "q95", sort.cell = TRUE)

DefaultAssay(pbmc) <- "ACTIVITY"
#FindMarkers(pbmc, ident.1 = "14") %>% head(20)

vec <- c("0" = "CD4_Naive",
         "1" = "CD8_Naive",
         "2" =  "CD4_Naive",
         "3" = "NKcell",
         "4" = "NKcell",
         "5" = "CD4_RTE",
         "6" = "CD8_cytotoxic",
         "7" = "NKcell",
         "8" = "Bcell_naive",
         "9" = "Monocytes_CD14",
         "10" = "CD8_RTE",
         "11" = "CD8_effectormemory",
         "12" = "Bcell_memory",
         "13" = "Monocytes_CD16",
         "14" = "pDCs",
         "15" = "Bcell_plasma")

pbmc$cl_anno <- vec[as.character(pbmc$seurat_clusters)]
#FindMarkers(pbmc, ident.1 = "6", ident.2 = "11") %>% head(20)
#FindMarkers(pbmc, ident.1 = "3") %>% head(50)


set.seed(7)
p1 <- DimPlot(pbmc, group.by = "cl_anno") +
  scale_color_manual(values = sample(jdb_palette("corona", 14))) +
  theme_void()
cowplot::ggsave2(p1, file = "../plots/PT3_PBMC_clusters-revision.pdf", width = 5, height = 3.5)

p2 <- FeaturePlot(object = pbmc, "heteroplasmy") +
  scale_color_viridis() + theme_void() + 
  theme(legend.position = "none")
cowplot::ggsave2(p2, file = "../plots/PT3_PBMC_heteroplasmy.pdf", width = 3.5, height = 3.5)



### Remake
pbmc2 <- readRDS("../../../pearson_large_data_files/output/PBMC_scATAC_3P1H-15FEB2021.rds")
old_umap <- data.frame(cell=rownames(pbmc2@meta.data), pbmc2@reductions$umap@cell.embeddings)
new_clustering <- data.frame(
  cell = gsub("-1", "-3", rownames(pbmc@meta.data)),
  pbmc@meta.data
)
mdf <- merge(old_umap, new_clustering, by = "cell")
set.seed(7)
x <- ggplot(mdf, aes(x = UMAP_1, y = UMAP_2, color = cl_anno)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = sample(jdb_palette("corona", 14))) +
  theme_void() + theme(legend.position = "none") + ggtitle("hm")
cowplot::ggsave2(x, file = "../plots/old_embedding_new_colors.pdf", width = 3.5, height = 3.5)

ggplot(mdf, aes(x = cl_anno, y = heteroplasmy)) + geom_boxplot()

pG <- ggplot(mdf %>% arrange(desc(X7del)), aes(x = UMAP_1, y = UMAP_2, color = X7del < 0.01)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("lightgrey", "red")) +
  theme_void() + theme(legend.position = "none") + ggtitle("hm")
cowplot::ggsave2(pG, file = "../plots/MDS_viz_pbmcs.pdf", width = 3.5, height = 3.5)
