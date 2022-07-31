library(Seurat)
library(Signac)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(GenomeInfoDb)
library(harmony)
library(BuenColors)
library(ggbeeswarm)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(viridis)
library(Matrix)

set.seed(1)

# Import ATAC data
summary_csv <- fread("../data/tcell-aggr-singlecell.csv.gz") %>%
  dplyr::filter(is__cell_barcode != "None")
summary_csv <- data.frame(summary_csv)
rownames(summary_csv) <- summary_csv[[1]]

df_4 <- fread("../data/deletion_heteroplasmy/PearsonCulture_Day21_del.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407"& version == "improved") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-4", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all)

# Combine everything
het_df <- rbind(df_4) %>% data.frame() # df_1
rownames(het_df) <- het_df[["barcode"]]
mdf <- merge(summary_csv, het_df,  by = "row.names")
rownames(mdf) <- mdf[["Row.names"]]
dim(mdf)
dim(summary_csv)
dim(het_df)
mdf <- mdf[complete.cases(mdf),]
mdf <- mdf %>% dplyr::filter(passed_filters > 1000 & reads_all >= 10)

# Now create the Seurat object
library(EnsDb.Hsapiens.v86)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "hg38"
seqlevelsStyle(annotations) <- 'UCSC'
annotations
annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))
annotations

# Set it up
mat <- Read10X_h5("../../../pearson_large_data_files/input/tcell-culture/tcell-aggr-filtered_peak_bc_matrix.h5")
fragments_file <-"../../../pearson_large_data_files/input/tcell-culture/tcell-aggr-fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = mat[,colnames(mat) %in% mdf$Row.names],
  sep = c(":", "-"),
  fragments = fragments_file,
  min.cells = 0,
  min.features = 0
)

pbmc2 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = mdf
)

head(pbmc2@meta.data)
dim(pbmc2)

# Filter
pbmc2$pct_reads_in_peaks <- pbmc2$pct_reads_in_peaks <- pbmc2$peak_region_fragments / pbmc2$passed_filters *100
pbmc2 <- subset(pbmc2, subset = reads_all > 10 & passed_filters > 1000 & pct_reads_in_peaks > 40)
dim(pbmc2)

# Add patient name too 
ptid <- substr(colnames(pbmc2), 18, 18)
pbmc2@meta.data$cultureID <- ptid
pbmc2@meta.data$timepoint <- case_when(ptid == "1" ~ "aPBMC", ptid  %in% c("2", "3") ~ "Day14",ptid == "4" ~ "Day21", TRUE ~ "Healthy")

table(pbmc2@meta.data$timepoint )


# Run linear dimension reduction
DefaultAssay(pbmc2) <- "peaks"
pbmc2 <- FindTopFeatures(pbmc2, min.cutoff = 'q50')

pbmc2 <- RunSVD(
  object = pbmc2,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

# Commands for a 2D embedding and graph-based clustering
set.seed(2022)
pbmc2 <- FindNeighbors(pbmc2, reduction = "lsi", dims = 2:30, graph.name = "LSI")
pbmc2 <- RunUMAP(object = pbmc2, reduction = 'lsi', dims = 2:30)
pbmc2 <- FindClusters(object = pbmc2, resolution = 0.5, algorithm = 3, graph.name = "LSI")
DimPlot(object = pbmc2,  label = TRUE, reduction = "umap") +
  scale_color_manual(values = jdb_palette("corona")) 

ggplot(pbmc2@meta.data, aes(x = heteroplasmy, color = seurat_clusters)) +
  stat_ecdf() +
  scale_color_manual(values = jdb_palette("corona")) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))

# add the gene information to the object
DefaultAssay(pbmc2) <- "peaks"
Annotation(pbmc2) <- annotations

# create a gene activity by cell matrix
gene.activities <- GeneActivity(pbmc2)

pbmc2[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc2 <- NormalizeData(
  object = pbmc2,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc2$nCount_RNA)
)

DefaultAssay(pbmc2) <- "RNA"
FeaturePlot(pbmc2, features = c("LEF1", "CD8B", "CD4", "CD8A","CXCL13", "CD86", "heteroplasmy"), reduction = "umap", max.cutoff = "q99",
            sort.cell = TRUE) &
  scale_color_viridis()
#saveRDS(pbmc2, file = "../../../pearson_large_data_files/output/invitro_tcell/Tcell_scATAC_culture-day21.rds")
