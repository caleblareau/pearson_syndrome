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

# Add heteroplasmy
df_1 <- fread("../../pbmc_scatac/data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>% data.frame() %>% 
  dplyr::filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all)

df_2 <- fread("../data/deletion_heteroplasmy/PearsonCulture_Day14v1_del.deletion_heteroplasmy.tsv") %>% data.frame() %>% 
  dplyr::filter(deletion == "del10381-15407" & version == "improved") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all)

df_3 <- fread("../data/deletion_heteroplasmy/PearsonCulture_Day14v2_del.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407"& version == "improved") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all)

df_4 <- fread("../data/deletion_heteroplasmy/PearsonCulture_Day21_del.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407"& version == "improved") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-4", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all)

# Combine everything
het_df <- rbind(df_2, df_3) %>% data.frame() # df_1
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

# Import ADT
import_kite_counts <- function(library, id){
  mtx <- fread(paste0("../data/pstag",library,"_fc/featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0("../data/pstag",library,"_fc/featurecounts.barcodes.txt.gz"), header = FALSE)[[1]], "-", id)
  colnames(matx) <- paste0(fread(paste0("../data/pstag",library,"_fc/featurecounts.genes.txt.gz"), header = FALSE)[[1]])
  return(t(matx))
}
tag_mat <- cbind(import_kite_counts("1", "2"), import_kite_counts("2", "3"))

# Add ADT, subset, and normalize
pbmc2[["ADT"]] <- CreateAssayObject(counts = tag_mat[,colnames(tag_mat) %in% colnames(pbmc2)], assay = "ADT")
pbmc2@meta.data$nADT_control = colSums(pbmc2[["ADT"]]@counts[grepl("control", rownames(pbmc2[["ADT"]]@counts)),])
pbmc2@meta.data$nADT_total <- colSums(pbmc2[["ADT"]]@counts)

pbmc2 <- subset(pbmc2, subset = pbmc2@meta.data$nADT_control/nADT_total < 0.02 & nADT_total > 150 & nADT_total < 1e4)
dim(pbmc2)

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
pbmc2 <- FindNeighbors(pbmc2, reduction = "lsi", dims = 2:20, graph.name = "LSI")
pbmc2 <- RunUMAP(object = pbmc2, reduction = 'lsi', dims = 2:20)
pbmc2 <- FindClusters(object = pbmc2, resolution = 1.2, algorithm = 3, graph.name = "LSI")
DimPlot(object = pbmc2,  label = TRUE, reduction = "umap")

pbmc2$caleb_cluster <- case_when(
  pbmc2$seurat_clusters %in% c(6) ~ "CD8/CD45RAhi",
  pbmc2$seurat_clusters %in% c(7,12) ~ "CD4/CD45RAhi",
  pbmc2$seurat_clusters %in% c(1) ~ "CD4/CD27hi/CD45RAhi",
  pbmc2$seurat_clusters %in% c(2,5,8,10,13) ~ "CD8/CD45RAlow" ,
  pbmc2$seurat_clusters %in% c(0,3,4,9,11) ~ "CD4/CD45RAlow" 
)

DimPlot(object = pbmc2,  label = TRUE, reduction = "umap", group.by = "caleb_cluster")

p2 <- ggplot(pbmc2@meta.data, aes(x = heteroplasmy, color = caleb_cluster)) +
  stat_ecdf() +
  scale_color_manual(values = jdb_palette("corona")) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))
cowplot::ggsave2(p2, file = paste0("../plots/ecdf_tcell_clusters_day14.pdf"), width = 2.9, height = 1.5)

# Process ADT data
DefaultAssay(pbmc2) <- 'ADT'
VariableFeatures(pbmc2) <- rownames(pbmc2[["ADT"]])
pbmc2 <- NormalizeData(pbmc2, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')
FeaturePlot(object = pbmc2, features = c("CD3D", "CD4", "CD8A", "CD55", "CD27",
                                         "CD45RA", "ITGA2", "ITGA4", "CD63", "CD45RO", "heteroplasmy"),
            sort.cell = TRUE, max.cutoff = "q95", reduction = "umap") & scale_color_viridis()

FeaturePlot(object = pbmc2, features = c("CD3D", "CD4", "CD8A", "CD55", "CD27",
                                         "CD45RA", "ITGA2", "ITGA4", "CD63", "CD45RO", "heteroplasmy"),
            sort.cell = FALSE, max.cutoff = "q95", reduction = "umap") & scale_color_viridis()

# Export plots
p1 <- FeaturePlot(object = pbmc2, features = c("CD8A"),
            sort.cell = TRUE, max.cutoff = "q95", reduction = "umap") + scale_color_viridis() +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1, file = "../plots/cd8a.png", width = 3, height = 3, dpi = 400)

p1 <- FeaturePlot(object = pbmc2, features = c("heteroplasmy"),
                  sort.cell = TRUE, max.cutoff = "q95", reduction = "umap") + scale_color_viridis() +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1, file = "../plots/heteroplasmy.png", width = 3, height = 3, dpi = 400)

p1 <- FeaturePlot(object = pbmc2, features = c("CD4"),
                  sort.cell = TRUE, max.cutoff = "q95", reduction = "umap") + scale_color_viridis() +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1, file = "../plots/cd4.png", width = 3, height = 3, dpi = 400)

p1 <- FeaturePlot(object = pbmc2, features = c("CD27"),
                  sort.cell = TRUE, max.cutoff = "q95", reduction = "umap") + scale_color_viridis() +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1, file = "../plots/cd27png", width = 3, height = 3, dpi = 400)

p1 <- FeaturePlot(object = pbmc2, features = c("CD45RA"),
                  sort.cell = TRUE, max.cutoff = "q95", reduction = "umap") + scale_color_viridis() +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1, file = "../plots/cd45ra.png", width = 3, height = 3, dpi = 400)

p1 <- FeaturePlot(object = pbmc2, features = c("CD45RA"),
                  sort.cell = TRUE, max.cutoff = "q95", reduction = "umap") + scale_color_viridis() +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1, file = "../plots/cd45ra.png", width = 3, height = 3, dpi = 400)

p1 <- DimPlot(object = pbmc2, group.by = "caleb_cluster",
                  reduction = "umap") + scale_color_manual(values = jdb_palette("corona"))+
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1, file = "../plots/clabels.png", width = 3, height = 3, dpi = 400)

#######


VlnPlot(pbmc2, features = c("CD4", "CD8A", "heteroplasmy", "CD45RA"))

DefaultAssay(pbmc2) <- "ADT"
pbmc2@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(rounded_heteroplasmy = round(mean(heteroplasmy))) %>%
  write.table(row.names = FALSE, sep = "\t", quote = FALSE)

ggplot(pbmc2@meta.data, aes(x = cultureID, y = heteroplasmy)) +
  geom_violin()

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
FeaturePlot(pbmc2, features = c("LEF1", "CD8A", "CD4"), reduction = "umap")
saveRDS(pbmc2, file = "../../../pearson_large_data_files/output/Tcell_scATAC_culture.rds")
