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
mdf <- merge(summary_csv, het_df, all.x = TRUE, by = "row.names")
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
pbmc2 <- subset(pbmc2, subset = mtDNA_depth > 10 & passed_filters > 1000 & pct_reads_in_peaks > 40)
dim(pbmc2)

# Add patient name too 
ptid <- substr(colnames(pbmc2), 18, 18)
pbmc2@meta.data$cultureID <- ptid
pbmc2@meta.data$timepoint <- case_when(ptid == "1" ~ "aPBMC", ptid  %in% c("2", "3") ~ "Day14",ptid == "4" ~ "Day21", TRUE ~ "Healthy")

table(pbmc2@meta.data$timepoint )

pbmc2 <- FindTopFeatures(pbmc2, min.cutoff = 'q50')

# Run linear dimension reduction (very similar to PCA)
pbmc2 <- RunSVD(
  object = pbmc2,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

# Commands for a 2D embedding and graph-based clustering
set.seed(2020)
pbmc2 <- RunUMAP(object = pbmc2, reduction = 'lsi', dims = 2:30)
DimPlot(object = pbmc2, group.by = "timepoint" )

DefaultAssay(pbmc2) <- "peaks"
pbmc2 <- FindNeighbors(pbmc2, reduction = "lsi")
pbmc2 <- FindClusters(object = pbmc2, resolution = 0.5)
DimPlot(object = pbmc2,  label = TRUE)

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

# Add ADT and normalize
pbmc2[["ADT"]] <- CreateAssayObject(counts = tag_mat[,colnames(tag_mat) %in% colnames(pbmc2)], assay = "ADT")
DefaultAssay(pbmc2) < "ADT"
pbmc2 <- NormalizeData(pbmc2, assay = "ADT", normalization.method = "CLR")
pbmc2 <- ScaleData(pbmc2, assay="ADT")

DefaultAssay(pbmc2) < "peaks"
FeaturePlot(object = pbmc2, features = "heteroplasmy" ) +
  scale_color_viridis()

set.seed(1)
cormat <- cor(pbmc2@meta.data$heteroplasmy,t(data.matrix(pbmc2@assays$ADT@data)),
              use = "pairwise.complete")
perm_cormat <- cor(sample(pbmc2@meta.data$heteroplasmy), t(data.matrix(pbmc2@assays$ADT@data)),
                   use = "pairwise.complete")

obs_df <- data.frame(
  gene = colnames(cormat),
  cor = cormat[1,]
) %>% dplyr::filter(!is.na(cor)) %>% arrange(desc(cor)) %>% mutate(rank = 1:n(), what = "observed") 
obs_df


FeaturePlot(object = pbmc2, features = c("CD3D", "CD4", "CD8A", "CD55", "CD27", "CD45RA", "ITGA2", "ITGA4", "CD63", "CD45RO"),
            sort.cell = TRUE, max.cutoff = "q95") & scale_color_viridis()

DefaultAssay(pbmc2) <- "ADT"
pbmc2@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(rounded_heteroplasmy = round(mean(heteroplasmy))) %>%
  write.table(row.names = FALSE, sep = "\t", quote = FALSE)

ggplot(pbmc2@meta.data, aes(x = cultureID, y = heteroplasmy)) +
  geom_violin()

# add the gene information to the object
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

# Pearson correlation to see factors most associated with activation / not
set.seed(1)
cormat <- cor(pbmc2@meta.data$heteroplasmy,t(data.matrix(pbmc2@assays$RNA@data)),
              use = "pairwise.complete")
perm_cormat <- cor(sample(pbmc2@meta.data$heteroplasmy), t(data.matrix(pbmc2@assays$RNA@data)),
                   use = "pairwise.complete")

obs_df <- data.frame(
  gene = colnames(cormat),
  cor = cormat[1,]
) %>% dplyr::filter(!is.na(cor)) %>% arrange(desc(cor)) %>% mutate(rank = 1:n(), what = "observed") 
obs_df

write.table(obs_df, file = "../output/heteroplasmy_Tcell_association_ranking.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

perm_df <- data.frame(
  gene = colnames(perm_cormat),
  cor = perm_cormat[1,]
) %>% dplyr::filter(!is.na(cor)) %>% arrange(desc(cor)) %>% mutate(rank = 1:n(), what = "permuted") 

p1 <- ggplot(obs_df, aes(x = rank, y = cor)) + 
  geom_point(size = 0.2) + scale_color_manual(values = c("black", "firebrick"))  + 
  geom_point(inherit.aes = FALSE, data = perm_df, aes(x = rank, y = cor), color = "lightgrey", size = 0.2) +
  labs(x = "Rank sorted genes", y = "Correlation") + 
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") + theme_void()
cowplot::ggsave2(p1, file = "../plots/Tcell_heteroplasmy_assoc.png", width = 3, height = 3, dpi = 500)

saveRDS(pbmc2, file = "../../../pearson_large_data_files/output/Tcell_scATAC_culture.rds")
