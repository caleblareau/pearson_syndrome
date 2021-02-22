library(Seurat)
library(Signac)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(GenomeInfoDb)
library(harmony)
library(BuenColors)
library(ggbeeswarm)
library(BSgenome.Hsapiens.UCSC.hg19)
library(EnsDb.Hsapiens.v75)
library(JASPAR2018)
library(TFBSTools)
library(motifmatchr)
library(matrixTests)
set.seed(1)

# Import ATAC data
summary_csv <- fread("../../../pearson_large_data_files/input/pbmcs_scatac/aggrs/P3H2_singlecell.csv") %>%
  dplyr::filter(cell_id != "None")
summary_csv <- data.frame(summary_csv)
rownames(summary_csv) <- summary_csv[[1]]

# Add heteroplasmy
# Add Heteroplasmy
df_BCI <- fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del6073-13095") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id)) %>%
  mutate(scale_heteroplasmy = scale(heteroplasmy))

df_CCF <- fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del8482-13447") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id)) %>%
  mutate(scale_heteroplasmy = scale(heteroplasmy))

df_PT3 <- fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id)) %>%
  mutate(scale_heteroplasmy = scale(heteroplasmy))

het_df <- rbind(df_BCI, df_CCF, df_PT3) %>% data.frame()
rownames(het_df) <- het_df[["barcode"]]
mdf <- merge(summary_csv, het_df, all.x = TRUE, by = "row.names")
rownames(mdf) <- mdf[["Row.names"]]
dim(mdf)
dim(summary_csv)
dim(het_df)

# Set it up
mat <- Read10X_h5("../../../pearson_large_data_files/input/pbmcs_scatac/aggrs/P3H2_filtered_peak_bc_matrix.h5")
fragments_file <- "../../../pearson_large_data_files/input/pbmcs_scatac/aggrs/P3H2_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = mat,
  sep = c(":", "-"),
  genome = 'hg19',
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

# import mtDNA depths to filter cells
bci <- fread("../data/depth/Pearson-PBMC-BCI_v12-mtMask_mgatk.depthTable.txt"); bci[["V1"]] <- gsub("-1", "-1",  bci[["V1"]])
ccf <- fread("../data/depth/Pearson-PBMC-CCF_v12-mtMask_mgatk.depthTable.txt"); ccf[["V1"]] <- gsub("-1", "-2",  ccf[["V1"]])
pt3 <- fread("../data/depth/Pearson-PBMC-PT3_v12-mtMask_mgatk.depthTable.txt"); pt3[["V1"]] <- gsub("-1", "-3",  pt3[["V1"]])

h9 <- fread("../data/depth/PBMC_H9_v12-mtMask_mgatk.depthTable.txt"); h9[["V1"]] <- gsub("-1", "-4",  h9[["V1"]])
h10 <- fread("../data/depth/PBMC_H10_v12-mtMask_mgatk.depthTable.txt"); h10[["V1"]] <- gsub("-1", "-5",  h10[["V1"]])

# Append mito depth to Seurat object
depth_vec <- c(bci$V2, ccf$V2, pt3$V2, h9$V2, h10$V2)
names(depth_vec) <- as.character(c(bci$V1, ccf$V1, pt3$V1, h9$V1, h10$V1))
pbmc2@meta.data$mtDNA_depth <- depth_vec[rownames(pbmc2@meta.data)]
summary(pbmc2@meta.data$mtDNA_depth)
pbmc2$pct_reads_in_peaks <- pbmc2$peak_region_fragments / pbmc2$passed_filters * 100
qplot(pbmc2@meta.data$pct_reads_in_peaks, log10(pbmc2@meta.data$passed_filters))

# Filter
pbmc2 <- subset(pbmc2, subset = mtDNA_depth > 10 & passed_filters > 1000 & pct_reads_in_peaks > 45)
dim(pbmc2)

# Add patient name too 
ptid <- substr(colnames(pbmc2), 18, 18)
pbmc2@meta.data$Patient_ID <- ptid
pbmc2@meta.data$Patient <- case_when(ptid == "1" ~ "BCI", ptid == "2" ~ "CCF",ptid == "3" ~ "PT3", TRUE ~ "Healthy")
pbmc2@meta.data$Disease <- case_when(ptid %in% c("1", "2", "3") ~ "Pearson", TRUE ~ "Healthy")

table(pbmc2@meta.data$Patient )

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
DimPlot(object = pbmc2, group.by = "Patient" )

pbmc2 <- RunHarmony(
  object = pbmc2,
  group.by.vars = 'Patient',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)
pbmc2 <- RunUMAP(pbmc2, dims = 2:30, reduction = 'harmony')
pbmc2 <- FindNeighbors(object = pbmc2, reduction = 'harmony', dims = 2:30) 
pbmc2 <-   FindClusters(object = pbmc2, verbose = FALSE, resolution = 0.5)
pbmc2$cluster_pheno <- paste0(pbmc2$seurat_clusters, "_", pbmc2$Disease)

corona = c("#1f77b4","#d62728","#2ca02c","#ff7f0e","#9467bd","#8c564b","#e377c2","#7f7f7f",
           "#bcbd22","#17becf","#ad494a","#e7ba52","#8ca252","#756bb1","#636363","#aec7e8",
           "#ff9896","#98df8a","#ffbb78","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7","#dbdb8d",
           "#9edae5","#e7969c","#e7cb94","#c7e9c0","#de9ed6","#d9d9d9")
DimPlot(object = pbmc2, label = TRUE)

col_vec <- c("#1f77b4","#756bb1", "#aec7e8",  "#2ca02c", "#7f7f7f", "#ff7f0e","#9467bd",
             "#8ca252", "#ffbb78", "#c7c7c7", "#e7ba52", "#8c564b"); names(col_vec) <-  as.character(0:11)

# Make general embedding plot for figure
ppng <- DimPlot(pbmc2, label = FALSE) +
  scale_color_manual(values = col_vec)+
  theme_void() + theme(legend.position = 'none') 
cowplot::ggsave2(ppng, file = "../plots/scATAC_plot_nice_colors.png", width = 5, height = 5, dpi = 500)

pbmc2$hetero2 <- ifelse(is.na(pbmc2$heteroplasmy), 100, pbmc2$heteroplasmy)
pbmc2$sh2 <- ifelse(is.na(pbmc2$scale_heteroplasmy), 100, pbmc2$scale_heteroplasmy)

phet_grid <- FeaturePlot(pbmc2,  features= c("hetero2"), split.by = "Patient") & 
  scale_color_gradientn( colours = viridisLite::viridis(256), limits = c(0,100), oob = scales::squish) &
  theme_void() & theme(legend.position = "none")
cowplot::ggsave2(phet_grid, file = "../plots/scATAC_plot_heteroplasmy.png", width = 20, height = 5.1, dpi = 500)

# Make violin plot
pbmc2@meta.data %>% dplyr::filter(Patient != "Healthy") -> plot_het
plot_het$plot_donor <- factor(plot_het$Patient, c("PT3", "CCF", "BCI"))
p1 <- plot_het %>% ggplot(aes(x = plot_donor, y = heteroplasmy)) +
  geom_violin(fill = NA) + coord_flip() +
  labs (x = "Patient", y = "% Heteroplasmy") +
  pretty_plot(fontsize = 8) + L_border()
cowplot::ggsave2(p1, file = "../plots/violins_heteroplasmy_nocluster.pdf", width = 3.5, height = 1.8)
plot_het %>% group_by(plot_donor) %>%
  summarize(count =n(), mean(heteroplasmy), median(heteroplasmy))

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"

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
boo <- pbmc2@meta.data$Disease == "Pearson" & (pbmc2@meta.data$seurat_clusters %in% c(0,1,2,4,6,7,8,9,10))
cormat <- cor(pbmc2@meta.data$scale_heteroplasmy[boo],t(data.matrix(pbmc2@assays$RNA@data[,boo])),
              use = "pairwise.complete")
perm_cormat <- cor(sample(pbmc2@meta.data$scale_heteroplasmy[boo]), t(data.matrix(pbmc2@assays$RNA@data[,boo])),
                   use = "pairwise.complete")

obs_df <- data.frame(
  gene = colnames(cormat),
  cor = cormat[1,]
) %>% dplyr::filter(!is.na(cor)) %>% arrange(desc(cor)) %>% mutate(rank = 1:n(), what = "observed") 
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

saveRDS(pbmc2, file = "../../../pearson_large_data_files/output/PBMC_scATAC_3P1H-15FEB2021.rds")
