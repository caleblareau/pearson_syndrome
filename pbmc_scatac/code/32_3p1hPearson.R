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

DimPlot(object = pbmc2, label = TRUE)
DimPlot(object = pbmc2, group.by = "Patient" )

pbmc2$hetero2 <- ifelse(is.na(pbmc2$heteroplasmy), 100, pbmc2$heteroplasmy)
pbmc2$sh2 <- ifelse(is.na(pbmc2$scale_heteroplasmy), 100, pbmc2$scale_heteroplasmy)

FeaturePlot(pbmc2,  features= c("sh2"), split.by = "Patient") & 
  scale_color_gradientn(colors = jdb_palette("brewer_spectra"))

boo <- pbmc2@meta.data$Disease == "Pearson" & (pbmc2@meta.data$seurat_clusters %in% c(0,1,2,4,6,7,8,9,10))


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
cormat <- cor(pbmc2@meta.data$scale_heteroplasmy[boo],t(data.matrix(pbmc2@assays$RNA@data[,boo])),
              use = "pairwise.complete")

# Filter to interesting genes
interesting_features <- c('NCAM1', 'CD4', 'CD8A', 'MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ', 'CD14', 'CCR5', 'CXCR6','FYCO1', 'KLF7')
FeaturePlot(pbmc2, features = interesting_features)

DefaultAssay(pbmc2) <- "RNA"
FindMarkers(pbmc2, ident.1 = "6")
FindMarkers(pbmc2, ident.1 = "4", ident.2 = "2")


if(FALSE){
  # Do azimuth label transfer
  reference <- SeuratDisk::LoadH5Seurat("../../../pearson_large_data_files/input/pbmc_multimodal.h5seurat")
  
  # Look at the typical seurat performance
  rna <- CreateSeuratObject(
    counts = reference@assays$SCT@counts,
    meta.data = reference@meta.data
  )
  rna <- NormalizeData(rna)
  rna <- FindVariableFeatures(rna, nfeatures = 3000)
  rna <- ScaleData(rna)
  rna <- RunPCA(rna, npcs = 30)
  rm(reference)
  
  DefaultAssay(pbmc2) <- "RNA"
  pbmc2 <- pbmc2 %>% FindVariableFeatures() %>% NormalizeData() %>% ScaleData()
  
  transfer.anchors <- FindTransferAnchors(reference = rna, query = pbmc2, features = VariableFeatures(object = rna), 
                                          reference.assay = "RNA", query.assay = "RNA", reduction = "cca")
  
  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$celltype.l2, 
                                       weight.reduction = pbmc2[["harmony"]], dims = 2:50)
  
  pbmc2$transfered_cluster <- celltype.predictions$predicted.id
  
  DimPlot(object = pbmc2, group.by = "transfered_cluster", label = TRUE )
}

#-------------------
# Look at motifsFeaturePlot(pbmc2,  features= c("pct_reads_in_peaks"))

DefaultAssay(pbmc2) <- 'peaks'

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = FALSE)
)
motif.matrix <- CreateMotifMatrix(
  features = pbmc2@assays$peaks@ranges,
  pwm = pwm,
  genome = 'hg19',
  use.counts = FALSE
)

pbmc2 <- AddMotifs(
  object = pbmc2,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  pfm = pwm
)
pbmc2 <- RunChromVAR(
  object = pbmc2,
  genome = BSgenome.Hsapiens.UCSC.hg19
)
Idents(pbmc2) <- "cluster_pheno"
DefaultAssay(pbmc2) <- 'RNA'
FindMarkers(pbmc2, ident.1 = "0_Pearson", ident.2 = "0_Healthy")


FeaturePlot(pbmc2, features = c("FOSB"),  split.by = "Patient") &
  scale_color_gradientn(colors = jdb_palette("solar_extra"))

sapply(ms, name) %>% sort()
# create a Motif object and add it to the assay


# Retain the actual motif matched positions, down to the base pair
DefaultAssay(pbmc2) <- 'peaks'
motif.positions <- matchMotifs(
  pwms = pwm,
  subject = granges(pbmc2),
  out = 'positions',
  genome = 'hg19'
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  positions = motif.positions,
  pwm = pwm
)

# Store all of this back in the Seurat object
pbmc2 <- SetAssayData(
  object = pbmc2,
  slot = 'motifs',
  new.data = motif
)


pbmc2 <- Footprint(
  object = pbmc2,
  motif.name = c("MA0139.1"),
  genome = BSgenome.Hsapiens.UCSC.hg19,
  in.peaks = TRUE
)
pbmcT <- subset(pbmc2, seurat_clusters %in% c(0))

PlotFootprint(pbmc2, features = c("MA0139.1"))


