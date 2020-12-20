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
summary_csv <- fread("../../../pearson_mtscatac_large_data_files/input/pbmcs_scatac/aggrs/pbmc_2pearson_2health_singlecell.csv") %>%
  filter(cell_id != "None")
summary_csv <- data.frame(summary_csv)
rownames(summary_csv) <- summary_csv[[1]]
mat <- Read10X_h5("../../../pearson_mtscatac_large_data_files/input/pbmcs_scatac/aggrs/pbmc_2pearson_2health_filtered_peak_bc_matrix.h5")
fragments_file <- "../../../pearson_mtscatac_large_data_files/input/pbmcs_scatac/aggrs/pbmc_2pearson_2health_fragments.tsv.gz"
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
  meta.data = summary_csv
)

head(pbmc2@meta.data)
dim(pbmc2)

# import mtDNA depths to filter cells
bci <- fread("../data/depth/Pearson-PBMC-BCI_v12-mtMask_mgatk.depthTable.txt"); bci[["V1"]] <- gsub("-1", "-1",  bci[["V1"]])
ccf <- fread("../data/depth/Pearson-PBMC-CCF_v12-mtMask_mgatk.depthTable.txt"); ccf[["V1"]] <- gsub("-1", "-2",  ccf[["V1"]])
h9 <- fread("../data/depth/PBMC_H9_v12-mtMask_mgatk.depthTable.txt"); h9[["V1"]] <- gsub("-1", "-3",  h9[["V1"]])
h10 <- fread("../data/depth/PBMC_H10_v12-mtMask_mgatk.depthTable.txt"); h10[["V1"]] <- gsub("-1", "-4",  h10[["V1"]])

# Append mito depth to Seurat object
depth_vec <- c(bci$V2, ccf$V2, h9$V2, h10$V2)
names(depth_vec) <- as.character(c(bci$V1, ccf$V1, h9$V1, h10$V1))
pbmc2@meta.data$mtDNA_depth <- depth_vec[rownames(pbmc2@meta.data)]
summary(pbmc2@meta.data$mtDNA_depth)
pbmc2$pct_reads_in_peaks <- pbmc2$peak_region_fragments / pbmc2$passed_filters * 100
qplot(pbmc2@meta.data$pct_reads_in_peaks, log10(pbmc2@meta.data$passed_filters))

# Filter
pbmc2 <- subset(pbmc2, subset = mtDNA_depth > 10 & passed_filters > 1000 & pct_reads_in_peaks)
dim(pbmc2)

# Add patient name too 
ptid <- substr(colnames(pbmc2), 18, 18)
pbmc2@meta.data$Patient_ID <- ptid
pbmc2@meta.data$Patient <- case_when(ptid == "1" ~ "BCI", ptid == "2" ~ "CCF", TRUE ~ "Healthy")
pbmc2@meta.data$Disease <- case_when(ptid %in% c("1", "2") ~ "Pearson", TRUE ~ "Healthy")

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
pbmc2 <- RunUMAP(pbmc2, dims = 2:20, reduction = 'harmony')
pbmc2 <- FindNeighbors(object = pbmc2, reduction = 'harmony', dims = 2:20) 
pbmc2 <-   FindClusters(object = pbmc2, verbose = FALSE, resolution = 0.5)
DimPlot(object = pbmc2, label = TRUE)
DimPlot(object = pbmc2, group.by = "Patient" )


# Look at motifs
DefaultAssay(pbmc) <- 'peaks'

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = FALSE)
)

#sapply(ms, name) %>% sort()
# create a Motif object and add it to the assay
motif.matrix <- CreateMotifMatrix(
  features = granges(pbmc2),
  pwm = pwm,
  genome = 'hg19',
  use.counts = FALSE
)

# Retain the actual motif matched positions, down to the base pair
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

Idents(pbmc2) <- "Patient"

pbmc2 <- Footprint(
  object = pbmc2,
  motif.name = c("MA1130.1", "MA1137.1"),
  genome = BSgenome.Hsapiens.UCSC.hg19,
  in.peaks = TRUE
)
Idents(pbmc2) <- "Disease"
PlotFootprint(pbmc2, features = c("MA1137.1"), )


pbmc2 <- AddMotifs(
  object = pbmc2,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  pfm = pwm
)
pbmc2 <- RunChromVAR(
  object = pbmc2,
  genome = BSgenome.Hsapiens.UCSC.hg19
)
df <- data.frame(D = ,)
colnames(df) <- c("Disease", unname(sapply(ms, name)))

rdf <- data.frame(
  row_wilcoxon_twosample( (pbmc2@assays$chromvar@data)["MA1130.1",pbmc2@meta.data$Disease== "Healthy"], (pbmc2@assays$chromvar@data)[,pbmc2@meta.data$Disease == "Pearson"]),
  motif = unname(sapply(ms, name)))

