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
make_combined_df <- function(x){
  hetdf <- fread(paste0("../data/namdc-metadata/",x,"_del.new.tsv")) %>%
    dplyr::filter(version == "improved")
  
  merge(
    hetdf, 
    fread(paste0("../data/reference_projections/namdc/",x,"_refmapped.csv.gz")),
    by.y = "cb", by.x = "cell_id"
  ) %>%dplyr::filter(reads_all >= 5) -> df_int
  
  merge(
    df_int, 
    fread(paste0("../data/reference_projections/namdc/",x,"_refmapped.csv.gz")),
    by.y = "cb", by.x = "cell_id"
  ) %>%dplyr::filter(reads_all >= 5) -> df
  rownames(df) <- df$cell_id
  df
}

df671 <- make_combined_df("NAMDC20008038671")

# Now create the Seurat object
library(EnsDb.Hsapiens.v86)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "hg38"
seqlevelsStyle(annotations) <- 'UCSC'
annotations
annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))
annotations

# Set it up
mat <- Read10X_h5("../../../pearson_large_data_files/input/pbmcs_scatac/fragments/NAMDC20008038671_hg38_raw_peak_bc_matrix.h5")[,rownames(df671)]
fragments_file <-"../../../pearson_large_data_files/input/pbmcs_scatac/fragments/NAMDC20008038671_hg38_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = mat,
  sep = c(":", "-"),
  fragments = fragments_file,
  min.cells = 0,
  min.features = 0
)

pbmc2 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = df671
)

head(pbmc2@meta.data)
dim(pbmc2)


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
cd4t <- subset(pbmc2, predicted.celltype.l1.x == "CD4 T")
cormat <- cor(cd4t@meta.data$heteroplasmy,t(data.matrix(cd4t@assays$RNA@data)),
              use = "pairwise.complete")
perm_cormat <- cor(sample(cd4t@meta.data$heteroplasmy), t(data.matrix(cd4t@assays$RNA@data)),
                   use = "pairwise.complete")

obs_df <- data.frame(
  gene = colnames(cormat),
  cor = cormat[1,]
) %>% dplyr::filter(!is.na(cor)) %>% arrange(desc(cor)) %>% mutate(rank = 1:n(), what = "observed") 
head(obs_df, 20)
tail(obs_df, 20)

library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
pbmc2 <- pbmc2[stringr::str_split_fixed(rownames(pbmc2),"-", 2)[,1] %in% paste0("chr", 1:22),]

pbmc2 <- AddMotifs(
  object = pbmc2,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
cv <- RunChromVAR(pbmc2)
