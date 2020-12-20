library(Seurat)
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
library(reticulate)

options(future.globals.maxSize = 4000 * 1024^2)
use_condaenv("base")
scanorama <- import('scanorama')

# Import
import_scRNAseq <- function(dir_base, name, pheno){
  data.dir <- paste0("../data/", dir_base)
  raw <- Read10X(data.dir = data.dir)
  colnames(raw) <- paste0(substr(colnames(raw), 1, 16), "-1")
  
  # import scrublet results
  singlets <- fread(paste0("../data/scrublet_out/", dir_base, ".scrub.tsv")) %>%
    data.frame() %>% dplyr::filter(!called) %>% pull(barcode)
  
  # Remove crazy high and low expressors
  n_feature_rna <- colSums(raw > 0)
  n_total_rna <- colSums(raw)
  pct_mito <- colSums(raw[grepl("^MT", rownames(raw)), ])/n_total_rna * 100
  qc_cells <- colnames(raw)[pct_mito < 10 & n_total_rna > 1000 & n_feature_rna > 500]
  
  # Filter for singlet cells and non-mito genes
  raw <- raw[!grepl("^MT", rownames(raw)), intersect(singlets, qc_cells)]
  raw <- CreateSeuratObject(counts = counts, project = "RNA")
  raw <- FindVariableFeatures(raw)
  raw <- NormalizeData(raw)
  raw <- ScaleData(raw)
}


# Import
pbmc1 <- import_scRNAseq("healthy_pbmc_8k_v2-remap", "H1", "Healthy")
pbmc2 <- import_scRNAseq("healthy_pbmc_4k_v2-remap", "H1", "Healthy")
pbmc3 <- import_scRNAseq("healthy_pbmc_5k_nextgem", "H2", "Healthy")
pbmc4 <- import_scRNAseq("healthy_pbmc_5k", "H2", "Healthy")

pBCI <- import_scRNAseq("pearson_bci", "pBCI", "Pearson")
pCCF <- import_scRNAseq("pearson_ccf", "pCCF", "Pearson")
pPT3 <- import_scRNAseq("pearson_mds", "pPT3", "Pearson")
features <- SelectIntegrationFeatures(object.list = list((pBCI), (pCCF), (pPT3)))
genes_list <- list(as.list(features),as.list(features),as.list(features))

# Integration.
integrated.data <- scanorama$integrate(list(t(data.matrix(pBCI)), t(data.matrix(pCCF)), t(data.matrix(pPT3))), genes_list)
so <- CreateSeuratObject(cbind(pBCI, pCCF, pPT3))
so[["sr"]] <- CreateDimReducObject(embeddings = do.call("rbind", integrated.data[[1]]), key = "sr_", assay = "RNA")
so <- RunUMAP(so, dims = 1:30,  reduction = "sr")
ggplot(data.frame(so@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point()

# Batch correction.
corrected.data <- scanorama$correct(datasets, genes_list, return_dense=TRUE)

# Integration and batch correction.
integrated.corrected.data <- scanorama$correct(list(data.matrix(pBCI), data.matrix(pCCF), data.matrix(pPT3)), genes_list,
                                               return_dimred=TRUE, return_dense=TRUE)

list(pBCI, pCCF, pPT3, pbmc1, pbmc2, pbmc3, pbmc4)