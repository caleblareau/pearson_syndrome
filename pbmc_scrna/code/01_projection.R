library(Seurat)
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
library(SeuratDisk)

options(future.globals.maxSize = 4000 * 1024^2)
reference <- LoadH5Seurat("../../../pearson_large_data_files/input/pbmc/pbmc_multimodal.h5seurat")

# Import
import_project_scRNAseq <- function(dir_base, name, pheno){
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
  length(qc_cells)
  
  # Filter for singlet cells and non-mito genes
  raw <- raw[!grepl("^MT", rownames(raw)), intersect(singlets, qc_cells)]
  raw <- CreateSeuratObject(counts = raw, project = "RNA")
  raw <- SCTransform(raw)
  anchors <- FindTransferAnchors(
    reference = reference,
    query = raw,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )
  projected <- MapQuery(
    anchorset = anchors,
    query = raw,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
  )
  df <- data.frame(
    projected@meta.data,
    projected@reductions$ref.umap@cell.embeddings)
  df$name <- name
  df$pheno <- pheno
  saveRDS(df, file = paste0("../output/seurat_projected_meta_", dir_base, ".rds"))
}


ph1 <- import_project_scRNAseq("pediatrichealthy_pbmc_31687", "P1", "Healthy")
ph2 <- import_project_scRNAseq("pediatrichealthy_pbmc_31697", "P2", "Healthy")

pBCI <- import_project_scRNAseq("pearson_bci", "pBCI", "Pearson")
pCCF <- import_project_scRNAseq("pearson_ccf", "pCCF", "Pearson")
pPT3 <- import_project_scRNAseq("pearson_mds", "pPT3", "Pearson")

pbmc1 <- import_project_scRNAseq("healthy_pbmc_8k_v2-remap", "H1", "Healthy")
pbmc2 <- import_project_scRNAseq("healthy_pbmc_4k_v2-remap", "H1", "Healthy")
pbmc3 <- import_project_scRNAseq("healthy_pbmc_5k_nextgem", "H2", "Healthy")
pbmc4 <- import_project_scRNAseq("healthy_pbmc_5k", "H2", "Healthy")

# Enrich me
lapply(list.files("../output", full.names = TRUE, pattern = "*seurat*"), readRDS) %>% rbindlist() %>%
  group_by(predicted.celltype.l2, name) %>% summarize(count = n()) %>%
  ungroup() %>% group_by(name) %>% mutate(prop = count/ sum(count) * 100) %>% data.frame()

