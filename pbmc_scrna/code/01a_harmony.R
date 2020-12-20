library(Seurat)
library(harmony)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
options(future.globals.maxSize = 4000 * 1024^2)

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
  raw <- CreateSeuratObject(counts = raw, project = "RNA")
  raw <- FindVariableFeatures(raw)
  raw <- NormalizeData(raw)
  raw <- ScaleData(raw)
  raw
}

# Import
pbmc1 <- import_scRNAseq("healthy_pbmc_8k_v2-remap", "H1", "Healthy")
pbmc2 <- import_scRNAseq("healthy_pbmc_4k_v2-remap", "H1", "Healthy")
pbmc3 <- import_scRNAseq("healthy_pbmc_5k_nextgem", "H2", "Healthy")
pbmc4 <- import_scRNAseq("healthy_pbmc_5k", "H2", "Healthy")

pBCI <- import_scRNAseq("pearson_bci", "pBCI", "Pearson")
pCCF <- import_scRNAseq("pearson_ccf", "pCCF", "Pearson")
pPT3 <- import_scRNAseq("pearson_mds", "pPT3", "Pearson")

merge_all <- merge(merge(pbmc1, pbmc2, add.cell.ids = c("h1", "h2")),
                 merge(merge(pBCI, pCCF, add.cell.ids = c("p1", "p2")), pPT3, add.cell.ids = c("", "p3")))
features <- SelectIntegrationFeatures(object.list = list(pBCI, pCCF, pPT3, pbmc1))
features <- pPT3@assays$RNA@var.features
merge_all <- ScaleData(merge_all, features = features, do.scale = FALSE)
merge_all <- RunPCA(merge_all, verbose = FALSE, features = features)
merge_all$cc <- substr(colnames(merge_all),1,2)
merge_all$cp <- substr(colnames(merge_all),1,3)

merge_all <- merge_all %>%  RunHarmony(c("cp"))

merge_all <- merge_all %>%  FindNeighbors(reduction = "harmony", dims = 1:30)  %>% RunUMAP(reduction = "harmony", dims = 1:30)
merge_all <- FindClusters(merge_all, resolution = 0.4) %>% identity()
DimPlot(merge_all, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, split.by = c('cp'), label = TRUE)

FindMarkers(merge_all, "11")
