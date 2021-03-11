library(BuenColors)

libs <- c("pearson_bci", "pearson_ccf", "pearson_mds", "healthy_pbmc_8k_v2-remap", "healthy_pbmc_4k_v2-remap", "healthy_pbmc_5k_nextgem", "healthy_pbmc_5k")

dir_base = "pearson_mds"
df <- readRDS(paste0("../output/seurat_projected_meta_", dir_base, ".rds"))
ggplot(df, aes(x = refUMAP_1, y = refUMAP_2, color = predicted.celltype.l2)) +
  geom_point() 

import_project_scRNAseq <- function(dir_base){
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
  df <- readRDS(paste0("../output/seurat_projected_meta_", dir_base, ".rds"))
  
  raw[["umapref"]] <- CreateDimReducObject(embeddings = data.matrix(df[,c("refUMAP_1", "refUMAP_2")]), 
                                       loadings = matrix(),
                                       key = "UM", assay = "RNA")
  raw$celltype <- df$predicted.celltype.l2
  raw$CXCL14exp <- raw@assays$SCT@data["CXCL14",]
  raw$celltype_short <- case_when(
    raw$celltype == "CD4 Naive" ~ "CD4naive",
    raw$celltype == "CD8 Naive" ~ "CD8naive",
    raw$celltype == "CD4 TCM" ~ "CD4EM",
    raw$celltype == "CD4 TEM" ~ "CD4EM",
    raw$celltype == "CD8 TCM" ~ "CD8EM",
    raw$celltype == "CD8 TEM" ~ "CD8EM", 
    TRUE ~ "other"
  )
  
  
}
  