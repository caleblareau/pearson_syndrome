library(dplyr)
library(Seurat)
library(data.table)
library(Matrix)

import_data <- function(dir_base, qn){
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
  raw <- raw[!grepl("^MT", rownames(raw)), intersect(singlets, qc_cells)]
  colnames(raw) <- gsub("-1", paste0("_", qn), colnames(raw))
  raw
}

healthy <- import_data("healthy_pbmc_4k_v2-remap", "healthy")
pearson <- import_data("pearson_mds", "pearson")
dim(pearson)
dim(healthy)

pos_df <- fread("../data/forinfercnv/infer_CNV_gene_annotations.tsv")

# Export
out_name <- "PBMCMDS"
base_dir <- paste0("../data/forinfercnv/", out_name)
dir.create(base_dir)

# Export cells file
cell_outfile = paste0(base_dir, "/cells.tsv")
data.frame(cells = c(colnames(pearson), colnames(healthy)), 
           status = c(rep("pearson",dim(pearson)[2]), rep("healthy",dim(healthy)[2]))) %>% 
  fwrite(file = cell_outfile, col.names = FALSE, sep = "\t")

# Export positions
gene_outfile = paste0(base_dir, "/gene_positions.tsv")

# Reorder chromosome names
chrs_order <- paste0(c(as.character(1:22)))
pos_df2 <- pos_df %>% filter(V1 %in% rownames(pearson) & V2 %in% chrs_order) 
pos_df2$V2 <- factor(pos_df2$V2, levels = chrs_order)

pos_df3 <- pos_df2 %>% arrange(V2, V3)
pos_df3 %>%
  fwrite(file = gene_outfile, col.names = FALSE, sep = "\t")

# Export counts
counts_outfile = paste0(base_dir, "/counts_mat.tsv")
cbind(pearson,healthy)[pos_df3$V1,] %>% data.matrix() %>% data.frame() %>%
  fwrite(file = counts_outfile, col.names = TRUE, sep = "\t", row.names = TRUE)
system(paste0("gzip ", counts_outfile))
out_name
