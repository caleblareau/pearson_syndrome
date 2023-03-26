library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)

# Import gene expression and set up object
m1 <- Read10X_h5("../data/PT1_rep1_mtDOGMA_gex.h5")
m2 <- Read10X_h5("../data/PT1_rep2_mtDOGMA_gex.h5")
colnames(m2) <- gsub("-1", "-2", colnames(m2))
full_mat <- cbind(m1, m2)
meta <- fread("../output/Pearson_DOGMA_full_meta_data.tsv") %>% data.frame()
m1 <- read.table("../data/qc/DOGMA_PT1_Rep1_per_barcode_metrics.csv.gz", sep = ",", header = TRUE)
m2 <- read.table("../data/qc/DOGMA_PT1_Rep2_per_barcode_metrics.csv.gz", sep = ",", header = TRUE)
rownames(m1) <- m1$barcode
rownames(m2) <- gsub("-1", "-2", m2$barcode)
meta2 <- rbind(m1, m2)
rownames(meta) <- meta$cb

# only keep "cells"
both_cells <- intersect(meta$cb, colnames(full_mat))
meta2d <- meta2[rownames(meta2) %in% both_cells,c("atac_mitochondrial_reads", "atac_raw_reads")]
meta2d$pct_mito_atac <- meta2d$atac_mitochondrial_reads/meta2d$atac_raw_reads*100

so <- CreateSeuratObject(
  full_mat[,both_cells], meta.data = cbind(meta[both_cells,], meta2d[both_cells,])
)
so[["pct_mito_rna"]] <- PercentageFeatureSet(so, pattern = "^MT-")

library(BuenColors)
p1 <- so@meta.data[,c("pct_mito_rna", "pct_mito_atac")] %>% reshape2::melt() %>%
  ggplot(aes(x = variable, y = value)) + 
  geom_boxplot(outlier.shape = NA, width = 0.7) + 
  pretty_plot(fontsize = 7) + L_border() +
  coord_cartesian(ylim = c(0, 40)) +
  labs(x = "modality", y = "% mitochondrial molecules") 
cowplot::ggsave2(p1, file = "../plots/modality_abundance.pdf", width = 1.5, height = 1.5)
