library(data.table)
library(Matrix)
library(dplyr)
library(Signac)

#  Import raw things
pbmc <- readRDS("../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")
md <- pbmc@meta.data
rna <- pbmc@assays$RNA@data
edgeR_diff <- readRDS("../../pbmc_scrna/output/5March-PearsonRNAseq-diffGE-edgeR.rds")
gene_superset <- intersect(unique(edgeR_diff$gene), rownames(rna))

# Process heteroplasmy
df_BCI <- fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
  filter(deletion == "del6073-13095") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id))

df_CCF <- fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
  filter(deletion == "del8482-13447") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id))

df_PT3 <- fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
  filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id))

het_df <- rbind(df_BCI, df_CCF, df_PT3)
bdf <- merge(md, het_df, by = "barcode") %>% group_by(patient) %>% mutate(scaled_heteroplasmy = scale(heteroplasmy, scale = FALSE))

compute_gene_score_assoc <- function(test_celltype){
  
  cells_df <- bdf %>% filter(predicted.id == test_celltype) 
  rna_mat <- data.matrix(rna[gene_superset,cells_df$barcode])
  het_vec <- cells_df %>% pull(scaled_heteroplasmy)
  
  lapply(1:dim(rna_mat)[1], function(idx){ 
    data.frame(
      Gene = rownames(rna_mat)[idx],
      Celltype = test_celltype,
      data.frame(matrix(summary(lm(rna_mat[idx,]~het_vec))$coefficients[2,], nrow = 1))[,c(1,3,4)]
    ) -> df
    colnames(df) <- c("gene", "celltype","atacEstimate", "atac_T_stat", "atac_pval")
    df
  }) %>% rbindlist() %>% data.frame() -> odf
  
  adf <- odf %>% arrange(atac_pval)
  adf$atac_adj_p <- p.adjust(adf$atac_pval)
  adf
}

test_celltypes <- names(table(md$predicted.id))[table(md$predicted.id) > 200]
lapply(test_celltypes, function(test_celltype){
  print(test_celltype)
  compute_gene_score_assoc(test_celltype)
}) %>% rbindlist() %>% data.frame() %>% arrange(atac_pval) -> all_assocs_ATAC_df

all_assocs_ATAC_df %>% group_by(celltype) %>% summarize(sum(atac_adj_p < 0.05, na.rm = TRUE))

mdf <- merge(all_assocs_ATAC_df, edgeR_diff, by = c("gene", "celltype"))

mdf <- mdf %>% arrange(atac_pval)


library(BuenColors)
ggplot(mdf, aes(x = logFC, y = atacEstimate)) +
  geom_point()

cor(mdf$atacEstimate, mdf$logFC)
mdf %>% filter(atacEstimate > 0 & logFC > 0)
mdf %>% filter(atacEstimate < 0 & logFC < 0)

mdf %>% filter(gene == "EOMES")
