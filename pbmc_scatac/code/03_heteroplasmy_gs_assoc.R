library(data.table)
library(Matrix)
library(dplyr)
library(Signac)

#  Import raw things
pbmc <- readRDS("../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")
md <- pbmc@meta.data
rna <- pbmc@assays$RNA@data
rna <- t(t(rna)/colSums(rna))*10000 # normalize
edgeR_diff <- readRDS("../../pbmc_scrna/output/5March-PearsonRNAseq-diffGE-edgeR.rds")
gene_superset <- intersect(unique(edgeR_diff$gene), rownames(rna))

# Process heteroplasmy
df_BCI <- fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del6073-13095") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id))

df_CCF <- fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del8482-13447") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id))

df_PT3 <- fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id))

het_df <- rbind(df_BCI, df_CCF, df_PT3)
bdf <- merge(md, het_df, by = "barcode") %>% group_by(patient) %>% mutate(scaled_heteroplasmy = scale(heteroplasmy, scale = FALSE))

# LOL
donors <- c("BCI", "CCF", "PT3")
compute_gene_score_assoc <- function(test_celltype){
  print(test_celltype)
  lapply(donors, function(pt_id_1){
    print(pt_id_1)
    cells_df <- bdf %>% dplyr::filter(predicted.id == test_celltype & Patient == pt_id_1) 
    ss_mat <- data.matrix(rna[,cells_df$barcode])
    het_vec <- cells_df %>% pull(heteroplasmy)
    
    cor_df <- data.frame(
      gene = rownames(ss_mat),
      celltype = test_celltype,
      patient = pt_id_1,
      cor_pearson = cor(t(ss_mat), het_vec,use = "pairwise.complete.obs"),
      cor_spearman = cor(t(ss_mat), het_vec,use = "pairwise.complete.obs", method = "spearman")
    )
    cor_df
  }) %>% rbindlist() -> adf
    adf
}

test_celltypes <- names(table(md$predicted.id))[table(md$predicted.id) > 75]
lapply(test_celltypes, function(test_celltype){
  print(test_celltype)
  compute_gene_score_assoc(test_celltype)
}) %>% rbindlist() %>% data.frame()  -> all_assocs_ATAC_df

all_assocs_ATAC_df %>% dplyr::filter(gene == "JUNB") %>% pull(cor_pearson) %>% summary()
all_assocs_ATAC_df %>% dplyr::filter(gene == "JUNB") %>% pull(cor_spearman) %>% summary()

all_assocs_ATAC_df %>% group_by(gene) %>% summarize(mc = mean(cor_pearson, na.rm = TRUE)) %>%
  arrange(desc(abs(mc))) %>% head(20)
