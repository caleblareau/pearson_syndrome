library(Seurat)
library(Signac)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BuenColors)
library(motifmatchr)
library(chromVARmotifs)
data("human_pwms_v2") #filtered collection of human motifs from cisBP database

source("../functions/import_counts_mat_10x.R")
set.seed(1)

# Import ATAC data
counts <- import_counts_mat_10x("../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_filtered_peak_bc_matrix/")
pbmc <- readRDS("../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")

# Import mtDNA data
df_BCI <- fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
  filter(deletion == "del6073-13095") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id))

df_CCF <- fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
  filter(deletion == "del8482-13447") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id))

df_PT3 <- fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
  filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id))

het_df <- rbind(df_BCI, df_CCF, df_PT3)
md <- pbmc@meta.data
bdf <- merge(md, het_df, by = "barcode") %>% group_by(patient) %>% mutate(scaled_heteroplasmy = scale(heteroplasmy, scale = FALSE))

# Format counts data for chromVAR
counts <- counts[,rownames(md)]
mdf <- data.frame(matrix(unlist(strsplit(rownames(counts), ":|-")), ncol = 3, byrow = TRUE))
mdf$chr <- as.character(mdf[,1])
mdf$start <- as.numeric(as.character(mdf[,2]))
mdf$end <- as.numeric(as.character(mdf[,3]))
gr <- makeGRangesFromDataFrame(mdf)

# Make a summarized experiment for CV
SE <- SummarizedExperiment(
  assays = list(counts = counts), 
  rowRanges = gr,
  colData = pbmc@meta.data
)

# Run chromVAR
SE <- filterPeaks(SE)
mm <- motifmatchr::matchMotifs(human_pwms_v2, SE, genome = BSgenome.Hsapiens.UCSC.hg19)
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
dev <- computeDeviations(object = SE,  annotations = mm)
devmat <- assays(dev)[["deviations"]][,bdf$barcode]
row.names(devmat) <- make.unique(unname(rowData(dev)$name))

# Look at the association between each TF and the 
compute_dev_assoc <- function(test_celltype){
  
  cells_df <- bdf %>% filter(predicted.id == test_celltype) 
  ss_mat <- data.matrix(devmat[,cells_df$barcode])
  het_vec <- cells_df %>% pull(scaled_heteroplasmy)
  
  lapply(1:dim(ss_mat)[1], function(idx){ 
    data.frame(
      TF = rownames(ss_mat)[idx],
      Celltype = test_celltype,
      cor_pearson = cor(ss_mat[idx,], het_vec,use = "pairwise.complete.obs"),
      cor_spearman = cor(ss_mat[idx,], het_vec,use = "pairwise.complete.obs", method = "spearman"),
      data.frame(matrix(summary(lm(ss_mat[idx,]~het_vec))$coefficients[2,], nrow = 1))[,c(1,3,4)]
    ) -> df
    colnames(df) <- c("gene", "celltype","cor_pearson", "cor_spearman","cvEstimate", "cvT_stat", "cv_pval")
    df
  }) %>% rbindlist() %>% data.frame() -> odf
  
  adf <- odf %>% arrange(cv_pval)
  adf$cv_pval_p <- p.adjust(adf$cv_pval)
  adf
}

# Compute per cell type
test_celltypes <- names(table(md$predicted.id))[table(md$predicted.id) > 200]
lapply(test_celltypes, function(test_celltype){
  print(test_celltype)
  compute_dev_assoc(test_celltype)
}) %>% rbindlist() %>% data.frame() %>% arrange(cv_pval) -> all_assocs_chromvar

saveRDS(all_assocs_chromvar, file = "../../../pearson_mtscatac_large_data_files/output/PBMC_scATAC_chromVAR_cors.rds")

# Import diff RNA-seq
rna_df <- readRDS("../../pbmc_scrna/output/5March-PearsonRNAseq-diffGE-edgeR.rds")
merge_go_df <- merge(rna_df, all_assocs_chromvar, by = c("celltype", "gene"))

"%ni%" <- Negate("%in%")
ggplot(merge_go_df, aes(x = logFC, y = cor_pearson, color = gene %in% c("FOS", "FOSB", "JUN", "JUNB"))) +
  geom_point() + facet_wrap(~celltype) + labs(x = "log FC WT / Pearson", y = "Deviation score correlation with heteroplasmy") +
  pretty_plot() + theme(legend.position = "bottom") + scale_color_manual(values = c("black", "red"))

merge_go_df %>% filter(gene %ni% c("FOS", "FOSB", "JUN", "JUNB")) %>% filter(logFC < -1) %>%
  filter(cor_spearman < -0.1)
merge_go_df %>% filter(gene %in% "ATF4")
