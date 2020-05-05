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

compute_dev_assoc <- function(test_celltype){
  
  cells_df <- bdf %>% filter(predicted.id == test_celltype) 
  ss_mat <- data.matrix(devmat[,cells_df$barcode])
  het_vec <- cells_df %>% pull(scaled_heteroplasmy)
  
  lapply(1:dim(ss_mat)[1], function(idx){ 
    data.frame(
      TF = rownames(ss_mat)[idx],
      Celltype = test_celltype,
      data.frame(matrix(summary(lm(ss_mat[idx,]~het_vec))$coefficients[2,], nrow = 1))[,c(1,3,4)]
    ) -> df
    colnames(df) <- c("TF", "celltype","cvEstimate", "cvT_stat", "cv_pval")
    df
  }) %>% rbindlist() %>% data.frame() -> odf
  
  adf <- odf %>% arrange(cv_pval)
  adf$cv_pval_p <- p.adjust(adf$cv_pval)
  adf
}

test_celltypes <- names(table(md$predicted.id))[table(md$predicted.id) > 200]
lapply(test_celltypes, function(test_celltype){
  print(test_celltype)
  compute_dev_assoc(test_celltype)
}) %>% rbindlist() %>% data.frame() %>% arrange(cv_pval) -> all_assocs_chromvar

all_assocs_chromvar %>% group_by(celltype) %>%
  summarize(sum(cv_pval_p < 0.01))

all_assocs_chromvar %>% group_by(celltype) %>%
  top_n(2, wt = -log10(cv_pval_p)) %>% data.frame()
