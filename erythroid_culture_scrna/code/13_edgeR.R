library(edgeR)
library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(stringr)

so <- readRDS("../../../pearson_large_data_files/output/invitro_erythroid/Invitro_erythroid_scRNA_seurat_object.rds")
so <- subset(so, module_score > 0.25)

"%ni%" <- Negate("%in%")

# Function to do differential gene expression with edgeR w/ covariates
run_edgeRQLFDetRate_CL_ery <- function(count, condt) {
  
  donor_id <- str_split_fixed(colnames(count), "_", 2)[,1]
  dge <- DGEList(count, group = condt)
  dge <- calcNormFactors(dge)
  
  # Adjust for total genes detected
  cdr <- scale(colMeans(count > 0))
  design <- model.matrix(~ cdr + condt) 
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit)
  tt <- topTags(qlf, n = Inf)
  df <- signif(tt@.Data[[1]], 3)
  df$gene <- rownames(df)
  
  # small function to pull CPM
  ps <- function(which_ones){
    rs <- rowSums(count[,which_ones])
    cpms <- round(rs/sum(rs) *1000000,1)[as.character(df$gene)]
    return(cpms)
  } 
  
  # Pull CPM values
  df$Pearson_cpm <- ps(condt == "Pearson")
  df$Healthy_cpm <- ps(condt == "Healthy")
  
  # Round
  df$logFC <- round(df$logFC,2)
  df$logCPM <- round(df$logCPM,2)
  df$signedTstat <- sign(df$logFC) * sqrt(round(df$F,2))
  df
  
}

# Filter out xy genes
gene_coords <- fread("../../pbmc_scrna/data/forinfercnv/infer_CNV_gene_annotations.tsv")
xy_genes <- gene_coords %>% filter(V2 %in% c("X", "Y")) %>% pull(V1)
counts <- (so@assays$RNA@counts)[!grepl("^RP|^MT-", rownames(so@assays$RNA@counts)),]
counts <- counts[!(rownames(counts) %in% xy_genes),]
meta_df <- so@meta.data

# Filer some more genes
rs <- rowSums(counts)
cpms <- round(rs/sum(rs) *1000000,1)
counts <- counts[cpms > 0.1,]
dim(counts)
ery_df6 <- run_edgeRQLFDetRate_CL_ery(counts[,meta_df$Day == "D6"],  meta_df$assignment[meta_df$Day == "D6"])
ery_df12 <- run_edgeRQLFDetRate_CL_ery(counts[,meta_df$Day == "D12"], meta_df$assignment[meta_df$Day == "D12"])

save(ery_df6, ery_df12, file = "../output/Erythroid_Person_edgeR_D6D12-DGE.rda")



