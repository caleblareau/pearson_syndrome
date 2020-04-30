library(edgeR)
library(Seurat)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(stringr)

"%ni%" <- Negate("%in%")

# Function to do differential gene expression with edgeR w/ covariates
run_edgeRQLFDetRate_CL <- function(count, condt) {
  
  donor_id <- str_split_fixed(colnames(count), "-", 2)[,1]
  dge <- DGEList(count, group = condt)
  dge <- calcNormFactors(dge)
  
  # adjust for sequencing technology (H2 is 10x v3), total genes detected
  cdr <- scale(colMeans(count > 0))
  design <- model.matrix(~ cdr + as.numeric(donor_id == "H2") + condt) 
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit)
  tt <- topTags(qlf, n = Inf)
  df <- tt@.Data[[1]]
  df$gene <- rownames(df)
  
  # small function to pull CPM
  ps <- function(which_ones){
    rs <- rowSums(count[,which_ones])
    cpms <- round(rs/sum(rs) *1000000,1)[as.character(df$gene)]
    return(cpms)
  } 
  
  # Pull CPM values
  df$Pearson_cpm <- ps(condt == "p")
  df$Healthy_cpm <- ps(condt == "H")
  
  # Pull donors
  df$H1v2 <- ps(donor_id == "H1")
  df$H2v3 <- ps(donor_id == "H2")
  df$pBCI <- ps(donor_id == "pBCI")
  df$pCCF <- ps(donor_id == "pCCF")
  df$pPT3 <- ps(donor_id == "pPT3")
  
  # Round
  df$logFC <- round(df$logFC,2)
  df$logCPM <- round(df$logCPM,2)
  df$signedTstat <- sign(df$logFC) * sqrt(round(df$F,2))
  
  df
  
}


# Import
so <- readRDS("../output/5March-PearsonRNAseq-integration.rds")
counts <- so@assays$RNA@counts
df <- so@meta.data
df$barcode <- rownames(df)

cts <- unique(sort(as.character(df$celltype)))
boo_ct <- table((sort(as.character(df$celltype)))) > 250 & cts %ni% c( "low count 1", "CCF cell 1")
cts_go <- cts[boo_ct]
list_of_de_mats <- lapply(cts_go, function(ct){
  print(ct)
  
  # Pull out celltype specific barcodes
  bcs <- df %>% filter(celltype == ct) %>% pull(barcode)
  counts_ct <- counts[,bcs]
  
  # Pull out barcodes
  boo_gene <- log10(rowSums(counts_ct)/sum(counts_ct)*1000000) > 0.5
  counts_go <- counts_ct[boo_gene,]
  cond <- substr(colnames(counts_go), 1,1)
  
  dfo <- run_edgeRQLFDetRate_CL(counts_go, cond)
  dfo$celltype <- ct
  dfo
})
melted_list <- rbindlist(list_of_de_mats) %>% arrange(FDR)
melted_list %>% filter(FDR < 1e-4) %>% arrange(desc(abs(logFC)))
saveRDS(melted_list, file = "../output/5March-PearsonRNAseq-diffGE-edgeR.rds")

ggplot(melted_list %>%  filter(celltype ==  "NK-cells"), aes(x = logFC, y = -log10(FDR + 1e-100))) +
  geom_point()
