library(data.table)
library(dplyr)
library(Matrix)
library(edgeR)
library(BuenColors)

# import RNA-seq data
mat <- fread("../data/rnaseq1_12July2019/matrix.mtx.gz", skip = 3) %>% data.matrix() 
sm <- sparseMatrix(i = mat[,1], j = mat[,2], x = mat[,3])
colnames(sm) <- as.character(read.table("../data/rnaseq1_12July2019/barcodes.tsv")[[1]])
rownames(sm) <- make.unique(as.character(read.table("../data/rnaseq1_12July2019/features.tsv")[[2]])[1:dim(sm)[1]])

# import donors and separate
ddf <- read.table("../output/12July_5pRNA_assignments.tsv", header = TRUE)
C1 <- sm[,ddf %>% filter(classify == "C1") %>% pull(barcode) ]
P1 <- sm[,ddf %>% filter(classify == "P1") %>% pull(barcode)]
P2 <- sm[,ddf %>% filter(classify == "P2") %>% pull(barcode)]

# Filter for CPM > 0.5
counts_all <- Matrix::rowSums(sm)
counts_cpm <- ((counts_all / sum(counts_all) ) * 1000000)
boo_all <- counts_all > 5

P1cpm_vec <- Matrix::rowSums(P1[boo_all,]) / sum(P1[boo_all,])  * 1000000
P2cpm_vec <- Matrix::rowSums(P2[boo_all,]) / sum(P2[boo_all,])  * 1000000
C1cpm_vec <- Matrix::rowSums(C1[boo_all,]) / sum(C1[boo_all,])  * 1000000


# Function to do differential gene expression with edgeR w/ covariates
run_edgeRQLFDetRate_CL <- function(A, B) {
  condt <- factor(c(rep("A", dim(A)[2]), rep("B", dim(B)[2])))
  count <- data.matrix(cbind(A[boo_all,], B[boo_all,]))
  dge <- DGEList(count, group = condt)
  dge <- calcNormFactors(dge)
  
  cdr <- scale(colMeans(count > 0))
  design <- model.matrix(~ cdr  + condt) 
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit)
  tt <- topTags(qlf, n = Inf)
  df <- tt@.Data[[1]]
  df$gene <- rownames(df)
  
  # Pull donors
  df$P1cpm <-  P1cpm_vec[rownames(df)]
  df$P2cpm <-  P2cpm_vec[rownames(df)]
  df$C1cpm <-  C1cpm_vec[rownames(df)]

  # Round
  df$logFC <- round(df$logFC,2)
  df$logCPM <- round(df$logCPM,2)
  df$signedTstat <- sign(df$logFC) * sqrt(round(df$F,2))
  
  df
  
}

results_P1_C1 <- run_edgeRQLFDetRate_CL(C1, P1)
results_P2_C1 <- run_edgeRQLFDetRate_CL(C1, P2)

saveRDS(results_P1_C1, file = "../output/results_P1_C1.rds")
saveRDS(results_P2_C1, file = "../output/results_P2_C1.rds")
