library(data.table)
library(Seurat)
library(Matrix)
library(dplyr)

# Import and munge
soupercelldf <- fread("../data/soupor_cell_assignments.tsv")
soupercelldf$barcode_aggr <- paste0(substr(soupercelldf$barcode, 1, 16), "-", soupercelldf$lane)
ddvec <- c("6", "6", "6", "6", "12", "12", "12", "12")
ccvec <- as.character(1:8)

# FIlter on the good stuf
soupercelldf <- soupercelldf %>% filter(assignment %in% c("Pearson", "Healthy") & stat == "singlet")
soupercelldf$Day <- ifelse(soupercelldf$lane %in% c("1", "2", "3", "4"), "D6", "D12")
meta_df <- data.frame(
  row.names = soupercelldf$barcode_aggr,
  soupercelldf
)

erythroid_gs <- fread("../data/erythroid_gene_set.txt", header = FALSE)[[1]]
process_me_sig <- function(idx){
  
  # Import Seurat matrix
  cc <- ccvec[idx]
  dd <- ddvec[idx]
  mat <- Read10X_h5(paste0("../data/PearsonMixErythroid_Day",dd,"_Channel",cc,".h5"))
  colnames(mat) <- gsub("-1", paste0("-", cc), colnames(mat))
  so <-  CreateSeuratObject(counts = mat[,colnames(mat) %in% soupercelldf[["barcode_aggr"]]],assay = "RNA", 
                            meta.data = meta_df)
  
  #QC
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so2 <- subset(so, subset = nFeature_RNA > 400 & nCount_RNA < 75000 & percent.mt < 25)
  so2 <- NormalizeData(so2)
  so2 <- AddModuleScore(so2, list("erythroid" = erythroid_gs), name = "erythroid_ms")   
  so2@meta.data$pct_erythroid <- colSums(so2@assays$RNA@counts[erythroid_gs,])/ colSums(so2@assays$RNA@counts)
  so2@meta.data$erythroid_ms1 <- ifelse(so2@meta.data$erythroid_ms1 > 10, 10, so2@meta.data$erythroid_ms1)
  
  so2@meta.data[,c("barcode_aggr", "erythroid_ms1","pct_erythroid")]
  
}
module_score_df <- lapply(1:8, process_me_sig) %>% rbindlist() %>% data.frame()
saveRDS(module_score_df, file = "../output/scRNA_erythroid_modulescore.rds")
