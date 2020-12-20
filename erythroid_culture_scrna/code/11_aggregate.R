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

process_me <- function(idx){
  cc <- ccvec[idx]
  dd <- ddvec[idx]
  mat <- Read10X_h5(paste0("../data/PearsonMixErythroid_Day",dd,"_Channel",cc,".h5"))
  colnames(mat) <- gsub("-1", paste0("-", cc), colnames(mat))
  mat
}
full_mat <- do.call("cbind", lapply( 1:8, process_me))
subset_mat <- full_mat[,soupercelldf[["barcode_aggr"]]]

meta_df <- data.frame(
  row.names = soupercelldf$barcode_aggr,
  soupercelldf
)
so <- CreateSeuratObject(
  subset_mat,
  project = "invitro_erythroid",
  assay = "RNA",
  meta.data = meta_df
)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
so2 <- subset(so, subset = nFeature_RNA > 400 & nCount_RNA < 75000 & percent.mt < 25)
so2@meta.data %>% group_by(Day, assignment) %>% summarize(count = n())
so2@meta.data %>% group_by(Day, assignment) %>% summarize(count = mean(nFeature_RNA))
saveRDS(so2, file = "../big_data/init_seurat_object.rds")
