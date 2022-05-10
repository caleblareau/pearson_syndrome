library(data.table)
library(Seurat)
library(Matrix)
library(dplyr)

# Import and munge
soupercelldf <- fread("../data/soupor_cell_assignments.tsv")
soupercelldf$barcode_aggr <- paste0(substr(soupercelldf$barcode, 1, 16), "-", soupercelldf$lane)
ddvec <- c("6", "6", "6", "6", "12", "12", "12", "12")
ccvec <- as.character(1:8)

# Filter on the good stuff
soupercelldf <- soupercelldf %>% dplyr::filter(assignment %in% c("Pearson", "Healthy") & stat == "singlet")
soupercelldf$Day <- ifelse(soupercelldf$lane %in% c("1", "2", "3", "4"), "D6", "D12")


if(FALSE){
  write.table(soupercelldf, file = "../output/Pearson_erythroid_scRNA_soupercell_assignment.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

process_me <- function(idx){
  cc <- ccvec[idx]
  dd <- ddvec[idx]
  mat <- Read10X_h5(paste0("../data/PearsonMixErythroid_Day",dd,"_Channel",cc,".h5"))
  colnames(mat) <- gsub("-1", paste0("-", cc), colnames(mat))
  mat
}
full_mat <- do.call("cbind", lapply( 1:8, process_me))
subset_mat <- full_mat[,soupercelldf[["barcode_aggr"]]]

meta_df_basic  <- data.frame(
  row.names = soupercelldf$barcode_aggr,
  soupercelldf
)

chr7 <- fread("../output/scRNA_MDS_annotations.tsv")
mds_vec <- chr7$MDS; names(mds_vec) <- paste0(substr(chr7$barcode, 1, 17), chr7$lane)
module <- readRDS("../output/scRNA_erythroid_modulescore.rds")
module_vec <- module$erythroid_ms1; names(module_vec) <- module[["barcode_aggr"]]

so <- CreateSeuratObject(
  subset_mat,
  project = "invitro_erythroid",
  assay = "RNA",
  meta.data = meta_df_basic
)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
so2 <- subset(so, subset = nFeature_RNA > 400 & nCount_RNA < 75000 & percent.mt < 25)
dim(so2)
so2@meta.data$MDS <- mds_vec[rownames(so2@meta.data)]
so2@meta.data$module_score <- module_vec[rownames(so2@meta.data)]

saveRDS(so2, file = "../../../pearson_large_data_files/output/invitro_erythroid/Invitro_erythroid_scRNA_seurat_object.rds")
