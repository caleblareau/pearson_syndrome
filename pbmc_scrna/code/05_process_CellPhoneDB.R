library(data.table)
library(Seurat)
library(dplyr)

so <- readRDS("../../../pearson_mtscatac_large_data_files/output/5March-PearsonRNAseq-integration.rds")
head(so@meta.data)
keep_celltypes <- c("CD14+ Monocytes","Naive CD4+ T-cells","Memory CD4+ T-cells","NK-like T-cells","Effector CD8+ T-cells",
                    "Dendritic cells","Naive CD8+ T-cells","CD16+ Monocytes","Activated B-cells","Gamma-Delta T-cell","pDC",
                    "Memory B-cells","NK-cells","CD56 bright NK-cell")
donors <- unique(so@meta.data$source)
meta2 <- data.frame(Cell = rownames(so@meta.data),
                    so@meta.data[,c('celltype', 'source')])

lapply(donors[2:5], function(pt){
  ss <- meta2 %>% filter(source == pt & celltype %in% keep_celltypes)
  
  # Process raw data
  dat <- data.frame(Gene = rownames(so@assays$RNA@data), 
                    data.matrix(so@assays$RNA@data[,ss$Cell]))
  write.table(dat, file = paste0("../../pbmc_scrna/output/cellphonedb_data/", pt, "_counts.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # Process meta data
  out_meta <- ss[,c("Cell", "celltype")]
  out_meta$Cell <- gsub("-", ".",out_meta$Cell)
  colnames(out_meta) <- c("Cell", "cell_type")
  write.table(out_meta, file = paste0("../../pbmc_scrna/output/cellphonedb_data/", pt, "_meta.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  pt
})