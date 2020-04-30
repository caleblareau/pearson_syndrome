library(data.table)
library(dplyr)
library(Seurat)
library(Matrix)
library(pheatmap)
library(matrixStats)

"%ni%" <- Negate("%in%")

# Pull data from previous steps
so <- readRDS("../output/5March-PearsonRNAseq-integration.rds")
counts <- so@assays$RNA@counts
df <- so@meta.data
df$barcode <- rownames(df)

# Define celltypes
cts <- unique(sort(as.character(df$celltype)))
boo_ct <- table((sort(as.character(df$celltype)))) > 250 & cts %ni% c( "low count 1", "CCF cell 1")
cts_go <- cts[boo_ct]
pts <- unique(df$source)

mat_list <- lapply(cts_go, function(ct){
  print(ct)
  sapply(pts, function(pt){
    
    # Pull out celltype specific barcodes
    bcs <- df %>% filter(celltype == ct & source == pt) %>% pull(barcode)
    counts_ct <- counts[,bcs]
    
    # Pull out barcodes
    round(rowSums(counts_ct)/sum(counts_ct)*1000000,1)
    
  }) -> mm
  colnames(mm) <- paste0(colnames(mm), "_",ct)
  mm
})

bm <- do.call("cbind", mat_list)

pdf("../plots/cluster_views.pdf", width = 5, height = 5)
lapply(1:13, function(idx){
  (pheatmap(cor(log2(bm[so@assays$integrated@var.features,1:5 + (idx-1)*5]+1))))
})
dev.off()