library(Seurat)
library(dplyr)

# import data
healthymd <- readRDS("../../../pearson_mtscatac_large_data_files/output/Healthy_scATAC_labelTransfer_HQcells.rds")@meta.data
pearsonmd <- readRDS("../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")@meta.data
frags_p <- fread("../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_fragments.tsv.gz")
frags_h1 <- fread("../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/PBMC_H9_v12-mtMask_fragments.tsv.gz")
frags_h2 <- fread("../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/PBMC_H10_v12-mtMask_fragments.tsv.gz")

frags_h1$V4 <- gsub("-1", "-4", frags_h1$V4)
frags_h2$V4 <- gsub("-1", "-5", frags_h2$V4)
frags_h <- rbind(frags_h1, frags_h2)

pearsonmd$predicted.id <- gsub("Activated B-cells", "Naive B-cells", as.character(pearsonmd$predicted.id))
tb <- table(pearsonmd$predicted.id)
celltypes <- names(tb)[tb >75 ]

table(pearsonmd$Patient)
lapply(celltypes, function(ct){
  
  ct2 <- gsub("[+]", "", gsub(" ", "_", ct))
  
  barcodes_h <- rownames(healthymd)[healthymd$predicted.id == ct] 
  barcodes_b <- pearsonmd %>% dplyr::filter(predicted.id == ct & Patient == "BCI") %>% pull(barcode)
  barcodes_c <- pearsonmd %>% dplyr::filter(predicted.id == ct & Patient == "CCF") %>% pull(barcode)
  barcodes_p <- pearsonmd %>% dplyr::filter(predicted.id == ct & Patient == "PT3") %>% pull(barcode)
  
  make_bw <- function(frags, bc, donor){
    cluster_gr <- frags %>% dplyr::filter(V4 %in% bc) %>%
      setnames(c("chr", "start", "end", "V4", "PCRn")) %>%
      makeGRangesFromDataFrame()
    
    reads_coverage <- coverage(cluster_gr)/length(cluster_gr)*1000000
    export.bw(reads_coverage, con = paste0("../../../pearson_mtscatac_large_data_files/output/pbmcs_bigwigs/PBMC_",donor, "-", as.character(ct2), ".bw"))
    donor
  }
  
  make_bw(frags_h, barcodes_h, "healthy")
  make_bw(frags_p, barcodes_b, "BCI")
  make_bw(frags_p, barcodes_c, "CCF")
  make_bw(frags_p, barcodes_p, "PT3")
  make_bw(frags_p, c(barcodes_p, barcodes_c, barcodes_c), "pearson")
  ct2
}) 



