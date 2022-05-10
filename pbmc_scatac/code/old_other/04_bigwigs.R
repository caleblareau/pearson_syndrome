library(rtracklayer)
library(GenomicRanges)
library(Signac)
library(Seurat)

frags <- fread("../../../pearson_large_data_files/input/pbmcs_scatac/aggrs/P3H2_fragments.tsv.gz")
pbmc2 <- readRDS("../../../pearson_large_data_files/output/PBMC_scATAC_3P1H-15FEB2021.rds")
cdf <- pbmc2@meta.data
if(FALSE){
  mv <- FindMarkers(pbmc2, "Pearson", "Healthy", group.by = "Disease", max.cells.per.ident = 5000)
  write.table(mv, file = "../output/FindMarkers_Seurat_topATAC.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
donors <- unique(cdf$Patient)
sapply(donors[2:4], function(donor){
  print(donor)
  possible_ids <- cdf %>% dplyr::filter(Patient == donor) %>% pull(barcode.x) %>% as.character()
  
  cluster_gr <- frags %>% dplyr::filter(V4 %in% possible_ids) %>%
    setnames(c("chr", "start", "end", "V4", "PCRn")) %>%
    makeGRangesFromDataFrame()
  
  reads_coverage <- coverage(cluster_gr)/length(cluster_gr)*1000000
  export.bw(reads_coverage, con = paste0("../../../pearson_large_data_files/output/allPBMC_", as.character(donor), ".bw"))
  donor
}) -> bulk2



