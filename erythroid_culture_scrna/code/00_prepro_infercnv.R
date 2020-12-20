library(dplyr)
library(Seurat)
library(data.table)
library(Matrix)

soupercelldf <- fread("../data/soupor_cell_assignments.tsv")
soupercelldf$barcode_aggr <- paste0(substr(soupercelldf$barcode, 1, 16), "-", soupercelldf$lane)
ddvec <- c("6", "6", "6", "6", "12", "12", "12", "12")
ccvec <- as.character(1:8)

process_inferCNV <- function(idx = 1){
  
  lane_pro <- as.character(idx)
  soupercelldf_lane <- soupercelldf %>% filter(lane == lane_pro & (stat == "singlet"))
  cc <- ccvec[idx]
  dd <- ddvec[idx]
  dat <- Read10X_h5(paste0("../data/PearsonMixErythroid_Day",dd,"_Channel",cc,".h5"))
  dim(dat)
  dat <- dat[,(colSums(dat) <  75000) & (colSums(dat >= 1) >= 400) & (colnames(dat) %in% soupercelldf_lane$barcode)]
  dim(dat)
  soupercelldf_lane2 <- soupercelldf_lane %>% filter(barcode %in% colnames(dat)) 
  pos_df <- fread("../data/for_infer_cnv/infer_CNV_gene_annotations.tsv")
  
  dat2 <- dat[pos_df[[1]],soupercelldf_lane2$barcode]
  colnames(dat2) <- gsub("-", ".", colnames(dat2))
  
  # Export
  out_name <- paste0("PearsonInVitro_out", lane_pro)
  base_dir <- paste0("../data/for_infer_cnv/processed/", out_name)
  dir.create(base_dir)
  
  # Export cells file
  cell_outfile = paste0(base_dir, "/cells.tsv")
  soupercelldf_lane2$barcode <- gsub("-", ".", soupercelldf_lane2$barcode)
  soupercelldf_lane2[,c("barcode", "assignment")] %>% 
    fwrite(file = cell_outfile, col.names = FALSE, sep = "\t")
  
  # Export positions
  gene_outfile = paste0(base_dir, "/gene_positions.tsv")
  
  # Reorder chromosome names
  chrs_order <- paste0(c(as.character(1:22)))
  pos_df2 <- pos_df %>% filter(V1 %in% rownames(dat) & V2 %in% chrs_order) 
  pos_df2$V2 <- factor(pos_df2$V2, levels = chrs_order)
  
  pos_df3 <- pos_df2 %>% arrange(V2, V3)
  pos_df3 %>%
    fwrite(file = gene_outfile, col.names = FALSE, sep = "\t")
  
  # Export counts
  counts_outfile = paste0(base_dir, "/counts_mat.tsv")
  (dat)[pos_df3$V1,] %>% data.matrix() %>% data.frame() %>%
    fwrite(file = counts_outfile, col.names = TRUE, sep = "\t", row.names = TRUE)
  system(paste0("gzip ", counts_outfile))
  out_name
}

process_inferCNV(1)
process_inferCNV(2)
process_inferCNV(3)
process_inferCNV(4)
process_inferCNV(5)
process_inferCNV(6)
process_inferCNV(7)
process_inferCNV(8)

