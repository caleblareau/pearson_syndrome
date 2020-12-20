library(data.table)
library(dplyr)

# Import data
lane = "1"
lapply(as.character(1:8), function(lane){
  celldf <- fread(paste0("../data/for_infer_cnv/processed/PearsonInVitro_out",lane,"/cells.tsv"), header = FALSE)
  gp <- fread(paste0("../data/for_infer_cnv/processed/PearsonInVitro_out",lane,"/gene_positions.tsv"))
  pearson_mat <- fread(paste0("../data/for_infer_cnv/processed/PearsonInVitro_out",lane,"/inferCNVoutput/", "infercnv.21_denoised.observations.txt")) 
  healthy_mat <- fread(paste0("../data/for_infer_cnv/processed/PearsonInVitro_out",lane,"/inferCNVoutput/", "infercnv.21_denoised.references.txt"))
  
  # Score the cells 
  genes_chr7_del <- gp %>% filter(V2 == "7" & V3 > 100000000) %>% pull(V1)
  cm_pearson <- colMeans(data.matrix(data.frame(pearson_mat[pearson_mat$V1 %in% genes_chr7_del,-1])))
  cm_healthy <- colMeans(data.matrix(data.frame(healthy_mat[healthy_mat$V1 %in% genes_chr7_del,-1])))
  healthy01 <- quantile(cm_healthy, 0.01)
  message(healthy01)
  message(table(cm_pearson < healthy01))
  bdf <- data.frame(
    barcode = gsub("[.]", "-", c(names(cm_pearson), names(cm_healthy))),
    chr7_del_ascore = c(cm_pearson, cm_healthy),
    what = c(rep("Pearson", length(cm_pearson)), rep("Healthy", length(cm_healthy))),
    MDS = c(cm_pearson < healthy01, rep(NA, length(cm_healthy))), 
    lane = lane
  )
  bdf
}) %>% rbindlist() %>% data.frame() -> odf
write.table(odf, file = "../output/scRNA_MDS_annotations.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)