library(Matrix)
library(dplyr)
library(data.table)
library(GenomicRanges)

ga_d06_1 <- readRDS("../../../pearson_mtscatac_large_data_files/output/invitro_erythroid/D6_1_gene_activities.rds")
ga_d06_2 <- readRDS("../../../pearson_mtscatac_large_data_files/output/invitro_erythroid/D6_2_gene_activities.rds")
ga_d12_1 <- readRDS("../../../pearson_mtscatac_large_data_files/output/invitro_erythroid/D12_1_gene_activities.rds")
ga_d12_2 <- readRDS("../../../pearson_mtscatac_large_data_files/output/invitro_erythroid/D12_2_gene_activities.rds")

genes <- intersect(rownames(ga_d12_2), intersect(rownames(ga_d12_1), intersect(rownames(ga_d06_1), rownames(ga_d06_2))))

process_gs_lib <- function(lib, mat){
  cells_dt <- fread(paste0("../output/cell_assignments_per_channel/Pearson_",lib,"_assign.tsv"))
  control_cells <- cells_dt %>% dplyr::filter(assign == "Control") %>% pull(barcode)
  pearson_cells <- cells_dt %>% dplyr::filter(assign == "Pearson") %>% pull(barcode)
  
  control_mat <- mat[genes,control_cells]
  pearson_mat <- mat[genes,pearson_cells]
  
  cpmdf <- data.frame(
    c = round(rowSums(control_mat)/sum(control_mat)*1000000,1),
    p = round(rowSums(pearson_mat)/sum(pearson_mat)*1000000,1)
  )
  colnames(cpmdf) <- paste0(c("Control", "Pearson"), "_", lib)
  cpmdf
}

gene_df <- cbind(
  gene = genes,
  process_gs_lib("D6_1",ga_d06_1),
  process_gs_lib("D6_2",ga_d06_2),
  process_gs_lib("D12_1",ga_d12_1),
  process_gs_lib("D12_2",ga_d12_2)
)
gene_df$overall_mean <- round(rowMeans(data.matrix(gene_df)[,2:9]), 2)

gene_df$D06_logFC <- round(log2((gene_df$Pearson_D6_1 + gene_df$Pearson_D6_2 + 1)/(gene_df$Control_D6_1 + gene_df$Control_D6_2 + 1)),2)
gene_df$D12_logFC <- round(log2((gene_df$Pearson_D12_1 + gene_df$Pearson_D12_2 + 1)/(gene_df$Control_D12_1 + gene_df$Control_D12_2 + 1)),2)

write.table(gene_df %>% arrange(desc(D12_logFC)) , file = "../output/Pearson_invitroErythroid_diff_geneScores.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#---
# Make bigwigs
#---

d12_1_frags <- fread("../../../pearson_mtscatac_large_data_files/input/invitro_ery_diff/fragments/Pearson_Healthy_D12_1_v12-mtMask.fragments.tsv.gz")
d12_2_frags <- fread("../../../pearson_mtscatac_large_data_files/input/invitro_ery_diff/fragments/Pearson_Healthy_D12_2_v12-mtMask.fragments.tsv.gz")
d6_1_frags <- fread("../../../pearson_mtscatac_large_data_files/input/invitro_ery_diff/fragments/Pearson_Healthy_D6_1_v12-mtMask.fragments.tsv.gz")
d6_2_frags <- fread("../../../pearson_mtscatac_large_data_files/input/invitro_ery_diff/fragments/Pearson_Healthy_D6_2_v12-mtMask.fragments.tsv.gz")


process_frags <- function(lib1, lib2, frags1, frags2, day){
  
  # Process cells
  cells_dt1 <- fread(paste0("../output/cell_assignments_per_channel/Pearson_",lib1,"_assign.tsv"))
  control_cells1 <- cells_dt1 %>% dplyr::filter(assign == "Control") %>% pull(barcode)
  pearson_cells1 <- cells_dt1 %>% dplyr::filter(assign == "Pearson") %>% pull(barcode)
  
  cells_dt2 <- fread(paste0("../output/cell_assignments_per_channel/Pearson_",lib2,"_assign.tsv"))
  control_cells2 <- cells_dt2 %>% dplyr::filter(assign == "Control") %>% pull(barcode)
  pearson_cells2 <- cells_dt2 %>% dplyr::filter(assign == "Pearson") %>% pull(barcode)
  
  # Process fragments/ bigwigs
  frags_c_gr <- makeGRangesFromDataFrame(rbind(frags1[frags1$V4 %in% control_cells1,], frags2[frags2$V4 %in% control_cells2,]), 
                                         seqnames.field = "V1", start.field = "V2", end.field = "V3")
  covc <- coverage(frags_c_gr)/length(frags_c_gr)*1000000
  rtracklayer::export.bw(covc, con = paste0("../../../pearson_mtscatac_large_data_files/output/invitro_erythroid/invitro_pearson_bigwig/Control_",day, ".bw"))
  
  frags_p_gr <- makeGRangesFromDataFrame(rbind(frags1[frags1$V4 %in% pearson_cells1,], frags2[frags2$V4 %in% pearson_cells2,]), 
                                         seqnames.field = "V1", start.field = "V2", end.field = "V3")
  covp <- coverage(frags_p_gr)/length(frags_p_gr)*1000000
  rtracklayer::export.bw(covp, con = paste0("../../../pearson_mtscatac_large_data_files/output/invitro_erythroid/invitro_pearson_bigwig/Pearson_",day, ".bw"))
}

process_frags("D12_1", "D12_2", d12_1_frags, d12_2_frags, "Day12")
process_frags("D6_1", "D6_2", d6_1_frags, d6_2_frags, "Day06")
