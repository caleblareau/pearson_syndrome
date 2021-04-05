library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(BuenColors)
library(magrittr)
library(Matrix)
library(CONICSmat)
library(data.table)

set.seed(1)

ddvec <- c("6", "6", "6", "6", "12", "12", "12", "12")
ccvec <- as.character(1:8)
mat_test <- Read10X_h5("../data/PearsonMixErythroid_Day12_Channel5.h5")
regions <- data.frame(
  Chrom = c(7,7,7),
  Start = c(0, 110000000,61528020),
  End = c(58169653, 159345973,159345973),
  Length = c(58169653, 499345973,97817953),
  row.names = c("7p","X7del", "7q")
)
gene_pos=getGenePositions(rownames(mat_test))

# Run through them all
process_chr7_del<- function(idx){
  
  output_file <- paste0("../output/conics_mat/Conics_", idx)
  cc <- ccvec[idx]
  dd <- ddvec[idx]
  mat <- Read10X_h5(paste0("../data/PearsonMixErythroid_Day",dd,"_Channel",cc,".h5"))
  
  soupercelldf <- fread("../data/soupor_cell_assignments.tsv")
  soupercelldf$barcode_aggr <- paste0(substr(soupercelldf$barcode, 1, 16), "-", soupercelldf$lane)
  ps_cells <- soupercelldf %>% dplyr::filter(assignment == "Pearson" & stat == "singlet" & lane == idx) %>% dplyr::pull(barcode)
  healthy_cells <- soupercelldf %>% dplyr::filter(assignment == "Healthy" & stat == "singlet" & lane == idx)%>% dplyr::pull(barcode)
  
  ps_mat <- mat[,colnames(mat) %in% ps_cells & colSums(mat > 0) > 500]
  healthy_mat <- mat[,colnames(mat) %in% healthy_cells & colSums(mat > 0) > 500]
  
  # Run CNVer for pearson
  ps_CPM <- t(t(ps_mat) /colSums(ps_mat)*1000000)
  ps_go_expr <- log2(ps_CPM/10+1)
  ps_go_expr=filterMatrix(ps_go_expr,gene_pos[,"hgnc_symbol"],minCells=5)
  ps_normFactor=calcNormFactors(ps_go_expr)
  ps_l=plotAll(ps_go_expr,ps_normFactor,regions,gene_pos, paste0(output_file, "_PS_CNVplots"))
  chr7_total_ps <- data.frame(round(ps_l,3), barcode = rownames(ps_l), ngenes = colSums(ps_mat > 0))
  write.table(chr7_total_ps, file = paste0(output_file, "_Pearson_stats.tsv"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  # Run CNVer for healthy
  healthy_CPM <- t(t(healthy_mat) /colSums(healthy_mat)*1000000)
  healthy_go_expr <- log2(healthy_CPM/10+1)
  healthy_go_expr=filterMatrix(healthy_go_expr,gene_pos[,"hgnc_symbol"],minCells=5)
  healthy_normFactor=calcNormFactors(healthy_go_expr)
  healthy_l=plotAll(healthy_go_expr,healthy_normFactor,regions,gene_pos, paste0(output_file, "_Healthy_CNVplots"))
  chr7_total_healthy <- data.frame(round(healthy_l,3), barcode = rownames(healthy_l), ngenes = colSums(healthy_mat > 0))
  write.table(chr7_total_healthy, file = paste0(output_file, "_Healthy_stats.tsv"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  # Run CNVer for pearson
  both_mat <- cbind(ps_mat, healthy_mat)
  bo_CPM <- t(t(both_mat) /colSums(both_mat)*1000000)
  bo_go_expr <- log2(bo_CPM/10+1)
  bo_go_expr=filterMatrix(bo_go_expr,gene_pos[,"hgnc_symbol"],minCells=5)
  bo_normFactor=calcNormFactors(bo_go_expr)
  bo_l=plotAll(bo_go_expr, bo_normFactor,regions,gene_pos, paste0(output_file, "_both_CNVplots"))
  chr7_total_bo <- data.frame(round(bo_l,3), barcode = rownames(bo_l), ngenes = colSums(both_mat > 0),
                              who = c(rep("Pearson", dim(ps_mat)[2]), rep("Healthy", dim(healthy_mat)[2])))
  write.table(chr7_total_bo, file = paste0(output_file, "_both_stats.tsv"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

lapply(1:8, process_chr7_del)

both8 <-fread("../output/conics_mat/Conics_8_both_stats.tsv")

ggplot(both8, aes(x = X7del, fill = who)) + 
  geom_histogram() + facet_wrap(~who)
