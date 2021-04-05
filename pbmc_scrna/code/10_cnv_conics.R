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

regions <- data.frame(
  Chrom = c(7,7,7),
  Start = c(0, 110000000,61528020),
  End = c(58169653, 159345973,159345973),
  Length = c(58169653, 499345973,97817953),
  row.names = c("7p","X7del", "7q")
)

output_file <- paste0("../output/Conics_PT3_PBMC")
mat <- Read10X(paste0("../data/pearson_mds/"))
gene_pos=getGenePositions(rownames(mat))

# Run CNVer for pearson
both_mat <- mat
bo_CPM <- t(t(both_mat) /colSums(both_mat)*1000000)
bo_go_expr <- log2(bo_CPM/10+1)
bo_go_expr=filterMatrix(bo_go_expr,gene_pos[,"hgnc_symbol"],minCells=5)
bo_normFactor=calcNormFactors(bo_go_expr)
bo_l=plotAll(bo_go_expr, bo_normFactor,regions,gene_pos, paste0(output_file, "_both_CNVplots"))
chr7_total_bo <- data.frame(round(bo_l,3), barcode = rownames(bo_l), ngenes = colSums(both_mat > 0))
write.table(chr7_total_bo, file = paste0(output_file, "_both_stats.tsv"),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

