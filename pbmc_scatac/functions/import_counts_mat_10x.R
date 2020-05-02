library(data.table)
library(Matrix)

# Accessory function to import the peaks x cells count matrix
import_counts_mat_10x <- function(directory){
  peaks <- fread(paste0(directory, "/peaks.bed.gz"), header = FALSE, col.names = c("chr", "start", "end"))
  peak_character <- paste0(peaks$chr, ":", as.character(peaks$start), "-", as.character(peaks$end))
  barcodes <- fread(paste0(directory, "/barcodes.tsv.gz"), header = FALSE)[[1]]
  mtx <- fread(paste0(directory, "/matrix.mtx.gz")[[1]], skip = 3, header = FALSE)
  
  # Assemble binary matrix
  mat <- Matrix::sparseMatrix(i = c(mtx[[1]],length(peaks)), 
                              j = c(mtx[[2]],length(barcodes)),
                              x = c(as.numeric(mtx[[3]]) > 0, 0))
  colnames(mat) <- barcodes
  rownames(mat) <- peak_character
  mat
}