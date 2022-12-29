library(dplyr)
library(data.table)
library(Matrix)

import_kite_counts <- function(library, id){
  mtx <- fread(paste0("../data/pstag",library,"_fc/featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0("../data/pstag",library,"_fc/featurecounts.barcodes.txt.gz"), header = FALSE)[[1]], "-", id)
  colnames(matx) <- paste0(fread(paste0("../data/pstag",library,"_fc/featurecounts.genes.txt.gz"), header = FALSE)[[1]])
  return(t(matx))
}
bcs <- fread("../data/tcell-aggr-singlecell.csv.gz") %>% filter(is__cell_barcode == 1) %>% pull(barcode)

m1 <- import_kite_counts("1", "2")
m2 <- import_kite_counts("2", "3")

write.table(data.frame(data.matrix(m1[,colnames(m1) %in% bcs])),
            file = "../output/tag-PT3-Day14v1.csv", row.names = TRUE, col.names = TRUE, sep = ",", quote = FALSE)
write.table(data.frame(data.matrix(m2[,colnames(m2) %in% bcs])),
            file = "../output/tag-PT3-Day14v2.csv", row.names = TRUE, col.names = TRUE, sep = ",", quote = FALSE)
