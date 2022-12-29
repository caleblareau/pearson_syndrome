library(Seurat)
library(data.table)
library(Signac)
library(dplyr)
library(Matrix)
library(harmony)
library(viridis)
library(mclust)
library(BuenColors)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)

# Import protein counts
import_kite_counts <- function(lib, nl){
  mtx <- fread(paste0("../data/protein_counts/",lib,"_counts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0("../data/protein_counts/",lib,"_counts.barcodes.txt.gz"), header = FALSE)[[1]], "-", nl)
  colnames(matx) <- (fread(paste0("../data/protein_counts/",lib,"_counts.genes.txt.gz"), header = FALSE)[[1]])
  return(t(matx))
}
library(dplyr)
bcs <- data.frame(data.table::fread('../output/Pearson_DOGMA_full_meta_data.tsv')) %>% dplyr::pull(cb)
l1a <- import_kite_counts("PT1_rep1_mtDOGMA_TSA_preAMP_S1", "1")
l1b <- import_kite_counts("PT1_rep1_mtDOGMA_TSA_sup_S3", "1")

l2a <- import_kite_counts("PT1_rep2_mtDOGMA_TSA_preAMP_S2", "2")
l2b <- import_kite_counts("PT1_rep2_mtDOGMA_TSA_sup_S4", "2")

write.table(data.frame(data.matrix(l1a[,colnames(l1a) %in% bcs])),
            file = "../output/PT1_rep1_mtDOGMA_TSA_preAMP.csv", row.names = TRUE, col.names = TRUE, sep = ",", quote = FALSE)
write.table(data.frame(data.matrix(l1b[,colnames(l1b) %in% bcs])),
            file = "../output/PT1_rep1_mtDOGMA_TSA_sup.csv", row.names = TRUE, col.names = TRUE, sep = ",", quote = FALSE)
write.table(data.frame(data.matrix(l2a[,colnames(l2a) %in% bcs])),
            file = "../output/PT1_rep2_mtDOGMA_TSA_preAMP.csv", row.names = TRUE, col.names = TRUE, sep = ",", quote = FALSE)
write.table(data.frame(data.matrix(l2b[,colnames(l2b) %in% bcs])),
            file = "../output/PT1_rep2_mtDOGMA_TSA_sup.csv", row.names = TRUE, col.names = TRUE, sep = ",", quote = FALSE)


