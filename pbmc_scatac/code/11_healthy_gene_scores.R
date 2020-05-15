library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(data.table)

dir <- "../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/"

gene.coords <- genes(EnsDb.Hsapiens.v75, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

process_healthy_signac_gs <- function(exp){
  cells <- colnames(paste0("../../../pearson_mtscatac_large_data_files/output/", exp, "-RIP.rds"))
  
  # create a gene by cell matrix
  gene.activities <- FeatureMatrix(
    fragments = paste0(dir, exp, "_v12-mtMask_fragments.tsv.gz"),
    features = genebodyandpromoter.coords,
    cells = cells,
    chunk = 20
  )
  
  gene.key <- genebodyandpromoter.coords$gene_name
  names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
  rownames(gene.activities) <- gene.key[rownames(gene.activities)]
  gene.activities
}

h9 <- process_healthy_signac_gs("PBMC_H9")
h10 <- process_healthy_signac_gs("PBMC_H10")

colnames(h9) <- gsub("-1", "-4", colnames(h9))
colnames(h10) <- gsub("-1", "-5", colnames(h10))

cg <- intersect(rownames(h9), rownames(h10))
gene.activities_combined <- cbind(h9[cg,], h10[cg,])

saveRDS(gene.activities_combined, file = "../../../pearson_mtscatac_large_data_files/output/Healthy_26April2020_geneActivityScores.rds")
