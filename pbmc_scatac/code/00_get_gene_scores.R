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

cells <- fread(paste0(dir, "pbmc_3donor_aggr_filtered_peak_bc_matrix/barcodes.tsv.gz"), header = FALSE)[[1]]

# create a gene by cell matrix
gene.activities <- FeatureMatrix(
  fragments = paste0(dir, "pbmc_3donor_aggr_fragments.tsv.gz"),
  features = genebodyandpromoter.coords,
  cells = cells,
  chunk = 20
)

gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

saveRDS(gene.activities, file = "../../../pearson_mtscatac_large_data_files/output/Perason_aggr_26April2020_geneActivityScores.rds")
