library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(data.table)

# Get coordinates
gene.coords <- genes(EnsDb.Hsapiens.v75, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

process_ga <- function(lib){
  
  # Get fragment files and cells
  frag_file <- paste0("../../../pearson_mtscatac_large_data_files/input/invitro_ery_diff/fragments/Pearson_Healthy_",lib,"_v12-mtMask.fragments.tsv.gz")
  scdf <- fread(paste0("../data/singlecell/Pearson_Healthy_",lib,"_singlecell.csv.gz"), header = TRUE) 
  cells <- scdf %>% dplyr::filter(cell_id != "None") %>% pull(barcode)
  
  # create a gene by cell matrix
  gene.activities <- FeatureMatrix(
    fragments = frag_file,
    features = genebodyandpromoter.coords,
    cells = cells,
    chunk = 20
  )
  
  gene.key <- genebodyandpromoter.coords$gene_name
  names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
  rownames(gene.activities) <- gene.key[rownames(gene.activities)]
  saveRDS(gene.activities, file = paste0("../../../pearson_mtscatac_large_data_files/output/invitro_erythroid/", lib, "_gene_activities.rds"))
}

process_ga("D6_1")
process_ga("D6_2")
process_ga("D12_1")
process_ga("D12_2")
