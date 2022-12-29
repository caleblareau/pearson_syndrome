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
l1a <- import_kite_counts("PT1_rep1_mtDOGMA_TSA_preAMP_S1", "1")
l1b <- import_kite_counts("PT1_rep1_mtDOGMA_TSA_sup_S3", "1")
common_l1s <- intersect(colnames(l1a), colnames(l1b))

l2a <- import_kite_counts("PT1_rep2_mtDOGMA_TSA_preAMP_S2", "2")
l2b <- import_kite_counts("PT1_rep2_mtDOGMA_TSA_sup_S4", "2")
common_l2s <- intersect(colnames(l2a), colnames(l2b))

new_set <- c(common_l1s,common_l2s)
full_protein_mat <- cbind(l1a[,common_l1s] + l1b[,common_l1s], l2a[,common_l2s] + l2b[,common_l2s])[,new_set]

# Import other attributes
# Import RNAseq
# Import scRNA-seq
process_scRNA  <- function(raw){
  
  # Remove crazy high and low expressors
  n_feature_rna <- colSums(raw > 0)
  n_total_rna <- colSums(raw)
  pct_mito <- colSums(raw[grepl("^MT", rownames(raw)), ])/n_total_rna * 100
  qc_cells <- colnames(raw)[pct_mito < 30 & n_total_rna > 1000 & n_feature_rna > 500]
  
  # Filter for singlet cells and non-mito genes
  raw <- raw[!grepl("^MT", rownames(raw)), qc_cells]
  raw
}

# Import the RNA
rep2_rna <- Read10X_h5("../data/PT1_rep2_mtDOGMA_gex.h5")
colnames(rep2_rna) <- gsub("-1", "-2", colnames(rep2_rna))
rep2_rna <- process_scRNA(rep2_rna)
rep1_rna <- process_scRNA(Read10X_h5("../data/PT1_rep1_mtDOGMA_gex.h5"))

# Import the meta data and create seurat object
rep1_qc <- fread("../data/qc/DOGMA_PT1_Rep1_per_barcode_metrics.csv.gz") %>% dplyr::filter(is_cell == 1) %>% data.frame()
rep2_qc <- fread("../data/qc/DOGMA_PT1_Rep2_per_barcode_metrics.csv.gz")%>% dplyr::filter(is_cell == 1) %>% data.frame()
rep2_qc$barcode <- gsub("-1", "-2", rep2_qc$barcode)
rownames(rep1_qc) <-rep1_qc$barcode; rownames(rep2_qc) <- rep2_qc$barcode

cellstate_meta <- merge(rbind(rep1_qc, rep2_qc),
                   data.frame(barcode = colnames(full_protein_mat), 
                              totalADT = colSums(full_protein_mat), 
                              totalCTRLadt = colSums(full_protein_mat[grepl("control", rownames(full_protein_mat)), ])), by = "barcode")

# Now append mtDNA meta
mito_meta <- fread("../output/Pearson_DOGMA_full_meta_data.tsv")[,c("cb", "orig.ident", "clip_het", "cov_het","coverage", "predicted.celltype.l2")]
full_meta <- merge(cellstate_meta, mito_meta,by.x = "barcode", by.y = "cb")
rownames(full_meta) <- full_meta$barcode
dim(full_meta)

# Do some visualization to determine thresholds
ggplot(full_meta, aes(x = totalCTRLadt, y = totalADT)) +
  geom_point() + scale_y_log10() + scale_x_log10() 
sum(full_meta$totalADT > 1e4 | full_meta$totalCTRLadt > 9 | full_meta$totalADT < 100)

# Set up Seurat object 
pbmc_pearson_dogma <- CreateSeuratObject(counts = cbind(rep1_rna,rep2_rna), 
                               meta.data = full_meta)
pbmc_pearson_dogma$pct_in_peaks <- pbmc_pearson_dogma$atac_peak_region_fragments / pbmc_pearson_dogma$atac_fragments *100
qplot(pbmc_pearson_dogma$pct_in_peaks)

# Prospectively subset on attributes
pbmc_pearson_dogma <- subset(pbmc_pearson_dogma,
                             subset = pct_in_peaks > 50 & totalCTRLadt < 10 & totalADT > 50 & nCount_RNA < (10^4.5))
dim(pbmc_pearson_dogma)

# Import ATAC-seq
peaks <- Read10X_h5("../../../pearson_large_data_files/input/pbmcs_dogma/dogma_pearson_filtered_feature_bc_matrix_aggr.h5")[["Peaks"]]

chrom_assay <- CreateChromatinAssay(
  counts = peaks[,colnames(pbmc_pearson_dogma)],
  sep = c(":", "-"),
  genome = NULL,
  fragments = '../../../pearson_large_data_files/input/pbmcs_dogma/atac_fragments.tsv.gz',
  min.cells = 0,
  min.features = 1
)
pbmc_pearson_dogma[["ADT"]] <- CreateAssayObject(full_protein_mat[,colnames(pbmc_pearson_dogma)])
pbmc_pearson_dogma[["peaks"]] <- chrom_assay

# Now do dimension reduction for each individually

# RNA Seurat stuff
DefaultAssay(pbmc_pearson_dogma) <- "RNA"
pbmc_pearson_dogma <- pbmc_pearson_dogma  %>% 
  NormalizeData() %>% ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(verbose = FALSE, assay = "RNA", reduction.name = "pca") 

# LSI dim reduction
DefaultAssay(pbmc_pearson_dogma) <- "peaks"
pbmc_pearson_dogma <- RunTFIDF(pbmc_pearson_dogma) %>% 
  FindTopFeatures( min.cutoff = 'q25') %>%
  RunSVD() 

# Do it for ADT
DefaultAssay(pbmc_pearson_dogma) <- "ADT"
pbmc_pearson_dogma <- pbmc_pearson_dogma  %>% 
  NormalizeData(assay = "ADT", normalization.method = "CLR", margin = 2) %>%
  ScaleData(assay = "ADT", do.scale = FALSE) %>%
  FindVariableFeatures(assay = "ADT") %>% 
  RunPCA(verbose = FALSE, assay = "ADT", reduction.name = 'apca') 

# Now run multimodal neighbors and embedding
pbmc_pearson_dogma <- FindMultiModalNeighbors(object = pbmc_pearson_dogma,
                                    reduction.list = list("pca", "lsi", "apca"),
                                    dims.list = list(1:30, 2:30, 1:20))
pbmc_pearson_dogma <- RunUMAP(pbmc_pearson_dogma, nn.name = "weighted.nn", reduction.name = "wnn.3.umap", reduction.key = "Uw3_" )
pbmc_pearson_dogma <- FindClusters(pbmc_pearson_dogma, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)

FeaturePlot(pbmc_pearson_dogma,
            features = c("peaks.weight","RNA.weight", "ADT.weight", "cov_het", "CD8A", "CD4", "CD20", "CD14", "CD45RA", "TCRv72", "TCRvd2"),
            min.cutoff = "q10", max.cutoff = "q90",
            reduction =  'wnn.3.umap',  pt.size = 0.1,
            by.col = FALSE) &
  scale_color_viridis()

FeaturePlot(pbmc_pearson_dogma,
            features = c("NCAM1", "CD19", "CD56"),
            min.cutoff = "q10", max.cutoff = "q90",
            reduction =  'wnn.3.umap',  pt.size = 0.1,
            by.col = FALSE) &
  scale_color_viridis()



DimPlot(pbmc_pearson_dogma, label = TRUE)
pO <- DimPlot(pbmc_pearson_dogma, reduction = "wnn.3.umap", label = FALSE) +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(pO, file = "../plots/umap_clusters_dogmawnn.png", width = 4, height = 4, dpi = 500)

pF <- FeaturePlot(pbmc_pearson_dogma, features = c( "CD4"), 
                          reduction =  'wnn.3.umap',  min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.01) +
  theme_void() + scale_color_viridis() + ggtitle("")+
  theme(legend.position = "none")
cowplot::ggsave2(pF, file = "../plots/CD4_expression.png", 
                 width = 4, height = 4, dpi = 500)

pH <- FeaturePlot(pbmc_pearson_dogma, features = c( "cov_het"), 
                  reduction =  'wnn.3.umap',   pt.size = 0.01) +
  theme_void() + scale_color_viridis() + ggtitle("")+
  theme(legend.position = "none")

cowplot::ggsave2(pH, file = "../plots/coverage_heteroplasmy.png", 
                 width = 4, height = 4, dpi = 500)

                        