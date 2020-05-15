library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(BuenColors)

peaks_pbmc_aggr <-diffloop::bedToGRanges("../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_filtered_peak_bc_matrix/peaks.bed.gz")

# function to get counts
getCountsFromFrags <- function(frag_gz_file,
                               peaks_gr,
                               barcodes){
  
  # Make GRanges of fragments that are solid for the cells that we care about
  frags_in <- data.table::fread(paste0(frag_gz_file)) %>% 
    data.frame() %>% filter(V4 %in% barcodes) 
  frags_valid <- c(
    GenomicRanges::makeGRangesFromDataFrame(frags_in, seqnames.field = "V1", start.field = "V2", end.field = "V2", keep.extra.columns = TRUE),
    GenomicRanges::makeGRangesFromDataFrame(frags_in, seqnames.field = "V1", start.field = "V3", end.field = "V3", keep.extra.columns = TRUE)
  )
  
  # Get a denominator, per cell
  denom <- table(GenomicRanges::mcols(frags_valid)$V4)
  barcodes_found <- names(denom)
  
  # Get the overlaps with peaks
  ovPEAK <- GenomicRanges::findOverlaps(peaks_gr, frags_valid)
  
  # Establish a numeric index for the barcodes for sparse matrix purposes
  id <- factor(as.character(GenomicRanges::mcols(frags_valid)$V4), levels = barcodes_found)
  
  # Make sparse matrix with counts with peaks by  unique barcode
  countdf <- data.frame(peaks = S4Vectors::queryHits(ovPEAK),
                        sample = as.numeric(id)[S4Vectors::subjectHits(ovPEAK)]) %>%
    dplyr::group_by(peaks,sample) %>% dplyr::summarise(count = n()) %>% data.matrix()
  
  m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks_gr)),
                            j = c(countdf[,2], length(barcodes_found)),
                            x = c(countdf[,3],0))
  colnames(m) <- barcodes_found
  
  # Make a polished colData
  colData <- data.frame(
    sample = barcodes_found,
    depth = as.numeric(denom),
    FRIP = Matrix::colSums(m)/as.numeric(denom)
  )
  # Make sure that the SE can be correctly constructed
  stopifnot(all(colData$sample == colnames(m)))
  
  # Make summarized Experiment
  SE <- SummarizedExperiment::SummarizedExperiment(
    rowRanges = peaks_gr,
    assays = list(counts = m),
    colData = colData
  )
  return(SE)
}

# Import experiment to make an RDS of the SE and a basic qc plot
importExperiment <- function(exp, peaks, frip_threshold){
  
  qcdf <- fread(paste0("../data/control_singlecell/", exp,"_v12-mtMask_singlecell.csv.gz"), header = TRUE, sep = ",") %>% 
    data.frame() %>% filter(cell_id != "None")
  bc <- as.character(qcdf$barcode)
  SE <- getCountsFromFrags(paste0("../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/", exp, "_v12-mtMask_fragments.tsv.gz"), peaks, bc)
  
  df <- data.frame(colData(SE))

  # Import mtDNA DF
  cov_mtDNA <- fread(paste0("../data/depth/",exp,"_v12-mtMask_mgatk.depthTable.txt"), header = FALSE, col.names = c("sample", "depth"))
  vec <- as.numeric(cov_mtDNA$depth); names(vec) <- as.character(cov_mtDNA$sample)
  colData(SE)$mtDNAcoverage <- vec[as.character(df$sample)] %>% unname()
  df$mtDNAcoverage <- vec[as.character(df$sample)] %>% unname()
  df$keep <- log10(df$depth) >= 3 & df$FRIP >= frip_threshold & df$mtDNAcoverage >= 20
  
  SE2 <- SE[, df$keep]
  saveRDS(SE2, file = paste0("../../../pearson_mtscatac_large_data_files/output/", exp, "-RIP.rds"))
  print(dim(SE2))
}

# Process healthy cells using existing peak set
importExperiment("PBMC_H10", peaks_pbmc_aggr, 0.25)
importExperiment("PBMC_H9", peaks_pbmc_aggr, 0.25)

