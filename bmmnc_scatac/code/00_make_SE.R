library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(ggrastr)
library(BuenColors)

peaks_bmmnc <- rowRanges(readRDS("../../../pearson_mtscatac_large_data_files/input/bmmnc/scATAC-Healthy-Hematopoiesis-191120.rds"))

# function to get counts
getCountsFromFrags <- function(frag_gz_file,
                               peaks_gr,
                               barcodes){
  
  # Make GRanges of fragments that are solid for the cells that we care about
  frags_valid <- data.table::fread(paste0(frag_gz_file)) %>% 
    data.frame() %>% filter(V4 %in% barcodes) %>%  # filter for barcodes in our search set
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
  
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
  
  qcdf <- fread(paste0("../data/", exp,"_v12_singlecell.csv"), header = TRUE, sep = ",") %>% 
    data.frame() %>% filter(cell_id != "None")
  bc <- as.character(qcdf$barcode)
  SE <- getCountsFromFrags(paste0("../../../pearson_mtscatac_large_data_files/input/bmmnc/", exp, "_v12_fragments.tsv.gz"), peaks, bc)
  
  df <- data.frame(colData(SE))
  
  # Import mtDNA DF
  cov_mtDNA <- data.frame(fread(paste0("../data/",exp,"_v12-mtMask_mgatk.depthTable.txt"),
                                header = FALSE, col.names = c("sample", "depth")))
  vec <- as.numeric(cov_mtDNA$depth); names(vec) <- as.character(cov_mtDNA$sample)
  colData(SE)$mtDNAcoverage <- vec[as.character(df$sample)] %>% unname()
  df$mtDNAcoverage <- vec[as.character(df$sample)] %>% unname()
  df$keep <- log10(df$depth) >= 3 & df$FRIP >= frip_threshold & df$mtDNAcoverage >= 20
  
  SE2 <- SE[, df$keep]
  saveRDS(SE2, file = paste0("../../../pearson_mtscatac_large_data_files/output/", exp, "_SummarizedExperiment.rds"))
  write.table(df, file = paste0("../output/barcode_qc/", exp, ".barcode_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  print(dim(SE2))
}

# from granja et al for projections
importExperiment("Pearson-BMMNC-PT3", peaks_bmmnc, 0.25)



