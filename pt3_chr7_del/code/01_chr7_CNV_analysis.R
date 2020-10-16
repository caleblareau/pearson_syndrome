library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(BuenColors)
library(dplyr)
library(Matrix)
library(CONICSmat)
library(data.table)

set.seed(1)

# Set up Signac/Seurat object
# Note: make sure at least Seurat v 3.2 and Signac v 1.0 is installed...
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

process_chr7_del<- function(fragment_file,single_cell_file,output_file){
  # Slower but convenient for the row parsing
  metadata <- read.csv(
    file = single_cell_file,
    header = TRUE,
    row.names = 1
  ) %>% dplyr::filter(cell_id != "None")
  
  # Filter for Channel 3
  metadata <- metadata[substr(rownames(metadata), 18,18)=="3",]
  
  # change to UCSC style since the data was mapped to hg19
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg19"
  mat <- t(data.matrix(metadata[,c(1,2,3)]))
  rownames(mat) <- c("chr1-1-2", "chr2-3-4", "chr3-1-3")
  
  mat <- mat[,colSums(mat) > 500]
  CA <- CreateChromatinAssay(
    counts = (mat),
    sep = c("-", "-"),
    genome = 'hg19',
    fragments = fragment_file,
    min.cells = 1,
    min.features = 1
  )
  
  # add the gene information to the object
  Annotation(CA) <- annotations
  
  pearson_CD34 <- CreateSeuratObject(
    counts = CA,
    assay = "peaks",
    meta.data = metadata
  )
  
  gene.activities <- GeneActivity(pearson_CD34)
  gene.activities <- gene.activities[,(colSums(gene.activities) > 500) & (substr(colnames(gene.activities), 18, 18) == "3")]
  
  CPM <- t(t(gene.activities) /colSums(gene.activities)*1000000)
  go_expr <- log2(CPM/10+1)
  gene_pos=getGenePositions(rownames(go_expr))
  
  regions <- data.frame(
    Chrom = c(7,7,7),
    Start = c(0, 110000000,61528020),
    End = c(58169653, 159345973,159345973),
    Length = c(58169653, 499345973,97817953),
    row.names = c("7p","X7del", "7q")
  )
  head(regions,n=5)
  go_expr=filterMatrix(go_expr,gene_pos[,"hgnc_symbol"],minCells=5)
  normFactor=calcNormFactors(go_expr)
  l=plotAll(go_expr,normFactor,regions,gene_pos, paste0(output_file, "_CNVplots"))
  
  consensus_barcodes <-rownames(l)
  pf <- data.table::fread(fragment_file)
  
  pf %>%  dplyr::filter(V4 %in% consensus_barcodes & V1 == "chr7") %>%
    mutate(del_region = V2 > 110000000) %>%
    group_by(V4, del_region) %>% summarize(count = n()) %>%
    group_by(V4) %>% mutate(hits = n()) %>% dplyr::filter(hits == 2) -> count_df
  
  chr7_total_df <- count_df %>% group_by(V4) %>% summarize(all = sum(count))
  chr7_total_del_df <- count_df  %>% dplyr::filter(del_region) 
  chr7_total_df$in_del <- chr7_total_del_df$count
  chr7_total_df$pct_in_del <- chr7_total_df$in_del / chr7_total_df$all*100
  
  chr7_total_df2 <- data.frame(chr7_total_df, l[chr7_total_df$V4,])
  write.table(chr7_total_df2, file = output_file,
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}
#process_chr7_del("Pearson-CD34-PT3_v12-mtMask_fragments.tsv.gz", "Pearson-CD34-PT3_v12-mtMask_singlecell.csv.gz", "Pearson-CD34-PT3.chr7DelQC.tsv")
#process_chr7_del("erythroid/Pearson_Healthy_D12_1_v12-mtMask.fragments.tsv.gz", "erythroid/Pearson_Healthy_D12_1_singlecell.csv.gz", "Pearson_Healthy_D12_1.chr7DelQC.tsv")
#process_chr7_del("erythroid/Pearson_Healthy_D12_2_v12-mtMask.fragments.tsv.gz", "erythroid/Pearson_Healthy_D12_2_singlecell.csv.gz", "Pearson_Healthy_D12_2.chr7DelQC.tsv")
#process_chr7_del("erythroid/Pearson_Healthy_D6_1_v12-mtMask.fragments.tsv.gz", "erythroid/Pearson_Healthy_D6_1_singlecell.csv.gz", "Pearson_Healthy_D6_1.chr7DelQC.tsv")
#process_chr7_del("erythroid/Pearson_Healthy_D6_2_v12-mtMask.fragments.tsv.gz", "erythroid/Pearson_Healthy_D6_2_singlecell.csv.gz", "Pearson_Healthy_D6_2.chr7DelQC.tsv")
#process_chr7_del("pbmcs/Pearson-BMMNC-PT3_v12_fragments.tsv.gz", "pbmcs/Pearson-BMMNC-PT3_v12_singlecell.csv", "Pearson-BMMNC.chr7DelQC.tsv")
#process_chr7_del("../pearson_asap/pearson_asap_fragments.tsv.gz", "../pearson_asap/singlecell.csv", "Pearson-ASAP.chr7DelQC.tsv")
process_chr7_del("data_in/pbmc_3donor_aggr_fragments.tsv.gz", "data_in/pbmc_3donor_aggr_singlecell.csv", "Pearson-PBMC.chr7DelQC.tsv")

d <- fread("Pearson-PBMC.chr7DelQC.tsv")
ggplot(d, aes(x = pct_in_del, y = X7del)) +
  geom_point(size = 0.2)
