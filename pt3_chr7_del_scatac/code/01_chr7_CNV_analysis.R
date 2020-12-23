library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(BuenColors)
library(magrittr)
library(Matrix)
library(CONICSmat)
library(data.table)

set.seed(1)

# Set up Signac/Seurat object
# Note: make sure at least Seurat v 3.2 and Signac v 1.0 is installed...
annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

getGenePositionsCL= function(gene_names,ensembl_version="dec2016.archive.ensembl.org",species="human"){
  if (species=="human"){
    ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", value = "ensembl_gene_id", host=ensembl_version)
    gene_positions <- biomaRt::getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'), values =gene_names, mart = ensembl)
  }

  gene_positions=gene_positions[!duplicated(gene_positions[,2]),]
  gene_positions[which(gene_positions[,3]=="X"),3]=23
  gene_positions[which(gene_positions[,3]=="Y"),3]=24
  gene_positions[which(gene_positions[,3]=="MT"),3]=0
  gene_positions[which(nchar(gene_positions[,3])>2),3]=0
  gene_positions=gene_positions[order(as.numeric(gene_positions[,3]),decreasing=F),]
  return(gene_positions)
}
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"
gp <- CONICSmat::getGenePositions(unique(annotations$gene_name))

process_chr7_del<- function(fragment_file,single_cell_file,output_file){
  # Slower but convenient for the row parsing
  metadata <- read.csv(
    file = single_cell_file,
    header = TRUE,
    row.names = 1
  ) %>% dplyr::filter(cell_id != "None")
  

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
  
  cells_in <- CreateSeuratObject(
    counts = CA,
    assay = "peaks",
    meta.data = metadata
  )
  
  gene.activities <- GeneActivity(cells_in)
  gene.activities <- gene.activities[,(colSums(gene.activities) > 500)]
  
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
    dplyr::mutate(del_region = V2 > 110000000) %>%
    dplyr::group_by(V4, del_region) %>% dplyr::summarize(count = n()) %>%
    dplyr::group_by(V4) %>% mutate(hits = n()) %>% dplyr::filter(hits == 2) -> count_df
  
  chr7_total_df <- count_df %>% dplyr::group_by(V4) %>% dplyr::summarize(all = sum(count))
  chr7_total_del_df <- count_df  %>% dplyr::filter(del_region) 
  chr7_total_df$in_del <- chr7_total_del_df$count
  chr7_total_df$pct_in_del <- chr7_total_df$in_del / chr7_total_df$all*100
  
  chr7_total_df2 <- data.frame(chr7_total_df, l[chr7_total_df$V4,])
  write.table(chr7_total_df2, file = output_file,
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}
#process_chr7_del("Pearson-CD34-PT3_v12-mtMask_fragments.tsv.gz", "Pearson-CD34-PT3_v12-mtMask_singlecell.csv.gz", "../output/Pearson-CD34-PT3.chr7DelQC.tsv")
#process_chr7_del("erythroid/Pearson_Healthy_D12_1_v12-mtMask.fragments.tsv.gz", "erythroid/Pearson_Healthy_D12_1_singlecell.csv.gz", "../output/Pearson_Healthy_D12_1.chr7DelQC.tsv")
#process_chr7_del("erythroid/Pearson_Healthy_D12_2_v12-mtMask.fragments.tsv.gz", "erythroid/Pearson_Healthy_D12_2_singlecell.csv.gz", "../output/Pearson_Healthy_D12_2.chr7DelQC.tsv")
#process_chr7_del("erythroid/Pearson_Healthy_D6_1_v12-mtMask.fragments.tsv.gz", "erythroid/Pearson_Healthy_D6_1_singlecell.csv.gz", "../output/Pearson_Healthy_D6_1.chr7DelQC.tsv")
#process_chr7_del("erythroid/Pearson_Healthy_D6_2_v12-mtMask.fragments.tsv.gz", "erythroid/Pearson_Healthy_D6_2_singlecell.csv.gz", "../output/Pearson_Healthy_D6_2.chr7DelQC.tsv")
#process_chr7_del("pbmcs/Pearson-BMMNC-PT3_v12_fragments.tsv.gz", "pbmcs/Pearson-BMMNC-PT3_v12_singlecell.csv", "../output/Pearson-BMMNC.chr7DelQC.tsv")
#process_chr7_del("../pearson_asap/pearson_asap_fragments.tsv.gz", "../pearson_asap/singlecell.csv", "../output/Pearson-ASAP.chr7DelQC.tsv")
#process_chr7_del("data_in/pbmc_3donor_aggr_fragments.tsv.gz", "data_in/pbmc_3donor_aggr_singlecell.csv", "../output/Pearson-PBMC.chr7DelQC.tsv")
process_chr7_del("../../../pearson_mtscatac_large_data_files/input/pbmcs_scatac/mtscatac_paper/PBMC_H9_v12-mtMask_fragments.tsv.gz",
                 "../../pbmc_scatac/data/control_singlecell/PBMC_H9_v12-mtMask_singlecell.csv", "../output/Healthy-PBMC_H9.chr7DelQC.tsv")
process_chr7_del("../../../pearson_mtscatac_large_data_files/input/pbmcs_scatac/mtscatac_paper/PBMC_H10_v12-mtMask_fragments.tsv.gz",
                 "../../pbmc_scatac/data/control_singlecell/PBMC_H10_v12-mtMask_singlecell.csv", "../output/Healthy-PBMC_H10.chr7DelQC.tsv")


