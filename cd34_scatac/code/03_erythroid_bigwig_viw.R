library(data.table)
library(dplyr)
library(matchingR)
library(GenomicRanges)
library(rtracklayer)

load("../output/CD34_umap_embedding_granja_proj3.rda")

pearson <- projection_df_pearson_pt %>% filter(celltype == "Pearson" & umap1 < -7.5)
healthy <- projection_df_control %>% filter(celltype == "Healthy" & umap1 < -7.5)

dist <- as.matrix(pdist::pdist(data.matrix(pearson[,c(2,3)]),
      data.matrix(healthy[,c(2,3)]))) * -1
colnames(dist) <- healthy$barcode
rownames(dist) <- pearson$barcode

pearson_matches <- galeShapley.collegeAdmissions(studentUtils = dist, collegeUtils = t(dist), slots = 1)

healthy$in_closest <- 1:dim(healthy)[1] %in% pearson_matches$matched.colleges[,1]
pearson$in_closest <- "Duh"
ggplot(rbind(healthy, pearson), aes(x = umap1, y = umap2, color = in_closest)) + 
  geom_point()

# Output matched healthy/pearson barcodes
healthy_bc <- healthy %>% filter(in_closest) %>% pull(barcode) %>% as.character()
pearson_bc <- pearson %>% pull(barcode) %>% as.character()

# Export healthy MEP bigwig
frags_healthy <- fread("../../../pearson_mtscatac_large_data_files/input/CD34/CD34_G10_v12-mtMask_fragments.tsv.gz")
healthy_gr <- frags_healthy %>% filter(V4 %in% healthy_bc) %>%
  setnames(c("chr", "start", "end", "V4", "PCRn")) %>%
  makeGRangesFromDataFrame()

healthy_coverage <- coverage(healthy_gr)/length(healthy_gr)*1000000
export.bw(healthy_coverage, con = paste0("../output/MEP_BW/Healthy_matched_MEP.bw"))

# Now export pearson
frags_pearson <- fread("../../../pearson_mtscatac_large_data_files/input/CD34/Pearson-CD34-PT3_v12-mtMask_fragments.tsv.gz")
pearson_gr <- frags_pearson %>% filter(V4 %in% pearson_bc) %>%
  setnames(c("chr", "start", "end", "V4", "PCRn")) %>%
  makeGRangesFromDataFrame()

pearson_coverage <- coverage(pearson_gr)/length(pearson_gr)*1000000
export.bw(pearson_coverage, con = paste0("../output/MEP_BW/Pearson_MEP.bw"))


healthy_counts <- assays(readRDS("../../../pearson_mtscatac_large_data_files/output/CD34_G10_SummarizedExperiment.rds"))[["counts"]][,healthy_bc]
pearson_counts <- assays(readRDS("../../../pearson_mtscatac_large_data_files/output/Pearson-CD34-PT3_SummarizedExperiment.rds"))[["counts"]][,pearson_bc]

d <- data.frame(rowRanges(readRDS("../../../pearson_mtscatac_large_data_files/output/CD34_G10_SummarizedExperiment.rds")),
           h = rowSums(healthy_counts), p = rowSums(pearson_counts))
d$h <- d$h / sum(d$h) * 1000000
d$p <- d$p / sum(d$p) * 1000000

d$diff <- log2((d$p + 0.1)/(d$h+0.1))
d$region <- paste0(d$seqnames,":", d$start, "-", d$end)
d  %>% filter(diff < - 6) %>% dplyr::select(c("seqnames", "start","end")) %>%
  write.table(row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
