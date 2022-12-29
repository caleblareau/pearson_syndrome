library(data.table)
library(dplyr)
library(BuenColors)
library(SummarizedExperiment)
library(Matrix)

estimate_coverage_heteroplasmy <- function(coord1, coord2, mat){
  idx <- 1:dim(mat)[1]
  in_del_boo <- (coord1 <= idx & coord2 >= idx)
  in_del <- colSums(mat[in_del_boo,])/sum(in_del_boo)
  out_del <- colSums(mat[!in_del_boo,])/sum(!in_del_boo)
  x <- in_del/out_del
  ifelse(x >1, 0, 1-x ) *100
}
se1 <- readRDS("../../../pearson_large_data_files/input/pbmcs_dogma/DOGMA_PT1_Rep1_mask_mgatk.rds")
het1 <- estimate_coverage_heteroplasmy(6073,13095, assays(se1)[["coverage"]])

se2 <- readRDS("../../../pearson_large_data_files/input/pbmcs_dogma/DOGMA_PT1_Rep2_mask_mgatk.rds")
het2 <- estimate_coverage_heteroplasmy(6073,13095,assays(se2)[["coverage"]])

df1 <- data.frame(
  cell = names(het1),
  cov_het = unname(het1), 
  coverage = se1$depth
)

df2 <- data.frame(
  cell = gsub("-1", "-2", names(het2)),
  cov_het = unname(het2), 
  coverage = se2$depth
)
dfboth <- rbind(df1, df2)

refmap<- rbind(
  fread("../data/DOGMA_PT1_Rep1_rna_refmapped.csv.gz"),
  fread("../data/DOGMA_PT1_Rep2_rna_refmapped.csv.gz") %>% mutate(cb = gsub("-1", "-2", cb))
)
mdf_initial<- merge(refmap, dfboth, by.x = "cb", by.y = "cell")
clipped_reads_df <- rbind(
  fread("../data/DOGMA_PT1_Rep1_delquant.deletion_heteroplasmy.tsv"),
  fread("../data/DOGMA_PT1_Rep2_delquant.deletion_heteroplasmy.tsv") %>% mutate(cell_id = gsub("-1", "-2", cell_id))) %>%
  dplyr::filter(version == "naive") 
colnames(clipped_reads_df) <- c("cell_id", "clip_het", colnames(clipped_reads_df)[3:8])
merge_df <- merge(mdf_initial, clipped_reads_df, by.x = "cb", by.y = "cell_id")
merge_df %>%
  group_by(predicted.celltype.l2) %>%
  dplyr::filter(reads_all > 5) %>%
  summarize(sum(clip_het < 0.01)/n()) %>% data.frame()

library(ggbeeswarm)
ggplot(merge_df %>% dplyr::filter(coverage > 10), aes(x = predicted.celltype.l2, y = cov_het)) +
  geom_quasirandom() + ggtitle("Clipped- read based heteroplasmy; min 10x reads spanning/supporting deletion (n = ~1000)")

ggplot(merge_df %>% dplyr::filter(coverage > 10), aes(x = cov_het, y = clip_het )) + geom_point() +
  labs(x = "Coverage heteroplasmy", y = "Clip heteroplasmy")

write.table(merge_df, file = "Pearson_DOGMA_full_meta_data.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
