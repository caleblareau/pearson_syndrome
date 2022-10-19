library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
"%ni%" <- Negate("%in%")


df2 <- data.frame(
  pos = c(16292, 15452, 4216, 4216, 11812, 1811, 14167),
  Pearson = c("C", "C", "T", "A", "A", "G", "T"),
  Control = c("T", "A", "C", "G", "G", "A", "C")
)

# Function to get the per barcode counts per donor by alt letter
extractme <- function(donor, letter, SE){
  idxx <- df2[which(df2[,donor] == letter), "pos"]
  Matrix::colSums(assays(SE)[[paste0(letter, "_counts_fw")]][idxx, ,drop = FALSE] + assays(SE)[[paste0(letter, "_counts_rev")]][idxx, ,drop = FALSE])
}

# Given a summarized experiment from mgatk, compute the essentials for ultimately determining contamination
process_SE_counts <- function(SE, lib){
  df3 <- data.frame(
    barcode = colnames(SE),
    control = extractme("Control", "A", SE) +  extractme("Control", "C", SE) +  extractme("Control", "G", SE) +  extractme("Control", "T", SE),
    pearson = extractme("Pearson", "A", SE) +  extractme("Pearson", "C", SE) +  extractme("Pearson", "G", SE) +  extractme("Pearson", "T", SE), 
    depth = colData(SE)$depth
  )
  df3 <-  df3 %>% mutate(control_perc = control/(pearson + control + 0.001)*100, pearson_perc = pearson/(control + pearson + 0.001)*100)  %>%
    mutate(minor_population = pmin(control_perc,pearson_perc))
  df3$assign <-  ifelse(df3$depth < 10, "Low_coverage",
                        ifelse( df3$control_perc > 95, "Control", 
                                ifelse(df3$pearson_perc > 95, "Pearson", 
                                       "Collision")))
  table(df3$assign)
  
  df3
}

estimate_coverage_heteroplasmy <- function(coord1, coord2, mat){
  idx <- 1:dim(mat)[1]
  in_del_boo <- (coord1 <= idx & coord2 >= idx)
  in_del <- colSums(mat[in_del_boo,])/sum(in_del_boo)
  out_del <- colSums(mat[!in_del_boo,])/sum(!in_del_boo)
  x <- in_del/out_del
  ifelse(x >1, 0, 1-x ) *100
}

se1 <- readRDS("HC_PT1_mix_CD3_CD28_Day14_rep1_ATAC_tenx.rds")
se2 <- readRDS("HC_PT1_mix_CD3_CD28_Day14_rep2_ATAC_tenx.rds")
df1 <- process_SE_counts(se1, "Rep1")
df2 <- process_SE_counts(se2, "Rep2")

df1$heteroplasmy <- estimate_coverage_heteroplasmy(6073,13095, assays(se1)[["coverage"]])
df2$heteroplasmy <- estimate_coverage_heteroplasmy(6073,13095, assays(se2)[["coverage"]])

rbind(df1,df2) %>%
  filter(assign == "Pearson") %>%
  ggplot(aes(x = heteroplasmy)) + stat_ecdf()

write.table(df1, file = paste0("Pearson_HC_PT1_mix_CD3_CD28_Day14_rep1_assign.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(df2, file = paste0("Pearson_HC_PT1_mix_CD3_CD28_Day14_rep2_assign.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


