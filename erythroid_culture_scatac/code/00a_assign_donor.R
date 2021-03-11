library(data.table)
library(Matrix)
library(SummarizedExperiment)

df2 <- data.frame(
  pos = c(16311, 8448, 961, 8618, 921, 16124, 189),
  Control = c("T", "T", "T", "C", "C", "C", "G"),
  Pearson = c("C", "C", "G", "T", "T", "T", "A")
)

# Function to get the per barcode counts per donor by alt letter
extractme <- function(donor, letter, SE){
  idxx <- df2[which(df2[,donor] == letter), "pos"]
  Matrix::colSums(assays(SE)[[paste0(letter, "_counts_fw")]][idxx, ,drop = FALSE] + assays(SE)[[paste0(letter, "_counts_rev")]][idxx, ,drop = FALSE])
}

# Given a summarized experiment from mgatk, compute the essentials for ultimately determining contamination
process_SE_counts <- function(lib){
  SE <- readRDS(paste0("../../../pearson_large_data_files/input/invitro_ery_scatac//mgatk_files/Pearson_Healthy_",lib,"_v12-mtMask_mgatk.rds"))
  df3 <- data.frame(
    barcode = colnames(SE),
    control = extractme("Control", "A", SE) +  extractme("Control", "C", SE) +  extractme("Control", "G", SE) +  extractme("Control", "T", SE),
    pearson = extractme("Pearson", "A", SE) +  extractme("Pearson", "C", SE) +  extractme("Pearson", "G", SE) +  extractme("Pearson", "T", SE), 
    depth = colData(SE)$depth
  )
  df3 <-  df3 %>% mutate(control_perc = control/(pearson + control + 0.001)*100, pearson_perc = pearson/(control + pearson + 0.001)*100)  %>%
    mutate(minor_population = pmin(control_perc,pearson_perc))
  df3$assign <-  ifelse(df3$depth < 10, "Low_coverage",
           ifelse( df3$control_perc > 98, "Control", 
                   ifelse(df3$pearson_perc > 98, "Pearson", 
                          "Collision")))
  table(df3$assign)
  write.table(df3, file = paste0("../output/cell_assignments_per_channel/Pearson_", lib, "_assign.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  df3
}

df6_1 <- process_SE_counts("D6_1")
df6_2 <- process_SE_counts("D6_2")

df12_1 <- process_SE_counts("D12_1")
df12_2 <- process_SE_counts("D12_2")

table(df6_1$assign)
table(df6_2$assign)
table(df12_1$assign)
table(df12_2$assign)

