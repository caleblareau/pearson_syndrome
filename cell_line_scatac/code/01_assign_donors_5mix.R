library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
"%ni%" <- Negate("%in%")

computeAFMutMatrix <- function(SE_go){
  cov <-  (assays(SE_go)[[paste0("A", "_counts_fw")]] + assays(SE_go)[[paste0("A", "_counts_rev")]]) +
    (assays(SE_go)[[paste0("C", "_counts_fw")]] + assays(SE_go)[[paste0("C", "_counts_rev")]]) +
    (assays(SE_go)[[paste0("G", "_counts_fw")]] + assays(SE_go)[[paste0("G", "_counts_rev")]]) +
    (assays(SE_go)[[paste0("T", "_counts_fw")]] + assays(SE_go)[[paste0("T", "_counts_rev")]]) + 0.0001
  ref_allele <- as.character(rowRanges(SE_go)$refAllele)
  
  getMutMatrix <- function(letter){
    mat <- (assays(SE_go)[[paste0(letter, "_counts_fw")]] + assays(SE_go)[[paste0(letter, "_counts_rev")]]) / cov
    rownames(mat) <- paste0(as.character(1:dim(mat)[1]), toupper(ref_allele), ">", letter)
    return(mat[toupper(ref_allele) != letter,])
  }
  
  rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
  
}

# Import AFs
SE <- readRDS("Pearson-Mix-5_v12-mtMask_mgatk/final/Pearson-Mix-5_v12-mtMask_mgatk.rds")
filt <- colData(SE)$depth >= 15 & colData(SE)$depth < 100
SE2 <- SE[,filt]
af <- computeAFMutMatrix( SE2 )
sum(af > 1)

x <- (Matrix::rowSums(af > 0.99) > 100) & (Matrix::rowSums(af < 0.005) > 100)
sum(x)
xx <- x & rownames(af) %ni% c("3107N>C", "3107N>A", "567A>C", "3107N>T")
afp <- af[xx,]

# Try something high dimensional 



dm <- t(data.matrix(afp))
classify <- case_when(
  dm[,"709G>A"] > 0.9 & dm[,"11251A>G"] > 0.9 ~ "aPD1",
  dm[,"3221A>G"] > 0.9 & dm[,"94G>A"] > 0.9 ~ "aPD2",
  dm[,"152T>C"] > 0.9 & dm[,"4646T>C"] > 0.9 & dm[,"3221A>G"] < 0.2 & dm[,"73A>G"] < 0.2 ~ "aPD3",
  
  dm[,"73A>G"] > 0.9 & dm[,"8411A>C"] > 0.9 & dm[,"3221A>G"] < 0.2 & dm[,"11251A>G"] < 0.2~ "C1",
  dm[,"8602T>C"] > 0.9 & dm[,"3010G>A"] > 0.9 & dm[,"73A>G"] < 0.2 ~ "C2",
  TRUE ~ "unassigned"
)
vv <- c("709G>A", "11251A>G", "3221A>G", "94G>A", "152T>C", "4646T>C", "73A>G",  "8411A>C", "8602T>C", "3010G>A")

assign_df <- data.frame(
  barcode = colnames(afp), fp_classify = classify,
  mean_cov = colData(SE2)$depth
) 

dm <- fread("../demuxlet/mix5/scSplit_result.csv")
merged_df <- merge(assign_df, dm, by.x = "barcode", by.y = "Barcode")
table(merged_df$fp_classify, merged_df$Cluster)

dm_vec <- c("aPD1", "aPD2", "aPD3", "C1", "C2", "unassigned")
names(dm_vec) <- c("SNG-0", "SNG-2", "SNG-1", "SNG-3", "SNG-5", "DBL-4")
merged_df$classify <- dm_vec[as.character(merged_df$Cluster)]
assign_df <- merged_df%>% arrange(fp_classify)

assign_vec <- c("aPD1" = "#2780FF", "aPD2" = "#960096", "aPD3" = "#000090",
                "C1" = "#333333", "C2" = "#999999", "unassigned" = "firebrick")

ha_col2 <- HeatmapAnnotation(df = data.frame(Assignment = assign_df$fp_classify),
                             col = list(Assignment = assign_vec))
ss <- c("4216T>C","16270C>T", "4646T>C",  "8411A>C","3010G>A")
pdf(paste0("output/Pearson_Mix_Homoplasmic5-all.pdf"), width=5, height=1.2)
hm <- Heatmap(data.matrix(afp[vv,as.character(assign_df$barcode)]), 
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 6),
              top_annotation=ha_col2,
              cluster_rows = FALSE, name = "AF",
              show_column_names = FALSE)
hm
dev.off()

png(paste0("output/Pearson_Mix_Homoplasmic5-all.png"), width=20, height=4.8, res = 1000, units =  "in")
hm <- Heatmap(data.matrix(afp[vv,as.character(assign_df$barcode)]), 
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 6),
              top_annotation=ha_col2,
              cluster_rows = FALSE, name = "AF",
              show_column_names = FALSE)
hm
dev.off()


