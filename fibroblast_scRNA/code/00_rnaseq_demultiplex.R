library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
"%ni%" <- Negate("%in%")

computeAFMutMatrix <- function(SE){
  cov <- assays(SE)[["coverage"]]+ 0.001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  
  getMutMatrix <- function(letter){
    mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
    rownames(mat) <- paste0(as.character(1:dim(mat)[1]), toupper(ref_allele), ">", letter)
    return(mat[toupper(ref_allele) != letter,])
  }
  
  rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
  
}

computeCovVec <- function(SE){
  cov <- assays(SE)[["coverage"]]+ 0.001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  
  getCovVec <- function(letter){
    vec <-rowMeans( cov)
    names(vec) <- paste0(as.character(1:length(vec)), toupper(ref_allele), ">", letter)
    return(vec[toupper(ref_allele) != letter])
  }
  
  c(getCovVec("A"), getCovVec("C"), getCovVec("G"), getCovVec("T"))
  
}


# Import AFs
SE <- readRDS("../data/Pearson-5p-mgatk.rds")
filt <- colData(SE)$depth > 3
SE2 <- SE[,filt]
af <- computeAFMutMatrix( SE2 )
cov_vec <- computeCovVec(SE2)

# Import mutations
muts <- read.table("../output/useful_SNPs.tsv")[,1] %>% as.character()
xx <- data.frame(cov = round(cov_vec[muts], 1), muts = muts) %>%
  filter(cov > 20) %>% pull(muts) %>% as.character()

afp <- af[c(xx),]

perc.rank <- function(x) trunc(rank(x))/length(x)
coverage_quantile <- perc.rank(colData(SE2)$depth)

ha_col <- HeatmapAnnotation(df = data.frame(CovRank = coverage_quantile),
                            col = list(CovRank = colorRamp2(c(0,1), c("white", "dodgerblue3"))))


dm <- t(data.matrix(afp))
classify <- case_when(
  dm[,"1888G>A"] > 0.95 & dm[,"8496T>C"] > 0.95 ~ "P1",
  dm[,"7768A>G"] > 0.95 ~ "P2",
  dm[,"7768A>G"] < 0.02 & dm[,"1888G>A"] < 0.02 & dm[,"8496T>C"] < 0.02 ~ "C1",
  TRUE ~ "unassigned"
)

assign_df <- data.frame(
  barcode = colnames(afp), classify,
  coverage_quantile ,
  mean_cov = colData(SE2)$depth
) %>% arrange(classify)
table(assign_df$classify)
assign_vec <- c("P1" = "firebrick", "P2" = "orange3", "C1" = "dodgerblue3", "unassigned" = "grey")

ha_col2 <- HeatmapAnnotation(df = data.frame(CovRank = assign_df$coverage_quantile, Assignment = assign_df$classify),
                             col = list(CovRank = colorRamp2(c(0,1), c("white", "dodgerblue3")),
                                        Assignment = assign_vec))
afp[afp >1] <- 1
pdf(paste0("../output/12July_5pRNA_assigned.pdf"), width=10, height=10)
hm <- Heatmap(data.matrix(afp[,as.character(assign_df$barcode)]), 
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              cluster_columns = FALSE,
              top_annotation=ha_col2,
              cluster_rows = TRUE, name = "AF",
              show_column_names = FALSE)
hm
dev.off()

write.table(assign_df, file = "../output/12July_5pRNA_assignments.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

