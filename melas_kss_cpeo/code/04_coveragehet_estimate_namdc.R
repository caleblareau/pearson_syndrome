library(data.table)
library(dplyr)
library(viridis)
library(BuenColors)
library(SummarizedExperiment)
library(Matrix)

# Simple function to estimate coverage-based heteroplasmy
estimate_coverage_heteroplasmy_mat<- function(coord1, coord2, mat){
  idx <- 1:dim(mat)[1]
  in_del_boo <- (coord1 <= idx & coord2 >= idx)
  in_del <- colSums(mat[in_del_boo,])/sum(in_del_boo)
  out_del <- colSums(mat[!in_del_boo,])/sum(!in_del_boo)
  x <- in_del/out_del
  ifelse(x >1, 0, 1-x ) *100
}


make_combined_df <- function(x){
  hetdf <- fread(paste0("../data/namdc-metadata/",x,"_del.new.tsv")) %>%
    filter(version == "improved")
  
  merge(
    hetdf, 
    fread(paste0("../data/reference_projections/namdc/",x,"_refmapped.csv.gz")),
    by.y = "cb", by.x = "cell_id"
  ) %>%filter(reads_all >= 5) -> df_int
  
  merge(
    df_int, 
    fread(paste0("../data/reference_projections/namdc/",x,"_refmapped.csv.gz")),
    by.y = "cb", by.x = "cell_id"
  ) %>%filter(reads_all >= 5) -> df
  df
}


df550 <- make_combined_df("NAMDC1412802550")
df671 <- make_combined_df("NAMDC20008038671") 
df507 <- make_combined_df("NAMDC18207019507")

cov550 <- assays(readRDS("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/NAMDC1412802550_hg38_mask_mgatk.rds"))[["coverage"]]
cov671 <- assays(readRDS("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/NAMDC20008038671_hg38_mask_mgatk.rds"))[["coverage"]]
cov507 <- assays(readRDS("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/NAMDC18207019507_hg38_mask_mgatk.rds"))[["coverage"]]

df550$coverage_heteroplasmy <- estimate_coverage_heteroplasmy_mat(8482,13447,cov550[,df550$cell_id]) 
df671$coverage_heteroplasmy <- estimate_coverage_heteroplasmy_mat(7982,15496,cov671[,df671$cell_id]) 
df507$coverage_heteroplasmy <- estimate_coverage_heteroplasmy_mat(8687,15103,cov507[,df507$cell_id]) 
df507$coverage_heteroplasmy2 <- estimate_coverage_heteroplasmy_mat(5924,6450,cov507[,df507$cell_id]) 

df507 %>% arrange(desc(coverage_heteroplasmy2)) %>% filter(reads_all > 10)

mdf_filt <- rbind(
  data.frame(df507[,c("heteroplasmy", "coverage_heteroplasmy")], donor = "d507"),
  data.frame(df550[,c("heteroplasmy", "coverage_heteroplasmy")], donor = "d550"),
  data.frame(df671[,c("heteroplasmy", "coverage_heteroplasmy")], donor = "d671")
  
)
qplot(df507$coverage_heteroplasmy, df507$coverage_heteroplasmy2)


mdf_filt %>% dplyr::filter(donor == "d507") -> aPD1
mdf_filt %>% dplyr::filter(donor == "d550") -> aPD2
mdf_filt %>% dplyr::filter(donor == "d671") -> aPD3

cor(aPD1$heteroplasmy, aPD1$coverage_heteroplasmy, method = "spearman")
cor(aPD2$heteroplasmy, aPD2$coverage_heteroplasmy,  method = "spearman")
cor(aPD3$heteroplasmy, aPD3$coverage_heteroplasmy, method = "spearman")

mdf_filt%>%
  ggplot(aes(x = heteroplasmy, y = coverage_heteroplasmy)) + 
  geom_point(size = 0.5) +
  facet_wrap(~donor) +
  pretty_plot(fontsize = 8) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) -> pG

cowplot::ggsave2(pG, file = "../output/compare_approaches_scatter_namdc.pdf", width = 6, height = 2)


