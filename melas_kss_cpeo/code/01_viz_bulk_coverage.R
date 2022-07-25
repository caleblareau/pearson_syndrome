library(data.table)
library(dplyr)
library(viridis)
library(BuenColors)
library(SummarizedExperiment)
library(Matrix)
library(RcppRoll)

cov550 <- assays(readRDS("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/NAMDC1412802550_hg38_mask_mgatk.rds"))[["coverage"]]
cov671 <- assays(readRDS("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/NAMDC20008038671_hg38_mask_mgatk.rds"))[["coverage"]]
cov507 <- assays(readRDS("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/NAMDC18207019507_hg38_mask_mgatk.rds"))[["coverage"]]

df <- data.frame(
  position = RcppRoll::roll_mean(1:16569, n = 5),
  KSS2 =  RcppRoll::roll_mean(rowMeans(cov550),n = 5),
  KSS1 = RcppRoll::roll_mean(rowMeans(cov671),n = 5),
  CPEO1 = RcppRoll::roll_mean(rowMeans(cov507),n = 5)
)

mdf <- reshape2::melt(df, id.vars = "position")
  
p1 <- ggplot(mdf, aes(x = position, y = value, color = variable)) +
  geom_line() + scale_y_log10(limits = c(5, 120)) +
  pretty_plot(fontsize = 7) + 
  L_border() +  theme(legend.position = "none")
cowplot::ggsave2(p1, filename = "../plots/NAMDC_coverage_pseudobulk.pdf", width = 2.3, height = 1)

#----------------------------------
# Setub for cells with heteroplasmy

make_combined_df <- function(x){
  hetdf <- fread(paste0("../data/namdc-metadata/",x,"_del.new.tsv")) %>%
    filter(version == "improved")
  merge(
    hetdf, 
    fread(paste0("../data/reference_projections/namdc/",x,"_refmapped.csv.gz")),
    by.y = "cb", by.x = "cell_id"
  ) %>%filter(reads_all >= 5) -> df
  df
}
df550 <- make_combined_df("NAMDC1412802550")
df671 <- make_combined_df("NAMDC20008038671") 
df507 <- make_combined_df("NAMDC18207019507")
sum(df550$reads_del)/sum(df550$reads_all)*100
sum(df671$reads_del)/sum(df671$reads_all)*100
sum(df507$reads_del)/sum(df507$reads_all)*100
df550 %>% filter(heteroplasmy > 0) %>% dim()
df671 %>% filter(heteroplasmy > 0) %>% dim()
df507 %>% filter(heteroplasmy > 0) %>% dim()

het_df <- data.frame(
  position = RcppRoll::roll_mean(1:16569, n = 10),
  KSS2 =  RcppRoll::roll_mean(rowMeans(cov550[,df550 %>% filter(heteroplasmy > 0) %>% pull(cell_id)]),n = 10),
  KSS1 = RcppRoll::roll_mean(rowMeans(cov671[,df671 %>% filter(heteroplasmy > 0) %>% pull(cell_id)]),n = 10),
  CPEO1 = RcppRoll::roll_mean(cov507[,df507 %>% filter(heteroplasmy > 0) %>% pull(cell_id)],n = 10)
)

p2 <- ggplot(reshape2::melt(het_df, id.vars = "position"), aes(x = position, y = value, color = variable)) +
  geom_line() + scale_y_log10(limits = c(5, 200)) +
  pretty_plot(fontsize = 7) + 
  facet_wrap(~variable,nrow = 3) +
  L_border() +  theme(legend.position = "none")
cowplot::ggsave2(p2, filename = "../plots/NAMDC_heteroplasmy_pseudobulk.pdf", width = 2.3, height = 1.8)

