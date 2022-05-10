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
