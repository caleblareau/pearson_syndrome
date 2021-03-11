library(dplyr)
library(BuenColors)
load("../output/Erythroid_Person_edgeR_D6D12-DGE.rda")

ery_mdf <- merge(ery_df6, ery_df12, by = "row.names") 

ery_mdf %>% dplyr::filter(logFC.x < -1.3 | logFC.y < -1.3)

df <- readRDS("../../../pearson_large_data_files/output/20Dec-PearsonRNAseq-diffGE-edgeR.rds")
df$Zstat <- qnorm((1E-250 + df$PValue)/2) * sign(df$logFC) * -1

df %>% group_by(gene) %>%
  summarize(SZ = sum(Zstat)/sqrt(n()), MP = mean(Pearson_cpm), MH = mean(Healthy_cpm)) %>%
  arrange(desc(SZ)) %>% mutate(rank = 1:n()) -> Rank_summarize_zscore_pbmc

mdf_other <- merge(Rank_summarize_zscore_pbmc, ery_df12, by = "gene") %>%
  dplyr::filter(MH > 10 | MP > 10) %>% dplyr::filter(Pearson_cpm > 10 | Healthy_cpm > 10)

mdf_other %>% dplyr::filter(signedTstat < -10 & SZ < -10) %>% dim()
mdf_other %>% dplyr::filter(signedTstat  > 10 & SZ > 10) %>% dim()
mdf_other %>% dplyr::filter(signedTstat > 10 & SZ < -10)%>% dim()
mdf_other %>% dplyr::filter(signedTstat  < -10 & SZ > 10) %>% dim()

ggplot(mdf_other, aes(x = signedTstat, y = SZ)) + 
  geom_point()
