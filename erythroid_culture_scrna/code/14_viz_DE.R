library(dplyr)
library(BuenColors)
load("../output/Erythroid_Person_edgeR_D6D12-DGE.rda")

ery_mdf <- merge(ery_df6, ery_df12, by = "row.names") 


metagenes <- c("SLC2A1", "SLC2A2", "SLC2A3", "SLC2A4","SLC2A5", "HK1", "HK2", "HK3", "GCK", "G6PC", "GPI", 
               "PFKM","PFKL", "FBP1", "FBP2", "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "GAPDHS", "PGK1", "PGK2", "PGAM1", "PGAM2",
               "ENO1", "ENO2", "ENO3", "PKM2", "PKLR", "LDHA", "LDHB", "DHC", "LDHAL6B", "PCK1", "GOT1", "MDH1", "MPC1", "MPC2", 
               "PC", "GOT2", "MDH2", "PDHA1", "PDHA2", "PDHB", "DLAT", "DLD", "PDHX")
ery_mdf %>% dplyr::filter(Row.names %in% metagenes) %>% data.frame() %>% arrange(desc(signedTstat.y))
  ggplot(aes(x = signedTstat.x , y = signedTstat.y)) + geom_point()

ery_mdf %>% dplyr::filter(logFC.x < -1.3 | logFC.y < -1.3)

df <- readRDS("../../pbmc_scrna/output/20Dec-PearsonRNAseq-diffGE-edgeR.rds")
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
