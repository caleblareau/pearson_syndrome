library(dplyr)
library(BuenColors)
load("../output/Erythroid_Person_edgeR_D6D12-DGE.rda")

ery_mdf <- rbind(ery_df6, ery_df12)
ery_mdf$Zstat <- qnorm((1E-250 + ery_mdf$PValue)/2) * sign(ery_mdf$logFC) * -1

ery_mdf%>% group_by(gene) %>%
  summarize(SZ = sum(Zstat)/sqrt(n()), MP = mean(Pearson_cpm), MH = mean(Healthy_cpm),
            mLogFC = mean(logFC)) %>%
  arrange(desc(SZ)) %>% mutate(rank = 1:n()) -> Rank_summarize_zscore_ery

df <- readRDS("../../pbmc_scrna/output/20Dec-PearsonRNAseq-diffGE-edgeR.rds")
df$Zstat <- qnorm((1E-250 + df$PValue)/2) * sign(df$logFC) * -1

df %>% group_by(gene) %>%
  summarize(SZ = sum(Zstat)/sqrt(n()), MP = mean(Pearson_cpm), MH = mean(Healthy_cpm), mLogFC = mean(logFC)) %>%
  arrange(desc(SZ)) %>% mutate(rank = 1:n()) -> Rank_summarize_zscore_pbmc

mdf_other <- merge(Rank_summarize_zscore_pbmc, Rank_summarize_zscore_ery, by = "gene") %>%
  dplyr::filter(MH.x > 10 | MP.x > 10 | MH.y > 10 | MP.y > 10)

mdf_other %>% dplyr::filter(mLogFC.x  < -0.2 & mLogFC.y < -0.2 & (MP.x/MH.x)<1 & (MP.y/MH.y)<1) %>% pull(gene) %>%
  data.frame() %>% write.table(quote = FALSE, row.names = FALSE)
mdf_other %>% dplyr::filter(mLogFC.x  > 0.2 & mLogFC.y  > 0.2 & (MP.x/MH.x)>1 & (MP.y/MH.y)>1)  %>% pull(gene) %>%
  data.frame() %>% write.table(quote = FALSE, row.names = FALSE)

ggplot(mdf_other, aes(x = signedTstat, y = SZ)) + 
  geom_point()
