library(dplyr)
library(BuenColors)
load("../output/Erythroid_Person_edgeR_D6D12-DGE.rda")

ery_mdf <- rbind(ery_df6, ery_df12)
ery_mdf$Zstat <- qnorm((1E-250 + ery_mdf$PValue)/2) * sign(ery_mdf$logFC) * -1

ery_mdf%>% group_by(gene) %>%
  summarize(SZ = sum(Zstat)/sqrt(n()), MP = mean(Pearson_cpm), MH = mean(Healthy_cpm),
            mLogFC = mean(logFC)) %>%
  arrange(desc(SZ)) %>% mutate(rank = 1:n()) -> Rank_summarize_zscore_ery
Rank_summarize_zscore_ery %>% dplyr::filter(gene %in% c("PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2"))
Rank_summarize_zscore_ery %>% dplyr::filter(gene %in% c("PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2"))


df <- readRDS("../../pbmc_scrna/output/20Dec-PearsonRNAseq-diffGE-edgeR.rds")
df$Zstat <- qnorm((1E-250 + df$PValue)/2) * sign(df$logFC) * -1

df %>% group_by(gene) %>%
  summarize(SZ = sum(Zstat)/sqrt(n()), MP = mean(Pearson_cpm), MH = mean(Healthy_cpm), mLogFC = mean(logFC)) %>%
  arrange(desc(SZ)) %>% mutate(rank = 1:n()) -> Rank_summarize_zscore_pbmc


p1 <- ggplot(Rank_summarize_zscore_ery, aes(x = rank, y = SZ)) + 
  geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Association rank", y = "Mean z-score")
cowplot::ggsave2(p1, file = "../plots/go_zscore.pdf", width =1.62, height = 2)

mdf_other <- merge(Rank_summarize_zscore_pbmc, Rank_summarize_zscore_ery, by = "gene") %>%
  dplyr::filter(MH.x > 10 | MP.x > 10 | MH.y > 10 | MP.y > 10)

mdf_other %>% dplyr::filter(mLogFC.x  < -0.2 & mLogFC.y < -0.2 & (MP.x/MH.x)<1 & (MP.y/MH.y)<1) %>% pull(gene) %>%
  data.frame() %>% write.table(quote = FALSE, row.names = FALSE)
mdf_other %>% dplyr::filter(mLogFC.x  > 0.2 & mLogFC.y  > 0.2 & (MP.x/MH.x)>1 & (MP.y/MH.y)>1)  %>% pull(gene) %>%
  data.frame() %>% write.table(quote = FALSE, row.names = FALSE)

ggplot(mdf_other, aes(x = signedTstat, y = SZ)) + 
  geom_point()
