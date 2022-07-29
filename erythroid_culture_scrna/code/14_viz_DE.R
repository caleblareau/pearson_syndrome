library(dplyr)
library(BuenColors)

# Gene sets
heme <- c("ALAD","COX10","CPOX","EARS2","EPRS","FECH","HMBS","PPOX","QARS","UROD")
cholesterol <- c("FDFT1","FDPS","GGPS1","HMGCR","HMGCS1","IDI2","LSS","MVD","MVK","PMVK","SQLE")

# Get scores
oxphos_pathway <- c("SNORD138","MIR4691","COX17","MIR7113","ATP5PD","ATP5MG","UQCR11","COX4I1",
                    "COX5B","COX6A1","COX6A2","COX6B1","COX6C","COX7A1","COX7A2","COX7B","COX7C",
                    "COX8A","COX11","COX15","UQCRQ","DMAC2L","SLC25A4","SLC25A5","SLC25A6","UQCR10",
                    "NDUFS7","MT-ATP6","MT-ATP8","MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ND1",
                    "MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT ND6","NDUFA1","NDUFA2","NDUFA3",
                    "NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFAB1","NDUFB1",
                    "NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
                    "NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFV1","NDUFS4","NDUFS5","NDUFS6",
                    "NDUFS8","NDUFV2","NDUFV3","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E",
                    "ATP5PB","ATP5MC1","ATP5MC2","ATP5MC3","ATP5ME","ATP5PF","ATP5PO","NDUFA12",
                    "SCO1","SDHA","SDHB","SDHC","SDHD","SURF1","UCP1","UCP2","UCP3","UQCRB","UQCRC1",
                    "UQCRC2","UQCRFS1","UQCRH","SLC25A14","COX7A2L","COX5A","ATP5IF1","SLC25A27","ATP5MF")


load("../output/Erythroid_Person_edgeR_D6D12-DGE.rda")

ery_mdf <- rbind(ery_df6, ery_df12)
ery_mdf$Zstat <- qnorm((1E-250 + ery_mdf$PValue)/2) * sign(ery_mdf$logFC) * -1

ery_mdf%>% group_by(gene) %>%
  summarize(SZ = sum(Zstat)/sqrt(n()), MP = mean(Pearson_cpm), MH = mean(Healthy_cpm),
            mLogFC = mean(logFC)) %>%
  arrange(desc(SZ)) %>% mutate(rank = 1:n()) -> Rank_summarize_zscore_ery


df <- readRDS("../../pbmc_scrna/output/24July2022-PearsonRNAseq-diffGE-edgeR.rds")
df$Zstat <- qnorm((1E-250 + df$PValue)/2) * sign(df$logFC) * -1

df %>% group_by(gene) %>%
  summarize(SZ = sum(Zstat)/sqrt(n()), MP = mean(Pearson_cpm), MH = mean(Healthy_cpm), mLogFC = mean(logFC)) %>%
  arrange(desc(SZ)) %>% mutate(rank = 1:n()) -> Rank_summarize_zscore_pbmc


p1 <- ggplot(Rank_summarize_zscore_ery, aes(x = rank, y = SZ)) + 
  geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Association rank", y = "Mean z-score")
p1

cowplot::ggsave2(p1, file = "../plots/go_zscore.pdf", width =1.62, height = 2)

mdf_other <- merge(Rank_summarize_zscore_pbmc, Rank_summarize_zscore_ery, by = "gene") %>%
  dplyr::filter(MH.x > 10 | MP.x > 10 | MH.y > 10 | MP.y > 10)

mdf_other %>% dplyr::filter(mLogFC.x  < -0.2 & mLogFC.y < -0.2 & (MP.x/MH.x)<1 & (MP.y/MH.y)<1) %>% pull(gene) %>%
  data.frame() %>% write.table(quote = FALSE, row.names = FALSE)
mdf_other %>% dplyr::filter(mLogFC.x  > 0.2 & mLogFC.y  > 0.2 & (MP.x/MH.x)>1 & (MP.y/MH.y)>1)  %>% pull(gene) %>%
  data.frame() %>% write.table(quote = FALSE, row.names = FALSE)

ggplot(mdf_other, aes(x = SZ.x, y = SZ.y)) + 
  geom_point()
