library(dplyr)
library(BuenColors)

# Gene sets
heme <- c("ALAD","COX10","CPOX","EARS2","EPRS","FECH","HMBS","PPOX","QARS","UROD")
cholesterol <- c("FDFT1","FDPS","GGPS1","HMGCR","HMGCS1","IDI2","LSS","MVD","MVK","PMVK","SQLE")
serine_glycine <- c("PHGDH", "PSAT1", "PSPH", "SHMT2")

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


####################

mdf <- merge(Rank_summarize_zscore_pbmc[,c(1,2)], Rank_summarize_zscore_ery[,c(1,2)], by = "gene", all = TRUE) 
colnames(mdf) <- c("gene", "Zscore_PBMC", "Zscore_ery")
mdf$Zscore_PBMC <- ifelse(is.na(mdf$Zscore_PBMC), 0, mdf$Zscore_PBMC)
mdf$Zscore_ery <- ifelse(is.na(mdf$Zscore_ery), 0, mdf$Zscore_ery)
cor(data.matrix((mdf %>% filter(gene %in% oxphos_pathway)))[,-1])
cor(data.matrix((mdf %>% filter(gene %in% c(heme, cholesterol, serine_glycine))))[,-1])

ggplot(mdf %>% filter(gene %in% oxphos_pathway), aes(x = Zscore_PBMC, y = Zscore_ery)) + 
  geom_point(size = 0.6) +
  geom_smooth(method='lm', formula= y~x) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "PBMC z-score", y = "Erythroid z-score") -> p1 
cowplot::ggsave2(p1, file = "../plots/oxphos.pdf", width = 2.1, height = 2.1)

ggplot(mdf %>% filter(gene %in% c(heme, cholesterol,serine_glycine)) %>%
         mutate(inheme = (gene %in% heme)*2 + (gene %in% serine_glycine)*1),
                aes(x = Zscore_PBMC, y = Zscore_ery, shape = as.character(inheme))) +   geom_point(size = 1) +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") +
  labs(x = "PBMC z-score", y = "Erythroid z-score") -> p2
cowplot::ggsave2(p2, file = "../plots/ery.pdf", width = 2.1, height = 2.1)

ggplot(mdf %>% filter(gene %in% c(heme, cholesterol,serine_glycine)) %>%
         mutate(inheme = (gene %in% heme)*2 + (gene %in% serine_glycine)*1),
                aes(x = Zscore_PBMC, y = Zscore_ery)) + 
         geom_point(size = 0.6) +
         geom_smooth(method='lm', formula= y~x) +
         
         pretty_plot(fontsize = 8) + L_border() +
         theme(legend.position = "none") +
         labs(x = "PBMC z-score", y = "Erythroid z-score") -> p21
       cowplot::ggsave2(p21, file = "../plots/ery2.pdf", width = 2.1, height = 2.1)
       
       