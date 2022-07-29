library(edgeR)
library(Seurat)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(stringr)

"%ni%" <- Negate("%in%")

###########
melted_list <- readRDS("../output/24July2022-PearsonRNAseq-diffGE-edgeR.rds")
melted_list$pearson_adultFC <- log2(((melted_list$Pearson_cpm+1)/(rowMeans(data.matrix(data.frame(melted_list[,c("HA1", "HA2")])))+1)))
melted_list$pearson_pediatricFC <- log2(((melted_list$Pearson_cpm+1)/(rowMeans(data.matrix(data.frame(melted_list[,c("Hped1", "Hped2")])))+1)))

ggplot(melted_list, aes(pearson_adultFC, pearson_pediatricFC)) +
  geom_point()

melted_list$Zstat <- qnorm((1E-250 + melted_list$PValue)/2) * sign(melted_list$logFC) * -1
melted_list %>% group_by(gene) %>%
  summarize(SZ = sum(Zstat)/sqrt(n()), 
            mLogFCpediatric = mean(pearson_pediatricFC),
            mLogFCadult = mean(pearson_adultFC)
  ) %>%
  arrange(desc(SZ)) %>% mutate(rank = 1:n()) -> Rank_summarize_zscore

# compute a bonferroni adjusted pvalue
qnorm(0.01/15000)

table((Rank_summarize_zscore$SZ) >10 )
table((Rank_summarize_zscore$SZ) < -10 )

p1 <- ggplot(Rank_summarize_zscore, aes(x = rank, y = SZ)) + 
  geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Association rank", y = "Mean z-score")
p1

cowplot::ggsave2(p1, file = "../plots/go_zscore.pdf", width =1.62, height = 2)

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


Rank_summarize_zscore %>%
  filter(gene %in% oxphos_pathway) %>% write.table()

Rank_summarize_zscore %>%
  filter(gene %in% c("FOS", "JUN", "JUNB", "FOSB")) %>% write.table()


write.table(Rank_summarize_zscore[,c(1,2)], 
            file = "../output/pbmc-total-assoc.rnk",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

Rank_summarize_zscore %>% filter(gene %in% c("UQCRB"))

Rank_summarize_zscore %>% 
  filter(mLogFCpediatric > 3 & mLogFCadult > 3) %>% data.frame()

Rank_summarize_zscore %>% 
  filter(mLogFCpediatric < -2.5 & mLogFCadult < -2.5) %>% data.frame()
