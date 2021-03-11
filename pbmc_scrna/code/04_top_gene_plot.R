library(data.table)
library(dplyr)
library(Seurat)
library(Matrix)
library(pheatmap)
library(matrixStats)

"%ni%" <- Negate("%in%")

# Pull data from previous steps
df <- readRDS("../output/20Dec-PearsonRNAseq-diffGE-edgeR.rds")
df$Zstat <- qnorm((1E-250 + df$PValue)/2) * sign(df$logFC) * -1
df %>% mutate(rat = (Pearson_cpm + 0.01) / (Healthy_cpm + 0.01)) %>% filter(rat > 100) %>% data.frame()

df %>% group_by(gene) %>%
  summarize(SZ = sum(Zstat)/sqrt(n()), MP = mean(Pearson_cpm), MH = mean(Healthy_cpm)) %>%
  arrange(desc(SZ)) %>% mutate(rank = 1:n()) -> Rank_summarize_zscore

metagenes <- c("SLC2A1", "SLC2A2", "SLC2A3", "SLC2A4","SLC2A5", "HK1", "HK2", "HK3", "GCK", "G6PC", "GPI", 
  "PFKM","PFKL", "FBP1", "FBP2", "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "GAPDHS", "PGK1", "PGK2", "PGAM1", "PGAM2",
  "ENO1", "ENO2", "ENO3", "PKM2", "PKLR", "LDHA", "LDHB", "DHC", "LDHAL6B", "PCK1", "GOT1", "MDH1", "MPC1", "MPC2", 
  "PC", "GOT2", "MDH2", "PDHA1", "PDHA2", "PDHB", "DLAT", "DLD", "PDHX")
Rank_summarize_zscore %>% dplyr::filter(gene %in% metagenes) %>% data.frame()


write.table(Rank_summarize_zscore[,c(1,2)], file = "../output/Rank_Zscore_combined.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
p1 <- ggplot(Rank_summarize_zscore, aes(x = rank, y = SZ)) + 
  geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Association rank", y = "Mean z-score")

Rank_summarize_zscore %>% dplyr::filter(gene %in% c("C12orf54", "CXCL14", "LINC01641"))
