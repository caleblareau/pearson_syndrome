library(data.table)
library(dplyr)
library(Seurat)
library(Matrix)
library(pheatmap)
library(matrixStats)

"%ni%" <- Negate("%in%")

# Pull data from previous steps
df <- readRDS("../../../pearson_large_data_files/output/20Dec-PearsonRNAseq-diffGE-edgeR.rds")
df$Zstat <- qnorm((1E-250 + df$PValue)/2) * sign(df$logFC) * -1
df %>% mutate(rat = (Pearson_cpm + 0.01) / (Healthy_cpm + 0.01)) %>% filter(rat > 100) %>% data.frame()

df %>% group_by(gene) %>%
  summarize(SZ = sum(Zstat)/sqrt(n()), MP = mean(Pearson_cpm), MH = mean(Healthy_cpm)) %>%
  arrange(desc(SZ)) %>% mutate(rank = 1:n()) -> Rank_summarize_zscore

write.table(Rank_summarize_zscore[,c(1,2)], file = "../output/Rank_Zscore_combined.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
p1 <- ggplot(Rank_summarize_zscore, aes(x = rank, y = SZ)) + 
  geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Association rank", y = "Mean z-score")

Rank_summarize_zscore %>% dplyr::filter(gene %in% c("C12orf54", "CXCL14", "LINC01641"))
