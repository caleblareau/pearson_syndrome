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

melted_list %>%
  filter(gene %in% c("UQCRB"))

write.table(Rank_summarize_zscore[,c(1,2)], 
            file = "../output/pbmc-total-assoc.rnk",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

Rank_summarize_zscore %>% filter(gene %in% c("UQCRB"))

Rank_summarize_zscore %>% 
  filter(mLogFCpediatric > 3 & mLogFCadult > 3) %>% data.frame()

Rank_summarize_zscore %>% 
  filter(mLogFCpediatric < -2.5 & mLogFCadult < -2.5) %>% data.frame()
