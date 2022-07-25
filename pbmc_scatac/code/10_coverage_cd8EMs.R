library(data.table)
library(dplyr)
library(viridis)
library(BuenColors)
library(SummarizedExperiment)
library(Matrix)

# Add Heteroplasmy
df_BCI <- merge(fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
                  dplyr::filter(deletion == "del6073-13095") ,
                fread(paste0("../data/reference_projection/BCI_refmapped.csv.gz")),
                by.y = "cb", by.x = "cell_id") %>%filter(reads_all >= 10) 

df_CCF <- merge(fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
                  dplyr::filter(deletion == "del8482-13447") ,
                fread(paste0("../data/reference_projection/CCF_refmapped.csv.gz")),
                by.y = "cb", by.x = "cell_id") %>%filter(reads_all >= 10) 

df_PT3 <- merge(fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
                  dplyr::filter(deletion == "del10381-15407") ,
                fread(paste0("../data/reference_projection/PT3_refmapped.csv.gz")),
                by.y = "cb", by.x = "cell_id") %>%filter(reads_all >= 10) 

# append permuted heteroplasmy
add_permuted_heteroplasmy <- function(df, celltype, donor){
  df <- df %>% filter(predicted.celltype.l2 == celltype)
  set.seed(1)
  df$permuted_heteroplasmy <- rbinom(n = length(df$reads_all), size = df$reads_all, prob = sum(df$reads_del)/sum(df$reads_all))/df$reads_all*100
  reshape2::melt(df[,c("heteroplasmy", "permuted_heteroplasmy")]) %>% 
    mutate(donor, celltype)
}

rdf <- rbind(
  add_permuted_heteroplasmy(df_BCI, "CD8 Naive", "PT1"),
  add_permuted_heteroplasmy(df_BCI, "CD8 TEM", "PT1"),
  add_permuted_heteroplasmy(df_BCI, "MAIT", "PT1"),
  add_permuted_heteroplasmy(df_CCF, "CD8 Naive", "PT2"),
  add_permuted_heteroplasmy(df_CCF, "CD8 TEM", "PT2"),
  add_permuted_heteroplasmy(df_CCF, "MAIT", "PT2"),
  add_permuted_heteroplasmy(df_PT3, "CD8 Naive", "PT3"),
  add_permuted_heteroplasmy(df_PT3, "CD8 TEM", "PT3"),
  add_permuted_heteroplasmy(df_PT3, "MAIT", "PT3")
)

library(ggbeeswarm)
p1 <- ggplot(rdf, aes(x = celltype, y = value, color = variable)) +
  geom_point(size = 0.1, position= position_jitterdodge()) + 
  facet_wrap(~donor) +
  scale_color_manual(values = c("black", "darkgrey")) +
  geom_violin(fill = NA,aes(color = variable),scale = "width") + 
  pretty_plot(fontsize = 7)  + theme(legend.position = "none") + 
  labs( x= "Celltype", y = "Heteroplasmy (%)")
cowplot::ggsave2(p1, file = "../plots/permuted_Tcell_coverage.pdf", width = 6.5, height = 2)


