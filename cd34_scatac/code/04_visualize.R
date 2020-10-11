library(data.table)
library(dplyr)
library(ggbeeswarm)
library(BuenColors)

pt3hetall <- fread("data/PT3_all_libraries.tsv") %>% data.frame()
pt3het <- pt3hetall %>% filter(cell_clust == "PBMC")
cov_df_pt3 <- fread("../Pearson_PBMC_CR-mtMask_mgatk/final/Pearson_PBMC_CR-mtMask_mgatk.depthTable.txt");  colnames(cov_df_pt3) <- c("cell_id", "coverage")
pt3 <- merge(merge(pt3het, cov_df_pt3), fread("projection_dfs/PT3_PBMC_final.tsv"), by.x = "cell_id", by.y = "barcode") %>% data.frame()
bci <- merge(fread("data/BCI_master.tsv"), fread("projection_dfs/BCI_PBMC_final.tsv"), by.x = "cell_id", by.y = "barcode")%>% data.frame()
ccf <- merge(fread("data/CCF_master.tsv"), fread("projection_dfs/CCF_PBMC_final.tsv"), by.x = "cell_id", by.y = "barcode")%>% data.frame()


table(pt3$cluster)
table(bci$cluster.y)
table(ccf$cluster.y)


pt3 %>% filter(coverage > 15) %>% group_by(cluster) %>%
  summarize(
    count =n(),
    mean_heteroplasmy = round(mean(corr_heteroplasmy), 2), 
    median_heteroplasmy = round(median(corr_heteroplasmy), 2)
  ) %>% data.frame() %>% write.table(sep = ",", row.names = FALSE, quote = FALSE)

bci %>% filter(coverage > 15) %>% group_by(cluster.y) %>%
  summarize(
    count =n(),
    mean_heteroplasmy = round(mean(heteroplasmy_corr), 2), 
    median_heteroplasmy = round(median(heteroplasmy_corr), 2)
  ) %>% data.frame() %>% write.table(sep = ",", row.names = FALSE, quote = FALSE)

ccf %>% filter(coverage > 15) %>% group_by(cluster.y) %>%
  summarize(
    count =n(),
    mean_heteroplasmy = round(mean(heteroplasmy_corr), 2), 
    median_heteroplasmy = round(median(heteroplasmy_corr), 2)
  ) %>% data.frame() %>% write.table(sep = ",", row.names = FALSE, quote = FALSE)


p1 <- ggplot(bci %>% filter(total > 15 ), aes(x = cluster.y, y = heteroplasmy_corr)) +
  geom_quasirandom() + geom_boxplot(color = "firebrick", fill = NA, outlier.shape = NA) +
  ggtitle("Patient BCI") + labs(x = "Annotated Celltype", y = "Deletion heteroplasmy") +
  pretty_plot() + L_border()

p2 <- ggplot(ccf %>% filter(total > 15 ), aes(x = cluster.y, y = heteroplasmy_corr)) +
  geom_quasirandom() + geom_boxplot(color = "firebrick", fill = NA, outlier.shape = NA) +
  ggtitle("Patient CCF") + labs(x = "Annotated Celltype", y = "Deletion heteroplasmy") +
  pretty_plot() + L_border()

p3 <- ggplot(pt3 %>% filter(coverage > 15 ), aes(x = cluster, y = corr_heteroplasmy)) +
  geom_quasirandom() + geom_boxplot(color = "firebrick", fill = NA, outlier.shape = NA) +
  ggtitle("Patient PT3 (Most Recent)") + labs(x = "Annotated Celltype", y = "Deletion heteroplasmy") +
  pretty_plot() + L_border()

cowplot::ggsave(cowplot::plot_grid(p1, p2, p3, nrow = 3), file = "distrbutions_classified.pdf", width = 10, height = 12)

