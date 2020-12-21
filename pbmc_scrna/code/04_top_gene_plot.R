library(data.table)
library(dplyr)
library(Seurat)
library(Matrix)
library(pheatmap)
library(matrixStats)

"%ni%" <- Negate("%in%")

# Pull data from previous steps
df <- readRDS("../../../pearson_mtscatac_large_data_files/output/20Dec-PearsonRNAseq-diffGE-edgeR.rds")
df$Zstat <- qnorm((1E-250 + df$PValue)/2) * sign(df$logFC) * -1

df %>% group_by(gene) %>%
  summarize(SZ = sum(Zstat)/sqrt(n())) %>%
  arrange(desc(SZ)) %>% mutate(rank = 1:n()) -> Rank_summarize_zscore

p1 <- ggplot(Rank_summarize_zscore, aes(x = rank, y = SZ)) + 
  geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Association rank", y = "Mean z-score")

p2 <- df[,c("gene", "celltype", "H1v2", "H2v3", "pBCI", "pCCF", "pPT3")] %>% filter(gene %in% c("JUN")) %>%
  reshape2::melt(id.vars = c("gene", "celltype")) %>% 
  ggplot(aes(x = variable, y = celltype, fill = log2(value))) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos")) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(limits = rev(levels(as.factor(df$celltype))), expand = c(0,0)) +
  labs(x = "", y = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") 


p3 <- df[,c("gene", "celltype", "H1v2", "H2v3", "pBCI", "pCCF", "pPT3")] %>% filter(gene %in% c("CXCL14")) %>%
  reshape2::melt(id.vars = c("gene", "celltype")) %>% 
  ggplot(aes(x = variable, y = celltype, fill = (value ))) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos")) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(limits = rev(levels(as.factor(df$celltype))), expand = c(0,0)) +
  labs(x = "", y = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") 

cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(0.9, 1.1, 1.1)), 
                 filename = "../plots/gene_analyses_top.pdf", width = 6, height = 1.9)
