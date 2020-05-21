library(fgsea)
library(msigdbr)
library(dplyr)
library(data.table)
library(BuenColors)

# Prepare gene set
ap1_genes <- fread("../misc/AP1_targets_Harmonizome.txt", header = FALSE)[[1]]
pathways <- list(ap1_genes)
names(pathways) <- c("AP1targets")

# Import genes
diff_df <- readRDS("../../../pearson_mtscatac_large_data_files/output/5March-PearsonRNAseq-diffGE-edgeR.rds")

lapply(unique(diff_df$celltype), function(ct){
  print(ct)
  ct_df <- diff_df %>% filter(celltype == ct)
  vec <- ct_df$logFC
  
  names(vec) <- as.character(ct_df$gene)
  
  fgsea_out <- fgseaMultilevel(pathways = pathways, 
                               stats = sort(vec, decreasing = TRUE),
                               minSize=15,
                               maxSize=2000)
 
  fgsea_out[,c("pathway", "pval", "padj", "NES")] %>% arrange(pval) %>% mutate(celltype = ct)
}) %>% rbindlist() %>% data.frame() -> enrich_out_all

enrich_out_all %>%
  mutate(log10padj = -1 * log10(padj)) %>%
  arrange(desc(log10padj)) -> plot_df

plot_bar <- ggplot(plot_df, aes(y = celltype, x = log10padj)) +
  geom_bar(stat = "identity", fill = jdb_palette("brewer_spectra")[2], color = "black") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_discrete(limits = rev((unique(plot_df$celltype))), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) + 
  theme(legend.position = "bottom") + geom_vline(xintercept = -1*log10(0.05), linetype = 2)

cowplot::ggsave2(plot_bar, file = "../plots/AP1_gsea.pdf", width = 2.5, height = 2)
