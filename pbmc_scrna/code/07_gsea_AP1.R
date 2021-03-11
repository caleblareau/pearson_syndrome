library(fgsea)
library(msigdbr)
library(dplyr)
library(data.table)
library(BuenColors)

# Prepare gene set
ap1_genes <- fread("../data/misc/AP1_targets_Harmonizome.txt", header = FALSE)[[1]]
s <- fread("../data/misc/senescence_up_genes_C2_fridman.tsv", header = FALSE)[[1]]
g <- fread("../../../glycolysis.txt", header = FALSE)[[1]]
o <- fread("../../../oxphos.txt", header = FALSE)[[1]]

pathways <- list(ap1_genes, s, g, o)
names(pathways) <- c("AP1targets", "senescence", "gly", "oxphos")

# Import genes
diff_df <- readRDS("../output/20Dec-PearsonRNAseq-diffGE-edgeR.rds")

lapply(unique(diff_df$celltype), function(ct){
  print(ct)
  ct_df <- diff_df %>% filter(celltype == ct)
  vec <- ct_df$logFC
  
  names(vec) <- as.character(ct_df$gene)
  
  fgsea_out <- fgseaMultilevel(pathways = pathways, 
                               stats = sort(vec, decreasing = TRUE),
                               minSize=15,
                               maxSize=2000)
  plotGseaTable(pathways = pathways, 
                stats = sort(vec, decreasing = TRUE),fgsea_out)
  
  fgsea_out[,c("pathway", "pval", "padj", "NES")] %>% arrange(pval) %>% mutate(celltype = ct)
}) %>% rbindlist() %>% data.frame() -> enrich_out_all

enrich_out_all %>%dplyr::filter(pathway == "gly")
enrich_out_all %>%dplyr::filter(pathway == "oxphos")



enrich_out_all %>%dplyr::filter(pathway == "senescence") %>%
  mutate(log10padj = -1 * log10(padj)) %>%
  arrange(desc(log10padj)) -> pb_sens

ggplot(pb_sens, aes(y = celltype, x = log10padj)) +
  geom_bar(stat = "identity", fill = jdb_palette("brewer_spectra")[2], color = "black") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_discrete(limits = rev((unique(pb_sens$celltype))), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) + 
  theme(legend.position = "bottom") + 
  geom_vline(xintercept = -1*log10(0.05), linetype = 2) -> plot_bar_senescence
plot_bar_senescence
enrich_out_all %>%dplyr::filter(pathway == "AP1targets") %>%
  mutate(log10padj = -1 * log10(padj)) %>%
  arrange(desc(log10padj)) -> pb_ap1 

ggplot(pb_ap1, aes(y = celltype, x = log10padj)) +
  geom_bar(stat = "identity", fill = jdb_palette("brewer_spectra")[2], color = "black") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_discrete(limits = rev((unique(pb_ap1$celltype))), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) + 
  theme(legend.position = "bottom") + 
  geom_vline(xintercept = -1*log10(0.05), linetype = 2) -> plot_bar_AP1targets


cowplot::ggsave2(cowplot::plot_grid(plot_bar_AP1targets, plot_bar_senescence, nrow = 1),
                 filename  = "../plots/sens_AP1_gsea.pdf", width = 5, height = 2)
