library(clusterProfiler)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)

# Import differential gene expression
diff_df <- readRDS("../../../pearson_large_data_files/output/pbmc/11July2022-PearsonRNAseq-diffGE-edgeR.rds")


tl <- bitr(diff_df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = FALSE) 
tl <- tl[!duplicated(tl),]; name_vec <- tl[[2]]; names(name_vec) <- tl[[1]]
diff_df$gene_id <- name_vec[as.character(diff_df$gene)]
diff_df <- diff_df[complete.cases(diff_df),]

# Compute per-celltype enrichments
lapply(unique(diff_df$celltype), function(ct_go){
  print(ct_go)
  ct <- diff_df %>% dplyr::filter(celltype == ct_go)
  geneList <- ct[["logFC"]]*-1*log10(ct$PValue); names(geneList) <- ct[["gene_id"]]; geneList <- sort(geneList, decreasing = TRUE)

  kk <-  gseMKEGG(geneList = geneList, # 
                   organism = 'hsa', pvalueCutoff = 1, minGSSize = 10)
  kk <- data.frame(kk)
  kk$celltype <- ct_go
  kk
}) %>% rbindlist() %>% data.frame()  -> all_enrich
all_enrich

pathways_sig <- all_enrich %>% mutate(z = abs(qnorm(pvalue/2)) * sign(NES)) %>%
  group_by(Description) %>% summarize(SZ = sum(z)/sqrt(n())) %>%
  mutate(pval = 2*pnorm(-abs(SZ))) %>% mutate(padj = p.adjust(pval)) %>% dplyr::filter(padj < 0.01) %>% 
  mutate(log_padj = -1*log10(padj)) %>% arrange(desc(log_padj))

d_con <- c("Cytochrome bc1 complex", "Cytochrome c oxidase", "F-type ATPase", "Respiratory complex I",
            "Gluconeogenesis", "Glycolysis")
names(d_con) <- c("Cytochrome bc1 complex", "Cytochrome c oxidase", "F-type ATPase, eukaryotes", 
           "NADH dehydrogenase (ubiquinone) 1 beta subcomplex", 
           "Gluconeogenesis, oxaloacetate => fructose-6P",
           "Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate")

enrich_out <- all_enrich %>% dplyr::filter(Description %in% pathways_sig$Description) %>%
  mutate(d_easy = d_con[Description])
enrich_out$d_easy <- factor(enrich_out$d_easy, levels = rev(unname(d_con)))
pathways_sig$d_easy <- factor(d_con[pathways_sig$Description], levels = rev(unname(d_con)))

p1 <- ggplot(enrich_out, aes(x = celltype, y = d_easy, color = NES, size = -log10(pvalue))) +
  geom_point() + scale_color_gradientn(colors = jdb_palette("brewer_spectra"),  limits = c(-2,2),  oob = scales::squish) +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "", y = "")  +
  theme(legend.position = "bottom")

p2 <- ggplot(pathways_sig, aes(x = log_padj, y = d_easy, fill = SZ > 0)) +
  geom_bar(stat = "identity", color = "black") + 
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(2,8)]) +
  pretty_plot(fontsize = 8) + L_border()  + labs(x = "-log10 adjusted p-value", y = "", fill = "")  +
  theme(legend.position = "bottom") + scale_x_continuous(expand = c(0,0), limits = c(0,18)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

cowplot::ggsave2(cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(2.5, 1)), 
                 filename = "../plots/KEGG_enrichment.pdf", width = 7.5, height = 3)
