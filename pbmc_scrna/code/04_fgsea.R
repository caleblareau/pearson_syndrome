library(fgsea)
library(msigdbr)
library(dplyr)
library(data.table)
library(BuenColors)

# Prepare gene set
m_df = msigdbr(species = "Homo sapiens", "H")
pathways <- split(m_df[["gene_symbol"]], m_df[["gs_name"]])
pathways[["AP1_TARGETS"]] <- fread("../data/misc/AP1_targets_Harmonizome.txt", header = FALSE)[[1]]

# Import genes
diff_df <- readRDS("../../../pearson_mtscatac_large_data_files/output/5March-PearsonRNAseq-diffGE-edgeR.rds")

lapply(unique(diff_df$celltype), function(ct){
  ct_df <- diff_df %>% filter(celltype == ct)
  vec <- ct_df$signedTstat

  names(vec) <- as.character(ct_df$gene)
  
  fgsea_out <- fgseaMultilevel(pathways = pathways, 
                     stats = sort(vec, decreasing = TRUE),
                     minSize=15,
                     maxSize=2000)
  if(FALSE){
    topPathways <- fgsea_out %>% filter(padj < 0.01) %>% pull(pathway)
    plotGseaTable(pathways[topPathways], vec, fgsea_out, 
                  gseaParam = 0.5)
    ct_df %>% filter(gene %in% pathways[["HALLMARK_HYPOXIA"]]) %>% head(100)
  }
  fgsea_out[,c("pathway", "pval", "padj", "NES")] %>% arrange(pval) %>% mutate(celltype = ct)
}) %>% rbindlist() %>% data.frame() -> enrich_out_all

enrich_out_all %>% filter(padj < 0.01) %>% arrange(NES) %>% tail(20)
enrich_out_all %>% group_by(celltype) %>% summarize(n = sum(padj < 0.1))

# Filter such that a pwayway has to be significant in at least 1 tissue
pathways_sig <- enrich_out_all %>% filter(padj < 0.01 ) %>% pull(pathway) %>% unique
enrich_out <- enrich_out_all %>% filter(pathway %in% pathways_sig)

# Rename and plot
ptwys <- gsub("HALLMARK_", "", as.character(enrich_out$pathway))
ptwy_levels <- rev(unique(sort(ptwys)))
enrich_out$pathway_order <- factor(ptwys, levels = ptwy_levels)
p1 <- ggplot(enrich_out, aes(x = celltype, y = pathway_order, color = NES, size = -log10(padj))) +
  geom_point() + scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  pretty_plot() + L_border() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "", y = "") 
ggsave(p1, file = "../plots/06March_GSEA.pdf", width = 7, height = 4)
