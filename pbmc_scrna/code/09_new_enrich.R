library(fgsea)
library(msigdbr)
library(dplyr)
library(data.table)
library(BuenColors)
library(clusterProfiler)

# Import genes
diff_df <- readRDS("../../../pearson_large_data_files/output/pbmc/11July2022-PearsonRNAseq-diffGE-edgeR.rds")
organism = "org.Hs.eg.db"

celltype_go = "CD14.Mono"

lapply(unique(diff_df$celltype), function(celltype_go){
  print(celltype_go)
  ct_df <- diff_df %>% filter(celltype == celltype_go) 
  universe <- ct_df %>% pull(gene)
  genes_up <- ct_df %>% filter(FDR < 0.05 & logFC > 0.25) %>% pull(gene)
  genes_down <- ct_df %>% filter(FDR < 0.05 & logFC < -0.25) %>% pull(gene)
  or <- function(df){
    sapply(1:dim(df)[1], function(idx){
      fg <- as.numeric(strsplit(df$GeneRatio[idx], "/",2)[[1]])
      bg <- as.numeric(strsplit(df$BgRatio[idx], "/",2)[[1]])
      (fg[1]/fg[2])/(bg[1]/bg[2])
    })
  }
  go_enrich_up <- enrichGO(gene = genes_up,
                           universe = universe,
                           OrgDb = organism, 
                           keyType = 'SYMBOL',
                           ont = "ALL",
                           qvalueCutoff = 0.05) %>% data.frame()
  go_enrich_up$logFC <- log2(or(go_enrich_up))
  go_enrich_down <- enrichGO(gene = genes_down,
                             universe = universe,
                             OrgDb = organism, 
                             keyType = 'SYMBOL',
                             ont = "ALL",
                             qvalueCutoff = 0.05)%>% data.frame()
  go_enrich_down$logFC <- -1*log2(or(go_enrich_down))
  
  both <- rbind(
    go_enrich_up[,c(1,2,3,10,11,8)] %>% 
      mutate(celltype = celltype_go, what = "up"),
    go_enrich_down[,c(1,2,3,10,11,8)]  %>% 
      mutate(celltype = celltype_go, what = "down")
  ) %>% arrange(qvalue)
  
}) %>% rbindlist() %>% data.frame() -> all_enrich
top_enrich <- all_enrich %>% arrange(qvalue) %>% filter(ONTOLOGY == "BP")
top_enrich[grepl("electron", top_enrich[["Description"]]),]
top_enrich[grepl("lipid", top_enrich[["Description"]]),]

top_enrich %>% filter(what == "up") %>% group_by(celltype) %>% slice_min(qvalue, n = 2) %>% data.frame() 
top_enrich %>% filter(what == "down") %>% group_by(celltype) %>% slice_min(qvalue, n = 2) %>% data.frame() 

top_enrich %>% filter(Count > 10) %>% arrange(desc(logFC*-log10(qvalue)))%>%
  mutate(rank = 1:n()) %>%
  ggplot(aes(x = rank, y = logFC * (-1)*log10(qvalue))) + geom_point()

         