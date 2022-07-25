library(data.table)
library(dplyr)
library(viridis)
library(BuenColors)

# Import Melas for one comparison
make_combined_df <- function(x){
  hetdf <- fread(paste0("../../melas_kss_cpeo/data/melas-metadata/",x,"_meta.tsv"))
  merge(
    hetdf, 
    fread(paste0("../../melas_kss_cpeo/data/reference_projections/melas//",x,"_refmapped.csv.gz")),
    by.y = "cb", by.x = "barcode"
  ) 
  
}
dfp9 <- make_combined_df("P9") 
dfp21 <- make_combined_df("P21")
dfp30 <- make_combined_df("P30")

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


prop0compute <- function(ssdf){
  sum(ssdf$heteroplasmy < 0.1) / dim(ssdf)[1]
 # ks.test(ssdf$heteroplasmy,ssdf$permuted_heteroplasmy)$p.value
  
}

ct <- unique(c(df_BCI$predicted.celltype.l2, df_CCF$predicted.celltype.l2, df_PT3$predicted.celltype.l2))
#ct <- c("MAIT", "NK", "B naive", "CD8 TEM", "CD4 TCM")
donors <- c("PT1", "PT2", "PT3", "xM9", "yM21",  "zM30")
donor_df_list <- list(df_BCI, df_CCF, df_PT3, dfp9, dfp21, dfp30)
lapply(ct, function(ctx){
  lapply(1:6, function(idx){
    subset_df1 <- donor_df_list[[idx]] %>% filter(predicted.celltype.l2 == ctx)
    data.frame(
      donor = donors[idx],
      celltype = ctx,
      n = dim(subset_df1)[1],
      prop0 = prop0compute(subset_df1)*100
    )
  }) %>% rbindlist() 
})%>% rbindlist() -> odf

df_BCI%>% group_by(predicted.celltype.l1, predicted.celltype.l2) %>%
  summarize(count =n()) %>% top_n(1,count)


# Define T cell colors
cells <- c("CD4 Naive", "CD4 TCM", "CD8 Naive", "CD8 TEM", "MAIT")
colors <- c("dodgerblue", "dodgerblue4", "pink", "firebrick", "purple4")

# Define other lineage colors
cells <- c(cells, c("CD14 Mono", "CD16 Mono", "cDC2", "pDC", "B naive", "B intermediate", "B memory", "NK", "NK_CD56bright"))
colors <- c(colors, rep("darkgrey", 4), rep("green4", 3), rep("orange2", 2))
names(colors) <- cells

px <- odf %>% arrange(desc(prop0)) %>%
  filter(celltype %in% names(colors)) %>%
  group_by(donor) %>% mutate(rank = 1:n()) %>%
  ggplot(aes(x = rank, y = prop0, color = celltype)) +
  geom_point() + facet_wrap(~donor) + scale_color_manual(values = colors) +
  pretty_plot(fontsize = 8) +
  labs(x = "Rank ordered cell types", y = "% cells with 0% heteroplasmy")
cowplot::ggsave2(px, file = "../plots/prop0_cells.pdf", width = 7, height = 3.6)
