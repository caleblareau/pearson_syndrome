library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)

all_mito_genes <- c("MT-ND1","MT-ND2","MT-CO1","MT-CO2","MT-ATP8","MT-ATP6" ,"MT-CO3" , "MT-ND3" , "MT-ND4L" ,"MT-ND4" , "MT-ND5",  "MT-ND6"  ,"MT-CYB")
bci <- Read10X("../data/pearson_bci/")
ccf <- Read10X("../data/pearson_ccf/")
pt3 <- Read10X("../data/pearson_mds/")
h1 <- Read10X("../data/healthy_pbmc_8k_v2-remap/")

# First look at total umis
import_anno_vec <- function(dir_base, short_id){
  idf <- readRDS(paste0("../output/seurat_projected_meta_", dir_base, ".rds"))
  vec <- idf$predicted.celltype.l1; names(vec) <- rownames(idf)
  (vec[vec %in% c("B", "CD4 T", "CD8 T", "Mono", "NK")])
}
anno_bci <- import_anno_vec("pearson_bci")
anno_ccf <- import_anno_vec("pearson_ccf")
anno_mds <- import_anno_vec("pearson_mds")
anno_h1 <- import_anno_vec("healthy_pbmc_8k_v2-remap")

rdf <- rbind(
  data.frame(pct_mt_UMIs = colSums(bci[all_mito_genes,])/colSums(bci)*100,celltype = anno_bci[colnames(bci)],donor = "PT1"),
  data.frame(pct_mt_UMIs = colSums(ccf[all_mito_genes,])/colSums(ccf)*100,celltype = anno_ccf[colnames(ccf)],donor = "PT2"),
  data.frame(pct_mt_UMIs = colSums(pt3[all_mito_genes,])/colSums(pt3)*100,celltype = anno_mds[colnames(pt3)],donor = "PT3"),
  data.frame(pct_mt_UMIs = colSums(h1[all_mito_genes,])/colSums(h1)*100,celltype = anno_h1[colnames(h1)],donor = "H1")
)
rdf <- rdf[complete.cases(rdf),]

rdf %>%
  ggplot(aes(x = celltype, color = donor, y = pct_mt_UMIs)) + 
  geom_boxplot()

# Subset to look at only deletion pct umis
genes_in_pt1 <- c("MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5")
genes_in_pt2 <- c("MT-ATP8", "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5")
genes_in_pt3 <- c("MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB")

make_pct_umis <- function(mat, genes_in_deletion, donor, deletion, anno_vec){
  data.frame(
    pct_mt_in_del_mito = colSums(mat[genes_in_deletion,])/colSums(mat[all_mito_genes,])*100,
    pct_umis_in_del = colSums(mat[genes_in_deletion,])/colSums(mat)*100,
    pct_umis_not_in_del = colSums(mat[all_mito_genes[!(all_mito_genes %in% genes_in_deletion)],])/colSums(mat)*100,
    donor, 
    deletion,
    celltype = anno_vec[colnames(mat)]
  )
}

pct_in_deletion_df <- rbind(
  make_pct_umis(bci, genes_in_pt1, "PT1", "PT1 del", anno_bci),
  #make_pct_umis(ccf, genes_in_pt1, "PT2", "delPT1", anno_ccf),
  #make_pct_umis(pt3, genes_in_pt1, "PT3", "delPT1", anno_mds),
  make_pct_umis(h1, genes_in_pt1, "H1", "PT1 del", anno_h1),
  
  #make_pct_umis(bci, genes_in_pt2, "PT1", "delPT2", anno_bci),
  make_pct_umis(ccf, genes_in_pt2, "PT2", "PT2 del", anno_ccf),
  #make_pct_umis(pt3, genes_in_pt2, "PT3", "delPT2", anno_mds),
  make_pct_umis(h1, genes_in_pt2, "H1", "PT2 del", anno_h1),
  
  #make_pct_umis(bci, genes_in_pt3, "PT1", "delPT3", anno_bci),
  #make_pct_umis(ccf, genes_in_pt3, "PT2", "delPT3", anno_ccf),
  make_pct_umis(pt3, genes_in_pt3, "PT3", "PT3 del", anno_mds),
  make_pct_umis(h1, genes_in_pt3, "H1", "PT3 del", anno_h1)
)
pct_in_deletion_df <- pct_in_deletion_df[complete.cases(pct_in_deletion_df),]

#make pvalues come out
pct_in_deletion_df %>%
  filter(deletion == "delPT2") %>%
  group_by(celltype) %>%
  do(w = wilcox.test(pct_umis_not_in_del~donor, data=., paired=FALSE)) %>% 
  summarise(celltype, Wilcox = w$p.value)

table(pct_in_deletion_df$donor)

p2 <- ggplot(pct_in_deletion_df, aes(x = celltype, y = pct_umis_not_in_del, color = donor)) +
  geom_boxplot(outlier.shape = NA) + facet_wrap(~deletion) + pretty_plot(fontsize = 8) +
  coord_cartesian(ylim = c(0,6)) + scale_color_manual(values  = c("darkgrey", "black", "dodgerblue4", "firebrick")) +
  labs(x = "Celltype", y = "% UMIs of MT genes outside deletion")  + theme(legend.position = "none")

p1 <- ggplot(pct_in_deletion_df, aes(x = celltype, y = pct_umis_in_del, color = donor)) +
  geom_boxplot(outlier.shape = NA) + facet_wrap(~deletion) + pretty_plot(fontsize = 8) +
  coord_cartesian(ylim = c(0,6)) + scale_color_manual(values  = c("darkgrey", "black", "dodgerblue4", "firebrick")) +
  labs(x = "Celltype", y = "% UMIs of MT genes inside deletion") + theme(legend.position = "none")

cowplot::ggsave2(cowplot::plot_grid(p1, p2, nrow = 1), 
                 file = "../plots/pct_UMIs_comparison.pdf", width = 7, height = 1.7)


