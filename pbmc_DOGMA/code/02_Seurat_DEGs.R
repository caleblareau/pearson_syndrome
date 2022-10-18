library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)

# Import gene expression and set up object
m1 <- Read10X_h5("../data/PT1_rep1_mtDOGMA_gex.h5")
m2 <- Read10X_h5("../data/PT1_rep2_mtDOGMA_gex.h5")
colnames(m2) <- gsub("-1", "-2", colnames(m2))
full_mat <- cbind(m1, m2)
meta <- fread("../output/Pearson_DOGMA_full_meta_data.tsv") %>% data.frame()
rownames(meta) <- meta$cb
both_cells <- intersect(meta$cb, colnames(full_mat))

so <- CreateSeuratObject(
  full_mat[,both_cells], meta.data = meta[both_cells,]
) %>% NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData() %>%
  RunPCA() 


# Now do correlations
Idents(object = so) <- "predicted.celltype.l2"

# Get scores
mito <- c("MT-ATP6","MT-ATP8","MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ND1",
          "MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6")

tcells <- so[,grepl("CD8 TEM",so$predicted.celltype.l2)]%>% NormalizeData 
cor_df <- data.frame(cor = cor(data.matrix(t(tcells@assays$RNA@data)), tcells$clip_het, use = "pairwise.complete")) %>%
  arrange(desc(cor)) 
cor_df$gene <- rownames(cor_df)
cor_df <- cor_df[complete.cases(cor_df),] 
cor_df$rank <- 1:dim(cor_df)[1]
cor_df$color <- case_when(
  cor_df$gene %in% mito ~ "mito", 
  TRUE ~ "other"
)

library(BuenColors)
p1 <- ggplot(cor_df, aes(x = rank, y = cor, color = color)) +
  geom_point(size = 0.5) + scale_color_manual(values = c("firebrick", "lightgrey")) +
  pretty_plot(fontsize = 7) + L_border() + theme(legend.position = "none") +
  labs(x = "Association rank", y = "RNA cor. with heteroplasmy")

cowplot::ggsave2(p1, file = "../plots/hetcor_cd8tem.pdf", width =1.5, height = 1.5)

####

do_binomial_test <- function(df, ct1, ct2){
  ct1v <- df %>% filter(predicted.celltype.l2 %in% ct1) %>%
    mutate(x = heteroplasmy < 0.01) %>% pull(x)
  ct2v <- df %>% filter(predicted.celltype.l2 %in% ct2) %>%
    mutate(x = heteroplasmy < 0.01) %>% pull(x)
  
  data.frame(
    ct1, ct2, 
    pvalue = prop.test(x = c(sum(ct1v), sum(ct2v)), n = c(length(ct1v), length(ct2v)))$p.value
  )
}

process_donor <- function(df){
  rbind(
    do_binomial_test(df, "CD8 TEM", "CD8 Naive"),
    do_binomial_test(df, "CD4 TCM", "CD4 Naive"),
    do_binomial_test(df,  "MAIT", c("CD8 Naive", "CD4 Naive", "CD8 TEM", "CD4 TEM"))
  ) 
}
keep <- c("MAIT", "CD4 Naive", "CD8 Naive", "CD8 TEM", "CD4 TCM")
dfgo <- so@meta.data %>% filter(predicted.celltype.l2 %in% keep & coverage >= 5)
dfgo$heteroplasmy <- dfgo$cov_het
process_donor(dfgo)
table(dfgo$predicted.celltype.l2)

library(ggbeeswarm)
so@meta.data %>% filter(predicted.celltype.l2 %in% keep) %>%
  filter(reads_all >= 5) %>%
  ggplot(aes(color = predicted.celltype.l2,  x =clip_het)) +
  stat_ecdf()

px <- dfgo %>%
  ggplot(aes(color = predicted.celltype.l2, x = heteroplasmy)) +
  stat_ecdf() + scale_color_manual(values = c("dodgerblue", "dodgerblue4", "pink", "firebrick", "purple4")) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))
cowplot::ggsave2(px, file = paste0("../plots/pt1_dogma_covhet.pdf"), width = 2.3, height = 1.5)



####

# Import protein counts
import_kite_counts <- function(lib, nl){
  mtx <- fread(paste0("../data/protein_counts/",lib,"_counts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0("../data/protein_counts/",lib,"_counts.barcodes.txt.gz"), header = FALSE)[[1]], "-", nl)
  colnames(matx) <- (fread(paste0("../data/protein_counts/",lib,"_counts.genes.txt.gz"), header = FALSE)[[1]])
  return(t(matx))
}
l1a <- import_kite_counts("PT1_rep1_mtDOGMA_TSA_preAMP_S1", "1")
l1b <- import_kite_counts("PT1_rep1_mtDOGMA_TSA_sup_S3", "1")
common_l1s <- intersect(colnames(l1a), colnames(l1b))

l2a <- import_kite_counts("PT1_rep2_mtDOGMA_TSA_preAMP_S2", "2")
l2b <- import_kite_counts("PT1_rep2_mtDOGMA_TSA_sup_S4", "2")
common_l2s <- intersect(colnames(l2a), colnames(l2b))

new_set <- intersect(colnames(so), c(common_l1s,common_l2s))
so_protein <- so[,new_set]
cmat <- cbind(l1a[,common_l1s] + l1b[,common_l1s], l2a[,common_l2s] + l2b[,common_l2s])[,new_set]
so_protein[["ADT"]] <- CreateAssayObject(cmat)
DefaultAssay(pbmc_lll) <- "ADT"
so_protein <- so_protein  %>% 
  NormalizeData(assay = "ADT", normalization.method = "CLR", margin = 2) %>%
  ScaleData(assay = "ADT", do.scale = FALSE)

##
tcells <- so_protein[,grepl("CD8 Naive",so_protein$predicted.celltype.l2)]
cor_df <- data.frame(cor = cor(data.matrix(t(tcells@assays$ADT@data)), tcells$clip_het, use = "pairwise.complete")) %>%
  arrange(desc(cor)) 
cor_df$gene <- rownames(cor_df)
cor_df <- cor_df[complete.cases(cor_df),] 
cor_df$rank <- 1:dim(cor_df)[1]
cor_df

qplot(tcells$cov_het, tcells@assays$ADT@data["B3GAT1",])
