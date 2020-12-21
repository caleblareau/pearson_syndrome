library(BuenColors)


libs <- c("pearson_bci", "pearson_ccf", "pearson_mds", "healthy_pbmc_8k_v2-remap", "healthy_pbmc_4k_v2-remap", "healthy_pbmc_5k_nextgem", "healthy_pbmc_5k")

dir_base = "pearson_ccf"
df <- readRDS(paste0("../output/seurat_projected_meta_", dir_base, ".rds"))
ggplot(df, aes(x = refUMAP_1, y = refUMAP_2, color = predicted.celltype.l2)) +
  geom_point() 
