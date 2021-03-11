library(Seurat)

process_counts_mean <- function(dir_base){
  df <- readRDS(paste0("../output/seurat_projected_meta_",dir_base,".rds"))
  df$barcode <- rownames(df)
  
  # Now import counts
  import_counts <- function(dir_base, short_id){
    data.dir <- paste0("../data/", dir_base)
    raw <- Read10X(data.dir = data.dir)
    colnames(raw) <- paste0( colnames(raw))
    raw
  }
  counts <- import_counts(dir_base, "")
  counts <- counts[grepl("^MT-", rownames(counts)),as.character(rownames(df))]
  data.frame(round(rowMeans(counts), 2))
}
count_df <- data.frame(
  H = process_counts_mean("healthy_pbmc_8k_v2-remap"),
  BCI = process_counts_mean("pearson_bci"),
  CCF = process_counts_mean("pearson_ccf"),
  PT3 = process_counts_mean("pearson_mds")
)
colnames(count_df) <- c("Healthy", "BCI", "CCF", "PT3")
count_df

coord <- c(3785, 4991, 10231, 10618, 11449, 13242, 14411, 15317, 6674, 7927, 9598, 8867, 8469)
names(coord) <- paste0("MT-", c("ND1", "ND2", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYB", "CO1", "CO2", "CO3", "ATP6", "ATP8"))
count_df_ordered <- reshape2::melt(data.matrix(count_df[names(sort(coord)),]))
count_df_ordered$Var2 <- factor(as.character(count_df_ordered$Var2), levels = c("PT3", "CCF", "BCI", "Healthy"))

# Make plot
p1 <- ggplot(count_df_ordered, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + scale_fill_gradientn(colors = jdb_palette("solar_blues")) + 
  pretty_plot(fontsize = 7) + L_border() + theme(legend.position = "bottom") +
  labs(x = "", y = "Donor", color = "Mean expression") +
  scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
p1
cowplot::ggsave2(p1, file = "../plots/mean_mtRNA_UMIs.pdf", width = 3, height = 1.85)
