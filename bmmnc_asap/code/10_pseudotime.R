library(data.table)
library(Seurat)
library(dplyr)
library(broom)
library(Matrix)
source("10a_PThelpers.R")
pearson_asap <- readRDS("../../../pearson_large_data_files/output/bone_marrow/pearson_asap_master_object.rds")
DimPlot(pearson_asap, label = TRUE)

# Semi supervised trajectory
ery_trajectory <- c('4', '13')

#Align single cells to Trajectory
df <- data.frame(row.names = colnames(pearson_asap),
                 x = pearson_asap@reductions$umap@cell.embeddings[,1],
                 y = pearson_asap@reductions$umap@cell.embeddings[,2],
                 Group = pearson_asap$seurat_clusters)
trajAligned1 <- alignTrajectory(df, ery_trajectory)

# simple function to add NAs to cells not in trajectory
augment <- function(df_all, df_t){
  need <- !(rownames(df_all) %in% rownames(df_t))
  dfa <- df_all[need,c("x", "y", "Group")]; dfa$ps <- NA
  colnames(dfa) <- c("x", "y", "Group", "pseudotime")
  total_df <- rbind(dfa, df_t)
  return(total_df)
}

df_ery <- augment(df, trajAligned1[[1]])
df_ery$pseudotime <- ifelse(df_ery$x > 0 | df_ery$y < -5, NA, df_ery$pseudotime)

# Visualize pseudotime
p2 <- ggplot(df_ery, aes(x,y,color=pseudotime)) + 
  geom_point(size = 0.5) + 
  viridis::scale_color_viridis(na.value = "lightgrey")  +
  pretty_plot(fontsize = 8) + theme_void() +
  theme(legend.position = "none")
cowplot::ggsave2(p2, file = "../plots/erythroid_invivo_pseudotime.png",
                 width = 5, height = 5, dpi = 500)

so_df <- data.frame(
  pearson_asap@meta.data[,c("heteroplasmy", "chr7", "reads_all")],
  (t(pearson_asap@assays$ADT@scale.data)[,c("CD71","CD72")])
)
mdf <- merge(df_ery, so_df, by = "row.names")
mdf$bins  <- cut(mdf$pseudotime, seq(0, 100, 10))

# visualize pseudotime
p1 <- ggplot(mdf[complete.cases(mdf),], aes(x = bins, y = heteroplasmy)) + 
  geom_boxplot(outlier.shape = NA) + pretty_plot(fontsize = 7) + 
  L_border() +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs( x = "Pseudotime bins", y = "Heteroplasmy (%)")
p1
cowplot::ggsave2(p1, file = "../plots/pseudotime_trajectory_heteroplasmy.pdf", width = 2.5, height = 1.6)

# Visualize mds pro
pbar <- mdf[complete.cases(mdf),] %>%
  group_by(bins) %>%
  summarize(mds_prop = mean(chr7 == "Monosomy7")) %>%
  ggplot(aes(x = bins, y = mds_prop*100, group = "one")) +
  geom_point() + geom_line() +
  #geom_bar(stat = "identity", color = "black", fill = "lightgrey") +
  pretty_plot(fontsize = 7) + L_border() + 
  labs( x = "Pseudotime bins", y = "% of cells with MDS deletion")
cowplot::ggsave2(pbar, file = "../plots/pseudotime_trajectory_mds.pdf", width = 3, height = 1.6)

mdf[complete.cases(mdf),] %>%
  group_by(bins) %>%
  summarize(mds_prop = mean(CD71)) %>%
  ggplot(aes(x = bins, y = mds_prop, group = "one")) +
  geom_point() + geom_line() +
  pretty_plot(fontsize = 7) + L_border() + 
  labs( x = "Pseudotime bins", y = "CD71expression")
