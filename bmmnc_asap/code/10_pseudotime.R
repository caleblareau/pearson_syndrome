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

#cowplot::ggsave2(p2, file = "../plots/erythroid_invivo_pseudotime.png",
#                 width = 5, height = 5, dpi = 500)

so_df <- data.frame(
  pearson_asap@meta.data[,c("heteroplasmy", "chr7", "reads_all")],
  (t(pearson_asap@assays$ADT@scale.data)[,c("CD71","CD235ab")])
)
mdf <- merge(df_ery, so_df, by = "row.names")
mdf$bins  <- cut(mdf$pseudotime, seq(0, 100, 10))

# visualize pseudotime
library(ggbeeswarm)
library(viridis)
p1 <- ggplot(mdf[complete.cases(mdf),], aes(x = bins, y = heteroplasmy, color = pseudotime)) + 
  geom_quasirandom(size = 0.2) + pretty_plot(fontsize = 7) + scale_color_viridis() + 
  stat_summary(
    geom = "point",
    fun = "median",
    col = "black",
    size = 1.5,
    shape = 15,
    fill = "black"
  ) +
  L_border() +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs( x = "Pseudotime bins", y = "Heteroplasmy (%)") + theme(legend.position = "none")
p1
cowplot::ggsave2(p1, file = "../plots/pseudotime_trajectory_heteroplasmy.pdf", width = 1.8, height = 1.6)

# Visualize mds pro
pbar <- mdf[complete.cases(mdf),] %>%
  group_by(bins) %>%
  summarize(mds_prop = mean(chr7 == "Monosomy7")*100, mds_sem = sqrt(var(chr7 == "Monosomy7"))/sqrt(n())*100) %>%
  ggplot(aes(x = bins, y = mds_prop, group = "one")) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mds_prop-mds_sem, ymax=mds_prop+mds_sem), width=.2) +
  pretty_plot(fontsize = 7) + L_border() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs( x = "Pseudotime bins", y = "% of cells with MDS deletion")
pbar

cowplot::ggsave2(pbar, file = "../plots/pseudotime_trajectory_mds.pdf", width = 1.8, height = 1.6)

pmito_coverage <- mdf[complete.cases(mdf),] %>%
  group_by(bins) %>%
  summarize(mean_cov = mean(reads_all), cov_sem = sqrt(var(reads_all))/sqrt(n())) %>%
  ggplot(aes(x = bins, y = mean_cov, group = "one")) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean_cov-cov_sem, ymax=mean_cov+cov_sem), width=.2) +
  pretty_plot(fontsize = 7) + L_border() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs( x = "Pseudotime bins", y = "mtDNA coverage")
pmito_coverage

cowplot::ggsave2(pmito_coverage, file = "../plots/pseudotime_trajectory_mito_coverage.pdf", width = 2.7, height = 2)

mdf[complete.cases(mdf),] %>%
  group_by(bins) %>%
  summarize(mds_prop = mean(CD235ab)) %>%
  ggplot(aes(x = bins, y = mds_prop, group = "one")) +
  geom_point() + geom_line() +
  pretty_plot(fontsize = 7) + L_border() + 
  labs( x = "Pseudotime bins", y = "CD71expression")



#####################################

# Order the cells and the gene / tags
ordered_cells <- mdf[!is.na(mdf$pseudotime),] %>% arrange((pseudotime))
adt_mat_ordered <- data.matrix(pearson_asap@assays$ADT@data)[,ordered_cells$Row.names]
markers_order <- c("CD71", "CD235a")
(markers_order)[!(markers_order %in% rownames(adt_mat_ordered))]
adt_mat_ordered <- adt_mat_ordered[markers_order,]

# Now group/order/smooth the cell states and make ArchR's nice heatmaps
# see: https://github.com/GreenleafLab/ArchR/blob/5f855d7b7ff3f57eb0d28f312c5ea5d373d27ebd/R/Trajectory.R
n_bins <- 100
groups_adt <- sapply(1:n_bins, function(idx){
  multiplier <- 100/n_bins
  bcs <- ordered_cells %>% dplyr::filter(pseudotime >= (idx-1)*multiplier & pseudotime < idx*multiplier) %>% pull(Row.names)
  rowMeans(adt_mat_ordered[,bcs], na.rm = TRUE)
})

groups_adt <- sapply(1:n_bins, function(idx){
  multiplier <- 100/n_bins
  bcs <- ordered_cells %>% dplyr::filter(pseudotime >= (idx-1)*multiplier & pseudotime < (idx*multiplier)) %>% pull(Row.names)
  rs <- rowSums(data.matrix(pearson_asap@assays$ADT@counts)[markers_order,bcs, drop = FALSE], na.rm = TRUE)
  log2(rs/sum(rs) *10000 + 1) # basic normalization
})

# Smooth and preprocess
smoothWindow = 16
smooth_groups_adt <- data.matrix((apply((groups_adt), 1, function(x) ArchR:::.centerRollMean((x), k = smoothWindow))))
smooth_groups_minmax_adt <- t(apply(smooth_groups_adt, 2, function(x)(x-min(x))/(max(x)-min(x))))

# Now make the plot
library(ComplexHeatmap)
pdf("../plots/erythroid_heatmap_ADT.pdf", width = 4, height = 1.5)
Heatmap(cbind(smooth_groups_minmax_adt[markers_order,]),
        col=as.character(jdb_palette("solar_rojos",type="continuous")),
        show_row_names = TRUE, 
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 0),
        show_column_names = FALSE)
dev.off()

# Now do GS
gs_mat_ordered <- data.matrix(pearson_asap@assays$ACTIVITY@counts)[,ordered_cells$Row.names]
groups_gs <- sapply(1:n_bins, function(idx){
  multiplier <- 100/n_bins
  bcs <- ordered_cells %>% dplyr::filter(pseudotime >= (idx-1)*multiplier & pseudotime < (idx*multiplier)) %>% pull(Row.names)
  print(length(bcs))
  rs <- rowSums(gs_mat_ordered[,bcs, drop = FALSE], na.rm = TRUE)
  log2(rs/sum(rs) *10000 + 1) # basic normalization
})
smooth_groups_gs <- data.matrix((apply((groups_gs), 1, function(x) ArchR:::.centerRollMean((x), k = smoothWindow))))
smooth_groups_minmax_gs <- t(apply(smooth_groups_gs, 2, function(x)(x-min(x))/(max(x)-min(x))))

# Make plot
genes_order <- c("HBB", "GATA1", "TMCC2", "HBA1", "HBA2", "ALAS2")
genes_order %in% rownames(groups_gs)
pdf("../plots/erythroid_heatmap_chromatin.pdf", width = 4, height = 1.5)
Heatmap(cbind(smooth_groups_minmax_gs[genes_order,]),
        col=as.character(jdb_palette("solar_blues",type="continuous")),
        show_row_names = TRUE, 
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 0),
        show_column_names = FALSE)
dev.off()

data.frame(cor  = cor(t(gs_mat_ordered), ordered_cells$pseudotime)) %>%
  arrange(desc(cor)) %>% head(50)
