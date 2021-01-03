library(data.table)
library(dplyr)
library(ggbeeswarm)
library(BuenColors)
library(ggrastr)
library(viridis)
load("../output/CD34_umap_embedding_granja_proj3.rda")

mds <- fread("../../pt3_chr7_del_scatac/output/Pearson-CD34-PT3.chr7DelQC.tsv")
del <- fread("../data/Pearson_CD34_PT3.deletion_heteroplasmy.tsv")  %>% dplyr::filter(version == "improved")

qplot(mds$pct_in_del, mds$X7del)
projection_df_pearson_pt_merge <- left_join(left_join(projection_df_pearson_pt, mds, by = c("barcode"= "V4")), 
                                            del, by = c("barcode" = "cell_id"))

ph1 <- projection_df_pearson_pt_merge %>%
  dplyr::filter(celltype == "Pearson" & reads_all >=15 & X7del <= 0.5) %>%
  ggplot(aes(x = heteroplasmy)) +
  geom_histogram(aes(y=..count../sum(..count..)), fill = "purple3") +
  scale_x_continuous(limits = c(-2,100)) +
  scale_y_continuous(limits = c(0,1)) +
  pretty_plot(fontsize = 7) + L_border() + 
  labs(x = "% Heteroplasmy", y = "Proportion")

ph2 <- projection_df_pearson_pt_merge %>%
  dplyr::filter(celltype == "Pearson" & reads_all >=15 & X7del > 0.5) %>%
  ggplot(aes(x = heteroplasmy)) +
  geom_histogram(aes(y=..count../sum(..count..)), fill = "black") +
  scale_x_continuous(limits = c(-2,100)) +
  scale_y_continuous(limits = c(0,1)) +
  pretty_plot(fontsize = 7) + L_border() + 
  labs(x = "% Heteroplasmy", y = "Proportion")

if(FALSE){
  cowplot:::ggsave2(cowplot::plot_grid(ph1, ph2,  nrow = 2),
                    filename = "../plots/histo2supplement.pdf", width = 1.8, height = 2.5)
}

pCDF <- projection_df_pearson_pt_merge %>%
  dplyr::filter(celltype == "Pearson" & reads_all >=15) %>%
  ggplot(aes(color = X7del > 0.5, x = heteroplasmy)) +
  stat_ecdf() +
  pretty_plot(fontsize = 7) + L_border() +
  scale_color_manual(values = c("purple3", "black")) +
  labs(x = "% Heteroplasmy", y = "Cumulative distribution") +
  theme(legend.position = "none")


if(FALSE){
  cowplot:::ggsave2(pCDF,
                    filename = "../plots/CDF_Heteroplasmy_Monosomy.pdf", width = 1.6, height = 1.6)
}

p1 <- ggplot(projection_df_C1[dim(projection_df_C1)[1]:1,], aes(x= umap1, y = umap2, color = celltype, label = celltype)) +
  geom_point_rast(size = 1, raster.dpi = 500) +
  labs(x = "UMAP1", y= "UMAP2", color = "C1 FACS ") +
  pretty_plot(fontsize = 6) + L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c(ejc_color_maps, "none" = "lightgrey"))

p2 <- ggplot(projection_df_pearson_pt[dim(projection_df_pearson_pt)[1]:1,], aes(x= umap1, y = umap2, color = celltype, label = celltype)) +
  geom_point_rast(size = 1, raster.dpi = 500) +
  labs(x = "UMAP1", y= "UMAP2", color = "Donor") +
  pretty_plot(fontsize = 6) + L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c("Pearson" = "firebrick", "base" = "lightgrey"))

p3 <- ggplot(projection_df_control[dim(projection_df_control)[1]:1,], aes(x= umap1, y = umap2, color = celltype, label = celltype)) +
  geom_point_rast(size = 1, raster.dpi = 500) +
  labs(x = "UMAP1", y= "UMAP2", color = "Donor") +
  pretty_plot(fontsize = 6 ) + L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c("Healthy" = "dodgerblue3", "base" = "lightgrey"))

if(FALSE){
  cowplot:::ggsave2(cowplot::plot_grid(p1, p2, p3,  nrow = 1),
                    filename = "../plots/Pearson_C1_initial_projection_CD34_viz.pdf", width = 5.4, height = 1.8)
}

projection_df_pearson_pt_merge$MDS <- ifelse(is.na(projection_df_pearson_pt_merge$X7del), "no", as.character(projection_df_pearson_pt_merge$X7del< 0.5))
pMDS <- ggplot(projection_df_pearson_pt_merge[dim(projection_df_pearson_pt_merge)[1]:1,], aes(x= umap1, y = umap2, color = MDS)) +
  geom_point_rast(size = 1, raster.dpi = 500) +
  labs(x = "UMAP1", y= "UMAP2", color = "Donor") +
  pretty_plot(fontsize = 7) + L_border() + 
  scale_color_manual(values = c("black", "lightgrey", "purple3")) +
  theme(legend.position = "none")

pHet <- ggplot(projection_df_pearson_pt_merge[dim(projection_df_pearson_pt_merge)[1]:1,], aes(x= umap1, y = umap2, color = heteroplasmy)) +
  geom_point_rast(size = 1, raster.dpi = 500) +
  labs(x = "UMAP1", y= "UMAP2", color = "Donor") +
  pretty_plot(fontsize = 7) + L_border() + 
  theme(legend.position = "non") +
  scale_color_viridis(na.value="lightgrey", end = max(projection_df_pearson_pt_merge$heteroplasmy, na.rm = TRUE)/100)
cowplot:::ggsave2(cowplot::plot_grid(pMDS, pHet, nrow = 1),
                  filename = "../plots/Pearson_MDS_projection_CD34_twoPanel.pdf", width = 3.6, height = 1.8)


