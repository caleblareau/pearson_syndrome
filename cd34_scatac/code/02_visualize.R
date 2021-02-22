library(data.table)
library(dplyr)
library(ggbeeswarm)
library(BuenColors)
library(viridis)
load("../output/CD34_umap_embedding_granja_proj3.rda")

mds <- fread("../../pt3_chr7_del_scatac/output/Pearson-CD34-PT3.chr7DelQC.tsv")
del <- fread("../data/Pearson_CD34_PT3.deletion_heteroplasmy.tsv")  %>% dplyr::filter(version == "improved")

qplot(mds$pct_in_del, mds$X7del)
projection_df_pearson_pt_merge <- left_join(left_join(projection_df_pearson_pt, mds, by = c("barcode"= "V4")), 
                                            del, by = c("barcode" = "cell_id"))
projection_df_pearson_pt_merge$MDS <- ifelse(is.na(projection_df_pearson_pt_merge$X7del), "no", as.character(projection_df_pearson_pt_merge$X7del< 0.5))

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

pMDS <- ggplot(projection_df_pearson_pt_merge[dim(projection_df_pearson_pt_merge)[1]:1,], aes(x= umap1, y = umap2, color = MDS)) +
  geom_point_rast(size = 1, raster.dpi = 500) +
  labs(x = "UMAP1", y= "UMAP2", color = "Donor") +
  pretty_plot(fontsize = 7) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "lightgrey", "red")) +
  theme(legend.position = "none")

pHet <- ggplot(projection_df_pearson_pt_merge[dim(projection_df_pearson_pt_merge)[1]:1,], aes(x= umap1, y = umap2, color = heteroplasmy)) +
  geom_point_rast(size = 1, raster.dpi = 500) +
  labs(x = "UMAP1", y= "UMAP2", color = "Donor") +
  pretty_plot(fontsize = 7) + L_border() + 
  theme(legend.position = "non") +
  scale_color_viridis(na.value="lightgrey", end = max(projection_df_pearson_pt_merge$heteroplasmy, na.rm = TRUE)/100)
cowplot:::ggsave2(cowplot::plot_grid(pMDS, pHet, nrow = 1),
                  filename = "../plots/Pearson_MDS_projection_CD34_twoPanel.pdf", width = 3.6, height = 1.8)



pPearsonmanual <- ggplot(projection_df_pearson_pt[dim(projection_df_pearson_pt)[1]:1,], aes(x= umap1, y = umap2, color = celltype, label = celltype)) +
  geom_point(size = 0.2) +
  theme_void() + theme(legend.position = "none") +
  scale_color_manual(values = c("Pearson" = "black", "base" = "lightgrey"))

pMDSmanual <- ggplot(projection_df_pearson_pt_merge[dim(projection_df_pearson_pt_merge)[1]:1,], aes(x= umap1, y = umap2, color = MDS)) +
  geom_point(size = 0.2) +
  theme_void() + theme(legend.position = "none") +
  scale_color_manual(values = c("dodgerblue3", "lightgrey", "red")) 

pHetmanual <- ggplot(projection_df_pearson_pt_merge[dim(projection_df_pearson_pt_merge)[1]:1,], aes(x= umap1, y = umap2, color = heteroplasmy)) +
  geom_point(size = 0.2) +
  theme_void() + theme(legend.position = "none") +
  scale_color_viridis(na.value="lightgrey", end = max(projection_df_pearson_pt_merge$heteroplasmy, na.rm = TRUE)/100)

cowplot:::ggsave2(pPearsonmanual,
                  filename = "../plots/Pearson_Projmanual.png", width = 3.6, height = 3.6, dpi = 300)

cowplot:::ggsave2(pMDSmanual,
                  filename = "../plots/Pearson_MDSmanual.png", width = 3.6, height = 3.6, dpi = 300)

cowplot:::ggsave2(pHetmanual,
                  filename = "../plots/Pearson_Hetmanual.png", width = 3.6, height = 3.6, dpi = 300)

ggplot(projection_df_pearson_pt[dim(projection_df_pearson_pt)[1]:1,], aes(x= umap1, y = umap2, color = celltype, label = celltype)) +
  geom_point(size = 0.2)


# Now assign the closet C1 cell
projection_df_C1_only <- projection_df_C1[projection_df_C1$celltype != "none",]
c1_umap_dist <- data.matrix(projection_df_C1_only[,c(2,3)]); rownames(c1_umap_dist) <- projection_df_C1_only[[1]]

# Assign closest cells-- pearson
projection_df_pearson_only <- projection_df_pearson_pt_merge[projection_df_pearson_pt_merge$celltype != "base",]
pearson_umap_dist <- data.matrix(projection_df_pearson_only[,c(2,3)])
d <- as.matrix(dist(rbind(c1_umap_dist,(pearson_umap_dist))))
pearson_assignments <- colnames(d)[max.col(-d[2075:dim(d)[1],1:2074])]

# Assign closest cells-- healthy
projection_df_control_only <- projection_df_control[projection_df_control$celltype != "base",]
control_umap_dist <- data.matrix(projection_df_control_only[,c(2,3)])
d <- as.matrix(dist(rbind(c1_umap_dist,(control_umap_dist))))
healthycontrol_assignments <- colnames(d)[max.col(-d[2075:dim(d)[1],1:2074])]


# MDS
data.frame(
  healthy = round(table(healthycontrol_assignments)/length(healthycontrol_assignments)*100,2),
  MDS = round(table(pearson_assignments[projection_df_pearson_pt_merge$MDS == "TRUE"])/length(pearson_assignments[projection_df_pearson_pt_merge$MDS == "TRUE"])*100,2),
  WT = round(table(pearson_assignments[projection_df_pearson_pt_merge$MDS == "FALSE"])/length(pearson_assignments[projection_df_pearson_pt_merge$MDS == "FALSE"])*100,2)
)[,c(1,2,4,6)] -> quant_df
colnames(quant_df) <- c("assign", "healthy", "MDS", "WT")
quant_dfp <- quant_df %>% reshape2::melt(id.vars = "assign")
pbar <- ggplot(quant_dfp, aes(x = variable, y = value, fill = assign)) + 
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = ejc_color_maps) +
  coord_flip() + pretty_plot(fontsize = 7) + L_border() + theme(legend.position = "none") + labs(x = "", y = "% of cells")

cowplot:::ggsave2(pbar,
                  filename = "../plots/Pearson_healthy_bars.pdf", width = 1.8, height = 1.4)

# Do a per-closest cell heteroplasmy
projection_df_pearson_only$assignment <- pearson_assignments

at_least_50 <- c("CLP", "LMPP", "CMP", "GMP", "MEP", "pDC")
pBoxhet <- ggplot(projection_df_pearson_only %>%
         dplyr::filter(assignment %in% at_least_50), aes(x = assignment, y = heteroplasmy, color = MDS)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Closet cell assignment", y = '% Heteroplasmy') +
  pretty_plot(fontsize = 7) + L_border() +
  scale_color_manual(values = c("dodgerblue3", "red")) 
cowplot::ggsave2(pBoxhet, file = "../plots/boxplot_CD34_heteroplasmy.pdf", width = 2.9, height = 1.8)

projection_df_pearson_only %>% dplyr::filter(assignment %in% at_least_50) %>%
  group_by(assignment, MDS) %>% summarize(median(heteroplasmy), count = n())

                                                                       