library(BuenColors)
library(data.table)
library(dplyr)
"%ni%" <- Negate("%in%")

df_BCI <- fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(reads_all>=20) %>% 
  dplyr::filter(deletion == "del6073-13095") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id)) %>%
  mutate(scale_heteroplasmy = scale(heteroplasmy))

df_CCF <- fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(reads_all>=20) %>% 
  dplyr::filter(deletion == "del8482-13447") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id)) %>%
  mutate(scale_heteroplasmy = scale(heteroplasmy))

df_PT3 <- fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id)) %>%
  dplyr::filter(reads_all>=20) %>% 
  mutate(scale_heteroplasmy = scale(heteroplasmy))

het_df <- rbind(df_BCI, df_CCF, df_PT3)

so <- readRDS("../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")
meta_data <- data.frame(
  so@meta.data,
  so@reductions$umap@cell.embeddings
)
bdf <- merge(meta_data, het_df, by = "barcode") %>% dplyr::filter(reads_all >= 20)

pi <- bdf$predicted.id
bdf$coarse <- case_when(
  grepl("B-cell", pi) ~ "Bcell",
  grepl("Monoc", pi) ~ "Monocyte",
  grepl("NK-ce", pi) ~ "NKcell",
  pi %in% c("pDC", "Dendritic cells") ~ "DC",
  grepl("T-cell", pi) ~ "Tcell",
  TRUE ~ "other"
)

library(ggbeeswarm)
bdf$predicted.id <- gsub("Activated B-cells", "Naive B-cells", bdf$predicted.id)
exclude_celltypes <- c("CCF cell 1", "IFN-activated T-cells", "Innate-like B-cells")
pbig <- ggplot(bdf %>% dplyr::filter(predicted.id %ni% exclude_celltypes), aes(x = predicted.id, y = heteroplasmy)) +
  geom_quasirandom(size = 0.2) + geom_boxplot(color = "dodgerblue2", outlier.shape = NA, fill = NA) + 
  facet_grid(Patient~coarse, scales = "free_x", space = "free_x") +
  pretty_plot(fontsize = 8) + labs (x = "", y = "Heteroplasmy") +
  stat_summary(fun=mean, geom="point", size=1, color="purple2", fill="purple2", shape = 15) 
cowplot::ggsave2(pbig, file = "../output/giant_grid.pdf", width = 7.8, height = 4)


p1 <- ggplot(bdf %>% dplyr::filter(coarse != "other" & Patient == "BCI"), aes(x = heteroplasmy, color = coarse)) +
  stat_ecdf() + pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") + scale_x_continuous(limits = c(0,100)) +
  labs(x = "Heteroplasmy", y = "Cumulative proportion", color = "celltype") +
  scale_color_manual(values =c ("#A65AC2", "#FFEC00", "#FF5A00", "#490C65", "#0081C9"))

p2 <- ggplot(bdf %>% dplyr::filter(coarse != "other" & Patient == "CCF"), aes(x = heteroplasmy, color = coarse)) +
  stat_ecdf() + pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") + scale_x_continuous(limits = c(0,100)) +
  labs(x = "Heteroplasmy", y = "Cumulative proportion", color = "celltype") +
  scale_color_manual(values =c ("#A65AC2", "#FFEC00", "#FF5A00", "#490C65", "#0081C9"))

p3 <- ggplot(bdf %>% dplyr::filter(coarse != "other" & Patient == "PT3"), aes(x = heteroplasmy, color = coarse)) +
  stat_ecdf() + pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") + scale_x_continuous(limits = c(0,100)) +
  labs(x = "Heteroplasmy", y = "Cumulative proportion", color = "celltype") +
  scale_color_manual(values =c ("#A65AC2", "#FFEC00", "#FF5A00", "#490C65", "#0081C9"))

cowplot::ggsave2(plot_grid(p1, p2, p3, nrow = 3), width = 1.1, height = 2.8, file = "../output/CDFs_3patients.pdf")

# Try different forms of smoothing
bdf$scale_plot_heteroplasmy <- ifelse(bdf$scale_heteroplasmy < -2, -2, bdf$scale_heteroplasmy)
bdf$scale_plot_heteroplasmy <- ifelse(bdf$scale_plot_heteroplasmy > 2, 2, bdf$scale_plot_heteroplasmy)

nn_smooth <- 20
nn_het <- matrix(bdf$scale_heteroplasmy[get.knn(bdf[,c("UMAP_1", "UMAP_2")], k = nn_smooth)$nn.index], ncol = nn_smooth)
bdf$smooth <- rowMeans(nn_het)
bdf$smooth_plot <-  ifelse(bdf$smooth < -2, -2, bdf$smooth)

p_het_scale <- ggplot(bdf, aes(x = UMAP_1, y = UMAP_2, color = smooth_plot)) + 
  geom_point(size = 0.5) + scale_color_gradientn(colors = jdb_palette("brewer_spectra")) +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(p_het_scale, file = "../output/scATAC_plot_scaled_heteroplasmy.png", width = 4, height = 4, dpi = 500)

exclude_celltypes <- c("CCF cell 1", "IFN-activated T-cells", "Innate-like B-cells")
ggplot(bdf %>% dplyr::filter(predicted.id %ni% exclude_celltypes), aes(x = predicted.id, y = heteroplasmy, color = Patient)) +
  geom_quasirandom() + geom_boxplot(color = "black", outlier.shape = NA, fill = NA) + facet_wrap(~Patient, nrow = 3) +
  stat_summary(fun.y=mean, geom="point", shape=2, size=2, color="black", fill="black") 

short_include <- c("NK-cells", "Gamma-Delta T-cell", "NK-like T-cells", "Effector CD8+ T-cells")
pbox <- ggplot(bdf %>% dplyr::filter(predicted.id %in% short_include), aes(x = predicted.id, y = heteroplasmy)) +
  geom_quasirandom(size = 0.2) + geom_boxplot(color = "dodgerblue2", outlier.shape = NA, fill = NA, width = 0.5) + 
  facet_wrap(~Patient, nrow = 3) + scale_x_discrete(limits = short_include) + pretty_plot(fontsize = 7) + 
  stat_summary(fun=mean, geom="point", size=1, color="green3", fill="green3", shape = 15)  +
  labs(x = "Annotated celltype", y = "% Heteroplasmy")
cowplot::ggsave2(pbox, file = "../output/boxplot_quasi_subset.pdf", height = 2.5, width = 3)

mwtest <- function(ctb, donor) {
  nk <-  bdf %>% dplyr::filter(predicted.id == "NK-cells" & Patient == donor) %>% pull(heteroplasmy)
  ct <-  bdf %>% dplyr::filter(predicted.id == ctb & Patient == donor) %>% pull(heteroplasmy)
  wilcox.test(ct,nk)$p.value
}

donor <- c("BCI", "CCF", "PT3")
cts <- c("Gamma-Delta T-cell", "NK-like T-cells", "Effector CD8+ T-cells")
lapply(donor, function(d){
  lapply(cts, function(cc){
   data.frame(donor = d, celtype = cc, pvalue =  mwtest(cc, d))
  }) %>% rbindlist() %>% data.frame() 
})%>% rbindlist() %>% data.frame()  -> all_pvalues_df
