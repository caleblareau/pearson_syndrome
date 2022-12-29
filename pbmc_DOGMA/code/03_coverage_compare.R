library(BuenColors)
library(dplyr)

dogma <- fread("../output/Pearson_DOGMA_full_meta_data.tsv")
mtscatac <- fread("../../pbmc_scatac/data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del6073-13095")
depth_t <- fread("../../pbmc_scatac/data/depth/Pearson-PBMC-BCI_v12-mtMask_mgatk.depthTable.txt")
vec <- depth_t[["V2"]]; names(vec) <- depth_t[["V1"]]
mtscatac$coverage <- vec[mtscatac$cell_id]

cdf <- data.frame(
  coverage = c(mtscatac$coverage, dogma$coverage),
  clip_coverage = c(mtscatac$reads_all, dogma$reads_all),
  what = c(rep("mtscatac", dim(mtscatac)[1]), rep("dogma", dim(dogma)[1]))
)

p1 <- ggplot(cdf %>%  dplyr::filter(coverage > 5), aes(x = coverage, color = what)) + 
  stat_ecdf() + scale_x_log10() +  pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  theme(legend.position = "none") + labs(x = "mtDNA coverage", y = "Cumulative fraction")

p2 <- ggplot(cdf %>%  dplyr::filter(clip_coverage > 5), aes(x = clip_coverage, color = what)) + 
  stat_ecdf() + scale_x_log10() + pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  theme(legend.position = "none") + labs(x = "# useful reads for clip heteroplasmy", y = "Cumulative fraction")
cowplot::ggsave2(cowplot::plot_grid(p1, p2), width = 3.8, height = 1.5, filename = "../plots/compare_tech.pdf")

cdf %>% filter(coverage > 5) %>%
  group_by(what) %>% summarize(median(coverage))

cdf %>% dplyr::filter(clip_coverage > 5) %>%
  group_by(what) %>% summarize(median(clip_coverage))
