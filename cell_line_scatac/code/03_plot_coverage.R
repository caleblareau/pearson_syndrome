library(SummarizedExperiment)
library(BuenColors)
library(zoo)

SE <- readRDS(paste0("../../../pearson_large_data_files/input/cellline_mix/Pearson-Mix-5_v12-mtMask_mgatk.rds"))
summary(SE$depth >= 20)
bdf <- merge(del_het, assign_df, by.x = "cell_id", by.y = "barcode")

resolution <- 5
pull_coverage <- function(cells){
  rollmean(rowMeans(assays(SE)[['coverage']][,cells]), resolution)
}

cov_df <- data.frame(
  pos = rollmean(1:16569,resolution),
  aPD1 = pull_coverage(bdf %>% dplyr::filter(fp_classify == "aPD1") %>% pull(cell_id) %>% unique()),
  aPD2 = pull_coverage(bdf %>% dplyr::filter(fp_classify == "aPD2") %>% pull(cell_id) %>% unique()),
  aPD3 = pull_coverage(bdf %>% dplyr::filter(fp_classify == "aPD3") %>% pull(cell_id) %>% unique()),
  C1 = pull_coverage(bdf %>% dplyr::filter(fp_classify == "C1") %>% pull(cell_id) %>% unique()),
  C2 = pull_coverage(bdf %>% dplyr::filter(fp_classify == "C2") %>% pull(cell_id) %>% unique())
)
color_vec <- c("aPD1" = "#2780FF", "aPD2" = "#960096", "aPD3" = "#000090",
               "C1" = "#333333", "C2" = "#999999", "unassigned" = "firebrick")
cov_df %>%  reshape2::melt(id.vars = "pos") %>% dplyr::filter(value > 5) %>%
  ggplot(aes(x = pos, y = value, color = variable)) +
  geom_line() + facet_grid(~variable) +
  scale_y_continuous(limits = c(1, 80)) +
  scale_color_manual(values = color_vec) +pretty_plot() +
  theme(legend.position = "none")  +
  labs(x = "Position on mtDNA chromosome", y = "Mean coverage / cell") -> pCov
cowplot::ggsave2(pCov, file = "../output/coverage_long_mix5.pdf", width = 7.2, height = 1.7)
