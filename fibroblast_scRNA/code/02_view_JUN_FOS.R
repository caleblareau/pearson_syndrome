library(data.table)
library(dplyr)
library(BuenColors)

diff <- readRDS("../output/results_P1_C1.rds")
diff[,c("gene", "P1cpm", "P2cpm", "C1cpm")] %>% filter(gene %in% c("JUN", "FOS", "FOSB", "JUNB")) %>%
  reshape2::melt()-> mdf

pjunfos <- ggplot(mdf, aes(x = variable, y = gene, fill = value)) +
  geom_tile() + scale_fill_gradientn(colors = jdb_palette("solar_blues")) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete( expand = c(0,0)) +
  labs(x = "", y = "", fill = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "bottom") 

cowplot::ggsave2(pjunfos, file = "../output/viz4_fibroblasts.pdf", width = 1.3, height = 2.0)

  