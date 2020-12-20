library(BuenColors)
library(dplyr)
library(data.table)

pbmc <- merge(fread("data/Pearson-PBMC.chr7DelQC.tsv"), 
              fread("data_in/pbmc_3donor_aggr_singlecell.csv"),
              by.x = "V4", by.y = "barcode")

cd34 <- merge(fread("data/Pearson-CD34-PT3.chr7DelQC.tsv"), 
              fread("data_in/Pearson-CD34-PT3_v12-mtMask_singlecell.csv.gz"),
              by.x = "V4", by.y = "barcode")

bmmnc <- merge(fread("data/Pearson-BMMNC.chr7DelQC.tsv"), 
              fread("data_in/Pearson-BMMNC-PT3_v12_singlecell.csv"),
              by.x = "V4", by.y = "barcode")

p1 <- ggplot(cd34, aes(x = pct_in_del)) +
  geom_histogram(fill = "lightgrey", color = "black", bins = 41) +
  aes(y=stat(count)/sum(stat(count))*100) + 
  pretty_plot(fontsize = 5) + L_border() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 11), breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(10, 50)) +
  geom_vline(xintercept = 26, color = "firebrick") +
  labs(x = "", y = "")

p2 <- ggplot(bmmnc, aes(x = pct_in_del)) +
  geom_histogram(fill = "lightgrey", color = "black", bins = 41) +
  aes(y=stat(count)/sum(stat(count))*100) + 
  pretty_plot(fontsize = 5) + L_border() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 11), breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(10, 50)) +
  geom_vline(xintercept = 26, color = "firebrick") +
  labs(x = "", y = "")


p3 <- ggplot(pbmc, aes(x = pct_in_del)) +
  geom_histogram(fill = "lightgrey", color = "black", bins = 41) +
  aes(y=stat(count)/sum(stat(count))*100) + 
  pretty_plot(fontsize = 5) + L_border() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 11), breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(10, 50)) +
  geom_vline(xintercept = 26, color = "firebrick") +
  labs(x = "", y = "")


cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, ncol = 1), 
                 file = "viz_histo.pdf", width = 1.8, height = 2)

sum(pbmc$pct_in_del < 26)/dim(pbmc)[1]
sum(cd34$pct_in_del < 26)/dim(cd34)[1]
sum(bmmnc$pct_in_del < 26)/dim(bmmnc)[1]
