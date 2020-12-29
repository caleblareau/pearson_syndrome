library(BuenColors)
library(dplyr)
library(data.table)

pbmc <- merge(fread("../output/Pearson-PBMC.chr7DelQC.tsv"), 
              fread("../../pbmc_scatac/data/"),
              by.x = "V4", by.y = "barcode")

cd34 <- merge(fread("../output/Pearson-CD34-PT3.chr7DelQC.tsv"), 
              fread("../../cd34_scatac/data/singlecell_sumstats/Pearson-CD34-PT3_v12-mtMask_singlecell.csv.gz"),
              by.x = "V4", by.y = "barcode")

bmmnc <- merge(fread("../output/Pearson-BMMNC.chr7DelQC.tsv"), 
              fread("../../bmmnc_scatac/data/Pearson-BMMNC-PT3_v12_singlecell.csv"),
              by.x = "V4", by.y = "barcode")

healthyC <- merge(fread("../output/Healthy-CD34_G10.chr7DelQC.tsv"), 
               fread("../../cd34_scatac/data/singlecell_sumstats/CD34_G10_v12-mtMask_singlecell.csv.gz"),
               by.x = "V4", by.y = "barcode")

healthyP <- merge(fread("../output/Healthy-PBMC_H10.chr7DelQC.tsv"), 
               fread("../../pbmc_scatac/data/control_singlecell/PBMC_H10_v12-mtMask_singlecell.csv"),
               by.x = "V4", by.y = "barcode")

p1 <- ggplot(cd34, aes(x = pct_in_del)) +
  geom_histogram(fill = "lightgrey", color = "black", bins = 41) +
  aes(y=stat(count)/sum(stat(count))*100) + 
  pretty_plot(fontsize = 5) + L_border() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 11), breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(10, 50)) +
  geom_vline(xintercept = 25, color = "firebrick") +
  labs(x = "", y = "")

p2 <- ggplot(bmmnc, aes(x = pct_in_del)) +
  geom_histogram(fill = "lightgrey", color = "black", bins = 41) +
  aes(y=stat(count)/sum(stat(count))*100) + 
  pretty_plot(fontsize = 5) + L_border() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 11), breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(10, 50)) +
  geom_vline(xintercept = 25, color = "firebrick") +
  labs(x = "", y = "")


p3 <- ggplot(pbmc, aes(x = pct_in_del)) +
  geom_histogram(fill = "lightgrey", color = "black", bins = 41) +
  aes(y=stat(count)/sum(stat(count))*100) + 
  pretty_plot(fontsize = 5) + L_border() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 11), breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(10, 50)) +
  geom_vline(xintercept = 25, color = "firebrick") +
  labs(x = "", y = "")

p4 <- ggplot(healthyC, aes(x = pct_in_del)) +
  geom_histogram(fill = "lightgrey", color = "black", bins = 41) +
  aes(y=stat(count)/sum(stat(count))*100) + 
  pretty_plot(fontsize = 5) + L_border() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 15), breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(10, 50)) +
  geom_vline(xintercept = 25, color = "firebrick") +
  labs(x = "", y = "")

p5 <- ggplot(healthyP, aes(x = pct_in_del)) +
  geom_histogram(fill = "lightgrey", color = "black", bins = 41) +
  aes(y=stat(count)/sum(stat(count))*100) + 
  pretty_plot(fontsize = 5) + L_border() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 11), breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(10, 50)) +
  geom_vline(xintercept = 25, color = "firebrick") +
  labs(x = "", y = "")

pbmc <- pbmc %>% filter(cell_id != "None")
cd34 <- cd34 %>% filter(cell_id != "None")
bmmnc <- bmmnc %>% filter(cell_id != "None")

healthyC <- healthyC %>% filter(cell_id != "None")
healthyP <- healthyP %>% filter(cell_id != "None")

cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, ncol = 1), 
                 file = "../output/visualize_histogram_pearson_maintext.pdf", width = 1.8, height = 2)

cowplot::ggsave2(cowplot::plot_grid(p4, p5, ncol = 1), 
                 file = "../output/visualize_histogram_healthy_supp.pdf", width = 1.8, height = 1.3)

sum(healthyC$pct_in_del < 25)/dim(healthyC)[1]
sum(healthyP$pct_in_del < 25)/dim(healthyP)[1]

sum(pbmc$pct_in_del < 25)/dim(pbmc)[1]
sum(cd34$pct_in_del < 25)/dim(cd34)[1]
sum(bmmnc$pct_in_del < 25)/dim(bmmnc)[1]
