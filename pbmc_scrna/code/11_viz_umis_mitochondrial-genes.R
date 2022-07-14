library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)

bci <- Read10X("../data/pearson_bci/")
ccf <- Read10X("../data/pearson_ccf/")
pt3 <- Read10X("../data/pearson_mds/")

rdf <- rbind(
  data.frame(
    barcode = colnames(bci),
    coverage = bci["MT-ND5",], modality = "scRNA", patient = "PT1"
  ),
  data.frame(
    barcode = colnames(ccf),
    coverage = ccf["MT-ND5",], modality = "scRNA", patient = "PT2"
  ),
  data.frame(
    barcode = colnames(pt3),
    coverage = pt3["MT-CYB",], modality = "scRNA", patient = "PT3"
  )
)

# import scATAC 
adf <- rbind(
  data.frame(fread("../../pbmc_scatac/data/depth/Pearson-PBMC-BCI_v12-mtMask_mgatk.depthTable.txt"), modality = "mtscATAC", patient = "PT1"),
  data.frame(fread("../../pbmc_scatac/data/depth/Pearson-PBMC-CCF_v12-mtMask_mgatk.depthTable.txt"), modality = "mtscATAC", patient = "PT2"),
  data.frame(fread("../../pbmc_scatac/data/depth/Pearson-PBMC-PT3_v12-mtMask_mgatk.depthTable.txt"), modality = "mtscATAC", patient = "PT3")
) %>% filter(V2 > 1)
colnames(adf) <- c("barcode", "coverage", "modality", "patient")
tdf <- rbind(adf, rdf)
tdf %>% group_by(patient, modality) %>%
  filter(coverage > 1) %>%
  summarize(mean(coverage), median(coverage))

p1 <- ggplot(tdf %>% filter(coverage > 1), aes(x = patient, y = coverage, color = modality)) +
  geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 200)) +
  labs(x = "Patient", color = "Modality", y = "Molecular coverage") +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  pretty_plot(fontsize = 8) + L_border() 
cowplot::ggsave2(p1, file = "../plots/boxplot-coverage.pdf", width = 3, height = 1.5)

# Plot heteroplasmy estimates from the two modalities
plothetdf <- data.frame(
  donor = c("PT1", "PT1", "PT2", "PT2", "PT3", "PT3"),
  modality = c("scRNA", "mtscATAC", "scRNA", "mtscATAC", "scRNA", "mtscATAC"),
  heteroplasmy = c(22.48, 41.2, 1.02, 66.0, 13.38, 53.6)
)

p1 <- ggplot(plothetdf, aes(x = donor, y = heteroplasmy, fill = modality)) +
  geom_bar(stat = "identity", position="dodge", color = "black") + 
  labs(x = "Patient", color = "Modality", y = "% Heteroplasmy") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(expand = c(0,0))
cowplot::ggsave2(p1, file = "../plots/heteroplasmy-estimates-pb-compare.pdf", width = 2, height = 1.5)
