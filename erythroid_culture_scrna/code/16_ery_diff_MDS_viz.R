library(edgeR)
library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(stringr)
library(BuenColors)
load( "../output/Erythroid_Person_edgeR_MDSstraify-D6D12-DGE.rda")
d6merge <- merge(ery_df6_mds, ery_df6_nomds, by = "row.names")
d12merge <- merge(ery_df12_mds, ery_df12_nomds, by = "row.names")

# Verify correlations are robust
cor(d12merge$signedTstat.x, d12merge$signedTstat.y)
cor(d6merge$signedTstat.x, d6merge$signedTstat.y)

# Gene sets
heme <- c("ALAD","COX10","CPOX","EARS2","EPRS","FECH","HMBS","PPOX","QARS","UROD")
cholesterol <- c("FDFT1","FDPS","GGPS1","HMGCR","HMGCS1","IDI2","LSS","MVD","MVK","PMVK","SQLE")
serine_glycine <- c("PHGDH", "PSAT1", "PSPH", "SHMT2")

# Get scores
oxphos_pathway <- c("SNORD138","MIR4691","COX17","MIR7113","ATP5PD","ATP5MG","UQCR11","COX4I1",
                    "COX5B","COX6A1","COX6A2","COX6B1","COX6C","COX7A1","COX7A2","COX7B","COX7C",
                    "COX8A","COX11","COX15","UQCRQ","DMAC2L","SLC25A4","SLC25A5","SLC25A6","UQCR10",
                    "NDUFS7","MT-ATP6","MT-ATP8","MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ND1",
                    "MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT ND6","NDUFA1","NDUFA2","NDUFA3",
                    "NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFAB1","NDUFB1",
                    "NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
                    "NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFV1","NDUFS4","NDUFS5","NDUFS6",
                    "NDUFS8","NDUFV2","NDUFV3","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E",
                    "ATP5PB","ATP5MC1","ATP5MC2","ATP5MC3","ATP5ME","ATP5PF","ATP5PO","NDUFA12",
                    "SCO1","SDHA","SDHB","SDHC","SDHD","SURF1","UCP1","UCP2","UCP3","UQCRB","UQCRC1",
                    "UQCRC2","UQCRFS1","UQCRH","SLC25A14","COX7A2L","COX5A","ATP5IF1","SLC25A27","ATP5MF")

d12merge$gene_class <- case_when(
  d12merge$Row.names  %in% oxphos_pathway ~ "OXPHOS",
  d12merge$Row.names  %in% cholesterol ~ "cholesterol",
  d12merge$Row.names  %in% heme ~ "heme",
  d12merge$Row.names  %in% serine_glycine ~ "serine_glycine",
  TRUE ~ "zother"
)

d12merge$alpha <- ifelse(d12merge$gene_class == "zother", 0.1, 1)
p2 <- ggplot(d12merge %>% arrange(desc(gene_class)), aes(x = signedTstat.x, y = signedTstat.y, color = gene_class, alpha = alpha)) +
  geom_point(size = 0.4) +  labs(x = "Pearson-MDS / Healthy statistic", y = "Pearson-no MDS / Healthy statistic") +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "green4", "purple2", "black")) +
  pretty_plot(fontsize= 8) + L_border() + theme(legend.position = "none")
cowplot::ggsave2(p2, file = "../plots/scatter_alpha_yn_mds.pdf", width = 2.2, height = 2.2)
