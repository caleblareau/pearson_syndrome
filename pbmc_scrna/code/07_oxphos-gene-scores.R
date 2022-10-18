library(Seurat)
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
library(SeuratDisk)
library(viridis)
options(future.globals.maxSize = 4000 * 1024^2)
reference <- LoadH5Seurat("../../../pearson_large_data_files/input/pbmc/pbmc_multimodal.h5seurat")

data.frame(
  ct = reference@meta.data$celltype.l2,
  t(reference@assays$SCT@data[c("STMN1", "TOP2A", "UBE2C", "MKI67", "PBK", "CENPF"),])
) %>% group_by(ct) %>%
  summarize(median(STMN1), median(TOP2A), median(UBE2C), median(MKI67), median(PBK), median(CENPF)) %>%data.frame()

# Get scores
oxphos_pathway <- c("SNORD138","MIR4691","COX17","MIR7113","ATP5PD","ATP5MG","UQCR11","COX4I1",
                    "COX5B","COX6A1","COX6A2","COX6B1","COX6C","COX7A1","COX7A2","COX7B","COX7C",
                    "COX8A","COX11","COX15","UQCRQ","DMAC2L","SLC25A4","SLC25A5","SLC25A6","UQCR10",
                    "NDUFS7","NDUFA1","NDUFA2","NDUFA3",
                    "NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFAB1","NDUFB1",
                    "NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
                    "NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFV1","NDUFS4","NDUFS5","NDUFS6",
                    "NDUFS8","NDUFV2","NDUFV3","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E",
                    "ATP5PB","ATP5MC1","ATP5MC2","ATP5MC3","ATP5ME","ATP5PF","ATP5PO","NDUFA12",
                    "SCO1","SDHA","SDHB","SDHC","SDHD","SURF1","UCP1","UCP2","UCP3","UQCRB","UQCRC1",
                    "UQCRC2","UQCRFS1","UQCRH","SLC25A14","COX7A2L","COX5A","ATP5IF1","SLC25A27","ATP5MF")
mito <- c("MT-ATP6","MT-ATP8","MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ND1",
          "MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6")
# Gene sets
heme <- c("ALAD","COX10","CPOX","EARS2","EPRS","FECH","HMBS","PPOX","QARS","UROD")
cholesterol <- c("FDFT1","FDPS","GGPS1","HMGCR","HMGCS1","IDI2","LSS","MVD","MVK","PMVK","SQLE")
glycine_serine <- c("PHGDH", "PSAT1", "PSPH", "SHMT2")
glycolysis <- c("ABCB6","ADORA2B","AGL","AGRN","AK3","AK4","AKR1A1","ALDH7A1","ALDH9A1","ALDOA",
                "ALDOB","ALG1","ANG","ANGPTL4","ANKZF1","ARPP19","ARTN","AURKA","B3GALT6","B3GAT1",
                "B3GAT3","B3GNT3","B4GALT1","B4GALT2","B4GALT4","B4GALT7","BIK","BPNT1","CACNA1H",
                "CAPN5","CASP6","CD44","CDK1","CENPA","CHPF","CHPF2","CHST1","CHST12","CHST2","CHST4",
                "CHST6","CITED2","CLDN3","CLDN9","CLN6","COG2","COL5A1","COPB2","CTH","CXCR4","CYB5A",
                "DCN","DDIT4","DEPDC1","DLD","DPYSL4","DSC2","ECD","EFNA3","EGFR","EGLN3","ELF3",
                "ENO1","ENO2","ERO1L","EXT1","EXT2","FAM162A","FBP2","FKBP4","FUT8","G6PD","GAL3ST1",
                "GALE","GALK1","GALK2","GAPDHS","GCLC","GFPT1","GLCE","GLRX","GMPPA","GMPPB","GNE",
                "GNPDA1","GOT1","GOT2","GPC1","GPC3","GPC4","GPR87","GUSB","GYS1","GYS2","HAX1",
                "HDLBP","HK2","HMMR","HOMER1","HS2ST1","HS6ST2","HSPA5","IDH1","IDUA","IER3",
                "IGFBP3","IL13RA1","IRS2","ISG20","KDELR3","KIF20A","KIF2A","LCT","LDHA","LDHC",
                "LHPP","LHX9","MDH1","MDH2","ME1","ME2","MED24","MERTK","MET","MIF","MIOX","MPI",
                "MXI1","NANP","NASP","NDST3","NDUFV3","NOL3","NSDHL","NT5E","P4HA1","P4HA2","PAM",
                "PAXIP1","PC","PDK3","PFKFB1","PFKP","PGAM1","PGAM2","PGK1","PGLS","PGM2","PHKA2",
                "PKM2","PKP2","PLOD1","PLOD2","PMM2","POLR3K","PPFIA4","PPIA","PPP2CB","PRPS1",
                "PSMC4","PYGB","PYGL","QSOX1","RARS","RBCK1","RPE","RRAGD","SAP30","SDC1","SDC2",
                "SDC3","SDHC","SLC16A3","SLC25A10","SLC25A13","SLC35A3","SLC37A4","SOD1","SOX9",
                "SPAG4","SRD5A3","STC1","STC2","STMN1","TALDO1","TFF3","TGFA","TGFBI","TKTL1","TPBG",
                "TPI1","TPST1","TSTA3","TXN","UGP2","VCAN","VEGFA","VLDLR","XYLT2","ZNF292")


reference <- AddModuleScore(reference, list(oxphos_pathway, heme, cholesterol, glycolysis, glycine_serine),
                            name = c("oxphos", "heme", "cholesterol", "glycolysis", "glycine_serine"))

pX <- FeaturePlot(reference, features = c("oxphos1"),raster = FALSE,
            min.cutoff = "q1", max.cutoff = "q95", reduction = "wnn.umap") + scale_color_viridis() +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(pX, file = "../plots/pbmc_umap_oxphos.png", width = 4, height = 4, dpi = 500)

library(forcats)
p1 <- ggplot(reference@meta.data %>%
         filter(grepl("CD4|CD8|MAIT", celltype.l2)), aes(x = fct_reorder(celltype.l2, oxphos1, median, .desc = TRUE), y = oxphos1)) +
  geom_violin() + geom_boxplot(fill = NA, color = "firebrick", outlier.shape = NA, width = 0.2) +
  pretty_plot(fontsize = 8) + L_border() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + labs(x = "T cell subset", y = "OXPHOS module score")

cowplot::ggsave2(p1, file = "../plots/tcell_oxphos_modules.pdf", width = 3, height = 2)



library(SeuratData)
#InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
bm <- AddModuleScore(bm, list(oxphos_pathway, heme, cholesterol, glycolysis, glycine_serine),
                     name = c("oxphos", "heme", "cholesterol", "glycolysis", "glycine_serine"))
bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "umap", 
              reduction.key = "umap_")
p1 <- DimPlot(bm, group.by = "celltype.l2",reduction = "umap", label = FALSE) +
   theme_void() + theme(legend.position = "none") + ggtitle("")
#cowplot::ggsave2(p1, file = "../plots/umap_clusters_BM.png", width = 4, height = 4, dpi = 500)

library(viridis)
FeaturePlot(bm, features = "glycolysis4", min.cutoff = "q1", max.cutoff = "q95") + scale_color_viridis()

p2 <- FeaturePlot(bm, features = c("glycolysis4"),
            min.cutoff = "q1", max.cutoff = "q95") + scale_color_viridis() +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p2, file = "../plots/umap_heme_glycolysis.png", width = 4, height = 4, dpi = 500)


