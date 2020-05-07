library(Seurat)
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
options(future.globals.maxSize = 4000 * 1024^2)

# Import
import_scRNAseq <- function(dir_base, name, pheno){
  data.dir <- paste0("../data/", dir_base)
  raw <- Read10X(data.dir = data.dir); colnames(raw) <- paste0(name, "-", colnames(raw))
  
  # import scrublet results
  singlets <- fread(paste0("../data/scrublet_out/", dir_base, ".scrub.tsv")) %>%
    data.frame() %>% dplyr::filter(!called) %>% pull(barcode)
  
  # Filter for singlet genes adn non-mitos
  raw <- raw[!grepl("^MT", rownames(raw)),paste0(name, "-", substr(singlets,1,16))]
  raw <- CreateSeuratObject(counts = raw,  min.cells = 3, min.features = 200); raw$source <- name; raw$pheno <- pheno
  #raw <- SCTransform(raw)
  raw <- NormalizeData(raw)
  raw <- FindVariableFeatures(raw)
  raw
  
}

# Import
pbmc1 <- import_scRNAseq("healthy_pbmc_8k_v2-remap", "H1", "Healthy")
pbmc2 <- import_scRNAseq("healthy_pbmc_4k_v2-remap", "H1", "Healthy")
pbmc3 <- import_scRNAseq("healthy_pbmc_5k_nextgem", "H2", "Healthy")
pbmc4 <- import_scRNAseq("healthy_pbmc_5k", "H2", "Healthy")

pBCI <- import_scRNAseq("pearson_bci", "pBCI", "Pearson")
pCCF <- import_scRNAseq("pearson_ccf", "pCCF", "Pearson")
pPT3 <- import_scRNAseq("pearson_mds", "pPT3", "Pearson")

features <- SelectIntegrationFeatures(object.list = list(pBCI, pCCF, pPT3, pbmc1, pbmc2, pbmc3, pbmc4))
pbmc.list <- lapply(X = list(pBCI, pCCF, pPT3, pbmc1, pbmc2, pbmc3, pbmc4), FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = pbmc.list, reduction = "rpca", 
                                  dims = 1:30)
pbmc.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Do dimension reduction
pbmc.integrated <- ScaleData(pbmc.integrated, verbose = FALSE)

pbmc.integrated <- RunPCA(pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:30, verbose = FALSE)

pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:30, verbose = FALSE)
pbmc.integrated <- FindClusters(pbmc.integrated, verbose = FALSE, resolution = 0.5)

FeaturePlot(pbmc.integrated, features = "nCount_RNA")

if(FALSE){
  # make figs for paper
  pA <- DimPlot(pbmc.integrated, reduction = "umap", split.by = "source", ncol = 5) +
    theme_void() + ggtitle("") + theme(legend.position  = "") +
    scale_color_manual(values = jdb_palette("corona"))
  ggsave(pA, file = "../plots/donor_split.png", width = 14, height = 3, dpi = 500)
  
  ppdf <- DimPlot(pbmc.integrated, label = TRUE, label.size = 1.5) + 
    pretty_plot(fontsize = 5) + theme(legend.position = "none") + L_border() +
    scale_color_manual(values = jdb_palette("corona"))
  cowplot::ggsave2(ppdf, file = "../plots/colors.pdf", width = 2, height = 2)
  
  ppng <-  DimPlot(pbmc.integrated, label = FALSE) +
    theme_void() + theme(legend.position = 'none') +
    scale_color_manual(values = jdb_palette("corona"))
  cowplot::ggsave2(ppng, file = "../plots/colors.png", width = 4, height = 4, dpi = 500)
  
}

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pbmc.integrated <- CellCycleScoring(pbmc.integrated, s.features = s.genes, g2m.features = g2m.genes)

FeaturePlot(pbmc.integrated, features = "G2M.Score")


new.cluster.ids <- c(
  '0'='Naive CD4+ T-cells',
  '1'='CD14+ Monocytes',
  '2'='Memory CD4+ T-cells',
  '3'='Naive CD8+ T-cells',
  '4'='Effector CD8+ T-cells',
  '5'='Naive B-cells',
  '6'="NK-cells",
  '7'='Memory B-cells',
  '8'='CD16+ Monocytes',
  '9' = "low count 1",
  '10' = "Gamma-Delta T-cell",
  '11'='Dendritic cells',
  '12'='NK-like T-cells',
  '13'='CCF cell 1',
  '14'='CD56 bright NK-cell',
  '15'='pDC',
  '16'='Platlets',
  '17'='IFN-activated T-cells',
  '18' = 'low count 2',
  '19' = 'G2M cells',
  '20'='CCF cell 2',
  '21'='Innate-like B-cells',
  '22'="CLEC9A+ DCs"
)

fm <- FindMarkers(pbmc.integrated, ident.1 = "CCF cell 2", min.pct = 0.3)

# Make plots for a supplement
if(FALSE){
  counts_df <- pbmc.integrated@meta.data %>% group_by(source, celltype) %>% summarize(n = n()) %>% ungroup() %>% 
    group_by(source) %>% mutate(prop = round(n/sum(n)*100,1)) %>% ungroup() %>%
    group_by(celltype) %>% mutate(meanP = mean(prop), sd = sd(prop), Z = (prop - meanP)/sd)
  
  ct_order <- new.cluster.ids <- c(
    'Naive CD4+ T-cells',
    'Memory CD4+ T-cells',
    'Naive CD8+ T-cells',
    'Effector CD8+ T-cells',
    "Gamma-Delta T-cell",
    'IFN-activated T-cells',
    'NK-like T-cells',
    "NK-cells",
    'CD56 bright NK-cell',
    'Naive B-cells',
    'Memory B-cells',
    'Innate-like B-cells',
    'CD14+ Monocytes',
    'CD16+ Monocytes',
    'pDC',
    'Dendritic cells',
    "CLEC9A+ DCs",
    'Platlets',
    'CCF cell 1',
    'CCF cell 2',
    'G2M cells',
    "low count 1",
    'low count 2'
  )
  
  counts_df$source <- factor(as.character(counts_df$source), levels = c("pPT3", "pCCF", "pBCI", "H2", "H1"))
  counts_df$celltype <- factor(as.character(counts_df$celltype), levels = ct_order)
  
  p1 <- ggplot(counts_df, aes(x = celltype, y = source, label = prop, color = Z)) + 
    geom_text(size = 1.25) + geom_tile(fill = NA, color = "black") + pretty_plot(fontsize = 6) + L_border() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "", y = "") +
    scale_color_gradientn(colors= c("dodgerblue", "black", "red")) +  theme(legend.position = "bottom") +
    scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
  cowplot::ggsave2(p1, file = "../plots/abudances.pdf", width = 3.5, height = 2.0)
  
  p_dots <- pbmc.integrated@meta.data %>% group_by(celltype) %>% summarize(medUmi = median(nFeature_RNA), medCount = median(nCount_RNA), count = n()) %>%
    ggplot() + geom_point(aes(x = medUmi, y = medCount, color = celltype, size = log10(count)*0.2)) +
    scale_y_log10() + scale_x_log10() +
    pretty_plot(fontsize = 6) + theme(legend.position = "none") + L_border() +
    scale_color_manual(values = jdb_palette("corona")) + labs(x = "median genes detected / cell", y = "median UMIs / cell") +
    scale_size(range = c(0, 1))
  cowplot::ggsave2(p_dots, file = "../plots/p_dots.pdf", width = 1.7, height = 1.7)
}

fpcl <- function(feature){
  FeaturePlot(pbmc.integrated, features = c(feature), min.cutoff = "q8", pt.size = 0.05) + 
    theme_void() +theme(legend.position = "none") + ggtitle("")
}
ggsave(plot_grid(fpcl("CD14"), fpcl("FCGR3A"), fpcl("FCER1A"),  fpcl("IL3RA"), fpcl("CLEC9A"), ncol = 5), 
       file = "../plots/Monocyte.png",width = 14, height = 3, dpi = 1000)

ggsave(plot_grid(fpcl("PPBP"),  fpcl("MS4A1"), fpcl("CD27"), fpcl("FCER2"), fpcl("MZB1"), ncol = 5), 
       file = "../plots/rare_Bcell.png", width = 14, height = 3, dpi = 500)

ggsave(plot_grid(fpcl("CCL5"),  fpcl("CD3D"), fpcl("GZMK"), fpcl("NCAM1"), fpcl("CD8A"), ncol = 5), 
       file = "../plots/NKcell.png", width = 14, height = 3, dpi = 500)

ggsave(plot_grid(fpcl("CCR7"),  fpcl("CD4"), fpcl("S100A4"), fpcl("TRDC"), fpcl("MX1"), ncol = 5), 
       file = "../plots/Tcell.png", width = 14, height = 3, dpi = 500)

ggsave(plot_grid(fpcl("G2M.Score"), ncol = 1), 
       file = "../plots/G2M_score.png", width = 2.8, height = 3, dpi = 500)

Idents(pbmc.integrated) <- new.cluster.ids[as.character(pbmc.integrated@meta.data$integrated_snn_res.0.5)]
pbmc.integrated$celltype <- Idents(pbmc.integrated)
saveRDS(pbmc.integrated, "../../../pearson_mtscatac_large_data_files/output/5March-PearsonRNAseq-integration.rds")