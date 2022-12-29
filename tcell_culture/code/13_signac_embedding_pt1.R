library(data.table)
library(BuenColors)
library(dplyr)
library(viridis)
library(Seurat)
library(Signac)

# Import and set up meta data
sc <- fread("../data/pt1_singlecell.csv.gz") %>%
  filter(is__cell_barcode == 1)
mdf <- merge(rbind(
  fread("../output/Pearson_HC_PT1_mix_CD3_CD28_Day14_rep1_assign.tsv"),
  fread("../output/Pearson_HC_PT1_mix_CD3_CD28_Day14_rep2_assign.tsv") %>%
    mutate(barcode = paste0(substr(barcode, 1, 16), "-2"))
), sc, by = "barcode") %>% mutate(frip = peak_region_fragments/passed_filters*100)
mdf <- mdf %>% filter(assign == "Pearson") %>% filter(frip > 30)
rownames(mdf) <- mdf$barcode

# Process mtscatac data
peaks <- Read10X_h5("../data/pt1_filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = peaks[,rownames(mdf)],
  sep = c(":", "-"),
  genome = NULL,
  fragments = '../../../pearson_large_data_files/input/tcell-culture/pt1_fragments.tsv.gz',
  min.cells = 1,
  min.features = 1
)
so <- CreateSeuratObject(chrom_assay, meta.data = mdf, assay = "peaks")

# LSI dim reduction
DefaultAssay(so) <- "peaks"
so <- RunTFIDF(so) %>% 
  FindTopFeatures( min.cutoff = 'q25') %>%
  RunSVD() %>% RunUMAP(dims = 2:25,reduction = "lsi") %>% FindNeighbors(dims = 2:25,reduction = "lsi")

DefaultAssay(so) <- "peaks"
so <- FindClusters(so, resolution = 0.5)
dim(so)
DimPlot(so, reduction = "umap")
FeaturePlot(so, features = "heteroplasmy") & scale_color_viridis()

# Now create the Seurat object
library(EnsDb.Hsapiens.v86)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "hg38"
seqlevelsStyle(annotations) <- 'UCSC'
annotations
annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))
annotations

# add the gene information to the object
DefaultAssay(so) <- "peaks"
Annotation(so) <- annotations

# create a gene activity by cell matrix
gene.activities <- GeneActivity(so)
so[['RNA']] <- CreateAssayObject(counts = gene.activities)
so <- NormalizeData(
  object = so,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(so$nCount_GA)
)

DefaultAssay(so) <- "RNA"
FeaturePlot(so, features = c("LEF1", "CD8B", "CD4", "CD8A","CXCL13", "CD86", "heteroplasmy"), reduction = "umap", max.cutoff = "q99",
            sort.cell = TRUE) &
  scale_color_viridis()

xdf <- data.frame(
  gene =rownames(so@assays$GA@data),
  value = cor(t(data.matrix(so@assays$GA@data)), so$heteroplasmy)
) %>% arrange(desc(value))
xdf <- xdf[complete.cases(xdf),]


so$caleb_cluster <- case_when(
  so$seurat_clusters %in% c(1) ~ "CD8 eff-like",
  so$seurat_clusters %in% c(0) ~ "Naive-like",
  so$seurat_clusters %in% c(2) ~ "CD4 eff-like" ,
)

p1 <- DimPlot(object = so, group.by = "caleb_cluster", pt.size = 2,
              reduction = "umap") +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1, file = "clabels-pt1.png", width = 3, height = 3, dpi = 400)

p2 <- ggplot(so@meta.data, aes(x = heteroplasmy, color = caleb_cluster)) +
  stat_ecdf() +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))
cowplot::ggsave2(p2, file = paste0("ecdf_tcell_clusters_PT1.pdf"), width = 2.35, height = 1.5)

p1 <- FeaturePlot(object = so, features = "CD86", sort.cell = TRUE, min.cutoff = "q10", max.cutoff = "q90",
              reduction = "umap") + scale_color_viridis() + 
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1, file = "PT1_GeneScores_CD86.png", width = 3, height = 3, dpi = 400)

p1 <- FeaturePlot(object = so, features = "heteroplasmy", sort.cell = TRUE,
                  reduction = "umap") + scale_color_viridis() + 
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(p1, file = "PT1_heteroplasmy.png", width = 3, height = 3, dpi = 400)


pn <- Nebulosa::plot_density(object = so, features = "CD86", 
                       reduction = "umap", size = 2) +
  theme_void() + theme(legend.position = "none") + ggtitle("")
cowplot::ggsave2(pn, file = "PT1_GeneScores_Nebulosa_CD86.png", width = 3, height = 3, dpi = 400)

                       