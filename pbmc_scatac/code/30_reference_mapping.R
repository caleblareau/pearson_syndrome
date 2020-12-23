library(data.table)
library(Matrix)
library(Seurat)
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
library(SeuratDisk)
options(future.globals.maxSize = 4000 * 1024^2)
reference <- LoadH5Seurat("../../../pearson_mtscatac_large_data_files/input/pbmc_multimodal.h5seurat")

raw <- readRDS("../../../pearson_mtscatac_large_data_files/output/Perason_aggr_26April2020_geneActivityScores.rds")
raw <- raw[,grepl("-3",colnames(raw))]
raw <- raw[,colSums(raw) > 1000]

# Import
assign_azimuth <- function(raw){
  raw <- CreateSeuratObject(counts = raw, project = "RNA")
  raw <- SCTransform(raw)
  anchors <- FindTransferAnchors(
    reference = reference,
    query = raw,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )
  projected <- MapQuery(
    anchorset = anchors,
    query = raw,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
  )
  df <- data.frame(
    projected@meta.data,
    projected@reductions$ref.umap@cell.embeddings)
  df
}
# Import
assign_seuratv3 <- function(raw){
  raw <- CreateSeuratObject(counts = raw, project = "RNA")
  raw <- SCTransform(raw)
  anchors <- FindTransferAnchors(
    reference = reference,
    query = raw,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )
  projected <- MapQuery(
    anchorset = anchors,
    query = raw,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
  )
  df <- data.frame(
    projected@meta.data,
    projected@reductions$ref.umap@cell.embeddings)
  df
}

# Look at the typical seurat performance
rna <- CreateSeuratObject(
  counts = reference@assays$SCT@counts,
  meta.data = reference@meta.data
)
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures = 3000)
rna <- ScaleData(rna)
rna <- RunPCA(rna, npcs = 30)

pbmc.atac <- CreateSeuratObject(
  counts = raw,
  assay = "ACTIVITY"
)
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- FindVariableFeatures(pbmc.atac)
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac)

transfer.anchors <- FindTransferAnchors(reference = rna, query = pbmc.atac, features = VariableFeatures(object = rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$celltype.l2, weight.reduction = "cca")
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)

azimuth_atac <- assign_azimuth(raw)
azimuth_atac$cca_annotation <- celltype.predictions$predicted.id

df_BCI <- fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del6073-13095") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id)) %>%
  mutate(scale_heteroplasmy = scale(heteroplasmy))

df_CCF <- fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del8482-13447") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id)) %>%
  mutate(scale_heteroplasmy = scale(heteroplasmy))

df_PT3 <- fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id)) %>%
  mutate(scale_heteroplasmy = scale(heteroplasmy))


tdf <- merge(merge(azimuth_atac, df_PT3, by.x = "row.names", by.y = "barcode"), 
             fread('../../pt3_chr7_del_scatac/output/Pearson-PBMC.chr7DelQC.tsv'), by.x = "Row.names", by.y = "V4")
tdf$MDS <- (tdf$X7del < 0.25)
library(ggbeeswarm)
tdf %>% filter(reads_all > 15) %>%
  ggplot(aes(x = predicted.celltype.l1, y = heteroplasmy, color = MDS)) +
  geom_boxplot( width =0.5, fill = NA)

ggplot(tdf, aes(x = refUMAP_1, y = refUMAP_2, color = heteroplasmy)) + geom_point() +
  scale_color_gradientn(colors = jdb_palette("brewer_spectra")) +
  facet_wrap(~MDS)

ggplot(azimuth_atac, aes(x = refUMAP_1, y = refUMAP_2, color = predicted.celltype.l2)) +
  geom_point()

ggplot(azimuth_atac, aes(x = refUMAP_1, y = refUMAP_2, color = predicted.celltype.l2)) +
  geom_point()

pp <- azimuth_atac %>% group_by(cca_annotation, predicted.celltype.l2) %>%
  summarize(count =n ()) %>% ungroup() %>% group_by(predicted.celltype.l2) %>% mutate(prop = count / sum(count)*100)

ggplot(pp, aes(x =cca_annotation, y = predicted.celltype.l2, fill = prop , label = count)) + 
  geom_tile() + geom_text() + scale_fill_gradientn(colors = jdb_palette("solar_rojos")) + pretty_plot() + L_border() +
  labs(x = "CCA Annotation", y = "Azimuth projection annotation") + theme(legend.position = "bottom")
