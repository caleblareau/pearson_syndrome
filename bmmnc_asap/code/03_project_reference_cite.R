library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(Matrix)

# Load the reference dataset
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
DefaultAssay(bm) <- 'ADT'

# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')
bm <- FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)
bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model = TRUE)
ref.obj <- bm; rm(bm)

#---------
# Import pearson asap data
#---------
import_deletion <- function(library_n){
  tab <- read.table(paste0("../data/mgatk_del/del_Pearson_BM_MNC_",library_n,".deletion_heteroplasmy.tsv"), header = TRUE)
  tab$cell_id <- gsub("-1", paste0("-", library_n),tab$cell_id )
  tab %>% filter(reads_all >= 10 & deletion == "del10381-15407" & version == "improved")
}

del_df <- rbind(import_deletion(1), import_deletion(2), import_deletion(3), import_deletion(4), import_deletion(5))
rownames(del_df) <- del_df$cell_id
counts <- Read10X_h5(filename = "../../../pearson_mtscatac_large_data_files/input/bonemarrow_mnc_asap/Pearson_ASAP_filtered_peak_bc_matrix.h5")

# Slower but convenient for the row parsing
metadata <- read.csv(
  file = "../data/Pearson_ASAP_singlecell.csv.gz",
  header = TRUE,
  row.names = 1
)

protein_mat <- readRDS("../data/Pearson_ASAP_allProteinCounts.rds")

# Do some quality control
# Not as stringent as it could be but reasonable

protein_QC <- data.frame(
  barcode = colnames(protein_mat),
  control_count = colSums(protein_mat[grepl("sotype", rownames(protein_mat)),]),
  total_count = colSums(protein_mat),
  pass_mito_qc = colnames(protein_mat) %in% del_df$cell_id
) %>% filter(pass_mito_qc)

# Set some QC thresholds
ggplot(protein_QC, aes(x = log10(total_count), y =control_count)) +
  geom_point()

ggplot(metadata %>% filter(cell_id != "None"), aes(x = log10(passed_filters), y = (peak_region_fragments)/passed_filters *100)) +
  geom_point()

protein_qc_barcodes <- protein_QC$barcode[protein_QC$total_count > 150 & protein_QC$control_count < 10]
mito_qc_barcodes <- protein_QC$barcode[protein_QC$total_count > 150 & protein_QC$control_count < 10]
atac_qc_barcodes <- metadata %>% filter(cell_id != "None") %>% filter(passed_filters > 10^3 & (peak_region_fragments/passed_filters) > 0.25) %>% rownames()

# Derive consensus barcodes based on QC of all 3 modalities 
consensus_barcodes <- intersect(protein_qc_barcodes, intersect(mito_qc_barcodes, atac_qc_barcodes))
length(consensus_barcodes) # over 20k which is nice  
counts_filt <- counts[,consensus_barcodes]
protein_mat_filt <- protein_mat[,consensus_barcodes]
del_df_filt <- del_df[consensus_barcodes,]
metadata_filt <- metadata[consensus_barcodes,]
full_meta_data <- cbind(metadata_filt,del_df_filt)

# Run Yuhan's transfer code
adt.matrix <- protein_mat_filt

rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD56(NCAM)")] <- "CD56"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD4-2")] <- "CD4"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD3-1")] <- "CD3"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD57_Recombinant" )] <- "CD57"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD38-1" )] <- "CD38"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "HLA-DR" )] <- "HLA.DR"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD278")] <-  "CD278-ICOS"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD11a/CD18")] <- "CD11a"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD127")] <- "CD127-IL7Ra"
adt.feature <- intersect(  rownames(ref.obj[["ADT"]]),  rownames(adt.matrix) )
setdiff(  rownames(ref.obj[["ADT"]]) ,  adt.feature)

# get supervised adt pca from reference
DefaultAssay(ref.obj) <- "ADT"
ref.obj <- RunSPCA(object = ref.obj, features = adt.feature, assay = "ADT", reduction.name = "sapca", reduction.key = "SAPC_", graph = "wsnn")
ref.obj <- ref.obj[adt.feature,]
ref.obj <- FindNeighbors(ref.obj)

# create another assay with matched ADTs
query <- CreateSeuratObject(counts =  adt.matrix[adt.feature,], assay = "ADT2")
DefaultAssay(query) <- "ADT2"
query <- NormalizeData(query, normalization.method = "CLR", margin = 2)
VariableFeatures(query ) <- adt.feature
query <- ScaleData(query) %>% RunPCA() %>% FindNeighbors()
anchors.transfer <- FindTransferAnchors(
  reference = ref.obj,
  query = query,
  reference.assay = "ADT",
  query.assay = 'ADT2', 
  reference.reduction = 'sapca',
  features = adt.feature,
  dims = 1:18,
  nn.method = "annoy",
  k.filter = NA,
  verbose = TRUE
)

### Find RPCA anchors, and replace the anchors in anchors.transfer
query.inte <- query
query.inte <- RunPCA(query.inte, features = adt.feature)
ref.obj.inte <- ref.obj
ref.obj.inte[["pca"]] <- ref.obj.inte[["sapca"]]
anchor.inte <- FindIntegrationAnchors(object.list = list(ref.obj.inte, query.inte),
                                      reference = c(1),
                                      anchor.features = adt.feature,
                                      dims = 1:18, 
                                      reduction = "rpca", 
                                      k.filter = NA)
# add integration anchors to transfer anchors
anchors.transfer@anchors <- anchor.inte@anchors[anchor.inte@anchors$dataset1 == 1 , 1:3]

query <- MapQuery(
  reference = ref.obj,
  query = query ,
  anchorset = anchors.transfer,
  refdata = list(celltype.l1 = "celltype.l1", 
                 celltype.l2 ="celltype.l2"),
  reference.reduction = "sapca",
  reduction.model = "wnn.umap"
)

table(query@meta.data$predicted.celltype.l1)
table(query@meta.data$predicted.celltype.l2)
query@meta.data$heteroplasmy <- del_df[rownames(query@meta.data),"heteroplasmy"]
DimPlot(query, reduction = "ref.umap", label = TRUE, group.by = "predicted.celltype.l2") 
library(BuenColors)
FeaturePlot(query, features = c("heteroplasmy"), reduction = "ref.umap") +
  scale_color_gradientn(colors=jdb_palette("solar_extra"))

FeaturePlot(query, features = c("predicted.celltype.l2.score"), reduction = "ref.umap") +
  scale_color_gradientn(colors=jdb_palette("solar_extra"))


query@meta.data %>% group_by(predicted.celltype.l2) %>%
  summarize(mean(heteroplasmy), median(heteroplasmy), count = n()) %>% data.frame()

            