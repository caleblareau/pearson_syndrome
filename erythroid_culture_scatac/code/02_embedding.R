library(Signac)
library(dplyr)
library(harmony)
library(Seurat)
library(BuenColors)

# Import scATAC
counts <- Read10X_h5(filename = "../../../pearson_large_data_files/input/erythroid-culture/atac/Pearson_Invitro_D6D12atac_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../../../pearson_large_data_files/input/erythroid-culture/atac//Pearson_Invitro_D6D12atac_singlecell.csv",
  header = TRUE,
  row.names = 1
)

import_mgatk_assignments <- function(lib, new_n, day){
  
  # Import the donor assignments
  dt <- fread(paste0("../output/cell_assignments_per_channel/Pearson_",lib,"_assign.tsv")) %>% data.frame()
  dt$day <- day
  dt$library <- lib
  
  dels <- fread(paste0("../data/dels/del_Pearson_Healthy_",lib,".deletion_heteroplasmy.tsv")) %>%
    dplyr::filter(version == "improved" & deletion == "del10381-15407")
  
  # Now import the chromosome 7 
  cnv.dt <- fread(paste0("../../pt3_chr7_del_scatac/output/Pearson_Healthy_",lib,".chr7DelQC.tsv")) %>% data.frame()
  dt2 <- merge(merge(dt, cnv.dt, by.x = "barcode", by.y = "V4"), dels, by.x = "barcode", by.y = "cell_id")
  dt2$barcode <- gsub("-1", paste0("-", as.numeric(new_n)), dt2$barcode)
  dt2[complete.cases(dt2),]
}

assign_df <- rbind(
  import_mgatk_assignments("D6_1", "1", "D06"),
  import_mgatk_assignments("D6_2", "2", "D06"),
  import_mgatk_assignments("D12_1", "3", "D12"),
  import_mgatk_assignments("D12_2", "4", "D12")
)

metadata2 <- merge(assign_df, metadata, by.x = "barcode", by.y = "row.names")
rownames(metadata2) <- metadata2$barcode

# Verify
stopifnot(sum(metadata2$cell_id == "None") == 0)
stopifnot(dim(metadata2)[1] == dim(assign_df)[1])
ery <- metadata2 %>% dplyr::filter(assign == "Pearson" & reads_all > 20)

#good sanity check
ggplot(metadata2 %>% dplyr::filter(assign == "Pearson"),
       aes(x = day, y = heteroplasmy)) + geom_violin()
pearson <- metadata2 %>% dplyr::filter(assign == "Pearson" & day == "D06")
smoothScatter( x = pearson$heteroplasmy, y = pearson$X7del, nbin = 328, 
               colramp = colorRampPalette(c("white", jdb_palette("solar_rojos"))), 
               nrpoints = 0, #xlim = c(10,50),
               ret.selection = FALSE, xlab="% Heteroplasmy", ylab="CONICs mixture component")

# Filter for worthwhile cells
metadata3 <- metadata2 %>% dplyr::filter(assign %in% c("Pearson", "Control"))
metadata3$MDS <- ifelse(metadata3$assign == "Control", "Control", ifelse(
  metadata3$X7del < 0.15, "del7", "WT"))

so <- CreateSeuratObject(
  counts = counts[,metadata3$barcode],
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata3
)
dim(so)
rm(counts)
so$pct_reads_in_peaks <- so$peak_region_fragments / so$passed_filters * 100

ggplot(so@meta.data, aes(x = log10(passed_filters + 1), y = pct_reads_in_peaks)) +
  geom_point() + facet_grid(library~assign)

# Now filter cells
so_filt <- subset(
  x = so,
  subset = assign %in% c("Control", "Pearson") & 
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 25 &
    passed_filters > 10^3
)
dim(so_filt)
dim(so)
rm(so)

# Do dimension reduction
so_filt <- RunTFIDF(so_filt)
so_filt <- FindTopFeatures(so_filt, min.cutoff = 'q50')
so_filt <- RunSVD(
  object = so_filt,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

so_filt <- RunUMAP(object = so_filt, reduction = 'lsi', dims = 2:30)
so_filt <- FindNeighbors(object = so_filt, reduction = 'lsi', dims = 2:30)
so_filt <- FindClusters(object = so_filt, verbose = FALSE, algorithm = 3)
DimPlot(object = so_filt, label = TRUE) + NoLegend()
DimPlot(object = so_filt, group.by = "day") 
DimPlot(object = so_filt, group.by = "assign") 

# Get rid of donor/pearson effects if we can with harmony
so_filt_hm <- RunHarmony(
  object = so_filt,
  group.by.vars = c('assign'),
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)

DefaultAssay(so_filt_hm) <- "peaks"
so_filt_hm <- RunUMAP(so_filt_hm, dims = 2:30, reduction = 'harmony')
so_filt_hm <- FindNeighbors(object = so_filt_hm, reduction = 'harmony', dims = 2:30)
so_filt_hm <- FindClusters(object = so_filt_hm, verbose = FALSE, algorithm = 3, resolution = 0.4)
DimPlot(so_filt_hm, label = TRUE)

DimPlot(object = so_filt_hm, group.by = "day") 
DimPlot(object = so_filt_hm, group.by = "assign") 
DimPlot(object = so_filt_hm, group.by = "MDS") +
  scale_color_manual(values = c("grey", "red", "black"))

# Import gene activity scores
ga1 <- readRDS("../../../pearson_large_data_files/output/invitro_erythroid/D6_1_gene_activities.rds")
ga2 <- readRDS("../../../pearson_large_data_files/output/invitro_erythroid/D6_2_gene_activities.rds")
ga3 <- readRDS("../../../pearson_large_data_files/output/invitro_erythroid/D12_1_gene_activities.rds")
ga4 <- readRDS("../../../pearson_large_data_files/output/invitro_erythroid/D12_2_gene_activities.rds")
colnames(ga1) <- gsub("-1", "-1", colnames(ga1)); colnames(ga2) <- gsub("-1", "-2", colnames(ga2));
colnames(ga3) <- gsub("-1", "-3", colnames(ga3)); colnames(ga4) <- gsub("-1", "-4", colnames(ga4))

# Harmonize ga 
commmon_ga <- intersect(intersect(intersect(rownames(ga1), rownames(ga2)), rownames(ga3)), rownames(ga4))
ga <- cbind(ga1[commmon_ga,], ga2[commmon_ga,], ga3[commmon_ga,], ga4[commmon_ga,])[,colnames(so_filt_hm)]

so_filt_hm[["ACTIVITY"]] <- CreateAssayObject(counts = ga)
DefaultAssay(so_filt_hm) <- "ACTIVITY"
so_filt_hm <- FindVariableFeatures(so_filt_hm)
so_filt_hm <- NormalizeData(so_filt_hm)

FeaturePlot(so_filt_hm, c("CXCL14", "ALAS2", "GATA1", "LYZ"), max.cutoff = "q90")

FindMarkers(so_filt_hm, "11") %>% head(20)
