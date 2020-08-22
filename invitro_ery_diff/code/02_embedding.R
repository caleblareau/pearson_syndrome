library(Signac)
library(dplyr)
library(harmony)

# Import data
counts <- Read10X_h5(filename = "../../../pearson_mtscatac_large_data_files/input/invitro_ery_diff/aggr/pearson_invitro_aggr_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../../../pearson_mtscatac_large_data_files/input/invitro_ery_diff/aggr/pearson_invitro_aggr_singlecell.csv",
  header = TRUE,
  row.names = 1
)

import_mgatk_assignments <- function(lib, new_n, day){
  dt <- fread(paste0("../output/cell_assignments_per_channel/Pearson_",lib,"_assign.tsv")) %>% data.frame()
  dt$day <- day
  dt$library <- lib
  dt$barcode <- gsub("-1", paste0("-", as.numeric(new_n)), dt$barcode)
  dt
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

so <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata2
)
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

so_filt_hm <- RunUMAP(so_filt_hm, dims = 2:30, reduction = 'harmony')
DimPlot(object = so_filt_hm, group.by = "day") 
DimPlot(object = so_filt_hm, group.by = "assign") 

import_del <- function(lib, new_n){
  dels <- fread(paste0("../data/dels/del_Pearson_Healthy_",lib,".deletion_heteroplasmy.tsv")) %>%
    dplyr::filter(version == "improved" & deletion == "del10381-15407")
  dels$barcode <- gsub("-1", paste0("-", as.numeric(new_n)), dels$cell_id)
  dels
}

del_df <- rbind(
  import_del("D6_1", "1"),
  import_del("D6_2", "2"),
  import_del("D12_1", "3"),
  import_del("D12_2", "4")
)

del_vec <- del_df$heteroplasmy; names(del_vec) <- del_df$barcode
so_filt_hm@meta.data$heteroplasmy <- ifelse(so_filt_hm$assign == "Control", NA, del_vec[colnames(so_filt_hm)])
FeaturePlot(so_filt_hm, feature = "heteroplasmy") +
  scale_color_gradientn(colors = jdb_palette("brewer_spectra"), na.value = "white")
