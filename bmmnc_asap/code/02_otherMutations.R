library(SummarizedExperiment)

source("../../global_functions/variant_calling.R")

import_channel <- function(channeln){
  se <- readRDS(paste0("../../../pearson_mtscatac_large_data_files/input/bonemarrow_mnc_asap/mgatk_files/Pearson-ASAP-BMMNC-c",channeln,"_mgatk.rds"))
  colnames(se) <- gsub("-1", paste0("-", as.character(channeln)), colnames(se))
  se
  }
SE <- cbind(import_channel("1"), import_channel("2"), import_channel("3"), import_channel("4"), import_channel("5"))

SEfilt <- SE[,SE$depth > 20]
se_called <- call_mutations_mgatk(SEfilt)
rm(SE); rm(SEfilt)

data.frame(rowData(se_called)) %>% 
  dplyr::filter(strand_correlation > 0.65 & log10(vmr) > -2 & n_cells_conf_detected > 4) %>% pull(variant) -> possible
saveRDS(data.matrix(assays(se_called)[["allele_frequency"]][possible,]), file = "../output/interesting_AFs.rds")
#pbmc <- readRDS("../../../pearson_mtscatac_large_data_files/output/asap/pearson_asap_master_object.rds")

af_df <- (data.frame(t(af) , barcode = colnames(af)))
df <- data.frame(pbmc@reductions$umap@cell.embeddings, cluster = pbmc@meta.data$seurat_clusters, 
                 barcode = rownames(pbmc@meta.data),  heteroplasmy = pbmc@meta.data$heteroplasmy)
mdf <- merge(af_df, df)
sort(cor(data.matrix(mdf))["heteroplasmy",])

ggplot(mdf, aes(x =   X14476G.A , y = heteroplasmy)) + geom_point()
