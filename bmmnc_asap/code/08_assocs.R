library(data.table)
library(Seurat)
library(dplyr)
library(broom)
library(Matrix)
mtDNA_variants <- readRDS("../output/interesting_AFs.rds")
pearson_asap <- readRDS( "../../../pearson_large_data_files/output/asap/pearson_asap_master_object.rds")

keep_cells <- intersect(colnames(mtDNA_variants), colnames(pearson_asap))
not_mds <- (pearson_asap@meta.data[keep_cells, "chr7"] == "Wildtype")
seurat_clusters <- as.factor(as.character(pearson_asap@meta.data[keep_cells, "seurat_clusters"]))
mtDNA_het <- mtDNA_variants[, keep_cells]
adts <- pearson_asap@assays$ADT@data[,keep_cells]
control_adts <- t(adts[grep("sotypeCtrl", rownames(adts)),])
colnames(control_adts) <- c("ct1", "ct2", "ct3", "ct4")
adts <-adts[!grepl("sotypeCtrl", rownames(adts)),]

# Perform all associations
lapply(1:dim(mtDNA_het)[1], function(i){
  var1 <- mtDNA_het[i,]
  df <- data.frame(
    var1, not_mds,control_adts,
    t(adts)
  ) %>% reshape2::melt(id.vars = c("var1", "not_mds", "ct1", "ct2", "ct3", "ct4"))
  ssdf <- df %>% group_by(variable) %>% do(tidy(lm(value ~ var1 + not_mds + ct1 + ct2 + ct3 + ct4 + seurat_clusters, data=.))) %>% dplyr::filter(term == "var1")
  ssdf$mut <- rownames(mtDNA_variants)[i]
  ssdf
}) %>% rbindlist() -> odf
odf$padj <- p.adjust(odf$p.value)
odf %>% arrange(padj) %>% head(20)

summary(mtDNA_het["9235T>C", not_mds] > 0.06)
summary(mtDNA_het["9235T>C", !not_mds] >0.06)

summary(mtDNA_het["9307G>A", not_mds] > 0.06)
summary(mtDNA_het["9307G>A", !not_mds] >0.06)
