library(Seurat)
library(dplyr)

# import data
healthymd <- readRDS("../../../pearson_mtscatac_large_data_files/output/Healthy_scATAC_labelTransfer_HQcells.rds")@meta.data
pearsonmd <- readRDS("../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")@meta.data
healthygs <- readRDS("../../../pearson_mtscatac_large_data_files/output/Healthy_26April2020_geneActivityScores.rds")
pearsongs <- readRDS("../../../pearson_mtscatac_large_data_files/output/Perason_aggr_26April2020_geneActivityScores.rds")

pearsonmd$predicted.id <- gsub("Activated B-cells", "Naive B-cells", as.character(pearsonmd$predicted.id))
tb <- table(pearsonmd$predicted.id)
celltypes <- names(tb)[tb >75 ]
genes_common <- intersect(rownames(pearsongs), rownames(healthygs))
pearsongs <- pearsongs[genes_common,]
healthygs <- healthygs[genes_common,]

table(pearsonmd$Patient)
lapply(celltypes, function(ct){
  
  cpmmer <- function(bc, mat){
    rs <- rowSums(mat[,colnames(mat) %in% bc])
    rs/sum(rs) *1000000
  }
  
  barcodes_h <- rownames(healthymd)[healthymd$predicted.id == ct] 
  barcodes_b <- pearsonmd %>% dplyr::filter(predicted.id == ct & Patient == "BCI") %>% pull(barcode)
  barcodes_c <- pearsonmd %>% dplyr::filter(predicted.id == ct & Patient == "CCF") %>% pull(barcode)
  barcodes_p <- pearsonmd %>% dplyr::filter(predicted.id == ct & Patient == "PT3") %>% pull(barcode)
  
  data.frame(
    celltype = ct,
    Healthy = cpmmer(barcodes_h,healthygs),
    BCI = cpmmer(barcodes_b,pearsongs),
    CCF = cpmmer(barcodes_c,pearsongs),
    PTI = cpmmer(barcodes_p,pearsongs), 
    gene = rownames(pearsongs)
  ) %>% reshape2::melt(id.vars = c("gene", "celltype"))
  
}) %>% rbindlist() %>% data.frame() -> all_cpm_melt

pB <- all_cpm_melt %>% dplyr::filter(gene == "JUN") %>%
  ggplot(aes(x = variable, y = celltype, fill = value)) + pretty_plot(fontsize = 7) + L_border()+
  geom_tile() + scale_y_discrete(limits = rev(levels(all_cpm_melt$celltype)), expand = c(0,0)) + 
  theme(legend.position = "bottom") + scale_fill_gradientn(colors = jdb_palette("brewer_violet")) +
  scale_x_discrete(expand = c(0,0)) 



cowplot::ggsave2(pB, file = "../output/JUN_GS.pdf", width = 1.7, height = 2.3)

