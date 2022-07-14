library(BuenColors)
library(dplyr)
library(data.table)

dir_base <- "pearson_mds"
rdf <- readRDS(paste0("../output/seurat_projected_meta_", dir_base, ".rds"))
mds <- fread("../output/PBMC_scRNA_MDS_annotations.tsv") %>%
  filter(what == "Pearson" ) %>% mutate(barcodeSimple = paste0(substr(barcode,1,16), "-1"))
mdf<- merge(rdf, mds, by.x = "row.names", by.y = "barcodeSimple")         

p1 <- ggplot(mdf %>% arrange((MDS)), aes(x = refUMAP_1, y = refUMAP_2, color = MDS)) +
  geom_point(size = 0.2) + scale_color_manual(values = c("lightgrey", "red")) +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/mds_pbmc_embedding.png", width = 3.4, height = 3.4, dpi = 500)
