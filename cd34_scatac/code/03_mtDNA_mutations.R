library(data.table)
library(dplyr)
library(ggbeeswarm)
library(BuenColors)
library(viridis)
source("../../global_functions/variant_calling.R")
load("../output/CD34_umap_embedding_granja_proj3.rda")
SE <- readRDS("../../../pearson_large_data_files/input/CD34/Pearson-CD34-PT3_mgatk.rds")
SEfilt <- SE[,SE$depth > 20]
se_called <- call_mutations_mgatk(SEfilt)
rm(SE); rm(SEfilt)

# Yolo
called <- data.frame(rowData(se_called)) %>% 
  dplyr::filter(strand_correlation > 0.65 & log10(vmr) > -2 & n_cells_conf_detected >=3 ) %>% pull(variant)

interesting <- unique(c("5557T>A", "13970G>A", "1719G>A", "7836T>C", "14476G>A", "12242A>G", "9450G>A", "9235T>C",  called))
  
mutdf <- data.frame(barcode = colnames(se_called), 
           t(data.matrix(assays(se_called)[["allele_frequency"]][interesting,])))
mdf <- merge(projection_df_pearson_pt, mutdf, by = "barcode")

plot_mutation <- function(mutation, cutoff, colorplot){
  mdf2 <- mdf[,c("umap1", "umap2", mutation)]
  colnames(mdf2) <- c("UMAP_1", "UMAP_2", "mutation")
  mdf2$colorme <- mdf2[[3]] > cutoff
  p1 <- ggplot(mdf2 %>% arrange((colorme)),
               aes(x =  UMAP_1, y = UMAP_2, color = colorme)) + geom_point(size = 1.5) +
     theme_void() + scale_color_manual(values = c("lightgrey", colorplot)) +
    theme(legend.position = "none")
  cowplot::ggsave2(p1, file = paste0("../plots/muts/",mutation,".png"), width = 4, height = 4, dpi = 300)
}


# MDS enriched mutations
plot_mutation("X9450G.A", 0.15, "red")
plot_mutation("X9235T.C", 0.05, "red")

plot_mutation("X1719G.A", 0.15, "red")
plot_mutation("X7836T.C", 0.15, "red")

# Non-MDS enriched mutations
plot_mutation("X14476G.A", 0.5, "dodgerblue3")
plot_mutation("X12242A.G", 0.25, "dodgerblue3")

# RTE
plot_mutation("X13970G.A", 0.25, "purple4")
plot_mutation("X5557T.A", 0.35, "purple4")
