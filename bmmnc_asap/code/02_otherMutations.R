library(SummarizedExperiment)
library(dplyr)
library(BuenColors)
library(viridis)
if(FALSE){
  source("../../global_functions/variant_calling.R")
  
  import_channel <- function(channeln){
    se <- readRDS(paste0("../../../pearson_large_data_files/input/bmmnc/mgatk//Pearson-ASAP-BMMNC-c",channeln,"_mgatk.rds"))
    colnames(se) <- gsub("-1", paste0("-", as.character(channeln)), colnames(se))
    se
  }
  SE <- cbind(import_channel("1"), import_channel("2"), import_channel("3"), import_channel("4"), import_channel("5"))
  
  SEfilt <- SE[,SE$depth > 20]
  se_called <- call_mutations_mgatk(SEfilt)
  rm(SE); rm(SEfilt)
  
  data.frame(rowData(se_called)) %>% 
    dplyr::filter(strand_correlation > 0.65 & log10(vmr) > -2 & n_cells_conf_detected >=3 ) %>% pull(variant) -> possible
  saveRDS(data.matrix(assays(se_called)[["allele_frequency"]][possible,]), file = "../output/interesting_AFs.rds")
}

# Overlap other point mutations with embedding
af <- readRDS("../output/interesting_AFs.rds")
pearson_asap <- readRDS( "../../../pearson_large_data_files/output/asap/pearson_asap_master_object.rds")


af_df <- (data.frame(t(af) , barcode = colnames(af)))
df <- data.frame(pearson_asap@reductions$umap@cell.embeddings, cluster = pearson_asap@meta.data$seurat_clusters, 
                 barcode = rownames(pearson_asap@meta.data),  heteroplasmy = pearson_asap@meta.data$heteroplasmy,
                 X7del = pearson_asap@meta.data$X7del,pct_in_del = pearson_asap@meta.data$pct_in_del,
                 reads_all = pearson_asap@meta.data$reads_all,
                 chr7 = pearson_asap@meta.data$chr7)
mdf <- merge(af_df, df, by = "barcode")

# Be more stringent on filtering
mdf <- mdf %>% dplyr::filter(X7del> 0.6 | X7del < 0.1)
mdf <- mdf %>% dplyr::filter(X7del> 0.6 | X7del < 0.1)

ggplot(mdf %>% arrange((X9450G.A)), aes(y = X7del, x = pct_in_del, color = X9450G.A)) +
  geom_point(size = 1) + 
  scale_color_viridis() 

var_meta_df_raw <- data.frame(
  variant = rownames(af),
  mean = rowMeans(af)*100
)

# Collect more info about variants
lapply(1:dim(var_meta_df_raw)[1], function(i){
  var_df <- data.frame(
    cluster = mdf2$cluster,
    variant = mdf2[,i+1]
  )
  tcells <- var_df$cluster %in% c(0,5,3,16,6,2)
  rtes <- var_df$cluster %in% c(5,6)
  
  data.frame(
    tcell_mean = mean(var_df$variant[tcells])*100,
    nont_mean = mean(var_df$variant[!tcells])*100,
    rte_mean = mean(var_df$variant[ rtes])*100,
    nonrte_mean = mean(var_df$variant[!rtes &tcells])*100,
    
    mds_mean = mean(var_df$variant[mdf2$chr7 == "Monosomy7"])*100,
    wt_mean = mean(var_df$variant[mdf2$chr7 == "Wildtype"])*100, 
    pct_mds = mean(var_df$variant[mdf2$chr7 == "Monosomy7"] > 0.1) *100,
    pct_wt = mean(var_df$variant[mdf2$chr7 == "Wildtype"] > 0.1) *100
  )
}) %>% data.table::rbindlist() %>% data.frame() -> vars_df
var_meta_df<- cbind(var_meta_df_raw,vars_df)

p1 <- ggplot(var_meta_df, aes(x = pct_mds, y = pct_wt, label = variant)) + 
  geom_text(size = 2) +
  pretty_plot(fontsize = 8) + L_border()  +
  scale_x_continuous(limits = c(0,4)) +
  scale_y_continuous(limits = c(0,4)) + labs(x = "% heteroplasmy MDS cells", y = "% heteroplasmy WT cells")
cowplot::ggsave2(p1, file = "../plots/mds_wt_variants_scatter.pdf", width = 2.2, height = 2.2)

plot_mutation <- function(mutation, cutoff, colorplot){
  mdf2 <- mdf[,c("UMAP_1", "UMAP_2", mutation, "chr7")]
  colnames(mdf2) <- c("UMAP_1", "UMAP_2", "mutation", "chr7")
  mdf2$colorme <- mdf2[[3]] > cutoff
  p1 <- ggplot(mdf2 %>% arrange((colorme)),
               aes(x =  UMAP_1, y = UMAP_2, color = colorme)) + geom_point(size = 1.5) +
    facet_wrap(~chr7) + theme_void() + scale_color_manual(values = c("lightgrey", colorplot)) +
    theme(legend.position = "none")
  cowplot::ggsave2(p1, file = paste0("../plots/muts/",mutation,".png"), width = 8, height = 4, dpi = 300)
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
