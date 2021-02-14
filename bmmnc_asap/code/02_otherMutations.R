library(SummarizedExperiment)
library(dplyr)
library(BuenColors)
library(viridis)
if(FALSE){
  source("../../global_functions/variant_calling.R")
  
  import_channel <- function(channeln){
    se <- readRDS(paste0("../../../pearson_mtscatac_data_files/input/bonemarrow_mnc_asap/mgatk_files/Pearson-ASAP-BMMNC-c",channeln,"_mgatk.rds"))
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
}

# Overlap other point mutaitons with embedding
af <- readRDS("../output/interesting_AFs.rds")
pearson_asap <- readRDS("../../../pearson_large_data_files/output/asap/pearson_asap_master_object.rds")

af_df <- (data.frame(t(af) , barcode = colnames(af)))
df <- data.frame(pearson_asap@reductions$umap@cell.embeddings, cluster = pearson_asap@meta.data$seurat_clusters, 
                 barcode = rownames(pearson_asap@meta.data),  heteroplasmy = pearson_asap@meta.data$heteroplasmy,
                 chr7 = pearson_asap@meta.data$chr7)
mdf <- merge(af_df, df)

var_meta_df_raw <- data.frame(
  variant = rownames(af),
  mean = rowMeans(af)*100
)

# Collect more info about variants
lapply(1:dim(var_meta_df_raw)[1], function(i){
  var_df <- data.frame(
    cluster = mdf$cluster,
    variant = mdf[,i+1]
  )
  tcells <- var_df$cluster %in% c(0,5,3,16,6,2)
  
  res.kw <- kruskal.test(variant ~ cluster, data = var_df)
  data.frame(kruskal_pvalue = (res.kw)[[3]], 
             tcell_mean = mean(var_df$variant[tcells])*100,
             nont_mean = mean(var_df$variant[!tcells])*100,
             mds_mean = mean(var_df$variant[mdf$chr7 == "Monosomy7"])*100,
             wt_mean = mean(var_df$variant[mdf$chr7 == "Wildtype"])*100)
}) %>% data.table::rbindlist() %>% data.frame() -> vars_df
var_meta_df<- cbind(var_meta_df_raw,vars_df)
var_meta_df$kruskal_pvalue_adj <- p.adjust(var_meta_df$kruskal_pvalue)
var_meta_df %>% mutate(rat = mds_mean/wt_mean) %>% arrange((rat))

plot_mutation <- function(variant){
  sub_df <- data.frame(
    UMAP1 = mdf$UMAP_1,
    UMAP2 = mdf$UMAP_2,
    mutation = mdf[,variant] * 100,
    chr7 = mdf$chr7
  )
  pm <- ggplot(sub_df %>% arrange((mutation)),
               aes(x =  UMAP1, y = UMAP2, color = mutation)) + geom_point()+ 
    pretty_plot() + L_border() + facet_wrap(~chr7)+
    scale_color_viridis() + 
    theme(legend.position = "bottom") + ggtitle(variant)
  
  cowplot::ggsave2(pm, file = paste0("../plots/all_muts/", variant, ".png"), 
                   width = 7, height = 4.5)
  variant
}
lapply(1:dim(var_meta_df)[1], function(i){
  colnames(mdf)[i+1] %>% plot_mutation() 
})

ggplot(var_meta_df, aes(x = log10(mds_mean/wt_mean), y = -log10(kruskal_pvalue_adj), label = variant)) + 
  geom_text() +
  pretty_plot() + L_border() +
  geom_vline(xintercept = 0, linetype = 2, color = "dodgerblue3")

ggplot(mdf %>% arrange((X3086T.C)),
       aes(x =  UMAP_1, y = UMAP_2, color = X3086T.C > 0.1)) + geom_point()
ggplot(mdf %>% arrange((X7836T.C)),
       aes(x =  UMAP_1, y = UMAP_2, color = X7836T.C > 0.1)) + geom_point()
ggplot(mdf %>% arrange((X13970G.A)),
       aes(x =  UMAP_1, y = UMAP_2, color = X13970G.A > 0.1)) + geom_point()
ggplot(mdf %>% arrange((X5043G.A)),
       aes(x =  UMAP_1, y = UMAP_2, color = X5043G.A > 0.1)) + geom_point()
ggplot(mdf %>% arrange((X14476G.A)),
       aes(x =  UMAP_1, y = UMAP_2, color = X14476G.A  > 0.1)) + geom_point()
