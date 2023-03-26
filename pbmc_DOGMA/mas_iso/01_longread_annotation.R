library(SummarizedExperiment)
library(Matrix)
library(data.table)
library(stringi)
library(dplyr)
library(BuenColors)

barcodes_channel1 <- substr(fread("../data/DOGMA_PT1_Rep1_mask_mgatk.depthTable.txt")[[1]], 1, 16)
barcodes_channel2 <- substr(fread("../data/DOGMA_PT1_Rep2_mask_mgatk.depthTable.txt")[[1]], 1, 16)

# estimate deletion heteroplasmy from dna
estimate_coverage_heteroplasmy <- function(coord1, coord2, mat){
  idx <- 1:dim(mat)[1]
  in_del_boo <- (coord1 <= idx & coord2 >= idx)
  in_del <- colSums(mat[in_del_boo,])/sum(in_del_boo)
  out_del <- colSums(mat[!in_del_boo,])/sum(!in_del_boo)
  x <- in_del/out_del
  ifelse(x >1, 0, 1-x ) *100
}
se1 <- readRDS("../../../pearson_large_data_files/input/pbmcs_dogma/DOGMA_PT1_Rep1_mask_mgatk.rds")
het1 <- estimate_coverage_heteroplasmy(6073,13095, assays(se1)[["coverage"]])

se2 <- readRDS("../../../pearson_large_data_files/input/pbmcs_dogma/DOGMA_PT1_Rep2_mask_mgatk.rds")
het2 <- estimate_coverage_heteroplasmy(6073,13095,assays(se2)[["coverage"]])

#######
df1 <- data.frame(
  cell = names(het1),
  cov_het = unname(het1), 
  coverage = se1$depth
)

df2 <- data.frame(
  cell = gsub("-1", "-2", names(het2)),
  cov_het = unname(het2), 
  coverage = se2$depth
)

# Import mas-iso-seq data
build_df_fwd<- function(what, where, rep){
  reads <- fread(paste0("data/",what,"_",where, "_", rep, ".txt.gz"))[["V10"]]
  parsed_bc <- substr(reads, 24, 24+15)
  parsed_umi <- substr(reads, 24+16, 24+16+11)
  data.frame(parsed_bc, parsed_umi, what,  where, nchar = nchar(reads))
}
rbind(
  build_df_fwd("intact1", "fwd", "1"),
  build_df_fwd("intact2", "fwd", "1"),
  build_df_fwd("deletion", "fwd", "1")
) -> raw_df

raw_df %>% group_by(what) %>% 
  summarize(median(nchar))

ggplot(raw_df %>% filter(nchar < 2500 & nchar > 100), aes(x = nchar, color = what)) +
  stat_ecdf() +
  labs (x = "cDNA length") +
  scale_color_manual(values = c("firebrick","dodgerblue4", "goldenrod")) +
  pretty_plot(fontsize = 7) + theme(legend.position = "none") + 
  labs(y = "Cumulative density") -> p1

cowplot::ggsave2(p1, file = "output/ecdf.pdf", width = 1.45, height = 1.45)


# Use this to estimate pseudobulk heteroplasmy
raw_df %>% group_by(what, parsed_bc, parsed_umi) %>% 
  summarize(n = n()) %>%
  group_by(what) %>% summarize(count = n())


plothetdf <- data.frame(
  donor = c("PT1", "PT1", "PT1"),
  modality = c("scRNA-seq", "mtscATAC-seq", "MAS-ISO-seq"),
  heteroplasmy = c(22.48, 41.2, 1031/(1031+1120+959)*100)
)

p1 <- ggplot(plothetdf, aes(x = modality, y = heteroplasmy,)) +
  geom_bar(stat = "identity", position="dodge", color = "black", fill = "lightgrey", width = 0.5) + 
  labs(x = "",  y = "% Heteroplasmy") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) 
p1
cowplot::ggsave2(p1, file = "output/heteroplasmy-estimates-masiso-compare.pdf", width = 1, height = 1.5)


# Now estimate at single cell leevel
raw_df %>% group_by(what, parsed_bc, parsed_umi) %>% 
  summarize(n = n()) %>%
  filter(parsed_bc %in% barcodes_channel1) %>%
  mutate(barcode = paste0(parsed_bc, "-1")) %>%
  mutate(new_what = ifelse(what == "deletion", "del", "wt")) %>%
  reshape2::dcast(barcode ~ new_what, fill = 0) -> mas_iso_ID

clipped_reads_df <- rbind(
  fread("../data/DOGMA_PT1_Rep1_delquant.deletion_heteroplasmy.tsv"),
  fread("../data/DOGMA_PT1_Rep2_delquant.deletion_heteroplasmy.tsv") %>% mutate(cell_id = gsub("-1", "-2", cell_id))) %>%
  dplyr::filter(version == "naive") 
colnames(clipped_reads_df) <- c("cell_id", "clip_het", colnames(clipped_reads_df)[3:8])

library(BuenColors)
library(ggbeeswarm)
merge(mas_iso_ID, df1, by.x = "barcode", by.y = "cell") %>%
  filter(!(del>0 & wt>0)) %>%
  mutate(anno = ifelse(del >0, "del", "wt")) -> plot_df
wilcox.test(cov_het~del, plot_df)

p1 <- plot_df %>% ggplot(aes(x = anno, y = cov_het)) + 
  geom_boxplot(outlier.shape = NA) +
  pretty_plot(fontsize = 7) + 
  L_border() + labs(x = "Fusion transcript detected", y = "DNA-based heteroplasmy")
p1
cowplot::ggsa


mdff <- merge(mas_iso_ID, df1, by.x = "barcode", by.y = "cell", all.y = TRUE)
mdff_filt <- mdff %>% filter(coverage > 10)
mdff_filt[complete.cases(mdff_filt),]
