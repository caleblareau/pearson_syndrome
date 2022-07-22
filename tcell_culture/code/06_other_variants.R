library(SummarizedExperiment)
library(Matrix)
library(dplyr)
source("../../global_functions/variant_calling.R")

set.seed(1)

# Add heteroplasmy
df_1 <- fread("../../pbmc_scatac/data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>% data.frame() %>% 
  dplyr::filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all) %>% filter(reads_all > 10)

# filter for tcells
tcells_pbmcs <- fread("../../pbmc_scatac/data/reference_projection/PT3_refmapped.csv.gz") %>% 
  dplyr::filter(predicted.celltype.l1 %in% c("CD4 T", "CD8 T", "other T")) %>% pull(cb)
df_1 <- df_1 %>% dplyr::filter(barcode %in% tcells_pbmcs)

df_2 <- fread("../data/deletion_heteroplasmy/PearsonCulture_Day14v1_del.deletion_heteroplasmy.tsv") %>% data.frame() %>% 
  dplyr::filter(deletion == "del10381-15407" & version == "improved") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all)%>% filter(reads_all > 10)

df_3 <- fread("../data/deletion_heteroplasmy/PearsonCulture_Day14v2_del.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407"& version == "improved") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all)%>% filter(reads_all > 10)

df_4 <- fread("../data/deletion_heteroplasmy/PearsonCulture_Day21_del.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407"& version == "improved") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-4", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all)%>% filter(reads_all > 10)


import_mgatk_se <- function(file, df, idx){
  se <- readRDS(file)
  colnames(se) <- gsub("1", idx, colnames(se))
  se[,sort(df$barcode)]
}

# Call variants per library and then combin
pbmc <-   import_mgatk_se("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/Pearson-PBMC-PT3_mgatk.rds", df_1, "1") %>% call_mutations_mgatk
day14 <-   import_mgatk_se("../../../pearson_large_data_files/input/tcell-culture/mgatk/PearsonCulture_Day14v1_mgatk.rds", df_2, "2")%>% call_mutations_mgatk
day21 <- import_mgatk_se("../../../pearson_large_data_files/input/tcell-culture/mgatk/PearsonCulture_Day21_mgatk.rds", df_4, "4")%>% call_mutations_mgatk

# Call variants
pbmc_vars <-  data.frame(rowData(pbmc)) %>%  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 20) %>% pull(variant)
day14_vars <-  data.frame(rowData(day14)) %>%  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 20) %>% pull(variant)
day21_vars <-  data.frame(rowData(day21)) %>%  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 20) %>% pull(variant)
combined_vars <- unique(c(pbmc_vars, day14_vars, day21_vars))

# remove variants that are problematic
combined_vars <- combined_vars[!(combined_vars %in% c("310T>C", "10382A>C"))]

# Mark transitions
transition <- c("C>T", "G>A", "A>G", "T>C")
table(substr(combined_vars,nchar(combined_vars)-2, nchar(combined_vars) )%in% transition)

afaf <- data.frame(
  combined_vars,
  pbmc_af = data.frame(rowData(pbmc))[combined_vars,"mean"],
  day14_af = data.frame(rowData(day14))[combined_vars,"mean"],
  day21_af = data.frame(rowData(day21))[combined_vars,"mean"]
)
write.table(afaf, file = "../output/afs_called_variants.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
library(BuenColors)
p1 <- ggplot(afaf, aes(x = pbmc_af*100, y = day14_af*100)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_y_log10(limits = c(5e-03, 10)) + scale_x_log10(limits = c(5e-03, 10)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x ="PBMC T cell heteroplasmy", y= "Day 14 heteroplasmy")
cowplot::ggsave2(p1, file = "../plots/scatter_frequency.pdf", width = 2, height = 1.8)


afaf %>% mutate(fc = log2(day14_af/pbmc_af)) %>% 
  arrange(desc(fc))

saveRDS(assays(pbmc)[["allele_frequency"]][combined_vars,],
        file = "../output/pbmc_called_variants.rds")
saveRDS(assays(day14)[["allele_frequency"]][combined_vars,],
        file = "../output/day14v1_called_variants.rds")
saveRDS(assays(day21)[["allele_frequency"]][combined_vars,],
        file = "../output/day21_called_variants.rds")

day14v2 <-   import_mgatk_se("../../../pearson_large_data_files/input/tcell-culture/mgatk/PearsonCulture_Day14v2_mgatk.rds", df_3, "3")%>% call_mutations_mgatk
saveRDS(assays(day14v2)[["allele_frequency"]][combined_vars,],
        file = "../output/day14v2_called_variants.rds")

