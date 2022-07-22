library(Seurat)
library(Signac)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(GenomeInfoDb)
library(harmony)
library(BuenColors)
library(ggbeeswarm)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(viridis)

set.seed(1)

# Import ATAC data
summary_csv <- fread("../data/tcell-aggr-singlecell.csv.gz") %>%
  dplyr::filter(is__cell_barcode != "None")
summary_csv <- data.frame(summary_csv)
rownames(summary_csv) <- summary_csv[[1]]

# Add heteroplasmy
df_1 <- fread("../../pbmc_scatac/data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>% data.frame() %>% 
  dplyr::filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all)

# filter for tcells
tcells_pbmcs <- fread("../../pbmc_scatac/data/reference_projection/PT3_refmapped.csv.gz") %>% 
  dplyr::filter(predicted.celltype.l1 %in% c("CD4 T", "CD8 T", "other T")) %>% pull(cb)
df_1 <- df_1 %>% dplyr::filter(barcode %in% tcells_pbmcs)

df_2 <- fread("../data/deletion_heteroplasmy/PearsonCulture_Day14v1_del.deletion_heteroplasmy.tsv") %>% data.frame() %>% 
  dplyr::filter(deletion == "del10381-15407" & version == "improved") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all)

df_3 <- fread("../data/deletion_heteroplasmy/PearsonCulture_Day14v2_del.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407"& version == "improved") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all)

df_4 <- fread("../data/deletion_heteroplasmy/PearsonCulture_Day21_del.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407"& version == "improved") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-4", cell_id)) %>%
  dplyr::select(barcode, heteroplasmy, reads_del, reads_wt, reads_all)


# Make ecdf plot
het_df_plot <- rbind(df_1, df_2, df_3, df_4) %>% data.frame() 
cc <- substr(het_df_plot$barcode, 18, 18)
het_df_plot$day <- case_when(
  cc == "1" ~ "aPBMC",
  cc %in% c("2", "3") ~ "Day14",
  cc == "4" ~ "Day21",
)

het_df_plot %>% dplyr::filter(reads_all > 10) %>%
  group_by(day) %>% summarize(count = n())

p1 <- ggplot(het_df_plot %>% dplyr::filter(reads_all > 10), aes(x = heteroplasmy, color = day)) +
  stat_ecdf() + scale_color_manual(values = c("orange", jdb_palette("brewer_spectra")[c(7,9)])) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100)) + theme(legend.position = "none")

cowplot::ggsave2(p1, file = paste0("../plots/ecdf_tcell.pdf"), width = 1.5, height = 1.5)
