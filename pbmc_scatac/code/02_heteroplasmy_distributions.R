library(BuenColors)
library(data.table)
library(dplyr)
"%ni%" <- Negate("%in%")

df_BCI <- fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
  filter(deletion == "del6073-13095") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id))

df_CCF <- fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
  filter(deletion == "del8482-13447") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id))

df_PT3 <- fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
  filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id))

het_df <- rbind(df_BCI, df_CCF, df_PT3)

meta_data <- readRDS("../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")@meta.data


bdf <- merge(meta_data, het_df, by = "barcode") %>% filter(reads_all > 20)

exclude_celltypes <- c("CCF cell 1", "IFN-activated T-cells", "Innate-like B-cells")
ggplot(bdf %>% filter(predicted.id %ni% exclude_celltypes), aes(x = predicted.id, y = heteroplasmy, color = Patient)) +
  geom_quasirandom() + geom_boxplot(color = "black", outlier.shape = NA, fill = NA) + facet_wrap(~Patient, nrow = 3)
