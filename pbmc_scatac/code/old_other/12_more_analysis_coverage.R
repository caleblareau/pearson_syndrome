library(dplyr)
library(data.table)
library(BuenColors)
set.seed(1)

df_BCI <- fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
  filter(reads_all > 10) %>%
  dplyr::filter(deletion == "del6073-13095") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id))

df_CCF <- fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
  filter(reads_all > 10) %>%
  dplyr::filter(deletion == "del8482-13447") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id)) 

df_PT3 <- fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
  filter(reads_all > 10) %>%
  dplyr::filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id)) 

ggplot(df_PT3, aes(x = reads_all, y = heteroplasmy)) +
  geom_point() + scale_x_log10()

ggplot(df_BCI, aes(x = reads_all, y = heteroplasmy)) +
  geom_point() 
