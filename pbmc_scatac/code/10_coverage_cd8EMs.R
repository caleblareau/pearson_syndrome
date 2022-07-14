library(data.table)
library(dplyr)
library(viridis)
library(BuenColors)
library(SummarizedExperiment)
library(Matrix)

# Add Heteroplasmy
df_BCI <- merge(fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
                  dplyr::filter(deletion == "del6073-13095") ,
                fread(paste0("../data/reference_projection/BCI_refmapped.csv.gz")),
                by.y = "cb", by.x = "cell_id") %>%filter(reads_all >= 10) 

df_CCF <- merge(fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
                  dplyr::filter(deletion == "del8482-13447") ,
                fread(paste0("../data/reference_projection/CCF_refmapped.csv.gz")),
                by.y = "cb", by.x = "cell_id") %>%filter(reads_all >= 10) 

df_PT3 <- merge(fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
                  dplyr::filter(deletion == "del10381-15407") ,
                fread(paste0("../data/reference_projection/PT3_refmapped.csv.gz")),
                by.y = "cb", by.x = "cell_id") %>%filter(reads_all >= 10) 

tcell_types <- c("MAIT", "CD8 TEM", "CD4 TCM", "CD4 Naive", "CD8 Naive")

bci_bc <- df_BCI %>%
  dplyr::filter(predicted.celltype.l2 == "CD8 TEM") %>%
  dplyr::filter(heteroplasmy < 0.01) %>% pull(cell_id)

ccf_bc <- df_CCF %>%
  dplyr::filter(predicted.celltype.l2 == "CD8 TEM") %>%
  dplyr::filter(heteroplasmy < 0.01) %>% pull(cell_id)

pt3_bc <- df_PT3 %>%
  dplyr::filter(predicted.celltype.l2 == "CD8 TEM") %>%
  dplyr::filter(heteroplasmy < 0.01) %>% pull(cell_id)

bci_mean_purified <- rowMeans(assays(readRDS("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/Pearson-PBMC-BCI_mgatk.rds"))[["coverage"]][,bci_bc])
ccf_mean_purified <- rowMeans(assays(readRDS("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/Pearson-PBMC-CCF_mgatk.rds"))[["coverage"]][,ccf_bc])
pt3_mean_purified <- rowMeans(assays(readRDS("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/Pearson-PBMC-PT3_mgatk.rds"))[["coverage"]][,pt3_bc])

bci_bc2 <- df_BCI %>%
  dplyr::filter(predicted.celltype.l2 == "CD8 TEM") %>%
  dplyr::filter(heteroplasmy > 30) %>% pull(cell_id)

ccf_bc2 <- df_CCF %>%
  dplyr::filter(predicted.celltype.l2 == "CD8 TEM") %>%
  dplyr::filter(heteroplasmy >30) %>% pull(cell_id)

pt3_bc2 <- df_PT3 %>%
  dplyr::filter(predicted.celltype.l2 == "CD8 TEM") %>%
  dplyr::filter(heteroplasmy > 30) %>% pull(cell_id)

bci_mean_there <- rowMeans(assays(readRDS("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/Pearson-PBMC-BCI_mgatk.rds"))[["coverage"]][,bci_bc2])
ccf_mean_there <- rowMeans(assays(readRDS("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/Pearson-PBMC-CCF_mgatk.rds"))[["coverage"]][,ccf_bc2])
pt3_mean_there <- rowMeans(assays(readRDS("../../../pearson_large_data_files/input/pbmcs_scatac/mgatk-output/Pearson-PBMC-PT3_mgatk.rds"))[["coverage"]][,pt3_bc2])

pb1 <- qplot(1:length(bci_mean_purified), bci_mean_purified)
pc1 <- qplot(1:length(ccf_mean_purified), ccf_mean_purified)
pp1 <- qplot(1:length(pt3_mean_purified), pt3_mean_purified)

pb2 <- qplot(1:length(bci_mean_there), bci_mean_there)
pc2 <- qplot(1:length(ccf_mean_there), ccf_mean_there)
pp2 <- qplot(1:length(pt3_mean_there), pt3_mean_there)

library(patchwork)
pb1 + pc1 + pp1 +
pb2 + pc2 + pp2

df_PT3 %>%
  group_by(predicted.celltype.l2) %>%
  summarize(mean_covearge = mean(reads_all), count = n()) %>%
  filter(count > 50)


