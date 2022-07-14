library(data.table)
library(dplyr)
library(viridis)
library(BuenColors)

reanno_class <- function(df){
  df$caleb_class <- case_when(df$predicted.celltype.l1 %in% c("CD4 T", "CD8 T", "other T") ~ "Tcell",
                              df$predicted.celltype.l1 %in% c("Mono") ~ "Mono",
                              df$predicted.celltype.l1 %in% c("DC") ~ "DC",
                              df$predicted.celltype.l1 %in% c("B") ~ "Bcell",
                              df$predicted.celltype.l1 %in% c("NK") ~ "NKcell",
                              TRUE ~ "meh"
  )
  df %>% filter(caleb_class != "meh")
}

# Add Heteroplasmy
df_BCI <- merge(fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
                  dplyr::filter(deletion == "del6073-13095") ,
                fread(paste0("../data/reference_projection/BCI_refmapped.csv.gz")),
                by.y = "cb", by.x = "cell_id") %>%filter(reads_all >= 10)  %>% reanno_class


df_CCF <- merge(fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
                  dplyr::filter(deletion == "del8482-13447") ,
                fread(paste0("../data/reference_projection/CCF_refmapped.csv.gz")),
                by.y = "cb", by.x = "cell_id") %>%filter(reads_all >= 10)  %>% reanno_class

df_PT3 <- merge(fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
                  dplyr::filter(deletion == "del10381-15407") ,
                fread(paste0("../data/reference_projection/PT3_refmapped.csv.gz")),
                by.y = "cb", by.x = "cell_id") %>%filter(reads_all >= 10)  %>% reanno_class




do_binomial_test <- function(df, ct1, ct2){
  ct1v <- df %>% filter(caleb_class %in% ct1) %>%
    mutate(x = heteroplasmy < 0.01) %>% pull(x)
  ct2v <- df %>% filter(caleb_class %in% ct2) %>%
    mutate(x = heteroplasmy < 0.01) %>% pull(x)
  
  data.frame(
    ct1, ct2, 
    pvalue = prop.test(x = c(sum(ct1v), sum(ct2v)), n = c(length(ct1v), length(ct2v)))$p.value
  )
}

process_donor <- function(df){
  rbind(
    do_binomial_test(df, "Tcell", "Bcell")  ) 
}

process_donor(df_BCI)
process_donor(df_CCF)
process_donor(df_PT3)


px <- df_BCI %>%
  ggplot(aes(color = caleb_class, x = heteroplasmy)) +
  stat_ecdf() + scale_color_manual(values = c("orange2", "green4", "dodgerblue3","purple2", "firebrick")) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))
cowplot::ggsave2(px, file = paste0("../plots/pBCI_heme_ecdf.pdf"), width = 2.3, height = 1.5)


px <- df_CCF%>%
  ggplot(aes(color = caleb_class, x = heteroplasmy)) +
  stat_ecdf() + scale_color_manual(values = c("orange2", "green4", "dodgerblue3","purple2", "firebrick")) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))
cowplot::ggsave2(px, file = paste0("../plots/pCCF_heme_ecdf.pdf"),width = 2.3, height = 1.5)

px <- df_PT3 %>%
  ggplot(aes(color = caleb_class, x = heteroplasmy)) +
  stat_ecdf() + scale_color_manual(values = c("orange2", "green4", "dodgerblue3","purple2", "firebrick")) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))
cowplot::ggsave2(px, file = paste0("../plots/pPT3_heme_ecdf.pdf"), width = 2.3, height = 1.5)
