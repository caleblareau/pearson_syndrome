library(data.table)
library(dplyr)
library(viridis)
library(BuenColors)

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

do_binomial_test <- function(df, ct1, ct2){
  ct1v <- df %>% filter(predicted.celltype.l2 %in% ct1) %>%
    mutate(x = heteroplasmy < 0.01) %>% pull(x)
  ct2v <- df %>% filter(predicted.celltype.l2 %in% ct2) %>%
    mutate(x = heteroplasmy < 0.01) %>% pull(x)
  
  data.frame(
    ct1, ct2, 
    pvalue = prop.test(x = c(sum(ct1v), sum(ct2v)), n = c(length(ct1v), length(ct2v)))$p.value
  )
}

process_donor <- function(df){
  rbind(
    do_binomial_test(df, "CD8 TEM", "CD8 Naive"),
    do_binomial_test(df, "CD4 TCM", "CD4 Naive"),
    do_binomial_test(df,  "MAIT", c("CD8 Naive", "CD4 Naive", "CD8 TEM", "CD4 TEM"))
  ) 
}

process_donor(df_BCI)
process_donor(df_CCF)
process_donor(df_PT3)


px <- df_BCI %>%
  filter(predicted.celltype.l2 %in% tcell_types) %>%
  ggplot(aes(color = predicted.celltype.l2, x = heteroplasmy)) +
  stat_ecdf() + scale_color_manual(values = c("dodgerblue", "dodgerblue4", "pink", "firebrick", "purple4")) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))

cowplot::ggsave2(px, file = paste0("../plots/pBCI_tcell_ecdf.pdf"), width = 2.3, height = 1.5)

px <- df_CCF%>%
  filter(predicted.celltype.l2 %in% tcell_types) %>%
  ggplot(aes(color = predicted.celltype.l2, x = heteroplasmy)) +
  stat_ecdf() + scale_color_manual(values = c("dodgerblue", "dodgerblue4", "pink", "firebrick", "purple4")) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))

cowplot::ggsave2(px, file = paste0("../plots/pCCF_tcell_ecdf.pdf"),width = 2.3, height = 1.5)

px <- df_PT3 %>%
  filter(predicted.celltype.l2 %in% tcell_types) %>%
  ggplot(aes(color = predicted.celltype.l2, x = heteroplasmy)) +
  stat_ecdf() + scale_color_manual(values = c("dodgerblue", "dodgerblue4", "pink", "firebrick", "purple4")) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))
cowplot::ggsave2(px, file = paste0("../plots/pPT3_tcell_ecdf.pdf"), width = 2.3, height = 1.5)

##########
pt1het <- ggplot(shuf(df_BCI), aes(x = refUMAP_1, y = refUMAP_2, color = heteroplasmy)) +
  geom_point(size = 1) + scale_color_viridis()  +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(pt1het, file = "../plots/bci_heteroplasmy.png", width = 5, height = 5.1, dpi = 500)

pt2het <- ggplot(df_CCF %>% arrange(desc(heteroplasmy)), aes(x = refUMAP_1, y = refUMAP_2, color = heteroplasmy)) +
  geom_point(size = 1) + scale_color_viridis()  +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(pt2het, file = "../plots/ccf_heteroplasmy.png", width = 5, height = 5.1, dpi = 500)

pt3het <- ggplot(shuf(df_PT3), aes(x = refUMAP_1, y = refUMAP_2, color = heteroplasmy)) +
  geom_point(size = 1) + scale_color_viridis()  +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(pt3het, file = "../plots/pt3_heteroplasmy.png", width = 5, height = 5.1, dpi = 500)

df_CCF %>%
  filter(predicted.celltype.l2 %in% tcell_types) %>%
  group_by(predicted.celltype.l2) %>% summarize(count = n())
