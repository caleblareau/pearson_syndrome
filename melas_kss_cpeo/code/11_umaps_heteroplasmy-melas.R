library(data.table)
library(dplyr)
library(viridis)
library(BuenColors)

make_combined_df <- function(x){
  hetdf <- fread(paste0("../data/melas-metadata/",x,"_meta.tsv"))
  merge(
    hetdf, 
    fread(paste0("../data/reference_projections/melas//",x,"_refmapped.csv.gz")),
    by.y = "cb", by.x = "barcode"
  ) 
  
}
dfp21 <- make_combined_df("P21")
dfp9 <- make_combined_df("P9") 
dfp30 <- make_combined_df("P30")


ggplot(dfp21 %>% arrange((heteroplasmy)), aes(x = refUMAP_1, y = refUMAP_2, color = heteroplasmy)) +
  geom_point() + scale_color_gradientn(colors = c("lightgrey", "firebrick"))  + pretty_plot() +
  theme_void() + theme(legend.position = "bottom") +labs(color = "")

ggplot(dfp9 %>% arrange((heteroplasmy)), aes(x = refUMAP_1, y = refUMAP_2, color = heteroplasmy)) +
  geom_point() + scale_color_gradientn(colors = c("lightgrey", "firebrick"))  + pretty_plot() +
  theme_void() + theme(legend.position = "bottom") +labs(color = "")

ggplot(dfp30 %>% arrange((heteroplasmy)), aes(x = refUMAP_1, y = refUMAP_2, color = heteroplasmy)) +
  geom_point() + scale_color_gradientn(colors = c("lightgrey", "firebrick"))  + pretty_plot() +
  theme_void() + theme(legend.position = "bottom") +labs(color = "")

dfp21 %>% group_by(predicted.celltype.l2) %>%
  summarize(count = n(), pct_het = mean(heteroplasmy < 0.01)*100) %>%
  arrange(desc(pct_het)) %>% filter(count > 50)

dfp9 %>% group_by(predicted.celltype.l2) %>%
  summarize(count = n(), pct_het = mean(heteroplasmy < 0.01)*100) %>%
  arrange(desc(pct_het)) %>% filter(count > 50)

dfp30 %>% group_by(predicted.celltype.l2) %>%
  summarize(count = n(), pct_het = mean(heteroplasmy < 0.01)*100) %>%
  arrange(desc(pct_het)) %>% filter(count > 50)

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
dfp9 %>%
  filter(predicted.celltype.l2 %in% tcell_types) %>%
  group_by(predicted.celltype.l2) %>% summarize(count = n())
process_donor(dfp9)

dfp21 %>%
  filter(predicted.celltype.l2 %in% tcell_types) %>%
  group_by(predicted.celltype.l2) %>% summarize(count = n())
process_donor(dfp21)

px <- dfp21 %>%
  filter(predicted.celltype.l2 %in% tcell_types) %>%
  ggplot(aes(color = predicted.celltype.l2, x = heteroplasmy*100)) +
  stat_ecdf() + scale_color_manual(values = c("dodgerblue", "dodgerblue4", "pink", "firebrick", "purple4")) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))
cowplot::ggsave2(px, file = paste0("../plots/melas_p21_tcell_ecdf.pdf"), width = 2.3, height = 1.5)

px <- dfp9 %>%
  filter(predicted.celltype.l2 %in% tcell_types) %>%
  ggplot(aes(color = predicted.celltype.l2, x = heteroplasmy*100)) +
  stat_ecdf() + scale_color_manual(values = c("dodgerblue", "dodgerblue4", "pink", "firebrick", "purple4")) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "% Heteroplasmy", y = "Cumulative fraction", color = "") +
  scale_x_continuous(limits = c(0, 100))
cowplot::ggsave2(px, file = paste0("../plots/melas_p9_tcell_ecdf.pdf"), width = 2.3, height = 1.5)
