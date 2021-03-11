library(data.table)
library(dplyr)
library(BuenColors)
library(Seurat)
library(Matrix)
library(ggbeeswarm)

# Import GTEX
SA <- fread("../data/bulk_gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
blood_ids <- SA[SA$SMTS == "Blood",] %>% pull(SAMPID)
gtex_exp <- fread("../data/bulk_gtex/GTEx_expression-4genes.tsv")
blood_gtex <- gtex_exp %>% filter(BARCODE %in% blood_ids)
blood_gtex$what <- ("GTEX")

# Import pearson
bci <- rowSums(Read10X(data.dir = "../data/pearson_bci"))
ccf <- rowSums(Read10X(data.dir = "../data/pearson_ccf"))
mds <- rowSums(Read10X(data.dir = "../data/pearson_mds"))

p1 <- (bci/sum(bci)*1000000)[c("CXCL14", "GPR85", "LINC01641", "C12orf54")]
p2 <- (ccf/sum(ccf)*1000000)[c("CXCL14", "GPR85", "LINC01641", "C12orf54")]
p3 <- (mds/sum(mds)*1000000)[c("CXCL14", "GPR85", "LINC01641", "C12orf54")]

pearson_df <- data.frame(
  BARCODE = c("BCI", "CCF", "PT3"),
  CXCL14 = c(p1[1], p2[1], p3[1]),
  GPR85 = c(p1[2], p2[2], p3[2]),
  LINC01641 = c(p1[3], p2[3], p3[3]),
  C12orf54 = c(p1[4], p2[4], p3[4]),
  what = "Pearson"
)
all_df <- rbind(pearson_df, blood_gtex)


p1 <- ggplot(all_df, aes(x = what, y = CXCL14)) +
  geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
  geom_quasirandom(data = pearson_df, color = "firebrick", size = 0.3) +
  pretty_plot(fontsize = 6) + L_border() + labs(x = "", y = "CXCL14 (per million)")

p2 <- ggplot(all_df, aes(x = what, y = GPR85)) +
  geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
  geom_quasirandom(data = pearson_df, color = "firebrick",  size = 0.3) +
  pretty_plot(fontsize = 6) + L_border() + labs(x = "", y = "GPR85 (per million)")

cowplot::ggsave2(cowplot::plot_grid(p1, p2, ncol = 1), 
                 file = "../plots/GTEx_GPR85_CXCL14.pdf", width = 1.2, height = 2.2)

p1 <- ggplot(all_df, aes(x = what, y =  C12orf54 )) +
  geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
  geom_quasirandom(data = pearson_df, color = "firebrick", size = 0.3) +
  pretty_plot(fontsize = 6) + L_border() + labs(x = "", y = "C12orf54 (per million)")

p2 <- ggplot(all_df, aes(x = what, y = LINC01641)) +
  geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
  geom_quasirandom(data = pearson_df, color = "firebrick",  size = 0.3) +
  pretty_plot(fontsize = 6) + L_border() + labs(x = "", y = "LINC01641 (per million)")

cowplot::ggsave2(cowplot::plot_grid(p1, p2, ncol = 1), 
                 file = "../plots/GTEx_other2.pdf", width = 1.2, height = 2.2)


