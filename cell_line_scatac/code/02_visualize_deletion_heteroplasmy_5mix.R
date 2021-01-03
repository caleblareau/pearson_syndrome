library(data.table)
library(ggbeeswarm)
library(BuenColors)

del_het <- fread("../data/tryall5Sample_dels.deletion_heteroplasmy.tsv")
assign_df <- fread("../output/firstpass_assignments.tsv")
bdf <- merge(del_het, assign_df, by.x = "cell_id", by.y = "barcode")

ggplot(bdf, aes(x = fp_classify, y = heteroplasmy, color = reads_del %in% c(1,2))) +
  geom_quasirandom() + facet_wrap(~deletion, ncol = 3) +
  pretty_plot() + L_border()

bdf %>% filter(fp_classify != "aPD1") %>% filter((heteroplasmy > 0) & (reads_del > 2) & (deletion == "del13157-15477")) %>%
  dplyr::filter(fp_classify != "unassigned") %>%
  pull(cell_id) -> hmm

bdf %>% filter(cell_id %in% hmm)
