library(data.table)
library(ggbeeswarm)
library(BuenColors)

del_het <- fread("Pearson-Mix-3_v12-mtMask_mgatk/final/threeSample_dels.deletion_heteroplasmy.tsv")
assign_df <- fread("output/firstpass_assignments_mix3.tsv")
bdf <- merge(del_het, assign_df, by.x = "cell_id", by.y = "barcode")

ggplot(bdf, aes(x = fp_classify, y = heteroplasmy, color = reads_del %in% c(1,2))) +
  geom_quasirandom() + facet_wrap(~deletion, ncol = 3) +
  pretty_plot() + L_border()

