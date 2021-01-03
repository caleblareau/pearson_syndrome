library(data.table)
library(ggbeeswarm)
library(BuenColors)

observed_patient_dels <- c("del13157-15477", "del9232-13413" ,"del8482-13445")
del_het <- fread("../data/tryall3Sample_dels.deletion_heteroplasmy.tsv")
assign_df <- fread("../output/firstpass_assignments_mix3.tsv")
bdf <- merge(del_het, assign_df, by.x = "cell_id", by.y = "barcode") %>%
  dplyr::filter(version == "improved" & deletion %in% observed_patient_dels)

#%>% dplyr::filter(fp_classify != "unassigned")
ggplot(bdf , aes(x = Cluster, y = heteroplasmy, color = reads_del %in% c(1,2))) +
  geom_quasirandom() + facet_wrap(~deletion, ncol = 3) +
  pretty_plot() + L_border()

