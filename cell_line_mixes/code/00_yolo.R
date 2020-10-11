library(data.table)
library(ggbeeswarm)
library(BuenColors)

"%ni%" <- Negate("%in%")
del_het <- fread("Pearson-Mix-5_v12-mtMask_mgatk/final/tryall5Sample_dels.deletion_heteroplasmy.tsv")
colnames(del_het) <- c("cell_id","heteroplasmy","reads_del","reads_wt","reads_all","bp_j","deletion", "what")

assign_df <- fread("output/firstpass_assignments.tsv")
bdf <- merge(del_het, assign_df, by.x = "cell_id", by.y = "barcode")

ggplot(bdf %>% filter(deletion %in% c("del13157-15477", "del9232-13413" ,"del8482-13445") & what == "improved"),
       aes(x = classify, y = heteroplasmy, color = (reads_del == 1))) +
  geom_quasirandom() + facet_grid(~deletion) + 
  pretty_plot() + L_border() + labs(x = "Donor annotation", y = "Heteroplasmy %", color = "# reads for del")




bdf %>% filter(deletion %ni% c("del13157-15477", "del9232-13413" ,"del8482-13445") & what == "test") %>%
  group_by(param) %>% summarize(count = sum(heteroplasmy > 0, na.rm = TRUE)) 

boo <- which((bdf[c(TRUE, FALSE),2] < 0.1) & (bdf[c(FALSE, TRUE),2] > 0.1))
head(bdf[c(TRUE, FALSE),][boo,])
head(bdf[c(FALSE, TRUE),][boo,])


bdf %>% group_by(deletion, fp_classify, what) %>% summarize(median(heteroplasmy), mean(heteroplasmy)) %>% data.frame()
