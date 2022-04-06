library(data.table)
library(ggbeeswarm)
library(BuenColors)
library(dplyr)

del_het <- fread("../data/all5Sample_dels.deletion_heteroplasmy.tsv")
assign_df <- fread("../output/firstpass_assignments_mix5.tsv")
bdf <- merge(del_het, assign_df, by.x = "cell_id", by.y = "barcode")
real_dels <- c("del13157-15477",  "del9232-13413","del8482-13445") 
bdf <- bdf %>% dplyr::filter(deletion %in% real_dels & fp_classify != "unassigned" & version == "improved")
bdf$deletion <- factor(as.character(bdf$deletion), real_dels)
bdf <- bdf %>% mutate(
  off_target = (!((fp_classify == "aPD1" & deletion == "del13157-15477") | 
                   (fp_classify == "aPD3" & deletion == "del8482-13445") |
                   (fp_classify == "aPD2" & deletion == "del9232-13413")) & heteroplasmy > 0)
)
p1 <- ggplot(bdf, aes(x = fp_classify, y = heteroplasmy, color = off_target)) +
  geom_quasirandom(size = 0.1) + facet_wrap(~deletion, ncol = 3) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = c("black", "firebrick")) +
  theme(legend.position = "none") + labs(x = "", y = "% Heteroplasmy")
#cowplot::ggsave2(p1, file = "../output/mix5_scatter.pdf", width = 3.5, height = 1.5)

bdf %>% dplyr::filter(!off_target) %>% group_by(fp_classify, deletion) %>% summarize(mean(heteroplasmy), median(heteroplasmy), count = n())
bdf %>% group_by(deletion) %>% summarize(sum(off_target, na.rm = TRUE), sum(reads_all))

#---
# Look at the 3 mix
#---
del_het <- fread("../data/all3Sample_dels.deletion_heteroplasmy.tsv")
assign_df <- fread("../output/firstpass_assignments_mix3.tsv")
bdf <- merge(del_het, assign_df, by.x = "cell_id", by.y = "barcode")
real_dels <- c("del13157-15477",  "del9232-13413","del8482-13445") 
bdf <- bdf %>% dplyr::filter(deletion %in% real_dels & fp_classify != "unassigned" & version == "improved")
bdf$deletion <- factor(as.character(bdf$deletion), real_dels)
bdf <- bdf %>% mutate(
  off_target = (!((fp_classify == "aPD1" & deletion == "del13157-15477") | 
                    (fp_classify == "aPD3" & deletion == "del8482-13445") |
                    (fp_classify == "aPD2" & deletion == "del9232-13413")) & heteroplasmy > 0)
)
p1 <- ggplot(bdf, aes(x = fp_classify, y = heteroplasmy, color = off_target)) +
  geom_quasirandom(size = 0.1) + facet_wrap(~deletion, ncol = 3) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = c("black", "firebrick")) +
  theme(legend.position = "none") + labs(x = "", y = "% Heteroplasmy")
#cowplot::ggsave2(p1, file = "../output/mix3_scatter.pdf", width = 3.5, height = 1.5)

bdf %>% dplyr::filter(!off_target) %>% group_by(fp_classify, deletion) %>% summarize(mean(heteroplasmy), median(heteroplasmy), count = n())
bdf %>% group_by(deletion) %>% summarize(sum(off_target, na.rm = TRUE), sum(reads_all))
