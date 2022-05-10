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

set.seed(42)

apd1  <- bdf %>% filter(fp_classify == "aPD1" & deletion == "del13157-15477")
apd1$permuted_heteroplasmy <- rbinom(n = length(apd1$reads_all), size = apd1$reads_all, prob = sum(apd1$reads_del)/sum(apd1$reads_all))/apd1$reads_all*100
apd2  <- bdf %>% filter(fp_classify == "aPD2" & deletion == "del9232-13413")
apd2$permuted_heteroplasmy <- rbinom(n = length(apd2$reads_all), size = apd2$reads_all, prob = sum(apd2$reads_del)/sum(apd2$reads_all))/apd2$reads_all*100
apd3  <- bdf %>% filter(fp_classify == "aPD3" & deletion == "del8482-13445")
apd3$permuted_heteroplasmy <- rbinom(n = length(apd3$reads_all), size = apd3$reads_all, prob = sum(apd3$reads_del)/sum(apd3$reads_all))/apd3$reads_all*100

var(apd1$heteroplasmy);var(apd1$permuted_heteroplasmy)
var(apd2$heteroplasmy);var(apd2$permuted_heteroplasmy)
var(apd3$heteroplasmy);var(apd3$permuted_heteroplasmy)

(var(apd1$heteroplasmy)-var(apd1$permuted_heteroplasmy))/var(apd1$heteroplasmy)
(var(apd2$heteroplasmy)-var(apd2$permuted_heteroplasmy))/var(apd2$heteroplasmy)
(var(apd3$heteroplasmy)-var(apd3$permuted_heteroplasmy))/var(apd3$heteroplasmy)

rdf <- rbind(
  reshape2::melt(apd1[,c("heteroplasmy", "permuted_heteroplasmy")]) %>% 
    mutate(deletion = "aPD1"),
  reshape2::melt(apd2[,c("heteroplasmy", "permuted_heteroplasmy")]) %>% 
    mutate(deletion = "aPD2"),
  reshape2::melt(apd3[,c("heteroplasmy", "permuted_heteroplasmy")]) %>% 
    mutate(deletion = "aPD3")
)

p1 <- ggplot(rdf, aes(x = deletion, y = value, color = variable)) +
  geom_violin(fill = NA) + scale_color_manual(values = c("black", "forestgreen")) +
  pretty_plot(fontsize = 7) + L_border() + theme(legend.position = "none") + 
  labs( x= "Deletion", y = "Heteroplasmy (%)")
cowplot::ggsave2(p1, file = "../output/permuted_VE.pdf", width = 1.8, height = 1.6)
