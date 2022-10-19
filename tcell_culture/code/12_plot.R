library(data.table)
library(BuenColors)
library(dplyr)
library(viridis)
sc <- fread("singlecell.csv") %>%
  filter(is__cell_barcode == 1)
mdf <- merge(rbind(
  fread("Pearson_HC_PT1_mix_CD3_CD28_Day14_rep1_assign.tsv"),
  fread("Pearson_HC_PT1_mix_CD3_CD28_Day14_rep2_assign.tsv") %>%
    mutate(barcode = paste0(substr(barcode, 1, 16), "-2"))
), sc, by = "barcode") %>% mutate(frip = peak_region_fragments/passed_filters*100)

ggplot(mdf, aes(x = log10(passed_filters), y = frip, color = assign))   + pretty_plot() + L_border() + 
  geom_point(size = 0.4) + scale_color_manual(values = jdb_palette("corona")) +
  facet_wrap(~assign)

ggplot(mdf %>% filter(assign == "Pearson"), aes(x = log10(passed_filters), y = frip, color = heteroplasmy))   + pretty_plot() + L_border() + 
  geom_point(size = 0.4) + scale_color_viridis()

mdf %>% filter(assign == "Pearson") %>% mutate(highFRIP = frip > 30)  %>%
  ggplot(aes(x=highFRIP, y = heteroplasmy)) + geom_violin()
