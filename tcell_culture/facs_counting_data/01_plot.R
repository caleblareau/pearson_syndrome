library(BuenColors)
library(dplyr)
library(data.table)

df <- fread("data_from_frank_josch.tsv") %>% data.frame()

p1 <- df %>%
  filter(experiment == "exp2") %>%
  ggplot(aes(x = Day, y = expansion, color = donor)) +
  geom_point() + geom_line() + scale_y_log10() +
  scale_color_manual(values = c("dodgerblue", "dodgerblue3", "dodgerblue4", "firebrick")) +
  pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "none") + labs(x = "Day", y = "Fold expansion")

p2 <- df[complete.cases(df),] %>%
  filter(experiment == "exp2") %>%
  ggplot(aes(x = Day, y = cd8_cd4_ratio,color = donor)) +
  geom_point() + geom_line() +
  scale_color_manual(values = c("dodgerblue", "dodgerblue3", "dodgerblue4", "firebrick")) +
  pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "none") + labs(x = "Day", y = "CD8+:CD4+ T cell ratio")
cowplot::ggsave2(p1, file = "expansion-fold-exp2.pdf", width = 1.8, height = 1.4)
cowplot::ggsave2(p2, file = "cd8_cd4_ratio-exp2.pdf", width = 1.8, height = 1.4)
