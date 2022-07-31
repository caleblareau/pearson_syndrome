library(BuenColors)
library(dplyr)
library(data.table)

df <- fread("../data/count_cd48_data_from_frank_josch.tsv") %>% data.frame()

p1 <- df %>%
  filter(experiment == "exp1") %>%
  ggplot(aes(x = Day+2, y = expansion, color = donor, shape = Cell_type)) +
  geom_point() + geom_line() + scale_y_log10() +
  scale_color_manual(values = c("dodgerblue",  "red")) +
  pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "none") + labs(x = "Day", y = "Fold expansion")
p1
p2 <- df[complete.cases(df),] %>%
  filter(experiment == "exp1") %>%
  ggplot(aes(x = Day+2, y = cd8_cd4_ratio,color = donor, shape = Cell_type)) +
  geom_point() + geom_line() +
  scale_color_manual(values = c("dodgerblue", "red")) +
  pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "none") + labs(x = "Day", y = "CD8+:CD4+ T cell ratio")
p2

cowplot::ggsave2(p1, file = "../plots/supplement-expansion-fold.pdf", width = 1.5, height = 1.4)
cowplot::ggsave2(p2, file = "../plots/supplement-cd8_cd4_ratio.pdf", width = 1.5, height = 1.4)


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
cowplot::ggsave2(p1, file = "../plots/expansion-fold-main.pdf", width = 1.8, height = 1.4)
cowplot::ggsave2(p2, file = "../plots/cd8_cd4_ratio-main.pdf", width = 1.8, height = 1.4)
