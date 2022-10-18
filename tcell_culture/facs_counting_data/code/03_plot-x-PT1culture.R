library(BuenColors)
library(dplyr)
library(data.table)

df <- fread("../data/josch-pt1-culture.txt") %>% data.frame()

p1 <- df %>%
  filter(!is.na(Expansion)) %>%
  ggplot(aes(x = Day, y = Expansion, color = donor)) +
  geom_point() + geom_line() + scale_y_log10() +
  scale_color_manual(values = c("dodgerblue",  "firebrick")) +
  pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "none") + labs(x = "Day", y = "Fold expansion")
p1
p2 <- df%>%
  filter(!is.na(CD8_CD4_ratio)) %>%
  ggplot(aes(x = Day+2, y = CD8_CD4_ratio,color = donor)) +
  geom_point() + geom_line() +
  scale_color_manual(values = c("dodgerblue", "firebrick")) +
  pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "none") + labs(x = "Day", y = "CD8+:CD4+ T cell ratio")
p2

cowplot::ggsave2(p1, file = "../plots/PT1-expansion-fold.pdf", width = 1.5, height = 1.4)
cowplot::ggsave2(p2, file = "../plots/PT1-cd8_cd4_ratio.pdf", width = 1.5, height = 1.4)

