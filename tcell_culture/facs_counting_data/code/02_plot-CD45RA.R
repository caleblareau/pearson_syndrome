library(BuenColors)
library(dplyr)
library(data.table)

df <- fread("../data/cd45ra_data_from_frank.tsv") %>% data.frame()

p1 <- df %>%
  filter(Marker == "CD4")%>%
  ggplot(aes(x = Day, y = Amount, color = Donor)) +
  geom_point() + geom_line() + 
  scale_color_manual(values = c("dodgerblue", "dodgerblue3", "dodgerblue4", "firebrick")) +
  pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "none") + labs(x = "Day",  y = "% of CD4, CD45RA positive")

p2 <- df %>%
  filter(Marker == "CD8") %>%
  ggplot(aes(x = Day, y = Amount,color = Donor)) +
  geom_point() + geom_line() +
  scale_color_manual(values = c("dodgerblue", "dodgerblue3", "dodgerblue4", "firebrick")) +
  pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "none") +  labs(x = "Day", y = "% of CD8, CD45RA positive")
cowplot::ggsave2(p1, file = "../plots/main-cd45ra-over-time-CD4.pdf", width = 1.8, height = 1.4)
cowplot::ggsave2(p2, file = "../plots/main-cd45ra-over-time-CD8.pdf", width = 1.8, height = 1.4)
