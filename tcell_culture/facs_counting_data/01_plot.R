library(BuenColors)
library(dplyr)
library(data.table)

df <- fread("data_from_frank_josch.tsv") %>% data.frame()

df %>%
  ggplot(aes(x = Day, y = expansion, shape = Cell_type, color = donor)) +
  geom_point() + geom_line() + scale_y_log10()

df[complete.cases(df),] %>%
  filter(donor %in% c("Do21", "Do23", "DoCL", "DoPS")) %>%
  ggplot(aes(x = Day, y = cd8_cd4_ratio, shape = Cell_type, color = donor)) +
  geom_point() + geom_line() 

