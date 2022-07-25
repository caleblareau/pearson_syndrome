library(data.table)
library(dplyr)
library(viridis)
library(BuenColors)

make_combined_df <- function(x){
  hetdf <- fread(paste0("../data/namdc-metadata/",x,"_del.new.tsv")) %>%
    filter(version == "improved")
  merge(
    hetdf, 
    fread(paste0("../data/reference_projections/namdc/",x,"_refmapped.csv.gz")),
    by.y = "cb", by.x = "cell_id"
  ) %>%filter(reads_all >= 5) -> df
  df
}
df550 <- make_combined_df("NAMDC1412802550")
df671 <- make_combined_df("NAMDC20008038671") 
df507 <- make_combined_df("NAMDC18207019507")


sum(df550$reads_del)/sum(df550$reads_all)*100
sum(df671$reads_del)/sum(df671$reads_all)*100
sum(df507$reads_del)/sum(df507$reads_all)*100

df550 %>% mutate(nzh = heteroplasmy > 0) %>%
  group_by(predicted.celltype.l2) %>%
  summarize(prop = sum(nzh)/length(nzh)*100, n = sum(nzh)) %>%
  filter(prop > 0) %>% arrange(desc(prop))

df671 %>% mutate(nzh = heteroplasmy > 0) %>%
  group_by(predicted.celltype.l2) %>%
  summarize( n = sum(nzh)) %>%
  filter(n > 0) %>% arrange(desc(n))

df507 %>% mutate(nzh = heteroplasmy > 0) %>%
  group_by(predicted.celltype.l2) %>%
  summarize(prop = sum(nzh)/length(nzh)*100, n = sum(nzh)) %>%
  filter(prop > 0) %>% arrange(desc(prop))



df671$heteroplasmy <- ifelse(df671$heteroplasmy  > 50, 50, df671$heteroplasmy )
p1 <- ggplot(df671 %>% arrange((heteroplasmy)), aes(x = refUMAP_1, y = refUMAP_2, color = heteroplasmy)) +
  geom_point(size = 3) + scale_color_gradientn(colors = c("lightgrey", "firebrick"))  + pretty_plot() +
  theme_void() + theme(legend.position = "none") +labs(color = "")

df550$heteroplasmy <- ifelse(df550$heteroplasmy  > 50, 50, df550$heteroplasmy )
p2 <- ggplot(df550 %>% arrange((heteroplasmy)), aes(x = refUMAP_1, y = refUMAP_2, color = heteroplasmy)) +
  geom_point(size = 3) + scale_color_gradientn(colors = c("lightgrey", "firebrick"))  + pretty_plot() +
  theme_void() + theme(legend.position = "none") +labs(color = "")

df507 <- make_combined_df("NAMDC18207019507")
df507$heteroplasmy <- ifelse(df507$heteroplasmy  > 50, 50, df507$heteroplasmy )
p3 <- ggplot(df507 %>% arrange((heteroplasmy)), aes(x = refUMAP_1, y = refUMAP_2, color = heteroplasmy)) +
  geom_point(size = 3) + scale_color_gradientn(colors = c("lightgrey", "firebrick"))  + pretty_plot() +
  theme_void() + theme(legend.position = "none") +labs(color = "")

#cowplot::ggsave2(p3, file = "../plots/namdc_df507.png", width = 8, height = 8, dpi = 400)
#cowplot::ggsave2(p1, file = "../plots/namdc_df671.png", width = 8, height = 8, dpi = 400)
#cowplot::ggsave2(p2, file = "../plots/namdc_df550.png", width = 8, height = 8, dpi = 400)

# one patient had a special additional deletion
make_combined_df_special <- function(x){
  hetdf <- fread(paste0("../data/namdc-metadata/",x,"-smallDeletion_del.new.tsv")) %>%
    filter(version == "improved")
  merge(
    hetdf, 
    fread(paste0("../data/reference_projections/namdc/",x,"_refmapped.csv.gz")),
    by.y = "cb", by.x = "cell_id"
  ) %>%filter(reads_all >= 5) -> df
  df
}

df507 <- make_combined_df_special("NAMDC18207019507")
df507$heteroplasmy <- ifelse(df507$heteroplasmy  > 50, 50, df507$heteroplasmy )
df507 %>% filter(heteroplasmy > 0) 
p4 <- ggplot(df507 %>% arrange((heteroplasmy)), aes(x = refUMAP_1, y = refUMAP_2, color = heteroplasmy)) +
  geom_point(size = 3) + scale_color_gradientn(colors = c("lightgrey", "firebrick"))  + pretty_plot() +
  theme_void() + theme(legend.position = "none") +labs(color = "")
cowplot::ggsave2(p4, file = "../plots/namdc_small-del.png", width = 8, height = 8, dpi = 400)

df671 %>% group_by(predicted.celltype.l2) %>%
  summarize(count = n(), pct_nonzero_het = mean(heteroplasmy > 0)*100) %>%
  arrange(desc(pct_nonzero_het)) %>% filter(count > 50)
