library(data.table)
library(dplyr)
library(BuenColors)

ps <- gsub("_del_find.SA_coverage_viz.pdf", "", list.files("pearson-data/", pattern = "*.SA_coverage_viz.pdf"))

lapply(ps, function(pearsonID){
  fread(paste0("pearson-data/",pearsonID,"_del_find.clip.tsv")) %>% 
    filter(SA > 5) %>% mutate(pearsonID)
}) %>% rbindlist() -> clip_df 

clip_df %>% group_by(position) %>%
  summarize(count = n(), total_SA = sum(SA), total_clip = sum(clip_count)) %>%
  filter(count  %in% c(1,2,3) & total_SA > 10 & total_clip > 10) %>% 
  arrange(desc(total_SA)) %>% data.frame() %>% head(50) %>%
  filter(position > 5000) %>% pull(position) -> positions

clip_df %>% group_by(position) %>%
  filter(position %in% c(5730:16023)) %>% # coordinates of heavy and light chain origin
  summarize(n = n(), sumSA = sum(SA), sumClip = sum(clip_count), pearsonID = pearsonID) %>% filter(n == 1) %>% arrange(desc(sumSA)) %>%
  group_by(pearsonID) -> top_clips


lapply(ps, function(pearsonID){
  fread(paste0("pearson-data/",pearsonID,"_del_find.SA.tsv")) %>%
    mutate(c1 = pmin(out1, out2), c2 = pmax(out1, out2)) %>%
    group_by(c1, c2) %>%
    summarize(count = n()) %>% mutate(donor = pearsonID)
}) %>% rbindlist() -> SA_df 

# Look for top interactions
p1 <- SA_df %>% group_by(donor) %>%
  filter(c1 %in% c(5730:16023) & c2 %in% c(5730:16023)) %>% # coordinates of heavy and light chain origin
  arrange(desc(count)) %>% mutate(rank = 1:n()) %>%
  top_n(8, wt = count) %>%
  ggplot(aes(x = as.character(rank), y = count)) +
  facet_wrap(~donor, scales = "free_y") +
  geom_point(size = 0.5) + pretty_plot(fontsize = 7) +
  L_border() + 
  labs(x = "Rank ordered deletion junctions", y = "count")
cowplot::ggsave2(p1, filename = "plots/Pearson_junctions.pdf", width = 3.5, height = 1.5)
