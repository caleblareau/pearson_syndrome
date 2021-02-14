library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(matrixStats)
library(ggrepel)
library(readxl)
library(gridExtra)
library(BuenColors)

# Author: Kopal Garg

# read in clip .tsv files and convert to a matrix
brain_files <- list.files("../data/brain_dels_SA/", pattern='*.tsv')
brain_fnames <- gsub(".tsv", "", brain_files )

brain_l <- lapply(paste0( brain_files), function(x){
  fread(paste0("../data/brain_dels_SA/",x))
})
names(brain_l) <- brain_fnames
brain <- brain_l %>% 
  bind_rows(., .id="ID")

# idxstats per chrM
reads=as.data.frame(read_excel('../data/meta/total_chrM_dedup_reads.xlsx', sheet='Brain', col_names = F)) %>% setNames(c('ID', 'reads'))
reads$ID=gsub(reads$ID, pattern='.qc.bam', replacement='')
brain = inner_join(brain, reads)

# Blacklist the recurrent ones
blacklist <- 
  brain %>% 
  group_by(position) %>%
  mutate(ratio = clip_count/reads) %>% # ratio = clip count / reads 
  mutate(z = scale(ratio)) %>% # z-scores of the ratio
  ungroup() %>% 
  group_by(ID) %>%
  mutate(pct=ntile(ratio,100)) %>% # as a percentile
  ungroup() %>% 
  group_by(position) %>%
  tally(pct) %>%
  filter(n>=length(brain_files)*95) %>% pull(position) # if present in > 95% of samples, add to blacklist (arbitrary)

brain_df <- brain %>% 
  group_by(ID) %>%
  mutate(ratio = clip_count/reads) %>% 
  mutate(z=scale(ratio)) %>% ungroup() %>%  
  mutate(blacklist = factor(ifelse(position %in% blacklist, TRUE,FALSE))) 

SA_file='../data/full_output/SRR7700727.SA.tsv'
id='SRR7700727'
SA= fread(SA_file) 

SA <- SA %>% 
  mutate(out3 = case_when(out2 > out1 ~ out1, out1 >= out2 ~ out2)) %>%
  mutate(out4 = case_when(out1 > out2 ~ out1, out2 >= out1 ~ out2)) %>%
  subset(., select=-c(out1,out2)) %>% setNames(c('main', 'SA')) %>% arrange(-main) %>% 
  group_by(SA, main) %>% add_tally() %>% 
  filter(n>1) %>% unique()

SA <- SA %>% filter(!(SA %in% blacklist) & !(main %in% blacklist)) %>%
  setNames(c('position', 'supp_A', 'weight'))

df <- brain_df %>% filter(ID==id) %>%
  mutate(z = case_when(position %in% c(1,16569) ~ 0, TRUE ~ z)) %>%
  arrange(-z) 

# summarise per bin for segment median plot on top of coverage
df2 = brain_df %>%filter(ID==id) %>% filter(blacklist==FALSE) %>% 
  mutate(pct = ntile(position, as.integer(16569/5))) %>% 
  group_by(pct) %>% 
  mutate(pct_median = median(coverage)) %>%
  mutate(df_median = median(df$coverage))

# remove the first few and last few positions 
df = df %>% filter(position %in% c(500:16000)) %>% 
  mutate(z_SA = scale(SA)) %>% 
  mutate(z_clip_count= scale(clip_count)) %>% 
  mutate(z_coverage= scale(coverage))

df_SA <- left_join(df,SA)

p1 <- df %>% ggplot(., aes(x = position, y = z_SA)) +
  geom_col(aes(fill=blacklist),width = 50,
           show.legend = FALSE) +
  scale_fill_manual(values=c("black", "lightgrey")) +
  geom_label_repel(
    data = (df %>% filter(position %in% c(500:16000)) %>% filter(blacklist!="TRUE") %>% filter(SA>0) %>% arrange(-z) )[1:6,],
    aes(label = position),
    size = 2, show.legend = FALSE) +
  pretty_plot(fontsize =  8) + L_border()+ 
  xlab("") +
  ylab("Deletion enrichment") + 
  ggtitle(id) + #y=z_SA
  geom_curve(aes(x = position, y = 0, xend = supp_A, yend = 0, alpha=weight), show.legend = FALSE,curvature=-0.5,data = df_SA %>% 
               filter(!is.na(weight)) %>% filter(weight>50))

p2 <- df2 %>% filter(position %in% c(500:16000)) %>% ggplot(., aes(x = position, y = pct_median)) + geom_line() +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "Position on mtDNA genome", y = "Coverage")

gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)
maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
grid.arrange(gp1,gp2)





