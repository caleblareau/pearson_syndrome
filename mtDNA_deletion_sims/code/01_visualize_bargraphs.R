library(gridExtra)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(BuenColors)

show_reason_window_plot <- function(bp1, bp2, stats_file_path){
  
  read_get_range <- function(file, bp1, bp2){
    read_tsv(file, col_names = c("start", "end", "lc", "rc", "clip_pos", "read_name", 'X6', 'X7')) %>%
      dplyr::select(-X6,-X7) %>%
      filter(start %in% (bp1-71):bp1 | end %in% bp2:(bp2+71))
  }
  
  files <- stats_file_path
  l <- lapply(files, function(x) read_get_range (x, bp1, bp2))
  names(l) <- gsub(".st.tsv", "", files)
  
  df <- bind_rows(l, .id = "id") %>%
    filter(clip_pos %in% c(0,bp1,bp2)) %>%
    mutate(clipped = factor(ifelse(clip_pos == 0, "unclipped", "clipped")),
           id = gsub(paste0("del_",bp1,"_",bp2,"."), "", id)) %>%
    separate(id, into = c("depth","h1", "h2"), sep = "\\.") %>%
    mutate(depth = as.numeric(depth),
           h1 = factor(h1),
           h2 = factor(h2)) %>%
    arrange(h2)
  
  df_count_start <- df %>%
    filter(start < bp1) %>%
    group_by(start, clipped,depth) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    mutate(normalized_count = count/depth) %>%
    mutate(start = start-(start[1]+71)) 
  
  df_count_end <- df %>%
    filter(end > bp2) %>%
    group_by(end, clipped,depth) %>%
    summarise(count  = n()) %>%
    ungroup() %>%
    mutate(normalized_count = count/depth) %>%
    mutate(end = end-(end[1]-1))
  
  
  p1=ggplot(df_count_start, aes(x = start-(df_count_start$start[1]+71), y = normalized_count,fill = clipped)) +
    geom_bar(stat = 'identity',position='fill', width= 2) +
    scale_fill_manual(values=c("#C1272D", "#2E3192")) + theme_minimal() +
    scale_y_continuous(labels = scales::percent_format()) + xlab('start') + ylab('percent') +
    scale_x_continuous(breaks = seq(-71,-1, by = 10)) 
  
  
  p2=ggplot(df_count_end, aes(x = end-(df_count_end$end[1]-1), y = normalized_count, fill = clipped)) +
    geom_bar(stat='identity',position='fill', width=2) +
    scale_fill_manual(values=c("#C1272D", "#2E3192")) + theme_minimal() +  
    scale_y_continuous(labels = scales::percent_format()) + xlab('end') + ylab('percent') + 
    scale_x_continuous(breaks = seq(1,71, by = 10)) 
  
  p12 =grid.arrange(p1,p2)
  ggsave(p12, file=paste0('8483_13445.pdf'), width=20,height=10)
  
  
}