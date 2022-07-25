library(Metrics)
library(viridisLite)
library(readr)
library(scales)
library(dplyr)
library(data.table)
library(BuenColors)

process_fp <- function(rl,extra){
  file_path <- paste0("../data_from_terra/df_",rl,".tsv.gz")
  merged<-as.data.frame(read_tsv(file_path))[c(1:5,8)]
  names(merged) <- c('cell_id', 'heteroplasmy', 'reads_del', 'reads_wt', 'reads_all', "del")
  
  fig1 <- merged
  merged <-fig1 %>%
    mutate(
      b_l = matrix(unlist(strsplit(cell_id,'_')), ncol=5,byrow=T)[,4],
      info = matrix(unlist(strsplit(cell_id,'_')), ncol=5,byrow=T)[,5],
      depth =  matrix(unlist(strsplit(info, '\\.')), ncol=5, byrow=T)[,2],
      b_r = matrix(unlist(strsplit(info,'\\.')), ncol=5,byrow=T)[,1],
      far=matrix(unlist(strsplit(cell_id,'_')), ncol=5,byrow=T)[,1],
      near_info=matrix(unlist(strsplit(cell_id,'_')), ncol=5,byrow=T)[,2],
      near=matrix(unlist(strsplit(near_info,'\\.')), ncol=2,byrow=T)[,1],
      del=paste0(b_l,'_',b_r),
      true_heteroplasmy=matrix(unlist(strsplit(info,'\\.')),ncol=5,byrow=T)[,4]) %>%
    subset(., select=-c(reads_del, reads_wt, reads_all,near_info,info)) %>%
    mutate(true_heteroplasmy = as.numeric(as.character(true_heteroplasmy)))
  
  rmse_info <- merged %>%
    group_by(del, far, near) %>%
    summarise(rmse_del_cond = rmse(heteroplasmy, true_heteroplasmy), meand = mean(depth)) %>%
    mutate(rmse_diff = rmse_del_cond - min(rmse_del_cond),
           near_far = paste0(far, "_", near))
  rmse_info <- rmse_info[order(as.numeric(rmse_info$far)),]
  rmse_info$rl <- paste0(extra, rl)
  rmse_info
}
mdf <- rbind(process_fp("50","a"),process_fp("72","b"),process_fp("100","c"))

get_rmse_diff <- function(del1){
  mdf %>%
    filter(rl == "b72") %>%
    filter(del==del1) %>%
    group_by(del) %>%
    slice_min(n = 1, order_by = rmse_del_cond) %>%
    pull(near_far) -> nf
  
  mdf %>%
    filter(del==del1) %>%
    filter(near_far %in% nf) %>% data.frame()
}

lapply(c("8427_15639", "9232_13413","8482_13445"), get_rmse_diff) %>% 
  rbindlist() -> difference_of_optimum_df


get_rmse_default <- function(del1){

  mdf %>%
    filter(del==del1) %>%
    filter(near_far %in% paste0("9_24")) %>% data.frame()
}

lapply(c("8427_15639", "9232_13413","8482_13445"), get_rmse_default) %>% 
  rbindlist()-> difference_of_default_df

px <- ggplot(difference_of_default_df, aes(x = del, y = rmse_diff, fill = as.character(rl))) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.5) +
  coord_flip() + labs(y = "Difference in RMSE (% Heteroplasmy)") +
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("darkgrey", "grey", "white")) +
  theme(legend.position = "none")
cowplot::ggsave2(px, file = "../output/RMSE_diff_default.pdf", width = 2, height = 1.6)
 