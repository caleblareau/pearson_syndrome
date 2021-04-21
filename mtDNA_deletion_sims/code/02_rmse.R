library(Metrics)
library(viridisLite)
library(readr)
library(scales)
library(dplyr)

file_path <- "../data_from_terra//summary_all_72.tsv.gz"
merged<-as.data.frame(read_tsv(file_path))[c(1:5,8)]
names(merged) <- c('cell_id', 'heteroplasmy', 'reads_del', 'reads_wt', 'reads_all', "del")

fig1 <- merged %>% dplyr::filter(del %in% c("13157_15477","9232_13413","8482_13445"))
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
  summarise(rmse_del_cond = rmse(heteroplasmy, true_heteroplasmy)) %>%
  mutate(rmse_scaled = scale(rmse_del_cond))
rmse_info <- rmse_info[order(as.numeric(rmse_info$far)),]
rmse_info %>% group_by(del) %>% arrange(rmse_del_cond) %>% slice_min(order_by = rmse_del_cond)
rmse_info %>% dplyr::filter(far == 9 & near == 24)

p1 <- ggplot(rmse_info, aes(x = (as.numeric(far)), y =(as.numeric(near)), fill = rmse_del_cond)) +
  geom_tile() +
  facet_wrap(~del) + 
  scale_fill_gradientn(colors = jdb_palette("calma_morado"),na.value = "white",  limits = (c(0,25)))+ xlab('Far') + ylab('Near') +
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "bottom") +labs(fill = "RMSE") +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

cowplot::ggsave2(p1, file = "../output/RMSE_pearson_del_Fig1.pdf", width = 5.7, height = 2.7)


file_path <- "../data_from_terra//summary_all_72.tsv.gz"
merged<-as.data.frame(read_tsv(file_path))[c(1:5,8)]
names(merged) <- c('cell_id', 'heteroplasmy', 'reads_del', 'reads_wt', 'reads_all', "del")

fig2 <- merged %>% dplyr::filter(del %in% c("6073_13095", "8482_13445", "10381_15407"))
merged <-fig2 %>%
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
  summarise(rmse_del_cond = rmse(heteroplasmy, true_heteroplasmy)) %>%
  mutate(rmse_scaled = scale(rmse_del_cond))
rmse_info <- rmse_info[order(as.numeric(rmse_info$far)),]
rmse_info %>% group_by(del) %>% arrange(rmse_del_cond) %>% slice_min(order_by = rmse_del_cond)
rmse_info %>% dplyr::filter(far == 9 & near == 24)

p1 <- ggplot(rmse_info, aes(x = (as.numeric(far)), y =(as.numeric(near)), fill = rmse_del_cond)) +
  geom_tile() +
  facet_wrap(~del) + 
  scale_fill_gradientn(colors = jdb_palette("calma_morado"),na.value = "white",  limits = (c(0,25)))+ xlab('Far') + ylab('Near') +
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "bottom") +labs(fill = "RMSE") +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

cowplot::ggsave2(p1, file = "../output/RMSE_pearson_del_primaryCells.pdf", width = 5.7, height = 2.7)



