library(data.table)
library(Matrix)
library(dplyr)

#  Import raw things
permuted_patients_mat <- readRDS("../data/PERMUTED_patients_spearman_mat.rds")
all_patients_mat <- readRDS("../data/all_patients_spearman_mat.rds")

get_sig_df<- function(ptid, idx){
  cutoff <- as.numeric(quantile((sort(permuted_patients_mat[[ptid]][,idx])), c(0.01,0.99)))
  
  data.frame(
    Patient = names(all_patients_mat)[ptid],
    Celltype = colnames(all_patients_mat[[ptid]])[idx],
    Correlation = all_patients_mat[[ptid]][,idx],
    Gene = make.unique(rownames(all_patients_mat[[ptid]]))
  ) %>% filter(Correlation < cutoff[1] | Correlation > cutoff[2])
}

lapply(1:3, function(ptid){
  lapply(1:13, function(idx){
    get_sig_df(ptid, idx)
  }) %>% rbindlist() %>% data.frame() 
})%>% rbindlist() %>% data.frame() %>% mutate(direction = Correlation > 0) -> everything_FDR01
table(everything_FDR01$direction)

sort(table(paste0(everything_FDR01$Gene, "_",everything_FDR01$direction)), decreasing = TRUE) %>% head(10)

everything_FDR01 %>% filter(Gene == "GAS7")
