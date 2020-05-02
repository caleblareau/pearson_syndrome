library(data.table)
library(dplyr)

estimate_chr7_cn <- function(depth_df, frag_file, out_file_pre){
  ff <- fread(frag_file)
  bcs <-depth_df %>% filter(V2 > 20) %>% pull(V1)
  ff %>% filter(V4 %in% bcs) %>%
    mutate(del_region = (V1 == "chr7") & (V2 > 110000000)) %>%
    group_by(V4, V1, del_region) %>% summarize(count = n()) -> count_df
  
  total_df <- count_df %>% group_by(V4) %>% summarize(all = sum(count))
  chr7_total_df <- count_df %>% filter(V1 == "chr7")  %>% group_by(V4) %>% summarize(all = sum(count))
  chr7_total_del_df <- count_df %>% filter(V1 == "chr7")  %>% filter(del_region) 
  
  # Make sure that all barcodes are represented
  total_df <- total_df %>% filter(V4 %in% chr7_total_del_df$V4)
  chr7_total_df <- chr7_total_df %>% filter(V4 %in% chr7_total_del_df$V4)
  
  data.frame(
    total_df, 
    chr7_total = chr7_total_df[[2]],
    chr7_total_del = chr7_total_del_df[[4]]
  ) %>% mutate(pct_reads_in_del = chr7_total_del/chr7_total*100) -> summary_df
  write.table(summary_df, file = paste0("../output/", out_file_pre, "_chr7del.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  out_file_pre
}
estimate_chr7_cn(
  fread("../../cd34_scatac/data/mgatk_depth/Pearson-CD34-PT3_v12-mtMask_mgatk.depthTable.txt", header = FALSE),
  "../../../pearson_mtscatac_large_data_files/input/CD34/Pearson-CD34-PT3_v12-mtMask_fragments.tsv.gz",
  "PT3_CD34"
)

estimate_chr7_cn(
  fread("../../bmmnc_scatac/data/Pearson-BMMNC-PT3_v12-mtMask_mgatk.depthTable.txt",header = FALSE),
  "../../../pearson_mtscatac_large_data_files/input/bmmnc/Pearson-BMMNC-PT3_v12_fragments.tsv.gz",
  "PT3_BMMNC"
)

estimate_chr7_cn(
  fread("../../pbmc_scatac/data/depth/Pearson-PBMC-PT3_v12-mtMask_mgatk.depthTable.txt", header = FALSE) %>%
    mutate(V1 = gsub(pattern = "-1", replacement = "-3", V1)),
  "../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_fragments.tsv.gz",
  "PT3_PBMC"
)
