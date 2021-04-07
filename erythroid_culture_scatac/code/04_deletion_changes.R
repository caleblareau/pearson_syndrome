library(data.table)
library(dplyr)
library(BuenColors)

process_one_library <- function(lib, day){
  cells <- fread(paste0("../output/cell_assignments_per_channel/Pearson_",lib,"_assign.tsv")) %>%
    dplyr::filter(assign == "Pearson")
  dels <- fread(paste0("../data/dels/del_Pearson_Healthy_",lib,".deletion_heteroplasmy.tsv")) %>%
    dplyr::filter(version == "improved" & deletion == "del10381-15407")
  mdf <- merge(cells, dels, by.x = "barcode", by.y = "cell_id")
  mdf$library <- lib
  mdf$day <- day
  
  mdf
}

all_df <- rbind(
  process_one_library("D6_1", "Day06"),
  process_one_library("D6_2", "Day06"),
  process_one_library("D12_1", "Day12"),
  process_one_library("D12_2", "Day12")
)

ggplot(all_df, aes(x = day, y = heteroplasmy, color = library)) + 
  geom_violin()
