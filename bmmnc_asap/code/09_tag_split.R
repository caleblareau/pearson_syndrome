library(dplyr)
protein <- readRDS("../data/Pearson_ASAP_allProteinCounts.rds")

lapply(1:5, function(i){
  ss <- data.matrix(protein[,grepl(as.character(i), colnames(protein))])
  colnames(ss) <- gsub(as.character(i), "1", colnames(ss))
  write.table(ss, row.names = TRUE, col.names = TRUE, sep = ",", quote = FALSE,
              file = paste0("../output/for_geo/Pearson-ASAP-BMMNC-tag-c", i, "-counts.csv"))
})