library(edgeR)
library(Seurat)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(stringr)

"%ni%" <- Negate("%in%")

# Function to do differential gene expression with edgeR w/ covariates
run_edgeRQLFDetRate_CL <- function(count, condt) {
  
  donor_id <- str_split_fixed(colnames(count), "_", 2)[,1]
  dge <- DGEList(count, group = condt)
  dge <- calcNormFactors(dge)
  
  # adjust for sequencing technology (H2 is 10x v3; Hped1 and Hped2 are muliome as well), as well as total genes detected, 
  cdr <- scale(colMeans(count > 0))
  design <- model.matrix(~ cdr + as.numeric(donor_id %in% c("H2a", "H2b")) + as.numeric(donor_id %in% c("Hped1", "Hped2")) + condt) 
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit)
  tt <- topTags(qlf, n = Inf)
  df <- signif(tt@.Data[[1]], 3)
  df$gene <- rownames(df)
  
  # small function to pull CPM
  ps <- function(which_ones){
    rs <- rowSums(count[,which_ones])
    cpms <- round(rs/sum(rs) *1000000,1)[as.character(df$gene)]
    return(cpms)
  } 
  
  # Pull CPM values
  df$Pearson_cpm <- ps(condt == "p")
  df$Healthy_cpm <- ps(condt == "H")
  
  # Pull donors
  df$H1v2 <- ps(donor_id %in% c("H1a", "H1b"))
  df$H2v3 <- ps(donor_id %in% c("H2a", "H2b"))
  df$Hped1 <- ps(donor_id %in% c("Hped1"))
  df$Hped2 <- ps(donor_id %in% c("Hped2"))
  df$pBCI <- ps(donor_id == "pBCI")
  df$pCCF <- ps(donor_id == "pCCF")
  df$pPT3 <- ps(donor_id == "pPT3")
  
  # Round
  df$logFC <- round(df$logFC,2)
  df$logCPM <- round(df$logCPM,2)
  df$signedTstat <- sign(df$logFC) * sqrt(round(df$F,2))
  df
  
}

#libs <- c("pearson_bci", "pearson_ccf", "pearson_mds", "healthy_pbmc_8k_v2-remap", "healthy_pbmc_4k_v2-remap", "healthy_pbmc_5k_nextgem", "healthy_pbmc_5k")
import_df <- function(dir_base, short_id){
  idf <- readRDS(paste0("../output/seurat_projected_meta_", dir_base, ".rds"))
  rownames(idf) <- paste0(short_id, "_", rownames(idf))
  idf
}

df <- rbind(
  import_df("pearson_bci", "pBCI"),
  import_df("pearson_ccf", "pCCF"),
  import_df("pearson_mds", "pPT3"),
  import_df("healthy_pbmc_8k_v2-remap", "H1a"),
  import_df("healthy_pbmc_4k_v2-remap", "H1b"),
  import_df("healthy_pbmc_5k_nextgem", "H2a"),
  import_df("healthy_pbmc_5k", "H2b"),
  import_df("pediatrichealthy_pbmc_31687", "Hped1"),
  import_df("pediatrichealthy_pbmc_31697", "Hped2")
  
)
df$celltype <- gsub(" ", ".", df$predicted.celltype.l2)
ctdf <- df

# Now import counts
import_counts <- function(dir_base, short_id){
  data.dir <- paste0("../data/", dir_base)
  gene_coords <- fread("../data/forinfercnv/infer_CNV_gene_annotations.tsv")
  xy_genes <- gene_coords %>% filter(V2 %in% c("X", "Y")) %>% pull(V1)
  raw <- Read10X(data.dir = data.dir)
  colnames(raw) <- paste0(short_id, "_", colnames(raw))
  
  raw[!(rownames(raw) %in% xy_genes),]
}

counts <- cbind(
  import_counts("pearson_bci", "pBCI"),
  import_counts("pearson_ccf", "pCCF"),
  import_counts("pearson_mds", "pPT3"),
  import_counts("healthy_pbmc_8k_v2-remap", "H1a"),
  import_counts("healthy_pbmc_4k_v2-remap", "H1b"),
  import_counts("healthy_pbmc_5k_nextgem", "H2a"),
  import_counts("healthy_pbmc_5k", "H2b"),
  import_counts("pediatrichealthy_pbmc_31687", "Hped1"),
  import_counts("pediatrichealthy_pbmc_31697", "Hped2")
)

counts <- counts[!grepl("^RP|^MT-", rownames(counts)),as.character(rownames(df))]
dim(counts)
dim(df)
df$barcode <- rownames(df)

cts <- unique(sort(as.character(df$celltype)))
boo_ct <- table((sort(as.character(df$celltype)))) >= 400 
cts_go <- cts[boo_ct]
length(cts_go)

library(ggbeeswarm)
p1 <- ctdf %>% group_by(celltype, name) %>% summarize(count = n())  %>% arrange((celltype)) %>% data.frame() %>%
  ungroup() %>% group_by(name) %>% mutate(prop = round(count / sum(count) *100, 2)) %>% data.frame() %>%
  dplyr::filter(celltype %in% cts_go) %>% 
  mutate(colorme = case_when(name %in% c("H1", "H2") ~ "HealthyAdult", 
                             name %in% c("P1", "P2") ~ "HealthyChild",
                             TRUE ~ "Pearson")) %>% 
  ggplot(aes(x = celltype, y = prop, color = colorme)) +
  geom_quasirandom(width = 0.3) +
  scale_color_manual(values = c("dodgerblue3", "green4","black")) + 
  pretty_plot(fontsize = 8)+ L_border() + theme(legend.position = "none") + 
  labs(x = "Azimuth labels", y = "% of cells") + 
  geom_vline(xintercept = 1:17 - 0.5, linetype = 2)
cowplot::ggsave2(p1, file = "../plots/azimuth_proportions.pdf", width = 6.5, height = 2)

list_of_de_mats <- lapply(cts_go, function(ct){
  print(ct)
  
  # Pull out celltype specific barcodes
  bcs <- df %>% filter(celltype == ct) %>% pull(barcode)
  counts_ct <- counts[,bcs]
  
  # Pull out barcodes
  boo_gene <- log10(rowSums(counts_ct)/sum(counts_ct)*1000000) > 0.5
  counts_go <- counts_ct[boo_gene,]
  cond <- substr(colnames(counts_go), 1,1)
  
  dfo <- run_edgeRQLFDetRate_CL(counts_go, cond)
  dfo$celltype <- ct
  dfo
})

melted_list <- rbindlist(list_of_de_mats) %>% arrange(FDR)
melted_list %>% arrange(desc(abs(logFC)))
saveRDS(melted_list, file = "../../../pearson_large_data_files/output/pbmc/11July2022-PearsonRNAseq-diffGE-edgeR.rds")


