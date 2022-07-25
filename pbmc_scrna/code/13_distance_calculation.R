library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(RANN)
library(data.table)
library(dplyr)
library(BuenColors)
library(stringr)
set.seed(1234)

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

# Create seurat graph and run PCA
so <- CreateSeuratObject(counts, meta.data = df)
so <- so %>% NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData() %>%
  RunPCA()
so$cellid <- rownames(so@meta.data)

# Get celltypes where we have good per-donor coverage
so@meta.data %>% group_by(celltype, name) %>%
  summarize(count = n()) %>%
  mutate(cb =count > 100) %>%
  ungroup() %>% group_by(celltype) %>% summarize(hit = sum(cb)) %>%
  filter(hit == 7)

subset_celltypes_to_test <- c("CD14.Mono", "CD16.Mono", "CD4.Naive", "CD4.TCM", "CD8.Naive", "NK")

estimate_distances_nn <- function(so, n_cells_test = 100, n_sims = 10){
  
  lapply(1:n_sims, function(sim_n){
    print(sim_n)
    lapply(subset_celltypes_to_test, function(celltype_x){ 
      
      # Pull data
      set.seed(sim_n)
      subset_df <- so@meta.data %>%
        mutate(mod_name =case_when(
          name %in% c("H1a", "H1b") ~ "H1",
          name %in% c("H2a", "H2b") ~ "H2",
          TRUE ~ name))%>%
        filter(celltype == celltype_x) %>%
        group_by(mod_name) %>% slice_sample(n = n_cells_test )
      cells <- subset_df %>% pull(cellid)
      
      # Compute pairwise distance for subset of cells
      d <- as.matrix(dist(so@reductions$pca@cell.embeddings[cells,1:30]) )
      colnames(d) <- make.unique(subset_df$mod_name)
      rownames(d) <- make.unique(subset_df$mod_name)
      dim(d)
      xy <- t(combn(colnames(d), 2))
      df_all_dist <- data.frame(xy, dist=d[xy]) 
      l1 <- str_split_fixed(df_all_dist$X1, "[.]", 2)[,1]
      l2 <- str_split_fixed(df_all_dist$X2, "[.]", 2)[,1]
      
      # annotate each distance
      pearson_ids <- c("pBCI", "pCCF", "pPT3")
      healthyadult_ids <- c("H1", "H2")
      pediatric_ids <- c("P1", "P2")
      same <- l1 == l2
      
      annotation <- case_when(
        same & (l1 %in% pearson_ids) ~ "aPearson_SameDonor",
        same & (l1 %in% healthyadult_ids) ~ "eAdult_SameDonor",
        same & (l1 %in% pediatric_ids) ~ "cPediatric_SameDonor",
        (l1 %in% pearson_ids) & (l2 %in% pearson_ids) ~ "bPearson_otherDonor",
        (l1 %in% healthyadult_ids) & (l2 %in% healthyadult_ids) ~ "fAdult_otherDonor",
        (l1 %in% pediatric_ids) & (l2 %in% pediatric_ids) ~ "dPediatric_otherDonor",
        
        ((l1 %in% pearson_ids) | (l2 %in% pearson_ids)) &  ((l1 %in% healthyadult_ids) | (l2 %in% healthyadult_ids)) ~ "zPearsonHealthy_between",
        ((l1 %in% healthyadult_ids) | (l2 %in% healthyadult_ids)) & ((l1 %in% pediatric_ids) | (l2 %in% pediatric_ids)) ~ "AdultPediatric_between",
        ((l1 %in% pediatric_ids) | (l2 %in% pediatric_ids)) & ((l1 %in% pearson_ids) | (l2 %in% pearson_ids)) ~ "PediatricPearson_between",
        
        TRUE ~ "hm"
      )
      
      out_df  <- data.frame(
        df_all_dist,
        annotation,
        sim_n = sim_n,
        celltype_x
      ) %>%
        filter(!(annotation %in% c("PediatricPearson_between", "AdultPediatric_between"))) # between distance is confounded by the technology it seems
      
    }) %>% rbindlist() %>% data.frame() -> all_distances_df_sim1
    all_distances_df_sim1
  }) %>% rbindlist() -> all_distances_df
  
  all_distances_df
}

distance_df <- estimate_distances_nn(so)

# visualize dinstances
p1 <- distance_df %>%
  ggplot(aes(x = celltype_x, y = dist, color = annotation)) +
  geom_boxplot(outlier.shape = NA) + 
  labs(x = "", y = "PCA distance", color = "") +
  pretty_plot() + L_border() +
  coord_cartesian(ylim = c(0, 40)) +
  scale_color_manual(values = c(jdb_palette("brewer_celsius")[c(2,1,6,7,8,9)], "black"))
cowplot::ggsave2(p1, file = "../plots/distances_scRNA.pdf", width = 7, height = 2.5)

median_df <- distance_df %>%
  group_by(annotation, celltype_x) %>% summarize(md = median(dist)) %>% data.frame()

# Pearson
round(median_df[7:12,3]/median_df[1:6,3],1)

# Healthy
round(median_df[19:24,3]/median_df[13:18,3],1)

# Pediatric
round(median_df[31:36,3]/median_df[25:30,3],1)

