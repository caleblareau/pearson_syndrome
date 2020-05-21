library(Seurat)
library(Signac)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BuenColors)
library(motifmatchr)
library(chromVARmotifs)
data("human_pwms_v2") #filtered collection of human motifs from cisBP database
"%ni%" <- Negate("%in%")
source("../functions/import_counts_mat_10x.R")
set.seed(1)

# Import ATAC data
#counts <- import_counts_mat_10x("../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_filtered_peak_bc_matrix/")
pbmc <- readRDS("../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")

# Import mtDNA data
df_BCI <- fread("../data/deletion_heteroplasmy/del_PBMC_BCI.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del6073-13095") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-1", cell_id))

df_CCF <- fread("../data/deletion_heteroplasmy/del_PBMC_CCF.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del8482-13447") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-2", cell_id))

df_PT3 <- fread("../data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv") %>%
  dplyr::filter(deletion == "del10381-15407") %>% mutate(barcode = gsub(pattern = "-1", replacement = "-3", cell_id))

het_df <- rbind(df_BCI, df_CCF, df_PT3)
md <- pbmc@meta.data
bdf <- merge(md, het_df, by = "barcode") %>% group_by(patient) %>% mutate(scaled_heteroplasmy = scale(heteroplasmy, scale = FALSE))

# Format counts data for chromVAR
counts <- counts[,rownames(md)]
mdf <- data.frame(matrix(unlist(strsplit(rownames(counts), ":|-")), ncol = 3, byrow = TRUE))
mdf$chr <- as.character(mdf[,1])
mdf$start <- as.numeric(as.character(mdf[,2]))
mdf$end <- as.numeric(as.character(mdf[,3]))
gr <- makeGRangesFromDataFrame(mdf)

# Make a summarized experiment for CV
SE <- SummarizedExperiment(
  assays = list(counts = counts), 
  rowRanges = gr,
  colData = pbmc@meta.data
)

# Run chromVAR
SE <- filterPeaks(SE)
mm <- motifmatchr::matchMotifs(human_pwms_v2, SE, genome = BSgenome.Hsapiens.UCSC.hg19)
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
dev <- computeDeviations(object = SE,  annotations = mm)
devmat <- assays(dev)[["deviations"]][,bdf$barcode]
row.names(devmat) <- make.unique(unname(rowData(dev)$name))
saveRDS(devmat, "../../../pearson_mtscatac_large_data_files/output/cv_deviations_pearson_pts_only.rds")

# Look at the association between each TF and the 
donors <- c("BCI", "CCF", "PT3")
compute_dev_assoc <- function(test_celltype, perm = 0){
  lapply(donors, function(pt_id_1){
    cells_df <- bdf %>% dplyr::filter(predicted.id == test_celltype & Patient == pt_id_1) 
    ss_mat <- data.matrix(devmat[,cells_df$barcode])
    het_vec <- cells_df %>% pull(heteroplasmy)
    if(perm != 0){
      set.seed(perm)
      het_vec <- sample(het_vec)
    }
    
    lapply(143, function(idx){  # specifically test JUN
      data.frame(
        TF = rownames(ss_mat)[idx],
        Celltype = test_celltype,
        cor_pearson = cor(ss_mat[idx,], het_vec,use = "pairwise.complete.obs"),
        cor_spearman = cor(ss_mat[idx,], het_vec,use = "pairwise.complete.obs", method = "spearman"),
        pt_id_1
      ) -> df
      colnames(df) <- c("gene", "celltype","cor_pearson", "cor_spearman", "patient")
      if(perm == 0){
        df$pvalue <- cor.test(ss_mat[idx,], het_vec,use = "pairwise.complete.obs", method = "spearman")$p.value
      }
      df
    }) %>% rbindlist() %>% data.frame() -> odf
    
    odf
  }) %>% rbindlist() -> adf
  adf
}

# Compute per cell type
test_celltypes <- names(table(md$predicted.id))[table(md$predicted.id) > 70]
lapply(test_celltypes, function(test_celltype){
  compute_dev_assoc(test_celltype)
}) %>% rbindlist() %>% data.frame() -> all_assocs_chromvar

lapply(1:100, function(pp){
  lapply(test_celltypes, function(test_celltype){
    compute_dev_assoc(test_celltype, perm = pp)
  }) %>% rbindlist() %>% data.frame() -> dd
  dd
}) %>% rbindlist() -> all_assocs_chromvar_permuted

perm <- all_assocs_chromvar_permuted %>% dplyr::filter(gene %in% c("JUN")) %>% pull(cor_spearman) 
obs <- all_assocs_chromvar %>% dplyr::filter(gene %in% c("JUN")) %>% pull(cor_spearman)
perm_mean <- all_assocs_chromvar_permuted %>% group_by(celltype, patient) %>% summarize(ss = mean(cor_spearman)) %>% pull(ss)
ks.test(perm, obs)
wilcox.test(perm, obs)

all_assocs_chromvar$padj <- p.adjust(all_assocs_chromvar$pvalue)
qplot(perm_mean, obs)

# Make the  token example plots
bci_df <- bdf %>% dplyr::filter(predicted.id == "Effector CD8+ T-cells" & Patient == "BCI") 
bci_df$JUN <- devmat[143,bci_df$barcode]
p1 <- ggplot(bci_df, aes(x = heteroplasmy, y = JUN)) + 
  geom_point(alpha = 0.5, size = 0.2) +  labs(x = "% Heteroplasmy", y = "JUN deviation score") +
  geom_smooth(method='lm', formula= y~x) + pretty_plot(fontsize = 7) + L_border()

ccf_df <- bdf %>% dplyr::filter(predicted.id == "Memory CD4+ T-cells" & Patient == "CCF") 
ccf_df$JUN <- devmat[143,ccf_df$barcode]
p2 <- ggplot(ccf_df, aes(x = heteroplasmy, y = JUN)) + 
  geom_point(alpha = 0.5, size = 0.2) +  labs(x = "% Heteroplasmy", y = "JUN deviation score") +
  geom_smooth(method='lm', formula= y~x) + pretty_plot(fontsize = 7) + L_border()

pt3_df <- bdf %>% dplyr::filter(predicted.id == "CD14+ Monocytes" & Patient == "PT3") 
pt3_df$JUN <- devmat[143,pt3_df$barcode]
p3 <- ggplot(pt3_df, aes(x = heteroplasmy, y = JUN)) + 
  geom_point(alpha = 0.5, size = 0.2) + labs(x = "% Heteroplasmy", y = "JUN deviation score") +
  geom_smooth(method='lm', formula= y~x) + pretty_plot(fontsize = 7) + L_border()

viz_df <- data.frame(
  val = c(obs, perm),
  what = c(rep("obs", length(obs)), rep("perm", length(perm)))
)

p4 <- ggplot(viz_df, aes(x = val, color = what)) + 
  geom_density() + 
  pretty_plot(fontsize = 7) + L_border() + labs(x = "All correlations", y = "density") +
  theme(legend.position = "none") + scale_color_manual(values = c("dodgerblue2", "black"))

cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, p4, nrow = 1), file = "../output/scATAC_4panels.pdf", width = 4.4, height = 1.1)

all_assocs_chromvar$celltype <- gsub("Activated B-cells", "Naive B-cells", all_assocs_chromvar$celltype)
pB <- ggplot(all_assocs_chromvar, aes( x= patient, y = celltype, fill = cor_spearman)) +
  geom_tile() + labs(x = "", y = "", fill = "") + pretty_plot(fontsize = 7) + L_border() +
  scale_fill_gradientn(colours = c("dodgerblue4", "white", "red"),
                       limits = c(-0.36, 0.36)) + 
  scale_y_discrete(limits = rev(sort(unique(all_assocs_chromvar$celltype))), expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(legend.position = "bottom")

cowplot::ggsave2(pB, file = "../output/correlation_grid_jun_heteroplasmy.pdf", width = 1.7, height = 2.3)



