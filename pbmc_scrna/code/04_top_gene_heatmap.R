library(Seurat)
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
options(future.globals.maxSize = 4000 * 1024^2)

edgeR_celltypes <- sort(unique(readRDS("../../../pearson_large_data_files/output/pbmc/11July2022-PearsonRNAseq-diffGE-edgeR.rds")$celltype))
libs <- c("pearson_bci", "pearson_ccf", "pearson_mds",
          "healthy_pbmc_8k_v2-remap", "healthy_pbmc_4k_v2-remap", "healthy_pbmc_5k_nextgem", "healthy_pbmc_5k",
          "pediatrichealthy_pbmc_31687", "pediatrichealthy_pbmc_31697")
import_df <- function(dir_base, short_id){
  idf <- readRDS(paste0("../output/seurat_projected_meta_", dir_base, ".rds"))
  rownames(idf) <- paste0(short_id, "_", rownames(idf))
  idf
}

md <- rbind(
  import_df("pearson_bci", "pBCI"),
  import_df("pearson_ccf", "pCCF"),
  import_df("pearson_mds", "pPT3"),
  import_df("healthy_pbmc_8k_v2-remap", "H1a"),
  import_df("healthy_pbmc_4k_v2-remap", "H1b"),
  import_df("healthy_pbmc_5k_nextgem", "H2a"),
  import_df("healthy_pbmc_5k", "H2b"),
  import_df("pediatrichealthy_pbmc_31687", "P1"),
  import_df("pediatrichealthy_pbmc_31697", "P2")
)

# Import
import_scRNAseq <- function(dir_base, name, pheno, n){
  data.dir <- paste0("../data/", dir_base)
  raw <- Read10X(data.dir = data.dir); colnames(raw) <- paste0(name, "_", substr(colnames(raw),1,16), "-1")
  
  # import scrublet results
  singlets <- fread(paste0("../data/scrublet_out/", dir_base, ".scrub.tsv")) %>%
    data.frame() %>% dplyr::filter(!called) %>% pull(barcode)
  
  # Filter for singlet genes and non-mitos
  raw <- raw[!grepl("^MT|^RP", rownames(raw)),paste0(name, "_", substr(singlets,1,16), "-1")]
  raw
}

# Import

pBCI <- import_scRNAseq("pearson_bci", "pBCI", "Pearson", "1")
pCCF <- import_scRNAseq("pearson_ccf", "pCCF", "Pearson", "2")
pPT3 <- import_scRNAseq("pearson_mds", "pPT3", "Pearson", "3")

pbmc1 <- import_scRNAseq("healthy_pbmc_8k_v2-remap", "H1a", "Healthy", "4")
pbmc2 <- import_scRNAseq("healthy_pbmc_4k_v2-remap", "H1b", "Healthy", "5")
pbmc3 <- import_scRNAseq("healthy_pbmc_5k_nextgem", "H2a", "Healthy", "6")
pbmc4 <- import_scRNAseq("healthy_pbmc_5k", "H2b", "Healthy", "7")
pbmc5 <- import_scRNAseq("pediatrichealthy_pbmc_31687", "P1", "Healthy", "8")
pbmc6 <- import_scRNAseq("pediatrichealthy_pbmc_31697", "P2", "Healthy", "9")

counts <- cbind(pbmc1, pbmc2, pbmc3, pbmc4, pbmc5, pbmc6, pBCI, pCCF, pPT3)

donors <- c("H1", "H2", "P1", "P2", "pBCI", "pCCF", "pPT3")
cts_go <- cts <-edgeR_celltypes

md$celltype <- gsub(" ", ".", md$predicted.celltype.l2)
lapply(donors, function(d){
  lapply(cts_go, function(ct){
    print(d)
    print(ct)
    bc <- rownames(md)[md$celltype == ct & md$name == d]
    vec <- rowSums(counts[,bc])
    cpm <- vec/sum(vec)*1000000
    data.frame(donor = d, celltype = ct, GPR85 = cpm["GPR85"], 
               CXCL14 = cpm["CXCL14"], PLEKHD1 = cpm["PLEKHD1"],
               JUN = cpm["JUN"], FOS = cpm["FOS"],
               CD69 = cpm["CD69"],
               JAK1 = cpm["JAK1"],
               ANXA1 = cpm["ANXA1"],
               HES1 = cpm["HES1"],
               NFAM1 = cpm["NFAM1"],
               
               JAML = cpm["JAML"],
               STAT3 = cpm["STAT3"],
               RGS1 = cpm["RGS1"],
               SLC2A3 = cpm["SLC2A3"],
               ALDOA = cpm["ALDOA"],
               LINC01641 = cpm["LINC01641"], C12orf54 = cpm["C12orf54"])
  }) %>% rbindlist() 
})%>% rbindlist() %>% data.frame() -> ddf


p_a <- ggplot(ddf, aes(x = donor, y = celltype, fill = log2(NFAM1 + 1))) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0,6), oob = scales::squish) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(limits = rev(levels(as.factor(ddf$celltype))), expand = c(0,0)) +
  labs(x = "", y = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") 
p_a
cowplot::ggsave2(p_a, file = "../plots/NFAM1_expression.pdf", width = 1.7, height = 2.3)


p_cd69 <- ggplot(ddf, aes(x = donor, y = celltype, fill = log2(CD69 + 1))) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(4,11), oob = scales::squish) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(limits = rev(levels(as.factor(ddf$celltype))), expand = c(0,0)) +
  labs(x = "", y = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") 
p_cd69

cowplot::ggsave2(p_cd69, file = "../plots/CD69_expression.pdf", width = 1.7, height = 2.3)


p_gpr85 <- ggplot(ddf, aes(x = donor, y = celltype, fill = log2(GPR85 + 1))) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0,6), oob = scales::squish) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(limits = rev(levels(as.factor(ddf$celltype))), expand = c(0,0)) +
  labs(x = "", y = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") 
p_gpr85

cowplot::ggsave2(p_gpr85, file = "../plots/GPR85_expression.pdf", width = 1.7, height = 2.3)

p_PLEKHD1 <- ggplot(ddf, aes(x = donor, y = celltype, fill = log2(PLEKHD1 + 1))) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0,6), oob = scales::squish) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(limits = rev(levels(as.factor(ddf$celltype))), expand = c(0,0)) +
  labs(x = "", y = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") 
p_PLEKHD1

p_C12orf54 <- ggplot(ddf, aes(x = donor, y = celltype, fill = log2(C12orf54 + 1))) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0,6), oob = scales::squish) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(limits = rev(levels(as.factor(ddf$celltype))), expand = c(0,0)) +
  labs(x = "", y = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") 
p_C12orf54
cowplot::ggsave2(p_C12orf54, file = "../plots/C12orf54_expression.pdf", width = 1.7, height = 2.3)


p_LINC01641 <- ggplot(ddf, aes(x = donor, y = celltype, fill = log2(LINC01641 + 1))) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0,6), oob = scales::squish) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(limits = rev(levels(as.factor(ddf$celltype))), expand = c(0,0)) +
  labs(x = "", y = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") 
p_LINC01641

cowplot::ggsave2(cowplot::plot_grid(p_C12orf54, p_LINC01641, p_PLEKHD1, nrow = 1), file = "../plots/supplment_expression.pdf", width = 5.1, height = 2.3)

p_cxcl14 <- ggplot(ddf, aes(x = donor, y = celltype, fill = log2(CXCL14 + 1))) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0,6), oob = scales::squish) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(limits = rev(levels(as.factor(ddf$celltype))), expand = c(0,0)) +
  labs(x = "", y = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") 
p_cxcl14

cowplot::ggsave2(p_cxcl14, file = "../plots/CXCL14_expression.pdf", width = 1.7, height = 2.3)


pJun <- ggplot(ddf, aes(x = donor, y = celltype, fill = log2(JUN +1))) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(6,14), oob = scales::squish) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(limits = rev(levels(as.factor(ddf$celltype))), expand = c(0,0)) +
  labs(x = "", y = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") 
pJun

pFos <- ggplot(ddf, aes(x = donor, y = celltype, fill = log2(FOS +1))) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(6,14), oob = scales::squish) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(limits = rev(levels(as.factor(ddf$celltype))), expand = c(0,0)) +
  labs(x = "", y = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") 
pFos

cowplot::ggsave2(cowplot::plot_grid(pJun, pFos, nrow = 1),
                 file = "../plots/JUNFOS_expression.pdf", width = 3.4, height = 2.3)
