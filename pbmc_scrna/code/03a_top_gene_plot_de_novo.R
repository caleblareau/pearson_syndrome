library(Seurat)
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
options(future.globals.maxSize = 4000 * 1024^2)

md <- readRDS("../../../pearson_mtscatac_large_data_files/output/5March-PearsonRNAseq-integration.rds")@meta.data

# Import
import_scRNAseq <- function(dir_base, name, pheno, n){
  data.dir <- paste0("../data/", dir_base)
  raw <- Read10X(data.dir = data.dir); colnames(raw) <- paste0(name, "-", substr(colnames(raw),1,16), "_", n)
  
  # import scrublet results
  singlets <- fread(paste0("../data/scrublet_out/", dir_base, ".scrub.tsv")) %>%
    data.frame() %>% dplyr::filter(!called) %>% pull(barcode)
  
  # Filter for singlet genes and non-mitos
  raw <- raw[!grepl("^MT|^RP", rownames(raw)),paste0(name, "-", substr(singlets,1,16), "_", n)]
  raw
}

# Import

pBCI <- import_scRNAseq("pearson_bci", "pBCI", "Pearson", "1")
pCCF <- import_scRNAseq("pearson_ccf", "pCCF", "Pearson", "2")
pPT3 <- import_scRNAseq("pearson_mds", "pPT3", "Pearson", "3")

pbmc1 <- import_scRNAseq("healthy_pbmc_8k_v2-remap", "H1", "Healthy", "4")
pbmc2 <- import_scRNAseq("healthy_pbmc_4k_v2-remap", "H1", "Healthy", "5")
pbmc3 <- import_scRNAseq("healthy_pbmc_5k_nextgem", "H2", "Healthy", "6")
pbmc4 <- import_scRNAseq("healthy_pbmc_5k", "H2", "Healthy", "7")

counts <- cbind(pbmc1, pbmc2, pbmc3, pbmc4, pBCI, pCCF, pPT3)

donors <- c("H1", "H2", "pBCI", "pCCF", "pPT3")
cts <- unique(sort(as.character(df$celltype)))
boo_ct <- table((sort(as.character(df$celltype)))) > 250 & cts %ni% c( "low count 1", "CCF cell 1")
cts_go <- cts[boo_ct]

lapply(donors, function(d){
  lapply(cts_go, function(ct){
    bc <- rownames(md)[md$celltype == ct & md$source == d]
    vec <- rowSums(counts[,bc])
    cpm <- vec/sum(vec)*1000000
    data.frame(donor = d, celltype = ct, value = cpm["GPR85"])
  }) %>% rbindlist() 
})%>% rbindlist() %>% data.frame() -> gpr85

p_gpr85 <- ggplot(gpr85, aes(x = donor, y = celltype, fill = log2(value + 1))) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0,6)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(limits = rev(levels(as.factor(gpr85$celltype))), expand = c(0,0)) +
  labs(x = "", y = "") + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") 

cowplot::ggsave2(p_gpr85, file = "../plots/GPR85_expression.pdf", width = 1.7, height = 2.3)

