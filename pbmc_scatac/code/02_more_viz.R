library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg19)
library(chromVAR)
library(JASPAR2018)
library(motifmatchr)
library(TFBSTools)

pbmc2 <- readRDS("../../../pearson_large_data_files/output/PBMC_scATAC_3P1H-15FEB2021.rds")

# Filter to interesting genes
interesting_features <- c('NCAM1', 'CD4', 'CD8A', 'MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1',
                          'LYZ', 'CD14', 'CCR5', 'CXCR6', "RORA", "EOMES")
FeaturePlot(pbmc2, features = interesting_features, max.cutoff = "q95")

DefaultAssay(pbmc2) <- "RNA"
FindMarkers(pbmc2, ident.1 = "6")
FindMarkers(pbmc2, ident.1 = "4", ident.2 = "2")

#-------------------

DefaultAssay(pbmc2) <- 'peaks'

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = FALSE)
)
motif.matrix <- CreateMotifMatrix(
  features = pbmc2@assays$peaks@ranges,
  pwm = pwm,
  genome = 'hg19',
  use.counts = FALSE
)

pbmc2 <- AddMotifs(
  object = pbmc2,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  pfm = pwm
)
pbmc2 <- RunChromVAR(
  object = pbmc2,
  genome = BSgenome.Hsapiens.UCSC.hg19
)
DefaultAssay(pbmc2) <- 'chromvar'
Idents(pbmc2) <- "seurat_clusters"
topdiff <- rownames(FindMarkers(pbmc2, 6) %>% head(20))
sapply(pwm, name)[c(topdiff)]

pearson_obj <- pbmc2[,pbmc2@meta.data$Patient != "Healthy"]
p_motif <- FeaturePlot(pearson_obj, features = c("MA0800.1"), max.cutoff = "q98") &
  scale_color_gradientn(colors = jdb_palette("solar_extra")) &
  theme_void() & theme(legend.position = "none")
cowplot::ggsave2(p_motif, file = "../plots/EOMES_motif_deviation_onlyPearson.png", width = 5, height = 5.1, dpi = 500)

p_motif2 <-  FeaturePlot(pearson_obj, features = c("MA0768.1")) &
  scale_color_gradientn(colors = jdb_palette("solar_extra")) &
  theme_void() & theme(legend.position = "none")
cowplot::ggsave2(p_motif2, file = "../plots/LEF1_motif_deviation_onlyPearson.png", width = 5, height = 5.1, dpi = 500)

p_motif3 <-  FeaturePlot(pearson_obj, features = c("MA0072.1"), max.cutoff = "q98") &
  scale_color_gradientn(colors = jdb_palette("solar_extra")) &
  theme_void() & theme(legend.position = "none")
cowplot::ggsave2(p_motif3, file = "../plots/Rora_motif_deviation_onlyPearson.png", width = 5, height = 5.1, dpi = 500)

set.seed(1)
boo <- pbmc2@meta.data$Disease == "Pearson" & (pbmc2@meta.data$seurat_clusters %in% c(0,1,2,4,6,7,8,9,10))
cormat <- cor(pbmc2@meta.data$scale_heteroplasmy[boo],t(data.matrix(pbmc2@assays$chromvar@data[,boo])),
              use = "pairwise.complete")

obs_df <- data.frame(
  TF = sapply(pwm, name),
  cor = cormat[1,]
) %>% dplyr::filter(!is.na(cor)) %>% arrange(desc(cor)) %>% mutate(rank = 1:n(), what = "observed") 

write.table(obs_df, file = "../output/heteroplasmy_Tcell_association_ranking_transcriptionfactoractivities.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

perm_cormat <- cor(sample(pbmc2@meta.data$scale_heteroplasmy[boo]), t(data.matrix(pbmc2@assays$chromvar@data[,boo])),
                   use = "pairwise.complete")

perm_df <- data.frame(
  gene = colnames(perm_cormat),
  cor = perm_cormat[1,]
) %>% dplyr::filter(!is.na(cor)) %>% arrange(desc(cor)) %>% mutate(rank = 1:n(), what = "permuted") 

p1 <- ggplot(obs_df, aes(x = rank, y = cor)) + 
  geom_point(size = 0.2) + scale_color_manual(values = c("black", "firebrick"))  + 
  geom_point(inherit.aes = FALSE, data = perm_df, aes(x = rank, y = cor), color = "lightgrey", size = 0.2) +
  labs(x = "Rank sorted TFs", y = "Correlation") + 
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") + theme_void()
cowplot::ggsave2(p1, file = "../plots/Tcell_heteroplasmy_assoc_TFs.png", width = 3, height = 3, dpi = 500)


#----
# Now do the JUN map
#---


vec_map <- sapply(pwm, name) %>% sort()
# create a Motif object and add it to the assay


# Retain the actual motif matched positions, down to the base pair
DefaultAssay(pbmc2) <- 'peaks'
motif.positions <- matchMotifs(
  pwms = pwm,
  subject = granges(pbmc2),
  out = 'positions',
  genome = 'hg19'
)

# Create a new Motif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  positions = motif.positions,
  pwm = pwm
)

# Store all of this back in the Seurat object
pbmc2 <- SetAssayData(
  object = pbmc2,
  slot = 'motifs',
  new.data = motif
)

pbmc2 <- Footprint(
  object = pbmc2,
  motif.name = c("MA0462.1"),
  genome = BSgenome.Hsapiens.UCSC.hg19,
  in.peaks = TRUE
)

Idents(pbmc2) <- pbmc2@meta.data$Disease
p1 <- PlotFootprint(pbmc2, features = c("MA0462.1")) &
  scale_color_manual(values = c("black", "firebrick"))
cowplot::ggsave2(p1, file = "../plots/JUN_PD_H.pdf", width = 5, height = 4)

Idents(pbmc2) <- pbmc2@meta.data$Patient
p1 <- PlotFootprint(pbmc2, features = c("MA0462.1")) &
  scale_color_manual(values = c("black", "firebrick", "grey", "dodgerblue2"))
cowplot::ggsave2(p1, file = "../plots/JUN_patient_footprint.pdf", width = 5, height = 4)


# Do the Tobias
mat <- pbmc2@assays$peaks@positionEnrichment[[1]][c(-26390, -26391),]
expected <- pbmc2@assays$peaks@positionEnrichment[[1]][c(26390),]
motif <-  pbmc2@assays$peaks@positionEnrichment[[1]][c(26391),]
healthy <- pbmc2@meta.data$Disease == "Healthy"

compute_stat_groups <- function(group_1,group_2){
  
  # Subset out values
  obsvec_g1 <- colMeans(mat[group_1,])
  obsvec_g2 <- colMeans(mat[group_2,])
  
  # Compute normalization values
  norm_value_expected <- mean(expected[c(25:75, 435:485)])
  norm_value_obs_g1 <- mean(obsvec_g1[c(25:75, 435:485)])
  norm_value_obs_g2 <- mean(obsvec_g2[c(25:75, 435:485)])
  
  normvec_g1 <- (obsvec_g1/norm_value_obs_g1) - (expected/norm_value_expected)
  normvec_g2 <- (obsvec_g2/norm_value_obs_g2) - (expected/norm_value_expected)
  
  stat = (mean(normvec_g1[normvec_g1 > 0]) - mean(normvec_g1[normvec_g1 < 0])) -
    (mean(normvec_g2[normvec_g2 > 0]) - mean(normvec_g2[normvec_g2 < 0]))
  
  stat
}

# Assume that group2 is a subset of not group 1
compute_stat <- function(group1, group2){
  obs_all_healthy_stat <- compute_stat_groups(group1, group2)
  
  sapply(1:100, function(i){
    set.seed(i)
    perm1 <- sample(group1)
    
    perm2 <- sample(which(!perm1), size = sum(group2))
    
    compute_stat_groups(perm1, perm2)
  }) -> perm_healthy_stats
  
  (obs_all_healthy_stat-mean(perm_healthy_stats))/sqrt(var(perm_healthy_stats))
}

compute_stat(healthy, !healthy)
compute_stat(healthy, pbmc2@meta.data$Patient == "BCI")
compute_stat(healthy, pbmc2@meta.data$Patient == "CCF")
compute_stat(healthy, pbmc2@meta.data$Patient == "PT3")
