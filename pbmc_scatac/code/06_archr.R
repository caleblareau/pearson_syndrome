library(ArchR)

addArchRGenome("hg19")
ArrowFiles <- createArrowFiles(
  inputFiles = "../../../pearson_mtscatac_large_data_files/input/pearson_3donor_pbmcs/pbmc_3donor_aggr_fragments.tsv.gz",
  sampleNames = "PBMC3donor",
  filterTSS = 1, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# Import reduced 
so <- readRDS("../../../pearson_mtscatac_large_data_files/output/scATAC_labelTransfer_HQcells.rds")

proj <- ArchRProject(
  ArrowFiles = "../../../pearson_mtscatac_large_data_files/output/PBMC3donor.arrow", 
  outputDirectory = "../../../pearson_mtscatac_large_data_files/output/Pearson_PBMC3donor",
  copyArrows = FALSE #This is recommened so that you maintain an unaltered copy for later usage.
)

# Subset cells to those desired
proj <- proj[paste0("PBMC3donor#", colnames(so)),]
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addImputeWeights(proj)

# Get gene score matrix
feat <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
mat <- assays(feat)[["GeneScoreMatrix"]] # replace mat 

# Now impute the values
imputed_gene_score_mat <- imputeMatrix(mat = data.matrix(mat), 
                                       imputeWeights = getImputeWeights(proj),
                                       logFile = "ArchRLogs/ArchR-addImputeWeights-498b101b0608-Date-2020-05-24_Time-13-32-08.log")
colnames(imputed_gene_score_mat) <- gsub("PBMC3donor#", "", colnames(imputed_gene_score_mat))
rownames(imputed_gene_score_mat) <- getFeatures(proj)
saveRDS(imputed_gene_score_mat, file = "../../../pearson_mtscatac_large_data_files/output/24May2020_ArchR_smoothed_activity_scores.rds")



