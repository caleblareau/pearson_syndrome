library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(Matrix)

estimate_coverage_heteroplasmy <- function(coord1, coord2, mat){
  idx <- 1:dim(mat)[1]
  in_del_boo <- (coord1 <= idx & coord2 >= idx)
  in_del <- colSums(mat[in_del_boo,])/sum(in_del_boo)
  out_del <- colSums(mat[!in_del_boo,])/sum(!in_del_boo)
  in_del/out_del
}

process_n <- function(n){
  SE <- readRDS(paste0("../../../pearson_large_data_files/input/cellline_mix/Pearson-Mix-",n,"_v12-mtMask_mgatk.rds"))
  summary(SE$depth >= 20)
  filt <- colData(SE)$depth >= 20 & colData(SE)$depth < 100
  assign_df <- fread(paste0("../output/firstpass_assignments_mix",n,".tsv"))
  mat <- assays(SE)[["coverage"]][,assign_df$barcode]
  het_df_cov <- data.frame(
    assign_df[,c(1,2)]
  )
  real_dels <- c("del13157-15477",  "del9232-13413","del8482-13445") 
  
  het_df_cov$del13157_15477 <- estimate_coverage_heteroplasmy(13157,15477,mat)
  het_df_cov$del9232_13413 <- estimate_coverage_heteroplasmy(9232,13413,mat)
  het_df_cov$del8482_13445 <- estimate_coverage_heteroplasmy(8482,13445,mat)
  het_df_cov_melt <- het_df_cov %>% reshape2::melt(id.vars = c("barcode", "fp_classify"))
  het_df_cov_melt %>% dplyr::filter(fp_classify != "unassigned")
}

five_df <- process_n("5")

five_df$value_transformed <- ifelse(five_df$value >1, 0, 1-five_df$value )

five_df <- five_df %>% mutate(
  off_target = (!((fp_classify == "aPD1" & variable == "del13157_15477") | 
                    (fp_classify == "aPD3" & variable == "del8482_13445") |
                    (fp_classify == "aPD2" & variable == "del9232_13413")) & value_transformed > 0)
)
five_df$deletion <- factor(as.character(five_df$variable), c("del13157_15477", "del9232_13413", "del8482_13445"))
p1 <- ggplot(five_df, aes(x = fp_classify, y = value_transformed, color = off_target)) +
  geom_quasirandom(size = 0.1) + facet_wrap(~deletion, ncol = 3) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = c("black", "firebrick")) +
  theme(legend.position = "none") + labs(x = "", y = "Coverage-based deletion heteroplasmy")
cowplot::ggsave2(p1, file = "../output/mix5_scatter_coverage_based.pdf", width = 3.5, height = 1.5)

# Do a side by side comparison
five_df$deletion <- gsub("_", "-", five_df$deletion)
del_het_clip <- fread("../data/all5Sample_dels.deletion_heteroplasmy.tsv") %>% dplyr::filter(version == "improved")
del_het_clip$barcode <- del_het_clip$cell_id
mdf <- merge(five_df, del_het_clip)
mdf_filt <- (mdf %>% dplyr::filter(
  (fp_classify == "aPD1" & variable == "del13157_15477") | 
    (fp_classify == "aPD3" & variable == "del8482_13445") |
    (fp_classify == "aPD2" & variable == "del9232_13413")))

mdf_filt%>%
  ggplot(aes(x = heteroplasmy, y = value_transformed)) + 
  geom_point(size = 0.5) +
  facet_wrap(~fp_classify) +
  pretty_plot(fontsize = 8) +
  geom_abline(intercept = 0, slope = 1/100, linetype = 2) -> pG
cowplot::ggsave2(pG, file = "../output/compare_approaches_scatter.pdf", width = 6, height = 2)

mdf_filt %>% dplyr::filter(fp_classify == "aPD1") -> aPD1
mdf_filt %>% dplyr::filter(fp_classify == "aPD2") -> aPD2
mdf_filt %>% dplyr::filter(fp_classify == "aPD3") -> aPD3

cor(aPD1$heteroplasmy, aPD1$value_transformed)
cor(aPD2$heteroplasmy, aPD2$value_transformed)
cor(aPD3$heteroplasmy, aPD3$value_transformed)


      
      