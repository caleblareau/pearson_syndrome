library(data.table)
library(Matrix)
library(dplyr)
library(lme4)
library(BuenColors)
library(pbapply)

#  Import raw things
df <- readRDS("../data/3Perason_master_integration_df.rds")
gs.mat <- readRDS("../data/filtered_3pbmcs_pearson.gene.activities.rds")
df_filt <- df %>% filter(coverage > 20)
gs.mat2 <- gs.mat[,df_filt$cell_id]
gs.mat3 = t(t(gs.mat2) / Matrix::colSums(gs.mat2))

gs.mat3@x = log1p(gs.mat3@x * 10000) # equivalent to adding a small pseudocount, but without making matrix dense
rm(gs.mat2); rm(gs.mat)
rownames(gs.mat3) <- make.unique(rownames(gs.mat3))

rowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

df_filt2 <- df_filt%>% filter(predicted.id %in% c("CD14+ Monocytes","CD16+ Monocytes "))
gs.mat4 <- gs.mat3[,as.character(df_filt2 %>% pull(cell_id))]

rv <- rowVars(gs.mat4)
variable_genes <- rv[rv >sort(rv, decreasing = TRUE)[5000]]

pblapply(names(variable_genes), function(gene){
  
  result = tryCatch({
    df_gene <- data.frame(
      gene = as.numeric(gs.mat4[gene,]),
      patient = df_filt2$patient,
      celltype = df_filt2$predicted.id,
      heteroplasmy = df_filt2$heteroplasmy / 100
    )
    mod <- lmer(data = df_gene, gene ~ heteroplasmy  + (1|patient))
    
    data.frame(
      gene = gene,
      t(data.frame(summary(mod)$coefficients[2,])),
      pvalue = car::Anova(mod)[,3], # remove 1 if making cell type a random effect
      row.names = gene
    )
    
  }, error = function(e) {
    data.frame(
      gene = gene,
      Estimate = 0, 
      Std..Error = 0,
      t.value = 0,
      pvalue = 1,
      row.names = gene
    )
  })
  result

}) %>% rbindlist() %>% data.frame() -> results_df
results_df$FDR <- p.adjust(results_df$pvalue)
#saveRDS(results_df, file = "../mixed_model_attempt_2re.rds")

ggplot(results_df, aes(x = Estimate, y = -log10(pvalue))) +
  geom_point()
results_df %>% filter(pvalue < 10e-5)
write.table(data.frame((results_df %>% filter(FDR < 0.01 & Estimate < 0) %>% arrange((FDR)))[,1]), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
