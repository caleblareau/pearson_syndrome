library(data.table)
"%ni%" <- Negate("%in%")
keep_vars <- c( "id_cp_interaction","interacting_pair","partner_a",
                "partner_b","gene_a","gene_b", "secreted", "receptor_a",
                "receptor_b", "annotation_strategy", "is_integrin")


p_bci <- fread("../output/cellphonedb_data/pBCI_cellphonedb_out/pvalues.txt") %>% 
  reshape2::melt(id.vars = keep_vars) %>% filter(value < 0.5) %>%
  mutate(duo = paste0(interacting_pair, "...", variable)) %>% pull(duo)

p_ccf <- fread("../output/cellphonedb_data/pCCF_cellphonedb_out/pvalues.txt") %>% 
  reshape2::melt(id.vars = keep_vars) %>% filter(value < 0.002) %>%
  mutate(duo = paste0(interacting_pair, "...", variable)) %>% pull(duo)

p_pt3 <- fread("../output/cellphonedb_data/pPT3_cellphonedb_out/pvalues.txt") %>% 
  reshape2::melt(id.vars = keep_vars) %>% filter(value < 0.002) %>%
  mutate(duo = paste0(interacting_pair, "...", variable)) %>% pull(duo)

p_h1 <- fread("../output/cellphonedb_data/H1_cellphonedb_out/pvalues.txt") %>% 
  reshape2::melt(id.vars = keep_vars) %>% filter(value < 0.002) %>%
  mutate(duo = paste0(interacting_pair, "...", variable)) %>% pull(duo)

p_h2 <- fread("../output/cellphonedb_data/H2_cellphonedb_out/pvalues.txt") %>% 
  reshape2::melt(id.vars = keep_vars) %>% filter(value < 0.002) %>%
  mutate(duo = paste0(interacting_pair, "...", variable)) %>% pull(duo)

both_healthy <- intersect(p_h1, p_h2)
either_healthy <- unique(c(p_h1, p_h2))

all_pearson <- intersect(intersect(p_bci, p_ccf), p_pt3)
either_pearson <- unique(c(p_bci, p_ccf, p_pt3))

sort(all_pearson[all_pearson %ni% both_healthy])
sort(both_healthy[both_healthy %ni% either_pearson])
