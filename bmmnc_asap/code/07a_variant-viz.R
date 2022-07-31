library(data.table)
library(dplyr)
library(BuenColors)

# Simple reverse complement function
reverse_complement <- function(s){
  chartr("ATGC","TACG",s)
}

get_enrich_mutation_df <- function(called_variants, ref_all){
  
  # Process 3 digit signature based on letters
  colnames(ref_all) <- c("pos", "ref")
  ref_all$ref <- toupper(ref_all$ref)
  l <- as.character(ref_all$ref)
  
  # Gs happen to be at the first and last position
  ref_all$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))
  
  # Remove Ns
  ref_all <- ref_all[!grepl("N", ref_all$three),]
  
  # Make every possible mutation
  ref_all_long <- rbind(ref_all,ref_all, ref_all,ref_all)
  ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = dim(ref_all)[1])
  ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt,]
  
  # add some meta data
  ref_all_long$variant <- paste0(as.character(ref_all_long$pos), ref_all_long$ref, ">", ref_all_long$alt)
  ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
  ref_all_long$change_rc <- reverse_complement(paste0(ref_all_long$ref, ref_all_long$alt))
  
  # A/G rich strand is "heavy" -- https://en.wikipedia.org/wiki/Heavy_strand
  table(ref_all$ref) # so the reference strand is light (more C/T)
  ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C","T"), "L", "H")
  
  # Change to C/T as ref allele
  ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
  ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", ref_all_long$three, ref_all_long$rc3)
  ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", ref_all_long$change, ref_all_long$change_rc)
  
  # Annotate with called variants
  ref_all_long$called <- ref_all_long$variant %in% called_variants
  
  # Compute changes in expected/observed
  total <- dim(ref_all_long)[1]
  total_called <- sum(ref_all_long$called)
  prop_df <- ref_all_long %>% group_by(three_plot, group_change, strand) %>%
    summarize(observed_prop_called = sum(called)/total_called, expected_prop = n()/total, n = n()) %>%
    mutate(fc_called = observed_prop_called/expected_prop)
  
  prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)
  prop_df 
}

# execute
ref_all <- fread("../../global_functions/chrM_refAllele.txt")
vars_clone <- readRDS("../output/interesting_AFs.rds")
prop_df <- get_enrich_mutation_df(rownames(vars_clone), ref_all)

# Plot the coordinate
pos <- substr(rownames(vars_clone) ,1, nchar(rownames(vars_clone))-3) %>%
  as.character() %>% as.numeric()

sum(pos >= 10381 & pos <= 15407)

fisher.test(
  matrix(c(8, 61, 15407-10381, 16569 - (15407 - 10381)), nrow = 2)
)

# Theme to remove any of the other plot riff raff
xxtheme <-   theme(
  axis.line = element_blank(),
  axis.ticks.y = element_blank(),   
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.margin = unit(c(0, 0, 0, 0), "cm"),
  plot.margin = unit(c(-0.35, -0.35, -0.35, -0.35), "cm")) 


data.frame(
  position = pos, 
  het = rowMeans(vars_clone)*100
) %>% arrange(desc(pos)) 
  
ggplot(aes(x = pos, y = het)) +
  geom_point() + expand_limits(y = c(-3, 1.1), x = c(1,16569)) +
  coord_polar(direction = 1) + labs(x = "", y = "heteroplasmy") + 
  theme(legend.position = "none") + xxtheme  -> P1
cowplot::ggsave2(P1, file = "../plots/het_var_view.pdf", width = 2, height = 2)


# Visualize the nucleotide bias
p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot(fontsize = 8) + L_border() + 
  theme(axis.title.x=element_blank(),
        axis.text.x =element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide (n = 69 mutations)", y = "Substitution Rate (Expected / Observed)")
cowplot::ggsave2(p1, file = "../plots/Pearson_69_mito_signature.pdf", width = 4, height = 2.4)

