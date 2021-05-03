library(data.table)
library(dplyr)
library(BuenColors)

# There are variable cutoffs used to assign the MDS / not
# based on how the CONICs caller emits the probabilities of the
# mixture components; see histograms for justification

import_mgatk_assignments <- function(lib, new_n, day){
  
  # Import the donor assignments
  dt <- fread(paste0("../output/cell_assignments_per_channel/Pearson_",lib,"_assign.tsv")) %>% 
    dplyr::filter(assign == "Pearson") %>%
    data.frame()
  dt$day <- day
  dt$library <- lib
  
  dels <- fread(paste0("../data/dels/del_Pearson_Healthy_",lib,".deletion_heteroplasmy.tsv")) %>%
    dplyr::filter(version == "improved" & deletion == "del10381-15407") %>% dplyr::filter(reads_all >= 20)
  
  # Now import the chromosome 7 
  cnv.dt <- fread(paste0("../../pt3_chr7_del_scatac/output/Pearson_Healthy_",lib,".chr7DelQC.tsv")) %>% data.frame()
  dt2 <- merge(merge(dt, cnv.dt, by.x = "barcode", by.y = "V4"), dels, by.x = "barcode", by.y = "cell_id")
  dt2$barcode <- gsub("-1", paste0("-", as.numeric(new_n)), dt2$barcode)
  dt2[complete.cases(dt2),]
}

ery_df <- rbind(
  import_mgatk_assignments("D6_1", "1", "D06"),
  import_mgatk_assignments("D6_2", "2", "D06"),
  import_mgatk_assignments("D12_1", "3", "D12"),
  import_mgatk_assignments("D12_2", "4", "D12")
)
ggplot( ery_df, aes(x = X7del)) + 
  geom_histogram()+
  facet_wrap(assign ~ day)

ery_df$MDS <- ery_df$X7del < 0.3
ery_df <- ery_df[,c("barcode", "day", "heteroplasmy", "X7del", "MDS")]
dim(ery_df)

if(FALSE){
  write.table(ery_df, file = "../output/Pearson_erythroid_scATAC_assignment_meta.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Process CD34
dels <- fread(paste0("../../cd34_scatac/data/Pearson_CD34_PT3.deletion_heteroplasmy.tsv")) %>%
  dplyr::filter(version == "improved" & deletion == "del10381-15407") %>% dplyr::filter(reads_all >= 20)
cnv.dt <- fread(paste0("../../pt3_chr7_del_scatac/output/Pearson-CD34-PT3.chr7DelQC.tsv")) %>% data.frame()
cd34 <- merge(cnv.dt, dels, by.x = "V4", by.y = "cell_id")
cd34 <- cd34[complete.cases(cd34),]
cd34$MDS <- cd34$X7del < 0.5
cd34$day <- "CD34"
cd34 <- cd34[,c("V4", "day", "heteroplasmy", "X7del", "MDS")]; colnames(cd34)[1] <- "barcode"

# Process BMMNC
dels <- fread(paste0("../../bmmnc_scatac/data/Pearson_BMMNC_PT3del.deletion_heteroplasmy.tsv")) %>%
  dplyr::filter(version == "improved" & deletion == "del10381-15407") %>% dplyr::filter(reads_all >= 20)
cnv.dt <- fread(paste0("../../pt3_chr7_del_scatac/output/Pearson-BMMNC.chr7DelQC.tsv")) %>% data.frame()
bmmnc <- merge(cnv.dt, dels, by.x = "V4", by.y = "cell_id")
bmmnc <- bmmnc[complete.cases(bmmnc),]
qplot(bmmnc$X7del, bins = 25)
bmmnc$MDS <- bmmnc$X7del < 0.3
bmmnc$day <- "BMMNC"
bmmnc <- bmmnc[,c("V4", "day", "heteroplasmy", "X7del", "MDS")]; colnames(bmmnc)[1] <- "barcode"

# Process PBMC
dels <- fread(paste0("../../pbmc_scatac/data/deletion_heteroplasmy/del_PBMC_PT3.deletion_heteroplasmy.tsv")) %>%
  dplyr::filter( deletion == "del10381-15407") %>% dplyr::filter(reads_all >= 20)
cnv.dt <- fread(paste0("../../pt3_chr7_del_scatac/output/Pearson-PBMC.chr7DelQC.tsv")) %>% data.frame()
dels$cell_id <- gsub("-1", "-3", dels$cell_id)
pbmc <- merge(cnv.dt, dels, by.x = "V4", by.y = "cell_id")
pbmc <- pbmc[complete.cases(pbmc),]
qplot(pbmc$X7del, bins = 25)
pbmc$MDS <- pbmc$X7del < 0.3
pbmc$day <- "PBMC"
pbmc <- pbmc[,c("V4", "day", "heteroplasmy", "X7del", "MDS")]; colnames(pbmc)[1] <- "barcode"

all_df <- rbind( ery_df, bmmnc,cd34)
p1 <- ggplot(all_df, aes(x = heteroplasmy, color = day)) +
   scale_x_continuous(limits = c(0,100))+
  stat_ecdf() + pretty_plot(fontsize = 7) + labs(x = "% Heteroplasmy", y = "Cumulative proportion", color = "Library") + L_border() +
  scale_color_manual(values = c("purple3", "dodgerblue2", "orange2", "firebrick"))
p1

cowplot::ggsave2(p1, file = "../plots/heteroplasmy_ecdf.pdf", height = 1.2, width = 2)

all_df <- all_df %>% mutate(nonzerohet=heteroplasmy > 0)
all_df$het_MDS <- paste0(all_df$nonzerohet, "_", all_df$MDS)
all_df %>% group_by(day,het_MDS) %>% summarize(count = n()) %>%
  group_by(day) %>% mutate(prop = count/sum(count)*100) %>%
  ggplot(aes(x = day, y = prop, fill = het_MDS)) +
  geom_bar(stat = "identity") + 
  labs(x = "Predicted MDS clone", y = "% Heteroplasmy > 0")

all_df %>% group_by(day,nonzerohet) %>% summarize(count = n()) %>%
  group_by(day) %>% mutate(prop = count/sum(count)*100) %>%
  ggplot(aes(x = day, y = prop, fill = nonzerohet)) +
  geom_bar(stat = "identity", color = "black") + 
  labs(x = "Cell source", y = "% Heteroplasmy > 0", fill = "") +
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = alpha(c(viridis(9)[1], viridis(9)[5]),0.5)) -> p1

all_df %>% group_by(day,MDS) %>% summarize(count = n()) %>%
  group_by(day) %>% mutate(prop = count/sum(count)*100) %>%
  ggplot(aes(x = day, y = prop, fill = MDS)) +
  geom_bar(stat = "identity", color = "black") + 
  labs(x = "Cell source", y = "% of cells", fill = "") +
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = alpha(c("dodgerblue3", "firebrick"), 0.5)) -> p2

cowplot::ggsave2(cowplot::plot_grid(p1, p2, nrow = 1), 
                 file = "../plots/bargraphs_mds_nonzerohet.pdf", width = 4.5, height = 1.1)

# Day 6 -> 12 fisher test mds
fisher.test(matrix(c(c(681,2248), c(534,1564)), 2, 2))

fisher.test(matrix(c(c(405,2524), c(179,1919)), 2, 2))
