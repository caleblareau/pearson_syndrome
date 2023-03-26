library(stringr)
library(Matrix)
library(data.table)
library(stringi)
library(dplyr)
library(BuenColors)
library(Rsubread)
library(seqinr)

# Run once to set table
reads <- fread("data/deletion_fwd_1.txt.gz")[["V10"]]
mid <- sapply(str_locate_all(reads, "AACGTTATGTAGCAGG"), function(x){
  ifelse(dim(x)[1] == 1, mean(x)-0.5, NA)
})
data.frame(
  mid,
  length = nchar(reads)
) -> rdf
rdf[complete.cases(rdf),] %>%
  filter(length < 2500 & length > 100) -> rdf_filt
rdf_filt %>%
  mutate(pct_before = mid/length*100, pct_after = (length-mid)/length*100) %>%
  arrange((mid)) %>%
  mutate(readid = 1:n(), start = 66, zero = 0, end = length - 65) %>%
  ggplot() +
  geom_rect( aes(xmin=zero, xmax=mid,
                 ymin=readid, ymax=readid),color = "goldenrod") +
  geom_rect( aes(xmin=mid, xmax=length,
                 ymin=readid, ymax=readid),color = "dodgerblue4") +
  pretty_plot(fontsize = 7) + labs(x = "cDNA molecule length", y = "MAS-ISO-seq reads") +
  scale_y_continuous(expand = c(0,0)) +   scale_x_continuous(expand = c(0,0)) -> p1
cowplot::ggsave2(p1, file = "output/read_viz.pdf", width = 2.4, height = 1.5)

