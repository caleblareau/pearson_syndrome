#!/usr/bin/env Rscript 

# load packages
require(data.table)
require(dplyr)
require(tidyr)
require(purrr)
require(optparse)
require(seqinr)

# generate mitochondrial chromosomes with the given deletions

option_list <- list(
    make_option(c("-b", "--breakpoints"), action="store",
                help=" input file breakpoints.tsv"),
    make_option(c("-o", "--outdir"), action="store",
                default="./", help="output directory"))

opt = parse_args(OptionParser(option_list=option_list))

#function to make bedfile from 2 breakpoints of chrM
to_bed <- function(breakpoint_l, breakpoint_r){
    chr=c("chrM", "chrM")
    start=c("1", breakpoint_r)
    end=c(breakpoint_l, 16569)
    out <- data.frame(chr, start, end)
    out
}


#import and process breakpoint data
breakpoints <- fread(opt$b) %>%
    data.frame(.) %>%
    dplyr::select(-V1) %>%
    filter(!is.na(V2)) %>%
    mutate(name = paste("del", V2, V3, sep="_")) %>%
    distinct(.)
names(breakpoints) <- c("bp1", "bp2", "name")

for(i in nrow(breakpoints)){
    out <- to_bed(breakpoints[i,1], breakpoints[i,2])
    file <- file.path(paste0("/broad/sankaranlab/jverboon/pearson/neat/bed", breakpoints[i,3], ".bed"))
    write.table(out, file, col.names = F, row.names = F, quote = F, sep = "\t")
}
