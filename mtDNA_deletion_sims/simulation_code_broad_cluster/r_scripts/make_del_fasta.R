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
    make_option(c("-f", "--fasta"), action="store",
                default="./", help="input wt fasta"),
    make_option(c("-r", "--read-length"), action="store",
                default="72", help="input wt fasta"))

opt = parse_args(OptionParser(option_list=option_list))
fasta <- opt$fasta

# function to generate said chromosomes
make_del <- function(chrM.fa, breakpointL, breakpointR){
    chrM <- system.file(chrM.fa, package = "seqinr")
    fa <- read.fasta(chrM.fa)[[1]] %>% as.character()
    fa_out <- fa[-((breakpointL + 1):(breakpointR -1))] %>% toupper(.)
    #fa_out <-c(fa_out, fa_out[1:(as.numeric(opt$r)-1)])
    name <- paste("chrM", "del", breakpointL, breakpointR, sep = "_")
    outfile <- paste0(dirname(chrM.fa), "/", name, ".fasta")
    write.fasta(fa_out, name, outfile)
}

#import and process breakpoint data
breakpoints <- fread(opt$b) %>%
    dplyr::select(-V1) %>%
    filter(!is.na(V2)) %>%
    mutate(name = paste("del", V2, V3, sep="_")) %>%
    distinct(.)
names(breakpoints) <- c("bp1", "bp2", "name")

#generate chromsomes
for(i in 1:nrow(breakpoints)){
    with(breakpoints, make_del(fasta, bp1[i], bp2[i]))
}


