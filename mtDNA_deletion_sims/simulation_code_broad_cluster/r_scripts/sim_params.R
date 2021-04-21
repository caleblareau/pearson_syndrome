#!/usr/bin/env Rscript 

# load packages
require(data.table)
require(dplyr)
require(tidyr)
require(purrr)
require(optparse)

# options
option_list <- list(
    make_option(c("-b", "--breakpoints"), action="store",
                help=" input file breakpoints.tsv"),
    make_option(c("-o", "--outdir"), action="store",
                default="./", help="out directory for simulated bams"),
    make_option(c("-p", "--heteroplasmy-increments"), action="store",
                type="numeric", default=25,
                help="number of heteroplasmies to sample"),
    make_option(c("-s", "--min-reads"), action="store",
                type="numeric", default=5000,
                help="min number of reads to simulate per condition"),
    make_option(c("-l", "--max-reads"), action="store",
                type="numeric", default=10000,
                help="max number of reads to simulate per condition"),
    make_option(c("-n", "--number-bams"), action="store",
                type="numeric", default=25,
                help="number of simulated bams"))

opt = parse_args(OptionParser(option_list=option_list))
sim_dir <- opt$outdir
dir.create(sim_dir)

# function to take input parameters + modified breakpoints datafame
# and output deletions and geneterate dataframe containing simulation
# parameters
generate_sim_params <- function(df, het_inc, min_r, max_r, number, outdir){
    set.seed(45)
    heteroplasmy <- sample(0:100, het_inc, replace = F)
    depth <- sample(min_r:max_r, number,  replace = F)
    deletions <- df$name
    l <- list(deletions, depth, heteroplasmy)
    out <- do.call(expand.grid, l)
    names(out) <- c("name", "reads", "h1")
    out <- out %>%
        mutate(h2 = 100 - h1) %>%
        left_join(., df) %>%
        mutate(name = paste(name, reads, h1, h2, "st", sep=".")) %>%
        mutate(bam = paste("chrM", "del", bp1, paste0(bp2, ".st.bam"), sep="_"))
    out
}

#read breakpoints file; keep unique breakpoints
breakpoints <- fread(opt$b) %>%
    dplyr::select(-V1) %>%
    filter(!is.na(V2)) %>%
    mutate(name = paste("del", V2, V3, sep="_")) %>%
    distinct(.)
names(breakpoints) <- c("bp1", "bp2", "name")

# make sim paramaters data.frame
out_params <- generate_sim_params(breakpoints, opt$`heteroplasmy-increments`,
                                  opt$`min-reads`, opt$`max-reads`,
                                  opt$`number-bams`, opt$outdir)

write.table(out_params,
            file = file.path(sim_dir, "sim_params.txt"), col.names = F,
            row.names = F, quote = F, sep ="\t")
