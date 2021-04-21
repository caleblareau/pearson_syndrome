#!/usr/bin/env Rscript 

# load packages
require(data.table)
require(dplyr)
require(tidyr)
require(purrr)
require(ggplot2)
require(optparse)
require(BuenColors)

#
option_list <- list(
    make_option(c("-i", "--in"), action="store",
                default="./", help="directory containg clip.tsv"),
    make_option(c("-o", "--outdir"), action="store",
                default="./", help="directory for output files"))

opt = parse_args(OptionParser(option_list=option_list))
dir_clip <- opt$i
out_dir <- opt$o
dir.create(out_dir)

# function to read in call .clip.tsv files from filepath
# contain these files and output a df
fread_clipped_pos <- function(files){
    lapply(files,
           function(x) fread(x) %>% 
               as.data.frame() %>%
               filter(clip_pos != 0) %>%
               mutate(sample = gsub(".clip.tsv", "", basename(x))) %>%
               mutate(norm = count/sum(count))) %>%
        bind_rows(., .id = "column_label") %>%
        mutate(sample = factor(sample))
}


# function to determine positions to blacklist based on presence
# in more than 75% of samples
# only makes sense if < 75% of samples have same breakpoints
blacklist_pos <- function(df){
    n_samp <- length(unique(df$sample))
    df %>%
        filter(norm > 0.005) %>%
        filter(!(clip_pos %in%  c(1, 16569))) %>%
        group_by(clip_pos) %>%
        summarise(count = n()) %>%
        ungroup(.) %>%
        filter(count > .75*n_samp) %>%
        .$clip_pos
}

# function to determine breakpoints based on being the 
# split candidates when start/end and 
determine_bp <- function(df, blacklist){
    name <- as.character(df$sample)[1]
    temp <- df %>%
        filter( !(clip_pos %in% blacklist) ) %>%
        arrange(-count) %>%
        head(., n = 4)
    
    # if the other top 4 positions are 3x less than 1st and last
    # of chromosome then unlikely to be deletion
    pos <- c(1, 16569)
    se <- temp %>% filter(clip_pos %in% pos) %>%
        .$count %>% sum(.)
    other <- temp %>%  filter( !(clip_pos %in% pos) ) %>%
        .$count %>% sum(.)
    
    #return breakpoint(highest, not start or end)
    #return breakpoint(highest, not start or end)
    with(temp %>% filter(!(clip_pos %in% pos)), 
         ifelse(se/other < 3,
                return(c(name, min(clip_pos), max(clip_pos))),
                return(c(name, NA, NA)))
         )
}

# function to plot clipped reads along mitochondrial chrom
# with and without blacklisted postions 
plot_chrM <- function(df, blacklist){
    temp <- df %>% dplyr::select(clip_pos, count) %>%
        mutate(no_blacklist = count) %>%
        dplyr::select(-count) %>%
        mutate(blacklist_positions = ifelse(clip_pos %in% c(blacklist, 1, 16569),
                               median(no_blacklist),
                               no_blacklist)) %>%
        gather(., treatment , counts, -clip_pos ) %>%
        mutate(treatment = factor(treatment))
    p <- ggplot(temp, aes(x = clip_pos, y = counts)) +
        geom_line() + pretty_plot()
    p + facet_grid(treatment ~ .)
}
    
#read in clip pos files
clip_files <- list.files(path = dir_clip, pattern = "*.clip.tsv", full.names = T)
clip_all <- fread_clipped_pos(clip_files)

# exploratory plot
p <- ggplot(data = clip_all, aes(x=clip_pos, y = norm,
                                 color = factor(sample))) +
    geom_line()

# determine blacklist sites
blacklist <- blacklist_pos(clip_all)

# empirically defined breakpoints
break_points <- lapply(unique(clip_all$sample),
                       function(x) clip_all %>% filter(sample == x) %>%
                           determine_bp(., blacklist)) %>%
    do.call(rbind, .) %>%
    data.frame(.)
write.table(break_points,
            file = file.path(out_dir, "breakpoints.tsv"), col.names = F,
            row.names = F, quote = F, sep ="\t")

# some plots for caleb to stare at
clip_names <- unique(clip_all$sample)
for(i in clip_names){
    df_temp <- clip_all %>% filter(sample == i)
    out_name <- paste0(i, "_clip_pos.pdf")
    file = file.path(out_dir, out_name)
    print(file)
    pdf(file = file,width=6,height=4)
    print(plot_chrM(df_temp, blacklist))
    dev.off()
}


