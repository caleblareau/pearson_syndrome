library(data.table)
library(BuenColors)
library(dplyr)
library(GenomicRanges)
library(IRanges)


# Import data once
del1 <- fread("chrM_del_13157_15477.st.tsv.gz")
del2 <- fread("chrM_del_9232_13413.st.tsv.gz")
del3 <- fread("chrM_del_8482_13445.st.tsv.gz")
wt <- fread("rCRS.st.tsv.gz")


# Simple function to estimate coverage-based heteroplasmy
estimate_coverage_heteroplasmy_vec <- function(coord1, coord2, vec){
  idx <- 1:length(vec)
  in_del_boo <- (coord1 <= idx & coord2 >= idx)
  in_del <- sum(vec[in_del_boo])/sum(in_del_boo)
  out_del <- sum(vec[!in_del_boo])/sum(!in_del_boo)
  x <- in_del/out_del
  ifelse(x >1, 0, 1-x ) *100
}


# Main function to estimate clipped heteroplasmy for comparison
compute_heteroplasmy_deletion <- function(dt1, dt2, breakpoint_l, breakpoint_r, window_far = 9,  window_near = 24){
  
  read_stats <- rbind(dt1, dt2)
  true_heteroplasmy = dim(dt1)[1]/(dim(dt1)[1] + dim(dt2)[1])
  coverage_vector <- (IRanges(start = read_stats$V1, read_stats$V2)) %>% coverage %>% as.numeric()
  
  coverage_estimate <- estimate_coverage_heteroplasmy_vec(
    breakpoint_l, breakpoint_r,
    coverage_vector
  )
  
  
  # Define parameters
  read_length = 72
  
  # Define windows based on breakpoints
  left_window <- (breakpoint_l-(read_length-window_far)):(breakpoint_l-window_near)
  right_window <- (breakpoint_r+window_near):(breakpoint_r+(read_length-window_far))
  
  # read in stats data and rename variables
  renames <- c("start", "end", "lc", "rc", "clip_pos", "read_name", "barcode", "n_clipped")
  read_stats <- read_stats %>% rename_all(~renames)
  read_stats[is.na(read_stats)] <- 0
  
  # clipped reads
  read_stats <- read_stats %>%
    filter(start %in% left_window | end %in% right_window) %>%
    filter(clip_pos %in% c((breakpoint_l), (breakpoint_r), 0))
  
  # calculate heteroplasmy; ifelse so each sequenced molecule only counts once
  out <- read_stats %>%
    group_by(read_name) %>%
    summarise(total = sum(lc) + sum(rc)) %>%
    mutate(total = ifelse(total > 0, 1, 0)) %>%
    summarise( clip_heteroplasmy = round(sum(100*total)/n(),2),
               reads_del = sum(total), reads_wt = n() - sum(total), reads_all = n()) %>%
    mutate(coverage_heteroplasmy = coverage_estimate, true_heteroplasmy = true_heteroplasmy * 100,
           coverage = mean(coverage_vector), deletion = paste0("del", as.character(breakpoint_l), "-", as.character(breakpoint_r)))  
  out
}

downsample_pe_reads <- function(dt, n, seed_value = 42){
  set.seed(seed_value)
  dt[dt$V6 %in% sample(unique(dt$V6), n),]
}


process_del_set <- function(deldf, bl, br, conditions_df, delid){
  lapply(1:7, function(index_id){
    print(index_id)
    lapply(1:100, function(seed_me){
      compute_heteroplasmy_deletion(
        downsample_pe_reads(deldf, conditions_df[index_id,1], seed_me),
        downsample_pe_reads(wt, conditions_df[index_id,2], seed_me),
        breakpoint_l = bl, breakpoint_r = br)
    }) %>% rbindlist() -> odf
    odf$index_id <- index_id
    odf
  }) %>% rbindlist()  -> out_df
  out_df$del <- delid
  out_df
}

conditions_df_specificity <- data.frame(
  ndelreads = 0,
  nwtreads = c(300,500, 700, 900, 1100, 1300, 1500, 2000)*4,
  ndelreads2 = c(300,500, 700, 900, 1100, 1300, 1500, 2000)*2
)

del1 <- fread("chrM_del_13157_15477.st.tsv.gz")
del2 <- fread("chrM_del_9232_13413.st.tsv.gz")
del3 <- fread("chrM_del_8482_13445.st.tsv.gz")

rrd <- rbind(
  process_del_set(del1, 13157,15477, conditions_df_specificity, "del1"),
  process_del_set(del1, 9232,13413, conditions_df_specificity, "del2"),
  process_del_set(del1, 8482,13445, conditions_df_specificity, "del3")
)

sem <- function(vec){
  sqrt(var(vec))/sqrt(length(vec))
}
rrd %>%
  group_by(index_id, del) %>%
  summarize(mean(coverage),
            mean(true_heteroplasmy),
            mean(clip_heteroplasmy), sem(clip_heteroplasmy),
            mean(coverage_heteroplasmy), sem(coverage_heteroplasmy),
  )

rrd %>%
  group_by(index_id, del) %>%
  summarize(mean(coverage),
            clip = mean(clip_heteroplasmy>=5)*100,  
            cov = mean(coverage_heteroplasmy>=5)*100,  
  ) %>% data.frame() %>%
  reshape2::melt(id.vars = c("index_id", "del", "mean.coverage.")) %>%
  ggplot(aes(x = mean.coverage., y = value, color = variable, shape = del)) +
  geom_point() + geom_line() +
  scale_color_manual(values = c("firebrick", "black")) +
  labs(x = "Simulated coverage", y = "% cells with >5% heteroplasmy (truth = 0)", shape = "", color = "") +
  pretty_plot(fontsize = 7) + L_border() + 
  theme(legend.position = "bottom")  -> p1
cowplot::ggsave2(p1, file = "precision_sim.pdf", width = 1.8, height = 2.2)
