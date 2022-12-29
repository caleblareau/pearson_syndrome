library(data.table)
library(BuenColors)
library(dplyr)
library(GenomicRanges)
library(IRanges)

# Rando funcitons
sem <- function(vec){
  sqrt(var(vec))/sqrt(length(vec))
}


# Import data once
del1 <- fread("../data/chrM_del_13157_15477.st.tsv.gz")
del2 <- fread("../data/chrM_del_9232_13413.st.tsv.gz")
del3 <- fread("../data/chrM_del_8482_13445.st.tsv.gz")
wt <- fread("../data/rCRS.st.tsv.gz")


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

# Master function to iteratively downsample and compute values
process_del_set <- function(deldf, bl, br, conditions_df, delid, n_resample = 100, window_far = 9,  window_near = 24){
  lapply(1:dim(conditions_df)[1], function(index_id){
    print(index_id)
    lapply(1:n_resample, function(seed_me){
      compute_heteroplasmy_deletion(
        downsample_pe_reads(deldf, conditions_df[index_id,1], seed_me),
        downsample_pe_reads(wt, conditions_df[index_id,2], seed_me),
        breakpoint_l = bl, breakpoint_r = br, window_far,  window_near)
    }) %>% rbindlist() -> odf
    odf$index_id <- index_id
    odf
  }) %>% rbindlist()  -> out_df
  out_df$del <- delid
  out_df
}

#----------------------
# Simulate recall of 0s
#----------------------
conditions_df_specificity <- data.frame(
  ndelreads = 0,
  nwtreads = c(300,500, 700, 900, 1100, 1300, 1500, 2000)*4,
  ndelreads2 = c(300,500, 700, 900, 1100, 1300, 1500, 2000)*2
)
rrd <- rbind(
  process_del_set(del1, 13157,15477, conditions_df_specificity, "del1", 100),
  process_del_set(del2, 9232,13413, conditions_df_specificity, "del2", 100),
  process_del_set(del3, 8482,13445, conditions_df_specificity, "del3", 100)
)

# Now visualize
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
cowplot::ggsave2(p1, file = "../output/precision_sim.pdf", width = 1.8, height = 2.2)

#----------------------
# Simulate sensitivity
#----------------------
conditions_df_sensitivity <- data.frame(
  ndelreads = round(c(300,500, 700, 900, 1100, 1300, 1500, 2000)*4/20),
  nwtreads = c(300,500, 700, 900, 1100, 1300, 1500, 2000)*4,
  ndelreads2 = c(300,500, 700, 900, 1100, 1300, 1500, 2000)*2
)
rrd_sens <- rbind(
  process_del_set(del1, 13157,15477, conditions_df_sensitivity, "del1", 100),
  process_del_set(del2, 9232,13413, conditions_df_sensitivity, "del2", 100),
  process_del_set(del3, 8482,13445, conditions_df_sensitivity, "del3", 100)
)

rrd_sens_10 <- rrd_sens

rrd_sens %>%
  group_by(index_id, del) %>%
  summarize(mean(coverage),
            clip = mean(clip_heteroplasmy>=1)*100,  
            cov = mean(coverage_heteroplasmy>=1)*100
  ) %>% data.frame() %>%
  reshape2::melt(id.vars = c("index_id", "del", "mean.coverage.")) %>%
  ggplot(aes(x = mean.coverage., y = value, color = variable, shape = del)) +
  geom_point() + geom_line() +
  scale_color_manual(values = c("firebrick", "black")) +
  labs(x = "Simulated coverage", y = "% cells with deletion detected", shape = "", color = "") +
  pretty_plot(fontsize = 7) + L_border() + 
  theme(legend.position = "none") + scale_y_continuous(limits = c(50, 100)) -> plot_sensistivity
cowplot::ggsave2(plot_sensistivity, file = "../output/sensitivity_plot.pdf", width = 2, height = 2)

#----------------------
# Simulate mean absolute error
#----------------------

# Process data

total_reads <- c(300,500, 700, 900, 1100, 1400, 2000, 3000, 4500, 6000, 9000, 12000)*4

gen_conditions_df2 <- function(prop_del){
  conditions_df <- data.frame(
    ndelreads = total_reads*prop_del,
    nwtreads = total_reads*(1-prop_del)
  )
  conditions_df
}
het_vec <- c(0,1:3*10)/100
me_conditions_df <- lapply(het_vec, function(x){
  gen_conditions_df2(x)
}) %>% rbindlist() %>% data.frame()

ns = 25 # total simulations
rrd_rmse_heteroplasmy <- rbind(
  process_del_set(del1, 13157,15477, me_conditions_df, "del1", n_resample = ns),
  process_del_set(del2, 9232,13413,  me_conditions_df, "del2", n_resample = ns),
  process_del_set(del3, 8482,13445,  me_conditions_df, "del3", n_resample = ns)
)

# Redo index
rrd_rmse_heteroplasmy$idx <- rep(rep(rep(1:length(total_reads), each = ns), length(het_vec)), 3) # 3 is number of deletions

rrd_rmse_heteroplasmy_summary <- rrd_rmse_heteroplasmy %>%
  group_by( del,idx) %>%
  summarize(mean_coverage = round(mean(coverage),0),
            clip = mean(abs(clip_heteroplasmy-true_heteroplasmy)),  
            th = mean(true_heteroplasmy),
            cov = mean(abs(coverage_heteroplasmy-true_heteroplasmy)),
            sem_clip = sem(clip_heteroplasmy),
            sem_cov = sem(coverage_heteroplasmy)
  ) %>% data.frame()
rrd_rmse_heteroplasmy_summary

plot_df2 <-data.frame(
  mean_coverage = c(rrd_rmse_heteroplasmy_summary$mean_coverage,rrd_rmse_heteroplasmy_summary$mean_coverage),
  deletion = c(rrd_rmse_heteroplasmy_summary$del, rrd_rmse_heteroplasmy_summary$del),
  what = c(rep("clip", dim(rrd_rmse_heteroplasmy_summary)[1]),rep("coverage", dim(rrd_rmse_heteroplasmy_summary)[1])),
  mean_abs_error = c(rrd_rmse_heteroplasmy_summary$clip, rrd_rmse_heteroplasmy_summary$cov),
  sem = c(rrd_rmse_heteroplasmy_summary$sem_clip, rrd_rmse_heteroplasmy_summary$sem_cov)
)


px <- plot_df2 %>%
  ggplot(aes(x = mean_coverage, y = mean_abs_error, color = what, fill = what)) +
  geom_point() + geom_line() +
  facet_wrap(~deletion) +
  scale_color_manual(values = c("firebrick", "black")) +
  scale_fill_manual(values = c("firebrick", "black")) +
  labs(x = "Simulated coverage", y = "Mean aboslute error (% heteroplasmy)", shape = "", color = "") +
  pretty_plot(fontsize = 7) +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 8)) +
  scale_x_log10() +
  geom_vline(xintercept = 50)

cowplot::ggsave2(px, file = "../output/error_heteroplsamy.pdf", width = 5, height = 2.4)


