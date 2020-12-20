library(infercnv)

process_infer_CNV <- function(idx){
  process_dir = paste0("../data/forinfercnv/PBMCMDS/")
  out_dir = paste0(process_dir, "inferCNVoutput")
  
  # create the infercnv object
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = paste0(process_dir, "counts_mat.tsv.gz"),
                                      annotations_file = paste0(process_dir, "cells.tsv"),
                                      delim="\t",
                                      gene_order_file = paste0(process_dir, "gene_positions.tsv"),
                                      ref_group_names= "healthy")
  
  # perform infercnv operations to reveal cnv signal
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff= 0.1,
                               out_dir = out_dir, 
                               cluster_by_groups=FALSE, 
                               plot_steps=FALSE,
                               denoise=TRUE,
                               HMM=FALSE
  )
  idx
}

process_infer_CNV("1")
