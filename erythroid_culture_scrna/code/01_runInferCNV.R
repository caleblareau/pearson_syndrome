library(infercnv)

process_infer_CNV <- function(idx){
  process_dir = paste0("../data/for_infer_cnv/processed/PearsonInVitro_out", idx, "/")
  out_dir = paste0(process_dir, "inferCNVoutput")
  
  # create the infercnv object
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = paste0(process_dir, "counts_mat.tsv.gz"),
                                      annotations_file = paste0(process_dir, "cells.tsv"),
                                      delim="\t",
                                      gene_order_file = paste0(process_dir, "gene_positions.tsv"),
                                      ref_group_names= "Healthy")
  
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
process_infer_CNV("5")
process_infer_CNV("2")
process_infer_CNV("3")
process_infer_CNV("4")
process_infer_CNV("6")
process_infer_CNV("7")
process_infer_CNV("8")

