#!/bin/bash

# submit jobs per cell_line

#$ -cwd
#$ -b y
#$ -l h_rt=01:00:00
#$ -N z.make_del_chr
#$ -t 1-1
#$ -tc 50
#$ -l h_vmem=8G 
#$ -binding linear:1

#get family name from file
Rscript ./r_scripts/make_del_fasta.R \
    -b ../outputs/breakpoints.tsv \
    -f ../sim_data/fasta/rCRS.fasta \
    -r 72

ls -1d ../sim_data/fasta/* | grep -v fa$ > ./fasta.txt
