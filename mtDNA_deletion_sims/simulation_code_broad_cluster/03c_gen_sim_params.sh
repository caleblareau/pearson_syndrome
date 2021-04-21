#!/bin/bash

# submit jobs per cell_line

#$ -cwd
#$ -b y
#$ -l h_rt=01:00:00
#$ -N z.generate_sim_paramts
#$ -t 1-1
#$ -tc 50
#$ -l h_vmem=8G 
#$ -binding linear:1

#generate simulation parameters
Rscript r_scripts/sim_params.R \
    -b ../outputs/breakpoints.tsv \
    -o . \
    -p 20 \
    -s 5000 \
    -l 10000 \
    -n 200
