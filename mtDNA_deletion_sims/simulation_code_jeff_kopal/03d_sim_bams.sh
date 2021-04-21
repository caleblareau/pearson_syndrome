#!/bin/bash

# submit jobs per cell_line

#$ -cwd
#$ -b y
#$ -l h_rt=01:00:00
#$ -N z.run_py
#$ -t 1-6
#$ -tc 1000
#$ -l h_vmem=8G 
#$ -binding linear:1

#get family name from file
params="sim_params.txt"
outdir="../sim_data"
outname=$(cat $params | awk -v ln=$SGE_TASK_ID "NR==ln"| cut -f1)
nreads=$(cat $params | awk -v ln=$SGE_TASK_ID "NR==ln"| cut -f2)
p1=$(cat $params | awk -v ln=$SGE_TASK_ID "NR==ln"| cut -f3)
p2=$(cat $params | awk -v ln=$SGE_TASK_ID "NR==ln"| cut -f4)
bam2=$(cat $params | awk -v ln=$SGE_TASK_ID "NR==ln"| cut -f7)
py="./py_scripts/simulate_mixture.py*"


python $py \
--bam1 ../sim_data/in_bams/rCRS.st.bam \
--bam2 ../sim_data/in_bams/$bam2 \
--p1 $p1 \
--p2 $p2 \
--ntotal $nreads \
--outname $outdir/$outname
