#!/bin/bash

# submit jobs per sim bam
# must have list of bams in dir
# create with find ../sim_data -maxdepth 1 -name '*bam' > bams.txt

#$ -cwd
#$ -b y
#$ -l h_rt=00:10:00
#$ -N z.run_py
#$ -t 1-6
#$ -tc 1000
#$ -l h_vmem=8G 
#$ -binding linear:1

# process bams generated from sim data to get clipped read stats
bams="bams.txt"
bam=$(cat $bams | awk -v ln=$SGE_TASK_ID "NR==ln")
echo $bam
tag="CB"
py="py_scripts/01_process_cell_reads.py"

python $py -i $bam -o $(echo $bam | sed 's/bam/tsv/g')

# clean up files
clip=$(echo $bam | sed 's/bam/clip.tsv/')
#rm $bam
#rm $bam".bai"
#rm $clip
