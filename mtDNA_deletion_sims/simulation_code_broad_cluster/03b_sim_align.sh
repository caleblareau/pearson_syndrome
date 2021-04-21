#!/bin/bash

# submit jobs per cell_line

#$ -cwd
#$ -b y
#$ -l h_rt=01:00:00
#$ -N z.make_sim_fastq
#$ -t 1-9
#$ -tc 50
#$ -l os=RedHat7
#$ -l h_vmem=2G 
#$ -binding linear:4
#$ -pe smp 4

#tools/params
source /broad/software/scripts/useuse
reuse BWA

wgsim="/broad/sankaranlab/jverboon/tools/wgsim/wgsim"
fastq_dir="../sim_data/fastq"
fastas="fasta.txt"
genome_fa="../sim_data/bwa/genome.fa"
fasta=$(cat $fastas | awk -v ln=$SGE_TASK_ID "NR==ln")
outname=$(basename $fasta | sed 's/.fasta//')
bamdir="../sim_data/in_bams"

# simulate data with wgsim
$wgsim $fasta \
    $fastq_dir/$outname"_1.fastq" \
    $fastq_dir/$outname"_2.fastq" \
    -d 500 \
    -s 50 \
    -1 72 \
    -2 72 \
    -R 0.05 \
    -S 1

bwa mem -t 4 -M $genome_fa $fastq_dir/$outname"_1.fastq" $fastq_dir/$outname"_2.fastq" \
    | samtools view -bS - \
    | samtools sort -@ 4 - -o $bamdir/$outname".st.bam"

