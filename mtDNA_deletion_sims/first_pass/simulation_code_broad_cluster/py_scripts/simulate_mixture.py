#!/usr/bin/env python

import sys
import re
import pysam
import os
from optparse import OptionParser
import random

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Mix two bams together in varying proportions"
opts = OptionParser(usage=usage)
opts.add_option("--bam1", help="First bam file")
opts.add_option("--bam2", help="Second bam file")
opts.add_option("--p1", help="Percent of reads in new bam to generate from first bam")
opts.add_option("--p2", help="Percent of reads in new bam to generate from second bam")
opts.add_option("--ntotal", help="Total number of reads to parse out (approximate / paired end)")
opts.add_option("--outname", help="otuput bam name")


options, arguments = opts.parse_args()

bam1in = options.bam1
bam2in = options.bam2
p1in = float(options.p1)
p2in = float(options.p2)
ntotal = float(options.ntotal)
base_out = options.outname

if(int(p1in + p2in) == 100):
	pass
else:
	sys.exit("Specify p1 / p2 that sums up to 100")

# Pull mtDNA read number
def get_n_mtDNA(bam_in):
	idxs1 = pysam.idxstats(bam_in).split("\n")
	n_mito_reads1 = 0
	for i in idxs1:
		if(i.split("\t")[0] == "chrM"):
			n_mito_reads1 = int(i.split("\t")[2])
	return(n_mito_reads1)

# Divide by 2 to account for paired end reads
n_mito_bam1 = int(get_n_mtDNA(bam1in)/2)
n_mito_bam2 = int(get_n_mtDNA(bam2in)/2)

# Determine the acceptance rate based on parameters
accept_rate_bam1 = p1in/100 * ntotal / n_mito_bam1
accept_rate_bam2 = p2in/100 * ntotal / n_mito_bam2

# Keep sets of reads for paired-end sequencing
reads_bam_1 = {"0"}
reads_bam_2 = {"0"}
reads_bam_1_reject = {"0"}
reads_bam_2_reject = {"0"}

#Establish an output file convention
basename_bam1 = os.path.basename(bam1in).split(".")[0]
basename_bam2 = os.path.basename(bam2in).split(".")[0]
#base_out = "mix_" + str(int(p1in)) + "_" + basename_bam1 + "_" + str(int(p2in)) + "_" + basename_bam2 
raw_out = base_out + ".raw.bam"
final_out = base_out + ".bam"

# Set up file handles
bam1handle = pysam.AlignmentFile(bam1in, "rb")
bam2handle = pysam.AlignmentFile(bam2in, "rb")
out = pysam.AlignmentFile(raw_out, "wb", template = bam1handle)

# Iterate through bam files and write out according to proportion
for read in bam1handle:
	read_name = read.query_name
	if(read_name in reads_bam_1 or (random.uniform(0, 1) < accept_rate_bam1)):
		if(not read_name in reads_bam_1_reject):
			out.write(read)
			reads_bam_1.add(read_name)
		else:
			reads_bam_1_reject.add(read_name)

for read in bam2handle:
	read_name = read.query_name
	if(read_name in reads_bam_2 or (random.uniform(0, 1) < accept_rate_bam2)):
		if(not read_name in reads_bam_2_reject):
			out.write(read)
			reads_bam_2.add(read_name)
		else:
			reads_bam_2_reject.add(read_name)

bam1handle.close()
bam2handle.close()
out.close()

# Final cleanup 
pysam.sort("-o", final_out, raw_out)
os.remove(raw_out)
pysam.index(final_out)

print("total reads (not divided by 2 for paired end):")
print(str(get_n_mtDNA(final_out)))

print("estimated heteroplasmy of bam 2:")
print(str(len(reads_bam_2) / (len(reads_bam_1) + len(reads_bam_2))))
