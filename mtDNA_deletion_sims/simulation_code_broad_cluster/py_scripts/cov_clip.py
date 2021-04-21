import re
import os
import pysam
import sys
from optparse import OptionParser
import numpy as np
import csv

# Parse out data
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process aligned bam files get split reads and coverasge"
opts = OptionParser(usage=usage)
opts.add_option("-i", "--input", help="Filename .bam file to read")

options, arguments = opts.parse_args()
inbam = options.input
bam_in = pysam.AlignmentFile(inbam, "rb")


def process_cigar_for_clip_position(cigar, tuple):

    pos = None

    # Case 1/2: start of read
    if(cigar[1] == "H" or cigar[2] == "H" or
       cigar[1] == "S" or cigar[2] == "S"):
        pos = tuple[0][1]

    # Case 2: end of read
    if(cigar[-1] == "H" or cigar[-1] == "S"):
        pos = tuple[-1][1]

    return(pos)
    

# loop and extract per position clipping counts
clip_pos_count = [0] * 16569

for read in bam_in.fetch('chrM'):
    seq = read.seq
    reverse = read.is_reverse
    #quality = read.query_qualities
    #align_qual_read = read.mapping_quality
    cigar_string = read.cigarstring
    positions = read.get_reference_positions()
    tuple = read.get_aligned_pairs(True)
    if(positions and cigar_string):
        clip_pos = process_cigar_for_clip_position(cigar_string, tuple)
        if clip_pos is not None:
            clip_pos_count[clip_pos] += 1

# get per base coverage; collapse to N bases not ber nucleotide
# and convert to list

cov =  bam_in.count_coverage('chrM', quality_threshold = 0, read_callback = "nofilter")
cov_out = np.add(np.add(cov[0],cov[1]), np.add(cov[2],cov[3])).tolist()

outfile_clip = re.sub(".bam$", ".clip.tsv", inbam)
with open(outfile_clip, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(cov_out,clip_pos_count))
